/*------------------------------------------------------------------------------
* preceph.c : precise ephemeris and clock functions
*
*          Copyright (C) 2007-2013 by T.TAKASU, All rights reserved.
*
* references :
*     [1] S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-c),
*         12 February, 2007
*     [2] J.Ray, W.Gurtner, RINEX Extensions to Handle Clock Information,
*         27 August, 1998
*     [3] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*     [4] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*         Space Technology Library, 2004
*
* version : $Revision: 1.1 $ $Date: 2014-08-20 12:56:44 $
* history : 2009/01/18 1.0  new
*           2009/01/31 1.1  fix bug on numerical error to read sp3a ephemeris
*           2009/05/15 1.2  support glonass,galileo,qzs
*           2009/12/11 1.3  support wild-card expansion of file path
*           2010/07/21 1.4  added api:
*                               eci2ecef(),sunmoonpos(),peph2pos(),satantoff(),
*                               readdcb()
*                           changed api:
*                               readsp3()
*                           deleted api:
*                               eph2posp()
*           2010/09/09 1.5  fix problem when precise clock outage
*           2011/01/23 1.6  support qzss satellite code
*           2011/09/12 1.7  fix problem on precise clock outage
*                           move sunmmonpos() to rtkcmn.c
*           2011/12/01 1.8  modify api readsp3()
*                           precede later ephemeris if ephemeris is NULL
*                           move eci2ecef() to rtkcmn.c
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

static const char rcsid[]="$Id: preceph.c,v 1.1 2014-08-20 12:56:44 privata Exp $";

#define SQR(x)      ((x)*(x))

#define NMAX        10              /* order of polynomial interpolation */
#define MAXDTE      900.0           /* max time difference to ephem time (s) */
#define EXTERR_CLK  1E-1            /* extrapolation error for clock (m/s) */ // PATCH
#define EXTERR_EPH  1E-3            /* extrapolation error for ephem (m/s^2) */ // PATCH

/* satellite code to satellite system ----------------------------------------*/
static int code2sys(char code)
{
    if (code=='G'||code==' ') return SYS_GPS;
    if (code=='R') return SYS_GLO;
    if (code=='E') return SYS_GAL; /* extension to sp3-c */
    if (code=='J') return SYS_QZS; /* extension to sp3-c */
    if (code=='C') return SYS_CMP; /* extension to sp3-c */
    return SYS_NONE;
}
/* read sp3 header -----------------------------------------------------------*/
static int readsp3h(FILE *fp, gtime_t *time, char *type, int *sats,
                    double *bfact, char *tsys)
{
    int i=0,j,k=0,ns=0,sys,prn;
    char buff[1024];
    
    trace(3,"readsp3h:\n");
    
    *tsys=0;
    bfact[0]=0.0;
    bfact[1]=0.0;
    while (1) {
        if (!fgets(buff,sizeof(buff),fp)) break;
	
	if (!strncmp(buff,"* ",2)) break;
	
        if (i==0) {
            *type=buff[2];
            if (str2time(buff,3,28,time)) return 0;
	    i++;
        }
	
	if (!strncmp(buff,"+ ",2)) {
	  if (!ns) ns=(int)str2num(buff,3,3);
            for (j=0;j<17&&k<ns;j++) {
              sys=code2sys(buff[9+3*j]);
              prn=(int)str2num(buff,10+3*j,2);
              if (k<MAXSAT) sats[k++]=satno(sys,prn);
            }
	}
	
	if (!strncmp(buff,"%c",2)) {
            if (!*tsys) {
	      strncpy(tsys,buff+9,3);
	      tsys[3]='\0';
	    }
	}
	
	if (!strncmp(buff,"%f",2)) {
	  if (!bfact[0] && !bfact[1]) {
            bfact[0]=str2num(buff, 3,10);
            bfact[1]=str2num(buff,14,12);
	    }
	}
	
    }
    
    return ns;
}
/* add precise ephemeris -----------------------------------------------------*/
static int addpeph(nav_t *nav, peph_t *peph)
{
    peph_t *nav_peph;
    
    if (nav->ne>=nav->nemax) {
        nav->nemax+=256;
        if (!(nav_peph=(peph_t *)realloc(nav->peph,sizeof(peph_t)*nav->nemax))) {
            trace(1,"readsp3b malloc error n=%d\n",nav->nemax);
            free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
            return 0;
        }
        nav->peph=nav_peph;
    }
    nav->peph[nav->ne++]=*peph;
    return 1;
}
/* read sp3 body -------------------------------------------------------------*/
static void readsp3b(FILE *fp, char type, int *sats, int ns, double *bfact,
                     char *tsys, int index, int opt, nav_t *nav)
{
    peph_t peph;
    gtime_t time;
    double val,std,base;
    int i,j,sat,sys,prn,n=ns*(type=='P'?1:2),pred_o,pred_c,v;
    char buff[1024];
    
    trace(3,"readsp3b: type=%c ns=%d index=%d opt=%d\n",type,ns,index,opt);
    
    while (fgets(buff,sizeof(buff),fp)) {
        
        if (!strncmp(buff,"EOF",3)) break;
        
        if (buff[0]!='*'||str2time(buff,3,28,&time)) {
            trace(2,"sp3 invalid epoch %31.31s\n",buff);
            continue;
        }
        if (!strcmp(tsys,"UTC")) time=utc2gpst(time); /* utc->gpst */
        peph.time =time;
        peph.index=index;
        
        for (i=0;i<MAXSAT;i++) for (j=0;j<4;j++) {
            peph.pos[i][j]=0.0;
            peph.std[i][j]=0.0f;
        }
        for (i=v=0;i<n&&fgets(buff,sizeof(buff),fp);i++) {
            
            if (strlen(buff)<4||buff[0]!='P') continue;
            
            sys=buff[1]==' '?SYS_GPS:code2sys(buff[1]);
            prn=(int)str2num(buff,2,2);
            if      (sys==SYS_SBS) prn+=100;
            else if (sys==SYS_QZS) prn+=192; /* extension to sp3-c */
            
            if (!(sat=satno(sys,prn))) continue;
            
            pred_c=strlen(buff)>=76&&buff[75]=='P';
            pred_o=strlen(buff)>=80&&buff[79]=='P';
            
            for (j=0;j<4;j++) {
                
                /* read option for predicted value */
                if (j< 3&&opt==1&& pred_o) continue;
                if (j< 3&&opt==2&&!pred_o) continue;
                if (j==3&&opt==1&& pred_c) continue;
                if (j==3&&opt==2&&!pred_c) continue;
                
                val=str2num(buff, 4+j*14,14);
                std=str2num(buff,61+j* 3,j<3?2:3);
                
                if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                    peph.pos[sat-1][j]=val*(j<3?1000.0:1E-6);
                    v=1; /* valid epoch */
                }
                if ((base=bfact[j<3?0:1])>0.0) {
                    peph.std[sat-1][j]=(float)(pow(base,std)*(j<3?1E-3:1E-12));
                }
            }
        }
        if (v) {
            if (!addpeph(nav,&peph)) return;
        }
    }
}
/* compare precise ephemeris -------------------------------------------------*/
static int cmppeph(const void *p1, const void *p2)
{
    peph_t *q1=(peph_t *)p1,*q2=(peph_t *)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}
/* combine precise ephemeris -------------------------------------------------*/
static void combpeph(nav_t *nav)
{
    int i,j,k,m;
    
    trace(3,"combpeph: ne=%d\n",nav->ne);
    
    qsort(nav->peph,nav->ne,sizeof(peph_t),cmppeph);
    
    for (i=0,j=1;j<nav->ne;j++) {
        
        if (fabs(timediff(nav->peph[i].time,nav->peph[j].time))<1E-9) {
            
            for (k=0;k<MAXSAT;k++) {
#if 0
                if (norm(nav->peph[j].pos[k],4)<=0.0) continue;
#endif
                for (m=0;m<4;m++) nav->peph[i].pos[k][m]=nav->peph[j].pos[k][m];
                for (m=0;m<4;m++) nav->peph[i].std[k][m]=nav->peph[j].std[k][m];
            }
        }
        else if (++i<j) nav->peph[i]=nav->peph[j];
    }
    nav->ne=i+1;
    
    trace(4,"combpeph: ne=%d\n",nav->ne);
}
/* read sp3 precise ephemeris file ---------------------------------------------
* read sp3 precise ephemeris/clock files and set them to navigation data
* args   : char   *file       I   sp3-c precise ephemeris file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*          int    opt         I   options (1: only observed, 2: only predicted)
* return : none
* notes  : see ref [1]
*          precise ephemeris is appended and combined
*          nav->peph and nav->ne must by properly initialized before calling the
*          function
*          only files with extensions of .sp3, .SP3, .eph* and .EPH* are read
*-----------------------------------------------------------------------------*/
extern void readsp3(const char *file, nav_t *nav, int opt)
{
    FILE *fp;
    gtime_t time={0};
    double bfact[2]={0};
    int i,j,n,ns,sats[MAXSAT]={0};
    char *efiles[MAXEXFILE],*ext,type=' ',tsys[4]="";
    
    trace(3,"readpephs: file=%s\n",file);
    
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return;
        }
    }
    /* expand wild card in file path */
    n=expath(file,efiles,MAXEXFILE);
    
    for (i=j=0;i<n;i++) {
        if (!(ext=strrchr(efiles[i],'.'))) continue;
        
        if (strcmp (ext,".sp3"  )&&strcmp (ext,".SP3"  )&&
            strncmp(ext,".eph",4)&&strncmp(ext,".EPH",4)) continue;
        
        if (!(fp=fopen(efiles[i],"r"))) {
            trace(2,"sp3 file open error %s\n",efiles[i]);
            continue;
        }
        /* read sp3 header */
        ns=readsp3h(fp,&time,&type,sats,bfact,tsys);
        
        /* read sp3 body */
        rewind(fp);
        readsp3b(fp,type,sats,ns,bfact,tsys,j++,opt,nav);
        
        fclose(fp);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);
    
    /* combine precise ephemeris */
    if (nav->ne>0) combpeph(nav);
}
/* read satellite antenna parameters -------------------------------------------
* read satellite antenna parameters
* args   : char   *file       I   antenna parameter file
*          gtime_t time       I   time
*          nav_t  *nav        IO  navigation data
* return : status (1:ok,0:error)
* notes  : only support antex format for the antenna parameter file
*-----------------------------------------------------------------------------*/
extern int readsap(const char *file, gtime_t time, nav_t *nav)
{
    pcvs_t pcvs={0};
    pcv_t pcv0={0},*pcv;
    int i;
    
    trace(3,"readsap : file=%s time=%s\n",file,time_str(time,0));
    
    if (!readpcv(file,&pcvs)) return 0;
    
    for (i=0;i<MAXSAT;i++) {
        pcv=searchpcv(i+1,"",time,&pcvs);
        nav->pcvs[i]=pcv?*pcv:pcv0;
    }
    free(pcvs.pcv);
    return 1;
}
/* read dcb parameters file --------------------------------------------------*/
static int readdcbf(const char *file, nav_t *nav)
{
    FILE *fp;
    double cbias;
    int sat,type=0;
    char buff[256];
    
    trace(3,"readdcbf: file=%s\n",file);
    
    if (!(fp=fopen(file,"r"))) {
        trace(2,"dcb parameters file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        
        if      (strstr(buff,"DIFFERENTIAL (P1-P2) CODE BIASES")) type=1;
        else if (strstr(buff,"DIFFERENTIAL (P1-C1) CODE BIASES")) type=2;
        else if (strstr(buff,"DIFFERENTIAL (P2-C2) CODE BIASES")) type=3;
        
        if (!type) continue;
        
        if (!(sat=satid2no(buff))||(cbias=str2num(buff,26,9))==0.0) continue;
        
        nav->cbias[sat-1][type-1]=cbias*1E-9*CLIGHT; /* ns -> m */
    }
    fclose(fp);
    
    return 1;
}
/* read dcb parameters ---------------------------------------------------------
* read differential code bias (dcb) parameters
* args   : char   *file       I   dcb parameters file (wild-card * expanded)
*          nav_t  *nav        IO  navigation data
* return : status (1:ok,0:error)
* notes  : currently only p1-c1 bias of code *.dcb file
*-----------------------------------------------------------------------------*/
extern int readdcb(const char *file, nav_t *nav)
{
    int i,j,n;
    char *efiles[MAXEXFILE]={0};
    
    trace(3,"readdcb : file=%s\n",file);
    
    for (i=0;i<MAXSAT;i++) for (j=0;j<3;j++) {
        nav->cbias[i][j]=0.0;
    }
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return 0;
        }
    }
    n=expath(file,efiles,MAXEXFILE);
    
    for (i=0;i<n;i++) {
        readdcbf(efiles[i],nav);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);
    
    return 1;
}
/* polynomial interpolation by Neville's algorithm ---------------------------*/
static double interppol(const double *x, double *y, int n)
{
    int i,j;
    
    for (j=1;j<n;j++) {
        for (i=0;i<n-j;i++) {
            y[i]=(x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
        }
    }
    return y[0];
}
/* satellite position by precise ephemeris -----------------------------------*/
static int pephpos(gtime_t time, int sat, const nav_t *nav, double *rs,
                   double *dts, double *vare, double *varc)
{
    double t[NMAX+1],p[3][NMAX+1],c[2],*pos,std=0.0,s[3],sinl,cosl;
    int i,j,k,index;
    
    trace(4,"pephpos : time=%s sat=%2d\n",time_str(time,3),sat);
    
    rs[0]=rs[1]=rs[2]=dts[0]=0.0;
    
    if (nav->ne<NMAX+1||
        timediff(time,nav->peph[0].time)<-MAXDTE||
        timediff(time,nav->peph[nav->ne-1].time)>MAXDTE) {
        trace(2,"no prec ephem %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    /* binary search */
    for (i=0,j=nav->ne-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->peph[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;
    
    /* polynomial interpolation for orbit */
    i=index-(NMAX+1)/2;
    if (i<0) i=0; else if (i+NMAX>=nav->ne) i=nav->ne-NMAX-1;
    
    for (j=0;j<=NMAX;j++) {
        t[j]=timediff(nav->peph[i+j].time,time);
        if (norm(nav->peph[i+j].pos[sat-1],3)<=0.0) {
            trace(2,"prec ephem outage %s sat=%2d\n",time_str(time,0),sat);
            return 0;
        }
    }
    for (j=0;j<=NMAX;j++) {
        pos=nav->peph[i+j].pos[sat-1];
#if 0
        p[0][j]=pos[0];
        p[1][j]=pos[1];
#else
        /* correciton for earh rotation ver.2.4.0 */
        sinl=sin(OMGE*t[j]);
        cosl=cos(OMGE*t[j]);
        p[0][j]=cosl*pos[0]-sinl*pos[1];
        p[1][j]=sinl*pos[0]+cosl*pos[1];
#endif
        p[2][j]=pos[2];
    }
    for (i=0;i<3;i++) {
        rs[i]=interppol(t,p[i],NMAX+1);
    }
    if (vare) {
        for (i=0;i<3;i++) s[i]=nav->peph[index].std[sat-1][i];
        std=norm(s,3);
        
        /* extrapolation error for orbit */
        if      (t[0   ]>0.0) std+=EXTERR_EPH*SQR(t[0   ])/2.0;
        else if (t[NMAX]<0.0) std+=EXTERR_EPH*SQR(t[NMAX])/2.0;
        *vare=SQR(std);
    }
    /* linear interpolation for clock */
    t[0]=timediff(time,nav->peph[index  ].time);
    t[1]=timediff(time,nav->peph[index+1].time);
    c[0]=nav->peph[index  ].pos[sat-1][3];
    c[1]=nav->peph[index+1].pos[sat-1][3];
    
    if (t[0]<=0.0) {
        if ((dts[0]=c[0])!=0.0) {
            std=nav->peph[index].std[sat-1][3]*CLIGHT-EXTERR_CLK*t[0];
        }
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])!=0.0) {
            std=nav->peph[index+1].std[sat-1][3]*CLIGHT+EXTERR_CLK*t[1];
        }
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->peph[index+i].std[sat-1][3]*CLIGHT+EXTERR_CLK*fabs(t[i]);    // PATCH
    }
    else {
        dts[0]=0.0;
    }
    if (varc) *varc=SQR(std);
    return 1;
}
/* satellite clock by precise clock ------------------------------------------*/
static int pephclk(gtime_t time, int sat, const nav_t *nav, double *dts,
                   double *varc)
{
    double t[2],c[2],std;
    int i,j,k,index;
    
    trace(4,"pephclk : time=%s sat=%2d\n",time_str(time,3),sat);
    
    if (nav->nc<2||
        timediff(time,nav->pclk[0].time)<-MAXDTE||
        timediff(time,nav->pclk[nav->nc-1].time)>MAXDTE) {
        trace(3,"no prec clock %s sat=%2d\n",time_str(time,0),sat);
        return 1;
    }
    /* binary search */
    for (i=0,j=nav->nc-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->pclk[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;
    
    /* linear interpolation for clock */
    t[0]=timediff(time,nav->pclk[index  ].time);
    t[1]=timediff(time,nav->pclk[index+1].time);
    c[0]=nav->pclk[index  ].clk[sat-1][0];
    c[1]=nav->pclk[index+1].clk[sat-1][0];
    
    if (t[0]<=0.0) {
        if ((dts[0]=c[0])==0.0) return 0;
        std=nav->pclk[index].std[sat-1][0]*CLIGHT-EXTERR_CLK*t[0];
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])==0.0) return 0;
        std=nav->pclk[index+1].std[sat-1][0]*CLIGHT+EXTERR_CLK*t[1];
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->pclk[index+i].std[sat-1][0]*CLIGHT+EXTERR_CLK*fabs(t[i]);    // PATCH
    }
    else {
        trace(3,"prec clock outage %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    if (varc) *varc=SQR(std);
    return 1;
}
/* satellite antenna phase center offset ---------------------------------------
* compute satellite antenna phase center offset in ecef
* args   : gtime_t time       I   time (gpst)
*          double *rs         I   satellite position and velocity (ecef)
*          int    sat         I   satellite number
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          pcv_t  *pcv        I   satellite antenna parameter
*          double *dant       I   satellite antenna phase center offset (ecef)
*                                 {dx,dy,dz} (m)
* return : none
*-----------------------------------------------------------------------------*/
extern void satantoff(gtime_t time, const double *rs, int sat, const pcv_t *pcv,
                      double *dant)
{
    double ex[3],ey[3],ez[3],es[3],r[3],rsun[3],gmst,erpv[5]={0};
    int i;
    
    trace(4,"satantoff: time=%s\n",time_str(time,3));
    
    /* sun position in ecef */
    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,&gmst);
    
    /* psi angle (PATCH) */
    double usun[3], vi[3], n[3], sinbeta, tanbeta, u[3], usat[3], sinmu, psi;
    double x0[3], y0[3], exn[3], eyn[3], zz[3], xx[3];
    for (i=0;i<3;i++) r[i]=-rs[i];
    if (!normv3(r,ez)) return;
    normv3(rsun,usun);
    vi[0]=rs[3]-OMGE*rs[1];
    vi[1]=rs[4]+OMGE*rs[0];
    vi[2]=rs[5];
    cross3(rs,vi,r);
    if (!normv3(r,n)) return;
    sinbeta=dot(n,usun,3);
    tanbeta=sinbeta/sqrt(1-sinbeta*sinbeta);
    cross3(usun,n,r);
    if (!normv3(r,u)) return;
    if (!normv3(rs,usat)) return;
    sinmu=dot(u,usat,3);
    psi=atan2(-tanbeta, sinmu);
    
    /* Beidou psi specific cases */
    int sys=satsys(sat, NULL);
    if (sys == SYS_CMP) {
      double inc=acos(n[2]);
      if (fabs(inc) < 0.1) {                             /* GEO */
        psi=0.0;
      }
      else {                                             /* IGSO & MEO */
        double beta=asin(sinbeta);
	if (fabs(beta*R2D)<4.0)                          /* Normal law */
	  psi=0.0;
      }
    }
    
    /* ex and ey */
    if (!normv3(vi,x0)) return;
    cross3(ez,x0,r);
    if (!normv3(r,y0)) return;
    for (i=0;i<3;i++) {
      ex[i] = x0[i]*cos(psi)+y0[i]*sin(psi);
      ey[i] = y0[i]*cos(psi)-x0[i]*sin(psi);
    }
    
    for (i=0;i<3;i++) { /* use L1 value */
        dant[i]=pcv->off[0][0]*ex[i]+pcv->off[0][1]*ey[i]+pcv->off[0][2]*ez[i];
    }
}
/* satellite position/clock by precise ephemeris/clock -------------------------
* compute satellite position/clock with precise ephemeris/clock
* args   : gtime_t time       I   time (gpst)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          int    opt         I   sat postion option
*                                 (0: center of mass, 1: antenna phase center)
*          double *rs         O   sat position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts        O   sat clock {bias,drift} (s|s/s)
*          double *var        IO  sat position and clock error variance (m)
*                                 (NULL: no output)
* return : status (1:ok,0:error or data outage)
* notes  : clock includes relativistic correction but does not contain code bias
*          before calling the function, nav->peph, nav->ne, nav->pclk and
*          nav->nc must be set by calling readsp3(), readrnx() or readrnxt()
*          if precise clocks are not set, clocks in sp3 are used instead
*-----------------------------------------------------------------------------*/
extern int peph2pos(gtime_t time, int sat, const nav_t *nav, int opt,
                    double *rs, double *dts, double *var)
{
    double rss[6],rst[3],dtss[1],dtst[1],dant[3]={0},vare=0.0,varc=0.0,tt=1E-3;
    int i;
    
    trace(4,"peph2pos: time=%s sat=%2d opt=%d\n",time_str(time,3),sat,opt);
    
    if (sat<=0||MAXSAT<sat) return 0;
    
    /* satellite position and clock bias */
    if (!pephpos(time,sat,nav,rss,dtss,&vare,&varc)||
        !pephclk(time,sat,nav,dtss,&varc)) return 0;
    
    time=timeadd(time,tt);
    if (!pephpos(time,sat,nav,rst,dtst,NULL,NULL)||
        !pephclk(time,sat,nav,dtst,NULL)) return 0;
    
    /* compute speed for antenna offset (PATCH) */
    for (i=0;i<3;i++) {
        rss[i+3]=(rst[i]-rss[i])/tt;
    }
    /* satellite antenna offset correction */
    if (opt) {
        satantoff(time,rss,sat,nav->pcvs+sat-1,dant);
    }
    for (i=0;i<3;i++) {
        rs[i  ]=rss[i]+dant[i];
        rs[i+3]=(rst[i]-rss[i])/tt;
    }
    /* relativistic effect correction */
    if (dtss[0]!=0.0) {
        dts[0]=dtss[0]-2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
        dts[1]=(dtst[0]-dtss[0])/tt;
    }
    else { /* no precise clock */
        dts[0]=dts[1]=0.0;
    }
//    if (var) *var=vare+varc;
    if (var) *var=varc;                 // (PATCH -> only clock)
    
    return 1;
}
