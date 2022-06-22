#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>
#include "rtklib.h"

using namespace std;

#define NB_FREQ 5

//////////////////////////////////////////////////////////////////////////////////////////////
static int biasstr2time(const char *s, int i, int n, gtime_t *t)
{
    double ep[6]={}; ep[1]=1.0; ep[2]=1.0;
    double day, sec;
    char str[256],*p=str;
    if (i<0||(int)strlen(s)<i||(int)sizeof(str)-1<i) return -1;
    for (s+=i;*s&&--n>=0;) *p++=*s++; *p='\0';
    if (sscanf(str,"%lf:%lf:%lf",ep,&day,&sec)<3)
        return -1;
    if (ep[0]<100.0) ep[0]+=ep[0]<80.0?2000.0:1900.0;
    *t=timeadd(epoch2time(ep), 86400.0*(day-1.0)+sec);
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Biases file functions (rinex format)
typedef struct {        /* precise bias type */
    gtime_t t1, t2;     /* time (GPST) */
    int index;          /* clock index for multiple files */
    int code;           /* bias code */
    int range;          /* range/phase bias */
    int sat;            /* satellite */
    double bias;        /* satellite bias (ns) */
} pbias_t;

typedef struct {        /* bias data type */
    int nb,nbmax;       /* number of biases */
    pbias_t *pbias;     /* biases */
} bias_t;


typedef struct {
  double code[MAXSAT][MAXCODE];
  double phase[MAXSAT][MAXCODE];
} pbiases_t;

typedef struct {
    gtime_t tmin, tmax;
    double dt;
    pbiases_t *pbiases;
} biases_t;

//////////////////////////////////////////////////////////////////////////////////////////////
static void generateBiases(bias_t *bias, biases_t *biases)
{
int i, nb, ii;
gtime_t tmin={}, tmax={};
double dt=0;
for (i=0;i<bias->nb;i++) {
  if (i==0) {
    tmin=bias->pbias[i].t1;
    tmax=bias->pbias[i].t2;
    dt=timediff(tmax,tmin);
    }
  if (timediff(bias->pbias[i].t1,tmin)<0.0) tmin=bias->pbias[i].t1;
  if (timediff(bias->pbias[i].t2,tmax)>0.0) tmax=bias->pbias[i].t2;
  if (timediff(bias->pbias[i].t2,bias->pbias[i].t1)<dt) dt=timediff(bias->pbias[i].t2,bias->pbias[i].t1);
  }
char st1[20], st2[20];
time2str(tmin, st1, 0);
time2str(tmax, st2, 0);

if (!dt) {
  biases->dt=0.0;
  biases->pbiases=NULL;
  return;
  }
biases->dt=dt;
biases->tmin=tmin;
biases->tmax=tmax;
nb=timediff(tmax,tmin)/dt;
biases->pbiases=(pbiases_t *)calloc(nb, sizeof(pbiases_t));
if (biases->pbiases==NULL) {
  fprintf(stderr, "Malloc error\n");
  exit(1);
  }

for (i=0;i<bias->nb;i++) {
  int i1=(int)(timediff(bias->pbias[i].t1,tmin)/dt);
  int i2=(int)(timediff(bias->pbias[i].t2,tmin)/dt);
  for (ii=i1;ii<i2;ii++) {
    if (bias->pbias[i].range)
      biases->pbiases[ii].code[bias->pbias[i].sat-1][bias->pbias[i].code]=bias->pbias[i].bias;
    else
      biases->pbiases[ii].phase[bias->pbias[i].sat-1][bias->pbias[i].code]=bias->pbias[i].bias;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
static int readrnxbias(FILE *fp, int index, bias_t *bias)
{
    static const char *sc[MAXCODE]={"1C", "1P", "1W", "1Y", "1M", "1N", "1S", "1L", "1E", "1A", "1B", "1X", "1Z",
                       "2C", "2D", "2S", "2L", "2X", "2P", "2W", "2Y", "2M", "2N", "5I", "5Q", "5X",
                       "7I", "7Q", "7X", "6A", "6B", "6C", "6X", "6Z", "6S", "6L", "8I", "8Q", "8X",
                       "2I", "2Q", "6I", "6Q", "3I", "3Q", "3X", "1I", "1Q", 
                       "5A", "5B", "5C", "9A", "9B", "9C", "9X"};

    pbias_t *nav_pbias;
    int i, code, range;
    char buff[200];
    
    if (!bias) return 0;
    
    while (fgets(buff,sizeof(buff),fp)) {
    
	if ((!strncmp(buff+1, "OSB", 3)) && (!strncmp(buff+65, "ns", 2))) {
	  int sat=satid2no(buff+11);
	  int range=(buff[25]=='C'?1:0);
	  for (i=0;i<MAXCODE;i++)
	    if (!strncmp(sc[i], buff+26, 2))
	      code=i;
	  gtime_t t1, t2;
	  char st1[20], st2[20];
	  if (biasstr2time(buff, 35, 14, &t1)) continue;
	  if (biasstr2time(buff, 50, 14, &t2)) continue;
	  time2str(t1, st1, 0);
	  time2str(t2, st2, 0);
	  double value=atof(buff+70);
          if (bias->nb>=bias->nbmax) {
            bias->nbmax+=1024;
            if (!(nav_pbias=(pbias_t *)realloc(bias->pbias,sizeof(pbias_t)*(bias->nbmax)))) {
                free(bias->pbias); bias->pbias=NULL; bias->nb=bias->nbmax=0;
                return -1;
            }
            bias->pbias=nav_pbias;
          }
          bias->nb++;
          bias->pbias[bias->nb-1].t1=t1;
          bias->pbias[bias->nb-1].t2=t2;
          bias->pbias[bias->nb-1].index=index;
          bias->pbias[bias->nb-1].sat=sat;
          bias->pbias[bias->nb-1].code=code;
          bias->pbias[bias->nb-1].range=range;
          bias->pbias[bias->nb-1].bias=value;
	}
    }
    return bias->nb>0;
}

//////////////////////////////////////////////////////////////////////////////////////////////
static int readrnxfilebias(const char *file, int index, bias_t *bias)
{
    FILE *fp;
    int cstat,stat;
    char tmpfile[1024];
    
    /* uncompress file */
    if ((cstat=uncompress(file,tmpfile))<0) {
        return 0;
    }
    if (!(fp=fopen(cstat?tmpfile:file,"r"))) {
        return 0;
    }
    /* read rinex file */
    stat=readrnxbias(fp,index,bias);
    
    fclose(fp);
    
    /* delete temporary file */
    if (cstat) remove(tmpfile);
    
    return stat;
}

//////////////////////////////////////////////////////////////////////////////////////////////
static int readrnxbias(const char *file, bias_t *bias)
{
    int i,n,index=0,stat=1;
    char *files[MAXEXFILE]={0};
    
    for (i=0;i<MAXEXFILE;i++) {
        if (!(files[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(files[i]); return 0;
        }
    }
    /* expand wild-card */
    n=expath(file,files,MAXEXFILE);
    
    /* read rinex clock files */
    for (i=0;i<n;i++) {
        if (readrnxfilebias(files[i],index++,bias)) {
            continue;
        }
        stat=0;
        break;
    }
    for (i=0;i<MAXEXFILE;i++) free(files[i]);
    
    if (!stat) return 0;
    
    return bias->nb;
}

//////////////////////////////////////////////////////////////////////////////////////////////
static int precbias(gtime_t time, int sat, const biases_t *biases, double code[MAXCODE], double phase[MAXCODE])
{
  int i,j;
    
  if (!biases->dt) return 0;
  if (timediff(time,biases->tmin)<0.0) return 0;
  if (timediff(time,biases->tmax)>0.0) return 0;
  i=(int)(timediff(time,biases->tmin)/biases->dt);
  for (j=0;j<MAXCODE;j++) {
    code[j] = biases->pbiases[i].code[sat-1][j];
    phase[j] = biases->pbiases[i].phase[sat-1][j];
    }

  return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////
static int readchannel(const char *file, int *channel)
{
  FILE *fp;
  int i, sat, chnl;
  char buff[50];

  if (!(fp = fopen(file, "r"))) {
    fprintf(stderr, "WARNING: Glonass channel file unavailable (-chan)\n");
    for (i=1;i<MAXPRNGLO;i++)
      channel[i-1] = -999;
    return 0;
  }

  while (fgets(buff, sizeof(buff), fp)) {

    if (!strncmp(buff, "%", 1))
      continue;

    sscanf(buff, "%d%d", &sat, &chnl);
    if ((sat < 0) || (sat > MAXPRNGLO)) {
      fprintf(stderr, "channel file error\n");
      exit(1);
    }
    channel[sat-1] = chnl;
  }
  fclose(fp);

  return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////
static void output(char *id, int day, double sec, double code[NB_FREQ], double phase[NB_FREQ], double doppler[NB_FREQ], int channel,
                   double rs[6], double dts, double iono, int ionoSource, double tropo, double codeBias[NB_FREQ], double phaseBias[NB_FREQ], int n1, int nw, double yaw, int discontinuity)
{
printf("%d %s %d %15.3lf \
%15.3lf %15.3lf %15.3lf %15.3lf %15.3lf \
%15.3lf %15.3lf %15.3lf %15.3lf %15.3lf \
%15.3lf %15.3lf %15.3lf %15.3lf %15.3lf \
%d \
%15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf \
%lf %d %lf %lf \
%lf %lf %lf %lf %lf \
%lf %lf %lf %lf %lf \
%d %d %d\n", 
      1, id, day, sec,
      code[0], code[1], code[2], code[3], code[4],
      phase[0], phase[1], phase[2], phase[3], phase[4],
      doppler[0], doppler[1], doppler[2], doppler[3], doppler[4],
      channel,
      rs[0], rs[1], rs[2], dts, rs[3], rs[4], rs[5],
      iono, ionoSource, tropo, yaw,
      codeBias[0], codeBias[1], codeBias[2], codeBias[3], codeBias[4],
      phaseBias[0], phaseBias[1], phaseBias[2], phaseBias[3], phaseBias[4],
      n1, nw, discontinuity);
}

//////////////////////////////////////////////////////////////////////////////////////////////
static void computeBiases(nav_t *pnav, int pcb, double cbias[MAXCODE], double pbias[MAXCODE], int sat, int p1, double codeBias[NB_FREQ], double phaseBias[NB_FREQ])
{
#define GPS_GAMMA (FREQ1/FREQ2*FREQ1/FREQ2)
#define GPS_LAM1 (CLIGHT/FREQ1)
#define GPS_LAM2 (CLIGHT/FREQ2)
int i;
for (i=0;i<NB_FREQ;i++)
  {
  codeBias[i]=0.0;
  phaseBias[i]=0.0;
  }
if (pcb)
  {
  for (i=0;i<MAXCODE;i++) {
    //!!!!!!! PATCH !!!!!!!!
    //if ((satsys(sat, NULL) == SYS_GAL) && ((i+1) == CODE_L1B) && cbias[i])        codeBias[0]=-CLIGHT/1E9*cbias[i];
    //if ((satsys(sat, NULL) == SYS_GAL) && ((i+1) == CODE_L5Q) && cbias[i])        codeBias[3]=-CLIGHT/1E9*cbias[i];
    //if ((satsys(sat, NULL) == SYS_GAL) && ((i+1) == CODE_L7Q) && cbias[i])        codeBias[4]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_GAL) && code2obs(i+1, NULL)[0] == '1' && cbias[i])        codeBias[0]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_GAL) && code2obs(i+1, NULL)[0] == '5' && cbias[i])        codeBias[3]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_GAL) && code2obs(i+1, NULL)[0] == '7' && cbias[i])        codeBias[4]=-CLIGHT/1E9*cbias[i];
    //!!!!!!! PATCH !!!!!!!!
    if ((satsys(sat, NULL) == SYS_GPS) && ((i+1) == CODE_L1C) && !p1 && cbias[i]) codeBias[0]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_GPS) && ((i+1) == CODE_L1W) && p1 && cbias[i])  codeBias[0]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_GLO) && ((i+1) == CODE_L1C) && cbias[i])        codeBias[0]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_GLO) && ((i+1) == CODE_L1W) && cbias[i])        codeBias[0]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_CMP) && ((i+1) == CODE_L2I) && cbias[i])        codeBias[0]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_GPS) && ((i+1) == CODE_L2X) && cbias[i])        codeBias[1]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_GLO) && ((i+1) == CODE_L2X) && cbias[i])        codeBias[1]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_GPS) && ((i+1) == CODE_L2W) && cbias[i])        codeBias[1]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_GLO) && ((i+1) == CODE_L2W) && cbias[i])        codeBias[1]=-CLIGHT/1E9*cbias[i];    
    if ((satsys(sat, NULL) == SYS_CMP) && ((i+1) == CODE_L6I) && cbias[i])        codeBias[2]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_GPS) && ((i+1) == CODE_L5I) && cbias[i])        codeBias[3]=-CLIGHT/1E9*cbias[i];
    if ((satsys(sat, NULL) == SYS_CMP) && ((i+1) == CODE_L7I) && cbias[i])        codeBias[4]=-CLIGHT/1E9*cbias[i];
    
    if ((satsys(sat, NULL) == SYS_GPS) && ((i+1) == CODE_L1C) && pbias[i]) phaseBias[0]=-pbias[i]*FREQ1/1E9;
    if ((satsys(sat, NULL) == SYS_GAL) && ((i+1) == CODE_L1B) && pbias[i]) phaseBias[0]=-pbias[i]*FREQ1/1E9;
    if ((satsys(sat, NULL) == SYS_CMP) && ((i+1) == CODE_L2I) && pbias[i]) phaseBias[0]=-pbias[i]*FREQ1_CMP/1E9;
    if ((satsys(sat, NULL) == SYS_GPS) && ((i+1) == CODE_L2W) && pbias[i]) phaseBias[1]=-pbias[i]*FREQ2/1E9;
    if ((satsys(sat, NULL) == SYS_CMP) && ((i+1) == CODE_L6I) && pbias[i]) phaseBias[2]=-pbias[i]*FREQ3_CMP/1E9;
    if ((satsys(sat, NULL) == SYS_GPS) && ((i+1) == CODE_L5I) && pbias[i]) phaseBias[3]=-pbias[i]*FREQ5/1E9;
    if ((satsys(sat, NULL) == SYS_GAL) && ((i+1) == CODE_L5Q) && pbias[i]) phaseBias[3]=-pbias[i]*FREQ5/1E9;
    if ((satsys(sat, NULL) == SYS_GAL) && ((i+1) == CODE_L7Q) && pbias[i]) phaseBias[4]=-pbias[i]*FREQ7/1E9;
    if ((satsys(sat, NULL) == SYS_CMP) && ((i+1) == CODE_L7I) && pbias[i]) phaseBias[4]=-pbias[i]*FREQ2_CMP/1E9;

    }
  }
  else
  {
  if (pnav->cbias[sat-1][0] && pnav->cbias[sat-1][1])
    {
    codeBias[0]=-pnav->cbias[sat-1][0]/(1.0-GPS_GAMMA);
    codeBias[1]=-GPS_GAMMA*pnav->cbias[sat-1][0]/(1.0-GPS_GAMMA);
    if (pnav->wlbias[sat-1])
      {
      double x=GPS_GAMMA*GPS_LAM1-GPS_LAM2;
      double alpha21=(GPS_LAM1-GPS_LAM2)/(GPS_LAM1+GPS_LAM2)/GPS_LAM1;
      double alpha22=(GPS_LAM1-GPS_LAM2)/(GPS_LAM1+GPS_LAM2)/GPS_LAM2;
      double y=pnav->wlbias[sat-1]-alpha21*codeBias[0]-alpha22*codeBias[1];
      phaseBias[0]=-GPS_LAM2*y/x;
      phaseBias[1]=-GPS_GAMMA*GPS_LAM1*y/x;
      }
    if (!p1)
      codeBias[0]-=pnav->cbias[sat-1][1];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
static int computeObs(obsd_t *data, double code[NB_FREQ], double phase[NB_FREQ], double doppler[NB_FREQ])
{
int i;
for (i=0;i<NB_FREQ;i++)
  {
  code[i]=0.0;
  phase[i]=0.0;
  doppler[i]=0.0;
  }
int p1=0;   // Cross-correlation
for (i=0;i<NFREQ+NEXOBS;i++)
  {
  if ((data->code[i] == CODE_L1C) || (data->code[i] == CODE_L1X) || (data->code[i] == CODE_L1I) || (data->code[i] == CODE_L1Q))
    {
    if (!p1 && data->P[i]) code[0]=data->P[i];
    if (!phase[0] && data->L[i]) phase[0]=data->L[i];
    if (!doppler[0] && data->D[i]) doppler[0]=data->D[i];
    }
  if ((data->code[i] == CODE_L1P) || (data->code[i] == CODE_L1W))
    {
    if (data->P[i])
      {
      p1=1;
      code[0]=data->P[i];
      }
    if (!phase[0] && data->L[i]) phase[0]=data->L[i];
    if (!doppler[0] && data->D[i]) doppler[0]=data->D[i];
    }
  if ((data->code[i] == CODE_L2P) || (data->code[i] == CODE_L2W))
    {
    code[1]=data->P[i];
    phase[1]=data->L[i];
    doppler[1]=data->D[i];
    }
  if ((data->code[i] == CODE_L6I) || (data->code[i] == CODE_L6Q))
    {
    code[2]=data->P[i];
    phase[2]=data->L[i];
    doppler[2]=data->D[i];
    }
  if ((data->code[i] == CODE_L5I) || (data->code[i] == CODE_L5Q) || (data->code[i] == CODE_L5X))
    {
    code[3]=data->P[i];
    phase[3]=data->L[i];
    doppler[3]=data->D[i];
    }
  if ((data->code[i] == CODE_L7I) || (data->code[i] == CODE_L7Q) || (data->code[i] == CODE_L7X))
//if ((data->code[i] == CODE_L8I) || (data->code[i] == CODE_L8Q) || (data->code[i] == CODE_L8X))
    {
    code[4]=data->P[i];
    phase[4]=data->L[i];
    doppler[4]=data->D[i];
    }
  }
return p1;
}

//////////////////////////////////////////////////////////////////////////////////////////////
static int typeBds(int sat, double rs[6])
{
#define GEO 0
#define IGSO 1
#define MEO 2
  int type = -999;
  static const double omegaEarth = 7292115.1467e-11;
  static const double gmWGS      = 398.6004418e12;
  double vix=(rs[3]-omegaEarth*rs[1]);
  double viy=(rs[4]+omegaEarth*rs[0]);
  double viz=(rs[5]);
  double nvi=sqrt(vix*vix+viy*viy+viz*viz);
  double nx=sqrt(rs[0]*rs[0]+rs[1]*rs[1]+rs[2]*rs[2]);
  double a = 1.0 / ( (2.0 / nx) - (nvi * nvi / gmWGS));
  double x = rs[1]*viz-rs[2]*viy;
  double y = rs[2]*vix-rs[0]*viz;
  double z = rs[0]*viy-rs[1]*vix;
  double inc=acos(z/sqrt(x*x+y*y+z*z));
  if (fabs(inc) < 0.1)
    type = GEO;
  else
  {
    if ( a > 35.0e6 )
      type = IGSO;
    else if ( a < 35.0e6 )
      type = MEO;
    else
      type = -1;
  } 
  return type;
}

//////////////////////////////////////////////////////////////////////////////////////////////
static void postprocess(const char *sobs, const char *ssp3, const char *sclk, const char *sbias, const char *satx, const char *sdcb, const char *ssbs, double xr[3], const char *schnl)
{
int j, imesiono;
static obs_t obs;
static nav_t nav;
static pcvs_t pcvs;
static bias_t bias;
static biases_t biases;
static sbs_t sbs;
static int tchannel[MAXPRNGLO] = { -999 };

memset(&obs, 0, sizeof(obs_t));
readrnx(sobs, 0, "", &obs, NULL, NULL);
memset(&nav, 0, sizeof(nav_t));
readsp3(ssp3, &nav, 0);

sbsreadmsg(ssbs, 0, &sbs);
imesiono=0;

readrnxc(sclk, &nav);
memset(&bias, 0, sizeof(bias_t));
readrnxbias(sbias, &bias);

generateBiases(&bias, &biases);
free(bias.pbias);

readpcv(satx, &pcvs);
readdcb(sdcb, &nav);

memset(tchannel, 0, MAXPRNGLO*sizeof(int));
readchannel(schnl, tchannel);

for (j=0;j<obs.n;j++)
  {
  if (obs.data[j].time.time)
    {
      double dtmax=1.0e9;

      while (1)
	{
	  if (imesiono >= sbs.n) break;
	  if (timediff(obs.data[j].time, gpst2time(sbs.msgs[imesiono].week, (double)sbs.msgs[imesiono].tow)) < 0.0) break;
	  sbsupdatecorr(&sbs.msgs[imesiono], &nav);
	  imesiono++;
	}
      double rs[6], dts[2], cbias[MAXCODE], pbias[MAXCODE];
      unsigned int day=(unsigned int)(obs.data[j].time.time/86400);
      unsigned int sec=(unsigned int)(obs.data[j].time.time-(time_t)day*86400);
      //double s = (double)sec+obs.data[j].time.sec;
      char id[4];
      pcv_t *pcv;
      satno2id(obs.data[j].sat, id);
      double code[NB_FREQ], phase[NB_FREQ], doppler[NB_FREQ];
      int p1=computeObs(&obs.data[j], code, phase, doppler);

      if (code[0]) {
	  pcv=searchpcv(obs.data[j].sat,"",obs.data[j].time,&pcvs);
	  
	  if (pcv){
	      gtime_t t;
	      nav.pcvs[obs.data[j].sat-1]=*pcv;
	      t=timeadd(obs.data[j].time,-code[0]/CLIGHT);	   
	      
	      if (peph2pos(t, obs.data[j].sat, &nav, 1, rs, dts, NULL)) {
		t=timeadd(t,-dts[0]);

		if (peph2pos(t, obs.data[j].sat, &nav, 1, rs, dts, NULL)) {
		  double codeBias[NB_FREQ], phaseBias[NB_FREQ];
		  int pcb=precbias(t, obs.data[j].sat, &biases, cbias, pbias);
		  computeBiases(&nav, pcb, cbias, pbias, obs.data[j].sat, p1, codeBias, phaseBias);
		  int channel=-999;
		  
		  if (*id == 'R')
		    channel=tchannel[atoi(id+1)-1];
		  
		  if (*id == 'C')
		    channel=typeBds(obs.data[j].sat, rs);
		  
		  int n1=0, nw=0;
		  if ((*id == 'G') || (*id == 'E') || (*id == 'C'))
		    nw=1;
		 
		  if (*id == 'G')
		    n1=1;
		  
		  double x = 1.0;
		  if ((*id == 'R') && (channel != -999))
		    x = FREQ1/(FREQ1_GLO+DFRQ1_GLO*channel);

		  if (*id == 'C')
		    x = FREQ1/FREQ1_CMP;
		  
		  double gamma = x*x;
		  double delay = 0.0;
		  if (xr[0] || xr[1] || xr[2]) {
		    double pos[3], e[3], azel[2], var;
		    ecef2pos(xr, pos);

		    if (rs[0] || rs[1] || rs[2]) {
		      geodist(rs, xr, e);
		      ecef2pos(xr, pos);
		      satazel(pos, e, azel);
		      int k = sbsioncorr(obs.data[j].time, &nav, pos, azel, &delay, &var);

		      if (!k) delay = 0.0;
		    }
		  }  
		  
		  double yaw = 0.;
		  int discontinuity = 0;
		  output(id, day+7305, sec, code, phase, doppler, channel, rs, dts[0]*CLIGHT, delay*gamma, 2, 0.0, codeBias, phaseBias, n1, nw, yaw, discontinuity);
		}
	      }
	  }
      }
    }
  }

freeobs(&obs);
freenav(&nav,0x08);
free(biases.pbiases);
free(pcvs.pcv);
}

//////////////////////////////////////////////////////////////////////////////////////////////
void ParseArgs(int argc, char *argv[], string &rinex, string &sp3, string &clk, string &bias, string &atx, string &dcb, string &sbs, double xr[3], string &chan)
{
  int iArg;
 
  if (argc<4 || argc>1000)
  {
    fprintf (stderr, "Usage : generateLowlevel\n"
             " - [E] -rinex : rinex_file\n"
             " - [E] -sp3 : sp3_file\n"
             " - [E] -clk : clk_file\n"
             " - [E] -bia : bia_file\n"
             " - [E] -atx : atx_file\n"
             " - [E] -dcb : dcb_file\n"
	     " - [E] -sbs : sbas_file xrover yrover zrover\n"
             " - [E] -chan : Glonass RFchannel_file\n"
             " - [S] > low_level_file\n");
    exit(1);
  }

  // Gestion des arguments.
  for (iArg=1; iArg<argc; iArg++)
  {
    if (strcmp(argv[iArg], "-rinex") == 0)
    {
      iArg++;
      rinex=(string)argv[iArg];
    }
    else if (strcmp(argv[iArg], "-chan") == 0)
    {
      iArg++;
      chan = (string)argv[iArg];
    }
    else if (strcmp(argv[iArg], "-sp3") == 0)
    {
      iArg++;
      sp3=(string)argv[iArg];
    }
    else if (strcmp(argv[iArg], "-clk") == 0)
    {
      iArg++;
      clk=(string)argv[iArg];
    }
    else if (strcmp(argv[iArg], "-bia") == 0)
    {
      iArg++;
      bias=(string)argv[iArg];
    }
    else if (strcmp(argv[iArg], "-atx") == 0)
    {
      iArg++;
      atx=(string)argv[iArg];
    }
    else if (strcmp(argv[iArg], "-dcb") == 0)
    {
      iArg++;
      dcb=(string)argv[iArg];
    }
    else if (strcmp(argv[iArg], "-sbs") == 0)
    {
      iArg++;
      sbs=(string)argv[iArg];
      iArg++;
      xr[0]=atof(argv[iArg]);
      iArg++;
      xr[1]=atof(argv[iArg]);
      iArg++;
      xr[2]=atof(argv[iArg]);
    }
    else
    {
      fprintf(stderr, "Unknown argument : '%s'.\n", argv[iArg]);
      exit(1);
    }
  } 
}

//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  string rinex, sp3, clk, bias, atx, dcb, sbs, chan;
  double xr[3];

  //ParseArgs(argc, argv, rinex, sp3, clk, bias, atx, dcb, chan);
  ParseArgs(argc, argv, rinex, sp3, clk, bias, atx, dcb, sbs, xr, chan);

  /* 
  if (argc < 8)
  {
    fprintf(stderr, "Usage: generateLowLevel rinex sp3 clk bia atx dcb chan\n");
    exit(1);
  }
  postprocess(rinex.c_str(), sp3.c_str(), clk.c_str(), bias.c_str(), atx.c_str(), dcb.c_str(), chan.c_str());
  */
  postprocess(rinex.c_str(), sp3.c_str(), clk.c_str(), bias.c_str(), atx.c_str(), dcb.c_str(), sbs.c_str(), xr, chan.c_str());
}
