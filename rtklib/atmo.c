/*------------------------------------------------------------------------------
* atmo.c : atmospheric data functions
*
*-----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <math.h>
#include "rtklib.h"

#define ROUND(x)    ((int)floor((x)+0.5))
#define THRESHOLD_MAX_DIST 200.0
#define MAX_STATION_NETWORK 20

static const char rcsid[]="$Id: atmo.c,v 1.1 2014-08-20 12:56:44 blota Exp $";

/* adjust weekly rollover of gps time ----------------------------------------*/
static void adjweek_cnesstec(rtcm_t *rtcm, double tow, gtime_t *time)
{
    double tow_p;
    int week;
    
    /* if no time, get cpu time */
    if (rtcm->time.time==0) *time=utc2gpst(timeget());
    tow_p=time2gpst(rtcm->time,&week);
    if      (tow<tow_p-302400.0) tow+=604800.0;
    else if (tow>tow_p+302400.0) tow-=604800.0;
    *time=gpst2time(week,tow);
}

extern int init_cnesstec(atmo_t *atmo)
{
    trace(3,"init_cnesstec:\n"); 
    atmo->nsta=0;
    atmo->cursta=0;
    atmo->staatmo=NULL;
    return 1;
}

extern void free_cnesstec(atmo_t *atmo)
{
    trace(3,"free_cnesstec:\n"); 
    free(atmo->staatmo);atmo->staatmo =NULL; atmo->nsta=0; atmo->cursta=0;
}
 
 
extern int update_stations_cnesstec(rtcm_t *rtcm, staatmo_t *staatmo)
{
  int i,j;
  for (i=0;i<rtcm->atmo.nsta;i++)
  {
    if ((sqrt(pow(rtcm->atmo.staatmo[i].x-staatmo->x,2)+pow(rtcm->atmo.staatmo[i].y-staatmo->y,2)+pow(rtcm->atmo.staatmo[i].z-staatmo->z,2))/1000.0) < 1.0)// dist in km
    { 
      //station update after decode
      rtcm->atmo.staatmo[i]=*staatmo;
      return i+1;
    }
  } 
  //new station to add if the for loop got no return
  staatmo_t *staatmo_realloc;    
  staatmo_realloc=(staatmo_t *)realloc(rtcm->atmo.staatmo, sizeof(staatmo_t)*(rtcm->atmo.nsta+1));
  if (staatmo_realloc == NULL)
  {
    trace(1,"update_stations: memory allocation error\n");
    return 0;
  } else {
    rtcm->atmo.staatmo=staatmo_realloc;
    rtcm->atmo.staatmo[rtcm->atmo.nsta]=*staatmo;      
    rtcm->atmo.nsta+=1;
    return rtcm->atmo.nsta;
  }
}

extern int atmocorr(atmo_t *atmo, gtime_t *time, double pos[6], staatmo_t *staatmo)
{ 
  int i, j, k, found, nsta_network=0, nsat_network, identical_sat[MAXPRNATMO]={0}, ambFixed_tot[MAXPRNATMO]={0};
  double sum_weight=0.0, tropo=0.0, covtropo=0.0, weight[MAX_STATION_NETWORK]={0.0},  dist[MAX_STATION_NETWORK]={0.0};
  staatmo_t staatmo_network[MAX_STATION_NETWORK];
  satatmo_t satatmo_network[MAXPRNATMO];
  for(k=0;k<MAXPRNATMO;k++)       
  {
    satatmo_network[k].sys=0;
    satatmo_network[k].prn=0;
    satatmo_network[k].ion=0.0;
    satatmo_network[k].covion=0.0;
    satatmo_network[k].ambFixed=0;
  }

  for (i=0;i<atmo->nsta;i++)
  {
dist[nsta_network]=sqrt(pow(atmo->staatmo[i].x-pos[0],2)+pow(atmo->staatmo[i].y-pos[1],2)+pow(atmo->staatmo[i].z-pos[2],2))/1000.0; // dist in km
    if(dist[nsta_network] < 1.0) //very close station
    {
      *staatmo=atmo->staatmo[i];
      trace(4,"atmo_corr: station is next to the rover (less than 1 kilometers)\n"); 
      return 1;
    }
    if (dist[nsta_network] < THRESHOLD_MAX_DIST) //station to take into account
    {
      if(nsta_network+1 >= MAX_STATION_NETWORK)
      {
        trace(3,"atmo_corr: too many station in the network\n"); 
	break;
      }
      staatmo_network[nsta_network]=atmo->staatmo[i];
      weight[nsta_network]=1/pow(dist[nsta_network],2); //partial weight
      sum_weight+=weight[nsta_network];  //sum of the partial weight
      nsta_network++; 
    }
  }
  if (nsta_network)
  {
    nsat_network=0;
    double penion=0.0;
    double penztd=0.0;
    for (i=0;i<nsta_network;i++)
    {
      penion=dist[i]/1000.0; //10km=1cm of covariance augmentation (iono)
      penztd=dist[i]/10000.0; //100km=1cm of covariance augmentation (ztd)
      weight[i]/=sum_weight; //total weight by station
      tropo+=staatmo_network[i].ztd*weight[i];
      covtropo+=(staatmo_network[i].covztd+penztd)*weight[i];
      for (j=0;j<staatmo_network[i].nsat;j++)
      {
        found=0;
        for (k=0;k<MAXPRNATMO;k++)       
        {
          if (staatmo_network[i].satatmo[j].sys==satatmo_network[k].sys && staatmo_network[i].satatmo[j].prn==satatmo_network[k].prn) 
	  {
	    satatmo_network[k].ion+=staatmo_network[i].satatmo[j].ion*weight[i];
            satatmo_network[k].covion+=(staatmo_network[i].satatmo[j].covion+penion)*weight[i]; 
	    identical_sat[k]++;	
	    ambFixed_tot[k]+=staatmo_network[i].satatmo[j].ambFixed;
	    found=1;
	  }       
        }
        if (!found) //new satellite
        {
          satatmo_network[nsat_network].sys=staatmo_network[i].satatmo[j].sys;
          satatmo_network[nsat_network].prn=staatmo_network[i].satatmo[j].prn;
       	  satatmo_network[nsat_network].ion=staatmo_network[i].satatmo[j].ion*weight[i];
          satatmo_network[nsat_network].covion=(staatmo_network[i].satatmo[j].covion+penion)*weight[i];
	  identical_sat[nsat_network]++;	 
	  ambFixed_tot[nsat_network]=staatmo_network[i].satatmo[j].ambFixed;
          nsat_network++;
        }           
      }
    }  
    for (k=0;k<MAXPRNATMO;k++) //ambFixed management      
      if (identical_sat[k])
        if((ambFixed_tot[k]/identical_sat[k])>0.8)
          satatmo_network[k].ambFixed=1;
          
    /* staatmo builder */
    staatmo->time.time=time->time;
    staatmo->time.sec=time->sec;    
    staatmo->x=pos[0];
    staatmo->y=pos[1];
    staatmo->z=pos[2];
    staatmo->ztd=tropo;
    staatmo->covztd=covtropo;
    staatmo->nsat=nsat_network;
    memcpy(staatmo->satatmo, satatmo_network, sizeof(satatmo_t)*MAXPRNATMO);
    /* staatmo builder */

    return 1;      
  }

  return 0;
}

extern int systemToBit(char system)
{
  int bit=0;
  switch(system){
    case 'G': bit=1; break;
    case 'R': bit=2; break;
    case 'E': bit=3; break;
    case 'C': bit=4; break;
    default: return bit;
  }
  return bit;
}

extern int bitToSys(int bit)
{
  int sys=0;
  switch(bit){
    case 1: sys=SYS_GPS; break;
    case 2: sys=SYS_GLO; break;
    case 3: sys=SYS_GAL; break;
    case 4: sys=SYS_CMP; break;
    default: return sys;
  }
  return sys;
}

/* Encode ATMO ascii */
extern int encode_cnesstec_ascii(rtcm_t *rtcm, char *enc)
{
    staatmo_t *staatmo;
    char field[100], id[4];
    int k;
   
    staatmo=&rtcm->atmo.staatmo[rtcm->atmo.cursta];
    *enc=0;
    unsigned int day=(unsigned int)(staatmo->time.time/86400);
    double sec=((double)staatmo->time.time-(double)day*86400.0)+staatmo->time.sec; 
    sprintf(field, "%d %.3f ", day+7305, sec); strcat(enc, field);
    sprintf(field, "%.3f %.3f %.3f ", staatmo->x, staatmo->y, staatmo->z); strcat(enc, field);
    sprintf(field, "%.3f %.3f ", staatmo->ztd, staatmo->covztd); strcat(enc, field);
    sprintf(field, "%d ", staatmo->nsat); strcat(enc, field);
    for (k=0;k<staatmo->nsat;k++) {
	  satno2id(satno(staatmo->satatmo[k].sys, staatmo->satatmo[k].prn), id);
      sprintf(field, "%s %.3f %.3f %d ", id,
			staatmo->satatmo[k].ion, staatmo->satatmo[k].covion, staatmo->satatmo[k].ambFixed);
      strcat(enc, field);
    }
    strcat(enc, "\n");

    return 1;
}

/* encode ATMO: atmospheric corrections -------------------------------------------*/
extern int encode_cnesstec(rtcm_t *rtcm, int sync)
{
    staatmo_t *staatmo;
    double tow;
    int i=24,j,k,epoch,week;
   
    setbitu(rtcm->buff,i,12,1272    ); i+=12; /* message number */
    staatmo=&rtcm->atmo.staatmo[rtcm->atmo.cursta];
    tow=time2gpst(staatmo->time,&week);
    epoch=ROUND(tow);
    setbitu(rtcm->buff,i,20,epoch); i+=20; /* gps epoch time */
    setbits(rtcm->buff,i,24,staatmo->x                 ); i+=24; /* x position*/
    setbits(rtcm->buff,i,24,staatmo->y                 ); i+=24; /* y position*/
    setbits(rtcm->buff,i,24,staatmo->z                 ); i+=24; /* z position*/
    setbitu(rtcm->buff,i,12,ROUND(staatmo->ztd/1E-3)   ); i+=12; /* zenith tropospheric delay */
    setbitu(rtcm->buff,i, 8,ROUND(staatmo->covztd/1E-3)); i+=8; /* zenith tropospheric delay covariance */
    setbitu(rtcm->buff,i, 8,staatmo->nsat              ); i+=8; /* number of sat */
    for (k=0;k<staatmo->nsat;k++) {	
	    setbitu(rtcm->buff,i, 2,staatmo->satatmo[k].sys               ); i+= 2; /* satellite system */	
        setbitu(rtcm->buff,i, 6,staatmo->satatmo[k].prn               ); i+= 6; /* satellite id */	
        setbits(rtcm->buff,i,16,ROUND(staatmo->satatmo[k].ion/1E-3)   ); i+=16; /* ionospheric delay */
        setbitu(rtcm->buff,i, 8,ROUND(staatmo->satatmo[k].covion/1E-3)); i+= 8; /* ionospheric delay covariance */
        setbitu(rtcm->buff,i, 2,staatmo->satatmo[k].ambFixed          ); i+= 2; /* ambiguity fixed */
    }
    rtcm->nbit=i;
    
    return 1;
}

/* decode ATMO: atmospheric corrections -------------------------------------------*/
extern int decode_cnesstec(rtcm_t *rtcm)
{
    double tow;
    int i=24,j,k,type,sat,system,sync;
    staatmo_t staatmo;
    
    
    type=getbitu(rtcm->buff,i,12); i+=12;
    tow=getbitu(rtcm->buff,i,20); i+=20;
    adjweek_cnesstec(rtcm,tow, &staatmo.time);

    staatmo.x=getbits(rtcm->buff,i,24); i+=24;
    staatmo.y=getbits(rtcm->buff,i,24); i+=24;
    staatmo.z=getbits(rtcm->buff,i,24); i+=24;
    staatmo.ztd=getbitu(rtcm->buff,i,12)*1E-3; i+=12;
    staatmo.covztd=getbitu(rtcm->buff,i,8)*1E-3; i+=8;
    staatmo.nsat=getbitu(rtcm->buff,i,8); i+=8;
    for (k=0;k<staatmo.nsat;k++) 
    {       	
      system=(int)getbitu(rtcm->buff,i,2); i+=2;
      staatmo.satatmo[k].sys=bitToSys(system);
      staatmo.satatmo[k].prn=(int)getbitu(rtcm->buff,i,6); i+=6;
      staatmo.satatmo[k].ion=getbits(rtcm->buff,i,16)*1E-3; i+=16;
      staatmo.satatmo[k].covion=getbitu(rtcm->buff,i,8)*1E-3; i+=8;
      staatmo.satatmo[k].ambFixed=getbitu(rtcm->buff,i,2); i+=2;
    }
    update_stations_cnesstec(rtcm, &staatmo);
    
    return 11;
}

/* decode ATMO: atmospheric corrections ascii -------------------------------------------*/
extern int decode_cnesstec_ascii(char *str, staatmo_t *staatmo)
{
char *token;
int i, day;
double sec;

token = strtok(str, " "); if (token == NULL) return 0;
day = atoi(token); token = strtok(NULL, " "); if (token == NULL) return 0;
sec = atof(token); token = strtok(NULL, " "); if (token == NULL) return 0;
staatmo->time.time=(double)(day-7305)*86400.0+floor(sec);
staatmo->time.sec=sec-floor(sec);
staatmo->x = atof(token); token = strtok(NULL, " "); if (token == NULL) return 0;
staatmo->y = atof(token); token = strtok(NULL, " "); if (token == NULL) return 0;
staatmo->z = atof(token); token = strtok(NULL, " "); if (token == NULL) return 0;
staatmo->ztd = atof(token); token = strtok(NULL, " "); if (token == NULL) return 0;
staatmo->covztd = atof(token); token = strtok(NULL, " "); if (token == NULL) return 0;
staatmo->nsat = atoi(token); token = strtok(NULL, " "); if (token == NULL) return 0;
for (i=0;i<staatmo->nsat;i++) 
  {
  staatmo->satatmo[i].sys = satsys(satid2no(token), &staatmo->satatmo[i].prn);
  token = strtok(NULL, " "); if (token == NULL) return 0;
  staatmo->satatmo[i].ion = atof(token); token = strtok(NULL, " "); if (token == NULL) return 0;
  staatmo->satatmo[i].covion = atof(token); token = strtok(NULL, " "); if (token == NULL) return 0;
  staatmo->satatmo[i].ambFixed = atoi(token); token = strtok(NULL, " "); if (token == NULL) return 0;
  }
return 1;
}

/* read atmo file -------------------------------------------*/
extern int read_cnesstec(FILE *fp, tab_staatmo_t *tab_staatmo)
{
int n;
char buff[10000];
staatmo_t staatmo;
    
tab_staatmo->n=0;
tab_staatmo->nr=0;
if (!tab_staatmo) return 0;
while (fgets(buff,sizeof(buff),fp))
  {
  if (decode_cnesstec_ascii(buff, &staatmo))
    {
    tab_staatmo->n++;
    if (!(tab_staatmo->staatmo=(staatmo_t *)realloc(tab_staatmo->staatmo,sizeof(staatmo_t)*tab_staatmo->n)))
      {
      free(tab_staatmo->staatmo);
      tab_staatmo->n=0;
      tab_staatmo->staatmo = NULL;
      return 0;
      }
    tab_staatmo->staatmo[tab_staatmo->n-1] = staatmo;
    }
  }
return tab_staatmo->n;
}

/* read atmo file -------------------------------------------*/
extern int update_stations_cnesstec_tab(rtcm_t *rtcm, tab_staatmo_t *tab_staatmo, gtime_t time)
{
staatmo_t staatmo;

while (1)
  {
  if (tab_staatmo->nr >= tab_staatmo->n) break;
  staatmo = tab_staatmo->staatmo[tab_staatmo->nr];
  if (timediff(staatmo.time, time) > 0.0) break;
  update_stations_cnesstec(rtcm, &staatmo);
  tab_staatmo->nr++;
  }

return 0;	
}

extern void free_tab_cnesstec(tab_staatmo_t *tab_staatmo)
{
if (tab_staatmo->staatmo != NULL) free(tab_staatmo->staatmo);
}
