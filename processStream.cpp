#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream> 
#include <unistd.h>

extern "C" {
#include "rtklib.h"
}
#include "rtrover.h"

static std::vector<rtcm_t> rtcm;
static std::vector<raw_t> raw;
static nav_t navIono;
static gtime_t vtec_time;
static vtec_t vtec;
static std::vector<double> msLast(MAXSAT,0.0), msCurrent(MAXSAT,0.0);
static double msCumul;
static std::vector<int> slot(MAXSAT,-999);
static int currentRover;
static int nbRover;
static bool verbose;

enum system { GPS, Glo, Gal, Bds, SystMax }; //same order as rtrover_utility.h
static std::vector<std::vector<rtrover_ionoCorr> > tabIonoCorr(SystMax); 

//////////////////////////////////////////////////////////////////////////////////////////////
void rtklib2bncTime(rtrover_time *t1, gtime_t *t2)
{
  int days=(unsigned int)(t2->time/86400);
  int sec=(unsigned int)(t2->time-(time_t)days*86400);
  t1->_mjd=days+40587;
  t1->_sec=(double)sec+t2->sec;
}

//////////////////////////////////////////////////////////////////////////////////////////////
int showmsg(char *format,...)
{
  va_list arg;
  va_start(arg,format); vfprintf(stderr,format,arg); va_end(arg);
  fprintf(stderr,*format?"\r":"\n");
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////
static char reverseSystem(const int syst)
{
  char system='G';
  switch (syst)
  {
    case GPS :
      system = 'G';
      break;
    case Glo :
      system = 'R';
      break;
    case Gal :
      system = 'E';
      break;
    case Bds :
      system = 'C';
      break;
   default :
      system = 'G';
  }
  return system;
}

//////////////////////////////////////////////////////////////////////
// Convert
static int convertSystem(const char system)
{
  int syst=GPS;
  switch (system)
  {
    case 'G' :
      syst = GPS;
      break;
    case 'R' :
      syst = Glo;
      break;
    case 'E' :
      syst = Gal;
      break;
    case 'C' :
      syst = Bds;
      break;
   default :
      syst = -1;
  }
  return syst;
}


//////////////////////////////////////////////////////////////////////////////////////////////
static int getMaxSatSyst (int syst)
{
  if (syst == GPS) return MAXPRNGPS;
  if (syst == Glo) return MAXPRNGLO;
  if (syst == Gal) return MAXPRNGAL;
  if (syst == Bds) return MAXPRNCMP;
}

//////////////////////////////////////////////////////////////////////////////////////////////
static void clearIonoCorr()
{
  for(int isyst=0; isyst<SystMax;isyst++)
  {
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)//reinit iono
    {
      tabIonoCorr[isyst][isat]._staName = NULL;
      tabIonoCorr[isyst][isat]._flag = 2;
      tabIonoCorr[isyst][isat]._satellite._system = reverseSystem(isyst);
      tabIonoCorr[isyst][isat]._satellite._number = isat+1;
      tabIonoCorr[isyst][isat]._value = 0.0;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
void convertDate(const int jj, const double ss, int *day, int *month, int *year, int *hour, int *minute, double *second)
{
  long i, j, k, y, m, d;
  double min;

  i = (long)jj + 712164;
  j = ((4*i-1) % 146097)/4;
  k = (((4*j+3) % 1461)+4)/4;
  y = 100 * ((4*i-1)/146097) + (4*j+3)/1461;
  m = (5*k-3)/153;
  d = (((5*k-3) % 153)+5)/5;

  if (m < 10)
  {
    *day = (int)d;
    *month = (int)(m+3);
    *year = (int)(y);
  }
  else
  {
    *day = (int)d;
    *month = (int)(m-9);
    *year = (int)(y+1);
  }

  *hour = (int)((long)ss / 3600);
  min = ss - 3600.0*(double)(*hour);
  *minute = (int)(min / 60);
  *second = min - 60.0*(double)(*minute);
}

//////////////////////////////////////////////////////////////////////////////////////////////
/*void settspan(gtime_t ts, gtime_t te)
{
}*/

//////////////////////////////////////////////////////////////////////////////////////////////
/*void settime(gtime_t time)
{
}*/

//////////////////////////////////////////////////////////////////////////////////////////////
//
// The configuration is the same for each rover
void readParam(char *confFileName, bool lowlevel, rtrover_opt *opt, rtrover_additionalOpt *addOpt)
{
  FILE *fichier;
  char mode[20], dum[200];
  static char _antexFileName[200];    // should stay static
  int use = -1;
  int nbSatFixAmb = 0;

  fichier=fopen(confFileName,"r");
  if (fichier != NULL)
  {
    fscanf(fichier, "%s%s\n", mode, dum);
    if (!strcmp(mode, "mode_PPP_AR"))
      opt->_mode=mode_PPP_AR;
    if (!strcmp(mode, "mode_PPP_DF"))
      opt->_mode=mode_PPP_DF;
    if (!strcmp(mode, "mode_PPP_SF"))
      opt->_mode=mode_PPP_SF;
    if (!strcmp(mode, "mode_SPP_DF"))
      opt->_mode=mode_SPP_DF;
    if (!strcmp(mode, "mode_SPP_SF"))
      opt->_mode=mode_SPP_SF;
    opt->_roverName = "";
    opt->_xyzAprRover[0] = 0.0;
    opt->_xyzAprRover[1] = 0.0;
    opt->_xyzAprRover[2] = 0.0;
    fscanf(fichier, "%s%s\n", _antexFileName, dum);
    opt->_antexFileName = _antexFileName;
    fscanf(fichier, "%d%d%d%s\n", &addOpt->_fixN1, &addOpt->_fixNw, &addOpt->_fixEx, dum);
    fscanf(fichier, "%d%s\n", &use, dum);//Gps
    addOpt->_useGPS = use ? true : false;
    fscanf(fichier, "%d%s\n", &use, dum);//Glonass
    opt->_useGlonass = use ? true : false;
    fscanf(fichier, "%d%s\n", &use, dum);//Galileo
    addOpt->_useGalileo = use ? true : false;
    fscanf(fichier, "%d%s\n", &use, dum);//Beidou
    addOpt->_useBeidou = use ? true : false;  
    fscanf(fichier, "%d%s\n", &addOpt->_correction, dum);
    fscanf(fichier, "%d%s\n", &addOpt->_reset, dum);
    fscanf(fichier, "%d%s\n", &addOpt->_outputVerbose, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_dt, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_maxAge, dum);
    fscanf(fichier, "%d%s\n", &addOpt->_minFix, dum);
    fscanf(fichier, "%d%s\n", &addOpt->_maxElim, dum);
    fscanf(fichier, "%d%s\n", &addOpt->_raim, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_thrMap, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigIniTro, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigModTro, dum);
    fscanf(fichier, "%d%s\n", &nbSatFixAmb, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_thrAmb, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigIniClkBias, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigModClkBias, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigIniIono, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigModIono, dum);
    fscanf(fichier, "%lf%lf%lf%s\n", &addOpt->_sigMeasIono[0], &addOpt->_sigMeasIono[1], &addOpt->_sigMeasIono[2], dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_thrMeasIono, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigMeasTropo, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_thrMeasTropo, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigIniPos, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigModPos, dum);
    fscanf(fichier, "%d%s\n", &addOpt->_dtMax, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_thrMeasCode, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_thrMeasPhase, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigMeasCodeGps, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigMeasPhaseGps, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigMeasCodeGlo, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigMeasPhaseGlo, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigMeasCodeGal, dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_sigMeasPhaseGal, dum);
    fscanf(fichier, "%lf%lf%lf%s\n", &addOpt->_sigMeasCodeBds[0], &addOpt->_sigMeasCodeBds[1], &addOpt->_sigMeasCodeBds[2], dum);
    fscanf(fichier, "%lf%lf%lf%s\n", &addOpt->_sigMeasPhaseBds[0], &addOpt->_sigMeasPhaseBds[1], &addOpt->_sigMeasPhaseBds[2], dum);
    fscanf(fichier, "%lf%s\n", &addOpt->_smooth, dum);
    fclose(fichier);
    addOpt->_lowlevel =lowlevel ? 1 : 0;
    addOpt->_nbSatFixAmb = (nbSatFixAmb >= 0) ? nbSatFixAmb : 0;
  }
  else
  {
    fprintf(stderr, "Error : can't open configuration file !\n");
    exit(1);
  }

  opt->_antNameRover="";
  opt->_logLevel=0;
  opt->_minobs=0;
}

//////////////////////////////////////////////////////////////////////////////////////////////
void readRover(char *roverFileName)
{
  std::ifstream file(roverFileName, std::ios::in);
  nbRover = 0;
 
  if(file)
  {
    std::string ligne;
    while(getline(file, ligne))
    {
      std::istringstream iss(ligne);
      std::string roverName;
      double x, y, z;
      iss >> roverName >> x >> y >> z;
      double aprPos[3]={x,y,z};
      rtrover_addRover(roverName, aprPos);
      nbRover++;
    }
    file.close();
  }
  else
  {
    fprintf(stderr, "Error : can't open rover file !\n");
    exit(1);
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////
void setDCB(char *dcbFileName)  // * wildcard
{
  rtrover_satDCBBiases satDCBBiases;
  nav_t nav_biases;

  readdcb(dcbFileName, &nav_biases);
  for (int i=0;i<MAXSAT;i++)
  {
    char id[4];
    satno2id(i+1, id);
    if ((id[0] == 'G') || (id[0] == 'R'))
    {
      satDCBBiases._satellite._system = id[0];
      satDCBBiases._satellite._number = atoi(id+1);
      satDCBBiases._time._mjd = 0;
      satDCBBiases._time._sec = 0.0;
      satDCBBiases._P1P2 = nav_biases.cbias[i][0];
      satDCBBiases._P1C1 = nav_biases.cbias[i][1];
      satDCBBiases._P2C2 = nav_biases.cbias[i][2];
      rtrover_putDCBBiases(1, &satDCBBiases);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
void decode_begin(char *confFileName, char *rover, char *dcbFileName, bool lowlevel)
{
  rtrover_opt opt;
  rtrover_additionalOpt addOpt;
  readParam(confFileName, lowlevel, &opt, &addOpt);
  readRover(rover);
  rtrover_initOptions(&opt);
  rtrover_setOptions(&opt);

  rtrover_setAdditionalOpt(&addOpt);

  setDCB(dcbFileName);

  navIono.seph=NULL;
  for(int isyst=0; isyst < SystMax; isyst ++)
    tabIonoCorr[isyst]=std::vector<rtrover_ionoCorr> (getMaxSatSyst(isyst));
  clearIonoCorr();

  gtime_t t0={0};
  vtec_time=t0;

  msCumul = 0.0;
  currentRover = 0;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////
void decode_end()
{
  for (int i=0;i<(int)rtcm.size();i++)
  {
    free_rtcm(&rtcm[i]);
  }
  for (int i=0;i<(int)raw.size();i++)
  {
    free_raw(&raw[i]);
  }
  rtrover_destroy();
}

//////////////////////////////////////////////////////////////////////////////////////////////
double vtecioncorr(gtime_t time, const double *pos, const double *possat)
{
#define MIN(x,y)    ((x)<(y)?(x):(y))
double E, A, lams, psipp, phipp;
int l, m, n;

if (timediff(time, vtec_time) > 300.0)
  return 0.0;

/* SSR2-VTEC-09 */
double r=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
double ux=pos[0]/r;
double uy=pos[1]/r;
double uz=pos[2]/r;

double lamr=atan2(uy,ux);                // longitude rover
double phir=atan(uz/sqrt(ux*ux+uy*uy));  // latitude rover

double nordx=-sin(phir)*cos(lamr);
double nordy=-sin(phir)*sin(lamr);
double nordz=cos(phir);
double estx=-sin(lamr);
double esty=cos(lamr);
double estz=0.0;

double vx=possat[0]-pos[0];
double vy=possat[1]-pos[1];
double vz=possat[2]-pos[2];
double nv=sqrt(vx*vx+vy*vy+vz*vz);
double cosp=(ux*vx+uy*vy+uz*vz)/nv;

E=M_PI/2.0-acos(cosp);
double x=ux*vx+uy*vy+uz*vz;
double projvx=vx-x*ux;
double projvy=vy-x*uy;
double projvz=vz-x*uz;
double projvn=sqrt(projvx*projvx+projvy*projvy+projvz*projvz);
double projn=(nordx*projvx+nordy*projvy+nordz*projvz)/projvn;
double proje=(estx*projvx+esty*projvy+estz*projvz)/projvn;
if (proje >= 0.0)
    A=acos(projn);
else
    A=2.0*M_PI-acos(projn);

psipp=M_PI/2.0-E-asin(r/(6370000.0+vtec.height)*cos(E));
phipp=asin(sin(phir)*cos(psipp)+cos(phir)*sin(psipp)*cos(A));

double lampp;
if (((phir>=0.0) && (tan(psipp)*cos(A)>tan(M_PI/2.0-phir))) || ((phir<0.0) && (-tan(psipp)*cos(A)>tan(M_PI/2.0+phir))))
  lampp=lamr+M_PI-asin(sin(psipp)*sin(A)/cos(phipp));
else
  lampp=lamr+asin(sin(psipp)*sin(A)/cos(phipp));

double t=(double)(time.time%86400);
lams=fmod(lampp+(t-50400.0)/43200.0*M_PI, 2.0*M_PI);
phipp=fmod(phipp+M_PI/2.0, 2.0*M_PI)-M_PI/2.0;
lams=fmod(lams+M_PI, 2.0*M_PI)-M_PI;

double xc = cos(phipp)*cos(lams);
double yc = cos(phipp)*sin(lams);
double zc = sin(phipp);

double tabcos[MAXVTEC+1]={0};
double tabsin[MAXVTEC+1]={0};
tabcos[0] = 1.0;
tabsin[0] = 0.0;
for (n=1;n<=vtec.nOrd;n++)
  {
  tabcos[n] = xc*tabcos[n-1] - yc*tabsin[n-1];
  tabsin[n] = yc*tabcos[n-1] + xc*tabsin[n-1];
  }

double C[MAXVTEC+1][MAXVTEC+1]={0};
double S[MAXVTEC+1][MAXVTEC+1]={0};
l=0;
for (m=0;m<=vtec.nOrd;m++)
  for (n=m;n<=vtec.nDeg;n++)
    C[n][m]=vtec.cosCoeffs[l++];
l=0;
for (m=1;m<=vtec.nOrd;m++)
  for (n=m;n<=vtec.nDeg;n++)
    S[n][m]=vtec.sinCoeffs[l++];
    
// Normalisation coefficients
double no[MAXVTEC+1][MAXVTEC+1]={0};
for (l=0;l<=vtec.nDeg;l++)
  {
  no[l][0] = sqrt(1.0/(2.0*l+1.0));
  for (m=1;m<=l;m++)
    {
    double p = 1.0;
    for (n=l-m+1.0;n<=l+m;n++)
        p*=n;
    no[l][m] = sqrt(p/(2.0*l+1.0)/2.0);
    }
  }

// Legendre function
double w[MAXVTEC+1][MAXVTEC+1]={0};
w[0][0] = 1;
w[1][1] = 1;
w[1][0] = zc;
for (n=2;n<=vtec.nDeg;n++)
  {
  w[n][n] = (double)(2*n-1)*w[n-1][n-1];
  w[n][n-1] = (double)(2*n-1)*zc*w[n-1][n-1];
  for (m=n-2;m>=0;m--)
    w[n][m] = (double)(2*n-1)/(double)(n-m)*zc*w[n-1][m]-(double)(n+m-1)/(double)(n-m)*w[n-2][m];
  }

double res=0.0;
for (n=0;n<=vtec.nDeg;n++)
  for (m=0;m<=MIN(n, vtec.nOrd);m++)
    res+=(C[n][m]*tabcos[m]+S[n][m]*tabsin[m])*w[n][m]/no[n][m];

double vtec=40.3e16*res/FREQ1/FREQ1;
if (vtec < 0.0) return 0.0;

return vtec/sin(E+psipp);
}

//////////////////////////////////////////////////////////////////////////////////////////////
void sbasIonoCorr(rtrover_satellite sat, obsd_t obs, double gamma)
{
  int jj, k;
  double delay, var, ss;
  double r[6], tropo[3], pos[3], possat[3], amb[6], iono[2], e[3], azel[2];
  rtrover_getRoverInformation(&jj, &ss, r, tropo);

  // SBAS Iono
  if (r[0] || r[1] || r[2])
  {
    ecef2pos(r, pos);
    rtrover_getSatelliteInformation(sat, possat, amb, iono);
    if (possat[0] || possat[1] || possat[2])
    {
      geodist(possat, r, e);
      ecef2pos(r, pos);
      satazel(pos, e, azel);
      k = sbsioncorr(obs.time, &navIono, pos, azel, &delay, &var);
      if (!k) delay = 0.0;
      //if (delay != 0.0) printf ("iono sbas\n");
      int isat=sat._number - 1;
      char syst=sat._system;
      //priority to iono from verbose mode
      rtrover_ionoCorr &ionoCorr = tabIonoCorr[convertSystem(syst)][isat];
      if ((delay != 0.0) && ((ionoCorr._value == 0.0) || (ionoCorr._flag == 3)))
      {
        ionoCorr._staName = NULL;
        ionoCorr._satellite = sat;
        ionoCorr._value = gamma * delay;
        ionoCorr._flag = 2;
      }
      //printf("sbas ion %d %lf\n", k, delay);
    }
  }

  // vtec iono
  if (r[0] || r[1] || r[2])
  {
    rtrover_getSatelliteInformation(sat, possat, amb, iono);
    if (possat[0] || possat[1] || possat[2])
    {
      delay = vtecioncorr(obs.time, r, possat);
      int isat=sat._number - 1;
      char syst=sat._system;
      //priority to iono from verbose mode
      rtrover_ionoCorr &ionoCorr = tabIonoCorr[convertSystem(syst)][isat];
      if ((delay != 0.0) && (ionoCorr._value == 0.0))
      {
        ionoCorr._staName = NULL;
        ionoCorr._satellite = sat;
        ionoCorr._value = gamma * delay;
        ionoCorr._flag = 3;
      }
      //printf("vtec ion %d %lf\n", k, delay);
    }
  }

  // sbas slow and fast
  double rs[3], dts=0;
  rs[0] = 0.0; rs[1] = 0.0; rs[2] = 0.0;
  k = sbssatcorr(obs.time, obs.sat, &navIono, rs, &dts, &var);
  if (k)
  {
    /*    
	  printf("sbas slow and fast %d %d %d %lf %lf %lf %lf\n",
	         obs.sat, navIono.sbssat.sat[obs.sat-1].lcorr.iode, navIono.sbssat.sat[obs.sat-1].fcorr.iodf,
		 rs[0], rs[1], rs[2], CLIGHT*dts);
    */
    rtrover_sbasCorr sbasCorr;
    sbasCorr._satellite = sat;
    rtklib2bncTime(&sbasCorr._time, &obs.time);
    sbasCorr._IODE = navIono.sbssat.sat[obs.sat-1].lcorr.iode;
    sbasCorr._orb[0] = rs[0];
    sbasCorr._orb[1] = rs[1];
    sbasCorr._orb[2] = rs[2];
    sbasCorr._clk = dts;
    rtrover_putSbasCorrs(1, &sbasCorr);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
void decode_stream_obs(obs_t *obs, int ubx)
{
  int i, j;
  int numSatRover=0;
  rtrover_satObs *satObsRover=NULL;

  msCurrent.assign(MAXSAT,0.0);

  if (obs->n)
  {
    satObsRover = new rtrover_satObs[obs->n];
    for (j=0;j<obs->n;j++)
    {
      double gamma=1.0; //PATCH
      char id[4];
      satno2id(obs->data[j].sat, id);
      if ((id[0] != 'G') && (id[0] != 'R') && (id[0] != 'E') && (id[0] != 'C'))
	continue;

      if (!obs->data[j].time.time)
        continue;
      satObsRover[numSatRover]._satellite._system = id[0];
      satObsRover[numSatRover]._satellite._number = atoi(id+1);
      rtklib2bncTime(&satObsRover[numSatRover]._time, &obs->data[j].time);
      satObsRover[numSatRover]._slotNumber = slot[obs->data[j].sat-1];
      satObsRover[numSatRover]._numObs =3;
      satObsRover[numSatRover]._obs = new rtrover_obs[satObsRover[numSatRover]._numObs];
      // Init
      for (i=0;i<satObsRover[numSatRover]._numObs;i++)
      {
        satObsRover[numSatRover]._obs[i]._rnxType2ch[0] = 'X';
        satObsRover[numSatRover]._obs[i]._rnxType2ch[1] = 'X';
        satObsRover[numSatRover]._obs[i]._code = 0.0;
        satObsRover[numSatRover]._obs[i]._codeValid = false;
        satObsRover[numSatRover]._obs[i]._phase = 0.0;
        satObsRover[numSatRover]._obs[i]._phaseValid = false;
        satObsRover[numSatRover]._obs[i]._doppler = 0.0;
        satObsRover[numSatRover]._obs[i]._dopplerValid = false;
      }
      
      for (i=0;i<NFREQ+NEXOBS;i++)      
      {
        double freq1 = FREQ1;
        double freq2 = FREQ2;
        double freq5 = FREQ5;
	double freq6 = FREQ3_CMP;//depends on the rtklib version
        double freq7 = FREQ7;
	//double freq1Bds = FREQ1_CMP;
	//double freq6Bds = FREQ3_CMP; //depends on the rtklib version
        if ((satObsRover[numSatRover]._slotNumber != -999) && (id[0] == 'R'))
        {
          gamma = (freq1 / (FREQ1_GLO+DFRQ1_GLO*satObsRover[numSatRover]._slotNumber)) * (freq1 / (FREQ1_GLO+DFRQ1_GLO*satObsRover[numSatRover]._slotNumber));
	  freq1 = FREQ1_GLO+DFRQ1_GLO*satObsRover[numSatRover]._slotNumber;
          freq2 = FREQ2_GLO+DFRQ2_GLO*satObsRover[numSatRover]._slotNumber;
        }
	if (id[0] == 'C')
	{
	  gamma = freq1 / FREQ1_CMP * freq1 / FREQ1_CMP;
	  freq1 = FREQ1_CMP;//depends on the rtklib version
	}

	if ((obs->data[j].code[i] == CODE_L1C) ||
	    (obs->data[j].code[i] == CODE_L1P) ||
	    (obs->data[j].code[i] == CODE_L1W) ||
	    (obs->data[j].code[i] == CODE_L1X) || 
	    (obs->data[j].code[i] == CODE_L1Q) || 
	    (obs->data[j].code[i] == CODE_L1I))
        {
	  satObsRover[numSatRover]._obs[0]._rnxType2ch[0] = '1';
	  switch (obs->data[j].code[i])
	  {
	    case CODE_L1C :
	      satObsRover[numSatRover]._obs[0]._rnxType2ch[1] = 'C';
	      break;
	    case CODE_L1P :
	      satObsRover[numSatRover]._obs[0]._rnxType2ch[1] = 'P';
	      break;
	    case CODE_L1W :
	      satObsRover[numSatRover]._obs[0]._rnxType2ch[1] = 'W';
	      break;
	    case CODE_L1X :
	    case CODE_L1Q :
	    case CODE_L1I :
	      satObsRover[numSatRover]._obs[0]._rnxType2ch[1] = 'X';
	      break;      
	  }
          satObsRover[numSatRover]._obs[0]._code = obs->data[j].P[i];
          satObsRover[numSatRover]._obs[0]._codeValid = true;
          satObsRover[numSatRover]._obs[0]._phase = obs->data[j].L[i];
          satObsRover[numSatRover]._obs[0]._phaseValid = true;
          satObsRover[numSatRover]._obs[0]._doppler = -obs->data[j].D[i]*CLIGHT/freq1;
          satObsRover[numSatRover]._obs[0]._dopplerValid = true;

          // ms
          msCurrent[obs->data[j].sat-1] = obs->data[j].P[i];
	  
	  sbasIonoCorr(satObsRover[numSatRover]._satellite, obs->data[j],gamma);
        }

        if ((obs->data[j].code[i] == CODE_L2X) ||
	    (obs->data[j].code[i] == CODE_L2Q) || (obs->data[j].code[i] == CODE_L2I) ||
	    (obs->data[j].code[i] == CODE_L2W) ||
	    (obs->data[j].code[i] == CODE_L2P))
        {
 	  if (id[0] == 'C') //some stations have the wrong BeiDou obs code
	  {
	    satObsRover[numSatRover]._obs[0]._rnxType2ch[0] = '1';
            satObsRover[numSatRover]._obs[0]._rnxType2ch[1] = 'X';
            satObsRover[numSatRover]._obs[0]._code = obs->data[j].P[i];
            satObsRover[numSatRover]._obs[0]._codeValid = true;
            satObsRover[numSatRover]._obs[0]._phase = obs->data[j].L[i];
            satObsRover[numSatRover]._obs[0]._phaseValid = true;
            satObsRover[numSatRover]._obs[0]._dopplerValid = true;
	    satObsRover[numSatRover]._obs[0]._doppler = -obs->data[j].D[i]*CLIGHT/freq1;
	  }
	  else
	  {
            satObsRover[numSatRover]._obs[1]._rnxType2ch[0] = '2';
	    switch (obs->data[j].code[i])
            {
              case CODE_L2X :
                satObsRover[numSatRover]._obs[1]._rnxType2ch[1] = 'X';
                break;
              case CODE_L2W :
                satObsRover[numSatRover]._obs[1]._rnxType2ch[1] = 'W';
                break;
              case CODE_L2P :
                satObsRover[numSatRover]._obs[1]._rnxType2ch[1] = 'P';
                break;      
	    }
            satObsRover[numSatRover]._obs[1]._code = obs->data[j].P[i];
            satObsRover[numSatRover]._obs[1]._codeValid = true;
            satObsRover[numSatRover]._obs[1]._phase = obs->data[j].L[i];
            satObsRover[numSatRover]._obs[1]._phaseValid = true;
            satObsRover[numSatRover]._obs[1]._dopplerValid = true;
	    satObsRover[numSatRover]._obs[1]._doppler = -obs->data[j].D[i]*CLIGHT/freq2;
	  }
        }
        if ((obs->data[j].code[i] == CODE_L5I) || (obs->data[j].code[i] == CODE_L5Q) || (obs->data[j].code[i] == CODE_L5X))
        {
          int ind = 2; // cas GPS
	  if (id[0] == 'E')
	    ind = 1;
	  satObsRover[numSatRover]._obs[ind]._rnxType2ch[0] = '5';
          satObsRover[numSatRover]._obs[ind]._rnxType2ch[1] = 'Q';
          satObsRover[numSatRover]._obs[ind]._code = obs->data[j].P[i];
          satObsRover[numSatRover]._obs[ind]._codeValid = true;
          satObsRover[numSatRover]._obs[ind]._phase = obs->data[j].L[i];
          satObsRover[numSatRover]._obs[ind]._phaseValid = true;
          satObsRover[numSatRover]._obs[ind]._dopplerValid = true;
          satObsRover[numSatRover]._obs[ind]._doppler = -obs->data[j].D[i]*CLIGHT/freq5;
        }
	if ((obs->data[j].code[i] == CODE_L6I) || (obs->data[j].code[i] == CODE_L6Q) || (obs->data[j].code[i] == CODE_L6X))
        {
	  satObsRover[numSatRover]._obs[1]._rnxType2ch[0] = '6';
          satObsRover[numSatRover]._obs[1]._rnxType2ch[1] = 'Q';
          satObsRover[numSatRover]._obs[1]._code = obs->data[j].P[i];
          satObsRover[numSatRover]._obs[1]._codeValid = true;
          satObsRover[numSatRover]._obs[1]._phase = obs->data[j].L[i];
          satObsRover[numSatRover]._obs[1]._phaseValid = true;
          satObsRover[numSatRover]._obs[1]._dopplerValid = true;
          satObsRover[numSatRover]._obs[1]._doppler = -obs->data[j].D[i]*CLIGHT/freq6;
        }
	if ((obs->data[j].code[i] == CODE_L7I) || (obs->data[j].code[i] == CODE_L7Q) || (obs->data[j].code[i] == CODE_L7X))
        {
	  satObsRover[numSatRover]._obs[2]._rnxType2ch[0] = '7';
          satObsRover[numSatRover]._obs[2]._rnxType2ch[1] = 'Q';
          satObsRover[numSatRover]._obs[2]._code = obs->data[j].P[i];
          satObsRover[numSatRover]._obs[2]._codeValid = true;
          satObsRover[numSatRover]._obs[2]._phase = obs->data[j].L[i];
          satObsRover[numSatRover]._obs[2]._phaseValid = true;
          satObsRover[numSatRover]._obs[2]._dopplerValid = true;
          satObsRover[numSatRover]._obs[2]._doppler = -obs->data[j].D[i]*CLIGHT/freq7;
        }	
      }
      numSatRover++;
    }
  }
  
  // sauts ms
  int n = 0;
  double x = 0.0;
  double deltaCumul = 0.0;
  for (i=0;i<MAXSAT;i++)
    if (msLast[i] && msCurrent[i])
    {
      x += msCurrent[i] - msLast[i];
      n++;
    }
  if (n)
  {
    x /= (double)n;
    deltaCumul = 1.0e-3*floor(x/CLIGHT/1.0e-3+0.5);
    msCumul += deltaCumul;
  }
  // fprintf(stderr, "cumul %lf %lf %lf\n", x, x/CLIGHT/1.0e-3+0.5, msCumul);
  if (ubx)
    for (i=0;i<numSatRover;i++)
    {
      if (satObsRover[i]._obs[0]._rnxType2ch[0] == '1')
      if ((satObsRover[i]._obs[0]._rnxType2ch[1] == 'C') /*|| (satObsRover[i]._obs[0]._rnxType2ch[1] == 'X')*/)
      {
        double freq1 = FREQ1;
        if ((satObsRover[i]._slotNumber != -999) && (satObsRover[i]._satellite._system == 'R'))
          freq1 = FREQ1_GLO+DFRQ1_GLO*satObsRover[i]._slotNumber;
        if (satObsRover[i]._obs[0]._phase)
          satObsRover[i]._obs[0]._phase += freq1*msCumul;
        if (satObsRover[i]._obs[0]._doppler)
          if (deltaCumul)
        satObsRover[i]._obs[0]._doppler = 0.0; 
      }
    }
  for (i=0;i<MAXSAT;i++)
    msLast[i] = msCurrent[i];

  rtrover_output output;
  
  // Put ionoCorr
  std::vector <rtrover_ionoCorr> newTabIonoCorr;
  for(int isyst=0;isyst<SystMax;isyst++)
  //for(int isyst=0;isyst<=GPS;isyst++) //only GPS, for now 
    for(int isat=0;isat<getMaxSatSyst(isyst);isat++) 
      newTabIonoCorr.push_back(tabIonoCorr[isyst][isat]);
  rtrover_putIonoCorrs(newTabIonoCorr.size(), &newTabIonoCorr[0]);
  // clear tabIonoCorr
  clearIonoCorr();

  rtrover_processEpoch(numSatRover, satObsRover, 0, NULL, &output);
  fputs(output._log, stdout);
  fflush(stdout);
  rtrover_freeOutput(&output);
  if (numSatRover)
  {
    for (j=0;j<numSatRover;j++)
      delete [] satObsRover[j]._obs;
    delete [] satObsRover;
  }

  if (verbose)
  {
    int jj;
    double ss, r[6], tropo[3];
    rtrover_getRoverInformation(&jj, &ss, r, tropo);
    if (jj)
    {
      int day, month, year, hour, minute;
      double second;
      convertDate(jj, ss,&day, &month, &year, &hour, &minute, &second);
      fprintf(stderr, "%d %02d-%02d-%02d %02d:%02d:%06.3lf %.3lf ", currentRover+1, year%100, month, day, hour, minute, second, tropo[0] + tropo[1]);
      for (int isyst=0; isyst<SystMax;isyst++)
      {
        for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
        {
          double possat[3], amb[6], iono[2];
          rtrover_satellite sat;
	  sat._system = reverseSystem(isyst);
          sat._number = isat+1;
          rtrover_getSatelliteInformation(sat, possat, amb, iono);
          if (amb[2] && iono[0])
          {           
	    if (amb[5] && (amb[5] < 0.01))             
	      fprintf(stderr, "%c%02d %.3lf ", sat._system, sat._number, iono[0]);
          }
        }
      }
      fprintf(stderr, "\n");
      fflush(stderr);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
void decode_stream_eph(nav_t *nav)
{
  
  for (int i=0;i<nav->n;i++)
  {
    char id[4];
    int isat=0;
    *id='\0';
    satno2id(nav->eph[i].sat, id);
    isat=atoi(id+1);
    
    if (satsys(nav->eph[i].sat,NULL)==SYS_GPS)
    {
      if ((isat >= MINPRNGPS) && (isat <= MAXPRNGPS))
      {
	rtrover_ephGPS eph;
	eph._satellite._system = 'G';
	eph._satellite._number = isat;	
	rtklib2bncTime(&eph._TOC, &nav->eph[i].toc);
	rtklib2bncTime(&eph._TOE, &nav->eph[i].toe);
	eph._IODE = nav->eph[i].iode;
	eph._IODC = nav->eph[i].iodc;
	eph._clock_bias = nav->eph[i].f0;  
	eph._clock_drift = nav->eph[i].f1;
	eph._clock_driftrate = nav->eph[i].f2;
	eph._Crs = nav->eph[i].crs;
	eph._Delta_n = nav->eph[i].deln;
	eph._M0 = nav->eph[i].M0;
	eph._Cuc = nav->eph[i].cuc;
	eph._e = nav->eph[i].e;
	eph._Cus = nav->eph[i].cus;
	eph._sqrt_A = sqrt(nav->eph[i].A);
	eph._Cic = nav->eph[i].cic;
	eph._OMEGA0 = nav->eph[i].OMG0;
	eph._Cis = nav->eph[i].cis;
	eph._i0 = nav->eph[i].i0;
	eph._Crc = nav->eph[i].crc;
	eph._omega = nav->eph[i].omg;
	eph._OMEGADOT = nav->eph[i].OMGd;
	eph._IDOT = nav->eph[i].idot;
	eph._TGD = nav->eph[i].tgd[0];
	eph._health = nav->eph[i].svh;
	rtrover_putGPSEphemeris(&eph);
      }
    }
    if (satsys(nav->eph[i].sat,NULL)==SYS_GAL)
    {  
      //only INAV eph Gal
      if ((isat >= MINPRNGAL) && (isat <= MAXPRNGAL) && (nav->eph[i].code == 1))
      {
	rtrover_ephGal eph;
	eph._satellite._system = 'E';
	eph._satellite._number = isat;
	rtklib2bncTime(&eph._TOC, &nav->eph[i].toc);
	rtklib2bncTime(&eph._TOE, &nav->eph[i].toe);
	eph._IODNav = nav->eph[i].iode;
	eph._clock_bias = nav->eph[i].f0;  
	eph._clock_drift = nav->eph[i].f1;
	eph._clock_driftrate = nav->eph[i].f2;
	eph._Crs = nav->eph[i].crs;
	eph._Delta_n = nav->eph[i].deln;
	eph._M0 = nav->eph[i].M0;
	eph._Cuc = nav->eph[i].cuc;
	eph._e = nav->eph[i].e;
	eph._Cus = nav->eph[i].cus;
	eph._sqrt_A = sqrt(nav->eph[i].A);
	eph._Cic = nav->eph[i].cic;
	eph._OMEGA0 = nav->eph[i].OMG0;
	eph._Cis = nav->eph[i].cis;
	eph._i0 = nav->eph[i].i0;
	eph._Crc = nav->eph[i].crc;
	eph._omega = nav->eph[i].omg;
	eph._OMEGADOT = nav->eph[i].OMGd;
	eph._IDOT = nav->eph[i].idot;
	//eph._BGD_1_5A = nav->eph[i].tgd[0];
	//eph._BGD_1_5B = nav->eph[i].tgd[1];
	eph._BGD_1_5A = nav->eph[i].tgd[1];
	eph._BGD_1_5B = nav->eph[i].tgd[2];
	eph._E5aHS = nav->eph[i].svh;
	eph._E5bHS = nav->eph[i].svh;
	rtrover_putGalEphemeris(&eph);
      }
    }
    if (satsys(nav->eph[i].sat,NULL)==SYS_CMP)
    {   
      if ((isat >= MINPRNCMP) && (isat <= MAXPRNCMP))
      {
	rtrover_ephBds eph;
	eph._satellite._system = 'C';
	eph._satellite._number = isat;
	rtklib2bncTime(&eph._TOC, &nav->eph[i].toc);
	rtklib2bncTime(&eph._TOE, &nav->eph[i].toe);
	eph._IODE = nav->eph[i].iode;
	eph._IODC = nav->eph[i].iodc;
	eph._clock_bias = nav->eph[i].f0;  
	eph._clock_drift = nav->eph[i].f1;
	eph._clock_driftrate = nav->eph[i].f2;
	eph._Crs = nav->eph[i].crs;
	eph._Delta_n = nav->eph[i].deln;
	eph._M0 = nav->eph[i].M0;
	eph._Cuc = nav->eph[i].cuc;
	eph._e = nav->eph[i].e;
	eph._Cus = nav->eph[i].cus;
	eph._sqrt_A = sqrt(nav->eph[i].A);
	eph._Cic = nav->eph[i].cic;
	eph._OMEGA0 = nav->eph[i].OMG0;
	eph._Cis = nav->eph[i].cis;
	eph._i0 = nav->eph[i].i0;
	eph._Crc = nav->eph[i].crc;
	eph._omega = nav->eph[i].omg;
	eph._OMEGADOT = nav->eph[i].OMGd;
	eph._IDOT = nav->eph[i].idot;
	eph._TGD_1_3 = nav->eph[i].tgd[0];
	eph._TGD_2_3 = nav->eph[i].tgd[1];
	eph._health = nav->eph[i].svh;
        eph._toes = nav->eph[i].toes;
	rtrover_putBdsEphemeris(&eph);
      }
    }
  }
  for (int i=0;i<nav->ng;i++)
  {
    char id[4];
    *id='\0';
    satno2id(nav->geph[i].sat, id);
    int isat=atoi(id+1);
    if ((isat >= MINPRNGLO) && (isat <= MAXPRNGLO) && (nav->geph[i].toe.time))
    {
      rtrover_ephGlo eph;
      eph._satellite._system = 'R';
      eph._satellite._number = isat;
      gtime_t toe_gpst = nav->geph[i].toe;
      gtime_t toe_utct = gpst2utc(toe_gpst);
      rtklib2bncTime(&eph._timeUTC, &toe_utct);
      eph._gps_utc = (int)timediff (toe_gpst, toe_utct);
      eph._E = nav->geph[i].age;
      eph._tau = nav->geph[i].taun;
      eph._gamma = nav->geph[i].gamn;
      eph._x_pos = nav->geph[i].pos[0]/1000.0;
      eph._x_velocity = nav->geph[i].vel[0]/1000.0;
      eph._x_acceleration = nav->geph[i].acc[0]/1000.0;
      eph._y_pos = nav->geph[i].pos[1]/1000.0;
      eph._y_velocity = nav->geph[i].vel[1]/1000.0;
      eph._y_acceleration = nav->geph[i].acc[1]/1000.0;
      eph._z_pos = nav->geph[i].pos[2]/1000.0;
      eph._z_velocity = nav->geph[i].vel[2]/1000.0;
      eph._z_acceleration = nav->geph[i].acc[2]/1000.0;
      eph._health = nav->geph[i].svh;
      eph._frequency_number = nav->geph[i].frq;
      slot[nav->geph[i].sat-1] = nav->geph[i].frq;
      rtrover_putGloEphemeris(&eph);
    }
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////
void putOrbCorrections(ssr_t *ssr)
{
  for (int i=0;i<MAXSAT;i++)
  {
    char id[4];
    satno2id(i+1, id);
    if (ssr[i].update && ssr[i].t0[0].time)
    {
      if ((id[0] != 'G') && (id[0] != 'R') && (id[0] != 'E') && (id[0] != 'C'))
        continue;
      rtrover_orbCorr corr;
      corr._satellite._system = id[0];
      corr._satellite._number = atoi(id+1);
      corr._iod = ssr[i].iode;
      if (id[0] == 'C')
        corr._iod = ssr[i].iodcrc;
      rtklib2bncTime(&corr._time, &ssr[i].t0[0]);
      for (int j=0;j<3;j++)
      {
        corr._rao[j] = ssr[i].deph[j];
        corr._dotRao[j] = ssr[i].ddeph[j];
        corr._dotDotRao[j] = 0.0;
      }
      rtrover_putOrbCorrections(1, &corr);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
void putClkCorrections(ssr_t *ssr)
{
  for (int i=0;i<MAXSAT;i++)
  {
    char id[4];
    satno2id(i+1, id);

    if (ssr[i].update && ssr[i].t0[1].time)
    {
      if ((id[0] != 'G') && (id[0] != 'R') && (id[0] != 'E') && (id[0] != 'C'))
        continue;
      rtrover_clkCorr corr;
      corr._satellite._system = id[0];
      corr._satellite._number = atoi(id+1);
      corr._iod = ssr[i].iode;
      if (id[0] == 'C')
        corr._iod = ssr[i].iodcrc;  
      rtklib2bncTime(&corr._time, &ssr[i].t0[1]);
      corr._dClk = ssr[i].dclk[0]/CLIGHT;
      corr._dotDClk = ssr[i].dclk[1]/CLIGHT;
      corr._dotDotDClk = ssr[i].dclk[2]/CLIGHT;
      rtrover_putClkCorrections(1, &corr);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
void putBiases(ssr_t *ssr)
{
  int i,j;

  for (i=0;i<MAXSAT;i++)
  {
    char id[4];
    satno2id(i+1, id);
    if (ssr[i].update && ssr[i].t0[4].time)
    {
      if ((id[0] != 'G') && (id[0] != 'R') && (id[0] != 'E') && (id[0] != 'C'))
        continue;
      rtrover_bias _biases[MAXCODE];
      rtrover_satBiases biases;
      biases._biases = _biases;
      biases._satellite._system = id[0];
      biases._satellite._number = atoi(id+1);
      rtklib2bncTime(&biases._time, &ssr[i].t0[4]);
      biases._numBiases = 0;
      for (j=0;j<MAXCODE;j++)
      {
        if (ssr[i].cbias[j])
        {
          if (j+1 == CODE_L1C)
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='1';
            _biases[biases._numBiases]._rnxType3ch[2]='C';
            _biases[biases._numBiases]._value=ssr[i].cbias[j];
            biases._numBiases++;
          }
          if (j+1 == CODE_L1W)
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='1';
            _biases[biases._numBiases]._rnxType3ch[2]='W';
            _biases[biases._numBiases]._value=ssr[i].cbias[j];
            biases._numBiases++;
          }
	  if (j+1 == CODE_L2P)
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='2';
            _biases[biases._numBiases]._rnxType3ch[2]='P';
            _biases[biases._numBiases]._value=ssr[i].cbias[j];
            biases._numBiases++;
          }
          if (j+1 == CODE_L2X)
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='2';
            _biases[biases._numBiases]._rnxType3ch[2]='X';
            _biases[biases._numBiases]._value=ssr[i].cbias[j];
            biases._numBiases++;
          }
          if (j+1 == CODE_L2W)
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='2';
            _biases[biases._numBiases]._rnxType3ch[2]='W';
            _biases[biases._numBiases]._value=ssr[i].cbias[j];
            biases._numBiases++;
          }
          if ((j+1 == CODE_L5I) || (j+1 == CODE_L5Q) || (j+1 == CODE_L5X))
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='5';
            _biases[biases._numBiases]._rnxType3ch[2]='Q';
            _biases[biases._numBiases]._value=ssr[i].cbias[j];
            biases._numBiases++;
          }
	  if ((j+1 == CODE_L1X) || (j+1 == CODE_L1I) || (j+1 == CODE_L1B) || (j+1 == CODE_L1Q) || (j+1 == CODE_L2I) || (j+1 == CODE_L2Q))
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='1';
            _biases[biases._numBiases]._rnxType3ch[2]='X';
            _biases[biases._numBiases]._value=ssr[i].cbias[j];
            biases._numBiases++;
          }
	  if ((j+1 == CODE_L6I) || (j+1 == CODE_L6Q))
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='6';
            _biases[biases._numBiases]._rnxType3ch[2]='Q';
            _biases[biases._numBiases]._value=ssr[i].cbias[j];
            biases._numBiases++;
          }
	  if ((j+1 == CODE_L7I) || (j+1 == CODE_L7Q) || (j+1 == CODE_L7X))
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='7';
            _biases[biases._numBiases]._rnxType3ch[2]='Q';
            _biases[biases._numBiases]._value=ssr[i].cbias[j];
            biases._numBiases++;
          }
        }
      }
      rtrover_putBiases(1, &biases);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
void putPhaseBiases(ssr_t *ssr)
{
  int i,j;

  for (i=0;i<MAXSAT;i++)
  {
    char id[4];
    satno2id(i+1, id);
    if (ssr[i].update && ssr[i].t0[4].time)
    {
      if ((id[0] != 'G') && (id[0] != 'R') && (id[0] != 'E') && (id[0] != 'C'))
        continue;
      rtrover_phaseBias _biases[MAXCODE];
      rtrover_satPhaseBiases biases;
      biases._biases = _biases;
      biases._satellite._system = id[0];
      biases._satellite._number = atoi(id+1);
      biases._yawAngleValue = ssr[i].yawAngle;
      biases._dispersive=false;
      biases._mwConsistency=false;
      rtklib2bncTime(&biases._time, &ssr[i].t0[4]);
      biases._numBiases = 0;
      for (j=0;j<MAXCODE;j++)
      {
        _biases[j]._discontinuityValue = 0;
        if (ssr[i].pbias[j])
        {
          if ((j+1 == CODE_L1C) || (j+1 == CODE_L1B) || (j+1 == CODE_L1X))  // Phase biases L1
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='1';
            _biases[biases._numBiases]._rnxType3ch[2]='Z';
            _biases[biases._numBiases]._value=ssr[i].pbias[j];
            _biases[biases._numBiases]._discontinuityValue=ssr[i].discontinuity[j];
	    _biases[biases._numBiases]._intIndicator=ssr[i].n1i[j]!=0;
            _biases[biases._numBiases]._wlIntIndicator=ssr[i].wli[j];
            biases._numBiases++;
          }
	  if ((j+1 == CODE_L2I)|| (j+1 == CODE_L1I)) // Phase biases Bds B1
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='1';
            _biases[biases._numBiases]._rnxType3ch[2]='Z';
            _biases[biases._numBiases]._value=ssr[i].pbias[j];
            _biases[biases._numBiases]._discontinuityValue=ssr[i].discontinuity[j];
	    _biases[biases._numBiases]._intIndicator=ssr[i].n1i[j]!=0;
            _biases[biases._numBiases]._wlIntIndicator=ssr[i].wli[j];
            biases._numBiases++;
          }
          if (j+1 == CODE_L2W)  // Phase biases L2
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='2';
            _biases[biases._numBiases]._rnxType3ch[2]='Z';
            _biases[biases._numBiases]._value=ssr[i].pbias[j];
            _biases[biases._numBiases]._discontinuityValue=ssr[i].discontinuity[j];
	    _biases[biases._numBiases]._intIndicator=ssr[i].n1i[j]!=0;
            _biases[biases._numBiases]._wlIntIndicator=ssr[i].wli[j];
            biases._numBiases++;
          }
          if ((j+1 == CODE_L5I) || (j+1 == CODE_L5Q) || (j+1 == CODE_L5X)) // Phase biases L5
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='5';
            _biases[biases._numBiases]._rnxType3ch[2]='Z';
            _biases[biases._numBiases]._value=ssr[i].pbias[j];
            _biases[biases._numBiases]._discontinuityValue=ssr[i].discontinuity[j];
	    _biases[biases._numBiases]._intIndicator=ssr[i].n1i[j]!=0;
            _biases[biases._numBiases]._wlIntIndicator=ssr[i].wli[j];
            biases._numBiases++;
          }
	  if ((j+1 == CODE_L6I) || (j+1 == CODE_L6Q)) // Phase biases L6 / B3
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='6';
            _biases[biases._numBiases]._rnxType3ch[2]='Z';
            _biases[biases._numBiases]._value=ssr[i].pbias[j];
            _biases[biases._numBiases]._discontinuityValue=ssr[i].discontinuity[j];
	    _biases[biases._numBiases]._intIndicator=ssr[i].n1i[j]!=0;
            _biases[biases._numBiases]._wlIntIndicator=ssr[i].wli[j];
            biases._numBiases++;
          }
	  if ((j+1 == CODE_L7I) || (j+1 == CODE_L7Q) || (j+1 == CODE_L7X)) // Phase biases L7 / B2
          {
            _biases[biases._numBiases]._rnxType3ch[0]='C';
            _biases[biases._numBiases]._rnxType3ch[1]='7';
            _biases[biases._numBiases]._rnxType3ch[2]='Z';
            _biases[biases._numBiases]._value=ssr[i].pbias[j];
            _biases[biases._numBiases]._discontinuityValue=ssr[i].discontinuity[j];
	    _biases[biases._numBiases]._intIndicator=ssr[i].n1i[j]!=0;
            _biases[biases._numBiases]._wlIntIndicator=ssr[i].wli[j];
            biases._numBiases++;
          }
          //fprintf(stderr, "PHASEP %d %d %lf\n", biases._satellite._number, j+1, ssr[i].pbias[j]);
        }
      }
      rtrover_putPhaseBiases(1, &biases);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
void decode_stream_ssr(ssr_t *ssr)
{
  putOrbCorrections(ssr);
  putClkCorrections(ssr);
  putBiases(ssr);
  putPhaseBiases(ssr);
}

//////////////////////////////////////////////////////////////////////////////////////////////
void decode_stream_sbsmsg(sbsmsg_t *sbsmsg)
{
  int type=getbitu(sbsmsg->msg,8,6);
  int i=sbsupdatecorr(sbsmsg, &navIono);
  /*
  if (type == 18)
    printf("sbsmsg %d %d\n", type, i);
  */
}

//////////////////////////////////////////////////////////////////////////////////////////////
void decode_tropo_iono(char *chaine)
{
  char *pChaine;
  //fprintf(stderr, "%s", sLigneFichier); 
  rtrover_tropoCorr tropoCorr;
  
  /////////////////////
  //read line
  /////////////////////
  
  //if errors
  pChaine = strtok(chaine, " "); if (!pChaine) return; //n 
  pChaine = strtok(NULL, " "); if (!pChaine) return; //date day
  pChaine = strtok(NULL, " "); if (!pChaine) return; //date second
  pChaine = strtok(NULL, " "); if (!pChaine) return; //tropo
  
  tropoCorr._value = atof(pChaine);
  //fprintf(stderr, "TROPO: %lf\n", tropoCorr._value); 

  pChaine = strtok (NULL, " "); if (!pChaine) return; 
  while (pChaine != NULL)
  {   
    //traitement	  
    if (!strncmp(pChaine,"G",1))
    {
      int isat;
      isat = atoi(1+pChaine)-1;
      tabIonoCorr[GPS][isat]._staName = NULL;
      tabIonoCorr[GPS][isat]._flag = 1;
      strncpy(&tabIonoCorr[GPS][isat]._satellite._system, pChaine, 1);
      tabIonoCorr[GPS][isat]._satellite._number = isat+1;
      pChaine = strtok (NULL, " "); if (!pChaine) return;
      tabIonoCorr[GPS][isat]._value = atof(pChaine);
      //fprintf(stderr, "IONO %c%d %lf\n", ionoCorr[isat]._satellite._system, ionoCorr[isat]._satellite._number, ionoCorr[isat]._value);
    }
    pChaine = strtok (NULL, " "); //fin de ligne si pChaine == null
  }
  rtrover_putTropoCorrs(&tropoCorr);
}

//////////////////////////////////////////////////////////////////////////////////////////////
void decode_char(int n, int format, int sec, unsigned char d)
{
  int j;
  int lastSize=0;
  static char sLigneFichier[200]="";

  //read tropo/iono corrections file
  if (format == 99)
  {   
    strcat(sLigneFichier, (const char*)&d); 
    if (d == '\n')
    {
      decode_tropo_iono(sLigneFichier);
      strcpy(sLigneFichier, ""); //raz buffer
    }
  }
  if (format == STRFMT_RTCM3)
  {
    if (n>=(int)rtcm.size())
    {      
      lastSize=rtcm.size();
      rtcm.resize(n+1);
      for(int i=lastSize; i<(int)rtcm.size() ; i++)
      {
        init_rtcm(&rtcm[i]);
        strcpy(rtcm[i].opt,"-EPHALL");
      }
    }
    if (!rtcm[n].time.time)
    {
      rtcm[n].time.time=sec;
      rtcm[n].time.sec=0.0;
    }
    j=input_rtcm3(&rtcm[n], (unsigned char)d);
    //if ((j == 1) && !n)
    if (j == 1)
      decode_stream_obs(&rtcm[n].obs, 0);
    if (j == 2)
      decode_stream_eph(&rtcm[n].nav);
    if (j == 10)
      decode_stream_ssr(rtcm[n].ssr);
    if (j == 20){
      vtec_time=rtcm[n].time;
      vtec=rtcm[n].vtec;
    }
  }
  if ((format == STRFMT_UBX) || (format == STRFMT_NVS))
  {
    if (n>=(int)raw.size())
    {
      lastSize=raw.size();
      raw.resize(n+1);
      for(int i=lastSize; i<(int)raw.size() ; i++)
      {
        init_raw(&raw[i]);
      }
    }
    /*
    if (!raw[n].time.time)
    {
      raw[n].time.time=sec;
      raw[n].time.sec=0.0;
    }
    */
    j=input_raw(&raw[n], format, (unsigned char)d);
    //if ((j == 1) && !n)
    if (j == 1)
      decode_stream_obs(&raw[n].obs, 1);
    if (j == 2)
      decode_stream_eph(&raw[n].nav);
    if (j == 3)
      decode_stream_sbsmsg(&raw[n].sbsmsg);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
void parseArgs(int argc, char *argv[], char *confFileName, char *roverfileName, char *dcbFileName, bool *lowlevel)
{
  int iArg = 1; 
  
  if (argc < 7 || argc > 9)
  {
    fprintf(stderr, "Usage: processStream -conf conf.txt -rover rover.txt -dcb DCB [-verbose] [-lowlevel]\n");
    exit(1);
  }
  
  while (iArg < argc)
  {
    if (strcmp(argv[iArg], "-conf") == 0)
    {
      iArg++;
      strcpy(confFileName, argv[iArg]);
    }
    else if (strcmp(argv[iArg], "-rover") == 0)
    {
      iArg++;
      strcpy(roverfileName, argv[iArg]);
    }
    else if (strcmp(argv[iArg], "-dcb") == 0)
    {
      iArg++;
      strcpy(dcbFileName, argv[iArg]);
    }
    else if (strcmp(argv[iArg], "-verbose") == 0)
    {
      verbose = true;
    }
    else if (strcmp(argv[iArg], "-lowlevel") == 0)
    {
      *lowlevel = true;
    }
    else
    {
      fprintf(stderr, "Usage: processStream -conf conf.txt -rover rover.txt -dcb DCB [-verbose] [-lowlevel]\n");
      exit(1);
    }
    iArg++;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{ 
  char confFileName[256];
  char roverfileName[256];  
  char dcbFileName[256];    
  bool lowlevel = false;
  
  verbose = false;
  
  parseArgs(argc, argv, confFileName, roverfileName, dcbFileName, &lowlevel);
  
  decode_begin(confFileName, roverfileName, dcbFileName, lowlevel);

  while (!feof(stdin))
  {
    char s[200], sh[150];
    int n, format, sec;
    *s='\0';
    fgets(s, sizeof s, stdin);
    if (*s)
    {
      sscanf(s, "%d%d%d%s", &n, &format, &sec, sh);
      if (n<1)
        continue;
      n--;
      // Set Current Rover
      if (n<nbRover)
      {
        currentRover=n;
        rtrover_setCurrentRover(currentRover);
      }
      // decode
      for (int i=0;i<(int)strlen(sh);i+=2)
      {
        unsigned char str[3];
        unsigned int d;
        strncpy((char *)str,sh+i,2);
        str[2]='\0';
        sscanf((char *)str,"%X", &d);	
        decode_char(n, format, sec, (unsigned char)d);
      }
    }
  }
  decode_end();

  return 1;
}
