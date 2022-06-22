#include <cstdio>
#include <cstring>
#include <unistd.h>
#include <fstream>
#include <sstream>

#include "rtrover_ppp.h"
#include "rtrover_frequency.h"
#include "rtrover_utility.h"

#define EPS 0.001
#define EPSDT (1.e-6)

class ROVER
{
  public:
    ROVER(std::string name, A_VECTOR3D aprPos);
    std::string _name;
    A_VECTOR3D _aprPos;
    double _tropoCorr;
    A_DATE _date;
    A_DATE _datePPP;
    A_PPP _ppp;
    std::vector<std::vector<A_MEASUREMENT> > _measurements;
    std::vector<std::vector<A_BIAS> > _bias;
    std::vector<std::vector<A_SSR_PARAMETER> > _positions;
    void clear();
};

static std::vector<ROVER> listRover;
static int outputVerbose;
static bool verbose;
static char cSyst[4]={'G','R','E','C'};

//////////////////////////////////////////////////////////////////////////////////////////////
// Class ROVER

ROVER::ROVER(std::string name, A_VECTOR3D aprPos)
{
  _name = name;
  _aprPos = aprPos;
  _tropoCorr = 0.0;
  _measurements = std::vector<std::vector<A_MEASUREMENT> > (SystMax);
  _bias = std::vector<std::vector<A_BIAS> > (SystMax);
  _positions = std::vector<std::vector<A_SSR_PARAMETER> >(SystMax);
  for (int isyst=0; isyst < SystMax; isyst++)
  {
    const int maxSatSyst = getMaxSatSyst(isyst);
    _measurements[isyst] = std::vector<A_MEASUREMENT>(maxSatSyst);
    _bias[isyst] = std::vector<A_BIAS>(maxSatSyst);
    _positions[isyst] = std::vector<A_SSR_PARAMETER>(maxSatSyst);
  }
}

void ROVER::clear()
{
  for (int isyst=0; isyst < SystMax; isyst++)
  {
    const int maxSatSyst = getMaxSatSyst(isyst);
    _measurements[isyst].clear();
    _measurements[isyst].resize(maxSatSyst);
    _bias[isyst].clear();
    _bias[isyst].resize(maxSatSyst);
    _positions[isyst].clear();
    _positions[isyst].resize(maxSatSyst);
  }
  _tropoCorr=0.0;
}


#define BNC_FORMAT
//#undef BNC_FORMAT

//////////////////////////////////////////////////////////////////////////////////////////////
// Write Output
static void writeOutput(A_RECEIVER_PRODUCT &recProduct, std::vector<std::vector<A_MEASUREMENT_PRODUCT> > measProduct, int numRover)
{
//14-07-17 14:51:58.000 FFMJ PPP  9  0  0  0  4053456.858 +- 8.703   617729.426 +- 7.618  4869396.227 +- 9.560 2.344938 + 0.006150 +- 0.499397
  A_DATE &date=listRover[numRover-1]._datePPP;

#ifdef BNC_FORMAT
  int _day, Mois, An, Heure, Minute;
  double _second;
  
  date.calendar(&_day, &Mois, &An, &Heure, &Minute, &_second);
  if(outputVerbose)
  {
    fprintf(stdout, "%02d-%02d-%02d %02d:%02d:%06.3lf %6.3f %6.3f %6.3f %6.3f %s PPP %2d %2d %2d %2d %12.3f +- %.3f %12.3f +- %.3f %12.3f +- %.3f %5.6f + %5.6f +- %.6f  %.4f\n",
	      An%100, Mois, _day, Heure, Minute, _second,
	      recProduct._HstaGPS, recProduct._HstaGlo,
	      recProduct._HstaGal, recProduct._HstaBds,
	      listRover[numRover-1]._name.c_str(),
	      recProduct._nbMes, recProduct._nbMesBloEx, recProduct._nbMesBloNw, recProduct._nbMesBloN1,
	      recProduct._x, recProduct._covX,
	      recProduct._y, recProduct._covY,
	      recProduct._z, recProduct._covZ,
	      recProduct._tropoV, recProduct._tropo, recProduct._covTropo, recProduct._integrity);
  }
  else
  {
    fprintf(stdout, "%02d-%02d-%02d %02d:%02d:%06.3lf %s PPP %2d %2d %2d %2d %12.3f +- %.3f %12.3f +- %.3f %12.3f +- %.3f %5.6f + %5.6f +- %.6f  %.4f\n",
	      An%100, Mois, _day, Heure, Minute, _second,
	      listRover[numRover-1]._name.c_str(),
	      recProduct._nbMes, recProduct._nbMesBloEx, recProduct._nbMesBloNw, recProduct._nbMesBloN1,
	      recProduct._x, recProduct._covX,
	      recProduct._y, recProduct._covY,
	      recProduct._z, recProduct._covZ,
	      recProduct._tropoV, recProduct._tropo, recProduct._covTropo, recProduct._integrity);
  }	      	      
  fflush (stdout);
#else  
  double xGps = recProduct._HstaGPS;
  if (xGps && (fabs(xGps) < EPS))
    xGps = EPS;
  double xGlo = recProduct._HstaGlo;
  if (xGlo && (fabs(xGlo) < EPS))
    xGlo = EPS;
  double xGal = recProduct._HstaGal;
  if (xGal && (fabs(xGal) < EPS))
    xGal = EPS;
  double xBds = recProduct._HstaBds;
  if (xBds && (fabs(xBds) < EPS))
    xBds = EPS;
  fprintf(stdout, "1 %d %.1f %s %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.3f %.6f %.4f %.4f %.4f %.6f %.6f %.6f %d %d %d %d\n",
          date.getDay(), date.getSecond(), listRover[numRover-1]._name.c_str(), 
	  xGps, recProduct._covHstaGPS, xGlo, recProduct._covHstaGlo,
	  xGal, recProduct._covHstaGal, xBds, recProduct._covHstaBds,
          recProduct._tropoV, recProduct._tropo, recProduct._covTropo,
          recProduct._x, recProduct._y, recProduct._z,
          recProduct._covX, recProduct._covY, recProduct._covZ,
          recProduct._nbMes, recProduct._nbMesBloEx, recProduct._nbMesBloNw, recProduct._nbMesBloN1);
  fflush (stdout);
#endif
  if (verbose) {
    fprintf(stderr, "%d %d %lf %.3lf ", numRover, listRover[numRover-1]._date.getDay(), listRover[numRover-1]._date.getSecond(), recProduct._tropoV + recProduct._tropo);
    for (int isyst=0; isyst<SystMax;isyst++)
      {
	for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
	  {
/*
    double a, ca;
    if (measProduct[isyst][isat]._ambP)
      {
      double a = measProduct[isyst][isat]._ambP; if (a && (fabs(a) < 0.001)) a=0.001;
      double ac = measProduct[isyst][isat]._covAmbP; if (ac && (fabs(ac) < 0.001)) ac=0.001;
      fprintf(stderr, "AMB %d %d %d %lf %lf %lf\n", numRover, isyst*37+isat+1, listRover[numRover-1]._date.getDay(), listRover[numRover-1]._date.getSecond(),
	  a, ac);
      }
    if (measProduct[isyst][isat]._iono)
      fprintf(stderr, "IONO %d %d %d %lf %lf %lf\n", numRover, isyst*37+isat+1, listRover[numRover-1]._date.getDay(), listRover[numRover-1]._date.getSecond(),
	  measProduct[isyst][isat]._iono, measProduct[isyst][isat]._covIono);
*/
	    if (measProduct[isyst][isat]._ambP && measProduct[isyst][isat]._iono)
	      {
		///double amb = measProduct[isyst][isat]._ambP;
		double covamb = measProduct[isyst][isat]._covAmbP;
		double iono = measProduct[isyst][isat]._iono;
		if (covamb && (covamb < 0.01))
		  fprintf(stderr, "%c%02d %.3lf %.3lf ", cSyst[isyst], +isat+1, iono, covamb);
	      }
	  }
      }
    fprintf(stderr, "\n");
    fflush(stderr);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////
// process all measurements for one date
static int processEpoch(int numRover)
{
  std::vector<std::vector<A_MEASUREMENT> > emptyMeas(SystMax);
  std::vector<std::vector<A_BIAS> > emptyBias(SystMax);
  std::vector<std::vector<A_SSR_PARAMETER> > emptyPos(SystMax);
  std::vector<std::vector<A_VECTOR3D> > emptyVel(SystMax);
  double dt;
  A_RECEIVER_PRODUCT recProduct;
  std::vector<std::vector<A_MEASUREMENT_PRODUCT> > measProduct(SystMax);
  
  ROVER &rover=listRover[numRover-1];
  
  //init vector
  for (int isyst=0; isyst<SystMax;isyst++)
  {
    const int maxSatSyst = getMaxSatSyst(isyst);
    emptyMeas[isyst] = std::vector<A_MEASUREMENT>(maxSatSyst);
    emptyBias[isyst] = std::vector<A_BIAS>(maxSatSyst);
    emptyPos[isyst] = std::vector<A_SSR_PARAMETER>(maxSatSyst);
    measProduct[isyst] = std::vector<A_MEASUREMENT_PRODUCT>(maxSatSyst);
  }
  
  if (rover._datePPP.getDay())
  {
    while (1) 
    {      
      dt = rover._datePPP + A_PPP::_settings._dt - rover._date;           
      if (dt >= 0.0 || (fabs(dt) <EPSDT)) 
      {	
        break;	  
      }
      rover._datePPP = rover._datePPP + A_PPP::_settings._dt;  
      if (rover._ppp.computePPP(emptyMeas, emptyBias, emptyPos, rover._tropoCorr,rover._datePPP, recProduct, measProduct, rover._aprPos, numRover-1))
        writeOutput(recProduct, measProduct, numRover);
    }
  }
  
  double fpas = rover._date.getSecond()/A_PPP::_settings._dt;
  double fracf = fabs(fpas - floor(fpas+0.5));
  if (((rover._datePPP.getDay()==0) || (fabs(rover._date - rover._datePPP - A_PPP::_settings._dt) < EPSDT)) && (fracf < EPSDBL)) 
  {
  
  	rover._datePPP = rover._date; 
  	if (rover._ppp.computePPP(rover._measurements, rover._bias, rover._positions, rover._tropoCorr, rover._datePPP, recProduct, measProduct, rover._aprPos, numRover-1))
  	{
    	writeOutput(recProduct, measProduct, numRover);
    	return 1;
  	}
  	return 0;
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// process one measurement
static void processMeasurement(int numRover, int isyst, int isat, A_DATE& date, A_MEASUREMENT& meas, A_BIAS& bias,
                       A_SSR_PARAMETER& pos, double tropo)
{
  static int lastRover = 1;

  if (numRover!=lastRover)
  {
    if (listRover[lastRover-1]._date.getDay())
    {
      processEpoch(lastRover);    
      listRover[lastRover-1].clear();
    }
  }
  if ((isat >= -1) && (isat < getMaxSatSyst(isyst)))
  {
    if ((numRover==lastRover) && (date.getSecond() != listRover[numRover-1]._date.getSecond()))
    {
      if (listRover[numRover-1]._date.getDay())
      {
        processEpoch(numRover);	
	listRover[numRover-1].clear();
      }
    }
    if ((isat != -1) && (isyst != -1))
    {
      listRover[numRover-1]._measurements[isyst][isat] = meas;
      listRover[numRover-1]._bias[isyst][isat] = bias;
      listRover[numRover-1]._positions[isyst][isat] = pos;
      listRover[numRover-1]._tropoCorr=tropo;
      listRover[numRover-1]._date.setDay(date.getDay());
      listRover[numRover-1]._date.setSecond(date.getSecond());
    }    
  }
  lastRover = numRover;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Read one measurement
static int readMeasurement(char *ligne, int *numRover, int *isyst, int *iSat, A_DATE *date, A_MEASUREMENT *meas, A_BIAS *bias, A_SSR_PARAMETER *pos, double *tropo)
{
  int iNb;
  int day;
  double second;
  char syst;
  char systSat[3];
  int n1i;

  if (*ligne == '%')
    return 0;

  iNb = sscanf(ligne,
  "%d%s%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%d%lf%lf%lf%lf%lf%lf%lf%lf%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%d%d%d", numRover, &systSat, &day, &second,
                      &meas->_code[F1], &meas->_code[F2], &meas->_code[F6], &meas->_code[F5], &meas->_code[F7],
                      &meas->_phase[F1], &meas->_phase[F2], &meas->_phase[F6], &meas->_phase[F5], &meas->_phase[F7],
                      &meas->_doppler[F1], &meas->_doppler[F2], &meas->_doppler[F6], &meas->_doppler[F5], &meas->_doppler[F7],&meas->_slot_typeBds,
                      &pos->_Xsp3, &pos->_Ysp3, &pos->_Zsp3, &pos->_Hsp3, &pos->_VXsp3, &pos->_VYsp3, &pos->_VZsp3, &pos->_iono,&pos->_typeIono,tropo,&pos->_yaw,
                      &bias->_code[F1], &bias->_code[F2], &bias->_code[F6], &bias->_code[F5], &bias->_code[F7],
                      &bias->_phase[F1], &bias->_phase[F2], &bias->_phase[F6], &bias->_phase[F5], &bias->_phase[F7],
                      &n1i,&bias->_wli, &bias->_discontinuity);

		  
  if (iNb != 44)
    return 0;
    
  iNb = sscanf(systSat, "%c%d", &syst, iSat);
  if (iNb != 2)
    return 0;
    
  bias->_n1i=n1i!=0;
    
  date->setDay(day);
  date->setSecond(second);
  (*iSat)--;
  *isyst=convertSystem(syst);
  
  for(int f=F1; f<FMax; f++)
    meas->_doppler[f] *= -FREQUENCY::LAMBDA((Frequency)f, *isyst, meas->_slot_typeBds);
  return 1;	
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Read rover file
static void readRover(char *roverFileName)
{
  std::ifstream file(roverFileName, std::ios::in);
 
  if(file)
  {
    std::string ligne;
    while(getline(file, ligne))
    {
      std::istringstream iss(ligne);
      std::string roverName;
      double x, y, z;
      iss >> roverName >> x >> y >> z;
      A_VECTOR3D aprPos(x,y,z);
      ROVER rover(roverName, aprPos);
      listRover.push_back(rover);
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
// Read conf file
static void readParam(char *confFile)
{
  FILE *file;
  char confMode[20], dum[200];
  double double1,double2,double3;
  int int1;
  int nbSatFixAmb = 0;
  char atxFileName[1024];

  file=fopen(confFile,"r");
  fscanf(file, "%s%s\n", confMode, dum);
  A_MODE &mode = A_PPP::_settings._mode;
  if (!strcmp(confMode, "mode_PPP_AR"))
  {  ///< precise point positining with ambiguity resolution
      mode._code = 1;
      mode._phase = 1;
      mode._freq[F1] = 1;
      mode._freq[F2] = 1;
      mode._freq[F5] = 1;
      mode._freq[F6] = 1;
      mode._freq[F7] = 1;
      mode._ambig = 1;
  }
  if (!strcmp(confMode, "mode_PPP_DF"))
  {///< dual-frequency precise point positioning
      mode._code = 1;
      mode._phase = 1;
      mode._freq[F1] = 1;
      mode._freq[F2] = 1;
      mode._freq[F5] = 1;
      mode._freq[F6] = 1;
      mode._freq[F7] = 1;
      mode._ambig = 0;
  }
  if (!strcmp(confMode, "mode_PPP_SF"))
  {///< single-frequency precise precise point positioning
      mode._code = 1;
      mode._phase = 1;
      mode._freq[F1] = 1;
      mode._freq[F2] = 0;
      mode._freq[F5] = 0;
      mode._freq[F6] = 0;
      mode._freq[F7] = 0;
      mode._ambig = 0;
  }
  if (!strcmp(confMode, "mode_SPP_DF"))
  {///< dual-frequency single point positioning
      mode._code = 1;
      mode._phase = 1;
      mode._freq[F1] = 1;
      mode._freq[F2] = 1;
      mode._freq[F5] = 1;
      mode._freq[F6] = 1;
      mode._freq[F7] = 1;
      mode._ambig = 0;
  }
  if (!strcmp(confMode, "mode_SPP_SF"))
  { ///< single-frequency single point positioning
      mode._code = 1;
      mode._phase = 1;
      mode._freq[F1] = 1;
      mode._freq[F2] = 0;
      mode._freq[F5] = 0;
      mode._freq[F6] = 0;
      mode._freq[F7] = 0;
      mode._ambig = 0;
  }
  fscanf(file, "%s%s\n", atxFileName, dum);
  fscanf(file, "%d%d%d%s\n", &A_FILTER::_filterSetting._fixN1, &A_FILTER::_filterSetting._fixNw, &A_FILTER::_filterSetting._fixEx, dum);
  fscanf(file, "%d%s\n", &A_FILTER::_filterSetting._gps._use, dum);
  fscanf(file, "%d%s\n", &A_FILTER::_filterSetting._glo._use, dum);
  fscanf(file, "%d%s\n", &A_FILTER::_filterSetting._gal._use, dum);
  fscanf(file, "%d%s\n", &A_FILTER::_filterSetting._bds._use, dum);
  fscanf(file, "%d%s\n", &int1, dum); //sbascorrection
  fscanf(file, "%d%s\n", &A_PPP::_settings._reset, dum);
  fscanf(file, "%d%s\n", &outputVerbose, dum);
  fscanf(file, "%lf%s\n", &A_PPP::_settings._dt, dum);
  fscanf(file, "%lf%s\n", &double1, dum); //maxAge
  fscanf(file, "%d%s\n", &A_FILTER::_filterSetting._nbMin, dum);
  fscanf(file, "%d%s\n", &A_FILTER::_filterSetting._maxElim, dum);
  fscanf(file, "%d%s\n", &A_FILTER::_filterSetting._raim, dum);
  fscanf(file, "%lf%s\n", &A_PPP::_settings._thrMap, dum);
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._sigIniTro = double1 * double1;
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._sigModTro = double1 * double1;
  fscanf(file, "%d%s\n", &nbSatFixAmb, dum);
  A_FILTER::_filterSetting._nbSatFixAmb = (nbSatFixAmb >= 0) ? nbSatFixAmb : 0;
  fscanf(file, "%lf%s\n", &A_FILTER::_filterSetting._thrAmb, dum);
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._sigIniBiasClk = double1 * double1;
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._sigModBiasClk = double1 * double1;
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._sigIniIono = double1 * double1;
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._sigModIono = double1 * double1;
  fscanf(file, "%lf%lf%lf%s\n", &double1,&double2,&double3, dum);
  A_FILTER::_filterSetting._sigMeasIono[0] = double1 * double1;
  A_FILTER::_filterSetting._sigMeasIono[1] = double2 * double2;
  A_FILTER::_filterSetting._sigMeasIono[2] = double3 * double3;
  fscanf(file, "%lf%s\n", &A_FILTER::_filterSetting._thrMeasIono, dum);
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._sigMeasTropo = double1 * double1;
  fscanf(file, "%lf%s\n", &A_FILTER::_filterSetting._thrMeasTropo, dum);
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._sigIniPos = double1 * double1;
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._sigModPos = double1 * double1;
  fscanf(file, "%d%s\n", &A_FILTER::_filterSetting._dtMax, dum);
  fscanf(file, "%lf%s\n", &A_FILTER::_filterSetting._thrMeasCode, dum);
  fscanf(file, "%lf%s\n", &A_FILTER::_filterSetting._thrMeasPhase, dum);
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._gps._sigCode[0] = double1 * double1;
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._gps._sigPhase[0] = double1 * double1;
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._glo._sigCode[0] = double1 * double1;
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._glo._sigPhase[0] = double1 * double1;
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._gal._sigCode[0] = double1 * double1;
  fscanf(file, "%lf%s\n", &double1, dum);
  A_FILTER::_filterSetting._gal._sigPhase[0] = double1 * double1;
  fscanf(file, "%lf%lf%lf%s\n", &double1, &double2, &double3, dum);
  A_FILTER::_filterSetting._bds._sigCode[0] = double1 * double1;
  A_FILTER::_filterSetting._bds._sigCode[1] = double2 * double2;
  A_FILTER::_filterSetting._bds._sigCode[2] = double3 * double3;
  fscanf(file, "%lf%lf%lf%s\n", &double1, &double2, &double3, dum);
  A_FILTER::_filterSetting._bds._sigPhase[0] = double1 * double1;
  A_FILTER::_filterSetting._bds._sigPhase[1] = double2 * double2;
  A_FILTER::_filterSetting._bds._sigPhase[2] = double3 * double3;
  fscanf(file, "%lf%s\n", &A_PPP::_settings._smooth, dum);
  fclose(file);
  
  A_PPP::_settings._map.loadATX(atxFileName);
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Parse processLowlevel args
static void ParseArgs(int argc, char *argv[], char * confFileName, char * roverFileName)
{
  if (argc < 5 || argc > 6)
  {
    fprintf (stderr, "Usage : processLowlevel\n"
             " - [E] < measurements.txt : lowlevel file\n"
             " - [E] -conf conf.txt : configuration file\n"
             " - [E] -rover rover.txt : rover file\n"
	     " - [E] [-verbose] : verbose file\n"
             " - [S] > output.sparse : output file\n");
    exit(1);
  }

  // Save args.
  for (int iArg=1; iArg<argc; iArg++)
  {
    if (strcmp(argv[iArg], "-conf") == 0)
    {
      iArg++;
      strcpy(confFileName, argv[iArg]);
    }
    else if (strcmp(argv[iArg], "-rover") == 0)
    {
      iArg++;
      strcpy(roverFileName, argv[iArg]);
    }
    else if (strcmp(argv[iArg], "-verbose") == 0)
    {
      verbose = true;
    }
    else
    {
      fprintf (stderr, "Usage : processLowlevel\n"
             " - [E] < measurements.txt : lowlevel file\n"
             " - [E] -conf conf.txt : configuration file\n"
             " - [E] -rover rover.txt : rover file\n"
	     " - [E] [-verbose] : verbose file\n"
             " - [S] > output.sparse : output file\n");
      exit(1);
    }
  } 
}

//////////////////////////////////////////////////////////////////////////////////////////////
//Main
int main(int argc, char *argv[])
{
  A_MEASUREMENT meas;
  A_BIAS bias;
  A_SSR_PARAMETER pos;
  int sat;
  int syst;
  int numRover;
  A_DATE date;
  double tropo;
  char confFileName[1024];
  char roverFileName[1024];
  char sLine[1024]; /* Read lines. */

  /* Read args and conf files */
  outputVerbose=0;
  verbose = false;
  ParseArgs(argc, argv,confFileName,roverFileName);
  readParam(confFileName);
  readRover(roverFileName);
   
  /* Process */
  while (1)
  {
    *sLine = 0;
    fgets (sLine, 1024, stdin);
    if (feof(stdin))
    {
      // End of File : if one date has been read, process of the measurements and stop
      if (!((long)date.getSecond() % 1L))
      {
     	A_DATE dateNext(date.getDay(),date.getSecond() + A_PPP::_settings._dt);
        processMeasurement(numRover, -1, -1, dateNext, meas, bias, pos, tropo);
      }
      break;
    }
    if (readMeasurement(sLine, &numRover, &syst, &sat, &date, &meas, &bias, &pos, &tropo))
    {
      if (!((long)date.getSecond() % 1L))
      {
        if ((date - listRover[numRover-1]._date) >= -0.001)
        {
          processMeasurement(numRover, syst, sat, date, meas, bias, pos, tropo);
	  date.setDay(0);
	  date.setSecond(0.0);
        }
      }
    }
    else
    {
      fprintf(stderr, "Error while reading measurement.\n");
      exit(1);
    }
    
  }

  return 0;
}
