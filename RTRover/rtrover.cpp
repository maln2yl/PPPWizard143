/************************************************************
Nom ......... : rtrover.cpp
Role ........ : main interface
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/

#include <cstdio>
#include <vector>
#include <cstring>

#include "rtrover.h"
#include "rtrover_rover.h"
#include "rtrover_ppp.h"
#include "rtrover_date.h"
#include "rtrover_broadcast.h"
#include "rtrover_vector3D.h"
#include "rtrover_utility.h"
#include "rtrover_frequency.h"

#define EPSDT (1.e-3)

// Array of ephemerides
static std::vector< std::vector<A_BROADCAST *> > tabEphPrec(SystMax);
static std::vector< std::vector<A_BROADCAST *> > tabEphLast(SystMax);

// Array of corrections
static std::vector<std::vector<rtrover_orbCorr> > tabOrbCorr(SystMax);
static std::vector<std::vector<rtrover_clkCorr> > tabClkCorr(SystMax);
static std::vector< std::vector<rtrover_sbasCorr> > tabSbasCorr(SystMax);

// Array of biases
static std::vector< std::vector<rtrover_satBiases> > tabSatBiases(SystMax);
static std::vector< std::vector<rtrover_satPhaseBiases> > tabSatPhaseBiases(SystMax);
static std::vector<std::vector<int> > tabSatBiasesPriority(SystMax);

// Global settings
static int correct, correctType, outputVerb, lowlevel;
static double maxAge;
A_PPP_SETTING A_PPP::_settings;
A_FILTER_SETTING A_FILTER::_filterSetting;

// Rover
static std::vector<A_ROVER> tabRover;
static int currentRover;

// Satellite position
static std::vector<std::vector<A_VECTOR3D> > tabSatPosition(SystMax);

//LOG
static FILE* log_rtrover;
static FILE* log_lowlevel;
static bool fileOpened;

///////////////////////////////////////////////////////////////////////////////////////////////////
// Initialization
static void initializeOrbCorr(rtrover_orbCorr& orbCorr)
{
  orbCorr._satellite._system='G';
  orbCorr._satellite._number=0;
  orbCorr._iod=0;
  orbCorr._time._mjd=0;
  orbCorr._time._sec=0.0;
  for(int i=0; i<3; i++)
  {
    orbCorr._rao[i]=0.0;
    orbCorr._dotRao[i]=0.0;
    orbCorr._dotDotRao[i]=0.0;
  }
}

static void initializeClkCorr(rtrover_clkCorr& clkCorr)
{
  clkCorr._satellite._system='G';
  clkCorr._satellite._number=0;
  clkCorr._iod=0;
  clkCorr._time._mjd=0;
  clkCorr._time._sec=0.0;
  clkCorr._dClk=0.0;
  clkCorr._dotDClk=0.0;
  clkCorr._dotDotDClk=0.0;
}

static void initializeSbasCorr(rtrover_sbasCorr& sbasCorr)
{
  sbasCorr._satellite._system='G';
  sbasCorr._satellite._number=0;
  sbasCorr._time._mjd=0;
  sbasCorr._time._sec=0.0;
  sbasCorr._IODE=0;
  for(int i=0; i<3; i++)
  {
    sbasCorr._orb[i]=0.0;
  }
  sbasCorr._clk=0.0;
}

static void initializeSatBiases(rtrover_satBiases& satBiases)
{
  satBiases._satellite._system='G';
  satBiases._satellite._number=0;
  satBiases._time._mjd=0;
  satBiases._time._sec=0.0;
  satBiases._numBiases=0;
  satBiases._biases=NULL;
}

static void initializePhaseBiases(rtrover_satPhaseBiases& phaseBiases)
{
  phaseBiases._satellite._system='G';
  phaseBiases._satellite._number=0;
  phaseBiases._time._mjd=0;
  phaseBiases._time._sec=0.0;
  phaseBiases._numBiases=0;
  phaseBiases._biases=NULL;
}

//////////////////////////////////////////////////////////////////////
// Convert
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
// TypeBds : Geo, IGSO or Meo
static int typeBds(const A_VECTOR3D xc, const A_VECTOR3D vv)
{
  int type = -1;
  if ((xc.getModule() < EPSDBL) || (vv.getModule() < EPSDBL))
    return type;

  static const double omegaEarth = 7292115.1467e-11;
  static const double gmWGS      = 398.6004418e12;
  A_VECTOR3D InVit, N;

  // Inertial velocity
  InVit.setX(vv.getX() - omegaEarth * xc.getY());
  InVit.setY(vv.getY() + omegaEarth * xc.getX());
  InVit.setZ(vv.getZ());
  N = xc ^ InVit;

  // Inclination
  const double inc = acos(N.getZ() / N.getModule());
  if (fabs(inc) < 0.1) {
    type = GEO;
  } else {

    const double v2 = InVit.getModule() * InVit.getModule();
    const double semiMajorAxis = 1.0 / ((2.0 / xc.getModule()) - (v2 / gmWGS));
    type = (semiMajorAxis >= 35.0e6) ? IGSO : MEO;
  } 

  return type;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Observations and biases
static void obsAndBiases(int isyst, int isat, rtrover_satObs *satObs, double codeFreq[FMax],
                         double phaseFreq[FMax], double dopplerFreq[FMax],double biaisCodeFreq[FMax],
			 double biaisPhaseFreq[FMax], bool n1i[FMax], int wli[FMax], int discontinuity[FMax], double *yawAngle)
{
  static char Pref[6]={'C','P','W','Q','X','B'};
  *yawAngle = 0.;
  for (int i=0;i<FMax;i++)
  {
    codeFreq[i]=0.0;
    phaseFreq[i]=0.0;
    dopplerFreq[i]=0.0;
    biaisCodeFreq[i]=0.0;
    biaisPhaseFreq[i]=0.0;
    n1i[i]=false;
    wli[i]=0;
    discontinuity[i] = 0;
  }
  
  rtrover_satBiases &satBiases = tabSatBiases[isyst][isat];
  rtrover_satPhaseBiases &satPhaseBias = tabSatPhaseBiases[isyst][isat];

  for (int k=5;k>=0;k--)
  {
    for (int f=0;f<FMax;f++)
    {
      for (int j=0;j<satObs->_numObs;j++)
      {
	char cf[FMax]={'1','2','5','6','7'};
	if ((Pref[k] == satObs->_obs[j]._rnxType2ch[1]) && (cf[f] == satObs->_obs[j]._rnxType2ch[0]))
	{
	  if (satObs->_obs[j]._codeValid)
		codeFreq[f] = satObs->_obs[j]._code;
	  if (satObs->_obs[j]._phaseValid)
		phaseFreq[f] = satObs->_obs[j]._phase;
	  if (satObs->_obs[j]._dopplerValid)
		dopplerFreq[f] = satObs->_obs[j]._doppler;

          biaisCodeFreq[f] = 0.0;
	  biaisPhaseFreq[f] = 0.0;
	  // PATCH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  if (isyst == 2)
	  {
	    if (satBiases._satellite._number)
	    {
              for (int l=0;l<satBiases._numBiases;l++)
	      {
		if ((satBiases._biases[l]._rnxType3ch[0] == 'C') &&
		  (satBiases._biases[l]._rnxType3ch[1] == satObs->_obs[j]._rnxType2ch[0]) /*&&
		  (satBiases._biases[l]._rnxType3ch[2] == satObs->_obs[j]._rnxType2ch[1])*/)
        	{
        	  biaisCodeFreq[f] = satBiases._biases[l]._value;
        	}
	      }
            }
	  }
	  // PATCH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  else
	  {
	    if (satBiases._satellite._number)
	    {
              for (int l=0;l<satBiases._numBiases;l++)
	      {
		if ((satBiases._biases[l]._rnxType3ch[0] == 'C') &&
		  (satBiases._biases[l]._rnxType3ch[1] == satObs->_obs[j]._rnxType2ch[0]) &&
		  (satBiases._biases[l]._rnxType3ch[2] == satObs->_obs[j]._rnxType2ch[1]))
        	{
        	  biaisCodeFreq[f] = satBiases._biases[l]._value;
        	}
	      }
            }
	  }
	  if (satPhaseBias._satellite._number)
          {
            *yawAngle = satPhaseBias._yawAngleValue;
	    for (int l=0;l<satPhaseBias._numBiases;l++)
            {
              if ((satPhaseBias._biases[l]._rnxType3ch[0] == 'C') &&
                  (satPhaseBias._biases[l]._rnxType3ch[1] == satObs->_obs[j]._rnxType2ch[0]) /*&&
                  (satPhaseBias._biases[l]._rnxType3ch[2] == 'Z')*/)
              {
		biaisPhaseFreq[f] = satPhaseBias._biases[l]._value;
		n1i[f] = satPhaseBias._biases[l]._intIndicator;
		wli[f] = satPhaseBias._biases[l]._wlIntIndicator;
		discontinuity[f] = satPhaseBias._biases[l]._discontinuityValue;
              }
            }
          }
	}
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////
void RSW_to_XYZ(A_VECTOR3D& rr, A_VECTOR3D& vv, A_VECTOR3D& rsw, A_VECTOR3D& xyz)
{
  A_VECTOR3D along=vv/vv.getModule();
  
  A_VECTOR3D cross=rr^vv;
  cross.doubleEtModule(1.0);
    
  A_VECTOR3D radial=along^cross;

  xyz=radial*rsw.getX()+along*rsw.getY()+cross*rsw.getZ();
}

////////////////////////////////////////////////////////////////////////////
// Apply ephemerides and clock corrections
bool applyCorr(const rtrover_orbCorr *orbCorr, rtrover_clkCorr *clkCorr, rtrover_time *tt, A_VECTOR3D& xc, double& tc, A_VECTOR3D& vv)
{
  double dtRao = diffDate(tt, &orbCorr->_time);

  // Position
  // --------
  
  A_VECTOR3D rao(orbCorr->_rao[0],orbCorr->_rao[1],orbCorr->_rao[2]);
  A_VECTOR3D dotRao(orbCorr->_dotRao[0], orbCorr->_dotRao[1], orbCorr->_dotRao[2]);
  A_VECTOR3D dotDotRao(orbCorr->_dotDotRao[0], orbCorr->_dotDotRao[1], orbCorr->_dotDotRao[2]);
  
  A_VECTOR3D raoHlp = rao + dotRao * dtRao + dotDotRao * 0.5 * dtRao * dtRao;

  if (raoHlp.getModule() > 20.0)
    return false;

  A_VECTOR3D dx;
  RSW_to_XYZ(xc, vv, raoHlp, dx);
  xc -= dx;

  // Velocity
  // --------
  A_VECTOR3D dotRaoHlp = dotRao + dotDotRao * dtRao ;

  A_VECTOR3D dv;
  RSW_to_XYZ(xc, vv, dotRaoHlp, dv);
  vv -= dv;

  // Clocks
  // ------
  double dtClk = diffDate(tt, &clkCorr->_time);
  //fprintf(stderr, "Clock %15.3lf %15.3lf\n", dtClk, CLIGHT*(clkCorr->_dClk + clkCorr->_dotDClk * dtClk + clkCorr->_dotDotDClk * dtClk * dtClk));

  tc += clkCorr->_dClk + clkCorr->_dotDClk * dtClk + clkCorr->_dotDotDClk * dtClk * dtClk;

  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Interpolated and corrected position
bool interpolatedCorrectedPosition(const rtrover_time *tt, const int isyst, const int isat, const double range, A_VECTOR3D& xc, double &tc, A_VECTOR3D& vv)
{
  rtrover_sbasCorr &sbasCorr = tabSbasCorr[isyst][isat];
  rtrover_orbCorr &orbCorr = tabOrbCorr[isyst][isat];
  rtrover_clkCorr &clkCorr = tabClkCorr[isyst][isat];
  A_BROADCAST *ephLast = tabEphLast[isyst][isat];
  A_BROADCAST *ephPrec = tabEphPrec[isyst][isat];

  if (correct)
  {
    if (correctType == SBAS)
    {
      if (!sbasCorr._satellite._number) return false;
      if (!sbasCorr._time._mjd && !sbasCorr._time._sec) return false;
      if (diffDate(tt, &sbasCorr._time) > 60.0) return false;
    }
    else
    {
      if (!orbCorr._satellite._number) return false;
      if (!clkCorr._satellite._number) return false;
      if (!orbCorr._time._mjd && !orbCorr._time._sec) return false;
      if (!clkCorr._time._mjd && !clkCorr._time._sec) return false;
      if (!clkCorr._iod && !orbCorr._iod) return false;
      if (diffDate(tt, &clkCorr._time) > 60.0) return false;
    }
  }


  if (tabEphLast[isyst][isat]!=NULL)
  {
    A_BROADCAST *eph=NULL;
    if (correct)
    {
      if (correctType == SBAS)
      {
	if (ephLast && (ephLast->getIODE() == sbasCorr._IODE))
          eph = ephLast;
        else if (ephPrec &&(ephPrec->getIODE() == sbasCorr._IODE))
          eph = ephPrec;
      }
      else
      {
	int clk_iod = 0;
	clk_iod = clkCorr._iod; 
	if (ephLast && (ephLast->getIODC() == clk_iod))
          eph = ephLast;
        else if (ephPrec && (ephPrec->getIODC() == clk_iod))
          eph = ephPrec;
      }
    }
    else
      eph = ephLast;

    if (eph)
    {
      rtrover_time tt0=*tt;
      tt0._sec -= range/CLIGHT;
      eph->computePosition(&tt0, xc, tc, vv, correct);
      if (!tc)
        return false;
      tt0._sec -= tc;
      eph->computePosition(&tt0, xc, tc, vv, correct);
      bool b=true;
      if (correct)
      {
        if (correctType == SBAS)
        {
          xc.setX(xc.getX() + sbasCorr._orb[0]);
          xc.setY(xc.getY() + sbasCorr._orb[1]);
          xc.setZ(xc.getZ() + sbasCorr._orb[2]);
          tc += sbasCorr._clk;
        }
        else
          b=applyCorr(&orbCorr, &clkCorr, &tt0, xc, tc, vv);
      } 
      
      switch (isyst)
      {
        case GPS :
	  PRINT_LOG(log_rtrover,"EPH GPS: %02d %d %15.3lf %15.3lf %15.3lf %15.3lf\n",isat+1, eph->getIODC(), xc.getX(), xc.getY(), xc.getZ(), CLIGHT*tc);
	  break;
        case Glo :
          PRINT_LOG(log_rtrover,"EPH GLO: %02d %15.3lf %15.3lf %15.3lf %15.3lf\n",isat+1, xc.getX(), xc.getY(), xc.getZ(), CLIGHT*tc);
          break;
        case Gal :
          PRINT_LOG(log_rtrover,"EPH GAL: %02d %15.3lf %15.3lf %15.3lf %15.3lf\n", isat+1, xc.getX(), xc.getY(), xc.getZ(), CLIGHT*tc);	 
          break;
	case Bds :
	  PRINT_LOG(log_rtrover,"EPH BDS: %02d %15.3lf %15.3lf %15.3lf %15.3lf\n", isat+1, xc.getX(), xc.getY(), xc.getZ(), CLIGHT*tc);	 
	  break;
        default :
          PRINT_LOG(log_rtrover,"EPH error !\n");	
      }
      
      return b;
    }
  }

  return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Epoch
static void processEpoch(std::vector<std::vector<rtrover_satObs> >& satObs, char *str_output)
{
  rtrover_time time;
  A_VECTOR3D xc, vv;
  double tc = 0.0;

 
  std::vector<std::vector<A_MEASUREMENT> > meas(SystMax);
  std::vector<std::vector<A_BIAS> > bias(SystMax);
  std::vector<std::vector<A_SSR_PARAMETER> > pos(SystMax);

  A_DATE date;
  time._mjd = 0;
  time._sec = 0.0;
  
  for (int isyst=0;isyst<SystMax;isyst++)
  {
    const int maxSatSyst = getMaxSatSyst(isyst);
    meas[isyst] = std::vector<A_MEASUREMENT>(maxSatSyst);
    bias[isyst] = std::vector<A_BIAS>(maxSatSyst);
    pos[isyst] = std::vector<A_SSR_PARAMETER>(maxSatSyst);
  }
  
  // Date
  for (int isyst=0;isyst<SystMax;isyst++)
  { 
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      if (satObs[isyst][isat]._satellite._number)
        time = satObs[isyst][isat]._time;
    }
  }
  if ((time._mjd == 0) && (time._sec == 0.0))
    return;
    
  PRINT_LOG(log_rtrover,"processEpoch: %d %lf\n", time._mjd, time._sec);

  // Time correction
  {
   
  if (tabRover[currentRover]._lastTime._mjd && tabRover[currentRover]._lastTime._sec)
  {
   
    for (double js=(tabRover[currentRover]._lastTime._mjd*86400.0+tabRover[currentRover]._lastTime._sec+A_PPP::_settings._dt);js<(time._mjd*86400.0+time._sec);js+=A_PPP::_settings._dt)
    {
      date.setDay((int) floor(js/86400)-33282.0);    
      date.setSecond(js-date.getDay()*86400.0);
      tabRover[currentRover]._ppp.computePPP(meas, bias, pos, 0.0, date, tabRover[currentRover]._stationProduct, tabRover[currentRover]._stationMeasurement, tabRover[currentRover]._roverAPrioriPosition, currentRover);
      tabRover[currentRover]._lastTime._mjd = (int) floor(js/86400);
      tabRover[currentRover]._lastTime._sec = js-tabRover[currentRover]._lastTime._mjd*86400.0;    		
    }
   }
  }

  // Observations and related biaises
  for (int isyst=0; isyst<SystMax;isyst++)
  {
    std::vector<rtrover_satObs> &satObsConst = satObs[isyst];
    std::vector<A_MEASUREMENT> &measConst = meas[isyst];
    std::vector<A_BIAS> &biasConst = bias[isyst];
    std::vector<A_VECTOR3D> &tabSatPositionConst = tabSatPosition[isyst];
    std::vector<rtrover_ionoCorr> &tabIonoCorrConst = tabRover[currentRover]._tabIonoCorr[isyst];
    std::vector<A_SSR_PARAMETER> &posConst = pos[isyst];

    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      rtrover_satObs &obsSat = satObsConst[isat];
      A_MEASUREMENT &measSat = measConst[isat];
      A_BIAS &biasSat = biasConst[isat];
      A_SSR_PARAMETER &posSat = posConst[isat];
      A_VECTOR3D &satPosition = tabSatPositionConst[isat];
      rtrover_ionoCorr &ionoCorr = tabIonoCorrConst[isat];
      if (obsSat._satellite._number)
      {
        double yawAngle = 0.;
	double codeFreq[FMax], phaseFreq[FMax], dopplerFreq[FMax], biaisCodeFreq[FMax], biaisPhaseFreq[FMax];
	bool n1i[FMax];
	int wli[FMax];
	int discontinuity[FMax];
	obsAndBiases(isyst,isat, &obsSat, codeFreq, phaseFreq, dopplerFreq, biaisCodeFreq, biaisPhaseFreq, n1i, wli, discontinuity, &yawAngle);
	if (interpolatedCorrectedPosition(&obsSat._time, isyst, isat, codeFreq[0], xc, tc, vv))
	{
	  measSat._code = std::vector<double>(codeFreq,codeFreq+FMax);
	  measSat._phase = std::vector<double>(phaseFreq,phaseFreq+FMax);
	  measSat._doppler = std::vector<double>(dopplerFreq,dopplerFreq+FMax);
	  measSat._slot_typeBds = obsSat._slotNumber;
	  if (isyst == Bds)
	    measSat._slot_typeBds = typeBds(xc, vv);
	  posSat._Xsp3 = xc.getX(); 
	  posSat._Ysp3 = xc.getY();
	  posSat._Zsp3 = xc.getZ();
	  posSat._Hsp3 = CLIGHT*tc;
	  satPosition = xc;
	  posSat._iono = ionoCorr._value;
	  posSat._typeIono = ionoCorr._flag;
	  posSat._yaw = yawAngle;
	  biasSat._code = std::vector<double>(biaisCodeFreq,biaisCodeFreq+FMax);
	  biasSat._phase = std::vector<double>(biaisPhaseFreq,biaisPhaseFreq+FMax);
	  posSat._VXsp3 = vv.getX(); 
	  posSat._VYsp3 = vv.getY();
	  posSat._VZsp3 = vv.getZ();

	  for (int f=0;f<FMax;f++)
	  {
	    if(wli[f] > 0)
	      biasSat._wli=1;
	    biasSat._n1i = biasSat._n1i | n1i[f];
	  }  
	  biasSat._discontinuity = discontinuity[0];

          if (lowlevel)
	  {
	    if (!fileOpened)
	    {
              log_lowlevel = fopen("lowlevel.txt", "a");
	      fileOpened = true;
            }
	    
	    fprintf(log_lowlevel,"%d %c%02d %d %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %d %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %15.3lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d\n", 
                currentRover+1, reverseSystem(isyst), isat+1, time._mjd-33282, time._sec,
		measSat._code[F1], measSat._code[F2], measSat._code[F6], measSat._code[F5], measSat._code[F7],
		measSat._phase[F1], measSat._phase[F2], measSat._phase[F6], measSat._phase[F5], measSat._phase[F7],
		measSat._doppler[F1], measSat._doppler[F2], measSat._doppler[F6], measSat._doppler[F5], measSat._doppler[F7], measSat._slot_typeBds,
		posSat._Xsp3, posSat._Ysp3, posSat._Zsp3, posSat._Hsp3, posSat._VXsp3, posSat._VYsp3, posSat._VZsp3, posSat._iono, posSat._typeIono, tabRover[currentRover]._tropoCorr._value, posSat._yaw, 
                biasSat._code[F1], biasSat._code[F2], biasSat._code[F6], biasSat._code[F5], biasSat._code[F7],
		biasSat._phase[F1], biasSat._phase[F2], biasSat._phase[F6], biasSat._phase[F5], biasSat._phase[F7],
		biasSat._n1i?1:0, biasSat._wli, biasSat._discontinuity);
         
          }
	}
      }
    }
  }
  
  double fpas = (time._mjd*86400.0+time._sec)/A_PPP::_settings._dt;
  double fracf = fabs(fpas - floor(fpas+0.5));
  
  if (( (tabRover[currentRover]._lastTime._mjd == 0) || (fabs(time._mjd*86400.0+time._sec - tabRover[currentRover]._lastTime._mjd*86400.0 - tabRover[currentRover]._lastTime._sec - A_PPP::_settings._dt) < EPSDT)) && (fracf < EPSDBL))
  {	  
  	tabRover[currentRover]._lastTime = time;
  	date.setDay(time._mjd-33282);
 	date.setSecond(time._sec);
  	if (tabRover[currentRover]._ppp.computePPP(meas, bias, pos, tabRover[currentRover]._tropoCorr._value, date,
      	    tabRover[currentRover]._stationProduct, tabRover[currentRover]._stationMeasurement, tabRover[currentRover]._roverAPrioriPosition, currentRover))
  	{
    	    char str[200];
            int day, month, year, hour, minute;
    	    double second;
    	    date.calendar(&day, &month, &year, &hour, &minute, &second);
    	    if (outputVerb)
    	    {
      		sprintf(str, "%02d-%02d-%02d %02d:%02d:%06.3lf %6.3f %6.3f %6.3f %6.3f %s PPP %2d %2d %2d %2d %12.3f +- %.3f %12.3f +- %.3f %12.3f +- %.3f %5.3f + %5.6f +- %5.6f  %.4f\n",
	    	   year%100, month, day, hour, minute, second,
	           tabRover[currentRover]._stationProduct._HstaGPS, tabRover[currentRover]._stationProduct._HstaGlo,
	    	   tabRover[currentRover]._stationProduct._HstaGal, tabRover[currentRover]._stationProduct._HstaBds, tabRover[currentRover]._roverName.c_str(),
	           tabRover[currentRover]._stationProduct._nbMes, tabRover[currentRover]._stationProduct._nbMesBloEx, tabRover[currentRover]._stationProduct._nbMesBloNw, tabRover[currentRover]._stationProduct._nbMesBloN1,
	           tabRover[currentRover]._stationProduct._x, tabRover[currentRover]._stationProduct._covX,
	           tabRover[currentRover]._stationProduct._y, tabRover[currentRover]._stationProduct._covY,
	           tabRover[currentRover]._stationProduct._z, tabRover[currentRover]._stationProduct._covZ,
	           tabRover[currentRover]._stationProduct._tropoV, tabRover[currentRover]._stationProduct._tropo, tabRover[currentRover]._stationProduct._covTropo, tabRover[currentRover]._stationProduct._integrity);
    	    }
    	    else
            {
      		sprintf(str, "%02d-%02d-%02d %02d:%02d:%06.3lf %s PPP %2d %2d %2d %2d %12.3f +- %.3f %12.3f +- %.3f %12.3f +- %.3f %5.3f + %5.6f +- %5.6f  %.4f\n",
		   year%100, month, day, hour, minute, second,
	           tabRover[currentRover]._roverName.c_str(),
	           tabRover[currentRover]._stationProduct._nbMes, tabRover[currentRover]._stationProduct._nbMesBloEx, tabRover[currentRover]._stationProduct._nbMesBloNw, tabRover[currentRover]._stationProduct._nbMesBloN1,
	           tabRover[currentRover]._stationProduct._x, tabRover[currentRover]._stationProduct._covX,
	           tabRover[currentRover]._stationProduct._y, tabRover[currentRover]._stationProduct._covY,
	           tabRover[currentRover]._stationProduct._z, tabRover[currentRover]._stationProduct._covZ,
	           tabRover[currentRover]._stationProduct._tropoV, tabRover[currentRover]._stationProduct._tropo, tabRover[currentRover]._stationProduct._covTropo, tabRover[currentRover]._stationProduct._integrity);
    	    }
    	strcat(str_output, str);
  	}
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Interface

RTROVER_API void rtrover_addRover(std::string roverName, double aprPos[3])
{
  A_ROVER rover(roverName,aprPos);
  tabRover.push_back(rover);
}

RTROVER_API void rtrover_setCurrentRover(const int n)
{
  currentRover = n;
}

RTROVER_API void rtrover_initObs(rtrover_obs* obs)
{
  PRINT_LOG(log_rtrover,"rtrover_initObs\n");
}

RTROVER_API void rtrover_initOptions(rtrover_opt* opt)
{
  log_rtrover = OPEN_LOG("rtrover.txt");
  fileOpened = false;
  
  PRINT_LOG(log_rtrover,"rtrover_initOptions\n");
  
  for (int isyst=0; isyst < SystMax ; isyst++)
  {
    const int maxSatSyst = getMaxSatSyst(isyst);
    tabEphPrec[isyst] = std::vector<A_BROADCAST *>(maxSatSyst);
    tabEphLast[isyst] = std::vector<A_BROADCAST *>(maxSatSyst);
    tabOrbCorr[isyst] = std::vector<rtrover_orbCorr>(maxSatSyst);
    tabClkCorr[isyst] = std::vector<rtrover_clkCorr>(maxSatSyst);
    tabSbasCorr[isyst] = std::vector<rtrover_sbasCorr>(maxSatSyst);
    tabSatBiases[isyst] = std::vector<rtrover_satBiases>(maxSatSyst);
    tabSatPhaseBiases[isyst] = std::vector<rtrover_satPhaseBiases>(maxSatSyst);
    tabSatBiasesPriority[isyst] = std::vector<int>(maxSatSyst,0);
    tabSatPosition[isyst] = std::vector<A_VECTOR3D>(maxSatSyst);
    
    for (int isat=0; isat < maxSatSyst ; isat++)
    {
      initializeOrbCorr(tabOrbCorr[isyst][isat]);
      initializeClkCorr(tabClkCorr[isyst][isat]);
      initializeSbasCorr(tabSbasCorr[isyst][isat]);
      initializeSatBiases(tabSatBiases[isyst][isat]);
      initializePhaseBiases(tabSatPhaseBiases[isyst][isat]);
    }

  }

  currentRover = 0;
  
}

RTROVER_API void rtrover_setOptions(const rtrover_opt* opt)
{
  PRINT_LOG(log_rtrover,"rtrover_setOptions\n");

  A_MODE &mode = A_PPP::_settings._mode;
  
  switch (opt->_mode)
  {
    case mode_PPP_DF:  ///< dual-frequency precise point positioning
      mode._code = 1;
      mode._phase = 1;
      mode._freq[F1] = 1;
      mode._freq[F2] = 1;
      mode._freq[F5] = 1;
      mode._freq[F6] = 1;
      mode._freq[F7] = 1;
      mode._ambig = 0;
      correct = 1;
      break;
    case mode_SPP_DF:  ///< dual-frequency single point positioning
      mode._code = 1;
      mode._phase = 1;
      mode._freq[F1] = 1;
      mode._freq[F2] = 1;
      mode._freq[F5] = 1;
      mode._freq[F6] = 1;
      mode._freq[F7] = 1;
      mode._ambig = 0;
      correct = 0;
      break;
    case mode_PPP_SF:  ///< single-frequency precise precise point positioning
      mode._code = 1;
      mode._phase = 1;
      mode._freq[F1] = 1;
      mode._freq[F2] = 0;
      mode._freq[F5] = 0;
      mode._freq[F6] = 0;
      mode._freq[F7] = 0;
      mode._ambig = 0;
      correct = 1;
      break;
    case mode_SPP_SF:  ///< single-frequency single point positioning
      mode._code = 1;
      mode._phase = 1;
      mode._freq[F1] = 1;
      mode._freq[F2] = 0;
      mode._freq[F5] = 0;
      mode._freq[F6] = 0;
      mode._freq[F7] = 0;
      mode._ambig = 0;
      correct = 0;
      break;
    case mode_PPP_AR:  ///< precise point positining with ambiguity resolution
      mode._code = 1;
      mode._phase = 1;
      mode._freq[F1] = 1;
      mode._freq[F2] = 1;
      mode._freq[F5] = 1;
      mode._freq[F6] = 1;
      mode._freq[F7] = 1;
      mode._ambig = 1;
      correct = 1;
      break;
    case mode_RTK:     ///< real-time kinematics
      mode._code = 0;
      mode._phase = 0;
      mode._freq[F1] = 0;
      mode._freq[F2] = 0;
      mode._freq[F5] = 0;
      mode._freq[F6] = 0;
      mode._freq[F7] = 0;
      mode._ambig = 0;
      correct = 0;
      break;
    case mode_PPP_FTTF: ///< fast time-to-fix ambiguity resolution PPP
      mode._code = 0;
      mode._phase = 0;
      mode._freq[F1] = 0;
      mode._freq[F2] = 0;
      mode._freq[F5] = 0;
      mode._freq[F6] = 0;
      mode._freq[F7] = 0;
      mode._ambig = 0;
      correct = 0;
      break;
  }

  A_FILTER::_filterSetting._glo._use = opt->_useGlonass;
  
  // When RTRover is called by BNC
  if (strlen(opt->_roverName) && (tabRover.size() == 0))
  {
    double aprPos[3];
    aprPos[0]=0.0;
    aprPos[1]=0.0;
    aprPos[2]=0.0;
    if ((opt->_xyzAprRover[0] != 0.0) || (opt->_xyzAprRover[1] != 0.0) || (opt->_xyzAprRover[2] != 0.0))
    {
      aprPos[0]=opt->_xyzAprRover[0];
      aprPos[1]=opt->_xyzAprRover[1];
      aprPos[2]=opt->_xyzAprRover[2];
    }
    rtrover_addRover(opt->_roverName,aprPos);
  }

  A_PPP::_settings._map.loadATX((char *)opt->_antexFileName);

  rtrover_additionalOpt addOpt;
  addOpt._dt = 1.0;
  addOpt._maxAge = 10.0;
  addOpt._maxElim = 3;
  addOpt._raim = 1;
  addOpt._minFix = 3600;
  addOpt._thrMap = 6.0;
  addOpt._sigIniTro = 0.5;
  addOpt._sigModTro = 0.000005;
  addOpt._nbSatFixAmb = 4;
  addOpt._thrAmb = 0.1;
  addOpt._sigIniClkBias = 0.0;
  addOpt._sigModClkBias = 0.001;
  addOpt._sigIniIono = 10.0;
  addOpt._sigModIono = 0.002;
  addOpt._sigIniPos = 10.0;
  addOpt._sigModPos = 10.0;
  addOpt._sigMeasIono[0]= 0.01;
  addOpt._sigMeasIono[1]= 0.2;
  addOpt._sigMeasIono[2]= 1.0;
  addOpt._thrMeasIono = 5.0;
  addOpt._sigMeasTropo = 0.1;
  addOpt._thrMeasTropo = 1.0;
  addOpt._dtMax = 300;
  addOpt._thrMeasCode = 10.0;
  addOpt._thrMeasPhase = 0.05;
  addOpt._sigMeasCodeGps = 1.0;
  addOpt._sigMeasPhaseGps = 0.01;
  addOpt._sigMeasCodeGlo = 5.0;
  addOpt._sigMeasPhaseGlo = 0.01;
  addOpt._sigMeasCodeGal = 5.0;
  addOpt._sigMeasPhaseGal = 0.01;
  addOpt._sigMeasCodeBds[0] = 5.0;
  addOpt._sigMeasPhaseBds[0] = 0.01;
  addOpt._sigMeasCodeBds[1] = 5.0;
  addOpt._sigMeasPhaseBds[1] = 0.01;
  addOpt._sigMeasCodeBds[2] = 5.0;
  addOpt._sigMeasPhaseBds[2] = 0.01;
  addOpt._smooth = 0.0;
  addOpt._correction = 0;
  addOpt._useGPS = true;
  addOpt._useGalileo = false;
  addOpt._useBeidou = false;
  addOpt._reset = 0;
  addOpt._outputVerbose = 0;
  addOpt._fixN1 = 1;
  addOpt._fixNw = 1;
  addOpt._fixEx = 1;
  addOpt._lowlevel = 0;
  rtrover_setAdditionalOpt(&addOpt);

  int use;
  
  char dum[50];
  FILE *fichier_rtrover=fopen("RTRover.ini","r");
  if (fichier_rtrover != NULL)
  {
    fscanf(fichier_rtrover, "%d%d%d%s\n", &addOpt._fixN1, &addOpt._fixNw, &addOpt._fixEx, dum);
    fscanf(fichier_rtrover, "%d%s\n", &use, dum); //useGps
    addOpt._useGPS = use ? true : false;
    fscanf(fichier_rtrover, "%d%s\n", &A_FILTER::_filterSetting._glo._use, dum);//useGlonass
    fscanf(fichier_rtrover, "%d%s\n", &use, dum); //useGalileo
    addOpt._useGalileo = use ? true : false;
    fscanf(fichier_rtrover, "%d%s\n", &use, dum); //useBeidou
    addOpt._useBeidou= use ? true : false;
    fscanf(fichier_rtrover, "%d%s\n", &addOpt._correction, dum);
    fscanf(fichier_rtrover, "%d%s\n", &addOpt._reset, dum);
    fscanf(fichier_rtrover, "%d%s\n", &addOpt._outputVerbose, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._dt, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._maxAge, dum);
    fscanf(fichier_rtrover, "%d%s\n", &addOpt._minFix, dum);
    fscanf(fichier_rtrover, "%d%s\n", &addOpt._maxElim, dum);
    fscanf(fichier_rtrover, "%d%s\n", &addOpt._raim, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._thrMap, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigIniTro, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigModTro, dum);
    fscanf(fichier_rtrover, "%d%s\n", &addOpt._nbSatFixAmb, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._thrAmb, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigIniClkBias, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigModClkBias, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigIniIono, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigModIono, dum);
    fscanf(fichier_rtrover, "%lf%lf%lf%s\n", &addOpt._sigMeasIono[0], &addOpt._sigMeasIono[1], &addOpt._sigMeasIono[2], dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._thrMeasIono, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigMeasTropo, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._thrMeasTropo, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigIniPos, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigModPos, dum);
    fscanf(fichier_rtrover, "%d%s\n", &addOpt._dtMax, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._thrMeasCode, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._thrMeasPhase, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigMeasCodeGps, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigMeasPhaseGps, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigMeasCodeGlo, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigMeasPhaseGlo, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigMeasCodeGal, dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._sigMeasPhaseGal, dum);
    fscanf(fichier_rtrover, "%lf%lf%lf%s\n", &addOpt._sigMeasCodeBds[0], &addOpt._sigMeasCodeBds[1], &addOpt._sigMeasCodeBds[2], dum);
    fscanf(fichier_rtrover, "%lf%lf%lf%s\n", &addOpt._sigMeasPhaseBds[0], &addOpt._sigMeasPhaseBds[1], &addOpt._sigMeasPhaseBds[2], dum);
    fscanf(fichier_rtrover, "%lf%s\n", &addOpt._smooth, dum);
    fclose(fichier_rtrover);
    rtrover_setAdditionalOpt(&addOpt);
  }
}

RTROVER_API void rtrover_setAdditionalOpt(const rtrover_additionalOpt* opt)
{
  PRINT_LOG(log_rtrover,"rtrover_setAdditionalOpt\n");
  
  for (int i = 0; i<3 ; i++)
  {
    A_FILTER::_filterSetting._gps._sigCode[i] = 0.0;
    A_FILTER::_filterSetting._gps._sigPhase[i] = 0.0;
    A_FILTER::_filterSetting._glo._sigCode[i] = 0.0;
    A_FILTER::_filterSetting._glo._sigPhase[i] = 0.0;
    A_FILTER::_filterSetting._gal._sigCode[i] = 0.0;
    A_FILTER::_filterSetting._gal._sigPhase[i] = 0.0;
    A_FILTER::_filterSetting._bds._sigCode[i] = 0.0;
    A_FILTER::_filterSetting._bds._sigPhase[i] = 0.0;
  }
 
  correctType = opt->_correction;
  maxAge = opt->_maxAge;
  if (maxAge > MAX_EPOCH*opt->_dt)
    maxAge = MAX_EPOCH*opt->_dt;
  outputVerb = opt->_outputVerbose;
  A_PPP::_settings._dt = opt->_dt;
  A_PPP::_settings._thrMap = opt->_thrMap;
  A_PPP::_settings._smooth = opt->_smooth;
  A_PPP::_settings._reset = opt->_reset;
//  A_PPP::_settings._triAR = opt->_triAR;
  A_FILTER::_filterSetting._gps._use = opt->_useGPS;
  A_FILTER::_filterSetting._gal._use = opt->_useGalileo;
  A_FILTER::_filterSetting._bds._use = opt->_useBeidou;
  A_FILTER::_filterSetting._nbMin = opt->_minFix;
  A_FILTER::_filterSetting._maxElim = opt->_maxElim;
  A_FILTER::_filterSetting._raim = opt->_raim;
  A_FILTER::_filterSetting._sigIniTro = opt->_sigIniTro * opt->_sigIniTro;
  A_FILTER::_filterSetting._sigModTro = opt->_sigModTro * opt->_sigModTro;
  A_FILTER::_filterSetting._nbSatFixAmb = (opt->_nbSatFixAmb >= 0) ? opt->_nbSatFixAmb : 0;
  A_FILTER::_filterSetting._thrAmb = opt->_thrAmb;
  A_FILTER::_filterSetting._sigIniBiasClk = opt->_sigIniClkBias * opt->_sigIniClkBias;
  A_FILTER::_filterSetting._sigModBiasClk = opt->_sigModClkBias * opt->_sigModClkBias;
  A_FILTER::_filterSetting._sigIniIono = opt->_sigIniIono * opt->_sigIniIono;
  A_FILTER::_filterSetting._sigModIono = opt->_sigModIono * opt->_sigModIono;
  A_FILTER::_filterSetting._sigIniPos = opt->_sigIniPos * opt->_sigIniPos;
  A_FILTER::_filterSetting._sigModPos = opt->_sigModPos * opt->_sigModPos;
  A_FILTER::_filterSetting._sigMeasIono[0] = opt->_sigMeasIono[0] * opt->_sigMeasIono[0];
  A_FILTER::_filterSetting._sigMeasIono[1] = opt->_sigMeasIono[1] * opt->_sigMeasIono[1];
  A_FILTER::_filterSetting._sigMeasIono[2] = opt->_sigMeasIono[2] * opt->_sigMeasIono[2];
  A_FILTER::_filterSetting._thrMeasIono = opt->_thrMeasIono;
  A_FILTER::_filterSetting._sigMeasTropo = opt->_sigMeasTropo * opt->_sigMeasTropo;
  A_FILTER::_filterSetting._thrMeasTropo = opt->_thrMeasTropo;
  A_FILTER::_filterSetting._thrMeasCode = opt->_thrMeasCode;
  A_FILTER::_filterSetting._thrMeasPhase = opt->_thrMeasPhase; 
  A_FILTER::_filterSetting._gps._sigCode[0] = opt->_sigMeasCodeGps * opt->_sigMeasCodeGps;
  A_FILTER::_filterSetting._gps._sigPhase[0] = opt->_sigMeasPhaseGps * opt->_sigMeasPhaseGps;
  A_FILTER::_filterSetting._glo._sigCode[0] = opt->_sigMeasCodeGlo * opt->_sigMeasCodeGlo;
  A_FILTER::_filterSetting._glo._sigPhase[0]  = opt->_sigMeasPhaseGlo * opt->_sigMeasPhaseGlo;
  A_FILTER::_filterSetting._gal._sigCode[0] = opt->_sigMeasCodeGal * opt->_sigMeasCodeGal;
  A_FILTER::_filterSetting._gal._sigPhase[0]  = opt->_sigMeasPhaseGal * opt->_sigMeasPhaseGal;
  A_FILTER::_filterSetting._bds._sigCode[0] = opt->_sigMeasCodeBds[0] * opt->_sigMeasCodeBds[0];
  A_FILTER::_filterSetting._bds._sigPhase[0]  = opt->_sigMeasPhaseBds[0] * opt->_sigMeasPhaseBds[0];
  A_FILTER::_filterSetting._bds._sigCode[1] = opt->_sigMeasCodeBds[1] * opt->_sigMeasCodeBds[1];
  A_FILTER::_filterSetting._bds._sigPhase[1]  = opt->_sigMeasPhaseBds[1] * opt->_sigMeasPhaseBds[1];
  A_FILTER::_filterSetting._bds._sigCode[2] = opt->_sigMeasCodeBds[2] * opt->_sigMeasCodeBds[2];
  A_FILTER::_filterSetting._bds._sigPhase[2]  = opt->_sigMeasPhaseBds[2] * opt->_sigMeasPhaseBds[2];
  A_FILTER::_filterSetting._dtMax = opt->_dtMax;  
//  A_FILTER::_filterSetting._fix = A_PPP::_settings._mode._ambig;
  A_FILTER::_filterSetting._fixN1 = opt->_fixN1;
  A_FILTER::_filterSetting._fixNw = opt->_fixNw;
  A_FILTER::_filterSetting._fixEx = opt->_fixEx;
  lowlevel= opt->_lowlevel;
}

RTROVER_API void rtrover_putGPSEphemeris(const rtrover_ephGPS* eph)
{
  PRINT_LOG(log_rtrover,"rtrover_putGPSEphemeris %d %lf %c%02d %d\n", eph->_TOC._mjd, eph->_TOC._sec, eph->_satellite._system, eph->_satellite._number, eph->_IODE);

  int isat=eph->_satellite._number-1;
  if ((isat != -1) && (isat < getMaxSatSyst(GPS)))
  {
    std::vector<A_BROADCAST*> &tabEphLastGPS = tabEphLast[GPS];
    std::vector<A_BROADCAST*> &tabEphPrecGPS = tabEphPrec[GPS];
    rtrover_satBiases &satBiasesGPS = tabSatBiases[GPS][isat];
    int &satBiasesPriority = tabSatBiasesPriority[GPS][isat];

    A_BROADCAST_GPS* ephLast = dynamic_cast<A_BROADCAST_GPS*>(tabEphLastGPS[isat]);
    if ( (ephLast==NULL) || diffDate(&eph->_TOC, &ephLast->getEph()._TOC) > 0.0) 
    {
      delete tabEphPrecGPS[isat];
      tabEphPrecGPS[isat] = ephLast;
      A_BROADCAST_GPS* brdc_gps=new A_BROADCAST_GPS(*eph);
      tabEphLastGPS[isat] = brdc_gps;

      if (satBiasesPriority <= 2)
      {      
	delete [] satBiasesGPS._biases;

	rtrover_satBiases satBiases;
	satBiases._satellite = eph->_satellite;
	satBiases._time = eph->_TOC;
	satBiases._numBiases = 1;

	rtrover_bias biases;
	biases._rnxType3ch[0] = 'C';
	biases._rnxType3ch[1] = '1';
	biases._rnxType3ch[2] = 'C';
	biases._value = -eph->_TGD*CLIGHT;
	satBiases._biases = &biases;

	satBiasesGPS = satBiases;
	satBiasesGPS._biases = new rtrover_bias[satBiases._numBiases];
          for(int j=0; j<satBiases._numBiases; j++)
            satBiasesGPS._biases[j] = satBiases._biases[j];
	satBiasesPriority = 2;
      }

    }
  }
}

RTROVER_API void rtrover_putGloEphemeris(const rtrover_ephGlo* eph)
{
  PRINT_LOG(log_rtrover,"rtrover_putGloEphemeris %d %lf %c%02d\n", eph->_timeUTC._mjd, eph->_timeUTC._sec, eph->_satellite._system, eph->_satellite._number);

  int isat=eph->_satellite._number-1;
  if ((isat != -1) && (isat < getMaxSatSyst(Glo)))
  {
    std::vector<A_BROADCAST*> &tabEphLastGlo = tabEphLast[Glo];
    std::vector<A_BROADCAST*> &tabEphPrecGlo = tabEphPrec[Glo];

    A_BROADCAST_Glo* ephLast = dynamic_cast<A_BROADCAST_Glo*>(tabEphLastGlo[isat]);
    if ((ephLast==NULL) || diffDate(&eph->_timeUTC, &ephLast->getEph()._timeUTC) > 0.0)
    {
      delete tabEphPrecGlo[isat];
      tabEphPrecGlo[isat] = ephLast;
      A_BROADCAST_Glo* brdc_glo=new A_BROADCAST_Glo(*eph);    
      tabEphLastGlo[isat] = brdc_glo;
    }
  }
}

RTROVER_API void rtrover_putGalEphemeris(const rtrover_ephGal* eph)
{
  PRINT_LOG(log_rtrover,"rtrover_putGalEphemeris %d %lf %c%02d %d\n", eph->_TOC._mjd, eph->_TOC._sec, eph->_satellite._system, eph->_satellite._number, eph->_IODNav);

  int isat=eph->_satellite._number-1;
  if ((isat != -1) && (isat < getMaxSatSyst(Gal)))
  {
    std::vector<A_BROADCAST*> &tabEphLastGal = tabEphLast[Gal];
    std::vector<A_BROADCAST*> &tabEphPrecGal = tabEphPrec[Gal];
    rtrover_satBiases &satBiasesGal = tabSatBiases[Gal][isat];
    int &satBiasesPriority = tabSatBiasesPriority[Gal][isat];

    A_BROADCAST_Gal* ephLast = dynamic_cast<A_BROADCAST_Gal*>(tabEphLastGal[isat]);
    if ( (ephLast==NULL) || diffDate(&eph->_TOC, &ephLast->getEph()._TOC) > 0.0) 
    {
      delete tabEphPrecGal[isat];
      tabEphPrecGal[isat] = ephLast;
      A_BROADCAST_Gal* brdc_gal=new A_BROADCAST_Gal(*eph);
      tabEphLastGal[isat] = brdc_gal;

      if (satBiasesPriority <= 2)
      {      
	delete [] satBiasesGal._biases;

	rtrover_satBiases satBiases;
	satBiases._satellite = eph->_satellite;
	satBiases._time = eph->_TOC;
	satBiases._numBiases = 3;
	
	rtrover_bias biases[satBiases._numBiases];
	biases[0]._rnxType3ch[0] = 'C';
	biases[0]._rnxType3ch[1] = '1';
	biases[0]._rnxType3ch[2] = 'C';
	biases[0]._value = - eph->_BGD_1_5A * CLIGHT;
	
	biases[1]._rnxType3ch[0] = 'C';
	biases[1]._rnxType3ch[1] = '5';
	biases[1]._rnxType3ch[2] = 'Q';
	biases[1]._value = - FREQUENCY::GAMMA(F5, Gal, 0) * eph->_BGD_1_5A * CLIGHT;
	
	biases[2]._rnxType3ch[0] = 'C';
	biases[2]._rnxType3ch[1] = '7';
	biases[2]._rnxType3ch[2] = 'Q';
	biases[2]._value =(( 1.0 - FREQUENCY::GAMMA(F7, Gal, 0)) *  eph->_BGD_1_5B - eph->_BGD_1_5A ) * CLIGHT;
	
	satBiases._biases = biases;

	satBiasesGal = satBiases;
	satBiasesGal._biases = new rtrover_bias[satBiases._numBiases];
          for(int j=0; j<satBiases._numBiases; j++)
            satBiasesGal._biases[j] = satBiases._biases[j];
	satBiasesPriority = 2;
      }

    }
  }
}

RTROVER_API void rtrover_putBdsEphemeris(const rtrover_ephBds* eph)
{
  PRINT_LOG(log_rtrover,"rtrover_putBdsEphemeris %d %lf %c%02d %d\n", eph->_TOC._mjd, eph->_TOC._sec, eph->_satellite._system, eph->_satellite._number, eph->_IODE);

  int isat=eph->_satellite._number-1;
  if ((isat != -1) && (isat < getMaxSatSyst(Bds)))
  {
    std::vector<A_BROADCAST*> &tabEphLastBds = tabEphLast[Bds];
    std::vector<A_BROADCAST*> &tabEphPrecBds = tabEphPrec[Bds];

    A_BROADCAST_Bds* ephLast = dynamic_cast<A_BROADCAST_Bds*>(tabEphLastBds[isat]);
    if ( (ephLast==NULL) || diffDate(&eph->_TOC, &ephLast->getEph()._TOC) > 0.0) 
    {
      delete tabEphPrecBds[isat];
      tabEphPrecBds[isat] = ephLast;
      A_BROADCAST_Bds* brdc_bds=new A_BROADCAST_Bds(*eph);
      tabEphLastBds[isat] = brdc_bds;
    }
  }
}

RTROVER_API void rtrover_putOrbCorrections(int numCorr, const rtrover_orbCorr* corr)
{
  FOR_LOG (int i=0;i<numCorr;i++)
  {   
    PRINT_LOG(log_rtrover,"rtrover_putOrbCorrections %d %lf %c%02d %d\n", corr[i]._time._mjd, corr[i]._time._sec,\
          corr[i]._satellite._system, corr[i]._satellite._number, corr[i]._iod);  
  }

  for (int i=0;i<numCorr;i++)
  {
    int isat = -1;
    int syst=convertSystem(corr[i]._satellite._system);
    if (corr[i]._satellite._number >= 1) 
      isat=corr[i]._satellite._number-1;
    if ((isat != -1) && (isat < getMaxSatSyst(syst)) && (syst != -1))
      tabOrbCorr[syst][isat] = corr[i];    
  }
}

RTROVER_API void rtrover_putClkCorrections(int numCorr, const rtrover_clkCorr* corr)
{
  FOR_LOG (int i=0;i<numCorr;i++)
  {
    PRINT_LOG(log_rtrover,"rtrover_putClkCorrections %d %lf %c%02d %d %lf %lf %lf\n", corr[i]._time._mjd, corr[i]._time._sec,\
          corr[i]._satellite._system, corr[i]._satellite._number, corr[i]._iod, corr[i]._dClk, corr[i]._dotDClk, corr[i]._dotDotDClk);
  }

  for (int i=0;i<numCorr;i++)
  {
    int isat = -1;
    int syst=convertSystem(corr[i]._satellite._system);
    if (corr[i]._satellite._number >= 1) 
      isat=corr[i]._satellite._number-1;
    if ((isat != -1) && (isat < getMaxSatSyst(syst)) && (syst != -1))
      tabClkCorr[syst][isat] = corr[i];
  }
}

RTROVER_API void rtrover_putBiases(int numBiases, const rtrover_satBiases* biases)
{
  FOR_LOG (int i=0;i<numBiases;i++)
  {
    PRINT_LOG(log_rtrover,"rtrover_putBiases %d %lf %c%02d %d: ", biases[i]._time._mjd, biases[i]._time._sec,\
            biases[i]._satellite._system, biases[i]._satellite._number, biases[i]._numBiases);
    FOR_LOG (int j=0;j<biases[i]._numBiases;j++)
    {
      PRINT_LOG(log_rtrover,"%c%c%c %15.3f ",
              biases[i]._biases[j]._rnxType3ch[0], biases[i]._biases[j]._rnxType3ch[1], biases[i]._biases[j]._rnxType3ch[2],
              biases[i]._biases[j]._value);
    }
    PRINT_LOG(log_rtrover,"\n");
  }

  for (int i=0;i<numBiases;i++)
  {
    int isat = -1;
    int syst=convertSystem(biases[i]._satellite._system);
    
    if (biases[i]._satellite._number >= 1) 
      isat=biases[i]._satellite._number-1;
      
    if ((isat != -1) && (isat < getMaxSatSyst(syst)) && (syst != -1))
    {
      int &satBiasesPriority = tabSatBiasesPriority[syst][isat];
      rtrover_satBiases &satBiases = tabSatBiases[syst][isat] ;
      
      if (satBiasesPriority <= 3)
      {
	delete [] satBiases._biases;
        satBiases = biases[i];
	satBiases._biases = new rtrover_bias[biases[i]._numBiases];
        for(int j=0; j<biases[i]._numBiases; j++)
          satBiases._biases[j] = biases[i]._biases[j];
        satBiasesPriority = 3;
      }
    }
  }
}

RTROVER_API void rtrover_putPhaseBiases(int numBiases, const rtrover_satPhaseBiases* biases)
{
  FOR_LOG (int i=0;i<numBiases;i++)
  {
    PRINT_LOG(log_rtrover,"rtrover_putPhaseBiases %d %lf %c%02d %8.3f %d: ", biases[i]._time._mjd, biases[i]._time._sec,\
          biases[i]._satellite._system, biases[i]._satellite._number, biases[i]._yawAngleValue, biases[i]._numBiases);
    FOR_LOG (int j=0;j<biases[i]._numBiases;j++)
    {
      PRINT_LOG(log_rtrover,"%c%c%c %8.3f %d %d %d %d  ",
            biases[i]._biases[j]._rnxType3ch[0], biases[i]._biases[j]._rnxType3ch[1], biases[i]._biases[j]._rnxType3ch[2],\
            biases[i]._biases[j]._value, biases[i]._biases[j]._discontinuityValue,//value du biais de phase en cycles
	    biases[i]._biases[j]._intIndicator,biases[i]._biases[j]._wlIntIndicator, biases[i]._biases[j]._discontinuityValue);
    }
    PRINT_LOG(log_rtrover,"\n");
  }

  for (int i=0;i<numBiases;i++)
  {
    int isat = -1;
    int syst=convertSystem(biases[i]._satellite._system);
    if (((syst == GPS) || (syst == Gal) || (syst == Bds)) && (biases[i]._satellite._number >= 1))
      isat=biases[i]._satellite._number-1;         
      
    if ((isat != -1) && (isat < getMaxSatSyst(syst)) && (syst != -1))
    {
      rtrover_satPhaseBiases &satPhaseBiases = tabSatPhaseBiases[syst][isat] ;
      
      delete [] satPhaseBiases._biases;
      satPhaseBiases = biases[i];

      //biais phase de m en cycles
      for (int j=0;j<biases[i]._numBiases;j++)
      {
        if (biases[i]._biases[j]._rnxType3ch[1] == '1')
          biases[i]._biases[j]._value /= FREQUENCY::LAMBDA(F1, syst, 0);
        if (biases[i]._biases[j]._rnxType3ch[1] == '2')
          biases[i]._biases[j]._value /= FREQUENCY::LAMBDA(F2, GPS, 0);
        if (biases[i]._biases[j]._rnxType3ch[1] == '5')
          biases[i]._biases[j]._value /= FREQUENCY::LAMBDA(F5, GPS, 0);
	if (biases[i]._biases[j]._rnxType3ch[1] == '6')
          biases[i]._biases[j]._value /= FREQUENCY::LAMBDA(F6, Bds, 0);
	if (biases[i]._biases[j]._rnxType3ch[1] == '7')
          biases[i]._biases[j]._value /= FREQUENCY::LAMBDA(F7, syst, 0);
      }
    
      satPhaseBiases._biases = new rtrover_phaseBias[biases[i]._numBiases];
      for(int j=0; j<biases[i]._numBiases; j++)
          satPhaseBiases._biases[j] = biases[i]._biases[j];
    }
  }
}


RTROVER_API void rtrover_putDCBBiases(int numBiases, const rtrover_satDCBBiases* biases)
{
  PRINT_LOG(log_rtrover,"rtrover_putDCBBiases\n");

  for (int i=0;i<numBiases;i++)
  {
    if ((biases[i]._satellite._system == 'G') && (biases[i]._satellite._number >= 1))
    {
      int isat=biases[i]._satellite._number-1;
      if ((isat != -1) && (isat < getMaxSatSyst(GPS)))
      {
	int &satBiasesPriority = tabSatBiasesPriority[GPS][isat];
	rtrover_satBiases &satBiasesGPS = tabSatBiases[GPS][isat] ;

	if (satBiasesPriority <= 1)
	{
	  double P1P2 = biases[i]._P1P2;
	  double P1C1 = biases[i]._P1C1;
	  double P2C2 = biases[i]._P2C2;
	  const static double a_L1_GPS = -FREQUENCY::GAMMA(F2, GPS, 0)/(1.0-FREQUENCY::GAMMA(F2, GPS, 0));
	  const static double a_L2_GPS = 1.0/(1.0-FREQUENCY::GAMMA(F2, GPS, 0));
	  delete [] satBiasesGPS._biases;

	  rtrover_satBiases satBiases;
	  satBiases._satellite = biases[i]._satellite;
	  satBiases._time = biases[i]._time;
	  satBiases._numBiases = 4;

	  rtrover_bias biases[4];
	  biases[0]._rnxType3ch[0] = 'C';
	  biases[0]._rnxType3ch[1] = '1';
	  biases[0]._rnxType3ch[2] = 'C';
	  biases[0]._value = 0.0;
	  if (P1P2 && P1C1)
            biases[0]._value = -a_L2_GPS*P1P2+P1C1;
	  biases[1]._rnxType3ch[0] = 'C';
	  biases[1]._rnxType3ch[1] = '1';
	  biases[1]._rnxType3ch[2] = 'W';
	  biases[1]._value = -a_L2_GPS*P1P2;
	  biases[2]._rnxType3ch[0] = 'C';
	  biases[2]._rnxType3ch[1] = '2';
	  biases[2]._rnxType3ch[2] = 'W';
	  biases[2]._value = a_L1_GPS*P1P2;
	  biases[3]._rnxType3ch[0] = 'C';
	  biases[3]._rnxType3ch[1] = '2';
	  biases[3]._rnxType3ch[2] = 'X';
	  biases[3]._value = 0.0;
	  if (P1P2 && P2C2)
            biases[3]._value = a_L1_GPS*P1P2+P2C2;

	  satBiases._biases = biases;

	  satBiasesGPS = satBiases;
	  satBiasesGPS._biases = new rtrover_bias[satBiases._numBiases];
          for(int j=0; j<satBiases._numBiases; j++)
            satBiasesGPS._biases[j] = satBiases._biases[j];
	  satBiasesPriority = 1;
	}
      }
    }
  }

}

RTROVER_API void rtrover_putIonoCorrs(int numCorr, const rtrover_ionoCorr* corr)
{
  FOR_LOG (int i=0;i<numCorr;i++)
  {
    PRINT_LOG(log_rtrover,"rtrover_putIonoCorrs %c%02d %lf\n", \
            corr[i]._satellite._system, corr[i]._satellite._number, corr[i]._value);
  }  
  
  for (int isyst=0;isyst<SystMax;isyst++)
  { 
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      tabRover[currentRover].initializeIonoCorr( tabRover[currentRover]._tabIonoCorr[isyst][isat]);
    }
  }
  
  for (int i=0;i<numCorr;i++)
  {
    int isat = -1;
    int syst=convertSystem(corr[i]._satellite._system);
    if (corr[i]._satellite._number >= 1) 
      isat=corr[i]._satellite._number-1;
    if ((isat != -1) && (isat < getMaxSatSyst(syst)) && (syst != -1))
      tabRover[currentRover]._tabIonoCorr[syst][isat] = corr[i];
  }
}

RTROVER_API void rtrover_putSbasCorrs(int numCorr, const rtrover_sbasCorr* corr)
{
  FOR_LOG (int i=0;i<numCorr;i++)
  {
    PRINT_LOG(log_rtrover,"rtrover_putSbasCorrs %c%02d %lf %lf %lf %lf\n", \
          corr[i]._satellite._system,
		  corr[i]._satellite._number, corr[i]._orb[0], corr[i]._satellite._number, corr[i]._orb[1], corr[i]._satellite._number, corr[i]._orb[2],
		  corr[i]._satellite._number, corr[i]._clk);
  }

  for (int i=0;i<numCorr;i++)
  {
    int isat = -1;
    int syst=convertSystem(corr[i]._satellite._system);
    if (corr[i]._satellite._number >= 1) 
      isat=corr[i]._satellite._number-1;
    if ((isat != -1) && (isat < getMaxSatSyst(syst)) && (syst != -1))
      tabSbasCorr[syst][isat] = corr[i];
  }
}

RTROVER_API void rtrover_putTropoCorrs(const rtrover_tropoCorr* corr)
{
  PRINT_LOG(log_rtrover,"rtrover_putTropoCorrs %lf\n", corr->_value);
  tabRover[currentRover]._tropoCorr._value = corr->_value;
}

RTROVER_API void rtrover_destroy()
{
  PRINT_LOG(log_rtrover,"rtrover_destroy\n");
    
  for (int isyst=0; isyst < SystMax ; isyst++)
  {
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      delete [] tabSatBiases[isyst][isat]._biases;
      delete [] tabSatPhaseBiases[isyst][isat]._biases;
      delete tabEphPrec[isyst][isat];
      delete tabEphLast[isyst][isat];  
    }
  }
  
  CLOSE_LOG(log_rtrover);
  if (lowlevel && fileOpened)
    fclose(log_lowlevel);  
}

RTROVER_API void rtrover_freeOutput(rtrover_output* output)
{
  PRINT_LOG(log_rtrover,"rtrover_freeOutput\n");
}

RTROVER_API 
void rtrover_processEpoch(int numSatRover,
                          const rtrover_satObs* satObsRover,
                          int numSatBase,
                          const rtrover_satObs* satObsBase,
                          rtrover_output* output
                          )
{
  static char str_output[MAX_EPOCH*200];

  FOR_LOG (int i=0;i<numSatRover;i++)
  {
    PRINT_LOG(log_rtrover,"rtrover_processEpoch Rover=%d, %d %lf %c%02d ", currentRover, satObsRover[i]._time._mjd, satObsRover[i]._time._sec, \
         satObsRover[i]._satellite._system, satObsRover[i]._satellite._number);
    FOR_LOG (int j=0;j<satObsRover[i]._numObs;j++)
    {
      PRINT_LOG(log_rtrover,"%c%c %15.3f %15.3f ",
              satObsRover[i]._obs[j]._rnxType2ch[0], satObsRover[i]._obs[j]._rnxType2ch[1],
              (satObsRover[i]._obs[j]._codeValid?satObsRover[i]._obs[j]._code:0.0),
              (satObsRover[i]._obs[j]._phaseValid?satObsRover[i]._obs[j]._phase:0.0));
    }
    PRINT_LOG(log_rtrover,"\n");
  }

  // Shift of array of observations
  if (tabRover[currentRover]._nbEpoch < MAX_EPOCH)
    tabRover[currentRover]._nbEpoch++;
  else
  {
    tabRover[currentRover].deleteTabSatObs(MAX_EPOCH-1);
  }
  for (int j=tabRover[currentRover]._nbEpoch-1;j>0;j--)
  {
    for (int isyst=0;isyst<SystMax;isyst++)
    { 
      for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
      {
        tabRover[currentRover]._tabSatObs[j][isyst][isat] = tabRover[currentRover]._tabSatObs[j-1][isyst][isat];
      }
    }
  }
  for (int isyst=0;isyst<SystMax;isyst++)
  { 
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      tabRover[currentRover].initializeSatObs(tabRover[currentRover]._tabSatObs[0][isyst][isat]);
    }
  }
  
  // Observations storage (index 0)
  tabRover[currentRover].newTabSatObs(0,numSatRover,satObsRover,tabEphLast[Glo]);

  // Process of array of epochs
  *str_output='\0';
  if (correct && (correctType == RTIGS))
  {
    for (int j=tabRover[currentRover]._nbEpoch-1;j>=0;j--)
    {
      double dto=0.0;
      double dtc=0.0;
      std::vector<std::vector<rtrover_satObs> > &tabSatObsEpoch = tabRover[currentRover]._tabSatObs[j];
      
      for (int isyst=0;isyst<SystMax;isyst++)
      { 
        if ((isyst == GPS && A_FILTER::_filterSetting._gps._use) || (isyst == Glo && A_FILTER::_filterSetting._glo._use) || (isyst == Gal && A_FILTER::_filterSetting._gal._use) || (isyst == Bds && A_FILTER::_filterSetting._bds._use))
        {
	  for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
          {
	    rtrover_satObs &satObs = tabSatObsEpoch[isyst][isat];
	    rtrover_clkCorr &clkCorr = tabClkCorr[isyst][isat];
	    if (satObs._satellite._number)
            {
	      double t = satObs._time._mjd * 86400.0 + satObs._time._sec;
	      if (t > dto)
		dto = t;
	    }
	    if (clkCorr._satellite._number)
            {
	      double t = clkCorr._time._mjd * 86400.0 + clkCorr._time._sec;
	      if (t > dtc)
		dtc = t;
	    }
          }
	}
      }
      if (!dto || !dtc)
      {
	tabRover[currentRover].deleteTabSatObs(j);
        tabRover[currentRover]._nbEpoch--;
        continue;
      }
      double dt=dto-dtc;
      if (dt < 0.0)
      {
	tabRover[currentRover].deleteTabSatObs(j);
        tabRover[currentRover]._nbEpoch--;
        continue;
      }
      if (dt < maxAge)
      {
	processEpoch(tabSatObsEpoch, str_output);
        tabRover[currentRover].deleteTabSatObs(j);
        tabRover[currentRover]._nbEpoch--;
        continue;
      }
      else
        break;
    }
  }
  else
  {
    for (int j=tabRover[currentRover]._nbEpoch-1;j>=0;j--)
    {
      std::vector<std::vector<rtrover_satObs> > &tabSatObsEpoch = tabRover[currentRover]._tabSatObs[j];
      processEpoch(tabSatObsEpoch, str_output);
      tabRover[currentRover].deleteTabSatObs(j);
      tabRover[currentRover]._nbEpoch--;
    }
  }

  output->_epoTime._mjd = 0;
  output->_epoTime._sec = 0.0;
  for (int i=0;i<3;i++)
    output->_xyzRover[i] = 0.0;
  for (int i=0;i<6;i++)
    output->_covMatrix[i] = 0.0;
  output->_log=str_output;
  output->_error = false;
}

/// Get rover position
RTROVER_API void rtrover_getRoverInformation(int *day, double *sec, double pos[6], double tropo[3])
{
  PRINT_LOG(log_rtrover,"rtrover_getRoverPosition\n");
  
  *day = 0;
  if (tabRover[currentRover]._lastTime._mjd)
    *day = tabRover[currentRover]._lastTime._mjd-33282;
  *sec = tabRover[currentRover]._lastTime._sec;
  pos[0] = tabRover[currentRover]._stationProduct._x;
  pos[1] = tabRover[currentRover]._stationProduct._y;
  pos[2] = tabRover[currentRover]._stationProduct._z;
  pos[3] = tabRover[currentRover]._stationProduct._covX;
  pos[4] = tabRover[currentRover]._stationProduct._covY;
  pos[5] = tabRover[currentRover]._stationProduct._covZ;
  tropo[0] = tabRover[currentRover]._stationProduct._tropoV;
  tropo[1] = tabRover[currentRover]._stationProduct._tropo;
  tropo[2] = tabRover[currentRover]._stationProduct._covTropo;
}

/// Get satellite position
RTROVER_API void rtrover_getSatelliteInformation(rtrover_satellite satellite, double pos[3], double amb[6], double iono[2])
{
  int isat = -1;

  pos[0] = 0.0;
  pos[1] = 0.0;
  pos[2] = 0.0;  
  int syst=convertSystem(satellite._system);
  if (satellite._number >= 1) 
    isat=satellite._number-1;
  if ((isat != -1) && (isat < getMaxSatSyst(syst)))
  {
    A_VECTOR3D &satPosition = tabSatPosition[syst][isat];
    pos[0] = satPosition.getX();
    pos[1] = satPosition.getY();
    pos[2] = satPosition.getZ();
    A_MEASUREMENT_PRODUCT &satMeasurement = tabRover[currentRover]._stationMeasurement[syst][isat];
    amb[0] = satMeasurement._ambEx;
    amb[1] = satMeasurement._ambA4;
    amb[2] = satMeasurement._ambP;
    amb[3] = satMeasurement._covAmbEx;
    amb[4] = satMeasurement._covAmbA4;
    amb[5] = satMeasurement._covAmbP;
    iono[0] = satMeasurement._iono;
    iono[1] = satMeasurement._covIono;
  }
}

