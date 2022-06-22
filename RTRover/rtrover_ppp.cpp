/************************************************************
Nom ......... : rtrover_ppp.cpp
Role ........ : ppp definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#include <cstdio>

#include "rtrover_vector3D.h"
#include "rtrover_ppp.h"
#include "rtrover_utility.h"
#include "rtrover_frequency.h"
#include "rtrover_filter.h"
#include "rtrover_givens_filter.h"
#include "rtrover_model.h"
#include "rtrover_date.h"
#include "rtrover_vector3D.h"

#define MAX_PARAM (3+SystMax)
#define COARSE_THR 50.0
#define COARSE_MAX 3

//////////////////////////////////////////////////////////////////////
static int givensFilter(std::vector<A_FILTER_MODELED_MEASUREMENT> &measurements, const double re_ResMax, const int maxElim, double *sol)
{
  double variance;
  double dResMax;
  int iMes_ResMax, iMes, iInc;
  int iter = 0;

  while (1)
  {
    int nb = 0;
    A_GIVENS_FILTER filter(MAX_PARAM);
    for (iMes=0; iMes<(int)measurements.size(); iMes++)
      if (measurements[iMes]._w)
      {
        filter.newMeasurement(measurements[iMes]);
        nb++;
      }
    iter++;

    if (iter > maxElim)
      break;
    if (nb < 4)
      break;
    filter.solution(sol, variance);
    if (!re_ResMax)
      return 1;
    iMes_ResMax = 0;
    dResMax = 0;
    for (iMes=0; iMes<(int)measurements.size(); iMes++)
    {
      if (measurements[iMes]._w)
      {
        double r = measurements[iMes]._res;
        for (iInc=0; iInc<(int)measurements[iMes]._der.size() ; iInc++)
        r -= sol[iInc]*measurements[iMes]._der[iInc];
        if (fabs(r) > dResMax)
        {
          iMes_ResMax = iMes;
          dResMax = fabs(r);
        }
      }
    }
    if (dResMax > re_ResMax)
    {
      filter.newMeasurement(measurements[iMes_ResMax]);
      measurements[iMes_ResMax]._w = 0.0;
      measurements[iMes_ResMax]._res = dResMax;
      continue;
    }
    return 1;
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
A_PPP::A_PPP()
{
  _model = std::vector<std::vector<A_MODEL> >(SystMax);
  for (int isyst=0; isyst<SystMax;isyst++)
    _model[isyst] = std::vector<A_MODEL>(getMaxSatSyst(isyst));
  
  _smoothedMeas = std::vector<A_SMOOTHED_MEASUREMENT>(FMax);
  
  _state = A_VECTOR3D();
   
  _bInit=0;

}

//////////////////////////////////////////////////////////////////////
void A_PPP::computeModel(const A_DATE &date, const A_VECTOR3D& MRec,
                         const std::vector<std::vector<A_MEASUREMENT> >& measurements,
			 const std::vector<std::vector<A_SSR_PARAMETER> >& positions,
                         std::vector<std::vector<A_MODELED_MEASUREMENT> >& mod)
{
  for (int isyst=0; isyst<SystMax;isyst++)
  {
    if ((isyst == GPS && A_FILTER::_filterSetting._gps._use) || (isyst == Glo && A_FILTER::_filterSetting._glo._use) || (isyst == Gal && A_FILTER::_filterSetting._gal._use) || (isyst == Bds && A_FILTER::_filterSetting._bds._use))
    {
      for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
      {
	mod[isyst][isat] = _model[isyst][isat].model(date, isyst, MRec, measurements[isyst][isat], positions[isyst][isat], _settings._map.getSatellites()[isyst][isat]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
void A_PPP::coarsePoint(const std::vector<std::vector<A_MEASUREMENT> >& measurements,
			const std::vector<std::vector<A_SSR_PARAMETER> >& positions,
			const A_DATE &date,
			A_VECTOR3D& point)
{
  double x;
  A_FILTER_MODELED_MEASUREMENT measurement;
  std::vector<A_FILTER_MODELED_MEASUREMENT> measurements_sat;
  double result[MAX_PARAM];
  A_VECTOR3D posRec;
  std::vector<std::vector<A_MODELED_MEASUREMENT> > mod(SystMax);
  std::vector<std::vector<A_MEASUREMENT> > measurementsBeidou;
  
  point=A_VECTOR3D();
  for (int isyst=0; isyst<SystMax;isyst++)
    mod[isyst] = std::vector<A_MODELED_MEASUREMENT>(getMaxSatSyst(isyst));
  
  for (int isyst=0; isyst<SystMax;isyst++)
  {
    if ((isyst == GPS && A_FILTER::_filterSetting._gps._use) || (isyst == Glo && A_FILTER::_filterSetting._glo._use) || (isyst == Gal && A_FILTER::_filterSetting._gal._use) || (isyst == Bds && A_FILTER::_filterSetting._bds._use))
    {
      for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
      {
	const A_SSR_PARAMETER &positionsSat=positions[isyst][isat];

	if (positionsSat._Xsp3 && positionsSat._Ysp3 && positionsSat._Zsp3)
	{
          A_VECTOR3D posSp3(positionsSat._Xsp3,positionsSat._Ysp3,positionsSat._Zsp3);
          posRec += posSp3;
	}
      }
    }
  }
  if (!posRec.getX() && !posRec.getY() && !posRec.getZ())
    return;

  x = posRec.getModule()/6.4e6;
  posRec /= x; 
  //fprintf(stderr, "POS %d %lf %lf %lf %lf\n", date.getDay(), date.getSecond(), posRec.getX(), posRec.getY(), posRec.getZ());  
  
  for (int i=0;i<5;i++)
  {
    computeModel(date, posRec, measurements, positions, mod);
    beidouCorrection(measurements, measurementsBeidou, mod);
    measurements_sat.clear();
    for (int isyst=0; isyst<SystMax;isyst++)
    {
      if ((isyst == GPS && A_FILTER::_filterSetting._gps._use) || (isyst == Glo && A_FILTER::_filterSetting._glo._use) || (isyst == Gal && A_FILTER::_filterSetting._gal._use) || (isyst == Bds && A_FILTER::_filterSetting._bds._use))
      {	
	for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
        {
          const A_SSR_PARAMETER &positionsSat=positions[isyst][isat];
          const A_MEASUREMENT &measurementsSat=measurementsBeidou[isyst][isat];
          const A_MODELED_MEASUREMENT &modSat=mod[isyst][isat];
  
          for (int f=F1; f<FMax; f++)
	  {
	    if (measurementsSat._code[f] && modSat._modCode[f] && positionsSat._Hsp3 && modSat._pos.getX() && modSat._pos.getY() && modSat._pos.getZ())
            {
              measurement._res = measurementsSat._code[f] - modSat._modCode[f] + positionsSat._Hsp3;
              measurement._w = 1.0;
              measurement._der.clear();
              measurement._der.resize(MAX_PARAM);
              measurement._der[0] = modSat._pos.getX();
              measurement._der[1] = modSat._pos.getY();
              measurement._der[2] = modSat._pos.getZ();
       	      measurement._der[3+isyst] = 1.0;
              measurements_sat.push_back(measurement);
            }
	  }
        }
      }
    }
    if (!givensFilter(measurements_sat, 0.0, COARSE_MAX, result))
      return;

    A_VECTOR3D result_vect(result[0],result[1],result[2]);
    if (result_vect.getModule() < 0.01)
    {
      if (!givensFilter(measurements_sat, COARSE_THR, COARSE_MAX, result))
	return;
      result_vect=A_VECTOR3D(result[0],result[1],result[2]);
      point=posRec+result_vect;
//  fprintf(stderr, "POS %d %lf %lf %lf %lf\n", date.getDay(), date.getSecond(), point.getX(), point.getY(), point.getZ());
      return;
    }
    
    posRec += result_vect;
  }
}

//////////////////////////////////////////////////////////////////////
static void bdsCodeCorr(int typeBds, double elevation, int freq, double &codeCorr)
{
  /* Elevation-dependent correction values for Beidou code measurements */
  //table for B1, B2 and B3
  static double elev_IGSO[3][10] = {{-0.55, -0.40, -0.34, -0.23, -0.15, -0.04, 0.09, 0.19, 0.27, 0.35},
                                    {-0.71, -0.36, -0.33, -0.19, -0.14, -0.03, 0.08, 0.17, 0.24, 0.33},
                                    {-0.27, -0.23, -0.21, -0.15, -0.11, -0.04, 0.05, 0.14, 0.19, 0.32}};
  static double elev_MEO[3][10] = {{-0.47, -0.38, -0.32, -0.23, -0.11, 0.06, 0.34, 0.69, 0.97, 1.05},
                                   {-0.40, -0.31, -0.26, -0.18, -0.06, 0.09, 0.28, 0.48, 0.64, 0.69},
                                   {-0.22, -0.15, -0.13, -0.10, -0.04, 0.05, 0.14, 0.27, 0.36, 0.47}};

   
  //Elevation in degrees
  int tElev;
  tElev = floor((elevation/10.0)+0.5);

  if (typeBds == IGSO)
  {
    switch (freq)
    {
      case F1 :
        codeCorr = elev_IGSO[0][tElev];
	break;
      case F6 :
        codeCorr = elev_IGSO[2][tElev];
	break;
      case F7 :
        codeCorr = elev_IGSO[1][tElev];
	break;     
      default :
        codeCorr = 0.0;	
    }
  }
  if (typeBds == MEO)
  {
    switch (freq)
    {
      case F1 :
        codeCorr = elev_MEO[0][tElev];
	break;
      case F6 :
        codeCorr = elev_MEO[2][tElev];
	break;
      case F7 :
        codeCorr = elev_MEO[1][tElev];
	break;     
      default :
        codeCorr = 0.0;	
    }
  }	 
}

//////////////////////////////////////////////////////////////////////
void A_PPP::smoothMeasurements(const A_DATE &date, const std::vector<std::vector<A_MEASUREMENT> >& raw_meas,
                                     std::vector<std::vector<A_MEASUREMENT> >& smooth_meas)
{
  smooth_meas = raw_meas;
  for(int f=F1;f<FMax;f++)
  {
    if(_settings._mode._freq[f])
    {
      _smoothedMeas[f].smoothMeasurements(date, raw_meas, _settings._smooth,(Frequency)f);
      for (int isyst=0; isyst<SystMax;isyst++)
        for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
          smooth_meas[isyst][isat]._code[f]=_smoothedMeas[f].getPl()[isyst][isat];
    }
  }	
}

//////////////////////////////////////////////////////////////////////
void A_PPP::beidouCorrection(const std::vector<std::vector<A_MEASUREMENT> >& meas,
                             std::vector<std::vector<A_MEASUREMENT> >& measBeidou,
			     const std::vector<std::vector<A_MODELED_MEASUREMENT> >& mod)
{
  measBeidou = meas;
  for (int isyst=0; isyst<SystMax;isyst++)
  {
    if ((isyst == Bds) && A_FILTER::_filterSetting._bds._use)
    {
      for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
      {
        const A_MEASUREMENT &mesSat=meas[isyst][isat];
        const A_MODELED_MEASUREMENT &modSat=mod[isyst][isat];
	double codeCorr = 0.0;
	for(int f=F1;f<FMax;f++)
	{
	  if (measBeidou[isyst][isat]._code[f] && modSat._elevation)
	  {
	    bdsCodeCorr(mesSat._slot_typeBds, modSat._elevation, f, codeCorr);
	    measBeidou[isyst][isat]._code[f] += codeCorr;
	  }
	}
      }
    }	
  }
}

//////////////////////////////////////////////////////////////////////
void A_PPP::computeResiduals(std::vector<std::vector<A_RESIDUAL> >& residuals,
                      const std::vector<std::vector<A_BIAS> >& bias,
                      const std::vector<std::vector<A_SSR_PARAMETER> >& positions,
                      const std::vector<std::vector<A_MEASUREMENT> >& mes,
                      const std::vector<std::vector<A_MODELED_MEASUREMENT> >& mod)
{

  A_MODE &mode=_settings._mode;
  
  for (int isyst=0; isyst<SystMax;isyst++)
  {
    if ((isyst == GPS && A_FILTER::_filterSetting._gps._use) || (isyst == Glo && A_FILTER::_filterSetting._glo._use) || (isyst == Gal && A_FILTER::_filterSetting._gal._use) || (isyst == Bds && A_FILTER::_filterSetting._bds._use))
    {
      
      for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
      {
        A_RESIDUAL &residualsSat=residuals[isyst][isat];
        const A_BIAS &biasSat=bias[isyst][isat];
        const A_SSR_PARAMETER &positionsSat=positions[isyst][isat];
        const A_MEASUREMENT &mesSat=mes[isyst][isat];
        const A_MODELED_MEASUREMENT &modSat=mod[isyst][isat];
	
        if ((modSat._mapping < _settings._thrMap) && positionsSat._Hsp3)
        {
          for (int f=F1; f<FMax; f++)
          {
	    if (mode._ambig && (biasSat._wli || biasSat._n1i))
	    {
              if (mode._freq[f] && mode._code && mesSat._code[f] && modSat._modCode[f] && biasSat._code[f])
                residualsSat._resCode[f] = mesSat._code[f] - modSat._modCode[f] + biasSat._code[f] + positionsSat._Hsp3;
              if (mode._freq[f] && mode._phase && mesSat._phase[f] && modSat._modPhase[f] && biasSat._phase[f])
                residualsSat._resPhase[f] = mesSat._phase[f] - modSat._modPhase[f] + biasSat._phase[f] + positionsSat._Hsp3/FREQUENCY::LAMBDA((Frequency)f, isyst, mesSat._slot_typeBds);
	    }
	    else
	    {
              if (mode._freq[f] && mode._code && mesSat._code[f] && modSat._modCode[f])
              {
                residualsSat._resCode[f] = mesSat._code[f] - modSat._modCode[f] + positionsSat._Hsp3;
                if (biasSat._code[f])
                  residualsSat._resCode[f] = residualsSat._resCode[f] + biasSat._code[f];
              }
              if (mode._freq[f] && mode._phase && mesSat._phase[f] && modSat._modPhase[f])
                residualsSat._resPhase[f] = mesSat._phase[f]-modSat._modPhase[f] + positionsSat._Hsp3/FREQUENCY::LAMBDA((Frequency)f, isyst, mesSat._slot_typeBds);
	    }
          }
          residualsSat._slot_typeBds = mesSat._slot_typeBds;
          residualsSat._iono = positionsSat._iono;
	  residualsSat._typeIono = positionsSat._typeIono;
          residualsSat._mapping = modSat._mapping;
          residualsSat._pos = modSat._pos;
	  residualsSat._n1i = false;
          residualsSat._wli = 0;
          residualsSat._ewli = false;
          residualsSat._discontinuity = 0;
	  if (mode._ambig)
	  {
	    residualsSat._n1i = biasSat._n1i;
            residualsSat._wli = biasSat._wli;
	    residualsSat._ewli = biasSat._wli;
	    residualsSat._discontinuity = biasSat._discontinuity;
	  }
        }
      }
    }    
  }
}

//////////////////////////////////////////////////////////////////////
double A_PPP::integrity(const std::vector<std::vector<A_RESIDUAL> >& residuals,
                 const A_RECEIVER_PRODUCT &recProduct,
		 const std::vector<std::vector<A_MEASUREMENT_PRODUCT> > &measProduct)
{
int nbMeas = 0;
int nbConst = 0;
double sqa= 0.0;

 for (int isyst=0; isyst<SystMax;isyst++) {
   bool c=false;

   for (int isat=0;isat<getMaxSatSyst(isyst);isat++) {
     double noise=0.0;
     if (isyst == GPS)
       noise = A_FILTER::_filterSetting._gps._sigCode[0]*residuals[isyst][isat]._mapping;
     else if (isyst == Glo)
       noise = A_FILTER::_filterSetting._glo._sigCode[0]*residuals[isyst][isat]._mapping;
     else if (isyst == Gal)
       noise = A_FILTER::_filterSetting._gal._sigCode[0]*residuals[isyst][isat]._mapping;
     else if (isyst == Bds) {
       if (residuals[isyst][isat]._slot_typeBds == GEO)
	 noise = A_FILTER::_filterSetting._bds._sigCode[0]*residuals[isyst][isat]._mapping;
       else if  (residuals[isyst][isat]._slot_typeBds == IGSO)
	 noise = A_FILTER::_filterSetting._bds._sigCode[1]*residuals[isyst][isat]._mapping;
       else if  (residuals[isyst][isat]._slot_typeBds == MEO)
	 noise = A_FILTER::_filterSetting._bds._sigCode[2]*residuals[isyst][isat]._mapping;
     }
    
    int n=0;
    double sqrl=0.0;
    
    for (int f=F1; f<FMax; f++) {
      if (measProduct[isyst][isat]._resCode[f] && !measProduct[isyst][isat]._elimCode[f] && noise) {
        sqrl += measProduct[isyst][isat]._resCode[f]*measProduct[isyst][isat]._resCode[f]/noise/noise;
        n++;
      }
    }
  
    if (n>=1) {
      sqa += sqrl/(double)n;
      nbMeas++;
      c=true;
    }
   }
   
   if (c)
     nbConst++;
 }

 double dop2=recProduct._covX*recProduct._covX+recProduct._covY*recProduct._covY+recProduct._covZ*recProduct._covZ;
 if (nbMeas <= 3+nbConst)
   return 0.0;
 double crit2=sqa/(double)(nbMeas-(3+nbConst))*dop2;
 // fprintf(stderr, "variance=%lf, dop2=%lf, nbMeas=%d, nbConst=%d, crit=%lf\n", sqa, dop2, nbMeas, nbConst, sqrt(crit2));
 return sqrt(crit2);
}

//////////////////////////////////////////////////////////////////////
int A_PPP::computePPP(const std::vector<std::vector<A_MEASUREMENT> >& measurements,
               const std::vector<std::vector<A_BIAS> >& biases,
               const std::vector<std::vector<A_SSR_PARAMETER> >& Positions,
	       const double tropo,
               const A_DATE &date,
               A_RECEIVER_PRODUCT& recProduct,
               std::vector<std::vector<A_MEASUREMENT_PRODUCT> >& measurementProduct,
	       const A_VECTOR3D& roverAPrioriPosition, const int rover)
{ 
  
  std::vector<std::vector<A_RESIDUAL> > residuals(SystMax);
  std::vector<std::vector<A_MODELED_MEASUREMENT> > mod(SystMax);
  A_VECTOR3D point;
  std::vector<std::vector<A_MEASUREMENT> > smoothMeas, smoothMeasBeidou;
  
  for (int isyst=0; isyst<SystMax;isyst++)
  {
    const int maxSatSyst = getMaxSatSyst(isyst);
    residuals[isyst] = std::vector<A_RESIDUAL>(maxSatSyst);
    mod[isyst] = std::vector<A_MODELED_MEASUREMENT>(maxSatSyst);
  }
  
  if(_settings._reset)
  {
    if (!((int)date.getSecond() % _settings._reset))
      _bInit = 0;
  }

  if (_bInit == 0)
  {
    _filter.initState(rover);
    _state.setX(0.0);
    _state.setY(0.0);
    _state.setZ(0.0);
    _bInit = 1;    
  }

  // State propagation
  if (_filter.propagateState(measurements, biases))
    {
    _bInit = 0;
    return 0;
    }

  // smoothing
  smoothMeasurements(date, measurements, smoothMeas);

  // Coarse point
  coarsePoint(smoothMeas, Positions, date, point);
  
  if (!point.getX() && !point.getY() && !point.getZ())
  {
    return 0;
  }

  // Current state initialization
  if (!_state.getX() && !_state.getY() && !_state.getZ())
  {
    if (!roverAPrioriPosition.getX() && !roverAPrioriPosition.getY() && !roverAPrioriPosition.getZ())
      _state = point;
    else      
      _state = roverAPrioriPosition;
  }

  // Model
  A_VECTOR3D MRec = point ;
#ifdef RESIDUALS
  MRec = roverAPrioriPosition;
#endif
  // Model
  computeModel(date, MRec, measurements, Positions, mod);

  // Beidou code orrection
  beidouCorrection(smoothMeas, smoothMeasBeidou, mod);

  // Residuals computation
  computeResiduals(residuals, biases, Positions, smoothMeasBeidou, mod);

  // State correction 
  if (_filter.correctState(date, _state, MRec, residuals, tropo, recProduct, measurementProduct))
    {
    _bInit = 0;
    return 0;
    }
  recProduct._integrity = integrity(residuals, recProduct, measurementProduct);
  
#ifdef RESIDUALS
  // Output residuals and PD for debug
  for (int isyst=0; isyst<SystMax;isyst++)
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
      {
      if (residuals[isyst][isat]._resCode[F1])
        {
        char c[SystMax]={'G','R','E','C'};
  	fprintf(stderr, "%c%02d %d %6.0lf ", c[isyst], isat+1, date.getDay(), date.getSecond());
  	for(int f=F1;f<FMax;f++)
  	  fprintf(stderr, "%lf ", residuals[isyst][isat]._resCode[f]);
  	for(int f=F1;f<FMax;f++)
  	  fprintf(stderr, "%lf ", residuals[isyst][isat]._resPhase[f]);
        A_VECTOR3D _pos=residuals[isyst][isat]._pos;
  	fprintf(stderr, "%lf %lf %lf %lf %lf %lf %lf\n", _pos.getX(), _pos.getY(), _pos.getZ(),
        	residuals[isyst][isat]._mapping, recProduct._tropoV, recProduct._tropo, residuals[isyst][isat]._iono);
        }
      }
#endif

  _state = A_VECTOR3D(recProduct._x, recProduct._y, recProduct._z);

  return 1;
}
