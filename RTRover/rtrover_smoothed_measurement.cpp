/************************************************************
Nom ......... : rtrover_smoothed_measurement.cpp
Role ........ : smoothed measurement definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#include "rtrover_smoothed_measurement.h"
#include "rtrover_utility.h"

//////////////////////////////////////////////////////////////////////
A_SMOOTHED_MEASUREMENT::A_SMOOTHED_MEASUREMENT()
{
  _dop = std::vector<std::vector<double> >(SystMax);
  _Pl = std::vector<std::vector<double> >(SystMax);
  for (int isyst=0; isyst<SystMax;isyst++)
  {
    _dop[isyst] = std::vector<double>(getMaxSatSyst(isyst),0.0);
    _Pl[isyst] = std::vector<double>(getMaxSatSyst(isyst),0.0);
  }
}

//////////////////////////////////////////////////////////////////////
void A_SMOOTHED_MEASUREMENT::smoothMeasurements(const A_DATE &date, const std::vector<std::vector<A_MEASUREMENT> >& mes, double alpha,int freq)
{
  double dt= date - _dateDop;
  for (int isyst=0; isyst<SystMax;isyst++)
  {   
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      const A_MEASUREMENT &mesSat = mes[isyst][isat];
      double &dopSat = _dop[isyst][isat];
      double &PlSat = _Pl[isyst][isat];
      double P = 0.0;
      double D = 0.0;

      P=mesSat._code[freq];
      D=mesSat._doppler[freq];
      
      if (dopSat && D &&P && PlSat && (dt < 10.0))
        PlSat = (1.0 - alpha) * P + alpha * (PlSat + dt * (dopSat + D)/2.0);
      else
        PlSat = P;
      
      dopSat = D;     
    }
  }
  _dateDop = date;
}
