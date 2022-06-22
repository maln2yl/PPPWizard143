/************************************************************
Nom ......... : rtrover_smoothed_measurement.h
Role ........ : smoothed measurement definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.2 1/30/2015
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef _RTROVER_SMOOTHED_MEASUREMENT
#define _RTROVER_SMOOTHED_MEASUREMENT

#include <vector>
#include <cstdlib>

#include "rtrover_date.h"
#include "rtrover_utility.h"


class A_SMOOTHED_MEASUREMENT
{
  public :
    A_SMOOTHED_MEASUREMENT();
    void smoothMeasurements(const A_DATE &date, const std::vector<std::vector<A_MEASUREMENT> >& mes, double alpha, int freq);
    std::vector<std::vector<double> >& getDop() { return _dop; }
    std::vector<std::vector<double> >& getPl() { return _Pl; }
  
  private :
    std::vector<std::vector<double> > _dop, _Pl;
    A_DATE _dateDop;
    
};

#endif
