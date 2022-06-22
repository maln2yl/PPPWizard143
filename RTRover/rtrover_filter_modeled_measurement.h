/************************************************************
Nom ......... : rtrover_filter_modeled_measurement.h
Role ........ : filter modeled measurement definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.0 9/30/2014
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef _RTROVER_FILTER_MODELED_MEASUREMENT
#define _RTROVER_FILTER_MODELED_MEASUREMENT

#include <vector>

class A_FILTER_MODELED_MEASUREMENT
{
  public :
    A_FILTER_MODELED_MEASUREMENT() : _res(0),_w(0),_der() {}
    A_FILTER_MODELED_MEASUREMENT(const double res, const double w, const std::vector<double> der) : _res(res),_w(w),_der(der) {} 
    double _res ; /* Measurement residual */
    double _w ; /* Measurement noise    */
    std::vector<double> _der; /* Partial derivatives  */
};

#endif
