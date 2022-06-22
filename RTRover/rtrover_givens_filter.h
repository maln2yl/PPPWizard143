/************************************************************
Nom ......... : rtrover_givens_filter.h
Role ........ : givens filter definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.0 9/30/2014
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef _RTROVER_GIVENS_FILTER
#define _RTROVER_GIVENS_FILTER

#include <vector>

#include "rtrover_filter_modeled_measurement.h"

class A_GIVENS_FILTER
{
  public:
    A_GIVENS_FILTER(const int size);
    void newMeasurement(const A_FILTER_MODELED_MEASUREMENT &y);
    void solution(double *sol, double& var);
    void variance(double *sigma);
  
  private:
    int _n;
    int _nb;
    std::vector<double> _d;
    std::vector<double> _r;
    double _sumRes;
};

#endif
