/************************************************************
Nom ......... : rtrover_kalman_filter.h
Role ........ : kalman filter definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.0 9/30/2014
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef _RTROVER_KALMAN_FILTER
#define _RTROVER_KALMAN_FILTER

#include <vector>

#include "rtrover_filter_modeled_measurement.h"

class A_KALMAN_FILTER
{
  public:
    A_KALMAN_FILTER() : _iN(0), _iT(0), _trU(), _trD(), _trS() {}
    int initialize(const int ie_Size);
    int updateDiagonalCovariance(const int ie_Index, const double re_Value);
    int recoverDiagonalCovariance(const int ie_Index, double &prs_Value);
    int addDiagonalCovariance(const int ie_Index, const double re_Value);
    int newMeasurementsElim(std::vector<A_FILTER_MODELED_MEASUREMENT> &poe_Measurements,
                            const double re_ResMax, std::vector<double> &pos_S);
                            int newMeasurementsElimComb(const int param_RejetMax,
                            std::vector<A_FILTER_MODELED_MEASUREMENT> &poe_Measurements, 
			    const double re_ResMax, std::vector<double> &pos_S);
    
  private:
    int     _iN, _iT;
    std::vector<double> _trU;
    std::vector<double> _trD;
    std::vector<double> _trS;

    int initializeSolution();
    int newMeasurement(const A_FILTER_MODELED_MEASUREMENT &poe_Measurement, const double re_Sign);
    int solution(std::vector<double> &pos_S);
    int newMeasurementsElimNumber( std::vector<A_FILTER_MODELED_MEASUREMENT> &poe_Measurements, 
                                   const double re_ResMax, const int ie_Number, std::vector<double> &pos_S);
};

#endif /* _RTROVER_KALMAN_FILTER */
