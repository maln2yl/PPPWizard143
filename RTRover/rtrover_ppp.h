/************************************************************
Nom ......... : rtrover_ppp.h
Role ........ : ppp definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef _RTROVER_PPP
#define _RTROVER_PPP

#include <vector>

#include "rtrover_map.h"
#include "rtrover_date.h"
#include "rtrover_model.h"
#include "rtrover_filter.h"
#include "rtrover_smoothed_measurement.h"
#include "rtrover_utility.h"

class A_MODE
{
  public:
    A_MODE() : _code(0), _phase(0), _freq(std::vector<int>(FMax,0)), _ambig(0) {}
    int _code, _phase;
    std::vector<int> _freq;
    int _ambig;
};

class A_PPP_SETTING
{
  public :
    A_PPP_SETTING() : _dt(0),_thrMap(0),_smooth(0),_reset(0),_map(),_mode() {}
    double _dt;
    double _thrMap;
    double _smooth;
    int _reset;
    A_MAP _map;
    A_MODE _mode;
};

// PPP
class A_PPP
{
  public :
    static A_PPP_SETTING _settings;
    A_PPP();
    std::vector<std::vector<A_MODEL> >& getModel() { return _model; }
    int computePPP(const std::vector<std::vector<A_MEASUREMENT> >& Mesures,
             const std::vector<std::vector<A_BIAS> >& Biais,
             const std::vector<std::vector<A_SSR_PARAMETER> >& Positions,
             const double tropo,
             const A_DATE &Date,
             A_RECEIVER_PRODUCT &stationProduct,
             std::vector<std::vector<A_MEASUREMENT_PRODUCT> >& measurementProduct,
             const A_VECTOR3D& roverAPrioriPosition, const int rover);

  private :
    std::vector<A_SMOOTHED_MEASUREMENT> _smoothedMeas;
    A_FILTER _filter;
    std::vector<std::vector<A_MODEL> > _model;
    A_VECTOR3D _state;
    int _bInit;
    
    void computeModel(const A_DATE &date, const A_VECTOR3D &MRec,
                      const std::vector<std::vector<A_MEASUREMENT> >& measurements,
		      const std::vector<std::vector<A_SSR_PARAMETER> >& positions,
                      std::vector<std::vector<A_MODELED_MEASUREMENT> >& mod);
    void beidouCorrection(const std::vector<std::vector<A_MEASUREMENT> >& meas,
                             std::vector<std::vector<A_MEASUREMENT> >& measBeidou,
			     const std::vector<std::vector<A_MODELED_MEASUREMENT> >& mod);
    void coarsePoint(const std::vector<std::vector<A_MEASUREMENT> >& measurements,
		     const std::vector<std::vector<A_SSR_PARAMETER> >& positions,
		     const A_DATE &date,
		     A_VECTOR3D& point);
    void smoothMeasurements(const A_DATE &date, const std::vector<std::vector<A_MEASUREMENT> >& raw_meas,
                            std::vector<std::vector<A_MEASUREMENT> >& smooth_meas);
    void computeResiduals(std::vector<std::vector<A_RESIDUAL> >& residuals,
             const std::vector<std::vector<A_BIAS> >& bias,
             const std::vector<std::vector<A_SSR_PARAMETER> >& positions,
             const std::vector<std::vector<A_MEASUREMENT> >& mes,
             const std::vector<std::vector<A_MODELED_MEASUREMENT> >& mod);



    double integrity(const std::vector<std::vector<A_RESIDUAL> >& residuals,
		     const A_RECEIVER_PRODUCT &recProduct,
		     const std::vector<std::vector<A_MEASUREMENT_PRODUCT> > &measProduct);
};

#endif
