/************************************************************
Nom ......... : rtrover_filter.h
Role ........ : filter definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.5 7/22/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef _RTROVER_FILTER
#define _RTROVER_FILTER

#include <cstdio>
#include <string>

#include "rtrover_vector3D.h"
#include "rtrover_kalman_filter.h"
#include "rtrover_utility.h"

typedef enum {EX, NW, N1} AmbiguityNature;
typedef enum {P1, P2, C5, C6, C7, L1, L2, L5, L6, L7, IONO, X, Y, Z, T, P, A4, Ex, JP, JA4, JEx, E} ParamNature;
typedef enum {pseudorange, carrierphase} Measure;

class A_RECEIVER_PRODUCT
{
  public :
    A_RECEIVER_PRODUCT() : _HstaGPS(0), _covHstaGPS(0), _HstaGlo(0), _covHstaGlo(0), _HstaGal(0), _covHstaGal(0),  _HstaBds(0), _covHstaBds(0),
                           _tropoV(0), _tropo(0), _covTropo(0), _x(0), _covX(0), _y(0), _covY(0), _z(0), _covZ(0), _nbMes(0), _nbMesBloNw(0), _nbMesBloN1(0) {}
    double _HstaGPS, _covHstaGPS, _HstaGlo, _covHstaGlo, _HstaGal, _covHstaGal, _HstaBds, _covHstaBds ;
    double _tropoV, _tropo, _covTropo;
    double _x,_covX;
    double _y, _covY;
    double _z, _covZ;
    int _nbMes, _nbMesBloEx, _nbMesBloNw, _nbMesBloN1;
    double _integrity;
};

class A_MEASUREMENT_PRODUCT
{
  public :
    A_MEASUREMENT_PRODUCT() : _resCode(std::vector<double>(FMax,0.0)),
                              _resPhase(std::vector<double>(FMax,0.0)),
			      _elimCode(std::vector<int>(FMax,0)),
                              _elimPhase(std::vector<int>(FMax,0)),
			      _ne(0), _nw(0), _n1(0),
			      _ambEx(0), _ambA4(0), _ambP(0),
			      _covAmbEx(0), _covAmbA4(0), _covAmbP(0),
			      _iono(0), _covIono(0) {}
    std::vector<double> _resCode, _resPhase;
    std::vector<int> _elimCode, _elimPhase;
    double _ne, _nw, _n1;
    double _ambEx, _ambA4, _ambP, _covAmbEx, _covAmbA4, _covAmbP;
    double _iono, _covIono;
};

class A_RESIDUAL
{ 
  public :
    A_RESIDUAL() : _resCode(std::vector<double>(FMax,0.0)), _resPhase(std::vector<double>(FMax,0.0)),
                   _slot_typeBds(0), _iono(0.0), _typeIono(0), _mapping(0.0), _pos(),
                  _n1i(false), _wli(0), _discontinuity(0) {}
    std::vector<double> _resCode, _resPhase;
    int _slot_typeBds;
    double _iono;
    int _typeIono;
    double _mapping;
    A_VECTOR3D _pos;
    bool _n1i;
    int _wli;
    bool _ewli;
    int _discontinuity;
};

class A_MEASUREMENT_SETTING
{
  public:
    A_MEASUREMENT_SETTING() : _use(0),_sigCode(std::vector<double>(3,0.0)),_sigPhase(std::vector<double>(3,0.0)) {}
    std::vector<double> _sigCode, _sigPhase;
    int _use;
};

class A_FILTER_SETTING
{
  public :    
    A_FILTER_SETTING() : _nbMin(0),_sigIniTro(0.0),_sigModTro(0.0),_sigIniBiasClk(0.0),_sigModBiasClk(0.0),
                         _nbSatFixAmb(4), _thrAmb(0.0),_sigMeasIono(std::vector<double>(3,0.0)),_thrMeasIono(0.0),_sigMeasTropo(0.0),_thrMeasTropo(0.0),
			 _gps(),_glo(),_gal(),_bds(),
			 _thrMeasCode(0.0),_thrMeasPhase(0.0),_sigIniIono(0.0),_sigModIono(0.0),
			 _sigIniPos(0.0),_sigModPos(0.0),_maxElim(0),_dtMax(0),_raim(0),
			 _fixN1(true), _fixNw(true), _fixEx(true) {}
    int _nbMin;
    double _sigIniTro;
    double _sigModTro;
    double _sigIniBiasClk;
    double _sigModBiasClk;
    int _nbSatFixAmb;
    double _thrAmb;
    std::vector<double> _sigMeasIono;
    double _thrMeasIono;
    double _sigMeasTropo;
    double _thrMeasTropo;
    A_MEASUREMENT_SETTING _gps;
    A_MEASUREMENT_SETTING _glo;
    A_MEASUREMENT_SETTING _gal;
    A_MEASUREMENT_SETTING _bds;
    double _thrMeasCode;
    double _thrMeasPhase;
    double _sigIniIono;
    double _sigModIono;
    double _sigIniPos;
    double _sigModPos;
    int _maxElim;
    int _dtMax;
    bool _fixN1, _fixNw, _fixEx; 
    int _raim;   
};

class A_LISTPARAM
{
  public :
    A_LISTPARAM();
    void setParam();
    int getIndex (const ParamNature param, const int syst=-1, const int ipass=-1) const;
    int getSize (const ParamNature param, const int syst=-1) const;
    int getNbInc() const { return _iNbInc; }
  private :
    std::vector< std::vector<int> > _iDebHP, _iDebHl;
    int _iDebT, _iDebX, _iDebY, _iDebZ, _iDebP, _iDebA4, _iDebEx, _iDebJP, _iDebJA4, _iDebJEx, _iDebe;
    std::vector< std::vector<int> > _iNbHP,_iNbHl;
    int _iNbT, _iNbX, _iNbY, _iNbZ, _iNbP, _iNbA4, _iNbEx, _iNbJP, _iNbJA4, _iNbJEx, _iNbe, _iNbInc;

};

class A_PASS_MANAGEMENT
{
  public:    
    static const int MaxPassSyst=12;
    int beginPass(int isyst, int isat);
    void endPass(int isyst, int ipass);
    int numPass(int isyst, int isat);
    void init();
    A_PASS_MANAGEMENT();
  private:
    std::vector< std::vector<int> > tabPassToSat;
    std::vector< std::vector<int> > tabSatToPass;
};

class A_FILTER
{
  public :    
    class A_FILTERLOG
    {
      public :
      FILE * _log;
      A_FILTERLOG() { _log =  OPEN_LOG("rtrover_log.txt"); }
      ~A_FILTERLOG() { CLOSE_LOG(_log); }
    };

    static A_FILTER_SETTING _filterSetting;
    A_FILTER();
    void initState(const int rover);
    int correctState(const A_DATE &date, const A_VECTOR3D& state, const A_VECTOR3D& point,
                         const std::vector <std::vector <A_RESIDUAL> >& residual,
			 const double tropo,
                         A_RECEIVER_PRODUCT &recProduct,
                         std::vector <std::vector <A_MEASUREMENT_PRODUCT> >& measurementProduct);
    int propagateState(const std::vector<std::vector<A_MEASUREMENT> >& measurements, const std::vector<std::vector<A_BIAS> >& biases);

  private :
    A_KALMAN_FILTER _filter;
    std::vector<std::vector<double> > _n, _a4, _ex;
    std::vector<double> _sol;
    std::vector<std::vector<int> > _nbMeasPas, _nbMesPas0, _discontinuityValue;
    A_KALMAN_FILTER _filterInt;
    std::vector<double> _solInt;
    A_LISTPARAM _params;
    A_PASS_MANAGEMENT _pass;
    static A_FILTERLOG _filterLog;
    int _rover;

    void fixAmbiguities(const AmbiguityNature an, const std::vector<std::vector<A_RESIDUAL> >& residuals, bool init);
    void fixJumps(const AmbiguityNature an, const std::vector<std::vector<A_RESIDUAL> >& residuals,
                        const std::vector<std::vector<int> > &elimL, std::vector<std::vector<int> > &fixed, int isyst, bool init);
    void eliminateMeasurement(const std::string strMeasure, std::vector<std::vector<double> >& measurementsRes,
                              std::vector<std::vector<int> >& measurementsElim, 
                              const std::vector<std::vector<A_RESIDUAL> >& residual,
                              std::vector<std::vector<A_FILTER_MODELED_MEASUREMENT> >& measurements,
                              const std::vector<std::vector<int> >& satMes,
                              const double thr_meas, std::vector<int>& reject);
    void computeMeasurementCode(const Frequency f, const std::vector<std::vector<A_RESIDUAL> >& residual,
                              std::vector<std::vector<A_FILTER_MODELED_MEASUREMENT> >& measurements,
                              std::vector<std::vector<int> >& satMes);
    void computeMeasurementPhase(const Frequency f, const std::vector<std::vector<A_RESIDUAL> >& residual,
                              std::vector<std::vector<A_FILTER_MODELED_MEASUREMENT> >& measurements,
                              std::vector<std::vector<int> >& satMes);		      
    void computeMeasurementIONO(const std::vector<std::vector<A_RESIDUAL> >& residual,
                              std::vector<std::vector<A_FILTER_MODELED_MEASUREMENT> >& measurements,
                              std::vector<std::vector<int> >& satMes);
    void computeMeasurementTROPO(const double tropo,
                                 double *measurementsRes,
				 int *measurementsElim);
    void initStateVectorAmb(int isyst, int isat, int ipass, double value, bool jumps);
};


#endif

