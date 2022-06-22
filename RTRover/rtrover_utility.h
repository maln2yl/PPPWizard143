/************************************************************
Nom ......... : rtrover_utility.cpp
Role ........ : common objects and fonctions
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef _RTROVER_UTILITY
#define _RTROVER_UTILITY

#include <vector>

#include "rtrover_interface.h"

#ifndef CLIGHT
#define CLIGHT (299792458.0)
#endif

#define EPSDBL (1.e-10)

#define MAX_GPS 32
#define MAX_GLO 24
#define MAX_GAL 36
#define MAX_BDS 35
#define MAX_EPOCH 600

#ifdef RTROVER_DEBUG
#define OPEN_LOG(filename) fopen(filename, "a")
#else
#define OPEN_LOG(filename) NULL
#endif

#ifdef RTROVER_DEBUG
#define CLOSE_LOG(file) fclose(file)
#else
#define CLOSE_LOG(file)
#endif

#ifdef RTROVER_DEBUG
#define PRINT_LOG(file, ...) fprintf(file,__VA_ARGS__)
#else
#define PRINT_LOG(file, ...)
#endif

#ifdef RTROVER_DEBUG
#define FOR_LOG for
#else
#define FOR_LOG(x)
#endif

#ifdef RTROVER_DEBUG
#define FFLUSH_LOG(file) fflush(file)
#else
#define FFLUSH_LOG(file)
#endif

enum system { GPS, Glo, Gal, Bds, SystMax };
typedef enum {F1, F2, F5, F6, F7, FMax} Frequency;
typedef enum {GEO, IGSO, MEO} BdsSat;

class A_SSR_PARAMETER
{
  public :
 A_SSR_PARAMETER() : _Xsp3(0.0), _Ysp3(0.0), _Zsp3(0.0), _Hsp3(0.0), _VXsp3(0.), _VYsp3(0.), _VZsp3(0.), _iono(0.0), _yaw(0.), _typeIono(0) {}
    double _Xsp3, _Ysp3, _Zsp3, _Hsp3, _VXsp3, _VYsp3, _VZsp3, _iono, _yaw;
    int _typeIono;
};

class A_MEASUREMENT
{
  public :
    A_MEASUREMENT() : _code(std::vector<double>(FMax,0.0)), _phase(std::vector<double>(FMax,0.0)), _doppler(std::vector<double>(FMax,0.0)), _slot_typeBds(0) {}
    std::vector<double>  _code, _phase, _doppler;    
    int _slot_typeBds;
};

class A_BIAS
{
  public:
    A_BIAS() : _code(std::vector<double>(FMax,0.0)), _phase(std::vector<double>(FMax,0.0)),
    _n1i(false), _wli(0), _discontinuity(0) {}
    std::vector<double> _code, _phase;
    bool _n1i;
    int _wli;
    int _discontinuity;
};

double diffDate(const rtrover_time *t1, const rtrover_time *t2);

int convertSystem(const char system);

inline int getMaxSatSyst (int syst)
{
  if (syst == GPS) return MAX_GPS;
  if (syst == Glo) return MAX_GLO;
  if (syst == Gal) return MAX_GAL;
  if (syst == Bds) return MAX_BDS;
};

#endif
