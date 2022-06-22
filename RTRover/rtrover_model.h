/************************************************************
Nom ......... : rtrover_model.h
Role ........ : model definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef _RTROVER_MODEL
#define _RTROVER_MODEL

#include "rtrover_vector3D.h"
#include "rtrover_date.h"
#include "rtrover_utility.h"

class A_SATELLITE_PHASE_MAP;

// Model
class A_MODELED_MEASUREMENT
{ 
  public :
    A_MODELED_MEASUREMENT():_modCode(std::vector<double>(FMax,0.0)), _modPhase(std::vector<double>(FMax,0.0)),_mapping(0.0),_elevation(0.0),_pos() {}
    std::vector<double>  _modCode, _modPhase; 
    double _mapping;
    double _elevation;
    A_VECTOR3D _pos;
};

class A_MODEL
{
  public:
    A_MODEL() : _wu(0) {}
    A_MODELED_MEASUREMENT model(const A_DATE &date,
                  const int isyst,
                  const A_VECTOR3D &posRec,
                  const A_MEASUREMENT &mes,
                  const A_SSR_PARAMETER &positions,
		  const A_SATELLITE_PHASE_MAP &satellites);
    static void geodesic(const double x, const double y, const double z, double *lambda, double *phi, double *height);
    static double tropDelay(const double cosPhi, const double h);
  
  private:
    double _wu; /* Windup state */
};

#endif
