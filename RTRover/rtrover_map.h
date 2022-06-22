/************************************************************
Nom ......... : rtrover_map.h
Role ........ : map definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef _RTROVER_MAP
#define _RTROVER_MAP

#include <vector>

#include "rtrover_vector3D.h"

class A_PHASE_MAP
{
  public:
    A_PHASE_MAP() : _zen1(0), _zen2(0), _dZen(0), _phasePattern() {}
    double _zen1, _zen2, _dZen;
    std::vector<double> _phasePattern;
};

// Map
class A_SATELLITE_PHASE_MAP
{
  public:
    A_SATELLITE_PHASE_MAP():_prn(0),_phaseMapL1x(0),_phaseMapL1y(0),_phaseMapL1z(0),_phaseMapL2x(0),_phaseMapL2y(0),_phaseMapL2z(0),_mapPhaseL1(),_mapPhaseL2() {}
    void phaseCenterVariation (const A_VECTOR3D &sta, const A_VECTOR3D &sat, double *x1, double *x2) const ;
    void clear();
    const int getPrn() const { return _prn; }
    const double getPhaseMapL1x() const { return _phaseMapL1x; }
    const double getPhaseMapL1y() const { return _phaseMapL1y; }
    const double getPhaseMapL1z() const { return _phaseMapL1z; }
    const double getPhaseMapL2x() const { return _phaseMapL2x; }
    const double getPhaseMapL2y() const { return _phaseMapL2y; }
    const double getPhaseMapL2z() const { return _phaseMapL2z; }
    A_PHASE_MAP& getMapPhaseL1() { return _mapPhaseL1; }
    A_PHASE_MAP& getMapPhaseL2() { return _mapPhaseL2; }
    void setPrn (const int prn) { _prn=prn; }
    void setPhaseMapL1x (const double phaseMapL1x) { _phaseMapL1x=phaseMapL1x; }
    void setPhaseMapL1y (const double phaseMapL1y) { _phaseMapL1y=phaseMapL1y; }
    void setPhaseMapL1z (const double phaseMapL1z) { _phaseMapL1z=phaseMapL1z; }
    void setPhaseMapL2x (const double phaseMapL2x) { _phaseMapL2x=phaseMapL2x; }
    void setPhaseMapL2y (const double phaseMapL2y) { _phaseMapL2y=phaseMapL2y; }
    void setPhaseMapL2z (const double phaseMapL2z) { _phaseMapL2z=phaseMapL2z; }    
    
  private:
    int _prn;
    double _phaseMapL1x, _phaseMapL1y, _phaseMapL1z;
    double _phaseMapL2x, _phaseMapL2y, _phaseMapL2z;
    A_PHASE_MAP _mapPhaseL1;
    A_PHASE_MAP _mapPhaseL2;
};

class A_MAP
{
  public:
    A_MAP();
    std::vector<std::vector<A_SATELLITE_PHASE_MAP> >& getSatellites() { return _satellites; }
    void loadATX(const char *fileName); 
    
  private:
    std::vector<std::vector<A_SATELLITE_PHASE_MAP> > _satellites;
};

#endif
