/************************************************************
Nom ......... : rtrover_broadcast.h
Role ........ : broadcast message definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef RTROVER_BROADCAST_H
#define RTROVER_BROADCAST_H

//#include "rtrover_interface.h"
#include "rtrover.h"
#include "rtrover_vector3D.h"

class A_BROADCAST
{
  public :
    virtual ~A_BROADCAST() {}
    virtual const int getIODE() const = 0;
    virtual const int getIODC() const = 0;
    virtual const int getSatelliteNumber() const = 0;
    virtual void computePosition(const rtrover_time *tt, A_VECTOR3D& xc, double& tc, A_VECTOR3D& vv, const int correction) = 0;
	

};

class A_BROADCAST_GPS : public A_BROADCAST
{
  public :
    A_BROADCAST_GPS(rtrover_ephGPS _eph) : _eph(_eph) {}
    A_BROADCAST_GPS();
    ~A_BROADCAST_GPS(){}
    const int getIODE() const { return _eph._IODE; }
    const int getIODC() const { return _eph._IODC; }
    const int getSatelliteNumber() const {return _eph._satellite._number;}
    void computePosition(const rtrover_time *tt, A_VECTOR3D& xc, double& tc, A_VECTOR3D& vv, const int correction);
    rtrover_ephGPS& getEph() { return _eph; }
	
  private :
    rtrover_ephGPS _eph;

};

class A_BROADCAST_Glo : public A_BROADCAST
{
  public :
    A_BROADCAST_Glo(rtrover_ephGlo _eph) : _eph(_eph) {}
    A_BROADCAST_Glo();
    ~A_BROADCAST_Glo(){}
    const int getIODE() const;	
    const int getIODC() const { return getIODE(); }
    const int getSatelliteNumber() const {return _eph._satellite._number;}
    void computePosition(const rtrover_time *tt, A_VECTOR3D& xc, double& tc, A_VECTOR3D& vv, const int correction);
    rtrover_ephGlo& getEph() { return _eph; }
	
  private :
    rtrover_ephGlo _eph;
    void gloDeriv(const double xv[6], const double acc[3], double va[6]);
    void rungeKutta4(const double yi[6], const double dx, const double acc[3], double yf[6]);
};

class A_BROADCAST_Gal : public A_BROADCAST
{
  public :
    A_BROADCAST_Gal(rtrover_ephGal _eph) : _eph(_eph) {}
    A_BROADCAST_Gal();
    ~A_BROADCAST_Gal(){}
    const int getIODE() const { return _eph._IODNav; }
    const int getIODC() const { return _eph._IODNav; }
    const int getSatelliteNumber() const {return _eph._satellite._number;}
    void computePosition(const rtrover_time *tt, A_VECTOR3D& xc, double& tc, A_VECTOR3D& vv, const int correction);
    rtrover_ephGal& getEph() { return _eph; }
	
  private :
    rtrover_ephGal _eph;

};

// RTCM3 BDS EPH encoding
//////////////////////////////////////////////////////////
// BNC ENCODING
#define BDSTOINT(type, value) static_cast<type>(round(value))
#define BDSADDBITS(a, b) {bitbuffer = (bitbuffer<<(a)) \
                       |(BDSTOINT(long long,b)&((1ULL<<a)-1)); \
                       numbits += (a); \
                       while(numbits >= 8) { \
                       buffer[size++] = bitbuffer>>(numbits-8);numbits -= 8;}}
#define BDSADDBITSFLOAT(a,b,c) {long long i = BDSTOINT(long long,(b)/(c)); \
                             BDSADDBITS(a,i)};

class A_BROADCAST_Bds : public A_BROADCAST
{
  public :
    A_BROADCAST_Bds(rtrover_ephBds _eph) : _eph(_eph) {}
    A_BROADCAST_Bds();
    ~A_BROADCAST_Bds(){}
    const int getIODE() const { return _eph._IODE; }
    const int getIODC() const;
    const int getSatelliteNumber() const {return _eph._satellite._number;}
    void computePosition(const rtrover_time *tt, A_VECTOR3D& xc, double& tc, A_VECTOR3D& vv, const int correction);
    rtrover_ephBds& getEph() { return _eph; }
	
  private :
    rtrover_ephBds _eph;

};

#endif
