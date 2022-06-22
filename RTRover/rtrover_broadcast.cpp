/************************************************************
Nom ......... : rtrover_broadcast.cpp
Role ........ : broadcast message definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#include <cmath>
#include <cstdio>

#include "rtrover_broadcast.h"
#include "rtrover_date.h"
#include "rtrover_utility.h"


//////////////////////////////////////////////////////////////////////
// Position GPS
A_BROADCAST_GPS::A_BROADCAST_GPS()
{
  _eph._TOC._mjd = 0;
  _eph._TOC._sec = 0.0; 
}

//////////////////////////////////////////////////////////////////////
void A_BROADCAST_GPS::computePosition(const rtrover_time *tt, A_VECTOR3D& xc, double& tc, A_VECTOR3D& vv, const int correction)
{
  static const double omegaEarth  = 7292115.1467e-11;
  static const double gmWGS       = 398.6005e12;
  static const double validityBrdc = 14400;

  xc=A_VECTOR3D();
  tc=0.0;
  vv=A_VECTOR3D();

  double a0 = _eph._sqrt_A * _eph._sqrt_A;
  if (a0 == 0) {
    return;
  }

  double n0 = sqrt(gmWGS/(a0*a0*a0));
  double tk = diffDate(tt, &_eph._TOE);
  if ( correction || (!correction && (tk<validityBrdc/2) && (tk>-validityBrdc/2)))
  {
    double n  = n0 + _eph._Delta_n;
    double M  = _eph._M0 + n*tk;
    double E  = M;
    double E_last;
    do {
      E_last = E;
      E = M + _eph._e*sin(E);
    } while ( fabs(E-E_last)*a0 > 0.001 );
    double v      = 2.0*atan( sqrt( (1.0 + _eph._e)/(1.0 - _eph._e) )*tan( E/2 ) );
    double u0     = v + _eph._omega;
    double sin2u0 = sin(2*u0);
    double cos2u0 = cos(2*u0);
    double r      = a0*(1 - _eph._e*cos(E)) + _eph._Crc*cos2u0 + _eph._Crs*sin2u0;
    double i      = _eph._i0 + _eph._IDOT*tk + _eph._Cic*cos2u0 + _eph._Cis*sin2u0;
    double u      = u0 + _eph._Cuc*cos2u0 + _eph._Cus*sin2u0;
    double xp     = r*cos(u);
    double yp     = r*sin(u);

    double deltat = (double)_eph._TOE._mjd-44244.0+_eph._TOE._sec/86400.0;
    int nweek = (int) floor(deltat/7.0);
    double TOEsec = (deltat - (nweek)*7.0)*86400.0;

    double OM     = _eph._OMEGA0 + (_eph._OMEGADOT - omegaEarth)*tk - omegaEarth*TOEsec;

    double sinom = sin(OM);
    double cosom = cos(OM);
    double sini  = sin(i);
    double cosi  = cos(i);
    xc.setX( xp * cosom - yp * cosi * sinom);
    xc.setY( xp * sinom + yp * cosi * cosom);
    xc.setZ( yp * sini);                 

    double diffT = diffDate(tt, &_eph._TOC);
    tc = _eph._clock_bias + _eph._clock_drift*diffT + _eph._clock_driftrate*diffT*diffT;

    // Velocity
    // --------
    double tanv2 = tan(v/2);
    double dEdM  = 1 / (1 - _eph._e*cos(E));
    double dotv  = sqrt((1.0 + _eph._e)/(1.0 - _eph._e)) / cos(E/2)/cos(E/2) / (1 + tanv2*tanv2) 
        	 * dEdM * n;
    double dotu  = dotv + (-_eph._Cuc*sin2u0 + _eph._Cus*cos2u0)*2*dotv;
    double dotom = _eph._OMEGADOT - omegaEarth;
    double doti  = _eph._IDOT + (-_eph._Cic*sin2u0 + _eph._Cis*cos2u0)*2*dotv;
    double dotr  = a0 * _eph._e*sin(E) * dEdM * n 
                  + (-_eph._Crc*sin2u0 + _eph._Crs*cos2u0)*2*dotv;
    double dotx  = dotr*cos(u) - r*sin(u)*dotu;
    double doty  = dotr*sin(u) + r*cos(u)*dotu;

    vv.setX( cosom   *dotx  - cosi*sinom   *doty      // dX / dr
             - xp*sinom*dotom - yp*cosi*cosom*dotom   // dX / dOMEGA
             + yp*sini*sinom*doti );                  // dX / di

    vv.setY( sinom   *dotx  + cosi*cosom   *doty
            + xp*cosom*dotom - yp*cosi*sinom*dotom
            - yp*sini*cosom*doti);

    vv.setZ( sini    *doty  + yp*cosi      *doti);

    // Relativistic Correction
    // -----------------------
    tc -= 2.0 * ( xc * vv ) / CLIGHT / CLIGHT;
  }
  else
  {
    tc = 0.0;
    fprintf(stderr,"Eph GPS too old ! : tk = %f\n", tk);
  }
}


//////////////////////////////////////////////////////////////////////
const int A_BROADCAST_Glo::getIODE() const
{
  A_DATE date(_eph._timeUTC._mjd-33282 , _eph._timeUTC._sec);
  date = date + (3*3600.0);
  return (int)(date.getSecond()/900);
}

A_BROADCAST_Glo::A_BROADCAST_Glo()
{
  _eph._timeUTC._mjd = 0;
  _eph._timeUTC._sec = 0.0; 
}

// Derivative of the state vector using a simple force model
//////////////////////////////////////////////////////////////////////
void A_BROADCAST_Glo::gloDeriv(const double xv[6], const double acc[3], double va[6])
{
  // State vector components
  // -----------------------
  double rr[3], vv[3];
  rr[0]=xv[0];rr[1]=xv[1];rr[2]=xv[2];
  vv[0]=xv[3];vv[1]=xv[4];vv[2]=xv[5];

  // Acceleration 
  // ------------
  static const double GM    = 398.60044e12;
  static const double AE    = 6378136.0;
  static const double OMEGA = 7292115.e-11;
  static const double C20   = -1082.6257e-6;

  double rho = sqrt(rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2]);
  double t1  = -GM/(rho*rho*rho);
  double t2  = 3.0/2.0 * C20 * (GM*AE*AE) / (rho*rho*rho*rho*rho);
  double t3  = OMEGA * OMEGA;
  double t4  = 2.0 * OMEGA;
  double z2  = rr[2] * rr[2];

  // Vector of derivatives
  // ---------------------
  va[0] = vv[0];
  va[1] = vv[1];
  va[2] = vv[2];
  va[3] = (t1 + t2*(1.0-5.0*z2/(rho*rho)) + t3) * rr[0] + t4*vv[1] + acc[0]; 
  va[4] = (t1 + t2*(1.0-5.0*z2/(rho*rho)) + t3) * rr[1] - t4*vv[0] + acc[1]; 
  va[5] = (t1 + t2*(3.0-5.0*z2/(rho*rho))     ) * rr[2]            + acc[2];
}


//////////////////////////////////////////////////////////////////////
// Fourth order Runge-Kutta numerical integrator for ODEs
void A_BROADCAST_Glo::rungeKutta4(const double yi[6], const double dx, const double acc[3], double yf[6])
{
  int i;
  double k1[6], k2[6], k3[6], k4[6], y1[6], zi[6];

  //  ColumnVector k1 = der(xi       , yi       , acc) * dx;
  for (i=0;i<6;i++) y1[i]= yi[i]; gloDeriv(y1, acc, zi); for (i=0;i<6;i++) k1[i] = dx*zi[i];

  //  ColumnVector k2 = der(xi+dx/2.0, yi+k1/2.0, acc) * dx;
  for (i=0;i<6;i++) y1[i]= yi[i]+k1[i]/2.0; gloDeriv(y1, acc, zi); for (i=0;i<6;i++) k2[i] = dx*zi[i];

  //  ColumnVector k3 = der(xi+dx/2.0, yi+k2/2.0, acc) * dx;
  for (i=0;i<6;i++) y1[i]= yi[i]+k2[i]/2.0; gloDeriv(y1, acc, zi); for (i=0;i<6;i++) k3[i] = dx*zi[i];

  // ColumnVector k4 = der(xi+dx    , yi+k3    , acc) * dx;
  for (i=0;i<6;i++) y1[i]= yi[i]+k3[i]; gloDeriv(y1, acc, zi); for (i=0;i<6;i++) k4[i] = dx*zi[i];

  //  ColumnVector yf = yi + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;
  for (i=0;i<6;i++) yf[i] = yi[i] + k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0;

}

//////////////////////////////////////////////////////////////////////
// Position Glonass
void A_BROADCAST_Glo::computePosition(const rtrover_time *tt, A_VECTOR3D& xc, double& tc, A_VECTOR3D& vv, const int correction)
{
  static const double nominalStep = 10.0;
  static const double validityBrdc = 15*60.0;

  xc=A_VECTOR3D();
  tc=0.0;
  vv=A_VECTOR3D();

  double dtPos = diffDate(tt, &_eph._timeUTC)-_eph._gps_utc;
  if ( correction || (!correction && (dtPos<validityBrdc) && (dtPos>-validityBrdc)))
  {
    int nSteps  = int(fabs(dtPos) / nominalStep) + 1;
    double step = dtPos / nSteps;

    double acc[3];
    acc[0] = _eph._x_acceleration * 1.e3;
    acc[1] = _eph._y_acceleration * 1.e3;
    acc[2] = _eph._z_acceleration * 1.e3;

    double yi[6];
    yi[0] = _eph._x_pos * 1.e3; 
    yi[1] = _eph._y_pos * 1.e3; 
    yi[2] = _eph._z_pos * 1.e3; 
    yi[3] = _eph._x_velocity * 1.e3; 
    yi[4] = _eph._y_velocity * 1.e3; 
    yi[5] = _eph._z_velocity * 1.e3; 
    //fprintf(stderr, "deb = %lf %lf %lf %lf %lf %lf\n", yi[0], yi[1], yi[2], yi[3], yi[4], yi[5]);

    for (int ii = 1; ii <= nSteps; ii++)
    { 
      double yf[6];
	  rungeKutta4(yi, step, acc, yf);
	  for (int i=0;i<6;i++)
	    yi[i] = yf[i];
    }

    xc.setX(yi[0]);
    xc.setY(yi[1]);
    xc.setZ(yi[2]);
    vv.setX(yi[3]);
    vv.setY(yi[4]);
    vv.setZ(yi[5]);
  //fprintf(stderr, "dtPos = %lf %lf %lf %lf %lf %lf %lf %lf\n", dtPos, xc.getX(), xc.getY(), xc.getZ(), tc, vv.getX(), vv.getY(), vv.getZ());

    double dtClk = diffDate(tt, &_eph._timeUTC)-_eph._gps_utc;
    tc = -_eph._tau + _eph._gamma * dtClk;
    //fprintf(stderr, "dtClk = %15.3lf %15.3lf %lg %15.3lf\n", dtClk, _eph._tau, _eph._gamma, tc*CLIGHT);
  }
  else
  {
    tc = 0.0;
    fprintf(stderr,"Eph Glo too old ! : dtPos = %f\n",dtPos);
  }
}

//////////////////////////////////////////////////////////////////////
// Position Galileo
A_BROADCAST_Gal::A_BROADCAST_Gal()
{
  _eph._TOC._mjd = 0;
  _eph._TOC._sec = 0.0; 
}

//////////////////////////////////////////////////////////////////////
void A_BROADCAST_Gal::computePosition(const rtrover_time *tt, A_VECTOR3D& xc, double& tc, A_VECTOR3D& vv, const int correction)
{
  static const double omegaEarth = 7292115.1467e-11;
  static const double gmWGS      = 398.6004418e12;
  static const double validityBrdc = 14400;

  xc=A_VECTOR3D();
  tc=0.0;
  vv=A_VECTOR3D();
  /*std::cout << "computePosition Gal " << tt->_mjd << " " << tt->_sec << std::endl;
  std::cout << _eph._TOC._mjd << " " << _eph._TOC._sec << std::endl;
  std::cout << _eph._TOE._mjd << " " << _eph._TOE._sec << std::endl;
  std::cout << _eph._IODNav << std::endl;
  std::cout << _eph._Crs << " " << _eph._Delta_n << " " << _eph._M0 << std::endl;
  std::cout << _eph._Cuc << " " << _eph._e << " " << _eph._Cus << std::endl;
  std::cout << _eph._sqrt_A << " " << _eph._Cic << " " << _eph._OMEGA0 << std::endl;
  std::cout << _eph._Cis << " " << _eph._i0 << " " << _eph._Crc << std::endl;
  std::cout << _eph._omega << " " << _eph._OMEGADOT << " " << _eph._IDOT << std::endl;*/

  double a0 = _eph._sqrt_A * _eph._sqrt_A;
  if (a0 == 0) {
    return;
  }

  double n0 = sqrt(gmWGS/(a0*a0*a0));
  double tk = diffDate(tt, &_eph._TOE);

  if ( correction || (!correction && (tk<validityBrdc) && (tk>0)))
  {
    double n  = n0 + _eph._Delta_n;
    double M  = _eph._M0 + n*tk;
    double E  = M;
    double E_last;
    do {
      E_last = E;
      E = M + _eph._e*sin(E);
    } while ( fabs(E-E_last)*a0 > 0.001 );
    double v      = 2.0*atan( sqrt( (1.0 + _eph._e)/(1.0 - _eph._e) )*tan( E/2 ) );
    double u0     = v + _eph._omega;
    double sin2u0 = sin(2*u0);
    double cos2u0 = cos(2*u0);
    double r      = a0*(1 - _eph._e*cos(E)) + _eph._Crc*cos2u0 + _eph._Crs*sin2u0;
    double i      = _eph._i0 + _eph._IDOT*tk + _eph._Cic*cos2u0 + _eph._Cis*sin2u0;
    double u      = u0 + _eph._Cuc*cos2u0 + _eph._Cus*sin2u0;
    double xp     = r*cos(u);
    double yp     = r*sin(u);

    double deltat = (double)_eph._TOE._mjd-51412.0+_eph._TOE._sec/86400.0;
    int nweek = (int) floor(deltat/7.0);
    double TOEsec = (deltat - (nweek)*7.0)*86400.0;

    double OM     = _eph._OMEGA0 + (_eph._OMEGADOT - omegaEarth)*tk - omegaEarth*TOEsec;

    double sinom = sin(OM);
    double cosom = cos(OM);
    double sini  = sin(i);
    double cosi  = cos(i);
    xc.setX( xp * cosom - yp * cosi * sinom);
    xc.setY( xp * sinom + yp * cosi * cosom);
    xc.setZ( yp * sini);                 

    double diffT = diffDate(tt, &_eph._TOC);
    tc = _eph._clock_bias + _eph._clock_drift*diffT + _eph._clock_driftrate*diffT*diffT;


    // Velocity
    // --------
    double tanv2 = tan(v/2);
    double dEdM  = 1 / (1 - _eph._e*cos(E));
    double dotv  = sqrt((1.0 + _eph._e)/(1.0 - _eph._e)) / cos(E/2)/cos(E/2) / (1 + tanv2*tanv2) 
        	 * dEdM * n;
    double dotu  = dotv + (-_eph._Cuc*sin2u0 + _eph._Cus*cos2u0)*2*dotv;
    double dotom = _eph._OMEGADOT - omegaEarth;
    double doti  = _eph._IDOT + (-_eph._Cic*sin2u0 + _eph._Cis*cos2u0)*2*dotv;
    double dotr  = a0 * _eph._e*sin(E) * dEdM * n 
                  + (-_eph._Crc*sin2u0 + _eph._Crs*cos2u0)*2*dotv;
    double dotx  = dotr*cos(u) - r*sin(u)*dotu;
    double doty  = dotr*sin(u) + r*cos(u)*dotu;

    vv.setX( cosom   *dotx  - cosi*sinom   *doty      // dX / dr
             - xp*sinom*dotom - yp*cosi*cosom*dotom   // dX / dOMEGA
             + yp*sini*sinom*doti );                  // dX / di

    vv.setY( sinom   *dotx  + cosi*cosom   *doty
            + xp*cosom*dotom - yp*cosi*sinom*dotom
            - yp*sini*cosom*doti);

    vv.setZ( sini    *doty  + yp*cosi      *doti);

    // Relativistic Correction
    // -----------------------
    //tc -= 2.0 * ( xc * vv ) / CLIGHT / CLIGHT;
    tc -= 4.442807309e-10 * _eph._e * sqrt(a0) *sin(E);
    //tc -= 4.442807633e-10 * _eph._e * sqrt(a0) *sin(E);
    //printf("EPH GAL brdc: %d %f %15.3lf %15.3lf %15.3lf %15.3lf\n",tt->_mjd, tt->_sec, xc.getX(), xc.getY(), xc.getZ(), CLIGHT*tc);	 
 
  }
  else
  {
    tc = 0.0;
    fprintf(stderr,"Eph Gal too old : tk = %f!\n",tk);
  }

}

//////////////////////////////////////////////////////////////////////
// Position BDS
A_BROADCAST_Bds::A_BROADCAST_Bds()
{
  _eph._TOC._mjd = 0;
  _eph._TOC._sec = 0.0; 
}

// Returns CRC24 (BNC fonction)
////////////////////////////////////////////////////////////////////////////
unsigned long CRC24(long size, const unsigned char *buf) {
  unsigned long crc = 0;
  int ii;
  while (size--) {
    crc ^= (*buf++) << (16);
    for(ii = 0; ii < 8; ii++) {
      crc <<= 1;
      if (crc & 0x1000000) {
        crc ^= 0x01864cfb;
      }
    }
  }
  return crc;
}
//////////////////////////////////////////////////////////////////////
const int A_BROADCAST_Bds::getIODC() const
{  
 /* unsigned char buffer[80];
  int size = 0;
  int numbits = 0;
  long long bitbuffer = 0;
  unsigned char *startbuffer = buffer;

  // BNC encoding
  BDSADDBITSFLOAT(14, _eph._IDOT, M_PI/static_cast<double>(1<<30)/static_cast<double>(1<<13))
  BDSADDBITSFLOAT(11, _eph._clock_driftrate, 1.0/static_cast<double>(1<<30)
      /static_cast<double>(1<<30)/static_cast<double>(1<<6))
  BDSADDBITSFLOAT(22, _eph._clock_drift, 1.0/static_cast<double>(1<<30)/static_cast<double>(1<<20))
  BDSADDBITSFLOAT(24, _eph._clock_bias, 1.0/static_cast<double>(1<<30)/static_cast<double>(1<<3))
  BDSADDBITSFLOAT(18, _eph._Crs, 1.0/static_cast<double>(1<<6))
  BDSADDBITSFLOAT(16, _eph._Delta_n, M_PI/static_cast<double>(1<<30)/static_cast<double>(1<<13))
  BDSADDBITSFLOAT(32, _eph._M0, M_PI/static_cast<double>(1<<30)/static_cast<double>(1<<1))
  BDSADDBITSFLOAT(18, _eph._Cuc, 1.0/static_cast<double>(1<<30)/static_cast<double>(1<<1))
  BDSADDBITSFLOAT(32, _eph._e, 1.0/static_cast<double>(1<<30)/static_cast<double>(1<<3))
  BDSADDBITSFLOAT(18, _eph._Cus, 1.0/static_cast<double>(1<<30)/static_cast<double>(1<<1))
  BDSADDBITSFLOAT(32, _eph._sqrt_A, 1.0/static_cast<double>(1<<19))
  BDSADDBITSFLOAT(18, _eph._Cic, 1.0/static_cast<double>(1<<30)/static_cast<double>(1<<1))
  BDSADDBITSFLOAT(32, _eph._OMEGA0, M_PI/static_cast<double>(1<<30)/static_cast<double>(1<<1))
  BDSADDBITSFLOAT(18, _eph._Cis, 1.0/static_cast<double>(1<<30)/static_cast<double>(1<<1))
  BDSADDBITSFLOAT(32, _eph._i0, M_PI/static_cast<double>(1<<30)/static_cast<double>(1<<1))
  BDSADDBITSFLOAT(18, _eph._Crc, 1.0/static_cast<double>(1<<6))
  BDSADDBITSFLOAT(32, _eph._omega, M_PI/static_cast<double>(1<<30)/static_cast<double>(1<<1))
  BDSADDBITSFLOAT(24, _eph._OMEGADOT, M_PI/static_cast<double>(1<<30)/static_cast<double>(1<<13))
  BDSADDBITS(5, 0)  // the last byte is filled by 0-bits to obtain a length of an integer multiple of 8

  return (int)CRC24(size, startbuffer); */
  return ((int)_eph._toes/720)%240; /* PATCH for the new IOD Beidou */
}

//////////////////////////////////////////////////////////////////////
void A_BROADCAST_Bds::computePosition(const rtrover_time *tt, A_VECTOR3D& xc, double& tc, A_VECTOR3D& vv, const int correction)
{
  static const double omegaEarth  = 7292115.0e-11;
  static const double gmWGS       = 398.6004418e12;
  static const double validityBrdc = 14400;
  static const double SIN_5 = -0.0871557427476582; /* sin(-5.0 deg) */
  static const double COS_5 = 0.9961946980917456; /* cos(-5.0 deg) */
  static const double NB_LEAP_SEC_BDS = 14.0; // delta time between BDS / UTC in 2006
  
  double OM, sinom, cosom, yi, dotom;
  double xg,yg,zg;  /* Beidou Geo satellite Positions */

  xc=A_VECTOR3D();
  tc=0.0;
  vv=A_VECTOR3D();
  /*std::cout << "computePosition Bds " << _eph._satellite._number << " " << tt->_mjd << " " << tt->_sec << std::endl;
  std::cout << _eph._TOC._mjd << " " << _eph._TOC._sec << std::endl;
  std::cout << _eph._TOE._mjd << " " << _eph._TOE._sec << std::endl;
  std::cout << _eph._IODE << std::endl;
  std::cout << _eph._Crs << " " << _eph._Delta_n << " " << _eph._M0 << std::endl;
  std::cout << _eph._Cuc << " " << _eph._e << " " << _eph._Cus << std::endl;
  std::cout << _eph._sqrt_A << " " << _eph._Cic << " " << _eph._OMEGA0 << std::endl;
  std::cout << _eph._Cis << " " << _eph._i0 << " " << _eph._Crc << std::endl;
  std::cout << _eph._omega << " " << _eph._OMEGADOT << " " << _eph._IDOT << std::endl;*/

  double a0 = _eph._sqrt_A * _eph._sqrt_A;
  if (a0 == 0) {
    return;
  }

  double n0 = sqrt(gmWGS/(a0*a0*a0));
  double tk = diffDate(tt, &_eph._TOE);
  if ( correction || (!correction && (tk<validityBrdc) && (tk>0)))
  {
    double n  = n0 + _eph._Delta_n;
    double M  = _eph._M0 + n*tk;
    double E  = M;
    double E_last;
    do {
      E_last = E;
      E = M + _eph._e*sin(E);
    } while ( fabs(E-E_last)*a0 > 0.001 );
    double v      = 2.0*atan( sqrt( (1.0 + _eph._e)/(1.0 - _eph._e) )*tan( E/2 ) );
    double u0     = v + _eph._omega;
    double sin2u0 = sin(2*u0);
    double cos2u0 = cos(2*u0);
    double r      = a0*(1 - _eph._e*cos(E)) + _eph._Crc*cos2u0 + _eph._Crs*sin2u0;
    double i      = _eph._i0 + _eph._IDOT*tk + _eph._Cic*cos2u0 + _eph._Cis*sin2u0;
    double u      = u0 + _eph._Cuc*cos2u0 + _eph._Cus*sin2u0;
    double xp     = r*cos(u);
    double yp     = r*sin(u);

    double deltat = (double)_eph._TOE._mjd-53736.0+_eph._TOE._sec/86400.0;
    int nweek = (int) floor(deltat/7.0);
    double TOEsec = (deltat - (nweek)*7.0)*86400.0 - NB_LEAP_SEC_BDS;
    
    double sini  = sin(i);
    double cosi  = cos(i);
    
    double omdt = omegaEarth * tk;
    double coso = cos(omdt);
    double sino = sin(omdt);

    if (_eph._satellite._number <= 5)
    {
      /* coordonnees terrestres */
      OM = _eph._OMEGA0 + _eph._OMEGADOT * tk - omegaEarth * TOEsec; /* longitude noeud ascendant */

      /* coordonnees dans un repere proche de l'ITRF */
      sinom = sin(OM);
      cosom = cos(OM);

      xg = xp * cosom - yp * sinom * cosi;
      yg = xp * sinom + yp * cosom * cosi;
      zg = yp * sini;

      yi = yg * COS_5 + zg * SIN_5;
      xc.setX( xg * coso + sino * yi);
      xc.setY( -xg * sino + coso * yi);
      xc.setZ( -yg * SIN_5 + zg * COS_5);

    }
    else
    {    
      OM     = _eph._OMEGA0 + (_eph._OMEGADOT - omegaEarth)*tk - omegaEarth*TOEsec;

      sinom = sin(OM);
      cosom = cos(OM);

      xc.setX( xp * cosom - yp * cosi * sinom);
      xc.setY( xp * sinom + yp * cosi * cosom);
      xc.setZ( yp * sini);
    }           

    double diffT = diffDate(tt, &_eph._TOC);
    tc = _eph._clock_bias + _eph._clock_drift*diffT + _eph._clock_driftrate*diffT*diffT;

    // Velocity
    // --------
    double tanv2 = tan(v/2);
    double dEdM  = 1 / (1 - _eph._e*cos(E));
    double dotv  = sqrt((1.0 + _eph._e)/(1.0 - _eph._e)) / cos(E/2)/cos(E/2) / (1 + tanv2*tanv2) 
        	 * dEdM * n;
    double dotu  = dotv + (-_eph._Cuc*sin2u0 + _eph._Cus*cos2u0)*2*dotv;

    if (_eph._satellite._number <= 5)
      dotom = _eph._OMEGADOT;
    else
      dotom = _eph._OMEGADOT - omegaEarth;
    
    double doti  = _eph._IDOT + (-_eph._Cic*sin2u0 + _eph._Cis*cos2u0)*2*dotv;
    double dotr  = a0 * _eph._e*sin(E) * dEdM * n 
                  + (-_eph._Crc*sin2u0 + _eph._Crs*cos2u0)*2*dotv;
    double dotx  = dotr*cos(u) - r*sin(u)*dotu;
    double doty  = dotr*sin(u) + r*cos(u)*dotu;
    
    double xpg = cosom   *dotx  - cosi*sinom   *doty      // dX / dr
               - xp*sinom*dotom - yp*cosi*cosom*dotom   // dX / dOMEGA
               + yp*sini*sinom*doti ;                  // dX / di

    double ypg = sinom   *dotx  + cosi*cosom   *doty
               + xp*cosom*dotom - yp*cosi*sinom*dotom
               - yp*sini*cosom*doti;

    double zpg = sini    *doty  + yp*cosi      *doti;  

    if (_eph._satellite._number <= 5)
    {
      double ypi = ypg * COS_5 + zpg * SIN_5;
      vv.setX( xpg * coso - xg * omegaEarth * sino + omegaEarth * coso * yi + sino * ypi );

      vv.setY( -xpg * sino - xg * omegaEarth * coso - omegaEarth * sino * yi + coso * ypi);

      vv.setZ( -ypg * SIN_5 + zpg * COS_5);
    }
    else
    {
      vv.setX( xpg );

      vv.setY( ypg);

      vv.setZ( zpg);
    }

    // Relativistic Correction
    // -----------------------
    //tc -= 2.0 * ( xc * vv ) / CLIGHT / CLIGHT;
    tc -= 4.442807633e-10 * _eph._e * sqrt(a0) *sin(E);
  }
  else
  {
    tc = 0.0;
    fprintf(stderr,"Eph Bds too old ! : tk = %f\n", tk);
  }
}
