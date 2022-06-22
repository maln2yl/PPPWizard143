/************************************************************
Nom ......... : rtrover_model.cpp
Role ........ : model definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#include <cmath>
#include <iostream>

#include "rtrover_vector3D.h"
#include "rtrover_utility.h"
#include "rtrover_frequency.h"
#include "rtrover_model.h"
#include "rtrover_date.h"
#include "rtrover_map.h"

static const double TWO_PI = 2. * M_PI;
static const double RAD2DEG = 180. / M_PI;

static const double omegaT = 7.2921151467E-5; /* (rd/s) */
static const double muGps = 3.986005e+14; /* (m**3/s**2) */

//////////////////////////////////////////////////////////////////////
static double modulo(const double x, const double a, const double   b)
{
  double rl_Modulo = 0.0 ;
  double rl_Base = 0.0 ;
  double rl_Offset = x - a ;
  rl_Base = b - a ;
  if ( ( rl_Offset > 0.0 ) && ( rl_Offset < rl_Base ) )
    rl_Modulo = rl_Offset ;
  else
    rl_Modulo =  rl_Offset - (rl_Base * floor ( rl_Offset / rl_Base ) ) ;

  return ( a + rl_Modulo ) ;
}

//////////////////////////////////////////////////////////////////////
/*------------------------------------------------------------------*/
/*
Sidereal time VEIS, or hour angle, is defined using the following formula:
Tsid = a1 + WT(t - t0)
t = tsystem - (SYSTEM - UTC)
t0 corresponds to a reference date (1/January/50).
There is a1=1.746647708617871rad  and WT = Earth_rotation_rate.
*/
/*------------------------------------------------------------------*/
static double siderealTime(const int Day, const double Sec, const double SystemUtc)

{
  double res;

  /* a1 offset of sidereal time at t0        */
  const double rl_A1 = 1.743078993205;
  /* 1,743078993205 : 14610 */

  /* t0 1/1/90 Julian Days */
  const long il_T0 = 14610;

  /* Offset in radian between a complete rotation and the Earth rotation in one days */
  const double rl_Eps = 1.72021795692135e-2;

  /* Days since t0                          */
  long il_DeltaDays = 0;

  /* Adjust seconds                         */
  double rl_AdjustSeconds = 0.0;

  /* Compute t - t0 in days */
  il_DeltaDays = Day - il_T0;

  /* Adjust seconds         */
  rl_AdjustSeconds = Sec - SystemUtc;

  /* Compute sideral time */
  res = rl_A1 + (rl_Eps * ((double) il_DeltaDays)) +
		     (omegaT * rl_AdjustSeconds);

  /* Position sideral time in the right quadrant */
  res = modulo(res, 0.0, TWO_PI);
  
  return res;

}

//////////////////////////////////////////////////////////////////////
static A_VECTOR3D moonPos(const int Day, const double Sec, const double SystemUtc)
{
  double   S, C, xv, yv, zv;
  double   rl_SiderealTime;
  double   rl_F = 0.0; 
  double   rl_D = 0.0; 
  double   rl_XLP = 0.0; 
  double   rl_CE = 0.0; 
  double   rl_SE = 0.0;  
  double   rl_CR = 0.0; 
  double   rl_SR = 0.0; 
  double   rl_Q = 0.0; 
  //double   rl_G = 0.0;
  double   rl_TE_Days = 0.0; 
  double   rl_TE_Seconds = 0.0; 
  double   rl_E = 0.0;
  double   rl_ROT = 0.0;
  double   rl_XL = 0.0; 
  double   rl_U = 0.0; 
  double   rl_DL = 0.0; 
  double   rl_B = 0.0; 
  double   rl_CU = 0.0; 
  double   rl_SU = 0.0; 
  double   rl_CB = 0.0; 
  double   rl_RX = 0.0; 
  double   rl_RY = 0.0; 
  double   rl_RZ = 0.0; 
  double   rl_DASR = 0.0; 
  double   rl_SB = 0.0; 
  double   rl_RL = 0.0;
  double   rl_S1 = 0.0;
  double   rl_S2 = 0.0;  
  double   rl_S3 = 0.0; 
  double   rl_S4 = 0.0; 
  double   rl_SXLP = 0.0;
  double   rl_C1 = 0.0; 
  double   rl_C2 = 0.0; 
  double   rl_C3 = 0.0; 
  double   rl_C4 = 0.0; 
  double   rl_CXLP = 0.0;
  double   rl_K1 = 0.0;                /* rl_XL - 2.0 * rl_D */
  double   rl_K2 = 0.0;                /* rl_XL + rl_F */
  double   rl_K3 = 0.0;                /* rl_XL - rl_F */
  double   rl_K4 = 0.0;                /* rl_F - 2.0 * rl_D */

  rl_TE_Seconds = ((Day * 86400.0) + Sec) - SystemUtc ;
  rl_TE_Days    = Day + (Sec - SystemUtc) / 86400.0 ;
  
  /* Simplified NEWCOMB theory (VEIS frame) */ 
  rl_ROT = 7.082201389e-12 * rl_TE_Seconds;
  //rl_G = 4.931445255 + (-8.203047484e-3 + 8.203047484e-7 * rl_TE_Days) ;
  rl_E = 4.091440975e-1 + (6.217909993e-5 - 6.217909993e-9  * rl_TE_Days);
  rl_F = 3.940394946 + (-2.308957241e+3 + 2.308957241e-1 * rl_TE_Days ) ;
  rl_XL = 3.236782911 + (- 2.280271271e+3 + 2.280271271e-1 * rl_TE_Days);
  rl_XLP = 2.338793558 + (- 1.720196511e+2 + 1.720196511e-2 * rl_TE_Days) ;
  rl_U = 1.192775464 + ( -2.299715112e+3 + 2.299715112e-1 * rl_TE_Days) ;
  rl_D = 2.057045056e-1 + (- 2.127687083e+3 + 2.127687083e-1 * rl_TE_Days) ;     
   
  rl_SR = sin(rl_ROT); rl_CR = cos(rl_ROT);
  rl_SE = sin(rl_E); rl_CE = cos(rl_E);

  rl_K1 = 3.236782911 
    -2.0 * 2.057045056e-1 
    + (-2.280271271e+3 + ((2.280271271e-1 - 2.0 * 2.127687083e-1 ) * rl_TE_Days) + 2.0 * 2.127687083e+3 )
    ; /* rl_XL - 2.0 * rl_D */

  rl_K2 = (3.236782911 + 3.940394946) 
    + 
    (-2.280271271e+3 - 2.308957241e+3 + (2.280271271e-1 + 2.308957241e-1) * rl_TE_Days)
    ; /* rl_XL + rl_F */

  rl_K3 = (3.236782911 - 3.940394946) 
    + 
    (-2.280271271e+3 + 2.308957241e+3 + (2.280271271e-1 - 2.308957241e-1) * rl_TE_Days)
    ; /* rl_XL - rl_F */

  rl_K4 = -2.0 * 2.057045056e-1 
    + 3.940394946 
    + ( 2.0 * 2.127687083e+3 - 2.308957241e+3 + (2.308957241e-1 - 2.0 * 2.127687083e-1) * rl_TE_Days )
    ; /* rl_F - 2.0 * rl_D */

  /* MOON */
  rl_S1 = sin(rl_XL); rl_C1 = cos(rl_XL);
  rl_S2 = sin(rl_K1); rl_C2 = cos(rl_K1);
  rl_S3 = sin(2.0 * rl_D); rl_C3 = cos(2.0 * rl_D);
  rl_S4 = sin(2.0 * rl_XL); rl_C4 = cos(2.0 * rl_XL);
  rl_SXLP = sin(rl_XLP); rl_CXLP = cos(rl_XLP);

  rl_DASR = 0.0545 * rl_C1 + 0.01002 * rl_C2 + 0.00825 * rl_C3 + 0.00297 * rl_C4 - 0.00012 * rl_CXLP;
  rl_DL = 0.10976 * rl_S1 - 0.02224 * rl_S2 + 0.01149 * rl_S3 + 0.00373 * rl_S4 - 0.00324 * rl_SXLP;
  rl_U = rl_U + rl_DL;
  rl_SU = sin(rl_U); rl_CU = cos(rl_U);
  
  rl_S1 = sin(rl_F); rl_C1 = cos(rl_F);
  rl_S2 = sin(rl_K2); rl_C2 = cos(rl_K2);
  rl_S3 = sin(rl_K3); rl_C3 = cos(rl_K3);
  rl_S4 = sin(rl_K4); rl_C4 = cos(rl_K4);
  rl_B = 0.0895 * rl_S1 + 0.00490 * rl_S2 + 0.00485 * rl_S3 - 0.00302 * rl_S4;
  rl_SB = sin(rl_B); rl_CB = cos(rl_B);
  
  rl_RX = rl_CU * rl_CB;
  rl_Q = rl_SU * rl_CB;
  rl_RY = rl_Q * rl_CE - rl_SB * rl_SE;
  rl_RZ = rl_SB * rl_CE + rl_Q * rl_SE;
  
  rl_Q = rl_RX * rl_CR + rl_RY * rl_SR;
  rl_RY = rl_RY * rl_CR - rl_RX * rl_SR;
  rl_RX = rl_Q;
  
  rl_RL = 384389000.30 / (1.0 + rl_DASR);
         
  xv = rl_RL * rl_RX;
  yv = rl_RL * rl_RY;
  zv = rl_RL * rl_RZ;

  rl_SiderealTime = siderealTime(Day, Sec, SystemUtc);

  S = sin(rl_SiderealTime);
  C = cos(rl_SiderealTime);

  return A_VECTOR3D(C*xv + S*yv, -S*xv + C*yv, zv);
}

//////////////////////////////////////////////////////////////////////
static A_VECTOR3D sunPos(const int Day, const double Sec, const double SystemUtc)
{
  double   S, C, xv, yv, zv;
  double   rl_SiderealTime;
  //double   rl_F = 0.0; 
  //double   rl_D = 0.0;
  double   rl_XLP = 0.0; 
  double   rl_CE = 0.0; 
  double   rl_SE = 0.0; 
  double   rl_CR = 0.0; 
  double   rl_SR = 0.0; 
  double   rl_Q = 0.0; 
  double   rl_G = 0.0; 
  double   rl_CL = 0.0; 
  double   rl_SL = 0.0; 
  double   rl_SX = 0.0; 
  double   rl_SY = 0.0; 
  double   rl_SZ = 0.0; 
  double   rl_TE_Days = 0.0; 
  double   rl_TE_Seconds = 0.0;  
  double   rl_RS = 0.0; 
  double   rl_E = 0.0; 
  double   rl_ROT = 0.0;
  //double   rl_XL = 0.0;
  //double   rl_U = 0.0; 
  double   rl_S1 = 0.0; 
  double   rl_S2 = 0.0; 
  double   rl_S3 = 0.0; 
  //double   rl_SXLP = 0.0;
  double   rl_C1 = 0.0; 
  double   rl_C2 = 0.0; 
  double   rl_C3 = 0.0; 
  double   rl_CXLP = 0.0;
  double   rl_K1 = 0.0;            /* rl_XLP + rl_G */
  double   rl_K2 = 0.0;            /* 2.0*rl_XLP + rl_G */

  rl_TE_Seconds = ((Day * 86400.0) + Sec) - SystemUtc ;
  rl_TE_Days    = Day + (Sec - SystemUtc) / 86400.0 ;
  
  /* Simplified NEWCOMB theory (VEIS frame) */ 
  rl_ROT = 7.082201389e-12 * rl_TE_Seconds;
  rl_G = 4.931445255 + ( -8.203047484e-3 + (8.203047484e-7 * rl_TE_Days) ) ;
  rl_E = 4.091440975e-1 + ( 6.217909993e-5 - (6.217909993e-9  * rl_TE_Days) ) ;
  //rl_F = 3.940394946 + (-2.308957241e+3 + (2.308957241e-1 * rl_TE_Days)) ;
  //rl_XL = 3.236782911 + (-2.280271271e+3 + (2.280271271e-1  * rl_TE_Days) ) ;
  rl_XLP = 2.338793558 + ( -1.720196511e+2 + (1.720196511e-2  * rl_TE_Days ) ) ;
  //rl_U = 1.192775464 + (-2.299715112e+3 + (2.299715112e-1  * rl_TE_Days ) ) ;
  //rl_D = 2.057045056e-1 + ( -2.127687083e+3 + (2.127687083e-1  * rl_TE_Days ) ) ;

  rl_K1 =  ((-8.203047484e-3 + (8.203047484e-7 * rl_TE_Days) ) 
	    + (2.338793558 + 4.931445255)) + 
           ( -1.720196511e+2 + (1.720196511e-2  * rl_TE_Days ) )  ;  /* rl_XLP + rl_G */ 

  rl_K2 =  ((-8.203047484e-3 + (8.203047484e-7 * rl_TE_Days) ) 
	   + ((2.0 * 2.338793558) + 4.931445255)) + 
           ( 2.0 * (-1.720196511e+2 + (1.720196511e-2  * rl_TE_Days )) )  ;  /* 2.0*rl_XLP + rl_G */ 

  rl_SR = sin(rl_ROT); rl_CR = cos(rl_ROT);
  rl_SE = sin(rl_E); rl_CE = cos(rl_E);

  /* SUN */
  rl_S1 = sin(rl_K1); rl_C1 = cos(rl_K1);
  rl_S2 = sin(rl_K2); rl_C2 = cos(rl_K2);
  rl_S3 = sin(rl_G); rl_C3 = cos(rl_G);

  rl_SL = 99972.0 * rl_S1 + (1671.0 * rl_S2 - 1678.0 * rl_S3 );
  rl_CL = 99972.0 * rl_C1 + (1671.0 * rl_C2 - 1678.0 * rl_C3 );

  rl_Q = sqrt(rl_CL * rl_CL + rl_SL * rl_SL);
  rl_CL = rl_CL / rl_Q;
  rl_SL = rl_SL / rl_Q;
  rl_SX = rl_CL;
  rl_SY = rl_SL * rl_CE;
  rl_SZ = rl_SL * rl_SE;
  rl_Q = rl_SX * rl_CR + rl_SY * rl_SR;
  rl_SY = rl_SY * rl_CR - rl_SX * rl_SR;
  rl_SX = rl_Q;

  //rl_SXLP = sin(rl_XLP);
  rl_CXLP = cos(rl_XLP);
  rl_RS = 0.149597870E12 / (1.0 + 0.016722 * rl_CXLP);

  xv = rl_RS * rl_SX;
  yv = rl_RS * rl_SY;
  zv = rl_RS * rl_SZ;

  rl_SiderealTime = siderealTime(Day, Sec, SystemUtc);

  S = sin(rl_SiderealTime);
  C = cos(rl_SiderealTime);

  return A_VECTOR3D(C*xv + S*yv, -S*xv + C*yv, zv);
}

//////////////////////////////////////////////////////////////////////
/* M_xyz->Ter = [X Y Z] */
static void attitudeMatrix_old(const int Day, const double Sec, A_VECTOR3D &X, A_VECTOR3D &Y, A_VECTOR3D &Z,
                                     const A_VECTOR3D &posSat)
{
  A_VECTOR3D posSun = sunPos(Day, Sec, 0);
  posSun.doubleEtModule(1.0);

  Z = posSat; Z.doubleEtModule(-1.0);
  Y = posSun ^ Z; Y.doubleEtModule(1.0);
  X = -Z ^ Y; X.doubleEtModule(1.0);

  if ((X * posSun) < 0.0)
  {
    X = -X;
    Y = -Y;
  }

}

//////////////////////////////////////////////////////////////////////
/* M_xyz->Ter = [X Y Z] */
static void attitudeMatrix(const int Day, const double Sec, const int slot_typeBds, A_VECTOR3D &X, A_VECTOR3D &Y, A_VECTOR3D &Z,
			   const A_VECTOR3D &posSat, const A_VECTOR3D &velSat, const double yawAngle)
{
  // Inertial velocity
  const A_VECTOR3D omegaVector(-omegaT * posSat.getY(), omegaT * posSat.getX(), 0.);
  const A_VECTOR3D inertialVel = velSat + omegaVector;

  // Calculation of attitude matrix
  Z = posSat;
  Z.doubleEtModule(-1.);

  A_VECTOR3D X0 = inertialVel;
  X0.doubleEtModule(1.);
  const A_VECTOR3D Y0 = Z ^ X0;

  // Beta angle (angle between orbital plane and direction of Sun)
  A_VECTOR3D N = posSat ^ inertialVel;
  N.doubleEtModule(1.);

  A_VECTOR3D posSun = sunPos(Day, Sec, 0);
  posSun.doubleEtModule(1.);

  const double sinbeta = N * posSun;
  const double betaDeg = asin(sinbeta) * RAD2DEG;

  // Psi angle
  double psi = yawAngle;

  // For BeiDou GEO, or BeiDou IGSO/MEO with beta angle < 4 deg
  if ((slot_typeBds == GEO) || ((fabs(betaDeg) < 4.) && ((slot_typeBds == IGSO) || (slot_typeBds == MEO)))) {
    psi = 0.;
  }

  if (fabs(psi) < EPSDBL) {

    // Calculation of psi angle
    A_VECTOR3D u = posSun ^ N;
    u.doubleEtModule(1.);

    A_VECTOR3D sat = posSat;
    sat.doubleEtModule(1.);

    const double tanbeta = sinbeta / sqrt(1. - sinbeta * sinbeta);
    const double sinmu = u * sat;
    psi = atan2(-tanbeta, sinmu);
    if (psi < 0.)
      psi += TWO_PI;
  }

  const double cospsi = cos(psi);
  const double sinpsi = sin(psi);
  X = X0 * cospsi + Y0 * sinpsi;
  Y = Y0 * cospsi - X0 * sinpsi;
}

//////////////////////////////////////////////////////////////////////
static double windup (const A_VECTOR3D &XRec, const A_VECTOR3D &YRec, const A_VECTOR3D &posRec,
                      const A_VECTOR3D &Xemi, const A_VECTOR3D &Yemi, const A_VECTOR3D &posemi)
{

  A_VECTOR3D DRecemi = posRec - posemi; DRecemi.doubleEtModule(1.0);

  /* Demi = Xemi - DRecemi*(Xemi.DRecemi) - DemiRec^Yemi */
  A_VECTOR3D Demi = Xemi - DRecemi*(Xemi * DRecemi) - (DRecemi ^ Yemi);
  Demi.doubleEtModule(1.0);

  /* DRec = XRec - DRecemi*(XRec.DRecemi) + DemiRec^YRec */
  A_VECTOR3D Drec = XRec - DRecemi*(XRec * DRecemi) + (DRecemi ^ YRec);
  Drec.doubleEtModule(1.0);

  /* sinX. */
  double dSinX = (Drec ^ Demi) * DRecemi;

  /* cosX. */
  double dCosX = Demi * Drec;

  /* windup angle */
  double dCorPhase = atan2 (dSinX, dCosX)/TWO_PI;
  return dCorPhase;
}

//////////////////////////////////////////////////////////////////////
static void geodesicFrame (A_VECTOR3D &X, A_VECTOR3D &Y, A_VECTOR3D &Z, const A_VECTOR3D &Hsta)
{
  A_VECTOR3D North(0.0, 0.0, 1.0);

  Z = Hsta;
  X = North ^ Z;
  X.doubleEtModule(1.0);
  Y = Z ^ X;
 }

//////////////////////////////////////////////////////////////////////
static double mapping(const double CosPhi)
{
  double SinPhi;
  SinPhi=sqrt(1.-CosPhi*CosPhi);
  return 1./(CosPhi+0.00143*SinPhi/(CosPhi+0.0455*SinPhi));
}

//////////////////////////////////////////////////////////////////////
/* Reference IERS 96, Borkowski, 1989                         */
void A_MODEL::geodesic(const double x, const double y, const double z, double *lambda, double *phi, double *height)
{
  double r, b, E, F, P, Q, D, st, s, v, G, t;
  static const double Re = 6378137.0;
  static const double Fl = 298.257222101;

  r = sqrt(x*x+y*y);
  b = fabs(Re - Re/Fl);
  if (z < 0.0) b = -b;
  E = ((z+b)*b/Re-Re)/r;
  F = ((z-b)*b/Re+Re)/r;
  P = (E*F + 1.0)*4.0/3.0;
  Q = (E*E - F*F)*2.0;
  D = P*P*P + Q*Q;
  if (D >= 0.0)
    {
      st = sqrt(D)+Q;
      s = fabs(exp(log(fabs(st))/3.0));
      if (st < 0.0) s = -s;
      v = P/s - s;
      v = -(Q+Q+v*v*v)/(3.0*P);
    }
  else
  v = 2.0*sqrt(-P)*cos(acos(Q/P/sqrt(-P))/3.0);
  G = 0.5*(E+sqrt(E*E+v));
  t = sqrt(G*G + (F-v*G)/(G+G-E)) - G;
  (*lambda) = atan2(y, x);
  (*phi) = atan((1.0-t*t)*Re/(2.0*b*t));
  (*height) = (r-Re*t)*cos(*phi) + (z-b)*sin(*phi);
}

//////////////////////////////////////////////////////////////////////
static A_VECTOR3D terrestrialTideObject(double GMA_GMT, const A_VECTOR3D &object, const A_VECTOR3D &sta)
{
  static const double H2 = 0.609;
  static const double L2 = 0.0852;
  static const double H3 = 0.292;
  static const double L3 = 0.015;

  double rsta, robject;
  double proj;
  double coef2, coef3, coef_object2, coef_sta2, coef_object3, coef_sta3;

  rsta = sta.getModule();
  robject = object.getModule();

  A_VECTOR3D usta = sta; usta.doubleEtModule(1.0);
  A_VECTOR3D uobject = object; uobject.doubleEtModule(1.0);

  proj = usta*uobject;

  coef2 = GMA_GMT*rsta*(rsta/robject)*(rsta/robject)*(rsta/robject);
  coef3 = coef2*(rsta/robject);
  coef_object2 = 3.0*L2*proj;
  coef_sta2 = 3.0*(H2/2.0-L2)*proj*proj - H2/2.0;
  coef_object3 = 3.0*L3/2.0*(5.0*proj*proj-1.0);
  coef_sta3 = 5.0/2.0*(H3-3.0*L3)*proj*proj*proj+3.0/2.0*(L3-H3)*proj;

// Order 3 set to 0
  A_VECTOR3D v = (uobject*coef_object2 + usta*coef_sta2)*coef2+(uobject*coef_object3 + usta*coef_sta3)*coef3*0.0;
  
  return v;
}

//////////////////////////////////////////////////////////////////////
static A_VECTOR3D terrestrialTideK1(const double ts, const A_VECTOR3D &sta)
{
  double rsta = sta.getModule();
  double cf = rsta*rsta*rsta;
  double coef_k1 = -(0.0253*sta.getZ()/cf)*(sta.getY()*cos(ts)+sta.getX()*sin(ts));
  return sta*coef_k1;
}


//////////////////////////////////////////////////////////////////////
static A_VECTOR3D terrestrialTide(const int Day, const double Sec,
                                     const double SystemUtc,
                                     const A_VECTOR3D &sta)
{
  static const double GML_GMT = 0.012300034;
  static const double GMS_GMT = 0.332946038E6;

  A_VECTOR3D moon=moonPos(Day, Sec, SystemUtc);
  A_VECTOR3D sun=sunPos(Day, Sec, SystemUtc);
  A_VECTOR3D imoon = terrestrialTideObject(GML_GMT, moon, sta);
  A_VECTOR3D isun = terrestrialTideObject(GMS_GMT, sun, sta);
  double ts = siderealTime(Day, Sec, SystemUtc);
  return imoon + isun + terrestrialTideK1(ts, sta);
}

//////////////////////////////////////////////////////////////////////
static double shapiro(const A_VECTOR3D &M, const A_VECTOR3D &Me)
{
  double r1 = M.getModule();
  double r2 = Me.getModule();
  A_VECTOR3D d = Me-M;
  double rho = d.getModule();
  return 2.0*muGps/(CLIGHT*CLIGHT)*log((r1+r2+rho)/(r1+r2-rho));
}

// Tropospheric Model (Saastamoinen)
//////////////////////////////////////////////////////////////////////
double A_MODEL::tropDelay(const double cosPhi, const double h) {

  double height = h;
  if (height < 0.0) height = 0.0;
  if (height > 5000.0) height = 5000.0;
  double pp =  1013.25 * pow(1.0 - 2.26e-5 * height, 5.225);
  double TT =  18.0 - height * 0.0065 + 273.15;
  double hh =  50.0 * exp(-6.396e-4 * height);
  double ee =  hh / 100.0 * exp(-37.2465 + 0.213166*TT - 0.000256908*TT*TT);

  double h_km = height / 1000.0;
  
  int    ii   = (int)(h_km + 1);
  double href = ii - 1;
  
  double bCor[7]; 
  bCor[0] = 1.156;
  bCor[1] = 1.006;
  bCor[2] = 0.874;
  bCor[3] = 0.757;
  bCor[4] = 0.654;
  bCor[5] = 0.563;
  bCor[6] = 0.472;
  double BB = bCor[ii-1] + (bCor[ii]-bCor[ii-1]) * (h_km - href);
  double r = (0.002277/cosPhi) * (pp + ((1255.0/TT)+0.05)*ee - BB*((1-cosPhi*cosPhi)/(cosPhi*cosPhi)));
  if (r < 0.0) r = 0.0;
  
  return r;
}

//////////////////////////////////////////////////////////////////////
A_MODELED_MEASUREMENT A_MODEL::model(const A_DATE &date,
                  const int isyst,
                  const A_VECTOR3D &posRec,
                  const A_MEASUREMENT &mes,
                  const A_SSR_PARAMETER &positions,
		  const A_SATELLITE_PHASE_MAP &satellites)
{
  double C0b;
  double cosPhi, tro;
  double sha;
  double Lon, Lat, Alt, Hx, Hy, Hz;
  A_MODELED_MEASUREMENT mod;

  if ((mes._code[F1]) && (posRec.getX() || posRec.getY() || posRec.getZ()))
  {
    A_VECTOR3D M(posRec);
    A_VECTOR3D Me(positions._Xsp3, positions._Ysp3, positions._Zsp3);

    /* Tide */
    M += terrestrialTide(date.getDay(), date.getSecond(), 0.0, M);

    /* Distance */
    A_VECTOR3D D = Me - M;
    C0b  = D.getModule();

    /* Sagnac */
    A_VECTOR3D Omega(0.0, 0.0, omegaT);
    A_VECTOR3D Vr = Omega ^ M;
    double Dist = -(Vr * D) /CLIGHT;

    /* Shapiro */
    sha = shapiro(M, Me);

    /* Elevation cosine */
    geodesic(M.getX(), M.getY(), M.getZ(),&Lon,&Lat,&Alt);
    Hx = cos(Lon)*cos(Lat);
    Hy = sin(Lon)*cos(Lat);
    Hz = sin(Lat);
    A_VECTOR3D H(Hx, Hy, Hz);
    cosPhi = (H * D) /C0b;

    /* Elevation test */
    if (cosPhi >= 0.0)
    {
      /* Troposphere */
      tro = tropDelay(cosPhi, Alt);

      /* Windup */
      A_VECTOR3D Xsat, Ysat, Zsat;
      A_VECTOR3D Xsta, Ysta, Zsta;
      A_VECTOR3D posSat(Me);
      A_VECTOR3D velSat(positions._VXsp3, positions._VYsp3, positions._VZsp3);
      A_VECTOR3D posSta(Hx, Hy, Hz);
      double corPhase;
      double yaw = positions._yaw * M_PI;
      attitudeMatrix(date.getDay(), date.getSecond(), mes._slot_typeBds, Xsat, Ysat, Zsat, posSat, velSat, yaw);
      geodesicFrame(Xsta, Ysta, Zsta, posSta);
      corPhase=windup(Xsta, Ysta, posSta, Xsat, Ysat, posSat);
      corPhase = corPhase + rint (_wu - corPhase);
      _wu = corPhase;

      /* Phase map */
      double mapPhase1, mapPhase2;
      satellites.phaseCenterVariation(M, Me, &mapPhase1, &mapPhase2);
      
      for(int f=F1; f<FMax;f++)
      {
        mod._modCode[f] = C0b + Dist + sha + tro;
	if (f == F1)
	  mod._modPhase[F1] = (C0b + Dist + sha + tro + mapPhase1)/FREQUENCY::LAMBDA(F1, isyst, mes._slot_typeBds) - corPhase;
	else
	  mod._modPhase[f] = (C0b + Dist + sha + tro + mapPhase2)/FREQUENCY::LAMBDA((Frequency)f, isyst, mes._slot_typeBds) - corPhase;
      }

      mod._mapping = mapping(cosPhi);
      mod._elevation=(asin(cosPhi)*180)/M_PI;
      mod._pos = -D/C0b;
    }
  }
  
  return mod;
}
