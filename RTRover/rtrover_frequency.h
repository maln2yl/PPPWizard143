/************************************************************
Nom ......... : rtrover_frequency.h
Role ........ : frequency plan
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef _RTROVER_FREQUENCY
#define _RTROVER_FREQUENCY

#include "rtrover_utility.h"

class FREQUENCY
{
#define GPS_F0 (10.23e+6)
#define GPS_F1 (154.0*GPS_F0)
#define GPS_F2 (120.0*GPS_F0)
#define GPS_F5 (115.0*GPS_F0)
#define GLO_F1B (1602e+6)
#define GLO_F2B (1246e+6)
#define GLO_F1O (0.5625e+6)
#define GLO_F2O (0.4375e+6)
#define GAL_F5B (118.0*GPS_F0)
#define BDS_F1 (152.6*GPS_F0)
#define BDS_F6 (124.0*GPS_F0)
  public :
    inline static double F(Frequency f, int isyst, int slot)
    {
    static double x[]={GPS_F1, GPS_F2, GPS_F5, 1.0, 1.0};
    static double a[]={GLO_F1B, GLO_F2B, 1.0, 1.0, 1.0};
    static double b[]={GLO_F1O, GLO_F2O, 1.0, 1.0};
    static double y[]={GPS_F1, 1.0, GPS_F5, 1.0, GAL_F5B};
    static double z[]={BDS_F1, 1.0, 1.0, BDS_F6, GAL_F5B};
    if (isyst == GPS) return x[f];
    if (isyst == Glo) return a[f]+b[f]*(double)slot;
    if (isyst == Gal) return y[f];
    if (isyst == Bds) return z[f];
    };
    
    inline static double LAMBDA(Frequency f, int isyst, int slot)
    {
    static double x[]={CLIGHT/GPS_F1, CLIGHT/GPS_F2, CLIGHT/GPS_F5, 1.0, 1.0};
    static double a[]={GLO_F1B, GLO_F2B, 1.0, 1.0, 1.0};
    static double b[]={GLO_F1O, GLO_F2O, 1.0, 1.0, 1.0};
    static double y[]={CLIGHT/GPS_F1, 1.0, CLIGHT/GPS_F5, 1.0, CLIGHT/GAL_F5B};
    static double z[]={CLIGHT/BDS_F1, 1.0, 1.0, CLIGHT/BDS_F6, CLIGHT/GAL_F5B};
    if (isyst == GPS) return x[f];
    if (isyst == Glo) return CLIGHT/(a[f]+b[f]*(double)slot);
    if (isyst == Gal) return y[f];
    if (isyst == Bds) return z[f];
    };
    
    inline static double GAMMA(Frequency f, int isyst, int slot)
    {
    static double x[]={1.0, GPS_F1/GPS_F2*GPS_F1/GPS_F2, GPS_F1/GPS_F5*GPS_F1/GPS_F5, 1.0, 1.0};
    static double a[]={GLO_F1B, GLO_F2B, 1.0, 1.0, 1.0};
    static double b[]={GLO_F1O, GLO_F2O, 1.0, 1.0, 1.0};
    static double y[]={1.0, 1.0, GPS_F1/GPS_F5*GPS_F1/GPS_F5, 1.0, GPS_F1/GAL_F5B*GPS_F1/GAL_F5B};
    static double z[]={1.0, 1.0, 1.0, BDS_F1/BDS_F6*BDS_F1/BDS_F6, BDS_F1/GAL_F5B*BDS_F1/GAL_F5B};
    if (isyst == GPS) return x[f];
    if (isyst == Glo)
    {
      double f1=a[0]+b[0]*(double)slot;
      double ff=a[f]+b[f]*(double)slot;
      return (f1/ff)*(f1/ff);
    }
    if (isyst == Gal) return y[f];
    if (isyst == Bds) return z[f];
    };
};

#endif
