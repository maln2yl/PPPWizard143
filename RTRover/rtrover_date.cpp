/************************************************************
Nom ......... : rtrover_date.cpp
Role ........ : date definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.0 9/30/2014
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#include "rtrover_date.h"

//////////////////////////////////////////////////////////////////////
A_DATE A_DATE::operator+ (const double d)
{
  int day;
  double sec; 

  day = (int) (d/86400.0);
  sec = d - (double) (86400.0*day);
  A_DATE date(_day + day,_second + sec);

  if ((date._second - 86400.0) >= 0.0)
  {
    date._second -= 86400.0;
    date._day++;
  }
  else if (date._second < 0.0)
  {
    date._second += 86400.0;
    date._day--;
  }

  return date;
}

//////////////////////////////////////////////////////////////////////
double A_DATE::operator- (const A_DATE &date) const
{
  int day;
  double sec;

  day = _day - date._day;
  sec = _second - date._second;

  return sec + (double)(day)*86400.0;
}

//////////////////////////////////////////////////////////////////////
void A_DATE::calendar(int *day, int *month, int *year, int *hour, int *minute, double *second)
{
  long i, j, k, y, m, d;
  double min;

  i = (long)_day + 712164;
  j = ((4*i-1) % 146097)/4;
  k = (((4*j+3) % 1461)+4)/4;
  y = 100 * ((4*i-1)/146097) + (4*j+3)/1461;
  m = (5*k-3)/153;
  d = (((5*k-3) % 153)+5)/5;

  if (m < 10)
  {
    *day = (int)d;
    *month = (int)(m+3);
    *year = (int)(y);
  }
  else
  {
    *day = (int)d;
    *month = (int)(m-9);
    *year = (int)(y+1);
  }

  *hour = (int)((long)_second / 3600);
  min = _second - 3600.0*(double)(*hour);
  *minute = (int)(min / 60);
  *second = min - 60.0*(double)(*minute);
}

//////////////////////////////////////////////////////////////////////
void A_DATE::fromCalendar(int day, int month, int year, int hour, int minute, double second)
{
  _day = (int)(367*(year-1950) - (7*(year+(month+9)/12))/4 + (275*month)/9 + day + 3381);
  _second = 3600.0*(double)hour + 60.0*(double)minute + second;
}
