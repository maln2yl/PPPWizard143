/************************************************************
Nom ......... : rtrover_givens_filter.cpp
Role ........ : givens filter definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.0 9/30/2014
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#include <cmath>
#include <stdio.h>

#include "rtrover_givens_filter.h"

#define TABR(I,J) ((I)*_n+(J))

//////////////////////////////////////////////////////////////////////
A_GIVENS_FILTER::A_GIVENS_FILTER(const int size)
{
  int i,j;
  
  _n = size;
  _d.clear();
  _d.resize(_n+1,0);
  _r.clear();
  _r.resize((_n+1)*(_n+1));
  for (j=0;j<=_n;j++)
  {
    for (i=0;i<=j-1;i++)
      _r[TABR(i,j)] = 0.0;
    _r[TABR(j,j)] = 1.0;
  }
  _sumRes = 0.0;
  _nb = 0;
}

//////////////////////////////////////////////////////////////////////
void A_GIVENS_FILTER::newMeasurement(const A_FILTER_MODELED_MEASUREMENT &y)
{
  int i,k;
  double delta0, deltap, dp, xi, xkp, cbar, sbar;
  double x[_n+1];

  for (i=0;i<_n;i++)
    x[i] = y._der[i];
  x[_n] = y._res;
  delta0 = y._w;
  _sumRes += y._w * y._res*y._res;
  _nb++;
  for (i=0;i<=_n;i++)
  {
    if (delta0 == 0.0)
      break;
    if (x[i] != 0.0)
    {
      dp = _d[i] + delta0 * x[i]*x[i];
      deltap = delta0 * (_d[i]/dp);
      cbar = _d[i]/dp;
      sbar = delta0 * (x[i]/dp);
      xi = x[i];
      x[i] = 0.0;
      _d[i] = dp;
      delta0 = deltap;
      for (k=i+1;k<=_n;k++)
      {
        xkp = x[k] - xi * _r[TABR(i,k)];
        _r[TABR(i,k)] = cbar * _r[TABR(i,k)] + sbar * x[k];
        x[k] = xkp;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
void A_GIVENS_FILTER::solution(double *sol, double& var)
{
  int i,j;
  double sol_t[_n];

  for (i=0;i<=_n-1;i++)
    sol_t[i] = _r[TABR(i,_n)];
  for (j=_n-1;j>=1;j--)
    for (i=j-1;i>=0;i--)
      sol_t[i] -= sol_t[j]*_r[TABR(i,j)];
  for (i=0;i<_n;i++)
    sol[i] = sol_t[i];
  var = _sumRes/(double)_nb;
}

//////////////////////////////////////////////////////////////////////
void A_GIVENS_FILTER::variance(double *sigma)
{
  int i,j,k;
  std::vector<double> _ri=_r;
  double s;

  for (j=1;j<=_n-1;j++)
    for (i=0;i<=j-1;i++)
      {
      s = _ri[TABR(i,j)];
      for (k=i+1;k<=j-1;k++)
        s -= _r[TABR(i,k)]*_r[TABR(k,j)];
      _ri[TABR(i,j)] = s;
      }
  for (i=0;i<=_n-2;i++)
    {
    s = 1.0/_d[i];
    for (j=i+1;j<=_n-1;j++)
      s += _ri[TABR(i,j)]*_ri[TABR(i,j)]/_d[j];
    sigma[i] = s;
    }
  sigma[_n-1] = 1.0/_d[_n-1];
}
