/************************************************************
Nom ......... : rtrover_kalman_filter.cpp
Role ........ : kalman filter definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.0 9/30/2014
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/

/************************************************************

OBJECTS FOR TYPE DEFINITION

************************************************************/
#include <cmath>

#include "rtrover_kalman_filter.h"

#define OK 0
#define ERREUR 1

/* Change index of upper triangular matrix to position in one dimensionnal array*/
/* start with (0,1) -> 0                                                        */
#define M_I(i,j) ( (i) + (( (j) * ((j) - 1)) / 2))

/* start with (1,2) -> 0                                                        */
#define M_IU(i,j) ( (i) - 1 + (( ((j) - 1) * ((j) - 2)) / 2))

/* start with (1,1) -> 0                                                        */
#define M_IM(i,j) ( (i) - 1 + (( (j) * ((j) - 1)) / 2))

/* Espilon to compare diagonal covariance with 0 in Kalman filter               */
static const double rg_Kalman_Epsilon =  1e-20;

/************************************************************

DECLARED PRIVATE SERVICES

************************************************************/

/************************************************************

PUBLIC SERVICES TO BE DEFINED

************************************************************/

#define MATHEMATICAL_EPSILON_ZERO   1.0E-15      
#define MATHEMATICAL_EPSILON_SQRT   1.0E-20

//////////////////////////////////////////////////////////////////////
/*static int equal_zero(
  const double re_ValueToTest,
  const double re_EpsilonZero
)
{
  return ( fabs ( re_ValueToTest ) < re_EpsilonZero ) ; 
}*/


//////////////////////////////////////////////////////////////////////
int A_KALMAN_FILTER::initialize(
  const int        ie_Size
)
{
  int il_I = 0; /* Loop index */
  int il_J = 0; /* Loop index */

  /*FDC End of Declarative Bloc */

  /* Function algorithm ... */
  int pis_Error = OK;

  /* Set the real size of Kalman filter */
  _iN = ie_Size;
  _iT = ((_iN * (_iN - 1)) >> 1);

  _trU.resize(_iT);
  _trD.resize(_iN);
  _trS.resize(_iN);

  /* Set U and D to 0 */
  for (il_I = _iN-1 ; il_I >= 0; il_I--)
  {
    /* Set D to 0 */
    _trD[il_I] = 0.;
  } /* End set D to 0 */

  /* Set the U matrix to 0 */
  il_J = M_IU(_iN - 1, _iN);
  for (il_I = il_J; il_I >= 0; il_I--)
  {
    _trU[il_I] = 0.;
  } /* End set U to 0; */

  return pis_Error;
}

//////////////////////////////////////////////////////////////////////
int A_KALMAN_FILTER::updateDiagonalCovariance(
  const int        ie_Index,
  const double     re_Value
)
{

  int il_I = 0;                              /* Index                  */
  int il_J = 0;                              /* Index                  */
  const int il_K = ie_Index - 1;             /* Index of covariance    */

  double rl_Pond = 0.0;                      /* Old covariance         */
  double rl_Temp = 0.0;                      /* Storage variable       */
  double rl_C   = 0.0;                       /* Cos term of gauss rotation*/
  double rl_S   = 0.0;                       /* Sin term of gauss rotation*/

  /* Allocation of objects */
  double trl_Mes[_iN];

  /* Function algorithm ... */
  int pis_Error = OK;

  /*FDC End of Declarative Bloc */

  /* Controle that index is in range of kalman size */
  if ( ie_Index < 1 || ie_Index > _iN)
  {
   return ERREUR;
  }

  /* Controle that covariance is positive */
  if (re_Value < 0.0)
  {
   return ERREUR;
  }

  /* Computation code is a translation from Mathlab code given in specification */
  if (ie_Index != 1)
  {
    /* Store the old diagonal value */
    rl_Pond = _trD[il_K];

    /* Store the old value of the K column of U */
    for (il_I = il_K - 1; il_I >= 0; il_I--)
    {
      /* Copy of ie_Index column of U to Mes */
      trl_Mes[il_I] = _trU[M_I(il_I, il_K)];

      /* New value of U is 0 */
      _trU[M_I(il_I, il_K)] = 0.0;
    }

    /* Scan the considered column of U */
    for (il_I = il_K - 1; il_I >= 0; il_I--)
    {

      /* Test if I component is not 0 */
      if (trl_Mes[il_I] != 0.0)
      {
	rl_Temp = _trD[il_I];

	_trD[il_I] += rl_Pond * trl_Mes[il_I] * trl_Mes[il_I];

	/* Compute rotation */
	rl_C = rl_Temp / _trD[il_I];
	rl_S = rl_Pond * trl_Mes[il_I] / _trD[il_I];
	rl_Pond *= rl_C;

	/* Apply rotation to U */
	for (il_J = il_I - 1; il_J >= 0; il_J--)
	{
	  rl_Temp = trl_Mes[il_J];
	  trl_Mes[il_J] -= trl_Mes[il_I] * _trU[M_I(il_J, il_I)];
	  _trU[M_I(il_J, il_I)] = rl_C * _trU[M_I(il_J, il_I)] + rl_S * rl_Temp;
	}

	trl_Mes[il_I] = 0.0;

      } /* end if trl_Mes[il_I] != 0.0 */

    } /* End for each element in column of U */
    
    /* Line K of U is set to 0 */
    for (il_I = il_K+1; il_I < _iN; il_I++)
    {
      _trU[M_I(il_K, il_I)] = 0.0;
    }

    /* Store the diagonal value */
    _trD[il_K] = re_Value;

  } /* end if not ie_Index == 1 */
  else
  {
    /* The first U line is set to 0 */
    for (il_I = _iN - 1; il_I > 0; il_I--)
    {
      _trU[M_I(0, il_I)] = 0.0;
    }
    _trD[0] = re_Value;
  } /* end if ie_Index == 1 */
  
  return pis_Error;
}

//////////////////////////////////////////////////////////////////////
int A_KALMAN_FILTER::recoverDiagonalCovariance(
  const int        ie_Index,
        double &   prs_Value
)
{

  int il_I = 0;                  /* Iterative index */
  const int il_K = ie_Index - 1; /* Term to recover */
  /*FDC End of Declarative Bloc */

  /* Function algorithm ... */
  int pis_Error = OK;

  /* Controle that index is in range of kalman size */
  if ( ie_Index < 1 || ie_Index > _iN)
  {
   return ERREUR;
  }

  /* Compute the result */
  prs_Value = _trD[il_K];

  for (il_I = _iN - 1; il_I > il_K; il_I--)
  {
    prs_Value += _trD[il_I] * 
		  _trU[M_I(il_K, il_I)] * _trU[M_I(il_K, il_I)];
  }
  
  return pis_Error;
}

//////////////////////////////////////////////////////////////////////
int A_KALMAN_FILTER::initializeSolution()
{

  int il_I = 0;             /* Iterative index wich scan solution item */
  /*FDC End of Declarative Bloc */

  /* Function algorithm ... */
  int pis_Error = OK;

  /* For each element in solution */
  for (il_I = _iN - 1; il_I >= 0; il_I--)
  {
    _trS[il_I] = 0.0;
  } /* End for each element in solution */

  return pis_Error;
}

//////////////////////////////////////////////////////////////////////
int A_KALMAN_FILTER::solution(
        std::vector<double> &pos_S
)
{

  int il_I = 0;             /* Iterative index wich scan solution item */
  /*FDC End of Declarative Bloc */

  /* Function algorithm ... */
  int pis_Error = OK;

  pos_S.resize(_iN);
          
  /* Return the internal solution S */
  /* For each element in solution */
  for (il_I = _iN - 1; il_I >= 0; il_I--)
  {
    pos_S[il_I] = _trS[il_I];
  }

  return pis_Error;
}

//////////////////////////////////////////////////////////////////////
int A_KALMAN_FILTER::newMeasurement(
  const A_FILTER_MODELED_MEASUREMENT &poe_Measurement,
  const double re_Sign
)
{
  int il_I = 0;                 /* Index             */
  int il_J = 0;                 /* Index             */
  int il_L = 0;                 /* Index             */

  double    rl_Trav5   = 0.0  ; /* Temporary storage */
  double    rl_UP      = 0.0  ; /* Temporary storage */
  double    rl_Wi      = 0.0  ; /* Variance of innovation */

  double trl_Trav1[_iN];
  double trl_Trav2[_iN];
  double trl_Trav3[_iN];
  double trl_Trav4[1+_iN];

  /*FDC End of Declarative Bloc */

  /* Function algorithm ... */
  int pis_Error = OK;

  if ((int)poe_Measurement._der.size() != _iN)
  {
   return ERREUR;
  }

  /* Initialisation of intermediate variable */
  il_L = 0;
  for (il_I = 0; il_I < _iN; il_I++)
  {
    trl_Trav1[il_I] = 0.0;
    trl_Trav3[il_I] = 0.0;
    trl_Trav4[il_I] = 0.0;

    for (il_J = 0; il_J < il_I; il_J++)
    {
      trl_Trav1[il_I] += _trU[il_L] * poe_Measurement._der[il_J];
      il_L++;
    }

    trl_Trav1[il_I] += poe_Measurement._der[il_I];
    trl_Trav2[il_I] = _trD[il_I] * trl_Trav1[il_I];
  } 

  trl_Trav4[_iN] = 0.0;
  trl_Trav4[0] = re_Sign*poe_Measurement._w;

  /* compute covariance matrix */
  il_L = 0;

  /* Stop if measurement is detected as invalid */
  for (il_J = 0; il_J < _iN && (pis_Error == OK); il_J++)
  {
    trl_Trav4[il_J + 1] = trl_Trav4[il_J] + _trD[il_J] * trl_Trav1[il_J] *  trl_Trav1[il_J];
    if (trl_Trav4[il_J + 1] == 0.0  || trl_Trav4[il_J] == 0.0)
    {
     return ERREUR;
    }
    else
    {
      _trD[il_J] *=  trl_Trav4[il_J] / trl_Trav4[il_J + 1] ;
  
      rl_Trav5 = trl_Trav1[il_J] / trl_Trav4[il_J];
  
      for (il_I = 0; il_I < il_J; il_I++)
      {
        rl_UP = _trU[il_L] - rl_Trav5 * trl_Trav3[il_I];
        trl_Trav3[il_I] += trl_Trav2[il_J] * _trU[il_L];
        _trU[il_L] = rl_UP;
        il_L++;
      }
	  
      trl_Trav3[il_J] += trl_Trav2[il_J];
    }
  }

  /* Compute the next part if no warning is thrown */
  if (pis_Error == OK)
  {
    /* Compute Kalman correction */
    rl_Trav5 = trl_Trav4[_iN];
    rl_Wi = rl_Trav5;
  
    /* Compute H . S */
    rl_Trav5 = 0.;
    for (il_I = _iN - 1; il_I >= 0; il_I--)
    {
      rl_Trav5 += poe_Measurement._der[il_I] * _trS[il_I];
    }

     rl_Trav5 -= poe_Measurement._res;
    /* Compute new solution */
    for (il_I = _iN - 1; il_I >= 0; il_I--)
    {
      _trS[il_I] -= rl_Trav5 * trl_Trav3[il_I] / rl_Wi;
    }
  
  }
  return pis_Error;
}

//////////////////////////////////////////////////////////////////////
int A_KALMAN_FILTER::newMeasurementsElim(
  std::vector<A_FILTER_MODELED_MEASUREMENT> &poe_Measurements,
  const double re_ResMax,
  std::vector<double> &pos_S
)

{
  int iErr = 0;
  int iMes = 0;
  int iInc = 0;
  double dResMax = 0.0;

  int pis_Error = OK;

  iErr = initializeSolution ();
  for (iMes = 0; iMes < (int)poe_Measurements.size(); iMes++)
    if (poe_Measurements[iMes]._w)
	  iErr = newMeasurement(poe_Measurements[iMes], 1.0);
    
  while (1)
  {
    iErr = solution(pos_S);
    if (!re_ResMax)
      break;
    dResMax = re_ResMax;
    int iMes_ResMax = -1;
    for (iMes=0; iMes<(int)poe_Measurements.size(); iMes++)
      {
      if (poe_Measurements[iMes]._w)
        {
	double r = poe_Measurements[iMes]._res;
        for (iInc=0; iInc<_iN; iInc++)
	   r -= pos_S[iInc]*poe_Measurements[iMes]._der[iInc];
        if (fabs(r) > dResMax)
          {
          iMes_ResMax = iMes;
          dResMax = fabs(r);
          }
        }
      }
    if (iMes_ResMax != -1)
      {
      iErr = newMeasurement(poe_Measurements[iMes_ResMax], -1.0);
      poe_Measurements[iMes_ResMax]._w = 0.0;
      poe_Measurements[iMes_ResMax]._res = dResMax;
      continue;
      }
    break;
  }

  for (iMes=0; iMes<(int)poe_Measurements.size(); iMes++)
    {
    if (poe_Measurements[iMes]._w)
      {
      for (iInc=0; iInc<_iN; iInc++)
        poe_Measurements[iMes]._res -= pos_S[iInc]*poe_Measurements[iMes]._der[iInc];
      }
    }

  return pis_Error;
}

//////////////////////////////////////////////////////////////////////
static void CombInit(int *a, int r)
{
  int i;
  for (i=0;i<r;i++)
    a[i] = i;
}

//////////////////////////////////////////////////////////////////////
static int CombGetNext (int *a, int r, int n)
{
  int i, j;
  if (!r)
    return 1;
  i = r - 1;
  while (1)
  {
    if (i < 0)
      return 1;
    if (a[i] == n - r +i)
      i--;
    else
      break;
  }
  a[i] = a[i] + 1;
  for (j = i + 1; j < r; j++)
    a[j] = a[i] + j - i;
  return 0;
}

//////////////////////////////////////////////////////////////////////
int A_KALMAN_FILTER::newMeasurementsElimNumber(
  std::vector<A_FILTER_MODELED_MEASUREMENT> &poe_Measurements,
  const double re_ResMax, const int ie_Number, std::vector<double> &pos_S)
{
  int i;
  int iMes = 0;
  int iInc = 0;
  double dResMin = 0.0;
  int iCombMes[ie_Number];
  int iCombMesMin[ie_Number];
  int pis_Error;
  A_KALMAN_FILTER backup;

  pis_Error = OK;
  backup = *this;
  dResMin = 1.0e10;
  CombInit(iCombMes, ie_Number);
  while (1)
  {
    *this = backup;
    initializeSolution ();
    for (iMes = 0; iMes < (int)poe_Measurements.size(); iMes++)
    {
      int next = 0;
      for (i=0;i<ie_Number;i++)
	if (iMes == iCombMes[i])
	  next = 1;
      if (!next)
	newMeasurement(poe_Measurements[iMes], 1.0);
    }
    solution(pos_S);
    double dRes = 0.0;
    int dep = 0;
    for (iMes=0; iMes<(int)poe_Measurements.size(); iMes++)
    {
      int next = 0;
      for (i=0;i<ie_Number;i++)
	if (iMes == iCombMes[i])
	  next = 1;
      if (!next)
      {
	double r = poe_Measurements[iMes]._res;
        for (iInc=0; iInc<_iN; iInc++)
	   r -= pos_S[iInc]*poe_Measurements[iMes]._der[iInc];
	dRes += r*r;
	if (re_ResMax && (fabs(r) > re_ResMax))
	  dep = 1;
      }
    }
    if ((dRes < dResMin) && !dep)
    {
      for (i=0;i<ie_Number;i++)
	iCombMesMin[i] = iCombMes[i];
      dResMin = dRes;
    }
    if (CombGetNext(iCombMes, ie_Number, poe_Measurements.size()))
      break;
  }
  if (dResMin > 1.0e9)
  {
    pis_Error = ERREUR;
    *this = backup;
    return pis_Error;
  }
  for (i=0;i<ie_Number;i++)
  {
    poe_Measurements[iCombMesMin[i]]._w = 0.0;
    poe_Measurements[iCombMesMin[i]]._res = dResMin;
  }
  *this = backup;
  initializeSolution ();
  for (iMes = 0; iMes < (int)poe_Measurements.size(); iMes++)
    if (poe_Measurements[iMes]._w)
      newMeasurement(poe_Measurements[iMes], 1.0);
  solution(pos_S);
  for (iMes=0; iMes<(int)poe_Measurements.size(); iMes++)
  {
    if (poe_Measurements[iMes]._w)
    {
      for (iInc=0; iInc<_iN; iInc++)
        poe_Measurements[iMes]._res -= pos_S[iInc]*poe_Measurements[iMes]._der[iInc];
    }
  }

  return pis_Error;
}

//////////////////////////////////////////////////////////////////////
int A_KALMAN_FILTER::newMeasurementsElimComb(const int param_RejetMax,
  std::vector<A_FILTER_MODELED_MEASUREMENT> &poe_Measurements, const double re_ResMax, std::vector<double> &pos_S)
{
  int iNumber = 0;
  int iNbMax = 0;
  int iMes = 0;

  if ((int)poe_Measurements.size() < param_RejetMax)
        iNbMax = poe_Measurements.size();
  else
	iNbMax = param_RejetMax; 
  for (iNumber=0;iNumber<iNbMax+1;iNumber++)
  {
    int iErr = newMeasurementsElimNumber(poe_Measurements, re_ResMax, iNumber, pos_S);
    if (iErr == OK) return OK;
  }
  for (iMes=0; iMes<(int)poe_Measurements.size(); iMes++)
    poe_Measurements[iMes]._w = 0;  
  initializeSolution ();
  solution(pos_S);
  return ERREUR;
}

//////////////////////////////////////////////////////////////////////
int A_KALMAN_FILTER::addDiagonalCovariance(
  const int        ie_Index,
  const double     re_Value
)
{
  int il_I = 0;                              /* Index                  */
  int il_J = 0;                              /* Index                  */
  const int il_K = ie_Index - 1;             /* Index of covariance    */
  const int il_N = _iN;                      /* Size of  covariance    */

  double rl_Pond = 0.0;                      /* Old covariance         */
  double rl_Temp = 0.0;                      /* Storage variable       */
  double rl_C   = 0.0;                       /* Cos term of gauss rotation*/
  double rl_S   = 0.0;                       /* Sin term of gauss rotation*/

  /* Allocation of objects */
  double trl_Mes[_iN];

  /* Function algorithm ... */
  int pis_Error = OK;

  /*FDC End of Declarative Bloc */

  /* Controle that index is in range of kalman size */
  if ( ie_Index < 1 || ie_Index > il_N)
  {
   return ERREUR;
  }

  /* Store the old diagonal value */
  rl_Pond = re_Value;

  /* Set mes to 0 */
  for (il_I = il_N - 1; il_I >= 0; il_I--)
  {
    trl_Mes[il_I] = 0;
  }
  trl_Mes[il_K] = 1;

  /* Integrate new variance in kalman covariance matrix */
  for (il_I = il_N - 1; il_I >= 0; il_I--)
  {

    /* Test if I component is not 0 */
    if (trl_Mes[il_I] != 0.0)
    {
      rl_Temp = _trD[il_I];

      _trD[il_I] += rl_Pond * trl_Mes[il_I] * trl_Mes[il_I];

      /* Test that the variance is positive */
      if (_trD[il_I] < 0)
      {
        return ERREUR;
      }
      /* Test the assimilation of variance */
      if (_trD[il_I] != 0.)
      {
        /* Compute rotation */
        rl_C = rl_Temp / _trD[il_I];
        rl_S = rl_Pond * trl_Mes[il_I] / _trD[il_I];
        rl_Pond *= rl_C;

        /* Apply rotation to U */
        for (il_J = il_I - 1; il_J >= 0; il_J--)
        {
          rl_Temp = trl_Mes[il_J];
          trl_Mes[il_J] -= trl_Mes[il_I] * _trU[M_I(il_J, il_I)];
          _trU[M_I(il_J, il_I)] = rl_C * _trU[M_I(il_J, il_I)] + rl_S * rl_Temp;
        }

        trl_Mes[il_I] = 0.0;

      } /* end if _trD[il_I] != 0.0 */
    } /* end if trl_Mes[il_I] != 0.0 */

  } /* End for each element in column of U */

  return pis_Error;
}
