/************************************************************
Nom ......... : rtrover_filter.cpp
Role ........ : filter definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>

#include "rtrover_utility.h"
#include "rtrover_kalman_filter.h"
#include "rtrover_model.h"
#include "rtrover_frequency.h"
#include "rtrover_filter.h"

//////////////////////////////////////////////////////////////////////
// Hard-coded settings
#define SIG_INI_CLK (1.0e8*1.0e8)
#define SIG_MOD_CLK (1.0e8*1.0e8)
#define SIG_INI_AMB (1.0e8*1.0e8)
#define EPS_FIX 1.0e-6
#define THR_AMB 0.25

A_FILTER::A_FILTERLOG A_FILTER::_filterLog;

//////////////////////////////////////////////////////////////////////
A_LISTPARAM::A_LISTPARAM()
{
  _iDebX = 0; _iNbX = 0;
  _iDebY = 0; _iNbY = 0;
  _iDebZ = 0; _iNbZ = 0;
  _iDebT = 0; _iNbT = 0;
  _iDebe = 0; _iNbe = 0;
  _iDebHP.clear();
  _iDebHl.clear();
  _iNbHP.clear();
  _iNbHl.clear();  
  _iDebHP.resize(SystMax);
  _iDebHl.resize(SystMax);
  _iNbHP.resize(SystMax);
  _iNbHl.resize(SystMax);
  for(int syst=0; syst<SystMax; syst++)
  {
   _iDebHP[syst].clear();
   _iDebHl[syst].clear();
   _iNbHP[syst].clear();
   _iNbHl[syst].clear();
   _iDebHP[syst].resize(FMax,0.0);
   _iDebHl[syst].resize(FMax,0.0);
   _iNbHP[syst].resize(FMax,0.0);
   _iNbHl[syst].resize(FMax,0.0);
  }
  _iDebP = 0; _iNbP = 0;
  _iDebA4 = 0; _iNbA4 = 0;
  _iDebEx = 0; _iNbEx = 0;
  _iDebJP = 0; _iNbJP = 0;
  _iDebJA4 = 0; _iNbJA4 = 0;
  _iDebJEx = 0; _iNbJEx = 0;
  _iNbInc = 0;

}
void A_LISTPARAM::setParam()
{
  int maxPass = (A_FILTER::_filterSetting._gps._use+A_FILTER::_filterSetting._glo._use+A_FILTER::_filterSetting._gal._use+A_FILTER::_filterSetting._bds._use) * A_PASS_MANAGEMENT::MaxPassSyst;
  _iDebX = 0; _iNbX = 1;
  _iDebY = _iDebX + _iNbX; _iNbY = 1;
  _iDebZ = _iDebY + _iNbY; _iNbZ = 1;
  _iDebT = _iDebZ + _iNbZ; _iNbT = 1;
  _iDebe = _iDebT + _iNbT; _iNbe = maxPass;
  _iDebHP.resize(SystMax);
  _iDebHl.resize(SystMax);
  _iNbHP.resize(SystMax);
  _iNbHl.resize(SystMax);
  for(int syst=0; syst<SystMax; syst++)
  {
   _iDebHP[syst].resize(FMax,0.0);
   _iDebHl[syst].resize(FMax,0.0);
   _iNbHP[syst].resize(FMax,0.0);
   _iNbHl[syst].resize(FMax,0.0);
  }
  _iDebHP[GPS][F1] = _iDebe + _iNbe; if (A_FILTER::_filterSetting._gps._use) _iNbHP[GPS][F1] = 1;
  _iDebHP[GPS][F2] = _iDebHP[GPS][F1] + _iNbHP[GPS][F1]; if (A_FILTER::_filterSetting._gps._use) _iNbHP[GPS][F2] = 1;
  _iDebHP[GPS][F5] = _iDebHP[GPS][F2] + _iNbHP[GPS][F2]; if (A_FILTER::_filterSetting._gps._use) _iNbHP[GPS][F5] = 1;
  _iDebHl[GPS][F1] = _iDebHP[GPS][F5] + _iNbHP[GPS][F5]; if (A_FILTER::_filterSetting._gps._use) _iNbHl[GPS][F1] = 1;
  _iDebHl[GPS][F2] = _iDebHl[GPS][F1] + _iNbHl[GPS][F1]; if (A_FILTER::_filterSetting._gps._use) _iNbHl[GPS][F2] = 1;
  _iDebHl[GPS][F5] = _iDebHl[GPS][F2] + _iNbHl[GPS][F2]; if (A_FILTER::_filterSetting._gps._use) _iNbHl[GPS][F5] = 1;
  _iDebHP[Glo][F1] = _iDebHl[GPS][F5] + _iNbHl[GPS][F5]; if (A_FILTER::_filterSetting._glo._use) _iNbHP[Glo][F1] = 1;
  _iDebHP[Glo][F2] = _iDebHP[Glo][F1] + _iNbHP[Glo][F1]; if (A_FILTER::_filterSetting._glo._use) _iNbHP[Glo][F2] = 1;
  _iDebHl[Glo][F1] = _iDebHP[Glo][F2] + _iNbHP[Glo][F2]; if (A_FILTER::_filterSetting._glo._use) _iNbHl[Glo][F1] = 1;
  _iDebHl[Glo][F2] = _iDebHl[Glo][F1] + _iNbHl[Glo][F1]; if (A_FILTER::_filterSetting._glo._use) _iNbHl[Glo][F2] = 1;
  _iDebHP[Gal][F1] = _iDebHl[Glo][F2] + _iNbHl[Glo][F2]; if (A_FILTER::_filterSetting._gal._use) _iNbHP[Gal][F1] = 1;
  _iDebHP[Gal][F5] = _iDebHP[Gal][F1] + _iNbHP[Gal][F1]; if (A_FILTER::_filterSetting._gal._use) _iNbHP[Gal][F5] = 1;
  _iDebHP[Gal][F7] = _iDebHP[Gal][F5] + _iNbHP[Gal][F5]; if (A_FILTER::_filterSetting._gal._use) _iNbHP[Gal][F7] = 1;
  _iDebHl[Gal][F1] = _iDebHP[Gal][F7] + _iNbHP[Gal][F7]; if (A_FILTER::_filterSetting._gal._use) _iNbHl[Gal][F1] = 1;
  _iDebHl[Gal][F5] = _iDebHl[Gal][F1] + _iNbHl[Gal][F1]; if (A_FILTER::_filterSetting._gal._use) _iNbHl[Gal][F5] = 1;
  _iDebHl[Gal][F7] = _iDebHl[Gal][F5] + _iNbHl[Gal][F5]; if (A_FILTER::_filterSetting._gal._use) _iNbHl[Gal][F7] = 1;
  _iDebHP[Bds][F1] = _iDebHl[Gal][F7] + _iNbHl[Gal][F7]; if (A_FILTER::_filterSetting._bds._use) _iNbHP[Bds][F1] = 1;
  _iDebHP[Bds][F6] = _iDebHP[Bds][F1] + _iNbHP[Bds][F1]; if (A_FILTER::_filterSetting._bds._use) _iNbHP[Bds][F6] = 1;
  _iDebHP[Bds][F7] = _iDebHP[Bds][F6] + _iNbHP[Bds][F6]; if (A_FILTER::_filterSetting._bds._use) _iNbHP[Bds][F7] = 1;
  _iDebHl[Bds][F1] = _iDebHP[Bds][F7] + _iNbHP[Bds][F7]; if (A_FILTER::_filterSetting._bds._use) _iNbHl[Bds][F1] = 1;
  _iDebHl[Bds][F6] = _iDebHl[Bds][F1] + _iNbHl[Bds][F1]; if (A_FILTER::_filterSetting._bds._use) _iNbHl[Bds][F6] = 1;
  _iDebHl[Bds][F7] = _iDebHl[Bds][F6] + _iNbHl[Bds][F6]; if (A_FILTER::_filterSetting._bds._use) _iNbHl[Bds][F7] = 1;
  _iDebP = _iNbT  + _iNbe + _iNbX + _iNbY + _iNbZ + \
           _iNbHP[GPS][F1] + _iNbHP[GPS][F2] + _iNbHP[GPS][F5] + \
           _iNbHl[GPS][F1] + _iNbHl[GPS][F2] + _iNbHl[GPS][F5] + \
           _iNbHP[Glo][F1] + _iNbHP[Glo][F2] + \
           _iNbHl[Glo][F1] + _iNbHl[Glo][F2] + \
           _iNbHP[Gal][F1] + _iNbHP[Gal][F5] + _iNbHP[Gal][F7] + \
           _iNbHl[Gal][F1] + _iNbHl[Gal][F5] + _iNbHl[Gal][F7] + \
           _iNbHP[Bds][F1] + _iNbHP[Bds][F6] + _iNbHP[Bds][F7] + \
           _iNbHl[Bds][F1] + _iNbHl[Bds][F6] + _iNbHl[Bds][F7];
  _iNbP = maxPass;
  _iDebA4 = _iDebP + _iNbP; _iNbA4 = maxPass;
  _iDebEx = _iDebA4 + _iNbA4; _iNbEx = maxPass;
  _iDebJP = _iDebEx + _iNbEx; _iNbJP = maxPass;
  _iDebJA4 = _iDebJP + _iNbJP; _iNbJA4 = maxPass;
  _iDebJEx = _iDebJA4 + _iNbJA4; _iNbJEx = maxPass;
  _iNbInc = _iNbHP[GPS][F1] + _iNbHP[GPS][F2] + _iNbHP[GPS][F5] + \
            _iNbHl[GPS][F1] + _iNbHl[GPS][F2] + _iNbHl[GPS][F5] + \
	    _iNbHP[Glo][F1] + _iNbHP[Glo][F2] + \
            _iNbHl[Glo][F1] + _iNbHl[Glo][F2] + \
	    _iNbHP[Gal][F1] + _iNbHP[Gal][F5] + _iNbHP[Gal][F7] + \
            _iNbHl[Gal][F1] + _iNbHl[Gal][F5] + _iNbHl[Gal][F7] + \
	    _iNbHP[Bds][F1] + _iNbHP[Bds][F6] + _iNbHP[Bds][F7] + \
            _iNbHl[Bds][F1] + _iNbHl[Bds][F6] + _iNbHl[Bds][F7] + \
	    _iNbT  + _iNbP + _iNbA4 + _iNbEx + _iNbJP + _iNbJA4 + \
	    _iNbJEx + _iNbe + _iNbX + _iNbY + _iNbZ;
}

//////////////////////////////////////////////////////////////////////
int A_LISTPARAM::getIndex(const ParamNature param, const int syst, const int pass) const
{
  int index = -1;
  int decalSyst = -1;
  int minConst= A_PASS_MANAGEMENT::MaxPassSyst + 1;
  
  if ((syst != -1) && (pass != -1))
  {
    switch(syst)
    {
      case GPS :
        decalSyst = 0;
        break;
      case Glo :
        if (A_FILTER::_filterSetting._gps._use)
          decalSyst = syst * (minConst-1);
        else
          decalSyst = 0;
        break;
      case Gal :
        if (A_FILTER::_filterSetting._gps._use && A_FILTER::_filterSetting._glo._use)
          decalSyst = syst * (minConst-1);
        else
        {
          int nbSyst=A_FILTER::_filterSetting._gps._use + A_FILTER::_filterSetting._glo._use;
          decalSyst = nbSyst * (minConst-1);
        }
        break;
      case Bds :
        if (A_FILTER::_filterSetting._gps._use && A_FILTER::_filterSetting._glo._use && A_FILTER::_filterSetting._gal._use)
          decalSyst = syst * (minConst-1);
        else
        {
          int nbSyst=A_FILTER::_filterSetting._gps._use + A_FILTER::_filterSetting._glo._use + A_FILTER::_filterSetting._gal._use;
          decalSyst = nbSyst * (minConst-1);
        }
        break;
      default :
        decalSyst = -1;
    }

  }

  switch(param)
  {
    case X :
      index = _iDebX;
      break;
    case Y :
      index = _iDebY;
      break;
    case Z :
      index = _iDebZ;
      break;
    case T :
      index = _iDebT;
      break;
    case E :
      if (decalSyst != -1)
        index = _iDebe + decalSyst + pass;
      else
        index = _iDebe;
      break;
    case P1 :
      index = _iDebHP[syst][F1];
      break;
    case P2 :
      index = _iDebHP[syst][F2];
      break;
    case C5 :
      index = _iDebHP[syst][F5];
      break;
    case C6 :
      index = _iDebHP[syst][F6];
      break;
    case C7 :
      index = _iDebHP[syst][F7];
      break;
    case L1 :
      index = _iDebHl[syst][F1];
      break;
    case L2 :
      index = _iDebHl[syst][F2];
      break;
    case L5 :
      index = _iDebHl[syst][F5];
      break;
    case L6 :
      index = _iDebHl[syst][F6];
      break;
    case L7 :
      index = _iDebHl[syst][F7];    
      break;
    case P :
      if (decalSyst != -1)
        index = _iDebP + decalSyst + pass;
      else
        index = _iDebP;
      break;
    case A4 :
      if (decalSyst != -1)
        index = _iDebA4 + decalSyst + pass;
      else
        index = _iDebA4;
      break;
    case Ex :
      if (decalSyst != -1)
        index = _iDebEx + decalSyst + pass;
      else
        index = _iDebEx;
      break;
    case JP :
      if (decalSyst != -1)
        index = _iDebJP + decalSyst + pass;
      else
        index = _iDebJP;
      break;
    case JA4 :
      if (decalSyst != -1)
        index = _iDebJA4 + decalSyst + pass;
      else
        index = _iDebJA4;
      break;
    case JEx :
      if (decalSyst != -1)
	index = _iDebJEx + decalSyst + pass;
      else
        index = _iDebJEx;
      break;
    default :
      index = -1;    
  }
  return index;
}

//////////////////////////////////////////////////////////////////////
int A_LISTPARAM::getSize(const ParamNature param, const int syst) const
{
  int size = -1;

  switch(param)
  {
    case X :
      size = _iNbX;
      break;
    case Y :
      size = _iNbY;
      break;
    case Z :
      size = _iNbZ;
      break;
    case T :
      size = _iNbT;
      break;
    case E :
      size = _iNbe;
      break;
    case P1 :
      size = _iNbHP[syst][F1];
      break;
    case P2 :
      size = _iNbHP[syst][F2];
      break;
    case C5 :
      size = _iNbHP[syst][F5];
      break;
    case C6 :
      size = _iNbHP[Bds][F6];
      break;
    case C7 :
      size = _iNbHP[syst][F7];
      break;
    case L1 :
      size = _iNbHl[syst][F1];
      break;
    case L2 :
      size = _iNbHl[syst][F2];
      break;
    case L5 :
      size = _iNbHl[syst][F5];
      break;
    case L6 :
      size = _iNbHl[Bds][F6];
      break;
    case L7 :     
      size = _iNbHl[syst][F7];
      break;
    case P :
      size = _iNbP;
      break;
    case A4 :
      size = _iNbA4;
      break;
    case Ex :
      size = _iNbEx;
      break;
    case JP :
      size = _iNbJP;
      break;
    case JA4 :
      size = _iNbJA4;
      break;
    case JEx :
      size = _iNbJEx;
      break;
    default :
      size = -1;    
  }
  return size;
}

//////////////////////////////////////////////////////////////////////
A_PASS_MANAGEMENT::A_PASS_MANAGEMENT()
{
  tabPassToSat = std::vector< std::vector<int> > (SystMax);
  tabSatToPass = std::vector< std::vector<int> > (SystMax);
  init();
}

//////////////////////////////////////////////////////////////////////
void A_PASS_MANAGEMENT::init()
{
  for (int isyst=0;isyst<SystMax;isyst++)
  {
    tabPassToSat[isyst] = std::vector<int>(MaxPassSyst,-1);
    tabSatToPass[isyst] = std::vector<int>(getMaxSatSyst(isyst),-1);
  }
}

//////////////////////////////////////////////////////////////////////
int A_PASS_MANAGEMENT::beginPass(int isyst, int isat)
{
  int npass = -1;
  for (int ipass=0;ipass<MaxPassSyst;ipass++)
  {
    if (tabPassToSat[isyst][ipass] == -1)
    {
      tabPassToSat[isyst][ipass] = isat;
      tabSatToPass[isyst][isat] = ipass;
      npass = ipass;
      break;
    }
  }
  return npass;
}

//////////////////////////////////////////////////////////////////////
void A_PASS_MANAGEMENT::endPass(int isyst, int ipass)
{
  int isat=tabPassToSat[isyst][ipass];
  if (isat != -1)
    tabSatToPass[isyst][isat] = -1;
  tabPassToSat[isyst][ipass] = -1;
}

//////////////////////////////////////////////////////////////////////
int A_PASS_MANAGEMENT::numPass(int isyst, int isat)
{
  return tabSatToPass[isyst][isat];
}

//////////////////////////////////////////////////////////////////////
void A_FILTER::fixAmbiguities(const AmbiguityNature an, const std::vector<std::vector<A_RESIDUAL> >& residuals, bool init)
{
  int BlocSta;
  int iSat_xMin, iPassMin;
  double y, z, dxMin, yry;

  for (int iSyst=0; iSyst<SystMax; iSyst++)
  {
    if ((iSyst == GPS && _filterSetting._gps._use) || (iSyst == Glo && _filterSetting._glo._use) || (iSyst == Gal && _filterSetting._gal._use) || (iSyst == Bds && _filterSetting._bds._use))    
    {
      ///////////////////////////////////////
      // Filter on the number os measurements
      int nbSat=0;
      for (int isat=0;isat<getMaxSatSyst(iSyst);isat++)
      {
        const A_RESIDUAL &residualsSat = residuals[iSyst][isat];
	bool process=((an == EX) && residualsSat._ewli) || ((an == NW) && residualsSat._wli) || ((an == N1) && residualsSat._n1i);
	double xL = 0.0;
	switch (iSyst)
	{
	  case GPS :
	    xL=(an==EX?residualsSat._resPhase[F5]:(an==NW?residualsSat._resPhase[F2]:residualsSat._resPhase[F1]));
	    break;
	  case Glo :
	    xL=0.0;
	    break;
	  case Gal :
	    xL=(an==EX?residualsSat._resPhase[F7]:(an==NW?residualsSat._resPhase[F5]:residualsSat._resPhase[F1]));
	    break;
	  case Bds :
	    xL=(an==EX?residualsSat._resPhase[F7]:(an==NW?residualsSat._resPhase[F6]:residualsSat._resPhase[F1]));
	    break;
	}
	if (process && xL)
	  nbSat++;
      }

      if (nbSat < _filterSetting._nbSatFixAmb)
	continue;
      ///////////////////////////////////////
      
      std::vector<int> bUsed(getMaxSatSyst(iSyst),0); 
      BlocSta = 0;
      for (int isat=0;isat<getMaxSatSyst(iSyst);isat++)
      {
	if (an==EX?_ex[iSyst][isat]:(an==NW?_a4[iSyst][isat]:_n[iSyst][isat]))
	  BlocSta = 1;
      }

      iSat_xMin = -1;
      iPassMin = -1;
      while (1)
      {
	dxMin = 100.0;
	for (int isat=0;isat<getMaxSatSyst(iSyst);isat++)
	{
	  bool process=((an == EX) && residuals[iSyst][isat]._ewli) || ((an == NW) && residuals[iSyst][isat]._wli) || ((an == N1) && residuals[iSyst][isat]._n1i);
	  if (process)
	  {
            int ipass=_pass.numPass(iSyst, isat);
            if ((ipass != -1) && (_nbMeasPas[iSyst][isat] >= _filterSetting._nbMin) && (an==N1?_a4[iSyst][isat]:true))
            //if ((ipass != -1) && (_nbMeasPas[iSyst][isat] >= _filterSetting._nbMin) && (an==N1?_a4[iSyst][isat]:(an==NW?_ex[iSyst][isat]:true)))
            {
	      double x, xL, xx, xy, xz, xmap;
	      bool b;
	      _filterInt.recoverDiagonalCovariance(1+_params.getIndex((an==EX?Ex:(an==NW?A4:P)),iSyst,ipass), x);
              x = sqrt(fabs(x));
	      if (init)
	      {
        	const A_RESIDUAL &residualsSat = residuals[iSyst][isat];
		xmap = residualsSat._mapping;
		switch (iSyst)
		{
		  case GPS :
	            xL=(an==EX?residualsSat._resPhase[F5]:(an==NW?residualsSat._resPhase[F2]:residualsSat._resPhase[F1]));
	            break;
		  case Glo :
	            xL=0.0;
	            break;
		  case Gal :
	            xL=(an==EX?residualsSat._resPhase[F7]:(an==NW?residualsSat._resPhase[F5]:residualsSat._resPhase[F1]));
	            break;
		  case Bds :
	            xL=(an==EX?residualsSat._resPhase[F7]:(an==NW?residualsSat._resPhase[F6]:residualsSat._resPhase[F1]));
	            break;
		}
        	b = xmap && xL && !BlocSta && (x < dxMin) && !bUsed[isat];      
	      }
	      else
	      {
        	y = _solInt[_params.getIndex((an==EX?Ex:(an==NW?A4:P)),iSyst,ipass)];
        	yry = fabs(y-rint(y));
        	b = x && (x < _filterSetting._thrAmb) && (yry < THR_AMB) &&
        	       (!(an==EX?_ex[iSyst][isat]:(an==NW?_a4[iSyst][isat]:_n[iSyst][isat]))) && (x <= dxMin) && !bUsed[isat] && BlocSta;
	      }
	      if (b)
	      {
        	iPassMin = ipass;
        	iSat_xMin = isat;
        	dxMin = x;
	      }
            }
	  }
	}
	if (dxMin == 100.0)
	  break;  
	bUsed[iSat_xMin] = 1;
	y = _solInt[_params.getIndex((an==EX?Ex:(an==NW?A4:P)),iSyst,iPassMin)];
	if (init)
	{
	  z = rint(y); if (!z) z = EPS_FIX;
	}
	else
	{
	  z = 0;
	  yry = fabs(y-rint(y));
	  if (yry < THR_AMB)
	  {
              z = rint(y); if (!z) z = EPS_FIX;
	  }
	}
	if (z)
	{
	  if (!BlocSta)
	  {
            PRINT_LOG(_filterLog._log,"Rover=%d : Receiver clock initialization\n", _rover);
	    FFLUSH_LOG(_filterLog._log);
	    switch (iSyst)
	    {
	      case GPS :
        	_filterInt.updateDiagonalCovariance(1+_params.getIndex(L1,iSyst), SIG_INI_CLK);
        	_filterInt.updateDiagonalCovariance(1+_params.getIndex(L2,iSyst), SIG_INI_CLK);
        	_filterInt.updateDiagonalCovariance(1+_params.getIndex(L5,iSyst), SIG_INI_CLK);
		break;
	      case Glo :
        	_filterInt.updateDiagonalCovariance(1+_params.getIndex(L1,iSyst), SIG_INI_CLK);
        	_filterInt.updateDiagonalCovariance(1+_params.getIndex(L2,iSyst), SIG_INI_CLK);
		break;
	      case Gal :
        	_filterInt.updateDiagonalCovariance(1+_params.getIndex(L1,iSyst), SIG_INI_CLK);
        	_filterInt.updateDiagonalCovariance(1+_params.getIndex(L5,iSyst), SIG_INI_CLK);
        	_filterInt.updateDiagonalCovariance(1+_params.getIndex(L7,iSyst), SIG_INI_CLK);
		break;
	      case Bds :
        	_filterInt.updateDiagonalCovariance(1+_params.getIndex(L1,iSyst), SIG_INI_CLK);
        	_filterInt.updateDiagonalCovariance(1+_params.getIndex(L6,iSyst), SIG_INI_CLK);
        	_filterInt.updateDiagonalCovariance(1+_params.getIndex(L7,iSyst), SIG_INI_CLK);
		break;
	    }
            BlocSta = 1;
	  }
	  if (an==EX)
	    PRINT_LOG(_filterLog._log,"Fix EX : Rover=%d, SAT=%1d,%02d, _N=%.0f, dxMin=%.03f, %d\n", _rover, iSyst, iSat_xMin+1, z, dxMin, init);
	  if (an==NW)
            PRINT_LOG(_filterLog._log,"Fix NW : Rover=%d, SAT=%1d,%02d, _N=%.0f, dxMin=%.03f, %d\n", _rover, iSyst, iSat_xMin+1, z, dxMin, init);
	  if (an==N1)
            PRINT_LOG(_filterLog._log,"Fix N1 : Rover=%d, SAT=%1d,%02d, _N=%.0f, dxMin=%.03f, %d\n", _rover, iSyst, iSat_xMin+1, z, dxMin, init);
	  FFLUSH_LOG(_filterLog._log);
	  double &v=(an==EX?_ex[iSyst][iSat_xMin]:(an==NW?_a4[iSyst][iSat_xMin]:_n[iSyst][iSat_xMin]));
	  v=z;
#define SIG_MES_AMB (0.0000001*0.0000001)
	  int iInc;
	  std::vector<double> _sol;
	  std::vector<A_FILTER_MODELED_MEASUREMENT> modeledMeasurement;
	  std::vector<double> der;
	  der.resize(_params.getNbInc(),0.0);  
	  der[(an==EX?_params.getIndex(Ex,iSyst,iPassMin):(an==NW?_params.getIndex(A4,iSyst,iPassMin):_params.getIndex(P,iSyst,iPassMin)))]=1.0;
	  A_FILTER_MODELED_MEASUREMENT a_modeledMeasurement(z,SIG_MES_AMB,der);
	  modeledMeasurement.push_back(a_modeledMeasurement);
	  for (iInc=0; iInc<_params.getNbInc(); iInc++)
            modeledMeasurement[0]._res -= _solInt[iInc]*modeledMeasurement[0]._der[iInc];
	  _filterInt.newMeasurementsElim(modeledMeasurement, 0.0, _sol);
	  for (iInc=0; iInc<_params.getNbInc(); iInc++)
            _solInt[iInc] += _sol[iInc];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
void A_FILTER::fixJumps(const AmbiguityNature an, const std::vector<std::vector<A_RESIDUAL> >& residuals,
                        const std::vector<std::vector<int> > &elimL, std::vector<std::vector<int> > &fixed, int isyst, bool init)
{
  std::vector<int> bUsed(getMaxSatSyst(isyst),0); 
  int iSat_xMin, iPassMin;
  double y, z, dxMin, yry;
  iSat_xMin = -1;
  iPassMin = -1;
  while (1)
  {
    dxMin = 100.0;
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      int ipass=_pass.numPass(isyst, isat);
      if (ipass != -1)
      {
        double x, xL;
        bool b;
        if (init)
        {
    	  const A_RESIDUAL &residualsSat = residuals[isyst][isat];
          x = residualsSat._mapping;
	  switch (isyst)
	  {
	    case GPS :
	      xL=(an==EX?residualsSat._resPhase[F5]:(an==NW?residualsSat._resPhase[F2]:residualsSat._resPhase[F1]));
	      break;
	    case Glo :
	      xL=(an==NW?residualsSat._resPhase[F2]:residualsSat._resPhase[F1]);
	      break;
	    case Gal :
	      xL=(an==EX?residualsSat._resPhase[F7]:(an==NW?residualsSat._resPhase[F5]:residualsSat._resPhase[F1]));
	      break;
	    case Bds :
	      xL=(an==EX?residualsSat._resPhase[F7]:(an==NW?residualsSat._resPhase[F6]:residualsSat._resPhase[F1]));
	      break;
	  }
    	  b = x && xL && (x < dxMin) && !bUsed[isat] && elimL[isyst][isat];	 
        }
        else
        {
    	  _filter.recoverDiagonalCovariance(1+_params.getIndex((an==EX?JEx:(an==NW?JA4:JP)),isyst,ipass), x);
    	  x = sqrt(fabs(x));
    	  y = _sol[_params.getIndex((an==EX?JEx:(an==NW?JA4:JP)),isyst,ipass)];
    	  yry = fabs(y-rint(y));
    	  b = (x > 0.001) && (yry < THR_AMB) && (x <= dxMin) && !bUsed[isat] && elimL[isyst][isat];
        }
        if (b)
        {
          iSat_xMin = isat;
          iPassMin = ipass;
          dxMin = x;
        }
      }
    }
    if (dxMin == 100.0)
      break;  
    bUsed[iSat_xMin] = 1;
    y = _sol[_params.getIndex((an==EX?JEx:(an==NW?JA4:JP)),isyst,iPassMin)];
    if (init)
    {
      z = rint(y); if (!z) z = EPS_FIX;
    }
    else
    {
      z = 0;
      yry = fabs(y-rint(y));
      if (yry < THR_AMB)
      {
          z = rint(y); if (!z) z = EPS_FIX;
      }
    }
    if (z)
    {
      if (an==EX)
	PRINT_LOG(_filterLog._log,"Fix Jump EX : Rover=%d, SAT=%1d,%02d, _N=%.0f, dxMin=%.03f, %d\n", _rover, isyst, iSat_xMin+1, z, dxMin, init);
      if (an==NW)
        PRINT_LOG(_filterLog._log,"Fix Jump NW : Rover=%d, SAT=%1d,%02d, _N=%.0f, dxMin=%.03f, %d\n", _rover, isyst, iSat_xMin+1, z, dxMin, init);
      if (an==N1)
        PRINT_LOG(_filterLog._log,"Fix Jump N1 : Rover=%d, SAT=%1d,%02d, _N=%.0f, dxMin=%.03f, %d\n", _rover, isyst, iSat_xMin+1, z, dxMin, init);
      FFLUSH_LOG(_filterLog._log);
      
      double &v=_sol[_params.getIndex((an==EX?JEx:(an==NW?JA4:JP)),isyst,iPassMin)];
      v=z;
#define SIG_MES_AMB (0.0000001*0.0000001)
      int iInc;
      std::vector<double> _s;
      std::vector<A_FILTER_MODELED_MEASUREMENT> modeledMeasurement;
      std::vector<double> der;
      der.resize(_params.getNbInc(),0.0);  
      der[(an==EX?_params.getIndex(JEx,isyst,iPassMin):(an==NW?_params.getIndex(JA4,isyst,iPassMin):_params.getIndex(JP,isyst,iPassMin)))]=1.0;
      A_FILTER_MODELED_MEASUREMENT a_modeledMeasurement(v,SIG_MES_AMB,der);
      modeledMeasurement.push_back(a_modeledMeasurement);
      for (iInc=0; iInc<_params.getNbInc(); iInc++)
        modeledMeasurement[0]._res -= _sol[iInc]*modeledMeasurement[0]._der[iInc];
      _filter.newMeasurementsElim(modeledMeasurement, 0.0, _s);
      for (iInc=0; iInc<_params.getNbInc(); iInc++)
        _sol[iInc] += _s[iInc];
      fixed[isyst][iSat_xMin]=1;
      if (init)
        break;
    }
  }
}

//////////////////////////////////////////////////////////////////////
void A_FILTER::computeMeasurementCode(const Frequency f, const std::vector<std::vector<A_RESIDUAL> >& residual,
                                    std::vector<std::vector<A_FILTER_MODELED_MEASUREMENT> >& measurements,
                                    std::vector<std::vector<int> >& satMes)
{
  
  A_FILTER_MODELED_MEASUREMENT modeled_measurement;
  
  for (int isyst=0; isyst<SystMax; isyst++)
  {
      measurements[isyst].clear();
      satMes[isyst].clear();
  }
  
  for (int isyst=0; isyst<SystMax; isyst++)
  {
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      const A_RESIDUAL &residualSat = residual[isyst][isat];
      int ipass=_pass.numPass(isyst, isat);
      if ((ipass != -1) && residualSat._resCode[(int)f] && residualSat._mapping && residualSat._pos.getX() && residualSat._pos.getY() && residualSat._pos.getZ() && _nbMeasPas[isyst][isat])
      {
	modeled_measurement._der.clear();
	modeled_measurement._der.resize(_params.getNbInc(),0.0);
	if (isyst == GPS)
	  modeled_measurement._w = _filterSetting._gps._sigCode[0]*residualSat._mapping;
	else if (isyst == Glo)
	  modeled_measurement._w = _filterSetting._glo._sigCode[0]*residualSat._mapping;
	else if (isyst == Gal)
	  modeled_measurement._w = _filterSetting._gal._sigCode[0]*residualSat._mapping;
	else if (isyst == Bds)
	{
	  if (residualSat._slot_typeBds == GEO)
	    modeled_measurement._w = _filterSetting._bds._sigCode[0]*residualSat._mapping;
	  else if  (residualSat._slot_typeBds == IGSO)
	    modeled_measurement._w = _filterSetting._bds._sigCode[1]*residualSat._mapping;
	  else if  (residualSat._slot_typeBds == MEO)
	    modeled_measurement._w = _filterSetting._bds._sigCode[2]*residualSat._mapping;
	}

	switch (f)
	{
	  case F1 :
	    modeled_measurement._der[_params.getIndex(P1,isyst)]=1.0;
	    modeled_measurement._der[_params.getIndex(E,isyst,ipass)]=1.0;
	    break;
	  case F2 :
	    modeled_measurement._der[_params.getIndex(P2,isyst)]=1.0;
	    modeled_measurement._der[_params.getIndex(P1,isyst)]=1.0;	
	    modeled_measurement._der[_params.getIndex(E,isyst,ipass)]=FREQUENCY::GAMMA(F2, isyst, residualSat._slot_typeBds);
	    break;
	  case F5 :
	    modeled_measurement._der[_params.getIndex(C5,isyst)]=1.0;
	    modeled_measurement._der[_params.getIndex(P1,isyst)]=1.0;	
	    modeled_measurement._der[_params.getIndex(E,isyst,ipass)]=FREQUENCY::GAMMA(F5, isyst, residualSat._slot_typeBds);
	    break;
	  case F6 :
	    modeled_measurement._der[_params.getIndex(C6,isyst)]=1.0;
	    modeled_measurement._der[_params.getIndex(P1,isyst)]=1.0;	
	    modeled_measurement._der[_params.getIndex(E,isyst,ipass)]=FREQUENCY::GAMMA(F6, isyst, residualSat._slot_typeBds);
	    break;
	  case F7 :
	    modeled_measurement._der[_params.getIndex(C7,isyst)]=1.0;
	    modeled_measurement._der[_params.getIndex(P1,isyst)]=1.0;
	    modeled_measurement._der[_params.getIndex(E,isyst,ipass)]=FREQUENCY::GAMMA(F7, isyst, residualSat._slot_typeBds);
	    break;
	}

	modeled_measurement._der[_params.getIndex(T)]=residualSat._mapping;
	modeled_measurement._der[_params.getIndex(X)]=residualSat._pos.getX();
	modeled_measurement._der[_params.getIndex(Y)]=residualSat._pos.getY();
	modeled_measurement._der[_params.getIndex(Z)]=residualSat._pos.getZ();
	modeled_measurement._res = residualSat._resCode[(int)f];

	measurements[isyst].push_back(modeled_measurement);
	satMes[isyst].push_back(isat);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
void A_FILTER::computeMeasurementPhase(const Frequency f, const std::vector<std::vector<A_RESIDUAL> >& residual,
                                    std::vector<std::vector<A_FILTER_MODELED_MEASUREMENT> >& measurements,
                                    std::vector<std::vector<int> >& satMes)
{
  
  A_FILTER_MODELED_MEASUREMENT modeled_measurement;
  
  for (int isyst=0; isyst<SystMax; isyst++)
  {
      measurements[isyst].clear();
      satMes[isyst].clear();
  }
  
  for (int isyst=0; isyst<SystMax; isyst++)
  {
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      const A_RESIDUAL &residualSat = residual[isyst][isat];
      int ipass=_pass.numPass(isyst, isat);
      if ((ipass != -1) && residualSat._resPhase[int(f)] && residualSat._mapping && residualSat._pos.getX() && residualSat._pos.getY() && residualSat._pos.getZ() && _nbMeasPas[isyst][isat])
      {
	modeled_measurement._der.clear();
	modeled_measurement._der.resize(_params.getNbInc(),0.0);
	if (isyst==GPS)
	  modeled_measurement._w = _filterSetting._gps._sigPhase[0]*residualSat._mapping;
	else if (isyst == Glo)
	  modeled_measurement._w = _filterSetting._glo._sigPhase[0]*residualSat._mapping;
	else if (isyst == Gal)
	  modeled_measurement._w = _filterSetting._gal._sigPhase[0]*residualSat._mapping;
	else if (isyst == Bds)
	{
	  if (residualSat._slot_typeBds == GEO)
	    modeled_measurement._w = _filterSetting._bds._sigPhase[0]*residualSat._mapping;
	  else if  (residualSat._slot_typeBds == IGSO)
	    modeled_measurement._w = _filterSetting._bds._sigPhase[1]*residualSat._mapping;
	  else if  (residualSat._slot_typeBds == MEO)
	    modeled_measurement._w = _filterSetting._bds._sigPhase[2]*residualSat._mapping;
	}


	switch (f)
	{
	  case F1 :
	    modeled_measurement._der[_params.getIndex(L1,isyst)]=1.0;
            modeled_measurement._der[_params.getIndex(P1,isyst)]=1.0;	
	    modeled_measurement._der[_params.getIndex(E,isyst,ipass)]=-1.0;
	    break;
	  case F2 :
	    modeled_measurement._der[_params.getIndex(L2,isyst)]=1.0;
	    modeled_measurement._der[_params.getIndex(P1,isyst)]=1.0;
	    modeled_measurement._der[_params.getIndex(E,isyst,ipass)]=-FREQUENCY::GAMMA(F2, isyst, residualSat._slot_typeBds);
	    if ((isyst == GPS) || (isyst == Glo))
	    {
	      modeled_measurement._der[_params.getIndex(A4,isyst,ipass)]=FREQUENCY::LAMBDA(F2, isyst, residualSat._slot_typeBds);
	      modeled_measurement._der[_params.getIndex(JA4,isyst,ipass)]=FREQUENCY::LAMBDA(F2, isyst, residualSat._slot_typeBds);
	    }
	    break;
	  case F5 :
	    modeled_measurement._der[_params.getIndex(L5,isyst)]=1.0;
	    modeled_measurement._der[_params.getIndex(P1,isyst)]=1.0;	
	    modeled_measurement._der[_params.getIndex(E,isyst,ipass)]=-FREQUENCY::GAMMA(F5, isyst, residualSat._slot_typeBds);
	    if ((isyst == Gal) || (isyst == GPS))
	    {
	      modeled_measurement._der[_params.getIndex(A4,isyst,ipass)]=FREQUENCY::LAMBDA(F5, isyst, residualSat._slot_typeBds);
	      modeled_measurement._der[_params.getIndex(JA4,isyst,ipass)]=FREQUENCY::LAMBDA(F5, isyst, residualSat._slot_typeBds);
	    }
	    if (isyst == GPS)
	    {
	      modeled_measurement._der[_params.getIndex(Ex,GPS,ipass)]=FREQUENCY::LAMBDA(F5, isyst, residualSat._slot_typeBds);
	      modeled_measurement._der[_params.getIndex(JEx,GPS,ipass)]=FREQUENCY::LAMBDA(F5, isyst, residualSat._slot_typeBds);
	    }
	    break;
	  case F6 :
	    modeled_measurement._der[_params.getIndex(L6,isyst)]=1.0;
	    modeled_measurement._der[_params.getIndex(P1,isyst)]=1.0;	
	    modeled_measurement._der[_params.getIndex(E,isyst,ipass)]=-FREQUENCY::GAMMA(F6, isyst, residualSat._slot_typeBds);
	    if (isyst == Bds)
	    {
	      modeled_measurement._der[_params.getIndex(A4,isyst,ipass)]=FREQUENCY::LAMBDA(F6, isyst, residualSat._slot_typeBds);
	      modeled_measurement._der[_params.getIndex(JA4,isyst,ipass)]=FREQUENCY::LAMBDA(F6, isyst, residualSat._slot_typeBds);
	    }
	    break;
	  case F7 :
	    modeled_measurement._der[_params.getIndex(L7,isyst)]=1.0;
	    modeled_measurement._der[_params.getIndex(P1,isyst)]=1.0;
	    modeled_measurement._der[_params.getIndex(E,isyst,ipass)]=-FREQUENCY::GAMMA(F7, isyst, residualSat._slot_typeBds);
	    if ((isyst == Gal) || (isyst == Bds))
	    {
	      modeled_measurement._der[_params.getIndex(A4,isyst,ipass)]=FREQUENCY::LAMBDA(F7, isyst, residualSat._slot_typeBds);
	      modeled_measurement._der[_params.getIndex(JA4,isyst,ipass)]=FREQUENCY::LAMBDA(F7, isyst, residualSat._slot_typeBds);
	      modeled_measurement._der[_params.getIndex(Ex,isyst,ipass)]=FREQUENCY::LAMBDA(F7, isyst, residualSat._slot_typeBds);
	      modeled_measurement._der[_params.getIndex(JEx,isyst,ipass)]=FREQUENCY::LAMBDA(F7, isyst, residualSat._slot_typeBds);
	    }
	    break;
	}
	modeled_measurement._der[_params.getIndex(P,isyst,ipass)]=FREQUENCY::LAMBDA(f, isyst, residualSat._slot_typeBds);
	modeled_measurement._der[_params.getIndex(JP,isyst,ipass)]=FREQUENCY::LAMBDA(f, isyst, residualSat._slot_typeBds);
	modeled_measurement._der[_params.getIndex(T)]=residualSat._mapping;
	modeled_measurement._der[_params.getIndex(X)]=residualSat._pos.getX();
	modeled_measurement._der[_params.getIndex(Y)]=residualSat._pos.getY();
	modeled_measurement._der[_params.getIndex(Z)]=residualSat._pos.getZ();
	modeled_measurement._res = FREQUENCY::LAMBDA(f, isyst, residualSat._slot_typeBds)*residualSat._resPhase[(int)f];

	if (modeled_measurement._w)
	{
	  measurements[isyst].push_back(modeled_measurement);
	  satMes[isyst].push_back(isat);
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
void A_FILTER::computeMeasurementIONO(const std::vector<std::vector<A_RESIDUAL> >& residual,
                                      std::vector<std::vector<A_FILTER_MODELED_MEASUREMENT> >& measurements,
                                      std::vector<std::vector<int> >& satMes)
{
 
  A_FILTER_MODELED_MEASUREMENT modeled_measurement;
  
  for (int isyst=0; isyst<SystMax; isyst++)
  {
    measurements[isyst].clear();
    satMes[isyst].clear();
  }

  for (int isyst=0; isyst<SystMax; isyst++)
  {
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      const A_RESIDUAL &residualSat=residual[isyst][isat];
      int ipass=_pass.numPass(isyst, isat);
      if ((ipass != -1) && residualSat._iono && _nbMeasPas[isyst][isat])
      {
	if(residualSat._typeIono  == 0)
	  modeled_measurement._w = _filterSetting._sigMeasIono[2];
	else	
	  modeled_measurement._w = _filterSetting._sigMeasIono[residualSat._typeIono - 1];

        modeled_measurement._der.clear();
        modeled_measurement._der.resize(_params.getNbInc(),0.0);
        modeled_measurement._der[_params.getIndex(E,isyst,ipass)]=1.0;
        modeled_measurement._res = residualSat._iono;

	if (modeled_measurement._w) {
	  measurements[isyst].push_back(modeled_measurement);
	  satMes[isyst].push_back(isat);
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
void A_FILTER::computeMeasurementTROPO(const double tropo, double *measurementsRes, int *measurementsElim)
{

  A_FILTER_MODELED_MEASUREMENT modeled_measurement;
  std::vector<A_FILTER_MODELED_MEASUREMENT> measurements;
  int iInc;
  std::vector<double> sol;
  *measurementsElim = 0;
  *measurementsRes = 0;

  if (!_filterSetting._sigMeasTropo)
    return;

  measurements.clear();

  modeled_measurement._w = _filterSetting._sigMeasTropo;
  modeled_measurement._der.clear();
  modeled_measurement._der.resize(_params.getNbInc(),0.0);
  modeled_measurement._der[_params.getIndex(T)]=1.0;
  modeled_measurement._res = tropo;

  measurements.push_back(modeled_measurement);

  for (iInc=0; iInc<_params.getNbInc(); iInc++)
    measurements[0]._res -= _sol[iInc]*measurements[0]._der[iInc];
  if (measurements[0]._w)
    *measurementsRes = measurements[0]._res;

  if (_filterSetting._raim)
    _filter.newMeasurementsElimComb(_filterSetting._maxElim, measurements, _filterSetting._thrMeasTropo, sol);
  else
    _filter.newMeasurementsElim(measurements, _filterSetting._thrMeasTropo, sol);
  for (iInc=0; iInc<_params.getNbInc(); iInc++)
    _sol[iInc] += sol[iInc];

  if (!measurements[0]._w)
  {
    PRINT_LOG(_filterLog._log, "TROPO rejection Rover=%d, residual : %8.2f\n", _rover, tropo);
    FFLUSH_LOG(_filterLog._log);
    *measurementsElim = 1;
  }
}

//////////////////////////////////////////////////////////////////////
void A_FILTER::eliminateMeasurement(const std::string strMeasure, std::vector<std::vector<double> >& measurementsRes, std::vector<std::vector<int> >& measurementsElim, 
                          const std::vector<std::vector<A_RESIDUAL> >& residual,
                          std::vector<std::vector<A_FILTER_MODELED_MEASUREMENT> >& measurements,
                          const std::vector<std::vector<int> >& satMes,
                          const double thr_meas,
                          std::vector<int>& reject)
{
  
  int iMes, iSat, iInc;
  std::vector<double> sol; 
  
  for (int isyst=0; isyst<SystMax; isyst++)
  {
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
    {
      reject[isyst]=0;
      std::vector<A_FILTER_MODELED_MEASUREMENT> &measurementsConst = measurements[isyst];
      const std::vector<int> &satMesConst = satMes[isyst];
      std::vector<double> &measurementsResConst = measurementsRes[isyst];
      std::vector<int> &measurementsElimConst = measurementsElim[isyst];
      for (iMes=0; iMes<(int)measurementsConst.size(); iMes++)
      {
	iSat = satMesConst[iMes];
	for (iInc=0; iInc<_params.getNbInc(); iInc++)
	  measurementsConst[iMes]._res -= _sol[iInc]*measurementsConst[iMes]._der[iInc];
	if (measurementsConst[iMes]._w)
	  measurementsResConst[iSat] = measurementsConst[iMes]._res;
      }

      if (_filterSetting._raim)
	_filter.newMeasurementsElimComb(_filterSetting._maxElim, measurementsConst, thr_meas, sol);
      else
	_filter.newMeasurementsElim(measurementsConst, thr_meas, sol);
      for (iInc=0; iInc<_params.getNbInc(); iInc++)
	_sol[iInc] += sol[iInc];

      for (iMes=0; iMes<(int)measurementsConst.size(); iMes++)
      {
	if (!measurementsConst[iMes]._w)
	{
          iSat = satMesConst[iMes];
          const A_RESIDUAL &residualSat = residual[isyst][iSat];
          PRINT_LOG(_filterLog._log,"%s rejection Rover=%d, SAT=%d,%d : %8.2f %8.2f %8.2f %8.2f %8.2f, Residual : %5.3f\n", strMeasure.c_str(),_rover, isyst,iSat+1, 
	             residualSat._resCode[F1],residualSat._resCode[F2],residualSat._resCode[F5],residualSat._resCode[F6],residualSat._resCode[F7], measurementsConst[iMes]._res);
          FFLUSH_LOG(_filterLog._log);

	  measurementsElimConst[iSat] = 1;
	  reject[isyst]++;
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
void A_FILTER::initStateVectorAmb(int isyst, int isat, int ipass, double value, bool jumps)
{
  _n[isyst][isat] = 0.0;
  _a4[isyst][isat] = 0.0;
  _ex[isyst][isat] = 0.0;
  _filter.updateDiagonalCovariance(1+_params.getIndex(P,isyst,ipass), value);
  _sol[_params.getIndex(P,isyst,ipass)] = 0.0;
  _filter.updateDiagonalCovariance(1+_params.getIndex(A4,isyst,ipass), value);
  _sol[_params.getIndex(A4,isyst,ipass)] = 0.0;
  _filter.updateDiagonalCovariance(1+_params.getIndex(Ex,isyst,ipass), value);
  _sol[_params.getIndex(Ex,isyst,ipass)] = 0.0;
  if (jumps)
  {
    _filter.updateDiagonalCovariance(1+_params.getIndex(JP,isyst,ipass), 0.0);
    _sol[_params.getIndex(JP,isyst,ipass)] = 0.0;
    _filter.updateDiagonalCovariance(1+_params.getIndex(JA4,isyst,ipass), 0.0);
    _sol[_params.getIndex(JA4,isyst,ipass)] = 0.0;
    _filter.updateDiagonalCovariance(1+_params.getIndex(JEx,isyst,ipass), 0.0);
    _sol[_params.getIndex(JEx,isyst,ipass)] = 0.0;
  }
}

static std::string selectMeasure(Measure measure,Frequency f)
{
  std::string stringMeas;
  switch (f)
  {
    case F1 :
      if (measure == pseudorange)
	stringMeas = "P1";
      else
	stringMeas = "L1";
      break;
    case F2 :
      if (measure == pseudorange)
	stringMeas = "P2";
      else
	stringMeas = "L2";
      break;
    case F5 :
      if (measure == pseudorange)
	stringMeas = "C5";
      else
	stringMeas = "L5";
      break;
    case F6 :
      if (measure == pseudorange)
	stringMeas = "C6";
      else
	stringMeas = "L6";
      break;
    case F7 :
      if (measure == pseudorange)
	stringMeas = "C7";
      else
	stringMeas = "L7";
      break;
  }
  return stringMeas;
}

//////////////////////////////////////////////////////////////////////
int A_FILTER::correctState(const A_DATE &date, const A_VECTOR3D& state, const A_VECTOR3D& point,
                          const std::vector<std::vector<A_RESIDUAL> >& residual,
			  const double tropo,
                          A_RECEIVER_PRODUCT &recProduct,
                          std::vector<std::vector<A_MEASUREMENT_PRODUCT> >& measurementProduct)
{

  
  int iInc;
  std::vector<std::vector<A_FILTER_MODELED_MEASUREMENT> > measurements(SystMax);
  std::vector<std::vector<int> > satMes(SystMax);
  std::vector<int> iReject(SystMax,0);
  int iRejectMax_GLO, iRejectMax_GPS, iRejectMax_GAL, iRejectMax_BDS;
  std::vector<std::vector<std::vector<int> > > elimCode(FMax, std::vector<std::vector<int> >(SystMax));
  std::vector<std::vector<std::vector<int> > > elimPhase(FMax, std::vector<std::vector<int> >(SystMax));
  std::vector<std::vector<int> > elimIONO(SystMax);
  std::vector<std::vector<std::vector<double> > > resCode(FMax, std::vector<std::vector<double> >(SystMax));
  std::vector<std::vector<std::vector<double> > > resPhase(FMax, std::vector<std::vector<double> >(SystMax));
  std::vector<std::vector<double> > resIONO(SystMax);
  double x;
  double resTropo;
  int elimTropo;
  std::string strMeasure;
   
  PRINT_LOG(_filterLog._log,"date %d %.1f\n", date.getDay(), date.getSecond());
  FFLUSH_LOG(_filterLog._log);
  //init vector
  for (int f=F1; f<FMax; f++)
  {
    for (int isyst=0; isyst<SystMax; isyst++)
    {
      if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
      {
        elimCode[f][isyst] = std::vector<int>(getMaxSatSyst(isyst),0);
        elimPhase[f][isyst] = std::vector<int>(getMaxSatSyst(isyst),0);
        resCode[f][isyst] = std::vector<double>(getMaxSatSyst(isyst),0.0);
        resPhase[f][isyst] = std::vector<double>(getMaxSatSyst(isyst),0.0);
      }
    }
  }
  
  for (int isyst=0; isyst<SystMax; isyst++)
  {
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
    {
      measurements[isyst] = std::vector<A_FILTER_MODELED_MEASUREMENT>(getMaxSatSyst(isyst));
      satMes[isyst] = std::vector<int>(getMaxSatSyst(isyst),0);
      elimIONO[isyst] = std::vector<int>(getMaxSatSyst(isyst),0);
      resIONO[isyst] = std::vector<double>(getMaxSatSyst(isyst),0.0);
    }
  }

  // Initial state (linearization around point)
  _sol[_params.getIndex(X)] = state.getX()-point.getX();
  _sol[_params.getIndex(Y)] = state.getY()-point.getY();
  _sol[_params.getIndex(Z)] = state.getZ()-point.getZ();


  for (int f=F1; f<FMax; f++)
  {
    strMeasure = selectMeasure(pseudorange,(Frequency) f);
    computeMeasurementCode((Frequency)f,residual,measurements,satMes);
    eliminateMeasurement(strMeasure, resCode[f], elimCode[f], residual, measurements, satMes, _filterSetting._thrMeasCode, iReject);
  }
  
  // Ionosphere (as a measurement)
  computeMeasurementIONO(residual,measurements,satMes);
  eliminateMeasurement("IONO",resIONO, elimIONO, residual, measurements, satMes, _filterSetting._thrMeasIono, iReject);

  // Troposphere (as a measurement)
  if (tropo)
  {
    double Longit,Latit,Altit;
    A_MODEL::geodesic(point.getX(), point.getY(), point.getZ(), &Longit, &Latit, &Altit);
    recProduct._tropoV = A_MODEL::tropDelay(1.0, Altit);
    computeMeasurementTROPO(tropo-recProduct._tropoV, &resTropo, &elimTropo);
  }

  // Initialization
  iRejectMax_GPS = 0;
  iRejectMax_GLO = 0;
  iRejectMax_GAL = 0;
  iRejectMax_BDS = 0;
  for (int f=F1; f<FMax; f++)
  {
    strMeasure = selectMeasure(carrierphase,(Frequency) f);
    computeMeasurementPhase((Frequency)f,residual,measurements,satMes);
    eliminateMeasurement(strMeasure, resPhase[f], elimPhase[f], residual, measurements, satMes, _filterSetting._thrMeasPhase, iReject);
    if (iReject[GPS] > iRejectMax_GPS)
      iRejectMax_GPS = iReject[GPS] ;
    if (iReject[Glo] > iRejectMax_GLO)
      iRejectMax_GLO = iReject[Glo] ;
    if (iReject[Gal] > iRejectMax_GAL)
      iRejectMax_GAL = iReject[Gal] ;
    if (iReject[Bds] > iRejectMax_BDS)
      iRejectMax_BDS = iReject[Bds] ;
  }

  // Full or partial reinitialization test
  bool jump[SystMax];
  bool nbJump = 0;
  for (int isyst=0; isyst<SystMax; isyst++)
  {
    jump[isyst] = false;
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
    for (int isat=0; isat<getMaxSatSyst(isyst); isat++)
    {
      if ( (iRejectMax_GPS>_filterSetting._maxElim && (resPhase[F1][isyst][isat] || resPhase[F2][isyst][isat] || resPhase[F5][isyst][isat]) && isyst==GPS)
        || (iRejectMax_GLO>_filterSetting._maxElim && (resPhase[F1][isyst][isat] || resPhase[F2][isyst][isat]) && isyst==Glo)
        || (iRejectMax_GAL>_filterSetting._maxElim && (resPhase[F1][isyst][isat] || resPhase[F5][isyst][isat] || resPhase[F7][isyst][isat]) && isyst==Gal)
        || (iRejectMax_BDS>_filterSetting._maxElim && (resPhase[F1][isyst][isat] || resPhase[F6][isyst][isat] || resPhase[F7][isyst][isat]) && isyst==Bds))
      {
        for (int f=F1; f<FMax; f++)
          elimPhase[f][isyst][isat] = 1;
        jump[isyst] = true;
      }
    }
    if (jump[isyst])
      nbJump++;
  }
  if (iRejectMax_GPS>_filterSetting._maxElim)
  {
    PRINT_LOG(_filterLog._log, "Rover=%d : too many rejects: GPS ambiguities supressions\n", _rover);
    FFLUSH_LOG(_filterLog._log);
  }
  if (iRejectMax_GLO>_filterSetting._maxElim)
  {
    PRINT_LOG(_filterLog._log, "Rover=%d : too many rejects: GLONASS ambiguities supressions\n", _rover);
    FFLUSH_LOG(_filterLog._log);
  }
  if (iRejectMax_GAL>_filterSetting._maxElim)
  {
    PRINT_LOG(_filterLog._log, "Rover=%d : too many rejects: GALILEO ambiguities supressions\n", _rover);
    FFLUSH_LOG(_filterLog._log);
  }
  if (iRejectMax_BDS>_filterSetting._maxElim)
  {
    PRINT_LOG(_filterLog._log, "Rover=%d : too many rejects: BEIDOU ambiguities supressions\n", _rover);
    FFLUSH_LOG(_filterLog._log);
  }
  if ( ((_filterSetting._gps._use == 0) || (iRejectMax_GPS>_filterSetting._maxElim))
    && ((_filterSetting._glo._use == 0) || (iRejectMax_GLO>_filterSetting._maxElim))
    && ((_filterSetting._gal._use == 0) || (iRejectMax_GAL>_filterSetting._maxElim)) 
    && ((_filterSetting._bds._use == 0) || (iRejectMax_BDS>_filterSetting._maxElim)))
  {
    PRINT_LOG(_filterLog._log,"Rover=%d : too many rejects: global ambiguities supressions\n", _rover);
    FFLUSH_LOG(_filterLog._log);
    for (iInc=_params.getIndex(X); iInc<_params.getIndex(X)+_params.getSize(X); iInc++)
      _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniPos);
    for (iInc=_params.getIndex(Y); iInc<_params.getIndex(Y)+_params.getSize(Y); iInc++)
      _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniPos);
    for (iInc=_params.getIndex(Z); iInc<_params.getIndex(Z)+_params.getSize(Z); iInc++)
      _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniPos);
  }
  
  // General phase jumps
  if (nbJump)
  {
    std::vector<std::vector<int> > elimL(SystMax);
    std::vector<std::vector<int> > fixedEx(SystMax);
    std::vector<std::vector<int> > fixedNw(SystMax);
    std::vector<std::vector<int> > fixedN1(SystMax);
    std::vector<std::vector<std::vector<int> > > elimJPhase(FMax,std::vector<std::vector<int> >(SystMax));
    //init vector
    for (int isyst=0; isyst<SystMax; isyst++)
    {
      if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
      {
        elimL[isyst] = std::vector<int>(getMaxSatSyst(isyst),0);
        fixedEx[isyst] = std::vector<int>(getMaxSatSyst(isyst),0);
        fixedNw[isyst] = std::vector<int>(getMaxSatSyst(isyst),0);
        fixedN1[isyst] = std::vector<int>(getMaxSatSyst(isyst),0);
      }
    }
    for (int f=F1; f<FMax; f++)
      for (int isyst=0; isyst<SystMax; isyst++)
        if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
          elimJPhase[f][isyst] = std::vector<int>(getMaxSatSyst(isyst),0);

    for (int isyst=0; isyst<SystMax; isyst++)
    {
      if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))    
        if (jump[isyst])
          for (int isat=0; isat<getMaxSatSyst(isyst); isat++)
          {
            int ipass=_pass.numPass(isyst, isat);
            if (ipass != -1)
            {
              if (elimPhase[F1][isyst][isat] || elimPhase[F2][isyst][isat] || elimPhase[F5][isyst][isat] || elimPhase[F6][isyst][isat] || elimPhase[F7][isyst][isat])
              {
                if (_filterSetting._fixN1) _filter.updateDiagonalCovariance(1+_params.getIndex(JP,isyst,ipass), SIG_INI_AMB);
                if (_filterSetting._fixNw) _filter.updateDiagonalCovariance(1+_params.getIndex(JA4,isyst,ipass), SIG_INI_AMB);
                if (_filterSetting._fixEx) _filter.updateDiagonalCovariance(1+_params.getIndex(JEx,isyst,ipass), SIG_INI_AMB);
                elimL[isyst][isat]=1;
              }
            }
          }  
    }
    for (int f=F1; f<FMax; f++)
    {
      strMeasure = selectMeasure(carrierphase,(Frequency) f);
      computeMeasurementPhase((Frequency)f,residual,measurements,satMes);
      eliminateMeasurement(strMeasure, resPhase[f], elimJPhase[f], residual, measurements, satMes, 0.0, iReject);
    }

    for (int isyst=0; isyst<SystMax; isyst++)
    {
      if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
      {
	if (_filterSetting._fixEx) fixJumps(EX, residual, elimL, fixedEx, isyst, true);
	if (_filterSetting._fixNw) fixJumps(NW, residual, elimL, fixedNw, isyst, true);
	if (_filterSetting._fixN1) fixJumps(N1, residual, elimL, fixedN1, isyst, true);
	if (_filterSetting._fixEx) fixJumps(EX, residual, elimL, fixedEx, isyst, false);  // Accumulation in fixed
	if (_filterSetting._fixNw) fixJumps(NW, residual, elimL, fixedNw, isyst, false);
	if (_filterSetting._fixN1) fixJumps(N1, residual, elimL, fixedN1, isyst, false);
	int notFixed=0;
	for (int isat=0; isat<getMaxSatSyst(isyst); isat++)
	{
          int ipass=_pass.numPass(isyst, isat);
          if (elimL[isyst][isat])
          {
	    bool bex, bnw, bn1;
	    bex=bnw=bn1=true;
	    if (isyst == GPS)
	    {
              if (_filterSetting._fixN1) bn1 = !resPhase[F1][isyst][isat] || (resPhase[F1][isyst][isat] && fixedN1[isyst][isat]);
              if (_filterSetting._fixNw) bnw = !resPhase[F2][isyst][isat] || (resPhase[F2][isyst][isat] && fixedNw[isyst][isat]);
              if (_filterSetting._fixEx) bex = !resPhase[F5][isyst][isat] || (resPhase[F5][isyst][isat] && fixedEx[isyst][isat]);
	    }
	    if (isyst == Glo)
	    {
              if (_filterSetting._fixN1) bn1 = !resPhase[F1][isyst][isat] || (resPhase[F1][isyst][isat] && fixedN1[isyst][isat]);
              if (_filterSetting._fixNw) bnw = !resPhase[F2][isyst][isat] || (resPhase[F2][isyst][isat] && fixedNw[isyst][isat]);
	    }
	    if (isyst == Gal)
	    {
              if (_filterSetting._fixN1) bn1 = !resPhase[F1][isyst][isat] || (resPhase[F1][isyst][isat] && fixedN1[isyst][isat]);
	      if (_filterSetting._fixNw) bnw = !resPhase[F5][isyst][isat] || (resPhase[F5][isyst][isat] && fixedNw[isyst][isat]);
	      if (_filterSetting._fixEx) bex = !resPhase[F7][isyst][isat] || (resPhase[F7][isyst][isat] && fixedEx[isyst][isat]);
	    }
	    if (isyst == Bds)
	    {
              if (_filterSetting._fixN1) bn1 = !resPhase[F1][isyst][isat] || (resPhase[F1][isyst][isat] && fixedN1[isyst][isat]);
	      if (_filterSetting._fixNw) bnw = !resPhase[F6][isyst][isat] || (resPhase[F6][isyst][isat] && fixedNw[isyst][isat]);
	      if (_filterSetting._fixEx) bex = !resPhase[F7][isyst][isat] || (resPhase[F7][isyst][isat] && fixedEx[isyst][isat]);
	    }
            if (!bn1 || !bnw || !bex)
            {
              PRINT_LOG(_filterLog._log, "Rover=%d : JUMP NOT FIXED %d %d\n", _rover, isyst, isat);
	      FFLUSH_LOG(_filterLog._log);
              initStateVectorAmb(isyst, isat, ipass, SIG_INI_AMB, false);
              notFixed++;
            }
          }
	}
	if (notFixed > _filterSetting._maxElim)	
	{
	return 1;
	}
	else
	{
	  if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
	  for (int isat=0; isat<getMaxSatSyst(isyst); isat++)
          {
            int ipass=_pass.numPass(isyst, isat);
            if (elimL[isyst][isat])
            {
	      bool bex, bnw, bn1;
	      bex=bnw=bn1=true;
	      if (isyst == GPS)
	      {
        	if (_filterSetting._fixN1) bn1 = !resPhase[F1][isyst][isat] || (resPhase[F1][isyst][isat] && fixedN1[isyst][isat]);
        	if (_filterSetting._fixNw) bnw = !resPhase[F2][isyst][isat] || (resPhase[F2][isyst][isat] && fixedNw[isyst][isat]);
        	if (_filterSetting._fixEx) bex = !resPhase[F5][isyst][isat] || (resPhase[F5][isyst][isat] && fixedEx[isyst][isat]);
	      }
	      if (isyst == Glo)
	      {
        	if (_filterSetting._fixN1) bn1 = !resPhase[F1][isyst][isat] || (resPhase[F1][isyst][isat] && fixedN1[isyst][isat]);
        	if (_filterSetting._fixNw) bnw = !resPhase[F2][isyst][isat] || (resPhase[F2][isyst][isat] && fixedNw[isyst][isat]);
	      }
	      if (isyst == Gal)
	      {
        	if (_filterSetting._fixN1) bn1 = !resPhase[F1][isyst][isat] || (resPhase[F1][isyst][isat] && fixedN1[isyst][isat]);
		if (_filterSetting._fixNw) bnw = !resPhase[F5][isyst][isat] || (resPhase[F5][isyst][isat] && fixedNw[isyst][isat]);
		if (_filterSetting._fixEx) bex = !resPhase[F7][isyst][isat] || (resPhase[F7][isyst][isat] && fixedEx[isyst][isat]);
	      }
	      if (isyst == Bds)
	      {
        	if (_filterSetting._fixN1) bn1 = !resPhase[F1][isyst][isat] || (resPhase[F1][isyst][isat] && fixedN1[isyst][isat]);
		if (_filterSetting._fixNw) bnw = !resPhase[F6][isyst][isat] || (resPhase[F6][isyst][isat] && fixedNw[isyst][isat]);
		if (_filterSetting._fixEx) bex = !resPhase[F7][isyst][isat] || (resPhase[F7][isyst][isat] && fixedEx[isyst][isat]);
	      }
	      if (bn1 && bnw && bex)
	      {
		elimPhase[F1][isyst][isat] = 0;
		elimPhase[F2][isyst][isat] = 0;
		elimPhase[F5][isyst][isat] = 0;
		elimPhase[F6][isyst][isat] = 0;
		elimPhase[F7][isyst][isat] = 0;
	      }
            }
	  }
	}
      }
    }
  }
  for (int isyst=0; isyst<SystMax; isyst++)
  {
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
    for (int isat=0; isat<getMaxSatSyst(isyst); isat++)
    {
      int ipass=_pass.numPass(isyst, isat);
      if (ipass != -1)
      {
        if (elimPhase[F1][isyst][isat] || elimPhase[F2][isyst][isat] || elimPhase[F5][isyst][isat] || elimPhase[F6][isyst][isat] || elimPhase[F7][isyst][isat])
        {
	  initStateVectorAmb(isyst, isat, ipass, SIG_INI_AMB, false);
	}
      }
    }  
  }

  for(int isyst=0; isyst<SystMax; isyst++)
  {
    _n[isyst].clear();
    _a4[isyst].clear();
    _ex[isyst].clear();
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
    {    
      _n[isyst].resize(getMaxSatSyst(isyst),0.0);
      _a4[isyst].resize(getMaxSatSyst(isyst),0.0);
      _ex[isyst].resize(getMaxSatSyst(isyst),0.0);
    }
  }
  _filterInt = _filter;
  for (iInc=0; iInc<_params.getNbInc(); iInc++)
    _solInt[iInc] = _sol[iInc];

  if (_filterSetting._fixEx) fixAmbiguities(EX, residual, true);
  if (_filterSetting._fixNw) fixAmbiguities(NW, residual, true);
  if (_filterSetting._fixN1) fixAmbiguities(N1, residual, true);
  if (_filterSetting._fixEx) fixAmbiguities(EX, residual, false);
  if (_filterSetting._fixNw) fixAmbiguities(NW, residual, false);
  if (_filterSetting._fixN1) fixAmbiguities(N1, residual, false);

  // Output
  if (_filterSetting._gps._use)
  {
    recProduct._HstaGPS = _solInt[_params.getIndex(P1,GPS)];
    _filterInt.recoverDiagonalCovariance(1+_params.getIndex(P1,GPS), x);
    recProduct._covHstaGPS = sqrt(fabs(x));
  }
  if (_filterSetting._glo._use)
  {
    recProduct._HstaGlo = _solInt[_params.getIndex(P1,Glo)];
    _filterInt.recoverDiagonalCovariance(1+_params.getIndex(P1,Glo), x);
    recProduct._covHstaGlo = sqrt(fabs(x));
  }
  if (_filterSetting._gal._use)
  {
    recProduct._HstaGal = _solInt[_params.getIndex(P1,Gal)];
    _filterInt.recoverDiagonalCovariance(1+_params.getIndex(P1,Gal), x);
    recProduct._covHstaGal = sqrt(fabs(x));
  }
  if (_filterSetting._bds._use)
  {
    recProduct._HstaBds = _solInt[_params.getIndex(P1,Bds)];
    _filterInt.recoverDiagonalCovariance(1+_params.getIndex(P1,Bds), x);
    recProduct._covHstaBds = sqrt(fabs(x));
  }
  recProduct._x = point.getX()+_solInt[_params.getIndex(X)];
  _filterInt.recoverDiagonalCovariance(1+_params.getIndex(X), x);
  recProduct._covX = sqrt(fabs(x));
  recProduct._y = point.getY()+_solInt[_params.getIndex(Y)];
  _filterInt.recoverDiagonalCovariance(1+_params.getIndex(Y), x);
  recProduct._covY = sqrt(fabs(x));
  recProduct._z = point.getZ()+_solInt[_params.getIndex(Z)];
  _filterInt.recoverDiagonalCovariance(1+_params.getIndex(Z), x);
  recProduct._covZ = sqrt(fabs(x));
  double Longit,Latit,Altit;
  A_MODEL::geodesic(recProduct._x, recProduct._y, recProduct._z, &Longit, &Latit, &Altit);
  recProduct._tropoV = A_MODEL::tropDelay(1.0, Altit);
  recProduct._tropo = _solInt[_params.getIndex(T)];
  _filterInt.recoverDiagonalCovariance(1+_params.getIndex(T), x);
  recProduct._covTropo = sqrt(fabs(x));
      
  recProduct._nbMes = 0;
  recProduct._nbMesBloEx = 0;
  recProduct._nbMesBloNw = 0;
  recProduct._nbMesBloN1 = 0;
  for (int isyst=0; isyst<SystMax; isyst++)
  {
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      if (_nbMeasPas[isyst][isat] && !_nbMesPas0[isyst][isat])
      {
	recProduct._nbMes++;
        if (_ex[isyst][isat])
          recProduct._nbMesBloEx++;
        if (_a4[isyst][isat])
          recProduct._nbMesBloNw++;
        if (_n[isyst][isat])
          recProduct._nbMesBloN1++;
      }
    }
  }

  for (int isyst=0; isyst<SystMax; isyst++)
  {
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      A_MEASUREMENT_PRODUCT &measProductSat = measurementProduct[isyst][isat];
      for (int f=F1; f<FMax; f++)
      {
        measProductSat._elimCode[f] = elimCode[f][isyst][isat];
	measProductSat._elimPhase[f] = elimPhase[f][isyst][isat];
	measProductSat._resCode[f] = resCode[f][isyst][isat];
	measProductSat._resPhase[f] = resPhase[f][isyst][isat];
      }
      measProductSat._ne = _ex[isyst][isat];
      measProductSat._nw = _a4[isyst][isat];
      measProductSat._n1 = _n[isyst][isat];

      measProductSat._ambEx = 0.0;
      measProductSat._ambA4 = 0.0;
      measProductSat._ambP = 0.0;
      measProductSat._iono = 0.0;
      measProductSat._covAmbEx = 0.0;
      measProductSat._covAmbA4 = 0.0;
      measProductSat._covAmbP = 0.0;
      measProductSat._covIono = 0.0;
      int ipass=_pass.numPass(isyst, isat);
      if (ipass != -1)
      {
        measProductSat._ambEx = _solInt[_params.getIndex(Ex,isyst,ipass)];
        measProductSat._ambA4 = _solInt[_params.getIndex(A4,isyst,ipass)];
        measProductSat._ambP = _solInt[_params.getIndex(P,isyst,ipass)];
        measProductSat._covAmbEx =sqrt(x);
        _filterInt.recoverDiagonalCovariance(1+_params.getIndex(A4,isyst,ipass), x);
        measProductSat._covAmbA4 =sqrt(x);
        _filterInt.recoverDiagonalCovariance(1+_params.getIndex(P,isyst,ipass), x);
        measProductSat._covAmbP =sqrt(x);
        measProductSat._iono = _solInt[_params.getIndex(E,isyst,ipass)];
        _filterInt.recoverDiagonalCovariance(1+_params.getIndex(E,isyst,ipass), x);
        measProductSat._covIono =sqrt(x);
      }
    }
  }

return 0;
}

//////////////////////////////////////////////////////////////////////
int A_FILTER::propagateState(const std::vector<std::vector<A_MEASUREMENT> >& measurements, const std::vector<std::vector<A_BIAS> >& biases)

{
  int iInc;
  std::vector<double> Q(_params.getNbInc(),0.0);
  
  // Code clock model noise
  for (iInc=_params.getIndex(P1,GPS); iInc<_params.getIndex(P1,GPS)+_params.getSize(P1,GPS); iInc++)
    Q[iInc] = SIG_MOD_CLK;
  for (iInc=_params.getIndex(P2,GPS); iInc<_params.getIndex(P2,GPS)+_params.getSize(P2,GPS); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(C5,GPS); iInc<_params.getIndex(C5,GPS)+_params.getSize(C5,GPS); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(P1,Glo); iInc<_params.getIndex(P1,Glo)+_params.getSize(P1,Glo); iInc++)
    Q[iInc] = SIG_MOD_CLK;
  for (iInc=_params.getIndex(P2,Glo); iInc<_params.getIndex(P2,Glo)+_params.getSize(P2,Glo); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(P1,Gal); iInc<_params.getIndex(P1,Gal)+_params.getSize(P1,Gal); iInc++)
    Q[iInc] = SIG_MOD_CLK;
  for (iInc=_params.getIndex(C5,Gal); iInc<_params.getIndex(C5,Gal)+_params.getSize(C5,Gal); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(C7,Gal); iInc<_params.getIndex(C7,Gal)+_params.getSize(C7,Gal); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(P1,Bds); iInc<_params.getIndex(P1,Bds)+_params.getSize(P1,Bds); iInc++)
    Q[iInc] = SIG_MOD_CLK;
  for (iInc=_params.getIndex(C6,Bds); iInc<_params.getIndex(C6,Bds)+_params.getSize(C6,Bds); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(C7,Bds); iInc<_params.getIndex(C7,Bds)+_params.getSize(C7,Bds); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;

  // Phase clock model noise
  for (iInc=_params.getIndex(L1,GPS); iInc<_params.getIndex(L1,GPS)+_params.getSize(L1,GPS); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(L2,GPS); iInc<_params.getIndex(L2,GPS)+_params.getSize(L2,GPS); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(L5,GPS); iInc<_params.getIndex(L5,GPS)+_params.getSize(L5,GPS); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(L1,Glo); iInc<_params.getIndex(L1,Glo)+_params.getSize(L1,Glo); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(L2,Glo); iInc<_params.getIndex(L2,Glo)+_params.getSize(L2,Glo); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(L1,Gal); iInc<_params.getIndex(L1,Gal)+_params.getSize(L1,Gal); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(L5,Gal); iInc<_params.getIndex(L5,Gal)+_params.getSize(L5,Gal); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(L7,Gal); iInc<_params.getIndex(L7,Gal)+_params.getSize(L7,Gal); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(L1,Bds); iInc<_params.getIndex(L1,Bds)+_params.getSize(L1,Bds); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(L6,Bds); iInc<_params.getIndex(L6,Bds)+_params.getSize(L6,Bds); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  for (iInc=_params.getIndex(L7,Bds); iInc<_params.getIndex(L7,Bds)+_params.getSize(L7,Bds); iInc++)
    Q[iInc] = _filterSetting._sigModBiasClk;
  
  // Tropo model noise
  for (iInc=_params.getIndex(T); iInc<_params.getIndex(T)+_params.getSize(T); iInc++)
    Q[iInc] = _filterSetting._sigModTro;

  // Iono model noise
  for (iInc=_params.getIndex(E); iInc<_params.getIndex(E)+_params.getSize(E); iInc++)
    Q[iInc] = _filterSetting._sigModIono;

  // Position model noise
  Q[_params.getIndex(X)] = _filterSetting._sigModPos;
  Q[_params.getIndex(Y)] = _filterSetting._sigModPos;
  Q[_params.getIndex(Z)] = _filterSetting._sigModPos;

  /* Covariance propagation */
  for (iInc=0; iInc<_params.getNbInc(); iInc++)
    if (Q[iInc])
      _filter.addDiagonalCovariance(1+iInc, Q[iInc]);

  /* pass transition */
  bool transition=false;
  for (int isyst=0; isyst<SystMax; isyst++)
  {      
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))    
    for (int isat=0; isat<getMaxSatSyst(isyst); isat++)
    {
      const A_MEASUREMENT & measurement = measurements[isyst][isat];
      bool meas=measurement._code[F1] || measurement._code[F2] || measurement._code[F5] || measurement._code[F6] || measurement._code[F7];
      int &nbMeasPasSat = _nbMeasPas[isyst][isat];
      int &nbMesPas0Sat = _nbMesPas0[isyst][isat];
      int &discontinuityValue = _discontinuityValue[isyst][isat];

      // Pass duration
      if (nbMeasPasSat)
      {
	if (meas)
	  nbMesPas0Sat=0;
	else
	  {
	  nbMesPas0Sat++;
	  }
      }
      else
      {
	if (meas)
          nbMeasPasSat=1;
      }

      // End pass
      if ( (nbMesPas0Sat > _filterSetting._dtMax) || ((discontinuityValue != biases[isyst][isat]._discontinuity) && (biases[isyst][isat]._discontinuity > 0)) )
      {
	transition=true;
	int ipass=_pass.numPass(isyst, isat);
	nbMeasPasSat = 0;
	if(meas)
	  nbMeasPasSat = 1;
	nbMesPas0Sat = 0;
	if ((discontinuityValue != biases[isyst][isat]._discontinuity) && (biases[isyst][isat]._discontinuity > 0))
	  discontinuityValue = biases[isyst][isat]._discontinuity;
//      if ((ipass != -1) && nbMeasPasSat)
	if (ipass != -1)
	{
	  _pass.endPass(isyst, ipass);
	  initStateVectorAmb(isyst, isat, ipass, 0.0, true);
          _filter.updateDiagonalCovariance(1+_params.getIndex(E,isyst,ipass), 0.0);
          _sol[_params.getIndex(E,isyst,ipass)] = 0.0;
	}
      }

      // Begin pass
      if (nbMeasPasSat == 1)
      {
	nbMeasPasSat = 1;
	nbMesPas0Sat = 0;
	int ipass=_pass.beginPass(isyst, isat);
	if (ipass != -1)
	{
	  initStateVectorAmb(isyst, isat, ipass, SIG_INI_AMB, true);
	  _filter.updateDiagonalCovariance(1+_params.getIndex(E,isyst,ipass), _filterSetting._sigIniIono);
	  _sol[_params.getIndex(E,isyst,ipass)] = 0.0;
	}
      }

      if (nbMeasPasSat)
	nbMeasPasSat++;
    }
  }

  // Reinit (on timeout)
  bool noPass=true;
  for (int isyst=0; isyst<SystMax; isyst++)
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))    
    for (int isat=0; isat<getMaxSatSyst(isyst); isat++)
      if (_nbMeasPas[isyst][isat])
        noPass=false;

  if (noPass && transition)
    return 1;

return 0;
}

//////////////////////////////////////////////////////////////////////
A_FILTER::A_FILTER()
{  
  _n = std::vector<std::vector<double> >(SystMax);
  _a4 = std::vector<std::vector<double> >(SystMax);
  _ex = std::vector<std::vector<double> >(SystMax);
  _nbMeasPas = std::vector<std::vector<int> >(SystMax);
  _nbMesPas0 = std::vector<std::vector<int> >(SystMax);
  _discontinuityValue = std::vector<std::vector<int> >(SystMax);

  for(int isyst=0; isyst<SystMax; isyst++)
  {
    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))    
    {
      const int maxSatSyst = getMaxSatSyst(isyst);
    _n[isyst] = std::vector<double>(maxSatSyst,0.0);
    _a4[isyst] = std::vector<double>(maxSatSyst,0.0);
    _ex[isyst] = std::vector<double>(maxSatSyst,0.0);
    _nbMeasPas[isyst] = std::vector<int>(maxSatSyst,0);
    _nbMesPas0[isyst] = std::vector<int>(maxSatSyst,0);
    _discontinuityValue[isyst] = std::vector<int>(maxSatSyst, 0);
    }
  }
  
  _rover = 0;
  _sol.clear();
  _solInt.clear();
}

//////////////////////////////////////////////////////////////////////
void A_FILTER::initState(const int rover)
{
  int iInc;
  
  _rover = rover;
  
  _params.setParam();
	     
  PRINT_LOG(_filterLog._log, "Rover=%d, _filter: %d parameters\n", _rover, _params.getNbInc());
  FFLUSH_LOG(_filterLog._log);
  
  _pass.init();
  _filter.initialize (_params.getNbInc());
  for (iInc=_params.getIndex(T); iInc<_params.getIndex(T)+_params.getSize(T); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniTro);

  for (iInc=_params.getIndex(P1,GPS); iInc<_params.getIndex(P1,GPS)+_params.getSize(P1,GPS); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, SIG_INI_CLK);
  for (iInc=_params.getIndex(P2,GPS); iInc<_params.getIndex(P2,GPS)+_params.getSize(P2,GPS); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(C5,GPS); iInc<_params.getIndex(C5,GPS)+_params.getSize(C5,GPS); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(L1,GPS); iInc<_params.getIndex(L1,GPS)+_params.getSize(L1,GPS); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(L2,GPS); iInc<_params.getIndex(L2,GPS)+_params.getSize(L2,GPS); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(L5,GPS); iInc<_params.getIndex(L5,GPS)+_params.getSize(L5,GPS); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);

  for (iInc=_params.getIndex(P1,Glo); iInc<_params.getIndex(P1,Glo)+_params.getSize(P1,Glo); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, SIG_INI_CLK);
  for (iInc=_params.getIndex(P2,Glo); iInc<_params.getIndex(P2,Glo)+_params.getSize(P2,Glo); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(L1,Glo); iInc<_params.getIndex(L1,Glo)+_params.getSize(L1,Glo); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(L2,Glo); iInc<_params.getIndex(L2,Glo)+_params.getSize(L2,Glo); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  
  for (iInc=_params.getIndex(P1,Gal); iInc<_params.getIndex(P1,Gal)+_params.getSize(P1,Gal); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, SIG_INI_CLK);
  for (iInc=_params.getIndex(C5,Gal); iInc<_params.getIndex(C5,Gal)+_params.getSize(C5,Gal); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(C7,Gal); iInc<_params.getIndex(C7,Gal)+_params.getSize(C7,Gal); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(L1,Gal); iInc<_params.getIndex(L1,Gal)+_params.getSize(L1,Gal); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(L5,Gal); iInc<_params.getIndex(L5,Gal)+_params.getSize(L5,Gal); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(L7,Gal); iInc<_params.getIndex(L7,Gal)+_params.getSize(L7,Gal); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
    
  for (iInc=_params.getIndex(P1,Bds); iInc<_params.getIndex(P1,Bds)+_params.getSize(P1,Bds); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, SIG_INI_CLK);
  for (iInc=_params.getIndex(C6,Bds); iInc<_params.getIndex(C6,Bds)+_params.getSize(C6,Bds); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(C7,Bds); iInc<_params.getIndex(C7,Bds)+_params.getSize(C7,Bds); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(L1,Bds); iInc<_params.getIndex(L1,Bds)+_params.getSize(L1,Bds); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(L6,Bds); iInc<_params.getIndex(L6,Bds)+_params.getSize(L6,Bds); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);
  for (iInc=_params.getIndex(L7,Bds); iInc<_params.getIndex(L7,Bds)+_params.getSize(L7,Bds); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniBiasClk);

  for (iInc=_params.getIndex(X); iInc<_params.getIndex(X)+_params.getSize(X); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniPos);
  for (iInc=_params.getIndex(Y); iInc<_params.getIndex(Y)+_params.getSize(Y); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniPos);
  for (iInc=_params.getIndex(Z); iInc<_params.getIndex(Z)+_params.getSize(Z); iInc++)
    _filter.updateDiagonalCovariance(1+iInc, _filterSetting._sigIniPos);

  for(int isyst=0; isyst<SystMax; isyst++)
  {    
    _n[isyst].clear();
    _a4[isyst].clear();    
    _ex[isyst].clear();   
    _nbMeasPas[isyst].clear();    
    _nbMesPas0[isyst].clear();
    _discontinuityValue[isyst].clear();

    if ((isyst == GPS && _filterSetting._gps._use) || (isyst == Glo && _filterSetting._glo._use) || (isyst == Gal && _filterSetting._gal._use) || (isyst == Bds && _filterSetting._bds._use))    
    {
      const int maxSatSyst = getMaxSatSyst(isyst);
      _n[isyst].resize(maxSatSyst,0.0);
      _a4[isyst].resize(maxSatSyst,0.0);
      _ex[isyst].resize(maxSatSyst,0.0);
      _nbMeasPas[isyst].resize(maxSatSyst,0);
      _nbMesPas0[isyst].resize(maxSatSyst,0);      
      _discontinuityValue[isyst].resize(maxSatSyst, 0);      
    }
  }
  _sol.clear();
  _sol.resize(_params.getNbInc(),0.0);
  _solInt.clear();
  _solInt.resize(_params.getNbInc(),0.0);
}
