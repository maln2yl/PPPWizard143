/************************************************************
Nom ......... : rtrover_rover.cpp
Role ........ : rover definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#include <string>

#include "rtrover_rover.h"
#include "rtrover_vector3D.h"
#include "rtrover_utility.h"
#include "rtrover_broadcast.h"
#include "rtrover_interface.h"

//////////////////////////////////////////////////////////////////////
A_ROVER::A_ROVER(std::string roverName,double aprPos[3])
{
  _roverName=roverName;
  _roverAPrioriPosition=A_VECTOR3D(aprPos[0],aprPos[1],aprPos[2]);
  _lastTime._mjd = 0;
  _lastTime._sec = 0.0;
  _nbEpoch = 0;
  _tabSatObs = std::vector<std::vector<std::vector<rtrover_satObs> > > (MAX_EPOCH,std::vector<std::vector<rtrover_satObs> >(SystMax));
  _tabIonoCorr = std::vector< std::vector<rtrover_ionoCorr> > (SystMax);
  _stationMeasurement = std::vector<std::vector<A_MEASUREMENT_PRODUCT> > (SystMax);
  
  for (int isyst=0; isyst < SystMax ; isyst++)
  {     
    _tabIonoCorr[isyst] = std::vector<rtrover_ionoCorr> (getMaxSatSyst(isyst));
    _stationMeasurement[isyst] = std::vector<A_MEASUREMENT_PRODUCT> (getMaxSatSyst(isyst));
    
    for (int isat=0; isat < getMaxSatSyst(isyst) ; isat++)
    {
      initializeIonoCorr(_tabIonoCorr[isyst][isat]);
    }
  }
  for (int iepoch=0; iepoch<MAX_EPOCH;iepoch++)
  {
    for (int isyst=0; isyst < SystMax ; isyst++)
    {     
      _tabSatObs[iepoch][isyst] = std::vector<rtrover_satObs>(getMaxSatSyst(isyst));
      for (int isat=0; isat < getMaxSatSyst(isyst) ; isat++)
      {
	  initializeSatObs(_tabSatObs[iepoch][isyst][isat]);
      }
    }
  }
  
  _tropoCorr._value = 0.0;
}

//////////////////////////////////////////////////////////////////////
A_ROVER::~A_ROVER()
{
  for (int isyst=0; isyst < SystMax ; isyst++)
  {
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      for (int j=0;j<_nbEpoch;j++)
        delete [] _tabSatObs[j][isyst][isat]._obs;
    }
  }
}

//////////////////////////////////////////////////////////////////////
void A_ROVER::initializeIonoCorr(rtrover_ionoCorr& ionoCorr)
{
  ionoCorr._staName=NULL;
  ionoCorr._satellite._system='G';
  ionoCorr._satellite._number=0;
  ionoCorr._value=0.0;
  ionoCorr._flag=2;
}

//////////////////////////////////////////////////////////////////////
void A_ROVER::initializeSatObs(rtrover_satObs& satObs)
{
  satObs._satellite._system='G';
  satObs._satellite._number=0;
  satObs._time._mjd=0;
  satObs._time._sec=0.0;
  satObs._slotNumber=0;
  satObs._numObs=0;
  satObs._obs=NULL;
}

//////////////////////////////////////////////////////////////////////
void A_ROVER::newTabSatObs(int indice, int numSatRover,const rtrover_satObs* satObsRover,
                           std::vector<A_BROADCAST*> tabEphLastGlo)
{
  for (int i=0;i<numSatRover;i++)
  {
    int isat = -1;
    std::vector<std::vector<rtrover_satObs> > &tabSatObsEpoch = _tabSatObs[0];
    if (((satObsRover[i]._satellite._system == 'G') || (satObsRover[i]._satellite._system == 'E') || (satObsRover[i]._satellite._system == 'C')) && (satObsRover[i]._satellite._number >= 1))
    {
      isat=satObsRover[i]._satellite._number-1;
      rtrover_satObs &satObs = tabSatObsEpoch[convertSystem(satObsRover[i]._satellite._system)][isat];
      satObs = satObsRover[i];
      satObs._obs = new rtrover_obs[satObsRover[i]._numObs];
      for(int j=0; j<satObsRover[i]._numObs; j++)
        satObs._obs[j] = satObsRover[i]._obs[j];
    }
    if ((satObsRover[i]._satellite._system == 'R') && (satObsRover[i]._satellite._number >= 1))
    {
      isat=satObsRover[i]._satellite._number-1;
      rtrover_satObs &satObs = tabSatObsEpoch[Glo][isat];
      if (tabEphLastGlo[isat] != NULL)
      {
        satObs  = satObsRover[i];
	satObs._obs = new rtrover_obs[satObsRover[i]._numObs];
        for(int j=0; j<satObsRover[i]._numObs; j++)
          satObs._obs[j] = satObsRover[i]._obs[j];
	A_BROADCAST_Glo* ephLast = dynamic_cast<A_BROADCAST_Glo*>(tabEphLastGlo[isat]);
	satObs._slotNumber = (int)ephLast->getEph()._frequency_number;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////// 
void A_ROVER::deleteTabSatObs(int indice)
{
  for (int isyst=0;isyst<SystMax;isyst++)
  { 
    for (int isat=0;isat<getMaxSatSyst(isyst);isat++)
    {
      delete [] _tabSatObs[indice][isyst][isat]._obs;
    }
  }
}
