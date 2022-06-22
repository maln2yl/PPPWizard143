/************************************************************
Nom ......... : rtrover_rover.h
Role ........ : rover definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef _RTROVER_ROVER
#define _RTROVER_ROVER

#include <string>

#include "rtrover_ppp.h"
#include "rtrover_vector3D.h"
#include "rtrover_interface.h"
#include "rtrover.h"
#include "rtrover_broadcast.h"

class A_ROVER
{
  public :
    A_ROVER(std::string roverName,double aprPos[3]);
    ~A_ROVER();
    void initializeSatObs(rtrover_satObs& satObs);
    void initializeIonoCorr(rtrover_ionoCorr& ionoCorr);
    void newTabSatObs(int indice, int numSatRover,const rtrover_satObs* satObsRover,
                      std::vector<A_BROADCAST*> tabEphLastGlo);
    void deleteTabSatObs(int indice);
    
    A_PPP _ppp;
    A_RECEIVER_PRODUCT _stationProduct;
    std::vector<std::vector<A_MEASUREMENT_PRODUCT> > _stationMeasurement;
    A_VECTOR3D _roverAPrioriPosition;
    std::vector<std::vector<std::vector<rtrover_satObs> > > _tabSatObs;
    std::vector< std::vector<rtrover_ionoCorr> > _tabIonoCorr;
    rtrover_tropoCorr _tropoCorr;
    rtrover_time _lastTime;
    int _nbEpoch;
    std::string _roverName;
};

#endif
