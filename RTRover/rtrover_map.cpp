/************************************************************
Nom ......... : rtrover_map.cpp
Role ........ : map definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "rtrover_vector3D.h"
#include "rtrover_map.h"
#include "rtrover_utility.h"

//////////////////////////////////////////////////////////////////////
static double computePhase(const double angle, const A_PHASE_MAP &map)
{

  int i1, i2;

  if (map._phasePattern.size() == 0.0)
    return 0.0;

  i1 = (int)((angle - map._zen1)/map._dZen);
  i2 = i1+1;

  if (i1 < 0)
    return map._phasePattern[0];
  if (i1 >= (int)map._phasePattern.size()-1)
    return map._phasePattern[map._phasePattern.size()-1];

  return ((map._phasePattern[i2]-map._phasePattern[i1])/
          ((map._zen1+i2*map._dZen)-(map._zen1+i1*map._dZen))) *
         (angle-(map._zen1+i1*map._dZen)) + map._phasePattern[i1];

}

//////////////////////////////////////////////////////////////////////
void A_SATELLITE_PHASE_MAP::phaseCenterVariation(const A_VECTOR3D &sta, const A_VECTOR3D &sat, double *x1, double *x2) const
{
  A_VECTOR3D V = sat - sta;
  double C0b  = V.getModule();
  double CosPhiSat = (sat*V)/sat.getModule()/C0b;
  double PhiSat = (acos(CosPhiSat)*180)/M_PI;
  *x1 = computePhase(PhiSat, _mapPhaseL1);
  *x2 = computePhase(PhiSat, _mapPhaseL2);
}

//////////////////////////////////////////////////////////////////////
void A_SATELLITE_PHASE_MAP::clear() 
{
  _prn=0;
  _phaseMapL1x=0;
  _phaseMapL1y=0;
  _phaseMapL1z=0;
  _phaseMapL2x=0;
  _phaseMapL2y=0;
  _phaseMapL2z=0;
  _mapPhaseL1._zen1 = 0;
  _mapPhaseL1._zen2 = 0;
  _mapPhaseL1._dZen = 0;
  _mapPhaseL1._phasePattern.clear();
  _mapPhaseL2._zen1 = 0;
  _mapPhaseL2._zen2 = 0;
  _mapPhaseL2._dZen = 0;
  _mapPhaseL2._phasePattern.clear();
}

//////////////////////////////////////////////////////////////////////
A_MAP::A_MAP()
{
  _satellites=std::vector<std::vector<A_SATELLITE_PHASE_MAP> >(SystMax);
  for (int isyst=0; isyst < SystMax ; isyst++)
    _satellites[isyst] = std::vector<A_SATELLITE_PHASE_MAP>(getMaxSatSyst(isyst));

}

//////////////////////////////////////////////////////////////////////
void A_MAP::loadATX(const char *fileName)
{
#define L_MAX 500
  FILE *fd;

  for (int isyst=0; isyst < SystMax ; isyst++)
  {
    _satellites[isyst].clear();
    _satellites[isyst].resize(getMaxSatSyst(isyst));
  }
    
  if (*fileName == 0)
    return;
  fd=fopen(fileName, "r");
  if (fd == NULL)
    return;

  char line[L_MAX+1]={"\n"};
  A_SATELLITE_PHASE_MAP satData;
  int syst = -1;
  int freq;

  freq = 0;
  while(!feof(fd))
  {
    *line = 0;
    if (fgets(line, L_MAX, fd) == NULL)
      continue;
    if (!strncmp(line+60, "START OF ANTENNA", 16))
    {
//      printf("START OF ANTENNA\n");
      satData.setPrn(0);
      freq = 0;
    }
    if (!strncmp(line+60, "START OF FREQUENCY", 18))
    {
      char a[20];
      sscanf(line, "%s", a);
      freq = atoi(a+1);
//      printf("START OF FREQUENCY %d\n", freq);
    }
    if (!strncmp(line+60, "END OF ANTENNA", 14))
    {
      if (satData.getPrn())
      {
	if (satData.getPrn()-1<getMaxSatSyst(syst))
          _satellites[syst][satData.getPrn()-1] = satData;
      }
      satData.clear();
      syst = -1;
      freq = 0;
//      printf("END OF ANTENNA\n");
    }
    if (!strncmp(line+60, "TYPE / SERIAL NO", 16))
    {
      char a1[20], a2[20], a3[20];
      sscanf(line, "%s%s%s", a1, a2, a3);
      if (!strcmp(a1, "BLOCK"))
      {
        if (!strcmp(a1, "BLOCK") && (*a3 == 'G'))
        {
          satData.setPrn(atoi(a3+1));
	  syst = GPS;
          //printf("TYPE / SERIAL NO: %s %s %s %d\n", a1, a2, a3, satData.getPrn());
        }
      }
      else if ((!strncmp(a1, "GLONASS",7) && (*a2 == 'R')) || (!strncmp(a1, "GALILEO",7) && (*a2 == 'E')) || (!strncmp(a1, "BEIDOU",6) && (*a2 == 'C')))
      {
        satData.setPrn(atoi(a2+1));
	syst = convertSystem(*a2);
        //printf("TYPE / SERIAL NO: %s %s %d\n", a1, a2, satData.getPrn());
      }
    }
    if (!strncmp(line+60, "ZEN1 / ZEN2 / DZEN", 18))
    {
      char a1[20], a2[20], a3[20];
      sscanf(line, "%s%s%s", a1, a2, a3);
      satData.getMapPhaseL1()._zen1 = atof(a1);
      satData.getMapPhaseL1()._zen2 = atof(a2);
      satData.getMapPhaseL1()._dZen = atof(a3);
      satData.getMapPhaseL1()._phasePattern.clear();
      if (satData.getMapPhaseL1()._dZen)
        satData.getMapPhaseL1()._phasePattern.resize((int)((satData.getMapPhaseL1()._zen2-satData.getMapPhaseL1()._zen1)/satData.getMapPhaseL1()._dZen)+1);
      satData.getMapPhaseL2()._zen1 = atof(a1);
      satData.getMapPhaseL2()._zen2 = atof(a2);
      satData.getMapPhaseL2()._dZen = atof(a3);
      satData.getMapPhaseL2()._phasePattern.resize(satData.getMapPhaseL1()._phasePattern.size());
//      printf("ZEN1 / ZEN2 / DZEN %s %s %s\n", a1, a2, a3);
    }
    if (!strncmp(line+60, "NORTH / EAST / UP", 17))
    {
      char a1[20], a2[20], a3[20];
      sscanf(line, "%s%s%s", a1, a2, a3);
      if (freq == 1)
      {
        satData.setPhaseMapL1x(atof(a1)/1000.0);
        satData.setPhaseMapL1y(atof(a2)/1000.0);
        satData.setPhaseMapL1z(atof(a3)/1000.0);
      }
      if ((freq == 2) || ((freq == 5) && (syst == Gal)) ||  ((freq == 7) && (syst == Bds)))
      {
        satData.setPhaseMapL2x(atof(a1)/1000.0);
        satData.setPhaseMapL2y(atof(a2)/1000.0);
        satData.setPhaseMapL2z(atof(a3)/1000.0);
      }
//      printf("NORTH / EAST / UP %s %s %s\n", a1, a2, a3);
    }
    if (!strncmp(line+3, "NOAZI", 5))
    {
      int i;
      char *pch;
      if (freq == 1)
      {
        pch = strtok (line+9," ");
        for (i=0;i<(int)satData.getMapPhaseL1()._phasePattern.size();i++)
        {
          satData.getMapPhaseL1()._phasePattern[i] = atof(pch)/1000.0;
          pch = strtok (NULL, " ");
        }
      }
      if (freq == 2)
      {
        pch = strtok (line+9," ");
        for (i=0;i<(int)satData.getMapPhaseL2()._phasePattern.size();i++)
        {
          satData.getMapPhaseL2()._phasePattern[i] = atof(pch)/1000.0;
          pch = strtok (NULL, " ");
        }
      }
//      printf("NOAZI\n");
    }

  }

  fclose(fd);
}

