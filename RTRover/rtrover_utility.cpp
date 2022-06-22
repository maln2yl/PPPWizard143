/************************************************************
Nom ......... : rtrover_utility.cpp
Role ........ : common fonctions
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#include "rtrover_interface.h"
#include "rtrover_utility.h"

//////////////////////////////////////////////////////////////////////
// Dates
double diffDate(const rtrover_time *t1, const rtrover_time *t2)
{
  return ((double)t1->_mjd-(double)t2->_mjd)*86400.0+(t1->_sec-t2->_sec);
}

//////////////////////////////////////////////////////////////////////
// Convert
int convertSystem(const char system)
{
  int syst=GPS;
  switch (system)
  {
    case 'G' :
      syst = GPS;
      break;
    case 'R' :
      syst = Glo;
      break;
    case 'E' :
      syst = Gal;
      break;
    case 'C' :
      syst = Bds;
      break;
   default :
      syst = -1;
  }
  return syst;
}

