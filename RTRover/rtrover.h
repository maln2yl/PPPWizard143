/************************************************************
Nom ......... : rtrover.h
Role ........ : main interface
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.3 2/15/2016
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef RTROVER_INTERFACE_EXT_H
#define RTROVER_INTERFACE_EXT_H

#include <string>

#include "rtrover_interface.h"

enum correctionType { RTIGS, SBAS };

/// Input additional options
struct rtrover_additionalOpt {
  double _dt;                    ///< mesurement interval
  int _minFix;                   ///< minimum step before AR
  int _maxElim;                  ///< RAIM maximum rejection
  int _raim;                     ///< advanced RAIM
  double _maxAge;                ///< maximum RTCM correction age
  double _thrMap;                ///< tropo mapping threshold (1/sin(ele))
  double _sigIniTro;             ///< tropo initial noise
  double _sigModTro;             ///< tropo model noise
  int _nbSatFixAmb;              ///< minimum satellite for ambiguity resolution
  double _thrAmb;                ///< ambiguity threshold for AR
  double _sigIniClkBias;         ///< initial clock bias noise
  double _sigModClkBias;         ///< mode clock bias noise
  double _sigIniIono;            ///< initial iono noise
  double _sigModIono;            ///< model iono noise
  double _sigIniPos;             ///< initial position noise
  double _sigModPos;             ///< model position noise
  double _sigMeasIono[3];        ///< iono measurement noise [precise ; SBAS ; global ]
  double _thrMeasIono;           ///< iono measurement rejection threshold
  double _sigMeasTropo;          ///< Tropo measurement noise
  double _thrMeasTropo;          ///< Tropo measurement rejection threshold
  double _sigMeasCodeGps;        ///< code GPS measurement noise
  double _sigMeasPhaseGps;       ///< phase GPS measurement noise
  double _sigMeasCodeGlo;        ///< code Glonass measurement noise
  double _sigMeasPhaseGlo;       ///< phase Glonass measurement noise
  double _sigMeasCodeGal;        ///< code Galileo measurement noise
  double _sigMeasPhaseGal;       ///< phase Galileo measurement noise
  double _sigMeasCodeBds[3];     ///< code Beidou measurement noise [ GEO; IGSO; MEO ]
  double _sigMeasPhaseBds[3];    ///< phase Beidou measurement noise [ GEO; IGSO; MEO ]
  double _thrMeasCode;           ///< code measurement rejection threshold
  double _thrMeasPhase;          ///< phase measurement rejection threshold
  double _smooth;                ///< doppler smoothing coefficient
  int _dtMax;         	         ///< maximum measurement gap
  int _correction;               ///< correction type (RTIGS or SBAS)
  bool _useGPS;                  ///< use GPS observations
  bool _useGalileo;              ///< use Galileo observations (in addition to GPS)
  bool _useBeidou;               ///< use BeiDou observations (in addition to GPS)
  int _fixN1, _fixNw, _fixEx;    ///< use tri-frequency ambiguity resolution
  int _reset;                    ///< Time between consecutive reset (for convergence tests)
  int _outputVerbose;            ///< verbose output
  int _lowlevel;                 ///< lowlevel
};

/// Galileo ephemeris (be careful to use only INAV)
struct rtrover_ephGal {
  rtrover_satellite _satellite;       ///< satellite 
  rtrover_time      _TOC;             ///< Reference Time of Clocks
  rtrover_time      _TOE;             ///< Reference Time of Ephemeris
  unsigned short    _IODNav;            ///< Issue of Data (Nav)
  double            _clock_bias;      ///< Clock correction [s]    
  double            _clock_drift;     ///< Clock correction rate [s/s]  
  double            _clock_driftrate; ///< Clock correction rate rate [s/s^2]
  double            _Crs;             ///< sine correction to radius [m]    
  double            _Delta_n;         ///< mean motion correction [rad/s]
  double            _M0;              ///< mean anomaly [rad]
  double            _Cuc;             ///< cosine correction to argument of latitude [rad]  
  double            _e;               ///< numerical eccentricity       
  double            _Cus;             ///< sine correction to argument of latitude [rad]
  double            _sqrt_A;          ///< semimajor axis square root [m^0.5]
  double            _Cic;             ///< cosine correction to inclination [rad]  
  double            _OMEGA0;          ///< longitude of the ascending node [rad]
  double            _Cis;             ///< sine correction to inclination [rad]  
  double            _i0;              ///< inclination angle [rad]  
  double            _Crc;             ///< cosine sine correction to radius [m]
  double            _omega;           ///< argument of perigee [rad]  
  double            _OMEGADOT;        ///< rate of right ascension [rad/s]
  double            _IDOT;            ///< rate of inclination angle [rad/s]
  double            _BGD_1_5A;        ///< differential group delay E1-E5A [s]
  double            _BGD_1_5B;        ///< differential group delay E1-E5B [s]
  int               _E5aHS;           ///< health flag
  int               _E5bHS;           ///< health flag
};

/// BeiDou ephemeris
struct rtrover_ephBds {
  rtrover_satellite _satellite;       ///< satellite 
  rtrover_time      _TOC;             ///< Reference Time of Clocks
  rtrover_time      _TOE;             ///< Reference Time of Ephemeris
  unsigned short    _IODE;            ///< Issue of Data (Nav)
  unsigned short    _IODC;            ///< Issue of Data (Clocks)
  double            _clock_bias;      ///< Clock correction [s]    
  double            _clock_drift;     ///< Clock correction rate [s/s]  
  double            _clock_driftrate; ///< Clock correction rate rate [s/s^2]
  double            _Crs;             ///< sine correction to radius [m]    
  double            _Delta_n;         ///< mean motion correction [rad/s]
  double            _M0;              ///< mean anomaly [rad]
  double            _Cuc;             ///< cosine correction to argument of latitude [rad]  
  double            _e;               ///< numerical eccentricity       
  double            _Cus;             ///< sine correction to argument of latitude [rad]
  double            _sqrt_A;          ///< semimajor axis square root [m^0.5]
  double            _Cic;             ///< cosine correction to inclination [rad]  
  double            _OMEGA0;          ///< longitude of the ascending node [rad]
  double            _Cis;             ///< - sine correction to inclination [rad]  
  double            _i0;              ///< inclination angle [rad]  
  double            _Crc;             ///< cosine sine correction to radius [m]
  double            _omega;           ///< argument of perigee [rad]  
  double            _OMEGADOT;        ///< rate of right ascension [rad/s]
  double            _IDOT;            ///< rate of inclination angle [rad/s]
  double            _TGD_1_3;         ///< differential group delay B1-B3 [s]
  double            _TGD_2_3;         ///< differential group delay B2-B3 [s]
  int               _health;          ///< health flag
  double            _toes;            ///< Toe in week [s]
};

/// Satellite-Specific Biases (DCB)
struct rtrover_satDCBBiases {
  rtrover_satellite _satellite;  ///< satellite 
  rtrover_time _time;            ///< bias reference time
  double _P1P2;                  ///< value of bias
  double _P1C1;                  ///< value of bias
  double _P2C2;                  ///< value of bias
};

/// Sbas Correction
struct rtrover_sbasCorr {
  rtrover_satellite _satellite;  ///< satellite 
  rtrover_time _time;            ///< bias reference time
  unsigned short _IODE;          ///< Issue of Data (Ephemeris)
  double _orb[3];                ///< value of the correction (position)
  double _clk;                   ///< clock correction
};

/// Tropospheric Correction
struct rtrover_tropoCorr {
  double _value;                 ///< value of the correction
};

RTROVER_API void rtrover_setAdditionalOpt(const rtrover_additionalOpt* opt);

/// Add (and internally store) Galileo Ephemeris
RTROVER_API void rtrover_putGalEphemeris(const rtrover_ephGal* eph);

/// Add (and internally store) BeiDou Ephemeris
RTROVER_API void rtrover_putBdsEphemeris(const rtrover_ephBds* eph);

/// Add (and internally store) DCB biases
RTROVER_API void rtrover_putDCBBiases(int numBiases, const rtrover_satDCBBiases* biases);

/// Add (and internally store) sbas corrections
RTROVER_API void rtrover_putSbasCorrs(int numCorr, const rtrover_sbasCorr* corr);

/// Add (and internally store) tropospheric corrections
RTROVER_API void rtrover_putTropoCorrs(const rtrover_tropoCorr* corr);

/// Get rover information
RTROVER_API void rtrover_getRoverInformation(int *day, double *sec, double pos[6], double tropo[3]);
/// Get satellite information
RTROVER_API void rtrover_getSatelliteInformation(rtrover_satellite satellite, double pos[3], double amb[6], double iono[2]);

/// Set Rover
RTROVER_API void rtrover_addRover(std::string roverName, double aprPos[3]);
RTROVER_API void rtrover_setCurrentRover(const int n);
#endif
