#ifndef RTROVER_INTERFACE_H
#define RTROVER_INTERFACE_H

/// \mainpage
///
/// RTRover (Real-Time Rover) library implements a GNSS client application
/// which can be used as
///   - SPP (Single-Point Positioning) client
///   - PPP (Precise-Point Positioning) client
///   - RTK (Real-Time Kinematic) client
///
/// It can process data stemming from various global positioning systems
/// (GPS NAVSTAR, Glonass, etc.) 
///
/// It is capable of processing
///   - dual frequency data
///   - single frequency data

/// \file rtrover_interface.h
/// 
/// \brief Interface for the RTRover library
///
/// This file contains the declaration of all structures and functions
/// necessary to run the application

/// Extern Definiton for DLL Support
#ifdef WIN32
#  include <windows.h>
#  ifdef RTRover_EXPORTS // automaticall added by cmake
#    define RTROVER_API EXTERN __declspec(dllexport)
#  else
#    ifdef RTRover_STATIC_LIB // added in CMakeLists.txt
#      define RTROVER_API EXTERN
#    else
#      define RTROVER_API EXTERN __declspec(dllimport)
#    endif
#  endif
#else
#  define RTROVER_API EXTERN
#endif

/// Extern Definiton for C++
#ifdef __cplusplus
#  define EXTERN extern "C"
#else
#  define EXTERN extern
#endif

/// Application Mode
enum rtrover_mode {
  mode_PPP_DF,  ///< dual-frequency precise point positioning
  mode_SPP_DF,  ///< dual-frequency single point positioning
  mode_PPP_SF,  ///< single-frequency precise precise point positioning
  mode_SPP_SF,  ///< single-frequency single point positioning
  mode_PPP_AR,  ///< precise point positining with ambiguity resolution
  mode_RTK,     ///< real-time kinematics
  mode_PPP_FTTF ///< fast time-to-fix ambiguity resolution PPP
};

/// Input options
struct rtrover_opt {
  rtrover_mode _mode;         ///< application mode
  const char* _roverName;     ///< name of the rover station
  const char* _baseName;      ///< name of the base station (if relevant)
  double _xyzAprRover[3];     ///< a priori coordinates of the rover station
  double _xyzAprBase[3];      ///< a priori coordinates of the base station
  double _neuEccRover[3];     ///< rover antenna eccentricity components
  double _neuEccBase[3];      ///< base antenna eccentricity components
  const char* _antNameRover;  ///< name of the rover antenna
  const char* _antNameBase;   ///< name of the base antenna
  const char* _antexFileName; ///< name of the ANTEX file
  int  _logLevel;             ///< level of details in log (0 ... no log, 1 ... normal log, 2 ... detailed log)
  int  _minobs;               ///< minimal number of observations
  bool _useGlonass;           ///< use Glonass observations (in addition to GPS)
};

/// Time
struct rtrover_time {
  unsigned int _mjd; ///< modified Julian date
  double       _sec; ///< seconds of day
};

/// Unique satellite identification
struct rtrover_satellite {
  char _system; ///< satellite system ('G' = GPS, 'R' = Glonass, 'E' = Galileo)
  int  _number; ///< satellite number (PRN if GPS, slot number if Glonass)
};

/// Output
struct rtrover_output {
  rtrover_time _epoTime;      ///< time of the processed epoch
  double       _xyzRover[3];  ///< resulting rover coordinates
  double       _covMatrix[6]; ///< covariance matrix of rover coordinates (upper triangle) 
  char*        _log;          ///< log message
  bool         _error;        ///< error flag
};

/// GPS ephemeris
struct rtrover_ephGPS {
  rtrover_satellite _satellite;       ///< satellite 
  rtrover_time      _TOC;             ///< Reference Time of Clocks
  rtrover_time      _TOE;             ///< Reference Time of Ephemeris
  unsigned short    _IODE;            ///< Issue of Data (Ephemeris)
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
  double            _Cis;             ///< sine correction to inclination [rad]  
  double            _i0;              ///< inclination angle [rad]  
  double            _Crc;             ///< cosine sine correction to radius [m]
  double            _omega;           ///< argument of perigee [rad]  
  double            _OMEGADOT;        ///< rate of right ascension [rad/s]
  double            _IDOT;            ///< rate of inclination angle [rad/s]
  double            _TGD;             ///< differential group delay L1-L2 [s]
  int               _health;          ///< health flag
};

/// Glonass ephemeris
struct rtrover_ephGlo {
  rtrover_satellite _satellite;        ///< satellite 
  rtrover_time      _timeUTC;          ///< Reference Time (UTC)
  int               _gps_utc;          ///< GPS - UTC [s]      
  double            _E;                ///< ephemeris age [days]   
  double            _tau;              ///< clock correction [s]      
  double            _gamma;            ///< clock correction rate [s/s]
  double            _x_pos;            ///< position x-component [km]     
  double            _x_velocity;       ///< velocity x-component [km/s]   
  double            _x_acceleration;   ///< acceleration x-component [km/s^2] 
  double            _y_pos;            ///< position y-component [km]         
  double            _y_velocity;       ///< velocity y-component [km/s]       
  double            _y_acceleration;   ///< acceleration y-component [km/s^2] 
  double            _z_pos;            ///< position z-component [km]         
  double            _z_velocity;       ///< velocity z-component [km/s]       
  double            _z_acceleration;   ///< acceleration z-component [km/s^2] 
  double            _health;           ///< health flag (0 = O.K.)
  double            _frequency_number; ///< frequency number (channel)
};

/// Single GNSS observation
struct rtrover_obs  {
  char   _rnxType2ch[2]; ///< RINEX version 3 observation type description (band and attribute)
  double _code;          ///< code (pseudorange) observation value
  bool   _codeValid;     ///< code validity flag
  double _phase;         ///< phase observation value
  bool   _phaseValid;    ///< phase validity flag
  double _doppler;       ///< doppler observation value
  bool   _dopplerValid;  ///< doppler validity flag
  double _snr;           ///< signal-to-noise value
  bool   _snrValid;      ///< signal-to-noise validity flag
  bool   _slip;          ///< cycle-slip flag
  int    _slipCounter;   ///< cycle-slip counter (negative value = undefined);
};

/// GNSS observations of a single satellite
struct rtrover_satObs {
  rtrover_satellite _satellite;  ///< satellite 
  rtrover_time      _time;       ///< observation time (according to receiver clock)
  int               _slotNumber; ///< slot number for Glonass satellites
  int               _numObs;     ///< number of observations
  rtrover_obs*      _obs;        ///< array of observations
};

/// Satellite Ephemeris Correction 
struct rtrover_orbCorr {
  rtrover_satellite _satellite;    ///< satellite 
  unsigned int    _iod;            ///< issue of data
  rtrover_time      _time;         ///< correction reference time
  double            _rao[3];       ///< radial, along-track, and out-of-plane correction components
  double            _dotRao[3];    ///< radial, along-track, and out-of-plane correction rate components
  double            _dotDotRao[3]; ///< radial, along-track, and out-of-plane correction rate rate components
};

/// Satellite Clock Correction 
struct rtrover_clkCorr {
  rtrover_satellite _satellite;  ///< satellite 
  unsigned int    _iod;        ///< issue of data
  rtrover_time      _time;       ///< correction reference time
  double            _dClk;       ///< clock correction 
  double            _dotDClk;    ///< clock correction rate
  double            _dotDotDClk; ///< clock correction rate rate
};

/// Single Bias
struct rtrover_bias {
  char   _rnxType3ch[3]; ///< bias description (RINEX v3 code or special)
  double _value;         ///< bias value
};

/// Satellite-Specific Biases
struct rtrover_satBiases {
  rtrover_satellite _satellite;  ///< satellite 
  rtrover_time      _time;       ///< bias reference time
  int               _numBiases;  ///< number of biases
  rtrover_bias*     _biases;     ///< array of biases
};

/// Single Phase Bias
struct rtrover_phaseBias {
  char   _rnxType3ch[3];         ///< bias description (RINEX v3 code or special)
  double _value;                 ///< bias value (m)
  int _discontinuityValue;       ///< Discontinuity counter
  bool _intIndicator;            ///< Integer indicator
  int _wlIntIndicator;           ///< WL-integer indicator
};

/// Satellite-Specific Phase Biases
struct rtrover_satPhaseBiases {
  rtrover_satellite _satellite;  ///< satellite 
  rtrover_time _time;            ///< bias reference time
  double _yawAngleValue;         ///< yaw (circles)
  bool _dispersive;              ///< Dispersive Bias Consistency Indicator
  bool _mwConsistency;           ///< MW Consistency Indicator
  int _numBiases;                ///< number of biases
  rtrover_phaseBias* _biases;    ///< array of biases
};

/// Ionospheric Correction
struct rtrover_ionoCorr {
  const char*       _staName;   ///< name (ID) of the reference station
  rtrover_satellite _satellite; ///< satellite 
  double            _value;     ///< value of the correction
  int               _flag;      ///< quality flag
};

/// Initialize rtrover_opt structure
RTROVER_API void rtrover_initOptions(rtrover_opt* opt);

/// Initialize rtrover_obs structure
RTROVER_API void rtrover_initObs(rtrover_obs* obs);

/// \brief Set (and internally store) all input options
///
/// This function must be called at the begining (before any other function),
/// it may be, however, called later changing the options on-the-fly.
RTROVER_API void rtrover_setOptions(const rtrover_opt* opt);

/// Add (and internally store) GPS Ephemeris
RTROVER_API void rtrover_putGPSEphemeris(const rtrover_ephGPS* eph);

/// Add (and internally store) Glonass Ephemeris
RTROVER_API void rtrover_putGloEphemeris(const rtrover_ephGlo* eph);

/// Add (and internally store) Ephemeris Corrections
RTROVER_API void rtrover_putOrbCorrections(int numCorr, const rtrover_orbCorr* corr);

/// Add (and internally store) Clock Corrections
RTROVER_API void rtrover_putClkCorrections(int numCorr, const rtrover_clkCorr* corr);

/// Add (and internally store) Satellite-Specific Biases
RTROVER_API void rtrover_putBiases(int numBiases, const rtrover_satBiases* biases);

/// Add (and internally store) Satellite-Specific Phase Biases
RTROVER_API void rtrover_putPhaseBiases(int numBiases, const rtrover_satPhaseBiases* biases);

/// Add (and internally store) ionospheric corrections
RTROVER_API void rtrover_putIonoCorrs(int numCorr, const rtrover_ionoCorr* corr);

/// Close processing, clean-up the memory
RTROVER_API void rtrover_destroy();

/// Free memory allocated in output structure
RTROVER_API void rtrover_freeOutput(rtrover_output* output);

/// \brief Process single epoch
///
/// This function call performs the actual processing of a single epoch of
/// data. The calling program must ensure that before calling this function
/// the input options are set, satellite ephemerides are available and 
/// (optionally) the orbital and clock corrections are added.
RTROVER_API 
void rtrover_processEpoch(int numSatRover,                   ///< number of satellites (rover)
                          const rtrover_satObs* satObsRover, ///< observations (rover)
                          int numSatBase,                    ///< number of satellites (base)
                          const rtrover_satObs* satObsBase,  ///< observations (base)
                          rtrover_output* output             ///< output
                          );

#endif
