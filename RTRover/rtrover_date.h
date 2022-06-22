/************************************************************
Nom ......... : rtrover_date.h
Role ........ : date definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.0 9/30/2014
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef _RTROVER_DATE
#define _RTROVER_DATE

class A_DATE
{ 
  public:
    A_DATE() : _day(0), _second(0.0) {}
    A_DATE(int day, double second) : _day(day), _second(second) {}
    
    const int getDay() const { return _day; }
    const double getSecond() const { return _second; }
    void setDay(const int day) { _day = day ; }
    void setSecond(const double second) { _second = second; }
    
    void calendar(int *day, int *month, int *year, int *hour, int *minute, double *second);
    void fromCalendar(const int day, const int month, const int year, const int hour, const int minute, const double second);
    A_DATE operator+(const double d);
    double operator- (const A_DATE &date) const;
    
  private :
    int _day;
    double _second;
};

#endif
