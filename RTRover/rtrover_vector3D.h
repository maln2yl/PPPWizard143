/************************************************************
Nom ......... : rtrover_vector3D.cpp
Role ........ : 3D vector definition and handling
Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
Version ..... : V1.0 9/30/2014
Licence ..... : see file LICENSE
 
CNES authorize you, free of charge, to circulate and
distribute for no charge, for purposes other than commercial,
the source and/or object code of COMPOSITE SOFTWARE on any
present and future support.
 
************************************************************/
#ifndef A_VECTOR3D_H
#define A_VECTOR3D_H

#include <cmath>

class A_VECTOR3D
{   
    public:
        A_VECTOR3D() : _x(0), _y(0), _z(0) {}
        A_VECTOR3D(const double x, const double y, const double z) : _x(x), _y(y), _z(z) {}
        A_VECTOR3D(const A_VECTOR3D& V) { _x = V._x; _y = V._y; _z = V._z; }

        const double getX() const { return _x ; }
	const double getY() const { return _y ; }
	const double getZ() const { return _z ; }
	void setX(const double x) { _x=x; }
	void setY(const double y) { _y=y; }
	void setZ(const double z) { _z=z; }
	
        double getModule() const { return sqrt(_x*_x + _y*_y + _z*_z); }
        void doubleEtModule(const double m) { double M = getModule(); _x = _x*m/M; _y = _y*m/M; _z = _z*m/M; }
        
        A_VECTOR3D operator+ (const A_VECTOR3D& V) const {return A_VECTOR3D(_x+V._x, _y+V._y, _z+V._z);}
        A_VECTOR3D operator- (const A_VECTOR3D& V) const {return A_VECTOR3D(_x-V._x, _y-V._y, _z-V._z);}
        A_VECTOR3D operator- () const {return A_VECTOR3D(-_x, -_y, -_z);}
        A_VECTOR3D operator* (const double R) const {return A_VECTOR3D(_x*R, _y*R, _z*R);}
        A_VECTOR3D operator/ (const double R) const {return A_VECTOR3D(_x/R, _y/R, _z/R);}
        double  operator* (const A_VECTOR3D& V) const {return _x*V._x + _y*V._y + _z*V._z;}
        A_VECTOR3D& operator+=(const A_VECTOR3D& V) {_x += V._x; _y += V._y; _z += V._z; return *this;}
        A_VECTOR3D& operator-=(const A_VECTOR3D& V) {_x -= V._x; _y -= V._y; _z -= V._z; return *this;}
        A_VECTOR3D& operator*=(const double R) {_x = _x * R; _y = _y * R; _z = _z * R; return *this;}
        A_VECTOR3D& operator/=(const double R) {_x = _x / R; _y = _y / R; _z = _z / R; return *this;}
        bool operator==(const A_VECTOR3D& V) const {return (_x == V._x && _y == V._y && _z == V._z);}
        bool operator!=(const A_VECTOR3D& V) const {return (_x != V._x || _y != V._y || _z != V._z);}
        
        //cross product
        A_VECTOR3D operator^  (const A_VECTOR3D& V) const {return A_VECTOR3D(_y*V._z - _z*V._y, _z*V._x - _x*V._z, _x*V._y - _y*V._x);}
        A_VECTOR3D& operator^= (const A_VECTOR3D& V) {_x = _y*V._z - _z*V._y; _y = _z*V._x - _x*V._z; _z = _x*V._y - _y*V._x; return *this;}
	
    private:
      double _x, _y, _z;
};

#endif
