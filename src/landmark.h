#ifndef NAV_LANDMARK
#define NAV_LANDMARK

#include "point.h"
#include "point3.h"

struct Landmark {
    Point xy;
    Point3 xyz;

    Landmark() {}
    Landmark(double _X, double _Y, double _Z, double _px, double _py) :
        xyz(Point3(_X, _Y, _Z)), xy(Point(_px, _py)) {}
};

std::ostream& operator<<(std::ostream &os, Landmark &landmark);

#endif