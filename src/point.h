#ifndef NAV_POINT
#define NAV_POINT

#include <cmath>

struct Point {
    double x = 0, y = 0;
    Point() {}
    Point(double _x, double _y) : x(_x), y(_y) {}
    inline double len() const { return sqrt(x * x + y * y); }
};

#endif