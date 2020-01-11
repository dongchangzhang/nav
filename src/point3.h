#ifndef NAV_POINT3
#define NAV_POINT3

#include <cmath>

struct Point3 {
    double x, y, z;
    Point3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    inline double len() {
        return sqrt(x * x + y * y + z * z);
    }
};

#endif