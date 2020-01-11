#ifndef NAV_POINT3
#define NAV_POINT3

#include <cmath>

struct Point3 {
    double x = 0, y = 0, z = 0;
    Point3() {}
    Point3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    inline double len() const {
        return sqrt(x * x + y * y + z * z);
    }
    inline double dist(const Point3 &p, bool mask_x=false, bool mask_y=false, bool mask_z=false) const {
        return sqrt(int(!mask_x) * (x - p.x) * (x - p.x) 
            + int(!mask_y) * (y - p.y) * (y - p.y) + int(!mask_z) * (z - p.z) * (z - p.z));
    }
};

#endif