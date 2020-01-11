#ifndef NAV_ORBIT
#define NAV_ORBIT

#include <cmath>
#include <iostream>

#include "point3.h"
#include "tool.h"

struct Orbit {
    double alpha = 0, base, now;
    Orbit(double _alpha, double _len_base) : alpha(_alpha), base(_len_base), now(_len_base) {
        if (base < 50) {
            std::cerr << "error params of orbit!" << std::endl;
            exit(-1);
        }
    }

    inline Point3 position(double beta) {
        double x = now * cos(radians(alpha)) * cos(radians(beta));
        double y = now * cos(radians(alpha)) * sin(radians(beta));
        double z = now * sin(radians(alpha));
        return Point3(x, y, z);
    }

    inline bool update(double forward) {
        double diff = now - forward;
        if (diff < 50) {
            return false;
        }
        now -= forward;
        return true;
    }
};

#endif