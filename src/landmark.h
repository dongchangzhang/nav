#ifndef NAV_LANDMARK
#define NAV_LANDMARK

struct Landmark {
    double X, Y, Z, px, py;
    Landmark() {}
    Landmark(double _X, double _Y, double _Z, double _px, double _py) 
        : X(_X), Y(_Y), Z(_Z), px(_px), py(_py) {}
};

std::ostream& operator<<(std::ostream &os, Landmark &landmark);

#endif