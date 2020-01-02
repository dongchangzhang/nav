#include <iostream>
#include "landmark.h"

std::ostream& operator<<(std::ostream &os, Landmark &landmark) {
    os << "(" << landmark.X << ", " << landmark.Y << ", " << landmark.Z << ")";
    return os;
}