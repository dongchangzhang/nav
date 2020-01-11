#include <iostream>
#include "landmark.h"

std::ostream& operator<<(std::ostream &os, Landmark &landmark) {
    os << "(" << landmark.xyz.x << ", " << landmark.xyz.y << ", " << landmark.xyz.z << ")";
    return os;
}