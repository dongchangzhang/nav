#include <iostream>
#include <random>
#include <ctime>
using namespace std;

int main() {
    default_random_engine e;
    e.seed((unsigned)time(0));
    normal_distribution<double> u(0, 0.034);
    double x = 0, y;
    for (int i = 0; i < 10000; ++i) {
        while (true) {
            y = fabs(u(e));
            if (y >= 0.05 and y <= 0.1 )
                break;
        }
        x += y;
    }
    cout << x / 10000 << endl;
    return 0;
}
