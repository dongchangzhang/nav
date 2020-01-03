// vhdsih

#include <cstdio>
#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>
#include <random>
#include <string>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>

#include "nav.h"
#include "tool.h"
#include "orbit.h"
#include "landmark.h"
#include "constants.h"

using namespace std;
using namespace Eigen;

bool run(int n, double _x, double _y, double _z, std::default_random_engine &engine) {
    double th = radians(2), wo = radians(2), ka = radians(2);
    double f = 5 * _MM, x0 = 6 * _UM, y0 = 6 * _UM;

    Matrix3d R;
    MatrixXd A, B, L, diff, error, loss;
    MatrixXd W, x, K, v;
    MatrixXd vars = MatrixXd::Zero(N_VARS, 1); // vars{Xs, Ys, Zs, th, wo, ka, f, x0, y0};
    MatrixXd vars_real, vars_base;

    vector<MatrixXd> errors;
    vector<MatrixXd> vars_solved_history;
    vector<Landmark> data, data_real, data_noise;

    vars << _x * _M, _y * _M, _z * _M, th, wo, ka, f, x0, y0;

    vars_real = vars;

    double len = sqrt(_x * _x + _y * _y + _z * _z);
    get_ideal_data(n, vars, data, len, engine);
    data_real = data;

    // add noise for vars and data
    add_noise(vars, data, engine);
    vars_base = vars;
    data_noise = data;

    vars_solved_history.push_back(vars);

    bool is_nan = false;
    for (int i = 0; i < MAX_ITERS; ++i) {
        gen_R_matrix(vars, R);
        get_matrix_A_and_B(n, data, vars, R, A, B);
        get_vector_L(n, data, vars, R, L);

        W = A * A.transpose();
        x = (B.transpose() * W * B).inverse() * B.transpose() * W * L;
        K = (A * A.transpose()).inverse() * (B * x - L);
        v = A.transpose() * K;

        // update params
        vars += x;

        // update data
        for (int k = 0; k < n; ++k) {
            data[k].X += v(k * 5 + 2, 0);
            data[k].Y += v(k * 5 + 3, 0);
            data[k].Z += v(k * 5 + 4, 0);
        }

        // calculate error
        loss = x.transpose() * x;

        if (loss(0, 0) > 1e3 ||isnan(loss(0, 0))) {
            return false;
        }
        diff = vars_real - vars;
        error = diff.transpose() * diff;

        errors.emplace_back(diff);
        vars_solved_history.push_back(vars);

        if (i % N_PRINT == 0)
            printf("step: %5d   loss: %10e   real_error: %10e\n", i, sqrt(loss(0, 0)), sqrt(error(0, 0)));

        // break rule
        if (sqrt(loss(0, 0)) < LIMIT)
            break;
    }

    dump(data_real, vars_real, vars_base, vars, errors, vars_solved_history);

    // // show result
    cout << "solved:" << endl;
    print(vars);
    // print(history_vars_base);
    // print(history_vars_real);
    return true;

}
int main() {
    int n = 10, n_nan = 0;
    Orbit orbit(60, 1010);
    Landmark pos;
    std::default_random_engine engine;
    engine.seed((unsigned) time(0));

    pos = orbit.position(10);
    while (orbit.update(10)) {
        for (int i = 0; i < 100; ++i) {
            pos = orbit.position(10);
            auto yes = run(10, pos.X, pos.Y, pos.Z, engine);
            if (!yes) ++n_nan;
        }
    }

    cout << "bad: " << n_nan << endl;
    return 0;
}