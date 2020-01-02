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

int run(int N, int _x, int _y, int _z) {
    double Xs = _x * _M, Ys = _y * _M, Zs = _z * _M, th = radians(2), wo = radians(2), ka = radians(2), f = 5 * _MM, x0 = 6 * _UM, y0 = 6 * _UM;

    // array<double, 9> vars{Xs, Ys, Zs, th, wo, ka, f, x0, y0};
    MatrixXd vars = MatrixXd::Zero(N_VARS, 1);
    vars << Xs, Ys, Zs, th, wo, ka, f, x0, y0;

    // save real params
    MatrixXd history_vars_real = vars;

    // data
    vector<Landmark> data, data_base;
    double len = sqrt(_x * _x + _y * _y + _z * _z);
    get_ideal_data(N, vars, data, len);

    // add error
    add_noise(vars, data);

    data_base = data;

    // save params-base
    MatrixXd history_vars_base = vars;

    Matrix3d R;
    MatrixXd A, B, L, error, error_matrix, loss;
    MatrixXd W, x, K, v;
    vector<MatrixXd> errors;
    vector<MatrixXd> vars_solved_history;

    vars_solved_history.push_back(vars);

    string label;
    cout << " -- start -- " << endl;

    int nan = 0;

    for (int i = 0; i < MAX_ITERS; ++i) {
        gen_R_matrix(vars, R);
        get_matrix_A_and_B(N, data, vars, R, A, B);
        get_vector_L(N, data, vars, R, L);

        W = A * A.transpose();
        x = (B.transpose() * W * B).inverse() * B.transpose() * W * L;
        K = (A * A.transpose()).inverse() * (B * x - L);
        v = A.transpose() * K;


        // update params
        vars += x;

        // update data
        for (int j = 0; j < N; ++j) {
            // data[j].px += v(j * 5 + 0, 0);
            // data[j].py += v(j * 5 + 1, 0);
            data[j].X += v(j * 5 + 2, 0);
            data[j].Y += v(j * 5 + 3, 0);
            data[j].Z += v(j * 5 + 4, 0);
        }

        // calculate error
        loss = x.transpose() * x;

        if (isnan(loss(0, 0))) {
            nan = 1;
            break;
        }
        error_matrix = history_vars_real - vars;
        error = error_matrix.transpose() * error_matrix;

        errors.emplace_back(error_matrix);
        vars_solved_history.push_back(vars);

        if (i % N_PRINT == 0)
            printf("step: %5d   loss: %10e   real_error: %10e\n", i, sqrt(loss(0, 0)), sqrt(error(0, 0)));

        // break rule
        if (sqrt(loss(0, 0)) < LIMIT)
            break;
    }

    dump(data_base, history_vars_real, history_vars_base, vars, errors, vars_solved_history);

    // // show result
    label = "solved";
    print(vars, label);
    // label = "base";
    // print(history_vars_base, label);
    // label = "real";
    // print(history_vars_real, label);
    return nan;

}

int main() {
    int N = 10, n_nan = 0;
    int x = 2, y = 2, z = 1000;
    Orbit orbit(45, 1010);
    Landmark pos;

    pos = orbit.position(10);
    while (orbit.update(10)) {
        for (int i = 0; i < 1000; ++i) {
            pos = orbit.position(10);
            x = pos.X;
            y = pos.Y;
            z = pos.Z;
            n_nan += run(N, 2 , 2, 1000);
        }
        break;
    }
    cout << "nan: " << n_nan << endl;
    cout << x << " " << y << " " << z << endl;
    return 0;
}