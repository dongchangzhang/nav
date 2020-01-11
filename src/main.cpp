// vhdsih
// 201912

#include <cstdio>
#include <cmath>
#include <ctime>
#include <chrono>
#include <iostream>
#include <fstream>
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

bool run(int n, const MatrixXd &vars, std::default_random_engine &engine) {
    Matrix3d R;
    MatrixXd A, B, L, W, x, K, v; 
    MatrixXd diff, error, loss, vars_real = vars, vars_base, vars_to_solved = vars;

    vector<MatrixXd> errors;
    vector<MatrixXd> vars_solved_history;
    vector<Landmark> data, data_real, data_noised;

    get_ideal_data(n, vars, data, engine);
    // backup 
    data_real = data;

    // add noise for vars and data
    add_noise(vars_to_solved, data, engine);

    // backup
    vars_base = vars_to_solved;
    data_noised = data;

    vars_solved_history.push_back(vars_to_solved);

    bool is_nan = false;
    for (int i = 0; i < MAX_ITERS; ++i) {
        gen_R_matrix(vars_to_solved, R);
        get_matrix_A_and_B(n, data, vars_to_solved, R, A, B);
        get_vector_L(n, data, vars_to_solved, R, L);

        W = A * A.transpose();
        x = (B.transpose() * W * B).inverse() * B.transpose() * W * L;
        K = (A * A.transpose()).inverse() * (B * x - L);
        v = A.transpose() * K;

        // update vars with dx
        vars_to_solved += x;

        // update data with dv
        for (int k = 0; k < n; ++k) {
            data[k].xyz.x += v(k * 5 + 2, 0);
            data[k].xyz.y += v(k * 5 + 3, 0);
            data[k].xyz.z += v(k * 5 + 4, 0);
        }

        // calculating loss
        loss = x.transpose() * x;

        // is a bad calculate procedure ?
        if (loss(0, 0) > 1e3 || isnan(loss(0, 0))) {
            cerr << "nan" << endl;
            return false;
        }

        // calculating real error (vars_real)
        diff = vars_real - vars_to_solved;
        //
        diff(3, 0) = degrees(diff(3, 0));
        diff(4, 0) = degrees(diff(4, 0));
        diff(5, 0) = degrees(diff(5, 0));

        diff(6, 0) *= 1e3;
        diff(7, 0) *= 1e6;
        diff(8, 0) *= 1e6;

        error = diff.transpose() * diff;

        // backup
        errors.emplace_back(diff);
        vars_solved_history.push_back(vars);

#ifdef LOGGING
        if (i % N_PRINT == 0)
            printf("step: %5d   loss: %10e   real_error: %10e\n", i, sqrt(loss(0, 0)), sqrt(error(0, 0)));
#endif

        // when to leave
        if (sqrt(loss(0, 0)) < LIMIT)
            break;
    }

#ifdef DUMP_DATA
    dump(data_real, vars_real, vars_base, vars_to_solved, errors, vars_solved_history);
#endif

    // show result
#ifdef LOGGING
    cout << "solved:" << endl;
    print(vars_to_solved);
#endif
    return true;

}

int main() {
    int idx = 0, n_nan = 0;
    Point3 position;
    Orbit orbit(60, 1010);

    std::default_random_engine engine;
    engine.seed((unsigned) time(0));
    std::normal_distribution<double> norm(0, 1.0/3.0); // u, stddev -> (-1, 1)

    MatrixXd vars = MatrixXd::Zero(N_VARS, 1); // vars {Xs, Ys, Zs, th, wo, ka, f, x0, y0};

#ifdef CALCU_TIME
    double len, _c;
    vector<pair<double, double>> times;
    chrono::steady_clock::time_point _a, _b;
#endif

    chrono::steady_clock::time_point time_begin = chrono::steady_clock::now();
    while (orbit.update(5)) {

#ifdef CALCU_TIME
        _a = chrono::steady_clock::now();
#endif
        position = orbit.position(10);
        // assign new vars
        vars(0, 0) = position.x;
        vars(1, 0) = position.y;
        vars(2, 0) = position.z;

        vars(3, 0) = radians(2 + norm(engine));
        vars(4, 0) = radians(2 + norm(engine));
        vars(5, 0) = radians(2 + norm(engine));

        vars(6, 0) = 7.6 * _MM;
        vars(7, 0) = 6 * _UM;
        vars(8, 0) = 6 * _UM;

        for (int i = 0; i < N_REPEAT; ++i) {
            ++idx;
            auto yes = run(10, vars, engine);
            if (!yes) ++n_nan;
        }
#ifdef CALCU_TIME

        _b = chrono::steady_clock::now();
        _c = chrono::duration<double, std::micro>(_b - _a).count() / N_REPEAT;
        len = sqrt(vars(0, 0) * vars(0, 0) + vars(1, 0) * vars(1, 0) + vars(2, 0) * vars(2, 0));
        times.emplace_back(pair<double, double>(len, _c));
#endif

#ifdef RUN_NOTICE
        if (idx % 200 == 0) cout << "." << endl;
#endif
    }
    chrono::steady_clock::time_point time_end = chrono::steady_clock::now();
    double us = chrono::duration<double, std::micro>(time_end - time_begin).count();
    printf("\n[INFO] Iteration count: %8d  |  Bad iter count: %8d  |  Time used: %10lf us\n", idx, n_nan, us);

#ifdef CALCU_TIME
    ofstream of("times.txt");
    for (auto t: times) {
        of << t.first << " " << t.second << endl;
    }
    of.close();
#endif
    
    return 0;
}