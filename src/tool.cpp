
#include "tool.h"
#include "nav.h"

#include <cstdio>
#include <random>
#include <iostream>
#include <fstream>

double radians(double degrees) {
    return degrees * PI / 180.0;
}

double degrees(double radians) {
    return radians * 180.0 / PI;
}

void get_ideal_data(int N, const Eigen::MatrixXd &vars, std::vector<Landmark> &data, std::default_random_engine &engine) {

    Eigen::Matrix3d R;
    gen_R_matrix(vars, R);
    double len = sqrt(vars(0, 0) * vars(0, 0) + vars(1, 0) * vars(1, 0) + vars(2, 0) * vars(2, 0));
    double boundary = len * tan(radians(22));
    std::normal_distribution<double> norm(0, 1/3.0); // u, stddev -> (-1, 1)

    // std::uniform_real_distribution<double> rand_xy(-boundary * 0.8, boundary * 0.8);
    // std::uniform_real_distribution<double> rand_z(-10, 10);

    double X, Y, Z, mX, mY, mZ, x, y;
    for (int i = 0; i < N; ++i) {
        X = norm(engine) * boundary;
        Y = norm(engine) * boundary;
        Z = norm(engine) * boundary;

        auto mean_XYZ = get_mean_XYZ(X, Y, Z, vars, R);

        mX = mean_XYZ[0];
        mY = mean_XYZ[1];
        mZ = mean_XYZ[2];

        x = vars(7, 0) - vars(6, 0) * mX / mZ;
        y = vars(8, 0) - vars(6, 0) * mY / mZ;
        data.emplace_back(Landmark(X, Y, Z, x, y));
    }
}

void add_noise(Eigen::MatrixXd &vars, std::vector<Landmark> &data, std::default_random_engine &engine) {
    // add error
    std::normal_distribution<double> norm(0, 1/3.0);

    // std::vector<double> vars_noise_unit{1 * _M, 1 * _M, 1 * _M, 
        // radians(5), radians(5), radians(5), 2 * _MM, 2 * _UM, 2 * _UM};
    std::vector<double> vars_noise_unit{
        NOISE_VARS_XYZ, NOISE_VARS_XYZ, NOISE_VARS_XYZ, 
        NOISE_VARS_ROTATE, NOISE_VARS_ROTATE, NOISE_VARS_ROTATE,
        NOISE_VARS_F, NOISE_VARS_X0Y0, NOISE_VARS_X0Y0
    };

    // for data
    for (auto &element: data) {
        element.X += norm(engine) * NOISE_DATA_XYZ;
        element.Y += norm(engine) * NOISE_DATA_XYZ;
        element.Z += norm(engine) * NOISE_DATA_XYZ;
    }

    // for vars
    for (int i = 0; i < N_VARS; ++i) {
        vars(i, 0) += vars_noise_unit[i] * norm(engine);
    }
}

void print(const Eigen::MatrixXd &vars) {

    auto _vars = vars;
    for (int k = 3; k <= 5; ++k) {
        _vars(k, 0) = degrees(vars(k, 0));
    }

    std::vector<std::string> titles{"Xs", "Ys", "Zs", "th", "wo", "ka", "f", "x0", "y0"};
    for (auto title: titles) {
        printf("%10s ", title.c_str());
    }
    printf("\n");

    for (int k = 0; k < 9; ++k) {
        printf("%10.5e ", _vars(k, 0));
    }
    printf("\n");
}

void dump(const std::vector<Landmark> &data, const Eigen::MatrixXd &real, const Eigen::MatrixXd &base, const Eigen::MatrixXd &solved, const std::vector<Eigen::MatrixXd> &errors, const std::vector<Eigen::MatrixXd> &vars_solved_history) {
    int id = 0, _id = -1;
    std::ifstream if_id(FID, std::ios::in);
    if (!if_id) {
        exit(-1);
    }
    while (!if_id.eof()) {
        if_id >> _id;
        id = std::max(id, _id);
    }
    id += 1;

#ifdef LOGGING
    std::cout << std::endl << "* saving the data into disk, the id is: " << id << std::endl;
#endif

    std::string vars_file = PVARS + std::to_string(id) + ".txt";
    std::string data_file = PDATA + std::to_string(id) + ".txt";
    std::string error_file = PERROR + std::to_string(id) + ".txt";
    std::string vars_history_file = PVARS_HISTORY + std::to_string(id) + ".txt";

    // add id
    std::ofstream of_id(FID, std::ios::app);
    of_id << id << std::endl;
    of_id.close();

    // write vars
    std::ofstream of_vars(vars_file);
    for (int k = 0; k < N_VARS; ++k) {
        of_vars << real(k, 0) << " ";
    }
    of_vars << std::endl;
    for (int k = 0; k < N_VARS; ++k) {
        of_vars << base(k, 0) << " ";
    }
    of_vars << std::endl;
    for (int k = 0; k < N_VARS; ++k) {
        of_vars << solved(k, 0) << " ";
    }
    of_vars << std::endl;
    of_vars.close();

    // write data
    std::ofstream of_data(data_file);
    for (auto landmark : data) {
        of_data << landmark.X << " " << landmark.Y << " " << landmark.Z << " " << landmark.px << " " << landmark.py
                << std::endl;
    }
    of_data.close();

    // write error
    std::ofstream of_error(error_file);
    for (auto error : errors) {
        for (int k = 0; k < N_VARS; ++k) {
            of_error << error(k, 0) << " ";
        }
        of_error << std::endl;
    }
    of_error.close();

    // write vars_solved_history
    std::ofstream of_vars_history(vars_history_file);
    for (auto var : vars_solved_history) {
        for (int k = 0; k < N_VARS; ++k) {
            of_vars_history << var(k, 0) << " ";
        }
        of_vars_history << std::endl;
    }
    of_vars_history.close();
}