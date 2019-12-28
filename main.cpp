#include <cstdio>
#include <cmath>
#include <ctime>
#include <array>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <zconf.h>

using namespace std;
using namespace Eigen;

const int N_VARS = 9;
const int MAX_ITERS = 1000;
const double PI = 3.141592653;
const double LIMIT = 1e-9;
const double unit_m = 1, unit_cm = 1e-2, unit_mm = 1e-3, unit_um = 1e-6;


const string id_file = "../records/id.txt";
const string vars_path = "../records/vars/";
const string data_path = "../records/data/";
const string error_path = "../records/errors/";

struct Landmark {
    double X, Y, Z, px, py;

    Landmark(double _X, double _Y, double _Z, double _px, double _py) : X(_X), Y(_Y), Z(_Z), px(_px), py(_py) {}
};

double radians(double degrees) {
    return degrees * PI / 180.0;
}

double degrees(double radians) {
    return radians * 180.0 / PI;
}

void gen_R_matrix(MatrixXd &vars, Matrix3d &R) {
    double cth = cos(vars(3, 0));
    double sth = sin(vars(3, 0));
    double cwo = cos(vars(4, 0));
    double swo = sin(vars(4, 0));
    double cka = cos(vars(5, 0));
    double ska = sin(vars(5, 0));

    R << cth * cka - sth * swo * ska, -cth * ska - sth * swo * cka, -sth * cwo,
            cwo * ska, cwo * cka, -swo,
            sth * cka + cth * swo * ska, -sth * ska + cth * swo * cka, cth * cwo;

}

vector<double> get_coefficient(double px, double py, double f, double x0, double y0,
                               double mZ, double th, double wo, double ka, Matrix3d &R) {
    double a1 = R(0, 0);
    double a2 = R(0, 1);
    double a3 = R(0, 2);
    double b1 = R(1, 0);
    double b2 = R(1, 1);
    double b3 = R(1, 2);
    double c1 = R(2, 0);
    double c2 = R(2, 1);
    double c3 = R(2, 2);

    double i_mZ = 1.0 / mZ;
    double a11 = i_mZ * (a1 * f - a3 * (px - x0));
    double a12 = i_mZ * (b1 * f + b3 * (px - x0));
    double a13 = i_mZ * (c1 * f + c3 * (px - x0));
    double a21 = i_mZ * (a2 * f + a3 * (py - y0));
    double a22 = i_mZ * (b2 * f + b3 * (py - y0));
    double a23 = i_mZ * (b3 * f + c3 * (py - y0));

    double a14 =
            (py - y0) * cos(wo) - ((px - x0) * ((px - x0) * cos(ka) - (py - y0) * sin(ka)) / f + f * cos(ka)) * cos(wo);
    // in word
    //double a14 =
    //       (py - y0) * sin(wo) - ((px - x0) * ((px - x0) * cos(ka) - (py - y0) * sin(ka)) / f + f * cos(ka)) * cos(wo);
    double a15 = -f * sin(ka) - (px - x0) * ((px - x0) * sin(ka) + (py - y0) * cos(ka)) / f;
    double a16 = py - y0;

    double a24 = -(px - x0) * sin(wo) -
                 ((py - y0) * ((px - x0) * cos(ka) - (py - y0) * sin(ka)) / f - f * sin(ka)) * cos(wo);
    double a25 = -f * cos(ka) - (py - y0) * ((px - x0) * sin(ka) + (py - y0) * cos(ka)) / f;
    double a26 = -(px - x0);

    double a17 = (px - x0) / f;
    double a27 = (py - y0) / f;

    double a18 = 1;
    double a28 = 0;
    double a19 = 0;
    double a29 = 1;

    return vector<double>{a11, a12, a13, a14, a15, a16, a17, a18, a19, a21, a22, a23, a24, a25, a26, a27, a28, a29};
}

vector<double> get_mean_XYZ(double X, double Y, double Z, MatrixXd &vars, Matrix3d &R) {
    double a1 = R(0, 0);
    double a2 = R(0, 1);
    double a3 = R(0, 2);
    double b1 = R(1, 0);
    double b2 = R(1, 1);
    double b3 = R(1, 2);
    double c1 = R(2, 0);
    double c2 = R(2, 1);
    double c3 = R(2, 2);

    double mX = a1 * (X - vars(0, 0)) + b1 * (Y - vars(1, 0)) + c1 * (Z - vars(2, 0));
    double mY = a2 * (X - vars(0, 0)) + b2 * (Y - vars(1, 0)) + c2 * (Z - vars(2, 0));
    double mZ = a3 * (X - vars(0, 0)) + b3 * (Y - vars(1, 0)) + c3 * (Z - vars(2, 0));

    return vector<double>{mX, mY, mZ};
}

void
get_matrixA_and_B(const int N, vector<Landmark> &data, MatrixXd &vars, Matrix3d &R, MatrixXd &A, MatrixXd &B) {
    A = MatrixXd::Zero(2 * N, 5 * N);
    B = MatrixXd::Zero(2 * N, N_VARS);
    double th = vars(3, 0), wo = vars(4, 0), ka = vars(5, 0), f = vars(6, 0), x0 = vars(7, 0), y0 = vars(8, 0);
    double  px, py, mZ;
    for (int i = 0; i < N; ++i) {
        px = data[i].px;
        py = data[i].py;
        auto mean_XYZ = get_mean_XYZ(data[i].X, data[i].Y, data[i].Z, vars, R);
        mZ = mean_XYZ[2];

        auto axx = get_coefficient(px, py, f, x0, y0, mZ, th, wo, ka, R);
        // axx index
        // a11 a12 a13 a14 a15 a16 a17 a18 a19 a21 a22 a23 a24 a25 a26 a27 a28 a29
        // 0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17
        A.block<1, 5>(2 * i, i * 5) << 1, 0, axx[0], axx[1], axx[2];
        A.block<1, 5>(2 * i + 1, i * 5) << 0, 1, axx[9], axx[10], axx[11];

        B.block<1, N_VARS>(2 * i, 0) << axx[0], axx[1], axx[2], axx[3], axx[4], axx[5], axx[6], axx[7], axx[8];
        B.block<1, N_VARS>(2 * i + 1, 0) << axx[9], axx[10], axx[11], axx[12], axx[13], axx[14], axx[15], axx[16], axx[17];
    }
}

void get_vector_L(int N, vector<Landmark> &data, MatrixXd &vars, Matrix3d &R, MatrixXd &L) {
    L = MatrixXd::Zero(2 * N, 1);
    double X, Y, Z, px, py, mX, mY, mZ, x, y;
    for (int i = 0; i < N; ++i) {
        X = data[i].X;
        Y = data[i].Y;
        Z = data[i].Z;
        px = data[i].px;
        py = data[i].py;

        auto mean_XYZ = get_mean_XYZ(X, Y, Z, vars, R);
        mX = mean_XYZ[0];
        mY = mean_XYZ[1];
        mZ = mean_XYZ[2];
        x = vars(7, 0) - vars(6, 0) * mX / mZ;
        y = vars(8, 0) - vars(6, 0) * mY / mZ;

        L(i * 2, 0) = px - x;
        L(i * 2 + 1, 0) = py - y;
    }
}

void get_ideal_data(int N, MatrixXd vars, vector<Landmark> &data) {
    Matrix3d R;
    gen_R_matrix(vars, R);
    default_random_engine engine;
    engine.seed((unsigned) time(0));
    // engine.seed(3);
    uniform_real_distribution<double> real_rand(-50, 50);

    double X, Y, Z, mX, mY, mZ, x, y;
    cout << " -- points -- " << endl;
    //vector<Landmark> tmp{Landmark(50, 50, 0, 0, 0), Landmark(-50, 50, 0, 0, 0),
    // Landmark(-50, -50, 0, 0, 0), Landmark(50, -50, 0, 0, 0), Landmark(25, 0, 10, 0, 0), Landmark(-25, 0, 10, 0, 0)};
    for (int i = 0; i < N; ++i) {
        X = real_rand(engine) * unit_m;
        Y = real_rand(engine) * unit_m;
        Z = real_rand(engine) * unit_m;

        printf("%10.5lf %10.5lf %10.5lf\n", X, Y, Z);

        auto mean_XYZ = get_mean_XYZ(X, Y, Z, vars, R);
        mX = mean_XYZ[0];
        mY = mean_XYZ[1];
        mZ = mean_XYZ[2];

        x = vars(7, 0) - vars(6, 0) * mX / mZ;
        y = vars(8, 0) - vars(6, 0) * mY / mZ;
        data.emplace_back(Landmark(X, Y, Z, x, y));
    }
    cout << " -- distance -- " << endl;
    for (int i = 0; i < data.size(); ++i) {
        for (int j = i + 1; j < data.size(); ++j) {
            double a = data[i].X - data[j].X;
            double b = data[i].Y - data[j].Y;
            double c = 0;
            double d = sqrt(a * a + b * b + c * c);
            printf("%10.4lf ", d);
        }
        cout << endl;
    }
}

bool add_noise(MatrixXd &vars, vector<Landmark> &data) {
    // add error
    default_random_engine engine;
    engine.seed((unsigned) time(0));
//    engine.seed(3);
    std::normal_distribution<double> norm(0, 1);

    // for data
    for (auto &element: data) {
        element.X += norm(engine) * 1 * unit_cm;
        element.Y += norm(engine) * 1 * unit_cm;
        element.Z += norm(engine) * 1 * unit_cm;
    }

    // for vars
    vars(0, 0) += 2 * norm(engine);
    vars(1, 0) += 2 * norm(engine);
    vars(2, 0) += 4 * norm(engine);

    vars(3, 0) += 5 * radians(1) * norm(engine);
    vars(4, 0) += 5 * radians(1) * norm(engine);
    vars(5, 0) += 5 * radians(1) * norm(engine);

    vars(6, 0) += 2 * unit_mm * norm(engine);

    vars(7, 0) += 2 * unit_um * norm(engine);
    vars(8, 0) += 2 * unit_um * norm(engine);

    return !(vars.rows() != N_VARS && vars.cols() != 1);
}

void print(MatrixXd &vars, string &label) {

    for (int k = 3; k <= 5; ++k) {
        vars(k, 0) = degrees(vars(k, 0));
    }

    cout << endl << " -- " << label << " -- " << endl;

    vector<string> titles{"Xs", "Ys", "Zs", "th", "wo", "ka", "f", "x0", "y0"};
    for (auto title: titles) {
        printf("%10s ", title.c_str());
    }
    printf("\n");

    for (int k = 0; k < 9; ++k) {
        printf("%10.5e ", vars(k, 0));
    }
    cout << endl;
    for (int k = 3; k <= 5; ++k) {
        vars(k, 0) = radians(vars(k, 0));
    }
}

void dump(vector<Landmark> &data, MatrixXd &real, MatrixXd &base, MatrixXd &solved, vector<MatrixXd> &errors) {
    //
    int id = 0, _id = -1;
    string tmp;
    ifstream if_id(id_file, ios::in);
    if (!if_id) {
        exit(-1);
    }
    while (!if_id.eof()) {
        if_id >> _id;
        id = max(id, _id);
    }
    id += 1;
    cout << endl << "data id:" << id << endl;

    string vars_file = vars_path + to_string(id) + ".txt";
    string data_file = data_path + to_string(id) + ".txt";
    string error_file = error_path + to_string(id) + ".txt";

    // add id
    ofstream of_id(id_file, ios::app);
    of_id << id << endl;
    of_id.close();

    // write vars
    ofstream of_vars(vars_file);
    for (int k = 0; k < N_VARS; ++k) {
        of_vars << real(k, 0) << " ";
    }
    of_vars << endl;
    for (int k = 0; k < N_VARS; ++k) {
        of_vars << base(k, 0) << " ";
    }
    of_vars << endl;
    for (int k = 0; k < N_VARS; ++k) {
        of_vars << solved(k, 0) << " ";
    }
    of_vars << endl;
    of_vars.close();

    // write data
    ofstream of_data(data_file);
    for (auto landmark : data) {
        of_data << landmark.X << " " << landmark.Y << " " << landmark.Z << " " << landmark.px << " " << landmark.py
                << endl;
    }
    of_data.close();

    // write error
    ofstream of_error(error_file);
    for (auto error : errors) {
        for (int k = 0; k < N_VARS; ++k) {
            of_error << error(k, 0) << " ";
        }
        of_error << endl;
    }
    of_error.close();
}

int main() {
    // the number of point
    int N = 20;

    double Xs = 2 * unit_m, Ys = 2 * unit_m, Zs = 1000 * unit_m,
            th = radians(2), wo = radians(2), ka = radians(2), f = 5 * unit_mm, x0 = 6 * unit_um, y0 = 6 * unit_um;

    // array<double, 9> vars{Xs, Ys, Zs, th, wo, ka, f, x0, y0};
    MatrixXd vars = MatrixXd::Zero(N_VARS, 1);
    vars << Xs, Ys, Zs, th, wo, ka, f, x0, y0;

    // save real params
    MatrixXd history_vars_real = vars;

    // data
    vector<Landmark> data;
    get_ideal_data(N, vars, data);

    // add error
    add_noise(vars, data);

    // save params-base
    MatrixXd history_vars_base = vars;

    Matrix3d R;
    MatrixXd A, B, L, error, error_matrix, loss;
    MatrixXd W, x, K, v;
    vector<MatrixXd> errors;

    string tmp = "areal";
    print(history_vars_real, tmp);
    tmp = "abase";
    print(history_vars_base, tmp);

    cout << " -- start -- " << endl;

    for (int i = 0; i < MAX_ITERS; ++i) {
        gen_R_matrix(vars, R);
        get_matrixA_and_B(N, data, vars, R, A, B);
        get_vector_L(N, data, vars, R, L);

        W = A * A.transpose();
        x = (B.transpose() * W * B).inverse() * B.transpose() * W * L;
        K = (A * A.transpose()).inverse() * (B * x - L);
        v = A.transpose() * K;

        // update params
        vars += x;

        // update data
        for (int j = 0; j < N; ++j) {
            data[j].X += v(j * 5 + 2, 0);
            data[j].Y += v(j * 5 + 3, 0);
            data[j].Z += v(j * 5 + 4, 0);
        }

        // calculate error
        loss = x.transpose() * x;
        error_matrix = history_vars_real - vars;
        error = error_matrix.transpose() * error_matrix;
        errors.emplace_back(error_matrix);

        if (i % 10 == 0)
            printf("step: %5d   loss: %10e   error: %10e\n", i, sqrt(loss(0, 0)), sqrt(error(0, 0)));

        // break rule
        if (sqrt(loss(0, 0)) < LIMIT)
            break;
    }

    dump(data, history_vars_real, history_vars_base, vars, errors);

    // show result
    string label = "solved";
    print(vars, label);
    label = "base";
    print(history_vars_base, label);
    label = "real";
    print(history_vars_real, label);

    return 0;
}

