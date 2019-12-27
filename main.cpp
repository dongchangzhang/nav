#include <cstdio>
#include <cmath>
#include <array>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>

using namespace std;
using namespace Eigen;
const double pi = 3.141592653;


struct Landmark {
    double X, Y, Z, px, py;

    Landmark(double _X, double _Y, double _Z, double _px, double _py) : X(_X), Y(_Y), Z(_Z), px(_px), py(_py) {}
};


Matrix3d gen_R_matrix(array<double, 9> &vars) {
    Matrix3d R;
    double cth = cos(vars[3]);
    double sth = sin(vars[3]);
    double cwo = cos(vars[4]);
    double swo = sin(vars[4]);
    double cka = cos(vars[5]);
    double ska = sin(vars[5]);

    R << cth * cka - sth * swo * ska, -cth * ska - sth * swo * cka, -sth * cwo,
            cwo * ska, cwo * cka, -swo,
            sth * cka + cth * swo * ska, -sth * ska + cth * swo * cka, cth * cwo;
    return R;
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

    vector<double> result{a11, a12, a13, a14, a15, a16, a17, a18, a19, a21, a22, a23, a24, a25, a26, a27, a28, a29};
    return result;
}

vector<double> get_mean_XYZ(double X, double Y, double Z, array<double, 9> &vars, Matrix3d &R) {
    double a1 = R(0, 0);
    double a2 = R(0, 1);
    double a3 = R(0, 2);
    double b1 = R(1, 0);
    double b2 = R(1, 1);
    double b3 = R(1, 2);
    double c1 = R(2, 0);
    double c2 = R(2, 1);
    double c3 = R(2, 2);

    double mX = a1 * (X - vars[0]) + b1 * (Y - vars[1]) + c1 * (Z - vars[2]);
    double mY = a2 * (X - vars[0]) + b2 * (Y - vars[1]) + c2 * (Z - vars[2]);
    double mZ = a3 * (X - vars[0]) + b3 * (Y - vars[1]) + c3 * (Z - vars[2]);

    return vector<double>{mX, mY, mZ};
}

void
get_matrixA_and_B(const int N, vector<Landmark> &data, array<double, 9> &vars, Matrix3d &R, MatrixXd &A, MatrixXd &B) {
    A = MatrixXd::Zero(2 * N, 5 * N);
    B = MatrixXd::Zero(2 * N, 9);
    double th = vars[3], wo = vars[4], ka = vars[5], f = vars[6], x0 = vars[7], y0 = vars[8];
    int k = 0;
    double X, Y, Z, px, py, mX, mY, mZ;
    double a11, a12, a13, a14, a15, a16, a17, a18, a19, a21, a22, a23, a24, a25, a26, a27, a28, a29;
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

        auto axx = get_coefficient(px, py, f, x0, y0, mZ, th, wo, ka, R);
        a11 = axx[0];
        a12 = axx[1];
        a13 = axx[2];
        a14 = axx[3];
        a15 = axx[4];
        a16 = axx[5];
        a17 = axx[6];
        a18 = axx[7];
        a19 = axx[8];
        a21 = axx[9];
        a22 = axx[10];
        a23 = axx[11];
        a24 = axx[12];
        a25 = axx[13];
        a26 = axx[14];
        a27 = axx[15];
        a28 = axx[16];
        a29 = axx[17];
        A.block<1, 5>(k, i * 5) << 1, 0, a11, a12, a13;

        A.block<1, 5>(k + 1, i * 5) << 0, 1, a21, a22, a23;
        B.block<1, 9>(k, 0) << a11, a12, a13, a14, a15, a16, a17, a18, a19;
        B.block<1, 9>(k + 1, 0) << a21, a22, a23, a24, a25, a26, a27, a28, a29;
        k = k + 2;
    }
}

void get_vector_L(int N, vector<Landmark> &data, array<double, 9> &vars, Matrix3d &R, MatrixXd &L) {
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
        x = vars[7] - vars[6] * mX / mZ;
        y = vars[8] - vars[6] * mY / mZ;

        L(i * 2, 0) = px - x;
        L(i * 2 + 1, 0) = py - y;
    }
}

void get_ideal_data(int N, array<double, 9> vars, vector<Landmark> &data) {
    auto R = gen_R_matrix(vars);
    default_random_engine engine;
    uniform_real_distribution<double> real_rand(-50, 50);

    double X, Y, Z, mX, mY, mZ, x, y;
    for (int i = 0; i < N; ++i) {
        X = real_rand(engine);
        Y = real_rand(engine);
        Z = real_rand(engine);

        auto mean_XYZ = get_mean_XYZ(X, Y, Z, vars, R);
        mX = mean_XYZ[0];
        mY = mean_XYZ[1];
        mZ = mean_XYZ[2];

        x = vars[7] - vars[6] * mX / mZ;
        y = vars[8] - vars[6] * mY / mZ;
        data.emplace_back(Landmark(X, Y, Z, x, y));

    }
}


int main() {
    double Xs = 2, Ys = 2, Zs = 1000, th = 2 * pi / 180.0, wo = 2 * pi / 180.0, ka =
            2 * pi / 180.0, f = 0.005, x0 = 0.000006, y0 = 0.000006;

    array<double, 9> vars{Xs, Ys, Zs, th, wo, ka, f, x0, y0};

    // real params
    array<double, 9> real = vars;

    // degrees to radians
    real[3] *= 180.0 / pi;
    real[4] *= 180.0 / pi;
    real[5] *= 180.0 / pi;

    // the number of point
    int N = 10;

    // data
    vector<Landmark> data;
    get_ideal_data(N, vars, data);

    // add error
    default_random_engine engine;
    std::normal_distribution<double> norm(0, 1);
    for (auto &elem: data) {
        elem.X += 0.01 * norm(engine);
        elem.Y += 0.01 * norm(engine);
        elem.Z += 0.01 * norm(engine);
    }

    vars[0] += 5 * norm(engine);
    vars[1] += 5 * norm(engine);
    vars[2] += 5 * norm(engine);
    vars[3] += 0.5 * norm(engine) * pi / 180.0;
    vars[4] += 0.5 * norm(engine) * pi / 180.0;
    vars[5] += 0.5 * norm(engine) * pi / 180.0;
    vars[6] = 0.007;
    vars[7] = 0.000008;
    vars[8] = 0.000001;

    // params-base
    auto base = vars;

    Matrix3d R;
    MatrixXd A, B, L;
    vector<double> error;

    for (int i = 0; i < 20000; ++i) {
        Matrix3d R = gen_R_matrix(vars);
        get_matrixA_and_B(N, data, vars, R, A, B);
        get_vector_L(N, data, vars, R, L);

        auto W = A * A.transpose();
        auto x = (B.transpose() * W * B).inverse() * B.transpose() * W * L;
        auto K = (A * A.transpose()).inverse() * (B * x - L);
        auto v = A.transpose() * K;

        // update params
        for (int k = 0; k < 9; ++k) {
            vars[k] += x(k, 0);
        }

        // update data
        for (int j = 0; j < N; ++j) {
            data[j].X += v(j * 5 + 2);
            data[j].Y += v(j * 5 + 3);
            data[j].Z += v(j * 5 + 4);
        }

        // cal error
        auto e = x.transpose() * x;
        if (i % 50 == 0)
            cout << "step: " << i << " error: " << sqrt(e(0, 0)) << endl;

        // rule
        if (sqrt(e(0, 0)) < 1e-10)
            break;

    }

    vars[3] *= 180.0 / pi;
    vars[4] *= 180.0 / pi;
    vars[5] *= 180.0 / pi;

    cout << "::Solved" << endl;
    for (auto var: vars) cout << var << " ";
    cout << endl;

    base[3] *= 180.0 / pi;
    base[4] *= 180.0 / pi;
    base[5] *= 180.0 / pi;
    cout << "::Base" << endl;
    for (auto var: base) cout << var << " ";
    cout << endl;
    cout << "::Real" << endl;
    for (auto var: real) cout << var << " ";
    cout << endl;

    return 0;
}

