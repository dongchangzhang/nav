#include <random>
#include <iostream>

#include "nav.h"
#include "tool.h"
#include "constants.h"


void gen_R_matrix(const Eigen::MatrixXd &vars, Eigen::Matrix3d &R) {
    double cth = cos(vars(3, 0));
    double sth = sin(vars(3, 0));
    double cwo = cos(vars(4, 0));
    double swo = sin(vars(4, 0));
    double cka = cos(vars(5, 0));
    double ska = sin(vars(5, 0));

    R << cth * cka - sth * swo * ska, -cth * ska - sth * swo * cka, -sth * cwo,
                           cwo * ska,                    cwo * cka, -swo,
         sth * cka + cth * swo * ska, -sth * ska + cth * swo * cka, cth * cwo;
}

std::vector<double> get_coefficient(double px, double py, double f, double x0, double y0, double mZ, double th, double wo, double ka, const Eigen::Matrix3d &R) {
    double a11 = (R(0, 0) * f + R(0, 2) * (px - x0)) / mZ;
    double a12 = (R(1, 0) * f + R(1, 2) * (px - x0)) / mZ;
    double a13 = (R(2, 0) * f + R(2, 2) * (px - x0)) / mZ;

    double a21 = (R(0, 1) * f + R(0, 2) * (py - y0)) / mZ;
    double a22 = (R(1, 1) * f + R(1, 2) * (py - y0)) / mZ;
    double a23 = (R(2, 1) * f + R(2, 2) * (py - y0)) / mZ;

    double a14 = (py - y0) * sin(wo) - (((px - x0) / f) * ((px - x0) * cos(ka) - (py - y0) * sin(ka)) + f * cos(ka)) * cos(wo);
    double a15 = -f * sin(ka) - ((px - x0) / f) * ((px - x0) * sin(ka) + (py - y0) * cos(ka));
    double a16 = py - y0;

    double a24 = -(px - x0) * sin(wo) - (((py - y0) / f) * ((px - x0) * cos(ka) - (py - y0) * sin(ka)) - f * sin(ka)) * cos(wo);
    double a25 = -f * cos(ka) - ((py - y0) / f) * ((px - x0) * sin(ka) + (py - y0) * cos(ka));
    double a26 = -(px - x0);

    double a17 = (px - x0) / f;
    double a27 = (py - y0) / f;

    double a18 = 1;
    double a28 = 0;
    double a19 = 0;
    double a29 = 1;

    return std::vector<double>{a11, a12, a13, a14, a15, a16, a17, a18, a19, a21, a22, a23, a24, a25, a26, a27, a28, a29};
}

std::vector<double> get_mean_XYZ(double X, double Y, double Z, const Eigen::MatrixXd &vars, const Eigen::Matrix3d &R) {
    //         a1 ->  R[0, 0]               b1 -> R[1, 0]                c1 -> R[2, 0]
    double mX = R(0, 0) * (X - vars(0, 0)) + R(1, 0) * (Y - vars(1, 0)) + R(2, 0) * (Z - vars(2, 0));
    //         a2 ->  R[0, 1]               b2 -> R[1, 1]                c2 -> R[2, 1]
    double mY = R(0, 1) * (X - vars(0, 0)) + R(1, 1) * (Y - vars(1, 0)) + R(2, 1) * (Z - vars(2, 0));
    //         a3 ->  R[0, 2]               b3 -> R[1, 2]                c3 -> R[2, 2]
    double mZ = R(0, 2) * (X - vars(0, 0)) + R(1, 2) * (Y - vars(1, 0)) + R(2, 2) * (Z - vars(2, 0));

    return std::vector<double>{mX, mY, mZ};
}

void get_matrix_A_and_B(const int N, const std::vector<Landmark> &data, const Eigen::MatrixXd &vars, const Eigen::Matrix3d &R, Eigen::MatrixXd &A, Eigen::MatrixXd &B) {
    A = Eigen::MatrixXd::Zero(2 * N, 5 * N);
    B = Eigen::MatrixXd::Zero(2 * N, N_VARS);
    double th = vars(3, 0), wo = vars(4, 0), ka = vars(5, 0), f = vars(6, 0), x0 = vars(7, 0), y0 = vars(8, 0);
    double px, py, mZ;
    std::vector<double> axx;
    for (int i = 0; i < N; ++i) {
        px = data[i].px;
        py = data[i].py;
        auto mean_XYZ = get_mean_XYZ(data[i].X, data[i].Y, data[i].Z, vars, R);
        mZ = mean_XYZ[2];

        axx = get_coefficient(px, py, f, x0, y0, mZ, th, wo, ka, R);
        // axx index:
        //     a11 a12 a13 a14 a15 a16 a17 a18 a19 a21 a22 a23 a24 a25 a26 a27 a28 a29
        //     0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17
        A.block<1, 5>(2 * i, i * 5) << 1, 0, axx[0], axx[1], axx[2];
        A.block<1, 5>(2 * i + 1, i * 5) << 0, 1, axx[9], axx[10], axx[11];

        B.block<1, N_VARS>(2 * i, 0) << axx[0], axx[1], axx[2], axx[3], axx[4], axx[5], axx[6], axx[7], axx[8];
        B.block<1, N_VARS>(2 * i + 1, 0)
                << axx[9], axx[10], axx[11], axx[12], axx[13], axx[14], axx[15], axx[16], axx[17];
    }
}

void get_vector_L(int N, const std::vector<Landmark> &data, const Eigen::MatrixXd &vars, const Eigen::Matrix3d &R, Eigen::MatrixXd &L) {
    L = Eigen::MatrixXd::Zero(2 * N, 1);
    std::vector<double> mean_XYZ;
    double px, py, mX, mY, mZ, x, y;
    for (int i = 0; i < N; ++i) {
        px = data[i].px;
        py = data[i].py;
        mean_XYZ = get_mean_XYZ(data[i].X, data[i].Y, data[i].Z, vars, R);
        mX = mean_XYZ[0];
        mY = mean_XYZ[1];
        mZ = mean_XYZ[2];
        x = vars(7, 0) - vars(6, 0) * mX / mZ;
        y = vars(8, 0) - vars(6, 0) * mY / mZ;

        L(i * 2, 0) = px - x;
        L(i * 2 + 1, 0) = py - y;
    }
}
