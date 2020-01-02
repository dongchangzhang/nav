#include <cstdio>
#include <cmath>
#include <array>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigen>

using namespace std;
using namespace Eigen;
using std::default_random_engine;
using std::uniform_int_distribution;
const double pi = 3.141592653;


struct Landmark {
    double X, Y, Z, px, py;
    Landmark(double _X, double _Y, double _Z, double _px, double _py) : X(_X), Y(_Y), Z(_Z), px(_px), py(_py) {}
};

struct Parameters {
    // 0   1   2   3   4   5   6   7  8
    // Xs, Ys, Zs, th, wo, ka, f, x0, y0
    // m   m   m   de  de  de  mm um  um
    array<double, 9> raw_params, params, error;
    double _m = 1000, _mm = 1, _um = 0.0001;

    Parameters() {
    }

    Parameters(double _Xs, double _Ys, double _Zs, double _th, double _wo, double _ka, double _f, double _x0, double _y0)
    : raw_params({_Xs, _Ys, _Zs, _th * pi / 180.0 , _wo * pi / 180.0, _ka * pi / 180.0, _f, _x0, _y0}) {
        params[0] = raw_params[0] * _m;
        params[1] = raw_params[1] * _m;
        params[2] = raw_params[2] * _m;

        params[3] = raw_params[3];
        params[4] = raw_params[4];
        params[5] = raw_params[5];

        params[6] = raw_params[6] * _mm;

        params[7] = raw_params[7] * _um;
        params[8] = raw_params[8] * _um;
    }

    void add_error() {
        std::random_device rd;
        std::default_random_engine rng {rd()};
        std::normal_distribution<double> norm(0,1);

        error[0] = norm(rng) * 2 * _m / 10;
        error[1] = norm(rng) * 2 * _m / 10;
        error[2] = norm(rng) * 2 * _m / 10;
        error[3] = norm(rng);
        error[4] = norm(rng);
        error[5] = norm(rng);
        error[6] = norm(rng) * _mm;
        error[7] = norm(rng) * _um;
        error[8] = norm(rng) * _um;

        for (int i = 0; i < 9; ++i) {
            params[i] += error[i];
        }
    }

    void update(MatrixXd &x) {
        for (int i = 0; i < 9; ++i) {
            params[i] += x(i, 0);
        }
    }
    void print() {
        for (int i = 0; i < 9; ++i) {
            cout << params[i] << " ";
        }
        cout << endl;
    }
};

Matrix3d gen_R_matrix(Parameters &params) {
    Matrix3d R;
    double th = params.params[3];
    double wo = params.params[4];
    double ka = params.params[5];

    double cth = cos(th);
    double sth = sin(th);
    double cwo = cos(wo);
    double swo = sin(wo);
    double cka = cos(ka);
    double ska = sin(ka);

    R << cth*cka-sth*swo*ska,  -cth*ska-sth*swo*cka, -sth*cwo,
         cwo*ska, cwo*cka, -swo,
         sth*cka+cth*swo*ska, -sth*ska+cth*swo*cka, cth*cwo;
    return R;
}

vector<double> get_coefficient(double px, double py, double mZ, Parameters &params) {
    Matrix3d R = gen_R_matrix(params);
    double th = params.params[3];
    double wo = params.params[4];
    double ka = params.params[5];
    double f = params.params[6];
    double x0 = params.params[7];
    double y0 = params.params[8];

    double a1 = R(0,0);
    double a2 = R(0,1);
    double a3 = R(0,2);
    double b1 = R(1,0);
    double b2 = R(1,1);
    double b3 = R(1,2);
    double c1 = R(2,0);
    double c2 = R(2,1);
    double c3 = R(2,2);

    double i_mZ = 1.0/mZ;
    double a11 = i_mZ*(a1*f-a3*(px-x0));
    double a12 = i_mZ*(b1*f+b3*(px-x0));
    double a13 = i_mZ*(c1*f+c3*(px-x0));
    double a21 = i_mZ*(a2*f+a3*(py-y0));
    double a22 = i_mZ*(b2*f+b3*(py-y0));
    double a23 = i_mZ*(b3*f+c3*(py-y0));

    double a14 = (py-y0)*cos(wo) - ((px-x0)*((px-x0)*cos(ka)-(py-y0)*sin(ka))/f+f*cos(ka))*cos(wo);
    double a15 = -f*sin(ka) - (px-x0)*((px-x0)*sin(ka)+(py-y0)*cos(ka))/f;
    double a16 = py-y0;

    double a24 = -(px-x0)*sin(wo) - ((py-y0)*((px-x0)*cos(ka)-(py-y0)*sin(ka))/f-f*sin(ka))*cos(wo);
    double a25 = -f*cos(ka) - (py-y0)*((px-x0)*sin(ka)+(py-y0)*cos(ka))/f;
    double a26 = -(px-x0);

    double a17 = (px-x0)/f;
    double a27 = (py-y0)/f;

    double a18 = 1;
    double a28 = 0;
    double a19 = 0;
    double a29 = 1;
    vector<double> result{a11,a12,a13,a14,a15,a16,a17,a18,a19,a21,a22,a23,a24,a25,a26,a27,a28,a29};
    return result;

}

vector<double> get_mean_XYZ(Matrix3d &R, double X, double Y, double Z, Parameters &params) {
    double Xs = params.params[0], Ys = params.params[1], Zs = params.params[2];
    double a1 = R(0,0);
    double a2 = R(0,1);
    double a3 = R(0,2);
    double b1 = R(1,0);
    double b2 = R(1,1);
    double b3 = R(1,2);
    double c1 = R(2,0);
    double c2 = R(2,1);
    double c3 = R(2,2);

    double mX = a1*(X-Xs) + b1*(Y-Ys) + c1*(Z-Zs);
    double mY = a2*(X-Xs) + b2*(Y-Ys) + c2*(Z-Zs);
    double mZ = a3*(X-Xs) + b3*(Y-Ys) + c3*(Z-Zs);
    return vector<double> {mX, mY, mZ};
}

void get_matrixA_and_B(const int N, vector<Landmark>& data, Parameters &params, MatrixXd &A, MatrixXd &B) {
    A = MatrixXd::Zero(2*N, 5*N);
    B = MatrixXd::Zero(2*N,9);
    Matrix3d R = gen_R_matrix(params);

    int k = 0;
    double X, Y, Z, px, py, mX, mY, mZ;
    double a11,a12,a13,a14,a15,a16,a17,a18,a19,a21,a22,a23,a24,a25,a26,a27,a28,a29;
    for (int i = 0; i < N; ++i) {
        X = data[i].X;
        Y = data[i].Y;
        Z = data[i].Z;
        px = data[i].px;
        py = data[i].py;

        auto mean_XYZ = get_mean_XYZ(R,X,Y,Z,params);
        mX = mean_XYZ[0];
        mY = mean_XYZ[1];
        mZ = mean_XYZ[2];

        auto axx = get_coefficient(px,py,mZ, params);
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

        A.block<1, 5>(k+1, i * 5) << 0,1,a21,a22,a23;
        B.block<1, 9>(k,0) << a11,a12,a13,a14,a15,a16,a17,a18,a19;
        B.block<1, 9>(k+1,0) << a21,a22,a23,a24,a25,a26,a27,a28,a29;
        k = k + 2;
    }
}

void get_vector_L(int N, Parameters &params, vector<Landmark>& data, MatrixXd &L) {
    L = MatrixXd::Zero(2 * N, 1);
    Matrix3d R = gen_R_matrix(params);
    double X, Y, Z, px, py, mX, mY, mZ, x, y;
    for (int i = 0; i < N; ++i) {
        X = data[i].X;
        Y = data[i].Y;
        Z = data[i].Z;
        px = data[i].px;
        py = data[i].py;

        auto mean_XYZ = get_mean_XYZ(R, X, Y, Z, params);
        mX = mean_XYZ[0];
        mY = mean_XYZ[1];
        mZ = mean_XYZ[2];

        x = params.params[7] - params.params[6] * mX / mZ;
        y = params.params[8] - params.params[6] * mY / mZ;

        L(i * 2, 0) = px - x;
        L(i * 2 + 1, 0) = py - y;
    }
}

void get_ideal_data(int N, Parameters &params, vector<Landmark> &data) {
    auto R = gen_R_matrix(params);
    MatrixXd V = MatrixXd::Zero(N,3);
    V << 50000,50000,0,
    -50000,50000,0,
    -50000,-50000,0,
    50000,-50000,0,
    25000,0,10000,
    -25000,0,10000;

    double X, Y, Z, mX, mY, mZ, x, y;
    for (int i = 0; i < N; ++i) {
        X = V(i, 0);
        Y = V(i, 1);
        Z = V(i, 2);

        auto mean_XYZ = get_mean_XYZ(R, X, Y, Z, params);
        mX = mean_XYZ[0];
        mY = mean_XYZ[1];
        mZ = mean_XYZ[2];

        x = params.params[7] - params.params[6] * mX / mZ;
        y = params.params[8] - params.params[6] * mY / mZ;
        data.emplace_back(Landmark(X, Y, Z, x, y));

    }
}


int main() {
    // x
    double Xs = 2, Ys = 2, Zs = 1000, th = 2, wo = 2, ka = 2, f = 5, x0 = 6, y0 = 6;
    Parameters params_real(Xs, Ys, Zs, th, wo, ka, f, x0, y0);
    Parameters params = params_real;

    params.add_error();

    // data
    int N = 6;
    vector<Landmark> data;
    get_ideal_data(N, params, data);

    default_random_engine gen;
    std::normal_distribution<double> dis(0,1);

    for (auto &elem: data) {
        elem.X += 0.01* dis(gen) * 100;
        elem.Y += 0.01 * dis(gen) * 100;
        elem.Z += 0.01 * dis(gen) * 100;
    }


    Matrix3d R;
    MatrixXd A, B, L;
//    MatrixXd X = MatrixXd::Zero(9, 1);
//    MatrixXd V = MatrixXd::Zero(5 * N, 1);

    vector<double> err;

    printf("Base %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", Xs, Ys, Zs, th * 180 / pi, wo * 180 / pi, ka * 180 / pi, f, x0, y0);
    for (int i = 0; i < 2000; ++i) {
        get_matrixA_and_B(N, data, params, A, B);
        get_vector_L(N,params, data, L);

        MatrixXd W = A * A.transpose();
        MatrixXd x = (B.transpose() * W * B).inverse() * B.transpose() * W * L;
        MatrixXd K = (A * A.transpose()).inverse() * (B * x - L);
        MatrixXd v = A.transpose() * K;

        params.update(x);

        for (int j = 0; j < N; ++j) {
            data[j].X += v(j * 5 + 2);
            data[j].Y += v(j * 5 + 3);
            data[j].Z += v(j * 5 + 4);
        }
//        cout << x << endl;
//        cout << " -------------------- * ------------------------ " << endl;
        auto e = x.transpose() * x;
        cout << "step: " << i << " error: " << sqrt(e(0, 0)) << endl;
        cout << " -------------------- * ------------------------ " << endl;

        if (sqrt(e (0, 0)) < 1e-8) break;

    }
    params_real.print();
    params.print();
//    printf("Anws %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", Xs, Ys, Zs, th * 180 / pi, wo * 180 / pi, ka * 180 / pi, f, x0, y0);
//    printf("Real %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", Xs0, Ys0, Zs0, th0, wo0, ka0, f0, x00, y00);
//    printf("%lf", fabs(th - th0));


    return 0;
}

