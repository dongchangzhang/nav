#include <cstdio>
#include <cmath>
#include <iostream>
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

Matrix3d gen_R_matrix(double th, double wo, double ka) {
    Matrix3d R;
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

vector<double> get_coefficient(double px, double py, double f, double x0, double y0, double mZ, double th, double wo, double ka, Matrix3d& R) {
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

vector<double> get_mean_XYZ(Matrix3d &R, double X, double Y, double Z, double Xs, double Ys, double Zs) {
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

void get_matrixA_and_B(const int N, vector<double>& VX, vector<double>& VY, vector<double>& VZ, vector<double>& Vpx, vector<double>& Vpy, double Xs, double Ys, double Zs, double f, double x0, double y0, double th, double wo, double ka, MatrixXd &A, MatrixXd &B) {
    A = MatrixXd::Zero(2*N, 5*N);
    B = MatrixXd::Zero(2*N,9);
    Matrix3d R = gen_R_matrix(th,wo,ka);
    int k = 0;
    double X, Y, Z, px, py, mX, mY, mZ;
    double a11,a12,a13,a14,a15,a16,a17,a18,a19,a21,a22,a23,a24,a25,a26,a27,a28,a29;
    for (int i = 0; i < N; ++i) {
        X = VX[i];
        Y = VY[i];
        Z = VZ[i];
        px = Vpx[i];
        py = Vpy[i];

        auto mean_XYZ = get_mean_XYZ(R,X,Y,Z,Xs,Ys,Zs);
        mX = mean_XYZ[0];
        mY = mean_XYZ[1];
        mZ = mean_XYZ[2];

        auto axx = get_coefficient(px,py,f,x0,y0,mZ,th,wo,ka,R);
        double a11,a12,a13,a14,a15,a16,a17,a18,a19,a21,a22,a23,a24,a25,a26,a27,a28,a29;
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

void get_vector_L(int N, Matrix3d &R, double Xs, double Ys, double Zs, double f, double x0, double y0, vector<double>& VX, vector<double> &VY, vector<double> &VZ, vector<double> &Vpx, vector<double> &Vpy, MatrixXd &L) {
    L = MatrixXd::Zero(2 * N, 1);
    double X, Y, Z, px, py, mX, mY, mZ, x, y;
    for (int i = 0; i < N; ++i) {
        X = VX[i];
        Y = VY[i];
        Z = VZ[i];
        px = Vpx[i];
        py = Vpy[i];

        auto mean_XYZ = get_mean_XYZ(R,X,Y,Z,Xs,Ys,Zs);
        mX = mean_XYZ[0];
        mY = mean_XYZ[1];
        mZ = mean_XYZ[2];
        x = x0 - f*mX/mZ;
        y = y0 - f*mY/mZ;

        L(i * 2, 0) = px - x;
        L(i * 2 + 1, 0) = py - y;
    }
}

void get_ideal_data(int N, double Xs, double Ys, double Zs, double th, double wo, double ka, double f, double x0, double y0,
        vector<double> &VX, vector<double> &VY, vector<double> &VZ, vector<double> &Vpx, vector<double> &Vpy) {
    auto R = gen_R_matrix(th,wo,ka);
    MatrixXd V = MatrixXd::Zero(N,3);
    V << 50,50,0,
    -50,50,0,
    -50,-50,0,
    50,-50,0,
    25,0,10,
    -25,0,10;

    double X, Y, Z, mX, mY, mZ, x, y;
    for (int i = 0; i < N; ++i) {
        X = V(i, 0);
        VX.push_back(X);
        Y = V(i, 1);
        VY.push_back(Y);
        Z = V(i, 2);
        VZ.push_back(Z);

        auto mean_XYZ = get_mean_XYZ(R,X,Y,Z,Xs,Ys,Zs);
        mX = mean_XYZ[0];
        mY = mean_XYZ[1];
        mZ = mean_XYZ[2];

        x = x0 - f*mX/mZ; Vpx.push_back(x);
        y = y0 - f*mY/mZ; Vpy.push_back(y);
    }


}



int main() {
    /*
    Matrix3d matrix = Matrix3d::Random();
    cout << matrix << endl;
    // transpose
    cout << " - * - " << endl;
    auto matrix_t = matrix.transpose();
    cout << matrix_t << endl;
    // inverse
    cout << " - * - " << endl;
    auto matrix_i = matrix.inverse();
    cout << matrix_i << endl;

    cout << " - * - " << endl;
    auto m = matrix_i * matrix_t;
    cout << m << endl;
     */
    double Xs = 2, Ys = 2, Zs = 1000;
    double th = 2 * pi / 180.0, wo = 2 * pi / 180.0, ka = 2 * pi / 180.0;
    double f = 0.5, x0 = 0.000006, y0 = 0.000006;

    double Xs0 = Xs, Ys0 = Ys, Zs0 = Zs;
    double th0 = th * 180 / pi, wo0 = wo * 180 / pi, ka0 = ka * 180 / pi;

    int N = 6;
    vector<double> VX, VY, VZ, Vpx, Vpy;
    get_ideal_data(N, Xs, Ys, Zs, th, wo, ka, f, x0, y0, VX, VY, VZ, Vpx, Vpy);

    default_random_engine gen;
    std::normal_distribution<double> dis(0,1);

    for (auto &elem: VX) {
        elem += 0.01* dis(gen);
    }
    for (auto &elem: VY) {
        elem += 0.01 * dis(gen);
    }
    for (auto &elem: VZ) {
        elem += 0.01 * dis(gen);
    }

    Matrix3d R;
    MatrixXd A, B, L;
    MatrixXd X = MatrixXd::Zero(9, 1);
    MatrixXd V = MatrixXd::Zero(5 * N, 1);
    Xs += 2 * dis(gen);
    Ys += 2 * dis(gen);
    Zs += 4 * dis(gen);
    th += 0.5 * dis(gen) * pi / 180.0;
    wo += 0.5 * dis(gen) * pi / 180.0;
    ka += 0.5 * dis(gen) * pi / 180.0;
    f = 0.7, x0 = 0.000008, y0 = 0.000001;


    vector<double> err;

    for (int i = 0; i < 2000; ++i) {
        get_matrixA_and_B(N, VX, VY, VZ, Vpx, Vpy, Xs, Ys, Zs, f, x0, y0, th, wo, ka, A, B);
        R = gen_R_matrix(th, wo, ka);
        get_vector_L(N, R, Xs, Ys, Zs, f, x0, y0, VX, VY, VZ, Vpx, Vpy, L);
        auto W = A * A.transpose();

        auto x = (B.transpose() * W * B).inverse() * B.transpose() * W * L;
        auto K = (A * A.transpose()).inverse() * (B * x - L);
        auto v = A.transpose() * K;

        Xs += x(0, 0);
        Ys += x(1, 0);
        Zs += x(2, 0);
        th += x(3, 0);
        wo += x(4, 0);
        ka += x(5, 0);
        f += x(6, 0);
        x0 += x(7, 0);
        y0 += x(8, 0);

        for (int j = 0; j < N; ++j) {
            VX[j] += v(j * 5 + 2);
            VY[j] += v(j * 5 + 3);
            VZ[j] += v(j * 5 + 4);
        }
        cout << x << endl;
        cout << " -------------------- * ------------------------ " << endl;
        auto e = x.transpose() * x;
        cout << "step: " << i << " error: " << sqrt(e(0, 0)) << endl;
        cout << " -------------------- * ------------------------ " << endl;

        if (sqrt(e (0, 0)) < 1e-8) break;

    }
    printf("Real %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", Xs, Ys, Zs, th * 180 / pi, wo * 180 / pi, ka * 180 / pi, f, x0, y0);
    printf("Anws %lf %lf %lf %lf %lf %lf\n", Xs0, Ys0, Zs0, th0, wo0, ka0);
    printf("%lf", fabs(th - th0));


    return 0;
}
