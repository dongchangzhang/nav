#ifndef NAV_NAV
#define NAV_NAV

#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>

#include "landmark.h"

void gen_R_matrix(Eigen::MatrixXd &vars, Eigen::Matrix3d &R);

std::vector<double> get_coefficient(double px, double py, double f, double x0, double y0, double mZ, double th, double wo, double ka, Eigen::Matrix3d &R);

std::vector<double> get_mean_XYZ(double X, double Y, double Z, Eigen::MatrixXd &vars, Eigen::Matrix3d &R);

void get_matrix_A_and_B(const int N, std::vector<Landmark> &data, Eigen::MatrixXd &vars, Eigen::Matrix3d &R, Eigen::MatrixXd &A, Eigen::MatrixXd &B);

void get_vector_L(int N, std::vector<Landmark> &data, Eigen::MatrixXd &vars, Eigen::Matrix3d &R, Eigen::MatrixXd &L);

#endif