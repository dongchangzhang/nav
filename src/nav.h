#ifndef NAV_NAV
#define NAV_NAV

#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>

#include "landmark.h"

void gen_R_matrix(const Eigen::MatrixXd &vars, Eigen::Matrix3d &R);

std::vector<double> get_coefficient(double px, double py, double f, double x0, double y0, double mZ, double th, double wo, double ka, const Eigen::Matrix3d &R);

std::vector<double> get_mean_XYZ(double X, double Y, double Z, const Eigen::MatrixXd &vars, const Eigen::Matrix3d &R);

void get_matrix_A_and_B(const int N, const std::vector<Landmark> &data, const Eigen::MatrixXd &vars, const Eigen::Matrix3d &R, Eigen::MatrixXd &A, Eigen::MatrixXd &B);

void get_vector_L(int N, const std::vector<Landmark> &data, const Eigen::MatrixXd &vars, const Eigen::Matrix3d &R, Eigen::MatrixXd &L);

#endif