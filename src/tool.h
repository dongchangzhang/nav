#ifndef NAV_TOOL
#define NAV_TOOL

#include <string>
#include <vector>
#include <random>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>

#include "landmark.h"
#include "constants.h"

double radians(double degrees);

double degrees(double radians);

void get_ideal_data(int N, const Eigen::MatrixXd &vars, std::vector<Landmark> &data, std::default_random_engine &engine);

void add_noise(Eigen::MatrixXd &vars, std::vector<Landmark> &data, std::default_random_engine &engine);

void print(const Eigen::MatrixXd &vars);

void dump(const std::vector<Landmark> &data, const Eigen::MatrixXd &real, const Eigen::MatrixXd &base, const Eigen::MatrixXd &solved, const std::vector<Eigen::MatrixXd> &errors, const std::vector<Eigen::MatrixXd> &vars_solved_history);

#endif