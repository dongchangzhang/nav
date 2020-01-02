#ifndef NAV_TOOL
#define NAV_TOOL

#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>

#include "landmark.h"
#include "constants.h"

double radians(double degrees);

double degrees(double radians);

void get_ideal_data(int N, Eigen::MatrixXd &vars, std::vector<Landmark> &data, double len);

void add_noise(Eigen::MatrixXd &vars, std::vector<Landmark> &data);

void print(Eigen::MatrixXd &vars, std::string &label);

void dump(std::vector<Landmark> &data, Eigen::MatrixXd &real, Eigen::MatrixXd &base, Eigen::MatrixXd &solved, std::vector<Eigen::MatrixXd> &errors, std::vector<Eigen::MatrixXd> &vars_solved_history);

#endif