//
// Created by vover on 11/22/17.
//

#ifndef POPOV_RESIDUALS_H
#define POPOV_RESIDUALS_H

#include <vector>

double residual_max (std::vector<double> V, std::vector<double> H, double nu, double X, double T, int N, int M, int j);
double residual_L2h (std::vector<double> V, std::vector<double> H, double nu, double X, double T, int N, int M, int j);
double residual_12 (std::vector<double> V, std::vector<double> H, double nu, double X, double T, int N, int M, int j);

#endif //POPOV_RESIDUALS_H
