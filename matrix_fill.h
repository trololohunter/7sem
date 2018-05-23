//
// Created by vover on 11/22/17.
//

#ifndef POPOV_MATRIX_FILL_H
#define POPOV_MATRIX_FILL_H

#include <vector>
#include "matrix.h"

void fill_matrix_0 (std::vector<double> &H, std::vector<double> &V, int n, int m, double x, double t);
void fill_matrix_0_1 (std::vector<double> &H, std::vector<double> &V, int n, int m, double x, double t);
void fill_matrix_0_2 (std::vector<double> &H, std::vector<double> &V, int n, int m, double x, double t);
void fill_matrix_H (matrix &H, std::vector<double> &Hpred, std::vector<double> &Vpred, int n, int m, double x, double t, int j);
void fill_matrix_V (matrix &V, std::vector<double> &Hpred, std::vector<double> &Vpred, std::vector<double> &Hnow, double mu, int n, int m, double x, double t, int j);
void resize_vectors (std::vector<double> &Hpred, std::vector<double> &Hnow, std::vector<double> &Vpred, std::vector<double> &Vnow, int m);

#endif //POPOV_MATRIX_FILL_H
