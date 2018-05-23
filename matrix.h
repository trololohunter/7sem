//
// Created by vover on 11/22/17.
//

#ifndef POPOV_MATRIX_H
#define POPOV_MATRIX_H

#include <vector>

class matrix {
public:
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> rhs;
    int n;

    matrix();
    ~matrix();
    matrix(int m);
    void resize(int m);
    int solver (std::vector<double> &x);


};


#endif //POPOV_MATRIX_H
