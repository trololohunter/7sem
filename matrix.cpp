//
// Created by vover on 11/22/17.
//

#include "matrix.h"

matrix::matrix()
{

}

matrix::matrix(int m)
{
    a.resize(m);
    b.resize(m-1);
    c.resize(m-1);
    rhs.resize(m);
    n = m;
}

matrix::~matrix()
{
    a.clear();
    b.clear();
    c.clear();
    rhs.clear();
}

void matrix::resize(int m)
{
    a.resize(m);
    b.resize(m-1);
    c.resize(m-1);
    rhs.resize(m);
    n = m;
}

int matrix::solver(std::vector<double>& x)
{
        double m;
        for (int i = 1; i < n; ++i)
        {
            m = c[i-1]/a[i-1];
            a[i] -= m*b[i-1];
            rhs[i] -= m*rhs[i-1];
        }

        x[n-1] = rhs[n-1]/a[n-1];

        for (int i = n - 2; i > -1; --i)
            x[i] = (rhs[i] - b[i]*x[i+1])/a[i];

        return 0;
}