//
// Created by vover on 11/22/17.
//

#include "matrix_fill.h"
#include "functions.h"

void resize_vectors (std::vector<double> &Hpred, std::vector<double> &Hnow, std::vector<double> &Vpred, std::vector<double> &Vnow, int m)
{
    Hpred.resize(m+1);
    Hnow.resize(m+1);
    Vpred.resize(m+1);
    Vnow.resize(m-1);
}

void fill_matrix_0 (std::vector<double> &H, std::vector<double> &V, int n, int m, double x, double t)
{
    int i = 0;
    double h = x/m;

    for (i = 0; i < m+1; ++i)
    {
        H[i] = ro(h*i,0);
        V[i] = u(h*i, 0);
    }
    V[0] = 0;
    V[m] = 0;

    return;
}

void fill_matrix_0_1 (std::vector<double> &H, std::vector<double> &V, int n, int m, double x, double t)
{
    int i = 0;
    double h = x / m;
    for (i = 0; i < m+1; ++i)
    {
        if (h * i < 4.5 || h*i > 5.5)
            H[i] = 1;
        else H[i] = 2;

        V[i] = 0;
    }

    return;
}

void fill_matrix_0_2 (std::vector<double> &H, std::vector<double> &V, int n, int m, double x, double t)
{
    int i = 0;
    double h = x / m;
    for (i = 0; i < m+1; ++i)
    {
        if (h * i < 4.5 || h*i > 5.5)
            V[i] = 0;
        else V[i] = 1;

        H[i] = 1;
    }

    return;
}

void fill_matrix_H (matrix &H, std::vector<double> &Hpred, std::vector<double> &Vpred,
                    int n, int m, double x, double t, int j)
{
    int i = 0;
    double h = x/m;
    double tau = t/n;
    double gamma = tau/h;

    H.a[0] = 1;
    H.b[0] = gamma * Vpred[1] / 2;
    H.rhs[0] = Hpred[0] - Hpred[0] * Vpred[1] * gamma/2
               + (Hpred[2] * Vpred[2] - 2 * Vpred[1] * Hpred[1]
                  - (Hpred[3] * Vpred[3] - 2 * Hpred[2] * Vpred[2] + Hpred[1] * Vpred[1]) / 2
                  + Hpred[0] * (Vpred[2] - 2 * Vpred[1] - (Vpred[1] - 2 * Vpred[2] + Vpred[3])/2)) * gamma/2
               + tau * f1(0, tau*j);

    for (i = 1; i < m; ++i)
    {
         H.a[i] = 1;
        H.b[i] = (Vpred[i] + Vpred[i+1]) * gamma/4;
        H.c[i-1] = (-1) * (Vpred[i] + Vpred[i-1]) * gamma/4;
        H.rhs[i] = Hpred[i] * (1 - (Vpred[i+1] - Vpred[i-1]) * gamma/4) + tau * f1(h*i, tau*j);
    }

    H.a[m] = 1;
    H.c[m-1] = (-1) * Vpred[m-1] * gamma / 2;
    H.rhs[m] = Hpred[m] + Hpred[m] * Vpred[m-1]  * gamma/2
               - (Hpred[m-2] * Vpred[m-2] - 2 * Hpred[m-1] * Vpred[m-1]
                  - (Hpred[m-3] * Vpred[m-3] - 2 * Hpred[m-2] * Vpred[m-2] + Hpred[m-1] * Vpred[m-1])/2
                  + Hpred[m] * (Vpred[m-2] - 2 * Vpred[m-1] - (Vpred[m-3] - 2*Vpred[m-2] + Vpred[m-1])/2)) * gamma/2
               + tau * f1(x, tau*j);

}
void fill_matrix_V (matrix &V, std::vector<double> &Hpred,
                    std::vector<double> &Vpred, std::vector<double> &Hnow,
                    double mu, int n, int m, double x, double t, int j)

{
    int i = 0;
    double h = x/m;
    double tau = t/n;
    double gamma = tau/h;

    for (i = 1; i < m; ++i)
    {
        V.a[i-1] = Hnow[i] + 2 * gamma * mu / h;
        if (i < m-1)
            V.b[i-1] = (Hnow[i + 1] * Vpred[i + 1] + Hnow[i] * Vpred[i]) * gamma / 3 - gamma * mu / h;
        if (i > 1)
            V.c[i-2] = (-1) * (Hnow[i] * Vpred[i] + Hnow[i-1] * Vpred[i-1]) * gamma / 3 - gamma * mu / h;
        V.rhs[i-1] = Hpred[i] * Vpred[i] - Vpred[i] * Vpred[i] * (Hnow[i+1] - Hnow[i-1]) * gamma / 6
                     - (p(Hnow[i+1]) - p(Hnow[i-1])) * gamma / 2 + tau * f2(h*i,j*tau,mu);
    }
    return;
}