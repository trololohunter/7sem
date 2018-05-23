//
// Created by vover on 11/22/17.
//

#include "residuals.h"
#include <cstdio>
#include <cmath>
#include "functions.h"

double residual_max (std::vector<double> V, std::vector<double> H, double nu, double X, double T, int N, int M, int j)
{
    double Vmax = 0, Hmax = 0;
    int i;
    double h = X/M;
    double tau = T/N;

    for (i = 0; i < M+1; ++i)
    {
        Hmax = (fabs(ro(i*h,j*tau) - H[i]) > Hmax) ? fabs(ro(i*h,j*tau) - H[i]) : Hmax;
        Vmax = (fabs(u(i*h,j*tau) - V[i]) > Vmax)  ? fabs(u(i*h,j*tau) - V[i])  : Vmax;

//        printf ("i == %d \t Vres == %e \t Hres == %e \n", i, fabs(u(i*h,j*tau) - Vpred[i]), fabs(ro(i*h,j*tau) - Hpred[i]));
//        printf ("i == %d \t Vf == %e \t Vmy == %e \n", i, u(i*h,j*tau), Vpred[i]);
//        printf ("i == %d \t Hf == %e \t Hmy== %e \n", i, ro(i*h,j*tau), Hpred[i]);
    }

    printf ("Vmax == %e \t Hmax == %e \n", Vmax, Hmax);
    return ((Hmax > Vmax) ? Hmax : Vmax);
}

double residual_L2h (std::vector<double> V, std::vector<double> H, double nu, double X, double T, int N, int M, int j)
{
    double Hres = 0;
    double Vres = 0;
    double h = X/M;
    double tau = T/N;
    int i = 0;

    Hres += 1/2 * h * (ro(0, j*tau) - H[0]) * (ro(0, j*tau) - H[0]);
    for (i = 1; i < M; ++i)
        Hres += h * (ro(i*h, j*tau) - H[i]) * (ro(i*h, j*tau) - H[i]);
    Hres += 1/2 * h * (ro(X, j*tau) - H[M]) * (ro(X, j*tau) - H[M]);

    //Vres += 1/2 * h * (u(0, j*tau) - V[i]) * (u(0, j*tau) - V[i]);
    for (i = 1; i < M; ++i)
        Vres += h * (u(i*h, j*tau) - V[i]) * (u(i*h, j*tau) - V[i]);
    //Vres += 1/2 * h * (u(X, j*tau) - V[i]) * (u(X, j*tau) - V[i]);

    printf ("Vmax == %e \t Hmax == %e \n", sqrt(Vres), sqrt(Hres));
    return ((Hres > Vres) ? sqrt(Hres) : sqrt(Vres));
}

double residual_12 (std::vector<double> V, std::vector<double> H, double nu, double X, double T, int N, int M, int j)
{
    double res = 0;

    double Hres = 0;
    double Vres = 0;
    double h = X/M;
    double tau = T/N;
    int i = 0;

    Hres += 1/2 * h * (ro(0, j*tau) - H[i]) * (ro(0, j*tau) - H[i]);
    for (i = 1; i < M; ++i)
        Hres += h * (ro(i*h, j*tau) - H[i]) * (ro(i*h, j*tau) - H[i]);
    Hres += 1/2 * h * (ro(X, j*tau) - H[i]) * (ro(X, j*tau) - H[i]);
    Hres += h * (ro(0, j*tau) - H[0]) * (ro(0, j*tau) - H[0]);
    for (i = 1; i < M; ++i)
        Hres += h * (ro(i*h, j*tau) - H[i]) * (ro(i*h, j*tau) - H[i]);

    Vres += 1/2 * h * (u(0, j*tau) - V[i]) * (u(0, j*tau) - V[i]);
    for (i = 1; i < M; ++i)
        Vres += h * (u(i*h, j*tau) - V[i]) * (u(i*h, j*tau) - V[i]);
    Vres += 1/2 * h * (u(X, j*tau) - V[i]) * (u(X, j*tau) - V[i]);
    Vres += h * (u(0, j*tau) - V[0]) * (u(0, j*tau) - V[0]);
    for (i = 1; i < M; ++i)
        Vres += h * (u(i*h, j*tau) - V[i]) * (u(i*h, j*tau) - V[i]);

    printf ("Vmax == %e \t Hmax == %e \n", sqrt(Vres), sqrt(Hres));
    return ((Hres > Vres) ? sqrt(Hres) : sqrt(Vres));

}