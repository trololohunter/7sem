//
// Created by vover on 11/22/17.
//

#include "functions.h"
#include <math.h>

double ro (double x, double t)
{
//    return 1;
    return (cos(M_PI * x / 10) + 1.5) * exp(t);
}

double ro_der_x (double x, double t)
{
    return (-1) * M_PI  * sin (M_PI * x /10) * exp(t) / 10;
}

double ro_der_t (double x, double t)
{
    return ro(x,t);
}

double u (double x, double t)
{
//    return 0;
    return sin(M_PI * x * x / 100) * cos(2 * M_PI * t);
}

double u_der_x (double x, double t)
{
    return M_PI * x * cos(2 * M_PI * t) * cos(M_PI * x * x/100) / 50 ;
}

double u_der_t (double x, double t)
{
    return (-1) * 2 * M_PI * sin(M_PI * t * 2) * sin (M_PI * x * x / 100);
}

double u_der_xx (double x, double t)
{
    return cos(2 * M_PI * t) * (M_PI * cos(M_PI * x * x / 100) / 50 - M_PI * M_PI * x * x * sin(M_PI * x * x / 100) / 2500);
}

double f (double x, double t)
{
    return 0;
}
double p (double x)
{
//    return x;
    return pow(x, 1.4);
}

double p_derivative (double x, double t)
{
//    return 0;
//    return ro_der_x(x,t);

    return 1.4 * pow(ro(x,t), 0.4) * ro_der_x(x,t);
}

double f1 (double x, double t)
{
//    return 0;
    return ro_der_t(x, t) + ro(x,t) * u_der_x(x,t) + u(x,t) * ro_der_x(x,t);
}
/*
double f2 (double x, double t, double mu)
{
    return 0;
    return ro(x,t) * u_der_t(x,t) + ro(x,t) * u(x,t) * u_der_x(x,t) + p_derivative(x,t) - mu * u_der_xx(x,t);
}
*/
double f2 (double x, double t, double mu)
{
//    return 0;
    return ro(x,t) * u_der_t(x,t) + u(x,t) * ro_der_t(x,t)
           + u(x,t) * u(x,t) * ro_der_x(x,t) + 2 * ro(x,t) * u(x,t) * u_der_x(x,t)
           + p_derivative(x,t) - mu * u_der_xx(x,t);
}
