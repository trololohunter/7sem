#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>

#include "matrix.h"
#include "residuals.h"
#include "matrix_fill.h"
#include "functions.h"

#define NN 100
#define MM 100
#define XX 10
#define TT 1
#define MU 1e-2
#define K 10000

void global_function (double mu, int n, int m, double x, double t);
double mass (std::vector<double> H, int m, double x);
double norm_L2h (std::vector<double> V, int M);

int main() {
    double mu = MU;
    double X = XX;
    double T = TT;
    int N = NN;
    int M = MM;

    while (1) {

        global_function(mu, N, M, X, T);

        if (scanf("%d", &M) != 1 || M < 1
            || scanf("%d", &N) != 1 || N < 1


                )
            break;
    }

    std::cout << "Hello, World!" << std::endl;
    return 0;
}

void global_function (double mu, int n, int m, double x, double t)
{

    double h = x/m;
    double tau = t/n;
    double mass_start;
    double mass_end;
    int i;
    int k = K;
    int j = 0;
    bool flag = true;
    double time = 0;

    matrix H(m+1);
    matrix V(m-1);

    std::vector<double> Hpred;
    std::vector<double> Hnow;
    std::vector<double> Vpred;
    std::vector<double> Vnow;

    resize_vectors (Hpred, Hnow, Vpred, Vnow, m);

//    while (k < n + 2 && k > 0) {
        FILE* ff1;
        FILE* ff2;
        FILE* d1;
        FILE* d2;
        ff1 = fopen("ff1.txt", "w");
        ff2 = fopen("ff2.txt", "w");
        d1 = fopen("dd1.txt", "w");
        d2 = fopen("dd2.txt", "w");


        fill_matrix_0(Hpred, Vpred, n, m, x, t);


//        fill_matrix_0_1(Hpred, Vpred, n, m, x, t);

//        fill_matrix_0_2(Hpred, Vpred, n, m, x, t);

//        mass_start = mass(Hpred,m,x);

//        while (flag) {

        for (j = 1; j < n + 1 && j < k; ++j) {

            fill_matrix_H(H, Hpred, Vpred, n, m, x, t, j);

            H.solver(Hnow);

            fill_matrix_V(V, Hpred, Vpred, Hnow, mu, n, m, x, t, j);

            V.solver(Vnow);

            for (i = 0; i < m+1; ++i)
            {
                Hpred[i] = Hnow[i];
            }
            for (i = 0; i < m - 1; ++i)
                Vpred[i + 1] = Vnow[i];
            Vpred[0] = 0;
            Vpred[m] = 0;
/*
            printf("step %d \t residual_max == %e \n", j, residual_max(Vpred, Hpred, mu, x, t, n, m, j));
            printf("step %d \t residual_L2h == %e \n", j, residual_L2h(Vpred, Hpred, mu, x, t, n, m, j));
            printf("step %d \t residual_12 == %e \n", j, residual_12(Vpred, Hpred, mu, x, t, n, m, j));
            printf("\n");

            int max_ind = 0;
            int min_ind = 0;
            for (i = 1; i < m+1; ++i)
                if (Hpred[i] > Hpred[max_ind]) max_ind = i;
                else
                    if (Hpred[i] < Hpred[min_ind]) min_ind = i;
            if ((norm_L2h(Vpred,m) < 1e-2 && ((Hpred[max_ind] - Hpred[min_ind]) < 1.2e-2))) flag = false;
            if ((rand() % 100) == 1) {
                printf ("time = %e \n", time);
                printf ("%e \n", (Hpred[max_ind] - Hpred[min_ind]));
                printf ("%e \n", norm_L2h(Vpred,m));
            }

     //       sleep(0.1);
            time += t/n;
*/
        }

        printf("step %d \t residual_max == %e \n", j, residual_max(Vpred, Hpred, mu, x, t, n, m, j-1));
        printf("step %d \t residual_L2h == %e \n", j, residual_L2h(Vpred, Hpred, mu, x, t, n, m, j-1));
        printf("step %d \t residual_12 == %e \n", j, residual_12(Vpred, Hpred, mu, x, t, n, m, j-1));
        printf("\n");

        for (i = 0; i < m + 1; ++i) {
            fprintf(ff1, "%lf \t %lf \n", h * i, ro(h * i, tau * (j - 1)));
            fprintf(ff2, "%lf \t %lf \n", h * i, Hpred[i]);
            fprintf(d1, "%lf \t %lf \n", h * i, u(h * i, tau * (j - 1)));
            fprintf(d2, "%lf \t %lf \n", h * i, Vpred[i]);
        }
/*
        mass_end = mass(Hpred,m,x);

        printf ("start mass = %lf \n", mass_start);
        printf ("end mass = %lf \n", mass_end);
        printf ("start mass - end mass = %lf \n", mass_start - mass_end);
        printf ("stabilization time = %e \n", time);
*/        fclose(ff1);
        fclose(ff2);
        fclose(d1);
        fclose(d2);
//        scanf("%d", &k);
//    }

    return;
}


double mass (std::vector<double> H, int m, double x)
{
    double res = 0;
    double h = x / m;

    res += h/2 * H[0];
    for(int i = 1; i < m; ++i)
        res += h * H[i];
    res += h/2 * H[m];

    return res;

}

double norm_L2h (std::vector<double> V, int M)
{
    double Vres = 0;
    int i = 0;

    for (i = 1; i < M; ++i)
        Vres += V[i] * V[i];
    return sqrt(Vres);
}