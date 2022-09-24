#include "matrix.h"
#include "integrals.h"
#include <stdio.h> 


int n = 4, N = n * n;


//отрезок для n мерных интегральных уравнений
double A = 0, B = 1;
double C = 0, D = 1;
double E = 0, F = 1;

//шаг для n мерных интегральных уравнений
double h1 = (B - A) / n;
double h2 = (D - C) / n;
double h3 = (F - E) / n;

inline double u(double y1, double y2, double y3) {     
    //return y1 * y1;
    //return y1 - pow(y1, 3) / 6.0;
    return y1 + y2 + y3;
}


void mk(double**& var, size_t type, size_t dim_s = 1) {
    //метод Коллокаций
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    size_t M = _msize(var) / sizeof(var[0]);
    int m = (int) sqrt(M);

    if(dim_s == 1)
        for (size_t i = 0; i < M; i++){
            double ksi = A + (i + 0.5) * h1;

            for (size_t j = 0; j < M; j++){
                double a = A + j * h1, b = A + (j + 1.0) * h1;
                
                var[i][j] = base_func(i, j) * (type - 1.0) - lambda * I_k(100, a, b, ksi);
            }
            var[i][M] = func(ksi);
        }
    else
        for (size_t i = 0; i < M; i++){
            short i1 = i / m, i2 = i % m;
            double ksi1 = A + (i1 + 0.5) * h1, ksi2 = A + (i2 + 0.5) * h2;

            for (size_t j = 0; j < M; j++) {
                short j1 = j / m, j2 = j % m;
                double a = A + j1 * h1, b = A + (j1 + 1.0) * h1, c = C + j2 * h2, d = C + (j2 + 1.0) * h2;

                var[i][j] = base_func(i, j) * (type - 1.0) - lambda * I_k(100, a, b, c, d, ksi1, ksi2);     
            }
            var[i][M] = func(ksi1, ksi2);
        }
}

void mg(double**& var, size_t type, size_t dim_s = 1) {
    //метод Галёркина
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    size_t M = _msize(var) / sizeof(var[0]);
    int m = (int)sqrt(M);

    if (dim_s == 1) {
        for (size_t i = 0; i < M; i++) {
            double a = A + i * h1, b = A + (i + 1.0) * h1;

            for (size_t j = 0; j < M; j++) {
                double c = A + j * h1, d = A + (j + 1.0) * h1;

                var[i][j] = h1 * base_func(i, j) * (type - 1.0) - lambda * In<double>(100, k, a, b, c, d);
            }
            var[i][M] = In<double>(100, a, b);
        }
    }
    else if (dim_s == 2) {
        for (size_t i = 0; i < M; i++) {
            short i1 = i / m, i2 = i % m;
            double a = A + i1 * h1, b = A + (i1 + 1.0) * h1;
            double c = C + i2 * h2, d = C + (i2 + 1.0) * h2;

            //printf("i=%d\n", i);
            for (size_t j = 0; j < M; j++) {
                short j1 = j / m, j2 = j % m;
                double e = A + j1 * h1, f = A + (j1 + 1.0) * h1;
                double g = C + j2 * h2, l = C + (j2 + 1.0) * h2;

                var[i][j] = h1 * h2 * base_func(i, j) * (type - 1.0) - lambda * In<double>(10, k, a, b, c, d, e, f, g, l);
            }
            var[i][M] = In<double>(100, func, a, b, c, d);
        }
    }
    else if (dim_s == 3) {
        int m = (int)pow(M, 1.0 / 3.0);

        for (size_t i = 0; i < M; i++) {
            short i1 = i % m, i2 = i / m - m * (i / m / m), i3 = i / m / m;

            double a = A + i1 * h1, b = a + h1;
            double c = C + i2 * h2, d = c + h2;
            double e = E + i3 * h3, f = e + h3;

            for (size_t j = 0; j < M; j++) {
                short j1 = j % m, j2 = j / m - m * (j / m / m), j3 = j / m / m;

                double g = A + j1 * h1, l = g + h1;
                double o = C + j2 * h2, p = o + h2;
                double q = E + j3 * h3, s = q + h3;

                var[i][j] = h1 * h2 * h3 * base_func(i, j) * (type - 1.0) - lambda * In<double>(1, k, a, b, c, d, e, f, g, l, o, p, q, s);
            }
            var[i][M] = In<double>(1, func, a, b, c, d, e, f);
        }
    }
}


int main1()
{
    return 0;
    double** a = createm<double>(n, n + 1.0);


    mg(a, 2);
    space(1);
    print(a);

    cout << gm(a) << "\n";


    space(1);


    for (size_t i = 0; i < n; i++)
    {
        //printf("%f %f\n", a[i][n], u(A + (i + 0.5) * h1));
    }
}

int main2() {
    return 0;
    double **a = createm<double>(N, N + 1.0), **a1 = createm<double>(N, N + 1.0);

    mg(a, 2.0, 2);
    
    mk(a1, 2.0, 2);

    //printf("mg\n");
    //print(a);
    //space(1);

    //printf("mk\n");
    //print(a1);
    //space(1);


    gm(a);


    gm(a1);

    space(1);

    printf("mg\n");
    for (size_t i = 0; i < N; i++)
    {
        if (i % n == 0 && i != 0) printf("\n");
        printf("%f ", a[i][N]);
        
    }
    space(1);

    printf("mk\n");
    for (size_t i = 0; i < N; i++)
    {
        /*short i1 = i / n, i2 = i % n;
        double ksi1 = A + (i1 + 0.5) * h1, ksi2 = A + (i2 + 0.5) * h2;

        printf("%f %f %f\n", ksi1, ksi2, a1[i][N]);*/
        
        if (i % n == 0 && i != 0) printf("\n");
        printf("%f ", a1[i][N]);
    }

}

int main() {
    N *= n;

    double** a = createm<double>(N, N + 1);

    mg(a, 2, 3);

    //print(a,"g");

    gm(a);
    double err_M = 0, err;
    for (size_t i = 0; i < N; i++)
    {
        short i1 = i % n, i2 = i / n - n * (i / n / n), i3 = i / n / n;
        err = fabs(u(A + (i1 + 0.5) * h1, C + (i2 + 0.5) * h2, E + (i3 + 0.5) * h3) - a[i][N]);
        if (err_M < err) err_M = err;
        printf("%f %f err=%f\n", u(A + (i1 + 0.5) * h1, C + (i2 + 0.5) * h2, E + (i3 + 0.5) * h3), a[i][N], err);
    }
    printf("\nMax err: %f", err_M);
    //max1=0.287156 max2=0.287156 sr1=0.085879 sr2=0.085879 n=2
    //max1=0.562298 max2=0.548107 sr1=0.099824 sr2=0.098840 n=3
    //max1=0.924069 max2=0.743103 sr1=0.356335 sr2=0.283103 n=4
    //max1=0.844089 max2=0.451763 sr1=0.106697 sr2=0.064411 n=5
    //max1=1.283065 max2=0.474908 sr1=0.376495 sr2=0.140924 n=6
    //max1=1.415884 max2=0.358216 sr1=0.387935 sr2=0.099369 n=7
    //max1=1.516391 max2=0.269375 sr1=0.396771 sr2=0.071402 n=8
    //max1=1.594911 max2=0.204893 sr1=0.402889 sr2=0.052463 n=9
    //max1=1.657861 max2=0.158189 sr1=0.408206 sr2=0.039484 n=10
}