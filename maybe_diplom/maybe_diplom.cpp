#include "matrix.h"
#include "integrals.h"

const double pi = 3.1415926, Eps = 0.0001;
const double  k0 = 1, k1 = 1.5 * k0;

const int n = 10, N = n * n;

double R = 5;
double lambda = 1;

//отрезок для одномерных интегральных уравнений
//double A = 0, B = 1;

//отрезок для двумерных интегральных уравнений
double A = 0, B = 1;
double C = 0, D = 1;

//шаг для одномерных интегральных уравнений
double h = (B - A) / (n-1);

//отрезок для двумерных интегральных уравнений
double h1 = (B - A) / n, h2 = (D - C) / n;


inline double u(double y1) {     
    return y1 * y1;
}


void mk(double**& var, size_t type, size_t dim_s = 1) {
    //метод Коллокаций
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    size_t M = _msize(var) / sizeof(var[0]);
    int m = (int) sqrt(M);

    if(dim_s == 1)
        for (size_t i = 0; i < M; i++)
        {
            double ksi = A + (i + 0.5) * h;

            for (size_t j = 0; j < M; j++){
            
                double a = A + j * h, b = A + (j + 1.0) * h;
                var[i][j] = base_func(i, j) * (type - 1.0) - lambda * I_k(100, a, b, ksi);
            }
            var[i][M] = f(ksi);
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
            var[i][M] = f(ksi1, ksi2);
        }
}

void mg(double**& var, size_t type, size_t dim_s = 1) {
    //метод Галёркина
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    size_t M = _msize(var) / sizeof(var[0]);

    if (dim_s == 1) {
        double(*function)(double, ...) = NULL;
        function = &k;

        for (size_t i = 0; i < M; i++) {

            double a = A + i * h, b = A + (i + 1.0) * h;

            for (size_t j = 0; j < M; j++) {

                double c = A + j * h, d = A + (j + 1.0) * h;
                var[i][j] = h * base_func(i, j) * (type - 1.0) - lambda * I(100, function, a, b, c, d);
            }
            var[i][M] = I(100, a, b);
        }
    }
    else {

    }
}


int main()
{
    double** a = createm(n, n + 1.0);

    //mk(a, 2, 2);
    mg(a, 2);

    cout << gm(a) << "\n";

    //print(a);

    space(1);

    /*for (size_t i = 0; i < N; i++)
    {
        short i1 = i / n, i2 = i % n;
        double ksi1 = A + (i1 + 0.5) * h1, ksi2 = A + (i2 + 0.5) * h2;
        
        printf("%f %f %f\n", ksi1, ksi2, a[i][N]);
    }*/

    for (size_t i = 0; i < n; i++)
    {
        printf("%f\n", a[i][n]);
    }
}


