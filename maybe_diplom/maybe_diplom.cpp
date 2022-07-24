#include "matrix.h"
#include "integrals.h"

const double pi = 3.1415926, Eps = 0.0001, k0 = 1, k1 = 1.5 * k0;
const int n = 100, N = n * n;
double R = 5;
double A = 0, B = 1, h = (B - A) / n,
C = 0, D = pi, h2 = (D - C) / n;
double lambda = 1;


inline double u(double y) {     
    //точное решение

    return y * y;
    //return 0;//C * 1.0 * (A * A - y * y);
}


void mk(double**& var, size_t type, size_t dim_s = 1) {
    //метод Коллокаций
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    for (size_t i = 0; i < n; i++)
    {
        double ksi = A + (i + 0.5) * h;

        for (size_t j = 0; j < n; j++){
            
            double a = A + j * h, b = A + (j + 1.0) * h;
            var[i][j] = base_func(i, j) * (type - 1.0) - lambda * I_k(100, a, b, ksi);
        }
        var[i][n] = f(ksi);
    }
}

void mg(double**& var, size_t type, size_t dim_s = 1) {
    //метод Галёркина
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    for (size_t i = 0; i < n; i++) {

        double a = A + i * h, b = A + (i + 1.0) * h;

        for (size_t j = 0; j < n; j++) {

            double c = A + j * h, d = A + (j + 1.0) * h;
            var[i][j] = h * base_func(i, j) * (type - 1.0) - lambda * I(100, a, b, c, d);
        }
        var[i][n] = I(100, a, b);
    }
}


int main()
{
    double** a = createm(n, n + 1.0);

    mg(a, 2);

    cout << gm(a) << "\n";

    print(a);

    space(2);

    for (size_t i = 0; i < n; i++)
    {
        printf("%f\n", a[i][n]);
    }
}


