#include "matrix.h"

const double pi = 3.1415926, Eps = 0.0001, k0 = 1, k1 = 1.5 * k0;
const int n = 10, N = n * n;
double R = 5;
double A = 0, B = 2 * pi, h1 = (B - A) / n,
C = 0, D = pi, h2 = (D - C) / n;


inline double u(double y) {     //точное решение

    //return y * y;
    return 0;//C * 1.0 * (A * A - y * y);
}

inline double k(double x, double y) {
    return x - y;
    //return 1.0 / (H * H + (x - y) * (x - y));
}

inline double f(double x) {
    return pow(x, 2) - 1.0 * x / 3.0 + 1.0 / 4.0;
    /*double one = (x - cc) / H;
    double two = (x + cc) / H;
    double three = sqrt(H * H + pow(x - cc, 2.0));
    double four = sqrt(H * H + pow(x + cc, 2.0));
    return one * atan(one) +
        two * atan(two) -
        2.0 * x / H * atan(x / H) +
        log((x * x + H * H) / (three * four));*/
}

inline double I(int N_i, double a, double b, double c, double d) {
    double h_i = (b - a) / N_i, h_j = (d - c) / N_i;
    double Sum = 0;
    for (int i = 0; i < N_i; i++) {
        for (int j = 0; j < N_i; j++) {
            Sum += k(a + i * h_i + h_i / 2.0, c + j * h_j + h_j / 2.0);
        }
    }
    return h_i * h_j * Sum;
}

inline double I(int N_i, double a, double b) {
    double h_int = (b - a) / N_i;
    double Sum = 0;
    for (int i = 0; i < N_i; i++) {
        Sum += f(a + i * h_int + h_int / 2.0);
    }
    return Sum * h_int;
}

void mg2(double**& var) {
    for (size_t i = 0; i < n; i++) {

        double a = A + i * h, b = A + (i + 1.0) * h;

        for (size_t j = 0; j < n; j++) {

            double c = A + j * h, d = A + (j + 1.0) * h;

            if (i == j) var[i][j] = h - I(100, a, b, c, d);
            if (i != j) var[i][j] = 0 - I(100, a, b, c, d);
        }
        var[i][n] = I(100, a, b);
    }
}

void mg1(double**& var) {
    for (size_t i = 0; i < n; i++) {

        double a = A + i * h1, b = A + (i + 1.0) * h1;

        for (size_t j = 0; j < n; j++)
        {
            double c = A + j * h, d = A + (j + 1.0) * h;

            var[i][j] = I(100, a, b, c, d);
        }
        var[i][n] = I(100, a, b);
    }
}

int main()
{
    double** a = createm(3,3, true), *a1 = createv(3,true);
    cout << "a=" << size(a) << " a1=" << size(a1) << "\n";
    print(a);
    a[1][0] = 1; a[0][1] = 1;
    print(a);
    cout << "a=" << size(a) << " a1=" << size(a1) << "\n";
    free(a);
    cout << size(a) << "\n";
}
