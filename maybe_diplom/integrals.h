#pragma once
#ifndef INTEGRALS_H
#define INTEGRALS_H
#include <math.h>

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

inline double base_func(int i, int j) {
    if (i == j) return 1;
    else return 0;
}



inline double I_k(int N_i, double a, double b, double ksi) {
    //одномерный интеграл для метода Коллокаций 
    double h_int = (b - a) / double(N_i);
    double Sum = 0.0;
    for (size_t i = 0; i < N_i; i++)
    {
        double l = (i + 0.5) * h_int + a;
        Sum += k(ksi, l) * h_int;
    }
    return Sum;
}

inline double I(int N_i, double a, double b, double c, double d) {
    //двумерный интеграл
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
    //одномерный интеграл
    double h_int = (b - a) / N_i;
    double Sum = 0;
    for (int i = 0; i < N_i; i++) {
        Sum += f(a + i * h_int + h_int / 2.0);
    }
    return Sum * h_int;
}

#endif INTEGRALS_H