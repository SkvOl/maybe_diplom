#pragma once
#ifndef INTEGRALS_H
#define INTEGRALS_H
#include <math.h>
#include <stdarg.h>

inline double k(double y1, ...) {
    double y2, x1, x2;
    va_list args;

    va_start(args, y1);
    y2 = va_arg(args, double);
    x1 = va_arg(args, double);
    x2 = va_arg(args, double);
    va_end(args);

    //return x1 + x2 - y1 * y2;
    return x1 - x2;
}

inline double f(double x1, ...) {
    va_list args;

    va_start(args, x1);
    double x2 = va_arg(args, double);
    va_end(args);

    //return 1.0 / 3.0;
    return pow(x1, 2) - 1.0 * x1 / 3.0 + 1.0 / 4.0;
}

inline double base_func(int i, int j) {
    return i == j ? 1 : 0;
}



inline double I_k(int N_i, double a, double b, double ksi) {
    //���������� �������� ��� ������ ���������� 
    double h_int = (b - a) / (N_i * 1.0);
    double Sum = 0.0;
    for (size_t i = 0; i < N_i; i++){
        double l = (i + 0.5) * h_int + a;
        Sum += k(ksi, l) * h_int;
    }
    return Sum;
}

inline double I_k(int N_i, double a, double b, double c, double d, double ksi1, double ksi2) {
    //��������� �������� ��� ������ ����������
    double h_i = (b - a) / N_i, h_j = (d - c) / N_i;
    double Sum = 0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            double l1 = a + (i + 0.5) * h_i, l2 = c + (j + 0.5) * h_j;
            Sum += k(l1, l2, ksi1, ksi2);
        }
    }
    return h_i * h_j * Sum;
}

inline double I(int N_i, double(*function)(double, ...), double a, double b, double c, double d) {
    //��������� ��������
    double h_i = (b - a) / N_i, h_j = (d - c) / N_i;
    double Sum = 0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            double l1 = a + (i + 0.5) * h_i, l2 = c + (j + 0.5) * h_j;
            Sum += (*function)(l1, l2);
        }
    }
    return h_i * h_j * Sum;
}

inline double I(int N_i, double a, double b) {
    //���������� ��������
    double h_int = (b - a) / (N_i * 1.0);
    double Sum = 0.0;
    for (size_t i = 0; i < N_i; i++) {
        double l = (i + 0.5) * h_int + a;
        Sum += f(l);
    }
    return Sum * h_int;
}

#endif INTEGRALS_H