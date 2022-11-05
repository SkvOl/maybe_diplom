#pragma once
#ifndef INTEGRALS_H
#define INTEGRALS_H
#include "constants.h"
#include "diff_geom.h"



//Для комплексных нужно переопределять 
inline complex<double> k_c(double x1, ...) {
    double x2, x3, y1, y2, y3;
    va_list args;

    va_start(args, x1);
    x2 = va_arg(args, double);
    x3 = va_arg(args, double);

    y1 = va_arg(args, double);
    y2 = va_arg(args, double);
    y3 = va_arg(args, double);
    va_end(args);

    double custom_eps = h1_obj <= h2_obj ? h1_obj : h2_obj;
    //double custom_eps = (h1 + h2) / 2.0;
    //double custom_eps = sqrt(h1 * h1 + h2 * h2);
    
    complex<double> i(0, 1), r = sqrt(pow(x1 - y1, 2) + pow(x2 - y2, 2) + pow(x3 - y3, 2));
    r.real() < (custom_eps / 100.0) ? r += (custom_eps / 100.0) : r = r;
    return exp(i * k0 * r) / 4.0 / pi / r;


    //complex<double> i(0, 1);
    //double r = sqrt(pow(x1 - y1, 2) + pow(x2 - y2, 2) + pow(x3 - y3, 2));
    //return (r > 1e-3 ? (exp(i * k0 * r) / 4.0 / pi / r) : 0);
}

inline complex<double> func_c(double x1, ...) {
    va_list args;
    va_start(args, x1);
    //complex<double> x2 = va_arg(args, complex<double>);
    va_end(args);

    complex<double> i_k_x(0, x1);
    return exp(i_k_x*k0);
}

inline void func_cv(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...) {
    double t1, t2, t3;
    complex<double>* var;

    va_list args;
    va_start(args, func_x3);
    t1 = va_arg(args, double);
    t2 = va_arg(args, double);
    t3 = va_arg(args, double);
    var = va_arg(args, complex<double>*);
    va_end(args);
    
    var[0] = 0;
    var[1] = exp(complex<double>(0, (*func_x1)(0, t1, t2, t3, NULL)) * k0);
    var[2] = 0;
}

inline void func_cv_tang(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...) {
    double t1, t2, t3;
    complex<double>* var;

    va_list args;
    va_start(args, func_x3);
    t1 = va_arg(args, double);
    t2 = va_arg(args, double);
    t3 = va_arg(args, double);
    var = va_arg(args, complex<double>*);
    va_end(args);


    double* normal = createv<double>(3);
    vector_n((*func_x1), (*func_x2), (*func_x3), t1, t2, t3, normal);
    double x1 = (*func_x1)(0, t1, t2, t3, NULL), x2 = (*func_x2)(0, t1, t2, t3, NULL), x3 = (*func_x3)(0, t1, t2, t3, NULL);

    var[0] = 0;
    var[1] = exp(complex<double>(0, x1) * k0);
    var[2] = 0;

    complex<double> tem = multv3(var, normal);

    var[0] -= normal[0] * tem;
    var[1] -= normal[1] * tem;
    var[2] -= normal[2] * tem;

    del(normal);
}

inline void permittivity_obj(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...) {
    double t1, t2, t3;
    complex<double>* var;

    va_list args;
    va_start(args, func_x3);
    t1 = va_arg(args, double);
    t2 = va_arg(args, double);
    t3 = va_arg(args, double);
    var = va_arg(args, complex<double>*);
    va_end(args);

    if ((*func_x2)(0, t1, t2, t3, NULL) >= 0) {
        var[0] = 0.2;
        var[1] = 0.2;
        var[2] = 0.2;
    }
    else {
        var[0] = 0.0;
        var[1] = 0.0;
        var[2] = 0.0;
    }
}


inline complex<double> In_c(int N_i, complex<double>(*function)(double, ...), double a, double b, double c, double d, double e, double f, double g, double h) {
    //четырёхмерный интеграл
    double h_i = (b - a) / N_i,
        h_j = (d - c) / N_i,
        h_k = (f - e) / N_i,
        h_z = (h - g) / N_i;

    complex<double> Sum=0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            for (size_t k = 0; k < N_i; k++) {
                for (size_t z = 0; z < N_i; z++) {

                    double l1 = a + (i + 0.5) * h_i,
                        l2 = c + (j + 0.5) * h_j,
                        l3 = e + (k + 0.5) * h_k,
                        l4 = g + (z + 0.5) * h_z;

                    if (l1 == l3 && l2 == l4) Sum += (*function)(l1, l2, l3 + 10.0 * mach_eps, l4 + 10.0 * mach_eps);
                    else Sum += (*function)(l1, l2, l3, l4);
                }
            }
        }
    }
    complex<double> h_i_c(h_i, 0), h_j_c(h_j, 0), h_k_c(h_k, 0), h_z_c(h_z, 0);
    return h_i_c * h_j_c * h_k_c * h_z_c * Sum;
}
inline complex<double> In_c(int N_i, complex<double>(*function)(double, ...), double a, double b, double c, double d) {
    //двумерный интеграл
    double h_i = (b - a) / N_i,
        h_j = (d - c) / N_i;

    complex<double> Sum(0, 0);
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            double l1 = a + (i + 0.5) * h_i,
                l2 = c + (j + 0.5) * h_j;
            Sum += (*function)(l1, l2);
        }
    }

    complex<double> h_i_c(h_i, 0), h_j_c(h_j, 0);
    return h_i_c * h_j_c * Sum;
}
//Для комплексных нужно переопределять 


template<typename type_k>
inline type_k k(type_k x1, ...) {
    type_k x2, x3, y1, y2, y3;
    va_list args;

    va_start(args, x1);
    x2 = va_arg(args, type_k);
    x3 = va_arg(args, type_k);

    y1 = va_arg(args, type_k);
    y2 = va_arg(args, type_k);
    y3 = va_arg(args, type_k);
    va_end(args);


    return x1 * x2 * x3 - y1 * y2 * y3;

    //return x1 + x2 - y1 * y2;

    //return x1 - y1;
}

template<typename type_func>
inline type_func func(type_func x1, ...) {
    va_list args;
    va_start(args, x1);
    type_func x2 = va_arg(args, type_func);
    type_func x3 = va_arg(args, type_func);
    va_end(args);

    //return sin(10.0 * x1 + 10.0 * x2) + (1.0 / 10000.0) * (-200.0 * x1 - 200.0 * x2 - 2.0) * sin(10.0) + (1.0 / 10000.0) * (100.0 * x1 + 100.0 * x2 - 99.0) * sin(20.0) + (1.0 / 500.0) * cos(10.0) - (1.0 / 500.0) * cos(20.0);
    //return (1.0 / 6.0) * (-6.0 * x1 - 6.0 * x2 + 3.0) * cos(1.0) + (1.0 / 2.0) * x1 + (3.0 / 2.0) * x2 - (1.0 / 2.0) * sin(1.0) - sin(x1) + 1.0 / 6.0;
    //return 1.0 / 3.0;
    //return pow(x1, 2) - 1.0 * x1 / 3.0 + 1.0 / 4.0;

    return -(x1 + x2 + x3) * (x1 * x2 * x3 - 9.0 / 8.0);
}

inline double base_func(int i, int j) {
    return i == j ? 1.0 : 0.0;
}

//inline void base_func(double* var,...) {
//    //size_t M = _msize(var) / sizeof(var[0]);
//    double it1, it2, it3, jt1, jt2, jt3;   
//    
//    va_list args;
//    va_start(args, var);
//    it1 = va_arg(args, double);
//    it2 = va_arg(args, double);
//    it3 = va_arg(args, double);
//
//    jt1 = va_arg(args, double);
//    jt2 = va_arg(args, double);
//    jt3 = va_arg(args, double);
//    va_end(args);
//
//    if (var[0] == 1) {
//        var[1] = 0.0;
//        var[2] = 0.0;
//        if ((jt1 >= it1 - h1_obj && jt1 <= it1) && (it2 == jt2)) {
//            var[0] = 1.0 + (jt1 - it1) / h1_obj;
//        }
//        else if ((jt1 <= it1 + h1_obj && jt1 >= it1) && it2 == jt2) {
//            var[0] = 1.0 - (jt1 - it1) / h1_obj;
//        }
//        else var[0] = 0.0;
//    }
//    else if (var[1] == 1) {
//        var[0] = 0.0;
//        var[2] = 0.0;
//        if ((it1 == jt1) && (jt2 >= it2 - h2_obj && jt2 <= it2)) {
//            var[1] = 1.0 + (jt2 - it2) / h2_obj;
//        }
//        else if ((it1 == jt1) && (jt2 <= it2 + h2_obj && jt2 >= it2)) {
//            var[1] = 1.0 - (jt2 - it2) / h2_obj;
//        }
//        else var[1] = 0.0;
//    }
//}

inline void base_func_rft(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...) {
    short down_index1, down_index2, up_index1, up_index2;
    double* var, t1, t2, t3;
    bool object;

    va_list args;
    va_start(args, func_x3);
    t1 = va_arg(args, double);
    t2 = va_arg(args, double);
    t3 = va_arg(args, double);
    down_index1 = va_arg(args, short);
    down_index2 = va_arg(args, short);
    up_index1 = va_arg(args, short);
    up_index2 = va_arg(args, short);
    object = va_arg(args, bool);
    var = va_arg(args, double*);
    va_end(args);
   
    
    double** J = createm<double>(3, 2), * v_vec = createv<double>(2);
    double A_base_func = object ? A_obj : A_screen, C_base_func = object ? C_obj : C_screen, h1_base_func = object ? h1_obj : h1_screen, h2_base_func = object ? h2_obj : h2_screen;
    

    Jacobi((*func_x1), (*func_x2), (*func_x3), t1, t2, 0, J);

    double d1, d2, d3, d4;
    if (up_index2 == 1) {
        d1 = A_base_func + h1_base_func * down_index1, d2 = d1 + 2.0 * h1_base_func, d3 = C_base_func + h2_base_func * down_index2, d4 = d3 + h2_base_func;
        v_vec[0] = d1 <= t1 && t1 <= d2 && d3 <= t2 && t2 <= d4 ? 1.0 - fabs(t1 - A_base_func - h1_base_func * (down_index1 + 1.0)) / h1_base_func / 1.0 : 0.0;
        v_vec[1] = 0.0;
    }
    else {
        d1 = A_base_func + h1_base_func * down_index1, d2 = d1 + h1_base_func, d3 = C_base_func + h2_base_func * down_index2, d4 = d3 + 2.0 * h2_base_func;
        v_vec[0] = 0.0;
        v_vec[1] = d1 <= t1 && t1 <= d2 && d3 <= t2 && t2 <= d4 ? 1.0 - fabs(t2 - C_base_func - h2_base_func * (down_index2 + 1.0)) / h2_base_func / 1.0 : 0.0;
    }


    mult(J, v_vec, var);

    del(J);
    del(v_vec);
}

inline void base_func(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...) {
    short down_index1, down_index2, down_index3;
    double* var, t1, t2, t3;
    bool object;

    va_list args;
    va_start(args, func_x3);
    t1 = va_arg(args, double);
    t2 = va_arg(args, double);
    t3 = va_arg(args, double);
    down_index1 = va_arg(args, short);
    down_index2 = va_arg(args, short);
    down_index3 = va_arg(args, short); 
    object = va_arg(args, bool);
    var = va_arg(args, double*);
    va_end(args);
    
    double** J = createm<double>(3, 2), * v_vec = createv<double>(2);
    double A_base_func = object ? A_obj : A_screen, 
           C_base_func = object ? C_obj : C_screen,
           E_base_func = object ? E_obj : E_screen,
           h1_base_func = object ? h1_obj : h1_screen, 
           h2_base_func = object ? h2_obj : h2_screen,
           h3_base_func = object ? h3_obj : h3_screen;

    Jacobi((*func_x1), (*func_x2), (*func_x3), t1, t2, t3, J);

    double d1, d2, d3, d4, d5, d6;
    
    d1 = A_base_func + h1_base_func * down_index1; d2 = d1 + 2.0 * h1_base_func; 
    d3 = C_base_func + h2_base_func * down_index2; d4 = d3 + h2_base_func;
    d5 = C_base_func + h2_base_func * down_index2; d6 = d3 + h2_base_func;

    v_vec[0] = 1.0;
    v_vec[1] = 1.0;
    v_vec[3] = 1.0;
   

    mult(J, v_vec, var);

    del(J);
    del(v_vec);
}


inline void grad(double(*function)(double, ...), double* res, double x1, double x2) {
    res[0] = ((*function)(x1 * (1 + mach_eps), x2) - (*function)(x1 * (1 - mach_eps), x2)) / 2.0 / x1 / mach_eps;
    res[1] = ((*function)(x1, x2 * (1 + mach_eps)) - (*function)(x1, x2 * (1 - mach_eps))) / 2.0 / x2 / mach_eps;
}

inline double div(double(*function)(double, ...), double x1, double x2) {
    double* res = createv<double>(2);

    grad((*function), res, x1, x2);
    return res[0] + res[1];
}



template<typename type_I_k1>
inline type_I_k1 I_k(int N_i, type_I_k1 a, type_I_k1 b, type_I_k1 ksi) {
    //одномерный интеграл для метода Коллокаций 
    type_I_k1 h_int = (b - a) / (N_i * 1.0);
    type_I_k1 Sum = 0.0;
    for (size_t i = 0; i < N_i; i++) {
        type_I_k1 l = (i + 0.5) * h_int + a;
        Sum += k(ksi, l) * h_int;
    }
    return Sum;
}

template<typename type_I_k2>
inline type_I_k2 I_k(int N_i, type_I_k2 a, type_I_k2 b, type_I_k2 c, type_I_k2 d, type_I_k2 ksi1, type_I_k2 ksi2) {
    //двумерный интеграл для метода Коллокаций
    type_I_k2 h_i = (b - a) / N_i, h_j = (d - c) / N_i;
    type_I_k2 Sum = 0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            type_I_k2 l1 = a + (i + 0.5) * h_i, l2 = c + (j + 0.5) * h_j;
            Sum += k(ksi1, ksi2, l1, l2);
        }
    }
    return h_i * h_j * Sum;
}

template<typename type_I1>
inline type_I1 In(int N_i, type_I1 a, type_I1 b) {
    //одномерный интеграл
    type_I1 h_int = (b - a) / (N_i * 1.0);
    type_I1 Sum = 0.0;
    for (size_t i = 0; i < N_i; i++) {
        type_I1 l = (i + 0.5) * h_int + a;
        Sum += func(l);
    }
    return Sum * h_int;
}

template<typename type_I2>
inline type_I2 In(int N_i, type_I2(*function)(type_I2, ...), type_I2 a, type_I2 b, type_I2 c, type_I2 d) {
    //двумерный интеграл
    type_I2 h_i = (b - a) / N_i,
        h_j = (d - c) / N_i;
    type_I2 Sum = 0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            type_I2 l1 = a + (i + 0.5) * h_i,
                l2 = c + (j + 0.5) * h_j;
            Sum += (*function)(l1, l2);
        }
    }
    return h_i * h_j * Sum;
}

template<typename type_I3>
inline type_I3 In(int N_i, type_I3(*function)(type_I3, ...), type_I3 a, type_I3 b, type_I3 c, type_I3 d, type_I3 e, type_I3 f) {
    //трёхмерный интеграл
    type_I3 h_i = (b - a) / N_i,
        h_j = (d - c) / N_i,
        h_k = (f - e) / N_i;

    type_I3 Sum = 0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            for (size_t k = 0; k < N_i; k++) {
                type_I3 l1 = a + (i + 0.5) * h_i,
                    l2 = c + (j + 0.5) * h_j,
                    l3 = e + (k + 0.5) * h_k;

                Sum += (*function)(l1, l2, l3);
            }
        }
    }
    return h_i * h_j * h_k * Sum;
}

template<typename type_I4>
inline type_I4 In(int N_i, type_I4(*function)(type_I4, ...), type_I4 a, type_I4 b, type_I4 c, type_I4 d, type_I4 e, type_I4 f, type_I4 g, type_I4 h) {
    //четырёхмерный интеграл
    type_I4 h_i = (b - a) / N_i,
        h_j = (d - c) / N_i,
        h_k = (f - e) / N_i,
        h_z = (h - g) / N_i;

    type_I4 Sum = 0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            for (size_t k = 0; k < N_i; k++) {
                for (size_t z = 0; z < N_i; z++) {

                    type_I4 l1 = a + (i + 0.5) * h_i,
                        l2 = c + (j + 0.5) * h_j,
                        l3 = e + (k + 0.5) * h_k,
                        l4 = g + (z + 0.5) * h_z;
                    Sum += (*function)(l1, l2, l3, l4);
                }
            }
        }
    }
    return h_i * h_j * h_k * h_z * Sum;
}

template<typename type_I6>
inline type_I6 In(int N_i, type_I6(*function)(type_I6, ...), type_I6 a, type_I6 b, type_I6 c, type_I6 d, type_I6 e, type_I6 f, type_I6 g, type_I6 h, type_I6 l, type_I6 o, type_I6 p, type_I6 q) {
    //шестимерый интеграл
    type_I6 h_i = (b - a) / N_i,
        h_j = (d - c) / N_i,
        h_k = (f - e) / N_i,
        h_z = (h - g) / N_i,
        h_v = (o - l) / N_i,
        h_w = (q - p) / N_i;

    type_I6 Sum = 0;
    for (size_t i = 0; i < N_i; i++) {
        for (size_t j = 0; j < N_i; j++) {
            for (size_t k = 0; k < N_i; k++) {
                for (size_t z = 0; z < N_i; z++) {
                    for (size_t v = 0; v < N_i; v++) {
                        for (size_t w = 0; w < N_i; w++) {

                            type_I6 l1 = a + (i + 0.5) * h_i,
                                l2 = c + (j + 0.5) * h_j,
                                l3 = e + (k + 0.5) * h_k,
                                l4 = g + (z + 0.5) * h_z,
                                l5 = l + (v + 0.5) * h_v,
                                l6 = p + (w + 0.5) * h_w;
                            Sum += (*function)(l1, l2, l3, l4, l5, l6);
                        }
                    }
                }
            }
        }
    }
    return h_i * h_j * h_k * h_z * h_v * h_w * Sum;
}

template<typename type_Is1>
inline type_Is1 Is(int N_i, type_Is1(*function)(type_Is1, ...), type_Is1 a, type_Is1 b) {
    //одномерный интеграл методом Симпсона
    type_Is1 h_i = (b - a) / N_i;
    type_Is1 Sum = 0.0;
    for (size_t i = 0; i < N_i; i++) {
        type_Is1 x1 = a + i * h_i, x2 = x1 + h_i;

        Sum += (x2 - x1) / 6.0 * ((*function)(x1) + 4.0 * (*function)(0.5 * (x1 + x2)) + (*function)(x2));
    }
    return Sum;
}

template<typename type_Is2>
inline type_Is2 Is(int N_i, type_Is2(*function)(type_Is2, ...), type_Is2 a, type_Is2 b, type_Is2 c, type_Is2 d) {
    //двумерный интеграл методом Симпсона
    type_Is2 h_i = (b - a) / N_i,
        h_j = (d - c) / N_i,
        Sum = 0.0;
    for (size_t i = 0; i < N_i; i++) {
        type_Is2 xi1 = a + i * h_i,
            xi2 = xi1 + h_i;
        for (size_t j = 0; j < N_i; j++) {
            type_Is2 xj1 = a + j * h_j,
                xj2 = xj1 + h_j;

            Sum += (xi2 - xi1) * (xj2 - xj1) / 6.0 * ((*function)(xi1, xj1) + 4.0 * (*function)(0.5 * (xi1 + xi2), 0.5 * (xj1 + xj2)) + (*function)(xi2, xj2));
        }
    }
    return Sum;
}

template<typename type_Is3>
inline type_Is3 Is(int N_i, type_Is3(*function)(type_Is3, ...), type_Is3 a, type_Is3 b, type_Is3 c, type_Is3 d, type_Is3 e, type_Is3 f) {
    //трёхмерный интеграл методом Симпсона
    type_Is3 h_i = (b - a) / N_i,
        h_j = (d - c) / N_i,
        h_k = (f - e) / N_i,
        Sum = 0.0;
    for (size_t i = 0; i < N_i; i++) {
        type_Is3 xi1 = a + i * h_i,
            xi2 = xi1 + h_i;

        for (size_t j = 0; j < N_i; j++) {
            type_Is3 xj1 = a + j * h_j,
                xj2 = xj1 + h_j;

            for (size_t k = 0; k < N_i; k++) {
                type_Is3 xk1 = a + k * h_k,
                    xk2 = xk1 + h_k;

                Sum += (xi2 - xi1) * (xj2 - xj1) * (xk2 - xk1) / 6.0 * ((*function)(xi1, xj1, xk1) + 4.0 * (*function)(0.5 * (xi1 + xi2), 0.5 * (xj1 + xj2), 0.5 * (xk1 + xk2)) + (*function)(xi2, xj2, xk2));
            }
        }
    }
    return Sum;
}

template<typename type_Is4>
inline type_Is4 Is(int N_i, type_Is4(*function)(type_Is4, ...), type_Is4 a, type_Is4 b, type_Is4 c, type_Is4 d, type_Is4 e, type_Is4 f, type_Is4 g, type_Is4 h) {
    //четырёхмерный интеграл методом Симпсона
    type_Is4 h_i = (b - a) / N_i,
        h_j = (d - c) / N_i,
        h_k = (f - e) / N_i,
        h_z = (h - g) / N_i,
        Sum = 0.0;
    for (size_t i = 0; i < N_i; i++) {
        type_Is4 xi1 = a + i * h_i,
            xi2 = xi1 + h_i;

        for (size_t j = 0; j < N_i; j++) {
            type_Is4 xj1 = a + j * h_j,
                xj2 = xj1 + h_j;

            for (size_t k = 0; k < N_i; k++) {
                type_Is4 xk1 = a + k * h_k,
                    xk2 = xk1 + h_k;

                for (size_t z = 0; z < N_i; z++) {
                    type_Is4 xz1 = a + z * h_z,
                        xz2 = xz1 + h_z;

                    Sum += (xi2 - xi1) * (xj2 - xj1) * (xk2 - xk1) * (xz2 - xz1) / 6.0 * ((*function)(xi1, xj1, xk1, xz1) + 4.0 * (*function)(0.5 * (xi1 + xi2), 0.5 * (xj1 + xj2), 0.5 * (xk1 + xk2), 0.5 * (xz1 + xz2)) + (*function)(xi2, xj2, xk2, xz2));
                }
            }
        }
    }
    return Sum;
}

template<typename type_Is6>
inline type_Is6 Is(int N_i, type_Is6(*function)(type_Is6, ...), type_Is6 a, type_Is6 b, type_Is6 c, type_Is6 d, type_Is6 e, type_Is6 f, type_Is6 g, type_Is6 h, type_Is6 l, type_Is6 o, type_Is6 p, type_Is6 q) {
    //шестимерный интеграл методом Симпсона
    type_Is6 h_i = (b - a) / N_i,
        h_j = (d - c) / N_i,
        h_k = (f - e) / N_i,
        h_z = (h - g) / N_i,
        h_v = (o - l) / N_i,
        h_w = (q - p) / N_i,
        Sum = 0.0;
    for (size_t i = 0; i < N_i; i++) {
        type_Is6 xi1 = a + i * h_i,
            xi2 = xi1 + h_i;

        for (size_t j = 0; j < N_i; j++) {
            type_Is6 xj1 = a + j * h_j,
                xj2 = xj1 + h_j;

            for (size_t k = 0; k < N_i; k++) {
                type_Is6 xk1 = a + k * h_k,
                    xk2 = xk1 + h_k;

                for (size_t z = 0; z < N_i; z++) {
                    type_Is6 xz1 = a + z * h_z,
                        xz2 = xz1 + h_z;

                    for (size_t v = 0; v < N_i; v++) {
                        type_Is6 xv1 = a + v * h_v,
                            xv2 = xv1 + h_v;

                        for (size_t w = 0; w < N_i; w++) {
                            type_Is6 xw1 = a + w * h_w,
                                xw2 = xw1 + h_w;

                            Sum += (xi2 - xi1) * (xj2 - xj1) * (xk2 - xk1) * (xz2 - xz1) * (xv2 - xv1) * (xw2 - xw1) / 6.0 * ((*function)(xi1, xj1, xk1, xz1, xv1, xw1) + 4.0 * (*function)(0.5 * (xi1 + xi2), 0.5 * (xj1 + xj2), 0.5 * (xk1 + xk2), 0.5 * (xz1 + xz2), 0.5 * (xv1 + xv2), 0.5 * (xw1 + xw2)) + (*function)(xi2, xj2, xk2, xz2, xv2, xw2));
                        }
                    }
                }
            }
        }
    }
    return Sum;
}

#endif INTEGRALS_H
