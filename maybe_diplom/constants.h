#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <math.h>

const int _n = 2, _N = 4 * _n * _n - 4 * _n;

const double mach_eps = sqrt(2.22045e-16);//1e-5;
const double pi = acos(-1);
const complex<double> k0(1, 0.01), k1 = 1.5 * k0;
//double k0 = 2 * pi;
const double R = 1.0;
const double lambda = 1.0;

//отрезок для n мерных интегральных уравнений
//const double A_obj = -pi / 2.0, B_obj = pi / 2.0;
//const double C_obj = 0.0, D_obj = 2.0 * pi;

const double A_obj = 0.0, B_obj = 2.0 * pi;
const double C_obj = 0.0, D_obj = pi;
const double E_obj = 0, F_obj = 1;


const double A_screen = 0.0, B_screen = 2.0*pi;
const double C_screen = 0.0, D_screen = pi ;

const double E_screen = 0.0, F_screen = 1.0;

//шаг  для n мерных интегральных уравнений
const double h1_obj = (B_obj - A_obj) / _n;
const double h2_obj = (D_obj - C_obj) / _n;
const double h3_obj = (F_obj - E_obj) / _n;

const double h1_screen = (B_screen - A_screen) / _n;
const double h2_screen = (D_screen - C_screen) / _n;
const double h3_screen = (F_screen - E_screen) / _n;

const double x1_start_obj = 0.0, x2_start_obj = 0.0, x3_start_obj = 0.0;
const double x1_start_screen = 5.0, x2_start_screen = 0.0, x3_start_screen = 0.0;

#endif CONSTANTS_H