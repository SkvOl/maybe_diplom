#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <math.h>

int _n = 2, _N = 2 * _n * _n - 2 * _n;

const double mach_eps = sqrt(2.22045e-16);
const double pi = acos(-1);
const complex<double> k0(1, 0.001), k1 = 1.5 * k0;
const double R = 1.0;
const double lambda = 1.0;

//отрезок для n мерных интегральных уравнений
double A = -pi / 2.0, B = pi / 2.0;
double C = 0.0, D = 2.0 * pi;
double E = 0, F = 1;

//шаг  для n мерных интегральных уравнений
double h1 = (B - A) / _n;
double h2 = (D - C) / _n;
double h3 = (F - E) / _n;

#endif CONSTANTS_H