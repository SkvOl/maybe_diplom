#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <math.h>

int _n = 2, _N = _n * _n;

//������� ��� n ������ ������������ ���������
double A = 0, B = 1;
double C = 0, D = 1;
double E = 0, F = 1;

//���  ��� n ������ ������������ ���������
double h1 = (B - A) / _n;
double h2 = (D - C) / _n;
double h3 = (F - E) / _n;

const double mach_eps = sqrt(2.22045e-16);
const double pi = acos(-1);
const complex<double> k0(1, 0.001), k1 = 1.5 * k0;
const double R = 1;
const double lambda = 1.0;

#endif CONSTANTS_H