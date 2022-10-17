#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <math.h>
#include <stdio.h> 

int n = 20, N = n * n;

const double mach_eps = sqrt(2.22045e-16);
const double pi = acos(-1);
const double k0 = 2.0 * pi / (3 * pow(10, 10)) * (5 * pow(10, 10)), k1 = 1.5 * k0;
//const double k0 = 0.01, k1 = 1.5 * k0;
//const double k0=1, k1 = 1.5 * k0;
const double R = 1.0;
const double lambda = 1.0;

//îòðåçîê äëÿ n ìåðíûõ èíòåãðàëüíûõ óðàâíåíèé
double A = 0, B = 5;
double C = 0, D = 5;
double E = 0, F = 5;

//øàã  äëÿ n ìåðíûõ èíòåãðàëüíûõ óðàâíåíèé
double h1 = (B - A) / n;
double h2 = (D - C) / n;
double h3 = (F - E) / n;

#endif CONSTANTS_H
