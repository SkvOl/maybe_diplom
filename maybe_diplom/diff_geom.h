#pragma once
#ifndef DIFF_GEOM_H
#define DIFF_GEOM_H
#include "matrix.h"
#include "integrals.h"

double x(double t1, double t2, double*& var, int k)
{
	switch (k) {
	case 0:return R * cos(t1) * cos(t2);
	case 1:
	{
		double* var1 = NULL;
		var[0] = (x(t1 * (1.0 + mach_eps), t2, var1, 0) - x(t1 * (1.0 - mach_eps), t2, var1, 0)) / 2.0 / t1 / mach_eps;
		var[1] = (x(t1, t2 * (1.0 + mach_eps), var1, 0) - x(t1, t2 * (1.0 - mach_eps), var1, 0)) / 2.0 / t2 / mach_eps;
		return 0;
	}
	}
}

double y(double t1, double t2, double*& var, int k)
{
	switch (k) {
	case 0:return R * cos(t1) * sin(t2);
	case 1:
	{
		double* var1 = NULL;
		var[0] = (y(t1 * (1.0 + mach_eps), t2, var1, 0) - y(t1 * (1.0 - mach_eps), t2, var1, 0)) / 2.0 / t1 / mach_eps;
		var[1] = (y(t1, t2 * (1.0 + mach_eps), var1, 0) - y(t1, t2 * (1.0 - mach_eps), var1, 0)) / 2.0 / t2 / mach_eps;

		return 0;
	}
	}
}

double z(double t1, double t2, double*& var, int k)
{
	switch (k) {
	case 0:return R * sin(t1);
	case 1:
	{
		double* var1 = NULL;
		var[0] = (z(t1 * (1.0 + mach_eps), t2, var1, 0) - z(t1 * (1.0 - mach_eps), t2, var1, 0)) / 2.0 / t1 / mach_eps;
		var[1] = (z(t1, t2 * (1.0 + mach_eps), var1, 0) - z(t1, t2 * (1.0 - mach_eps), var1, 0)) / 2.0 / t2 / mach_eps;
		return 0;
	}
	}
}

void vector_n(double(*function_x)(double, double, double*&, int), double(*function_y)(double, double, double*&, int), double(*function_z)(double, double, double*&, int), double t1, double t2, double*& var)
{
	double* ksi_x = createv<double>(2), * ksi_y = createv<double>(2), * ksi_z = createv<double>(2);
	(*function_x)(t1, t2, ksi_x, 1); (*function_y)(t1, t2, ksi_y, 1); (*function_z)(t1, t2, ksi_z, 1);
	var[0] = (ksi_y[0] * ksi_z[1] - ksi_z[0] * ksi_y[1]) / sqrt(pow(ksi_y[0] * ksi_z[1] - ksi_z[0] * ksi_y[1], 2) + pow(ksi_z[0] * ksi_x[1] - ksi_x[0] * ksi_z[1], 2) + pow(ksi_x[0] * ksi_y[1] - ksi_y[0] * ksi_x[1], 2));
	var[1] = (ksi_z[0] * ksi_x[1] - ksi_x[0] * ksi_z[1]) / sqrt(pow(ksi_y[0] * ksi_z[1] - ksi_z[0] * ksi_y[1], 2) + pow(ksi_z[0] * ksi_x[1] - ksi_x[0] * ksi_z[1], 2) + pow(ksi_x[0] * ksi_y[1] - ksi_y[0] * ksi_x[1], 2));
	var[2] = (ksi_x[0] * ksi_y[1] - ksi_y[0] * ksi_x[1]) / sqrt(pow(ksi_y[0] * ksi_z[1] - ksi_z[0] * ksi_y[1], 2) + pow(ksi_z[0] * ksi_x[1] - ksi_x[0] * ksi_z[1], 2) + pow(ksi_x[0] * ksi_y[1] - ksi_y[0] * ksi_x[1], 2));
	del<double>(ksi_x); del<double>(ksi_y); del<double>(ksi_z);
}

double det_g(double(*function_x)(double, double, double*&, int), double(*function_y)(double, double, double*&, int), double(*function_z)(double, double, double*&, int), double t1, double t2, double**& var)
{
	double* ksi_x = createv<double>(2), * ksi_y = createv<double>(2), * ksi_z = createv<double>(2);
	(*function_x)(t1, t2, ksi_x, 1); (*function_y)(t1, t2, ksi_y, 1); (*function_z)(t1, t2, ksi_z, 1);
	var[0][0] = ksi_x[0] * ksi_x[0] + ksi_y[0] * ksi_y[0] + ksi_z[0] * ksi_z[0];
	var[0][1] = ksi_x[0] * ksi_x[1] + ksi_y[0] * ksi_y[1] + ksi_z[0] * ksi_z[1];
	var[1][0] = var[0][1];
	var[1][1] = ksi_x[1] * ksi_x[1] + ksi_y[1] * ksi_y[1] + ksi_z[1] * ksi_z[1];
	del<double>(ksi_x); del<double>(ksi_y); del<double>(ksi_z);
	return var[0][0] * var[1][1] - var[0][1] * var[1][0];
}

double I(int N_i, double(*function_x)(double, double, double*&, int), double(*function_y)(double, double, double*&, int), double(*function_z)(double, double, double*&, int), double** var, double t1_a, double t1_b, double t2_c, double t2_d)
{
	double h_i = (t1_b - t1_a) / N_i,
		h_j = (t2_d - t2_c) / N_i;
	double Sum = 0;
	for (size_t i = 0; i < N_i; i++) {
		for (size_t j = 0; j < N_i; j++) {
			Sum += sqrt(det_g(x, y, z, t1_a + (i + 0.5) * h_i, t2_c + (j + 0.5) * h_j, var));
		}
	}
	return h_i * h_j * Sum;
}

#endif DIFF_GEOM_H