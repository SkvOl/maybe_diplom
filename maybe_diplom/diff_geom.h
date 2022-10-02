#pragma once
#ifndef DIFF_GEOM_H
#define DIFF_GEOM_H
#include "matrix.h"
#include "integrals.h"

inline void base_func(double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double t1, double t2, int down_index1, int down_index2, int up_index1, int up_index2, double** var);
inline complex<double> k_c(double x1, ...);


double x1_screen(double t1, double t2, double* var, int k)
{
	switch (k) {
	case 0:return R * cos(t1) * cos(t2);
	case 1:
	{
		double* var1 = NULL;
		var[0] = (x1_screen(t1 * (1.0 + mach_eps), t2, var1, 0) - x1_screen(t1 * (1.0 - mach_eps), t2, var1, 0)) / 2.0 / (t1 + mach_eps) / mach_eps;
		var[1] = (x1_screen(t1, t2 * (1.0 + mach_eps), var1, 0) - x1_screen(t1, t2 * (1.0 - mach_eps), var1, 0)) / 2.0 / (t2 + mach_eps) / mach_eps;
		return 0;
	}
	}
}

double x2_screen(double t1, double t2, double* var, int k)
{
	switch (k) {
	case 0:return R * cos(t1) * sin(t2);
	case 1:
	{
		double* var1 = NULL;
		var[0] = (x2_screen(t1 * (1.0 + mach_eps), t2, var1, 0) - x2_screen(t1 * (1.0 - mach_eps), t2, var1, 0)) / 2.0 / (t1 + mach_eps) / mach_eps;
		var[1] = (x2_screen(t1, t2 * (1.0 + mach_eps), var1, 0) - x2_screen(t1, t2 * (1.0 - mach_eps), var1, 0)) / 2.0 / (t2 + mach_eps) / mach_eps;

		return 0;
	}
	}
}

double x3_screen(double t1, double t2, double* var, int k)
{
	switch (k) {
	case 0:return R * sin(t1);
	case 1:
	{
		double* var1 = NULL;
		var[0] = (x3_screen(t1 * (1.0 + mach_eps), t2, var1, 0) - x3_screen(t1 * (1.0 - mach_eps), t2, var1, 0)) / 2.0 / (t1 + mach_eps) / mach_eps;
		var[1] = (x3_screen(t1, t2 * (1.0 + mach_eps), var1, 0) - x3_screen(t1, t2 * (1.0 - mach_eps), var1, 0)) / 2.0 / (t2 + mach_eps) / mach_eps;
		return 0;
	}
	}
}

void vector_n(double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double t1, double t2, double* var)
{
	double* ksi_x1 = createv<double>(2), * ksi_x2 = createv<double>(2), * ksi_x3 = createv<double>(2);
	(*function_x1)(t1, t2, ksi_x1, 1); (*function_x2)(t1, t2, ksi_x2, 1); (*function_x3)(t1, t2, ksi_x3, 1);
	var[0] = (ksi_x2[0] * ksi_x3[1] - ksi_x3[0] * ksi_x2[1]) / sqrt(pow(ksi_x2[0] * ksi_x3[1] - ksi_x3[0] * ksi_x2[1], 2) + pow(ksi_x3[0] * ksi_x1[1] - ksi_x1[0] * ksi_x3[1], 2) + pow(ksi_x1[0] * ksi_x2[1] - ksi_x2[0] * ksi_x1[1], 2));
	var[1] = (ksi_x3[0] * ksi_x1[1] - ksi_x1[0] * ksi_x3[1]) / sqrt(pow(ksi_x2[0] * ksi_x3[1] - ksi_x3[0] * ksi_x2[1], 2) + pow(ksi_x3[0] * ksi_x1[1] - ksi_x1[0] * ksi_x3[1], 2) + pow(ksi_x1[0] * ksi_x2[1] - ksi_x2[0] * ksi_x1[1], 2));
	var[2] = (ksi_x1[0] * ksi_x2[1] - ksi_x2[0] * ksi_x1[1]) / sqrt(pow(ksi_x2[0] * ksi_x3[1] - ksi_x3[0] * ksi_x2[1], 2) + pow(ksi_x3[0] * ksi_x1[1] - ksi_x1[0] * ksi_x3[1], 2) + pow(ksi_x1[0] * ksi_x2[1] - ksi_x2[0] * ksi_x1[1], 2));
	del(ksi_x1); del(ksi_x2); del(ksi_x3);
}

double det_g(double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double t1, double t2, double** var)
{
	double* ksi_x1 = createv<double>(2), * ksi_x2 = createv<double>(2), * ksi_x3 = createv<double>(2);
	(*function_x1)(t1, t2, ksi_x1, 1); (*function_x2)(t1, t2, ksi_x2, 1); (*function_x3)(t1, t2, ksi_x3, 1);
	var[0][0] = ksi_x1[0] * ksi_x1[0] + ksi_x2[0] * ksi_x2[0] + ksi_x3[0] * ksi_x3[0];
	var[0][1] = ksi_x1[0] * ksi_x1[1] + ksi_x2[0] * ksi_x2[1] + ksi_x3[0] * ksi_x3[1];
	var[1][0] = var[0][1];
	var[1][1] = ksi_x1[1] * ksi_x1[1] + ksi_x2[1] * ksi_x2[1] + ksi_x3[1] * ksi_x3[1];
	del(ksi_x1); del(ksi_x2); del(ksi_x3);
	return var[0][0] * var[1][1] - var[0][1] * var[1][0];
}

double det_g_reverse(double** var1, double** var2) {
	double coef = 1.0 / (var1[0][0] * var1[1][1] - var1[0][1] * var1[1][0]);
	var2[0][0] = var1[1][1] / coef;
	var2[0][1] = -var1[0][1] / coef;
	var2[1][0] = -var1[1][0] / coef;
	var2[1][1] = var1[0][0] / coef;

	return var2[0][0] * var2[1][1] - var2[0][1] * var2[1][0];
}

void Jacobian(double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double t1, double t2, double** var) {
	double* ksi_x1 = createv<double>(2), * ksi_x2 = createv<double>(2), * ksi_x3 = createv<double>(2);
	(*function_x1)(t1, t2, ksi_x1, 1); (*function_x2)(t1, t2, ksi_x2, 1); (*function_x3)(t1, t2, ksi_x3, 1);
	var[0][0] = ksi_x1[0]; var[0][1] = ksi_x1[1];
	var[1][0] = ksi_x2[0]; var[1][1] = ksi_x2[1];
	var[2][0] = ksi_x3[0]; var[2][1] = ksi_x3[1];
	del(ksi_x1); del(ksi_x2); del(ksi_x3);
}

template<typename type_A_v>
inline void A_v(int N_i, int k, type_A_v(*function)(double, ...), double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double** var1, type_A_v** var2, double t1_a, double t1_b, double t2_c, double t2_d, double t1, double t2, int down_index1, int down_index2, int up_index2)
{
	var2[0][0] = 0.0;
	var2[1][0] = 0.0;
	var2[2][0] = 0.0;

	double h_i = (t1_b - t1_a) / N_i,
		h_j = (t2_d - t2_c) / N_i;

	double** v = createm<double>(3, 1);
	type_A_v Sum = 0;

	switch (k) {
	case 0: 
	{
		for (size_t i = 0; i < N_i; i++) {
			for (size_t j = 0; j < N_i; j++) {
				double l1 = t1_a + (i + 0.5) * h_i,
					l2 = t2_c + (j + 0.5) * h_j;

				base_func((*function_x1), (*function_x2), (*function_x3), l1, l2, down_index1, down_index2, 0, up_index2, v);
				type_A_v ker = (*function)((*function_x1)(t1, t2, NULL, 0), (*function_x2)(t1, t2, NULL, 0), (*function_x3)(t1, t2, NULL, 0), (*function_x1)(l1, l2, NULL, 0), (*function_x2)(l1, l2, NULL, 0), (*function_x3)(l1, l2, NULL, 0));
				double sq_det = sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, var1));

				var2[0][0] += ker * v[0][0] * sq_det;
				var2[1][0] += ker * v[1][0] * sq_det;
				var2[2][0] += ker * v[1][0] * sq_det;
			}
		}
		var2[0][0] *= h_i * h_j;
		var2[1][0] *= h_i * h_j;
		var2[2][0] *= h_i * h_j;
	}
	case 1:
	{
		for (size_t i = 0; i < N_i; i++) {
			for (size_t j = 0; j < N_i; j++) {
				double l1 = t1_a + (i + 0.5) * h_i,
					l2 = t2_c + (j + 0.5) * h_j;

				base_func((*function_x1), (*function_x2), (*function_x3), l1, l2, down_index1, down_index2, 0, up_index2, v);
				type_A_v ker1 = (*function)((*function_x1)(t1 * (1.0 + mach_eps), t2, NULL, 0), (*function_x2)(t1 * (1.0 + mach_eps), t2, NULL, 0), (*function_x3)(t1 * (1.0 + mach_eps), t2, NULL, 0), (*function_x1)(l1, l2, NULL, 0), (*function_x2)(l1, l2, NULL, 0), (*function_x3)(l1, l2, NULL, 0));
				type_A_v ker2 = (*function)((*function_x1)(t1 * (1.0 - mach_eps), t2, NULL, 0), (*function_x2)(t1 * (1.0 - mach_eps), t2, NULL, 0), (*function_x3)(t1 * (1.0 - mach_eps), t2, NULL, 0), (*function_x1)(l1, l2, NULL, 0), (*function_x2)(l1, l2, NULL, 0), (*function_x3)(l1, l2, NULL, 0));
				double sq_det = sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, var1));

				var2[0][0] += (ker1 - ker2) * v[0][0] * sq_det;
				var2[1][0] += (ker1 - ker2) * v[1][0] * sq_det;
				var2[2][0] += (ker1 - ker2) * v[1][0] * sq_det;
			}
		}
		var2[0][0] *= h_i * h_j / 2.0 / (t1 + mach_eps) / mach_eps;
		var2[1][0] *= h_i * h_j / 2.0 / (t1 + mach_eps) / mach_eps;
		var2[2][0] *= h_i * h_j / 2.0 / (t1 + mach_eps) / mach_eps;
	}
	case 2: 
	{
		for (size_t i = 0; i < N_i; i++) {
			for (size_t j = 0; j < N_i; j++) {
				double l1 = t1_a + (i + 0.5) * h_i,
					l2 = t2_c + (j + 0.5) * h_j;

				base_func((*function_x1), (*function_x2), (*function_x3), l1, l2, down_index1, down_index2, 0, up_index2, v);
				type_A_v ker1 = (*function)((*function_x1)(t1, t2 * (1.0 + mach_eps), NULL, 0), (*function_x2)(t1, t2 * (1.0 + mach_eps), NULL, 0), (*function_x3)(t1, t2 * (1.0 + mach_eps), NULL, 0), (*function_x1)(l1, l2, NULL, 0), (*function_x2)(l1, l2, NULL, 0), (*function_x3)(l1, l2, NULL, 0));
				type_A_v ker2 = (*function)((*function_x1)(t1, t2 * (1.0 - mach_eps), NULL, 0), (*function_x2)(t1, t2 * (1.0 - mach_eps), NULL, 0), (*function_x3)(t1, t2 * (1.0 - mach_eps), NULL, 0), (*function_x1)(l1, l2, NULL, 0), (*function_x2)(l1, l2, NULL, 0), (*function_x3)(l1, l2, NULL, 0));
				double sq_det = sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, var1));
				
				var2[0][0] += (ker1 - ker2) * v[0][0] * sq_det;
				var2[1][0] += (ker1 - ker2) * v[1][0] * sq_det;
				var2[2][0] += (ker1 - ker2) * v[1][0] * sq_det;
			}
		}
		var2[0][0] *= h_i * h_j / 2.0 / (t2 + mach_eps) / mach_eps;
		var2[1][0] *= h_i * h_j / 2.0 / (t2 + mach_eps) / mach_eps;
		var2[2][0] *= h_i * h_j / 2.0 / (t2 + mach_eps) / mach_eps;
	}
	}
	del(v);
}

inline void d_x(short type, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double t1, double t2, double** var) {
	if (type == 1) {
		double* ksi_x1 = createv<double>(2), * ksi_x2 = createv<double>(2), * ksi_x3 = createv<double>(2);
		(*function_x1)(t1, t2, ksi_x1, 1); (*function_x2)(t1, t2, ksi_x2, 1); (*function_x3)(t1, t2, ksi_x3, 1);
		var[0][0] = ksi_x1[0];
		var[1][0] = ksi_x2[0];
		var[2][0] = ksi_x3[0];
		del(ksi_x1); del(ksi_x2); del(ksi_x3);
	}
	else if (type == 2) {
		double* ksi_x1 = createv<double>(2), * ksi_x2 = createv<double>(2), * ksi_x3 = createv<double>(2);
		(*function_x1)(t1, t2, ksi_x1, 1); (*function_x2)(t1, t2, ksi_x2, 1); (*function_x3)(t1, t2, ksi_x3, 1);
		var[0][0] = ksi_x1[1];
		var[1][0] = ksi_x2[1];
		var[2][0] = ksi_x3[1];
		del(ksi_x1); del(ksi_x2); del(ksi_x3);
	}
}

inline double I(int N_i, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double** var, double t1_a, double t1_b, double t2_c, double t2_d)
{
	double h_i = (t1_b - t1_a) / N_i,
		h_j = (t2_d - t2_c) / N_i;
	double Sum = 0;
	for (size_t i = 0; i < N_i; i++) {
		for (size_t j = 0; j < N_i; j++) {
			double l1 = t1_a + (i + 0.5) * h_i,
				   l2 = t2_c + (j + 0.5) * h_j;
			Sum += sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, var));
		}
	}
	return h_i * h_j * Sum;
}

inline complex<double> S(int N_i, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double a, double b, double c, double d, double e, double f, double g, double h, int down_index1, int down_index2, int up_index2) {
	double h_i = (b - a) / N_i,
		h_j = (d - c) / N_i;
	complex<double> Sum = 0;

	for (size_t i = 0; i < N_i; i++) {
		for (size_t j = 0; j < N_i; j++) {
			double l1 = a + (i + 0.5) * h_i,
				l2 = c + (j + 0.5) * h_j;

			double** tensor = createm<double>(2, 2), ** tensor_reverse = createm<double>(2, 2);

			double sq_det = det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, tensor);
			det_g_reverse(tensor, tensor_reverse);

			for (size_t mu = 0; mu < 2; mu++)
			{
				for (size_t nu = 0; nu < 2; nu++)
				{
					double** private_diff_x = createm<double>(3, 1);
					complex<double>** Sum3 = createm<complex<double>>(3, 1, true);
					d_x(nu + 1, (*function_x1), (*function_x2), (*function_x3), l1, l2, private_diff_x);

					for (size_t alpha = 0; alpha < 3; alpha++)
					{
						for (size_t beta = 0; beta < 3; beta++)
						{
							double** private_diff_x2 = createm<double>(3, 1);
							complex<double>** A_v_res1 = createm<complex<double>>(3, 1), ** proizv = createm<complex<double>>(3, 1, true);
							A_v(3, alpha + 1, k_c, x1_screen, x2_screen, x3_screen, tensor, A_v_res1, e, f, g, h, l1, l2, down_index1, down_index2, up_index2);
							d_x(beta + 1, (*function_x1), (*function_x2), (*function_x3), l1, l2, private_diff_x2);
						}
					}

					del(private_diff_x);
				}
			}
			del(tensor); del(tensor_reverse);
		}
	}
}
#endif DIFF_GEOM_H