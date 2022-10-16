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

inline void Jacobian(double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double t1, double t2, double** var) {
	double* ksi_x1 = createv<double>(2), * ksi_x2 = createv<double>(2), * ksi_x3 = createv<double>(2);
	(*function_x1)(t1, t2, ksi_x1, 1); (*function_x2)(t1, t2, ksi_x2, 1); (*function_x3)(t1, t2, ksi_x3, 1);
	var[0][0] = ksi_x1[0]; var[0][1] = ksi_x1[1];
	var[1][0] = ksi_x2[0]; var[1][1] = ksi_x2[1];
	var[2][0] = ksi_x3[0]; var[2][1] = ksi_x3[1];
	del(ksi_x1); del(ksi_x2); del(ksi_x3);
}

template<typename type_A_v>
inline void A_v(int N_i, int k, type_A_v(*function)(double, ...), double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double** var1, type_A_v** var2,double t1, double t2, int down_index1, int down_index2, int up_index2)
{
	var2[0][0] = 0.0;
	var2[1][0] = 0.0;
	var2[2][0] = 0.0;

	double e = A + down_index1 * h1, f = e + h1;
	double g = C + down_index2 * h2, l = g + h2;

	double h_i = (f - e) / N_i,
		   h_j = (l - g) / N_i;

	double** v = createm<double>(3, 1);


	switch (k) {
	case 0: {
		for (size_t i = 0; i < N_i; i++) {
			for (size_t j = 0; j < N_i; j++) {
				double l1 = e + (i + 0.5) * h_i,
					   l2 = g + (j + 0.5) * h_j;

				base_func((*function_x1), (*function_x2), (*function_x3), l1, l2, down_index1, down_index2, 0, up_index2, v);
				type_A_v ker = (*function)((*function_x1)(t1, t2, NULL, 0), (*function_x2)(t1, t2, NULL, 0), (*function_x3)(t1, t2, NULL, 0), (*function_x1)(l1, l2, NULL, 0), (*function_x2)(l1, l2, NULL, 0), (*function_x3)(l1, l2, NULL, 0));
				double sq_det = sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, var1));

				var2[0][0] += ker * v[0][0] * sq_det;
				var2[1][0] += ker * v[1][0] * sq_det;
				var2[2][0] += ker * v[2][0] * sq_det;
			}
		}
		var2[0][0] *= h_i * h_j;
		var2[1][0] *= h_i * h_j;
		var2[2][0] *= h_i * h_j;
		break;
	}
	case 1: {
		for (size_t i = 0; i < N_i; i++) {
			for (size_t j = 0; j < N_i; j++) {
				double l1 = e + (i + 0.5) * h_i,
					l2 = g + (j + 0.5) * h_j;

				base_func((*function_x1), (*function_x2), (*function_x3), l1, l2, down_index1, down_index2, 0, up_index2, v);
				type_A_v ker1 = (*function)((*function_x1)(t1 * (1.0 + mach_eps), t2, NULL, 0), (*function_x2)(t1 * (1.0 + mach_eps), t2, NULL, 0), (*function_x3)(t1 * (1.0 + mach_eps), t2, NULL, 0), (*function_x1)(l1, l2, NULL, 0), (*function_x2)(l1, l2, NULL, 0), (*function_x3)(l1, l2, NULL, 0));
				type_A_v ker2 = (*function)((*function_x1)(t1 * (1.0 - mach_eps), t2, NULL, 0), (*function_x2)(t1 * (1.0 - mach_eps), t2, NULL, 0), (*function_x3)(t1 * (1.0 - mach_eps), t2, NULL, 0), (*function_x1)(l1, l2, NULL, 0), (*function_x2)(l1, l2, NULL, 0), (*function_x3)(l1, l2, NULL, 0));
				double sq_det = sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, var1));

				var2[0][0] += (ker1 - ker2) * v[0][0] * sq_det;
				var2[1][0] += (ker1 - ker2) * v[1][0] * sq_det;
				var2[2][0] += (ker1 - ker2) * v[2][0] * sq_det;
			}
		}
		var2[0][0] *= h_i * h_j / 2.0 / (t1 + mach_eps) / mach_eps;
		var2[1][0] *= h_i * h_j / 2.0 / (t1 + mach_eps) / mach_eps;
		var2[2][0] *= h_i * h_j / 2.0 / (t1 + mach_eps) / mach_eps;
		break;
	}
	case 2: {
		for (size_t i = 0; i < N_i; i++) {
			for (size_t j = 0; j < N_i; j++) {
				double l1 = e + (i + 0.5) * h_i,
					l2 = g + (j + 0.5) * h_j;

				base_func((*function_x1), (*function_x2), (*function_x3), l1, l2, down_index1, down_index2, 0, up_index2, v);
				type_A_v ker1 = (*function)((*function_x1)(t1, t2 * (1.0 + mach_eps), NULL, 0), (*function_x2)(t1, t2 * (1.0 + mach_eps), NULL, 0), (*function_x3)(t1, t2 * (1.0 + mach_eps), NULL, 0), (*function_x1)(l1, l2, NULL, 0), (*function_x2)(l1, l2, NULL, 0), (*function_x3)(l1, l2, NULL, 0));
				type_A_v ker2 = (*function)((*function_x1)(t1, t2 * (1.0 - mach_eps), NULL, 0), (*function_x2)(t1, t2 * (1.0 - mach_eps), NULL, 0), (*function_x3)(t1, t2 * (1.0 - mach_eps), NULL, 0), (*function_x1)(l1, l2, NULL, 0), (*function_x2)(l1, l2, NULL, 0), (*function_x3)(l1, l2, NULL, 0));
				double sq_det = sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, var1));

				var2[0][0] += (ker1 - ker2) * v[0][0] * sq_det;
				var2[1][0] += (ker1 - ker2) * v[1][0] * sq_det;
				var2[2][0] += (ker1 - ker2) * v[2][0] * sq_det;
			}
		}
		var2[0][0] *= h_i * h_j / 2.0 / (t2 + mach_eps) / mach_eps;
		var2[1][0] *= h_i * h_j / 2.0 / (t2 + mach_eps) / mach_eps;
		var2[2][0] *= h_i * h_j / 2.0 / (t2 + mach_eps) / mach_eps;
		break;
	}
	}
	del(v);
}

inline void d_x(short type, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double t1, double t2, double* var) {
	if (type == 1) {
		double* ksi_x1 = createv<double>(2), * ksi_x2 = createv<double>(2), * ksi_x3 = createv<double>(2);
		(*function_x1)(t1, t2, ksi_x1, 1); (*function_x2)(t1, t2, ksi_x2, 1); (*function_x3)(t1, t2, ksi_x3, 1);
		var[0] = ksi_x1[0];
		var[1] = ksi_x2[0];
		var[2] = ksi_x3[0];
		del(ksi_x1); del(ksi_x2); del(ksi_x3);
	}
	else if (type == 2) {
		double* ksi_x1 = createv<double>(2), * ksi_x2 = createv<double>(2), * ksi_x3 = createv<double>(2);
		(*function_x1)(t1, t2, ksi_x1, 1); (*function_x2)(t1, t2, ksi_x2, 1); (*function_x3)(t1, t2, ksi_x3, 1);
		var[0] = ksi_x1[1];
		var[1] = ksi_x2[1];
		var[2] = ksi_x3[1];
		del(ksi_x1); del(ksi_x2); del(ksi_x3);
	}
}

inline double I(int N_i, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double** var, double a, double b, double c, double d)
{
	double h_i = (b - a) / N_i,
		   h_j = (d - c) / N_i;
	
	double Sum = 0;
	for (size_t i = 0; i < N_i; i++) {
		for (size_t j = 0; j < N_i; j++) {
			double l1 = a + (i + 0.5) * h_i,
				   l2 = c + (j + 0.5) * h_j;
			Sum += sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, var));
		}
	}
	return h_i * h_j * Sum;
}

template<typename type_div>
inline type_div div(short k, double** tensor, double **tensor_reverse, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double t1, double t2, int down_index1, int down_index2, int up_index2) {
	size_t M = _msize(tensor_reverse) / sizeof(tensor_reverse[0]);
	size_t N = _msize(tensor_reverse[0]) / sizeof(tensor_reverse[0][0]);

	type_div Sum = 0;

	switch (k) {
	case 0: {
		det_g((*function_x1), (*function_x2), (*function_x3), t1, t2, tensor);
		det_g_reverse(tensor, tensor_reverse);
		type_div** A_v_res;
		double* x;
		for (size_t mu = 0; mu < M; mu++)
		{
			for (size_t nu = 0; nu < N; nu++)
			{
				A_v_res = createm<type_div>(3, 1);
				x = createv<double>(3);

				A_v(2, mu + 1, k_c, (*function_x1), (*function_x2), (*function_x3), tensor, A_v_res, t1, t2, down_index1, down_index2, up_index2);
				d_x(nu + 1, (*function_x1), (*function_x2), (*function_x3), t1, t2, x);

				Sum += multv2(A_v_res, x) * tensor_reverse[mu][nu];

				del(A_v_res); del(x);
			}
		}
		break;
	}
	case 1: {
		double** tensor2 = createm<double>(M, N), ** tensor_reverse2 = createm<double>(M, N);

		det_g((*function_x1), (*function_x2), (*function_x3), t1 * (1.0 + mach_eps), t2, tensor);
		det_g((*function_x1), (*function_x2), (*function_x3), t1 * (1.0 - mach_eps), t2, tensor2);
		det_g_reverse(tensor, tensor_reverse);
		det_g_reverse(tensor2, tensor_reverse2);
		
		type_div** A_v_res1, ** A_v_res2;
		double* x1, *x2;
		for (size_t mu = 0; mu < M; mu++)
		{
			for (size_t nu = 0; nu < N; nu++)
			{
				A_v_res1 = createm<type_div>(3, 1); A_v_res2 = createm<type_div>(3, 1);
				x1 = createv<double>(3); x2 = createv<double>(3);

				A_v(2, mu + 1, k_c, (*function_x1), (*function_x2), (*function_x3), tensor, A_v_res1, t1 * (1.0 + mach_eps), t2, down_index1, down_index2, up_index2);
				A_v(2, mu + 1, k_c, (*function_x1), (*function_x2), (*function_x3), tensor2, A_v_res2, t1 * (1.0 - mach_eps), t2, down_index1, down_index2, up_index2);
				d_x(nu + 1, (*function_x1), (*function_x2), (*function_x3), t1 * (1.0 + mach_eps), t2, x1);
				d_x(nu + 1, (*function_x1), (*function_x2), (*function_x3), t1 * (1.0 - mach_eps), t2, x2);
				
				Sum += (multv2(A_v_res1, x1) * tensor_reverse[mu][nu] - multv2(A_v_res2, x2) * tensor_reverse2[mu][nu]);

				del(A_v_res1); del(x1); 
				del(A_v_res2); del(x2);
			}
		}
		del(tensor2); del(tensor_reverse2);
		Sum /= ((t1 + mach_eps) * mach_eps);
		break;
	}
	case 2: {
		double** tensor2 = createm<double>(M, N), ** tensor_reverse2 = createm<double>(M, N);

		det_g((*function_x1), (*function_x2), (*function_x3), t1, t2 * (1.0 + mach_eps), tensor);
		det_g((*function_x1), (*function_x2), (*function_x3), t1, t2 * (1.0 - mach_eps), tensor2);
		det_g_reverse(tensor, tensor_reverse);
		det_g_reverse(tensor2, tensor_reverse2);
		
		type_div** A_v_res1, ** A_v_res2;
		double* x1, * x2;
		for (size_t mu = 0; mu < M; mu++)
		{
			for (size_t nu = 0; nu < N; nu++)
			{
				A_v_res1 = createm<type_div>(3, 1); A_v_res2 = createm<type_div>(3, 1);
				x1 = createv<double>(3); x2 = createv<double>(3);

				A_v(2, mu + 1, k_c, (*function_x1), (*function_x2), (*function_x3), tensor, A_v_res1, t1, t2 * (1.0 + mach_eps), down_index1, down_index2, up_index2);
				A_v(2, mu + 1, k_c, (*function_x1), (*function_x2), (*function_x3), tensor2, A_v_res2, t1, t2 * (1.0 - mach_eps), down_index1, down_index2, up_index2);
				d_x(nu + 1, (*function_x1), (*function_x2), (*function_x3), t1, t2 * (1.0 + mach_eps), x1);
				d_x(nu + 1, (*function_x1), (*function_x2), (*function_x3), t1, t2 * (1.0 - mach_eps), x2);

				Sum += (multv2(A_v_res1, x1) * tensor_reverse[mu][nu] - multv2(A_v_res2, x2) * tensor_reverse2[mu][nu]);


				del(A_v_res1); del(x1);
				del(A_v_res2); del(x2);
			}
		}
		del(tensor2); del(tensor_reverse2);
		Sum /= ((t2 + mach_eps) * mach_eps);
		break;
	}
	}
	return Sum;
}

template<typename type_grad>
inline void grad(double** tensor, double** tensor_reverse, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double t1, double t2, int down_index1, int down_index2, int up_index2, type_grad* var) {
	size_t M = _msize(tensor_reverse) / sizeof(tensor_reverse[0]);
	size_t N = _msize(tensor_reverse[0]) / sizeof(tensor_reverse[0][0]);
	
	det_g((*function_x1), (*function_x2), (*function_x3), t1, t2, tensor);
	det_g_reverse(tensor, tensor_reverse);

	type_grad div_item;
	double* x;
	for (size_t mu = 0; mu < M; mu++)
	{
		for (size_t nu = 0; nu < N; nu++)
		{
			x = createv<double>(3);
			d_x(nu + 1, (*function_x1), (*function_x2), (*function_x3), t1, t2, x);
			div_item = div<type_grad>(mu + 1, tensor, tensor_reverse, (*function_x1), (*function_x2), (*function_x3), t1, t2, down_index1, down_index2, up_index2);
			
			var[0] += tensor_reverse[mu][nu] * div_item * x[0];
			var[1] += tensor_reverse[mu][nu] * div_item * x[1];
			var[2] += tensor_reverse[mu][nu] * div_item * x[2];
			del(x);
		}
	}

}

template<typename type_S>
inline type_S S(int N_i, double** tensor, double** tensor_reverse, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), int k, int l, int i1, int i2, int j1, int j2) {
	double a = A + i1 * h1, b = a + h1;
	double c = C + i2 * h2, d = c + h2;


	double h_i = (b - a) / N_i,
		   h_j = (d - c) / N_i;

	type_S Sum = 0;
	type_S* grad_res, ** A_v_res, * Sum_ker = createv<type_S>(3);
	double** v;
	for (size_t i = 0; i < N_i; i++) {
		for (size_t j = 0; j < N_i; j++) {
			double l1 = a + (i + 0.5) * h_i,
				   l2 = c + (j + 0.5) * h_j;
			//cout << l1 << " " << l2 << "\n";
			grad_res = createv<type_S>(3, true);
			A_v_res = createm<type_S>(3, 1);
			v = createm<double>(3, 1);

			grad(tensor, tensor_reverse, (*function_x1), (*function_x2), (*function_x3), l1, l2, j1, j2, l, grad_res);
			A_v(N_i, 0, k_c, (*function_x1), (*function_x2), (*function_x3), tensor, A_v_res, l1, l2, j1, j2, l);
			base_func((*function_x1), (*function_x2), (*function_x3), l1, l2, i1, i2, 0, k, v);
			
			Sum_ker[0] = grad_res[0] + k0 * k0 * A_v_res[0][0];
			Sum_ker[1] = grad_res[1] + k0 * k0 * A_v_res[1][0];
			Sum_ker[2] = grad_res[2] + k0 * k0 * A_v_res[2][0];

			Sum += multv1<type_S>(v, Sum_ker) * sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, tensor));
			del(v); del(A_v_res); del(grad_res);
		}
	}
	return h_i * h_j * Sum;
}

template<typename type_f>
inline type_f f_vec(int N_i, double** tensor, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), int k, int i1, int i2) {
	double a = A + i1 * h1, b = a + h1;
	double c = C + i2 * h2, d = c + h2;

	double h_i = (b - a) / N_i,
		   h_j = (d - c) / N_i;

	type_f Sum = 0;
	double** v;
	type_f* E0;
	for (size_t i = 0; i < N_i; i++) {
		for (size_t j = 0; j < N_i; j++) {
			double l1 = a + (i + 0.5) * h_i,
				   l2 = c + (j + 0.5) * h_j;
			v = createm<double>(3, 1);
			E0 = createv<type_f>(3);
			base_func((*function_x1), (*function_x2), (*function_x3), l1, l2, i1, i2, 0, k, v);
			func_cv((*function_x1)(l1, l2, NULL, 0), E0);	

			cout << "v:\n";
			print(v);

			Sum += multv1<type_f>(v, E0) * sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, tensor));
			del(v); del(E0);
		}
	}
	return h_i * h_j * Sum;
}

#endif DIFF_GEOM_H