#pragma once
#ifndef DIFF_GEOM_H
#define DIFF_GEOM_H
#include "matrix.h"
#include "integrals.h"

inline void base_func(double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double t1, double t2, int down_index1, int down_index2, int up_index1, int up_index2, double* var, bool object);
inline complex<double> k_c(double x1, ...);


inline double x1_obj(double t1, double t2, double* var, int k)
{
	switch (k) {
	//case 0: return x1_start_obj + R * cos(t1) * cos(t2);
	case 0: return x1_start_obj + R * cos(t1) * sin(t2);
	case 1:
	{
		double* var1 = NULL;
		var[0] = (x1_obj(t1 + mach_eps, t2, var1, 0) - x1_obj(t1 - mach_eps, t2, var1, 0)) / 2.0 / mach_eps;
		var[1] = (x1_obj(t1, t2 + mach_eps, var1, 0) - x1_obj(t1, t2 - mach_eps, var1, 0)) / 2.0 / mach_eps;

		return 0;
	}
	}
}

inline double x2_obj(double t1, double t2, double* var, int k)
{
	switch (k) {
	//case 0:return x2_start_obj + R * cos(t1) * sin(t2);
	case 0: return x2_start_obj + R * sin(t1) * sin(t2);
	case 1:
	{
		double* var1 = NULL;
		var[0] = (x2_obj(t1 + mach_eps, t2, var1, 0) - x2_obj(t1 - mach_eps, t2, var1, 0)) / 2.0 / mach_eps;
		var[1] = (x2_obj(t1, t2 + mach_eps, var1, 0) - x2_obj(t1, t2 - mach_eps, var1, 0)) / 2.0 / mach_eps;

		return 0;
	}
	}
}

inline double x3_obj(double t1, double t2, double* var, int k)
{
	switch (k) {
	//case 0:return x3_start_obj + R * sin(t1);
	case 0:return x3_start_obj + R * cos(t2);
	case 1:
	{
		double* var1 = NULL;
		var[0] = (x3_obj(t1 + mach_eps, t2, var1, 0) - x3_obj(t1 - mach_eps, t2, var1, 0)) / 2.0 / mach_eps;
		var[1] = (x3_obj(t1, t2 + mach_eps, var1, 0) - x3_obj(t1, t2 - mach_eps, var1, 0)) / 2.0 / mach_eps;

		return 0;
	}
	}
}




inline double x1_screen(double t1, double t2, double* var, int k)
{
	switch (k) {
	case 0: return x1_start_screen + R * cos(t1) * sin(t2);
	case 1:
	{
		double* var1 = NULL;
		var[0] = (x1_screen(t1 + mach_eps, t2, var1, 0) - x1_screen(t1 - mach_eps, t2, var1, 0)) / 2.0 / mach_eps;
		var[1] = (x1_screen(t1, t2 + mach_eps, var1, 0) - x1_screen(t1, t2 - mach_eps, var1, 0)) / 2.0 / mach_eps;

		return 0;
	}
	}
}

inline double x2_screen(double t1, double t2, double* var, int k)
{
	switch (k) {
	case 0: return x2_start_screen + R * sin(t1) * sin(t2);
	case 1:
	{
		double* var1 = NULL;
		var[0] = (x2_screen(t1 + mach_eps, t2, var1, 0) - x2_screen(t1 - mach_eps, t2, var1, 0)) / 2.0 / mach_eps;
		var[1] = (x2_screen(t1, t2 + mach_eps, var1, 0) - x2_screen(t1, t2 - mach_eps, var1, 0)) / 2.0 / mach_eps;

		return 0;
	}
	}
}

inline double x3_screen(double t1, double t2, double* var, int k)
{
	switch (k) {
	case 0:return x3_start_screen + R * cos(t2);
	case 1:
	{
		double* var1 = NULL;
		var[0] = (x3_screen(t1 + mach_eps, t2, var1, 0) - x3_screen(t1 - mach_eps, t2, var1, 0)) / 2.0 / mach_eps;
		var[1] = (x3_screen(t1, t2 + mach_eps, var1, 0) - x3_screen(t1, t2 - mach_eps, var1, 0)) / 2.0 / mach_eps;

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
	double coef = var1[0][0] * var1[1][1] - var1[0][1] * var1[1][0];
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
inline void A_v(int N_i, int k, type_A_v(*function)(double, ...), double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double** var1, type_A_v* var2,double t1, double t2, int down_index1, int down_index2, int up_index2, bool object)
{
	var2[0] = 0.0;
	var2[1] = 0.0;
	var2[2] = 0.0;

	double e = 0, f = 0;
	double g = 0, l = 0;

	if (object) {
		e = A_obj + down_index1 * h1_obj; f = up_index2 == 1 ? e + 2.0 * h1_obj : e + h1_obj;
		g = C_obj + down_index2 * h2_obj; l = up_index2 == 2 ? g + 2.0 * h2_obj : g + h2_obj;
	}
	else {
		e = A_screen + down_index1 * h1_screen; f = up_index2 == 1 ? e + 2.0 * h1_screen : e + h1_screen;
		g = C_screen + down_index2 * h2_screen; l = up_index2 == 2 ? g + 2.0 * h2_screen : g + h2_screen;
	}


	double h_i = (f - e) / N_i,
		   h_j = (l - g) / N_i;

	double* v = createv<double>(3);

	switch (k) {
	case 0: {
		for (size_t i = 0; i < N_i; i++) {
			for (size_t j = 0; j < N_i; j++) {
				double l1 = e + (i + 0.65) * h_i,
					   l2 = g + (j + 0.65) * h_j;
				
				base_func((*function_x1), (*function_x2), (*function_x3), l1, l2, down_index1, down_index2, 0, up_index2, v, object);
				type_A_v ker = (*function)((*function_x1)(t1, t2, NULL, 0), (*function_x2)(t1, t2, NULL, 0), (*function_x3)(t1, t2, NULL, 0), (*function_x1)(l1, l2, NULL, 0), (*function_x2)(l1, l2, NULL, 0), (*function_x3)(l1, l2, NULL, 0));
				double sq_det = sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, var1));

				var2[0] += ker * v[0] * sq_det;
				var2[1] += ker * v[1] * sq_det;
				var2[2] += ker * v[2] * sq_det;
			}
		}
		var2[0] *= h_i * h_j;
		var2[1] *= h_i * h_j;
		var2[2] *= h_i * h_j;

		break;
	}
	case 1: {
		complex<double>* var21 = createv<type_A_v>(3), *var22 = createv<type_A_v>(3);
		A_v(N_i, 0, (*function), (*function_x1), (*function_x2), (*function_x3), var1, var21, t1 + mach_eps, t2, down_index1, down_index2, up_index2, object);
		A_v(N_i, 0, (*function), (*function_x1), (*function_x2), (*function_x3), var1, var22, t1 - mach_eps, t2, down_index1, down_index2, up_index2, object);

		var2[0] = (var21[0] - var22[0]) / 2.0 / mach_eps;
		var2[1] = (var21[1] - var22[1]) / 2.0 / mach_eps;
		var2[2] = (var21[2] - var22[2]) / 2.0 / mach_eps;
		
		del(var21); del(var22);
		break;
	}
	case 2: {
		type_A_v* var21 = createv<type_A_v>(3), * var22 = createv<type_A_v>(3);
		A_v(N_i, 0, (*function), (*function_x1), (*function_x2), (*function_x3), var1, var21, t1, t2 + mach_eps, down_index1, down_index2, up_index2, object);
		A_v(N_i, 0, (*function), (*function_x1), (*function_x2), (*function_x3), var1, var22, t1, t2 - mach_eps, down_index1, down_index2, up_index2, object);

		var2[0] = (var21[0] - var22[0]) / 2.0 / mach_eps;
		var2[1] = (var21[1] - var22[1]) / 2.0 / mach_eps;
		var2[2] = (var21[2] - var22[2]) / 2.0 / mach_eps;

		del(var21); del(var22);
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
inline type_div div(short k, double** tensor, double **tensor_reverse, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double t1, double t2, int down_index1, int down_index2, int up_index2, bool object) {
	size_t M = _msize(tensor_reverse) / sizeof(tensor_reverse[0]);
	size_t N = _msize(tensor_reverse[0]) / sizeof(tensor_reverse[0][0]);

	type_div Sum = 0;

	switch (k) {
	case 0: {
		det_g((*function_x1), (*function_x2), (*function_x3), t1, t2, tensor);
		det_g_reverse(tensor, tensor_reverse);
		type_div* A_v_res;
		double* x;
		for (size_t mu = 0; mu < M; mu++)
		{
			for (size_t nu = 0; nu < N; nu++)
			{
				A_v_res = createv<type_div>(3);
				x = createv<double>(3);

				A_v(4, mu + 1, k_c, (*function_x1), (*function_x2), (*function_x3), tensor, A_v_res, t1, t2, down_index1, down_index2, up_index2, object);
				d_x(nu + 1, (*function_x1), (*function_x2), (*function_x3), t1, t2, x);

				Sum += multv3(A_v_res, x) * tensor_reverse[mu][nu];


				del(A_v_res); del(x);
			}
		}
		break;
	}
	case 1: {
		double** tensor2 = createm<double>(M, N), ** tensor_reverse2 = createm<double>(M, N);
		double** tensor22 = createm<double>(M, N), ** tensor_reverse22 = createm<double>(M, N);
		Sum = div<type_div>(0, tensor2, tensor_reverse2, (*function_x1), (*function_x2), (*function_x3), t1 + mach_eps, t2, down_index1, down_index2, up_index2, object);
		Sum -= div<type_div>(0, tensor22, tensor_reverse22, (*function_x1), (*function_x2), (*function_x3), t1 - mach_eps, t2, down_index1, down_index2, up_index2, object);

		Sum /= (2.0 * mach_eps);
		del(tensor22); del(tensor_reverse22);
		
		break;
	}
	case 2: {
		double** tensor2 = createm<double>(M, N), ** tensor_reverse2 = createm<double>(M, N);
		double** tensor22 = createm<double>(M, N), ** tensor_reverse22 = createm<double>(M, N);
		Sum = div<type_div>(0, tensor2, tensor_reverse2, (*function_x1), (*function_x2), (*function_x3), t1, t2 + mach_eps, down_index1, down_index2, up_index2, object);
		Sum -= div<type_div>(0, tensor22, tensor_reverse22, (*function_x1), (*function_x2), (*function_x3), t1, t2 - mach_eps, down_index1, down_index2, up_index2, object);
		Sum /= (2.0 * mach_eps);
		del(tensor22); del(tensor_reverse22);
		break;
	}
	}
	return Sum;
}

template<typename type_grad>
inline void grad(double** tensor, double** tensor_reverse, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double t1, double t2, int down_index1, int down_index2, int up_index2, type_grad* var, bool object) {
	size_t M = _msize(tensor_reverse) / sizeof(tensor_reverse[0]);
	size_t N = _msize(tensor_reverse[0]) / sizeof(tensor_reverse[0][0]);
	
	det_g((*function_x1), (*function_x2), (*function_x3), t1, t2, tensor);
	det_g_reverse(tensor, tensor_reverse);
	
	type_grad div_item;
	double* x;
	for (size_t mu = 0; mu < M; mu++)
	{
		div_item = div<type_grad>(mu + 1, tensor, tensor_reverse, (*function_x1), (*function_x2), (*function_x3), t1, t2, down_index1, down_index2, up_index2, object);
		for (size_t nu = 0; nu < N; nu++)
		{
			x = createv<double>(3);
			d_x(nu + 1, (*function_x1), (*function_x2), (*function_x3), t1, t2, x);
			
			var[0] += tensor_reverse[mu][nu] * div_item * x[0];
			var[1] += tensor_reverse[mu][nu] * div_item * x[1];
			var[2] += tensor_reverse[mu][nu] * div_item * x[2];

			del(x);
		}	
	}
}

template<typename type_S_screen>
inline type_S_screen S_screen(int N_i, double** tensor, double** tensor_reverse, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), int k, int l, int i1, int i2, int j1, int j2) {
	double a = 0, b = 0;
	double c = 0, d = 0;
	
	a = A_screen + i1 * h1_screen; b = k == 1 ? a + 2.0 * h1_screen : a + h1_screen;
	c = C_screen + i2 * h2_screen; d = k == 2 ? c + 2.0 * h2_screen : c + h2_screen;
	

	double h_i = (b - a) / N_i,
		   h_j = (d - c) / N_i;

	type_S_screen Sum = 0;
	type_S_screen * grad_res, * A_v_res, * Sum_ker = createv<type_S_screen>(3);
	double * v, * normal;
	
	for (size_t i = 0; i < N_i; i++) {
		for (size_t j = 0; j < N_i; j++) {
			double l1 = a + (i + 0.5) * h_i,
				l2 = c + (j + 0.5) * h_j;

			grad_res = createv<type_S_screen>(3, true);
			A_v_res = createv<type_S_screen>(3);

			v = createv<double>(3);
			normal = createv<double>(3);

			grad(tensor, tensor_reverse, (*function_x1), (*function_x2), (*function_x3), l1, l2, j1, j2, l, grad_res, false);
			A_v(N_i, 0, k_c, (*function_x1), (*function_x2), (*function_x3), tensor, A_v_res, l1, l2, j1, j2, l, false);
			
			
			base_func((*function_x1), (*function_x2), (*function_x3), l1, l2, i1, i2, 0, k, v, false);
			vector_n((*function_x1), (*function_x2), (*function_x3), l1, l2, normal);

			Sum_ker[0] = grad_res[0] + k0 * k0 * A_v_res[0];
			Sum_ker[1] = grad_res[1] + k0 * k0 * A_v_res[1];
			Sum_ker[2] = grad_res[2] + k0 * k0 * A_v_res[2];
			
			
			complex<double> tem = multv3(Sum_ker, normal);

			Sum_ker[0] -= normal[0] * tem;
			Sum_ker[1] -= normal[1] * tem;
			Sum_ker[2] -= normal[2] * tem;

			

			Sum += multv3<type_S_screen>(Sum_ker,v) * sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, tensor));
			del(v); del(normal); del(A_v_res); del(grad_res);
		}
	}
	return h_i * h_j * Sum;
}

template<typename type_S_obj>
inline type_S_obj S_obj(int N_i, double** tensor, double** tensor_reverse, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), int k, int l, int i1, int i2, int j1, int j2) {
	double a = 0, b = 0;
	double c = 0, d = 0;

	a = A_obj + i1 * h1_obj; b = k == 1 ? a + 2.0 * h1_obj : a + h1_obj;
	c = C_obj + i2 * h2_obj; d = k == 2 ? c + 2.0 * h2_obj : c + h2_obj;


	double h_i = (b - a) / N_i,
		h_j = (d - c) / N_i;

	type_S_obj Sum = 0;
	type_S_obj* grad_res, * A_v_res, * Sum_ker = createv<type_S_obj>(3);
	double* v;
	for (size_t i = 0; i < N_i; i++) {
		for (size_t j = 0; j < N_i; j++) {
			double l1 = a + (i + 0.5) * h_i,
				l2 = c + (j + 0.5) * h_j;

			grad_res = createv<type_S_obj>(3, true);
			A_v_res = createv<type_S_obj>(3);
			v = createv<double>(3);

			grad(tensor, tensor_reverse, (*function_x1), (*function_x2), (*function_x3), l1, l2, j1, j2, l, grad_res, true);

			A_v(N_i, 0, k_c, (*function_x1), (*function_x2), (*function_x3), tensor, A_v_res, l1, l2, j1, j2, l, true);
			base_func((*function_x1), (*function_x2), (*function_x3), l1, l2, i1, i2, 0, k, v, true);

			Sum_ker[0] = grad_res[0] + k0 * k0 * A_v_res[0];
			Sum_ker[1] = grad_res[1] + k0 * k0 * A_v_res[1];
			Sum_ker[2] = grad_res[2] + k0 * k0 * A_v_res[2];

			Sum += multv3<type_S_obj>(Sum_ker, v) * sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, tensor));


			del(v); del(A_v_res); del(grad_res);
		}
	}
	return h_i * h_j * Sum;
}

template<typename type_f>
inline type_f f_vec(int N_i,void(*function)(double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), double, double, complex<double>*), double** tensor, double(*function_x1)(double, double, double*, int), double(*function_x2)(double, double, double*, int), double(*function_x3)(double, double, double*, int), int k, int i1, int i2, bool object) {
	double a = 0, b = 0;
	double c = 0, d = 0;

	if (object) {
		a = A_obj + i1 * h1_obj; b = k == 1 ? a + 2.0 * h1_obj : a + h1_obj;
		c = C_obj + i2 * h2_obj; d = k == 2 ? c + 2.0 * h2_obj : c + h2_obj;
	}
	else {
		a = A_screen + i1 * h1_screen; b = k == 1 ? a + 2.0 * h1_screen : a + h1_screen;
		c = C_screen + i2 * h2_screen; d = k == 2 ? c + 2.0 * h2_screen : c + h2_screen;
	}

	double h_i = (b - a) / N_i,
		   h_j = (d - c) / N_i;

	type_f Sum = 0;
	double* v;
	type_f* E0;
	for (size_t i = 0; i < N_i; i++) {
		for (size_t j = 0; j < N_i; j++) {
			double l1 = a + (i + 0.5) * h_i,
				   l2 = c + (j + 0.5) * h_j;
			v = createv<double>(3);
			E0 = createv<type_f>(3);
			base_func((*function_x1), (*function_x2), (*function_x3), l1, l2, i1, i2, 0, k, v, object);
			(*function)((*function_x1), (*function_x2), (*function_x3), l1, l2, E0);

			//cout << multv1<type_f>(v, E0) << "\n";
			
			Sum += multv3<type_f>(E0, v) * sqrt(det_g((*function_x1), (*function_x2), (*function_x3), l1, l2, tensor));
			
			del(v); del(E0);
		}
	}
	return h_i * h_j * Sum;
}

#endif DIFF_GEOM_H