#pragma once
#ifndef DIFF_GEOM_H
#define DIFF_GEOM_H
#include "matrix.h"
#include "integrals.h"
#include "constants.h"

inline void base_func_rft(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...);
inline void base_func(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...);
inline complex<double> k_c(double x1, ...);


inline double x1_obj(int k, ...)
{
	double t1, t2, t3, * var;

	va_list args;
	va_start(args, k);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	var = va_arg(args, double*);
	va_end(args);

	switch (k) {
		//case 0: return x1_start_obj + R * cos(t1) * cos(t2);
	case 0: return x1_start_obj + R * cos(t1) * sin(t2);
	case 1:
	{
		var[0] = (x1_obj(0, t1 + mach_eps, t2, t3, NULL) - x1_obj(0, t1 - mach_eps, t2, t3, NULL)) / 2.0 / mach_eps;
		var[1] = (x1_obj(0, t1, t2 + mach_eps, t3, NULL) - x1_obj(0, t1, t2 - mach_eps, t3, NULL)) / 2.0 / mach_eps;
		var[2] = (x1_obj(0, t1, t2, t3 + mach_eps, NULL) - x1_obj(0, t1, t2, t3 - mach_eps, NULL)) / 2.0 / mach_eps;

		return 0;
	}
	}
}

inline double x2_obj(int k, ...)
{
	double t1, t2, t3, * var;

	va_list args;
	va_start(args, k);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	var = va_arg(args, double*);
	va_end(args);

	switch (k) {
		//case 0:return x2_start_obj + R * cos(t1) * sin(t2);
	case 0: return x2_start_obj + R * sin(t1) * sin(t2);
	case 1:
	{
		var[0] = (x2_obj(0, t1 + mach_eps, t2, t3, NULL) - x2_obj(0, t1 - mach_eps, t2, t3, NULL)) / 2.0 / mach_eps;
		var[1] = (x2_obj(0, t1, t2 + mach_eps, t3, NULL) - x2_obj(0, t1, t2 - mach_eps, t3, NULL)) / 2.0 / mach_eps;
		var[2] = (x2_obj(0, t1, t2, t3 + mach_eps, NULL) - x2_obj(0, t1, t2, t3 - mach_eps, NULL)) / 2.0 / mach_eps;

		return 0;
	}
	}
}

inline double x3_obj(int k, ...)
{
	double t1, t2, t3, * var;

	va_list args;
	va_start(args, k);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	var = va_arg(args, double*);
	va_end(args);

	switch (k) {
		//case 0:return x3_start_obj + R * sin(t1);
	case 0:return x3_start_obj + R * cos(t2);
	case 1:
	{
		var[0] = (x3_obj(0, t1 + mach_eps, t2, t3, NULL) - x3_obj(0, t1 - mach_eps, t2, t3, NULL)) / 2.0 / mach_eps;
		var[1] = (x3_obj(0, t1, t2 + mach_eps, t3, NULL) - x3_obj(0, t1, t2 - mach_eps, t3, NULL)) / 2.0 / mach_eps;
		var[2] = (x3_obj(0, t1, t2, t3 + mach_eps, NULL) - x3_obj(0, t1, t2, t3 - mach_eps, NULL)) / 2.0 / mach_eps;

		return 0;
	}
	}
}




inline double x1_screen(int type, ...)
{
	double t1, t2, t3, * var;

	va_list args;
	va_start(args, type);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	var = va_arg(args, double*);
	va_end(args);

	switch (type) {
	case 0: return x1_start_screen + R * cos(t1) * sin(t2);
	case 1:
	{
		var[0] = (x1_screen(0, t1 + mach_eps, t2, t3, NULL) - x1_screen(0, t1 - mach_eps, t2, t3, NULL)) / 2.0 / mach_eps;
		var[1] = (x1_screen(0, t1, t2 + mach_eps, t3, NULL) - x1_screen(0, t1, t2 - mach_eps, t3, NULL)) / 2.0 / mach_eps;
		var[2] = (x1_screen(0, t1, t2, t3 + mach_eps, NULL) - x1_screen(0, t1, t2, t3 - mach_eps, NULL)) / 2.0 / mach_eps;

		return 0;
	}
	}
}

inline double x2_screen(int type, ...)
{
	double t1, t2, t3, * var;

	va_list args;
	va_start(args, type);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	var = va_arg(args, double*);
	va_end(args);

	switch (type) {
	case 0: return x2_start_screen + R * sin(t1) * sin(t2);
	case 1:
	{
		var[0] = (x2_screen(0, t1 + mach_eps, t2, t3, NULL) - x2_screen(0, t1 - mach_eps, t2, t3, NULL)) / 2.0 / mach_eps;
		var[1] = (x2_screen(0, t1, t2 + mach_eps, t3, NULL) - x2_screen(0, t1, t2 - mach_eps, t3, NULL)) / 2.0 / mach_eps;
		var[2] = (x2_screen(0, t1, t2, t3 + mach_eps, NULL) - x2_screen(0, t1, t2, t3 - mach_eps, NULL)) / 2.0 / mach_eps;

		return 0;
	}
	}
}

inline double x3_screen(int type, ...)
{
	double t1, t2, t3, * var;

	va_list args;
	va_start(args, type);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	var = va_arg(args, double*);
	va_end(args);

	switch (type) {
	case 0: return x3_start_screen + R * cos(t2);
	case 1:
	{

		var[0] = (x3_screen(0, t1 + mach_eps, t2, t3, NULL) - x3_screen(0, t1 - mach_eps, t2, t3, NULL)) / 2.0 / mach_eps;
		var[1] = (x3_screen(0, t1, t2 + mach_eps, t3, NULL) - x3_screen(0, t1, t2 - mach_eps, t3, NULL)) / 2.0 / mach_eps;
		var[2] = (x3_screen(0, t1, t2, t3 + mach_eps, NULL) - x3_screen(0, t1, t2, t3 - mach_eps, NULL)) / 2.0 / mach_eps;

		return 0;
	}
	}
}



void vector_n(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...)
{
	double t1, t2, t3, * var;
	va_list args;
	va_start(args, func_x3);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	var = va_arg(args, double*);
	va_end(args);

	double* ksi_x1 = createv<double>(3), * ksi_x2 = createv<double>(3), * ksi_x3 = createv<double>(3);
	(*func_x1)(1, t1, t2, t3, ksi_x1); (*func_x2)(1, t1, t2, t3, ksi_x2); (*func_x3)(1, t1, t2, t3, ksi_x3);
	
	var[0] = (ksi_x2[0] * ksi_x3[1] - ksi_x3[0] * ksi_x2[1]) / sqrt(pow(ksi_x2[0] * ksi_x3[1] - ksi_x3[0] * ksi_x2[1], 2) + pow(ksi_x3[0] * ksi_x1[1] - ksi_x1[0] * ksi_x3[1], 2) + pow(ksi_x1[0] * ksi_x2[1] - ksi_x2[0] * ksi_x1[1], 2));
	var[1] = (ksi_x3[0] * ksi_x1[1] - ksi_x1[0] * ksi_x3[1]) / sqrt(pow(ksi_x2[0] * ksi_x3[1] - ksi_x3[0] * ksi_x2[1], 2) + pow(ksi_x3[0] * ksi_x1[1] - ksi_x1[0] * ksi_x3[1], 2) + pow(ksi_x1[0] * ksi_x2[1] - ksi_x2[0] * ksi_x1[1], 2));
	var[2] = (ksi_x1[0] * ksi_x2[1] - ksi_x2[0] * ksi_x1[1]) / sqrt(pow(ksi_x2[0] * ksi_x3[1] - ksi_x3[0] * ksi_x2[1], 2) + pow(ksi_x3[0] * ksi_x1[1] - ksi_x1[0] * ksi_x3[1], 2) + pow(ksi_x1[0] * ksi_x2[1] - ksi_x2[0] * ksi_x1[1], 2));
	del(ksi_x1); del(ksi_x2); del(ksi_x3);
}

double det_g(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...)
{
	double t1, t2, t3, ** var;
	va_list args;
	va_start(args, func_x3);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	var = va_arg(args, double**);
	va_end(args);
	
	double* ksi_x1 = createv<double>(3), * ksi_x2 = createv<double>(3), * ksi_x3 = createv<double>(3);
	(*func_x1)(1, t1, t2, t3, ksi_x1); (*func_x2)(1, t1, t2, t3, ksi_x2); (*func_x3)(1, t1, t2, t3, ksi_x3);
	
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

inline void Jacobi(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...) {
	double t1, t2, t3, **var;
	va_list args;
	va_start(args, func_x3);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	var = va_arg(args, double**);
	va_end(args);
	
	double* ksi_x1 = createv<double>(3), * ksi_x2 = createv<double>(3), * ksi_x3 = createv<double>(3);
	(*func_x1)(1, t1, t2, t3, ksi_x1); (*func_x2)(1, t1, t2, t3, ksi_x2); (*func_x3)(1, t1, t2, t3, ksi_x3); 
	
	
	var[0][0] = ksi_x1[0]; var[0][1] = ksi_x1[1];// var[0][2] = ksi_x1[2];
	var[1][0] = ksi_x2[0]; var[1][1] = ksi_x2[1];// var[1][2] = ksi_x2[2];
	var[2][0] = ksi_x3[0]; var[2][1] = ksi_x3[1];// var[2][2] = ksi_x3[2];
	
	del(ksi_x1); del(ksi_x2); del(ksi_x3);
	
}

template<typename type_A_v>
inline void A_v(type_A_v(*func_ker)(double, ...), double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...)
{
	short N_i, type, down_index1, down_index2, up_index2;
	double** tensor, t1, t2, t3;
	type_A_v* var;
	bool object;

	va_list args;
	va_start(args, func_x3);
	N_i = va_arg(args, short);
	type = va_arg(args, short);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	down_index1 = va_arg(args, short);
	down_index2 = va_arg(args, short);
	up_index2 = va_arg(args, short);
	object = va_arg(args, bool);
	tensor = va_arg(args, double**);
	var = va_arg(args, type_A_v*);
	va_end(args);



	var[0] = 0.0;
	var[1] = 0.0;
	var[2] = 0.0;

	double e = 0.0, f = 0.0;
	double g = 0.0, l = 0.0;

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

	switch (type) {
	case 0: {
		for (size_t i = 0; i < N_i; i++) {
			for (size_t j = 0; j < N_i; j++) {
				double l1 = e + (i + 0.65) * h_i,
					   l2 = g + (j + 0.65) * h_j;

				base_func_rft((*func_x1), (*func_x2), (*func_x3), l1, l2, t3, down_index1, down_index2, 0, up_index2, object, v);
				type_A_v ker = (*func_ker)((*func_x1)(0, t1, t2, NULL), (*func_x2)(0, t1, t2, NULL), (*func_x3)(0, t1, t2, NULL), (*func_x1)(0, l1, l2, NULL), (*func_x2)(0, l1, l2, NULL), (*func_x3)(0, l1, l2, NULL));
				double sq_det = sqrt(det_g((*func_x1), (*func_x2), (*func_x3), l1, l2, 0, tensor));

				var[0] += ker * v[0] * sq_det;
				var[1] += ker * v[1] * sq_det;
				var[2] += ker * v[2] * sq_det;
			}
		}
		var[0] *= h_i * h_j;
		var[1] *= h_i * h_j;
		var[2] *= h_i * h_j;

		break;
	}
	case 1: {
		type_A_v* var2 = createv<type_A_v>(3), * var1 = createv<type_A_v>(3);
		A_v((*func_ker), (*func_x1), (*func_x2), (*func_x3), N_i, 0, t1 + mach_eps, t2, t3, down_index1, down_index2, up_index2, object, tensor, var2);
		A_v((*func_ker), (*func_x1), (*func_x2), (*func_x3), N_i, 0, t1 - mach_eps, t2, t3, down_index1, down_index2, up_index2, object, tensor, var1);


		var[0] = (var2[0] - var1[0]) / 2.0 / mach_eps;
		var[1] = (var2[1] - var1[1]) / 2.0 / mach_eps;
		var[2] = (var2[2] - var1[2]) / 2.0 / mach_eps;

		del(var2); del(var1);
		break;
	}
	case 2: {
		type_A_v* var2 = createv<type_A_v>(3), * var1 = createv<type_A_v>(3);
		A_v((*func_ker), (*func_x1), (*func_x2), (*func_x3), N_i, 0, t1, t2 + mach_eps, t3, down_index1, down_index2, up_index2, object, tensor, var2);
		A_v((*func_ker), (*func_x1), (*func_x2), (*func_x3), N_i, 0, t1, t2 - mach_eps, t3, down_index1, down_index2, up_index2, object, tensor, var1);


		var[0] = (var2[0] - var1[0]) / 2.0 / mach_eps;
		var[1] = (var2[1] - var1[1]) / 2.0 / mach_eps;
		var[2] = (var2[2] - var1[2]) / 2.0 / mach_eps;

		del(var2); del(var1);
		break;
	}
	case 3: {
		type_A_v* var2 = createv<type_A_v>(3), * var1 = createv<type_A_v>(3);
		A_v((*func_ker), (*func_x1), (*func_x2), (*func_x3), N_i, 0, t1, t2, t3 + mach_eps, down_index1, down_index2, up_index2, object, tensor, var2);
		A_v((*func_ker), (*func_x1), (*func_x2), (*func_x3), N_i, 0, t1, t2, t3 - mach_eps, down_index1, down_index2, up_index2, object, tensor, var1);


		var[0] = (var2[0] - var1[0]) / 2.0 / mach_eps;
		var[1] = (var2[1] - var1[1]) / 2.0 / mach_eps;
		var[2] = (var2[2] - var1[2]) / 2.0 / mach_eps;

		del(var2); del(var1);
		break;
	}
	}
	del(v);
}

inline void d_x(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...) {
	short type;
	double t1, t2, t3, *var;
	
	va_list args;
	va_start(args, func_x3);
	type = va_arg(args, short);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	var = va_arg(args, double*);
	va_end(args);
	
	switch (type) {
	case 1: {
		double* ksi_x1 = createv<double>(3), * ksi_x2 = createv<double>(3), * ksi_x3 = createv<double>(3);
		(*func_x1)(1, t1, t2, t3, ksi_x1); (*func_x2)(1, t1, t2, t3, ksi_x2); (*func_x3)(1, t1, t2, t3, ksi_x3);

		var[0] = ksi_x1[0];
		var[1] = ksi_x2[0];
		var[2] = ksi_x3[0];
		del(ksi_x1); del(ksi_x2); del(ksi_x3);
		break;
	}
	case 2: {
		double* ksi_x1 = createv<double>(3), * ksi_x2 = createv<double>(3), * ksi_x3 = createv<double>(3);
		(*func_x1)(1, t1, t2, t3, ksi_x1); (*func_x2)(1, t1, t2, t3, ksi_x2); (*func_x3)(1, t1, t2, t3, ksi_x3);
		var[0] = ksi_x1[1];
		var[1] = ksi_x2[1];
		var[2] = ksi_x3[1];
		del(ksi_x1); del(ksi_x2); del(ksi_x3);
		break;
	}
	case 3: {
		double* ksi_x1 = createv<double>(3), * ksi_x2 = createv<double>(3), * ksi_x3 = createv<double>(3);
		(*func_x1)(1, t1, t2, t3, ksi_x1); (*func_x2)(1, t1, t2, t3, ksi_x2); (*func_x3)(1, t1, t2, t3, ksi_x3);
		var[0] = ksi_x1[2];
		var[1] = ksi_x2[2];
		var[2] = ksi_x3[2];
		del(ksi_x1); del(ksi_x2); del(ksi_x3);
		break;
	}
	}
}

template<typename type_div>
inline type_div div(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...) {
	short type, down_index1, down_index2, up_index2;
	double** tensor, ** tensor_reverse, t1, t2, t3;
	bool object;
	
	va_list args;
	va_start(args, func_x3);
	type = va_arg(args, short);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	down_index1 = va_arg(args, short);
	down_index2 = va_arg(args, short);
	up_index2 = va_arg(args, short);
	object = va_arg(args, bool);
	tensor = va_arg(args, double**);
	tensor_reverse = va_arg(args, double**);
	va_end(args);

	size_t M = _msize(tensor_reverse) / sizeof(tensor_reverse[0]);
	size_t N = _msize(tensor_reverse[0]) / sizeof(tensor_reverse[0][0]);

	type_div Sum = 0;

	switch (type) {
	case 0: {
		
		det_g((*func_x1), (*func_x2), (*func_x3), t1, t2, t3, tensor);
		det_g_reverse(tensor, tensor_reverse);
		type_div* A_v_res;
		double* x;
		for (size_t mu = 0; mu < M; mu++)
		{
			for (size_t nu = 0; nu < N; nu++)
			{
				A_v_res = createv<type_div>(3);
				x = createv<double>(3);

				A_v(k_c, (*func_x1), (*func_x2), (*func_x3), 4, mu + 1, t1, t2, t3, down_index1, down_index2, up_index2, object, tensor, A_v_res);
				d_x((*func_x1), (*func_x2), (*func_x3), nu + 1, t1, t2, t3, x);

				Sum += multv3(A_v_res, x) * tensor_reverse[mu][nu];

				del(A_v_res); del(x);
			}
		}
		break;
	}
	case 1: {
		double** tensor1 = createm<double>(M, N), ** tensor_reverse1 = createm<double>(M, N);
		double** tensor2 = createm<double>(M, N), ** tensor_reverse2 = createm<double>(M, N);

		Sum = div<type_div>((*func_x1), (*func_x2), (*func_x3), 0, t1 + mach_eps, t2, t3, down_index1, down_index2, up_index2, object, tensor1, tensor_reverse1);
		Sum -= div<type_div>((*func_x1), (*func_x2), (*func_x3), 0, t1 - mach_eps, t2, t3, down_index1, down_index2, up_index2, object, tensor2, tensor_reverse2);

		Sum /= (2.0 * mach_eps);
		del(tensor1); del(tensor_reverse1);
		del(tensor2); del(tensor_reverse2);
		break;
	}
	case 2: {
		double** tensor1 = createm<double>(M, N), ** tensor_reverse1 = createm<double>(M, N);
		double** tensor2 = createm<double>(M, N), ** tensor_reverse2 = createm<double>(M, N);

		Sum = div<type_div>((*func_x1), (*func_x2), (*func_x3), 0, t1, t2 + mach_eps, t3, down_index1, down_index2, up_index2, object, tensor1, tensor_reverse1);
		Sum -= div<type_div>((*func_x1), (*func_x2), (*func_x3), 0, t1, t2 - mach_eps, t3, down_index1, down_index2, up_index2, object, tensor2, tensor_reverse2);

		Sum /= (2.0 * mach_eps);
		del(tensor1); del(tensor_reverse1);
		del(tensor2); del(tensor_reverse2);
		break;
	}
	case 3: {
		double** tensor1 = createm<double>(M, N), ** tensor_reverse1 = createm<double>(M, N);
		double** tensor2 = createm<double>(M, N), ** tensor_reverse2 = createm<double>(M, N);

		Sum = div<type_div>((*func_x1), (*func_x2), (*func_x3), 0, t1, t2, t3 + mach_eps, down_index1, down_index2, up_index2, object, tensor1, tensor_reverse1);
		Sum -= div<type_div>((*func_x1), (*func_x2), (*func_x3), 0, t1, t2, t3 - mach_eps, down_index1, down_index2, up_index2, object, tensor2, tensor_reverse2);

		Sum /= (2.0 * mach_eps);
		del(tensor1); del(tensor_reverse1);
		del(tensor2); del(tensor_reverse2);
		break;
	}
	}
	return Sum;
}

template<typename type_grad>
inline void grad(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...) {
	short down_index1, down_index2, up_index2;
	double** tensor, ** tensor_reverse, t1, t2, t3;
	bool object;
	type_grad* var;

	va_list args;
	va_start(args, func_x3);
	t1 = va_arg(args, double);
	t2 = va_arg(args, double);
	t3 = va_arg(args, double);
	down_index1 = va_arg(args, short);
	down_index2 = va_arg(args, short);
	up_index2 = va_arg(args, short);
	object = va_arg(args, bool);
	tensor = va_arg(args, double**);
	tensor_reverse = va_arg(args, double**);
	var = va_arg(args, type_grad*);
	va_end(args);
	
	size_t M = _msize(tensor_reverse) / sizeof(tensor_reverse[0]);
	size_t N = _msize(tensor_reverse[0]) / sizeof(tensor_reverse[0][0]);
	
	det_g((*func_x1), (*func_x2), (*func_x3), t1, t2, t3, tensor);
	det_g_reverse(tensor, tensor_reverse);
	
	type_grad div_item;
	double* x;
	for (size_t mu = 0; mu < M; mu++)
	{
		
		div_item = div<type_grad>((*func_x1), (*func_x2), (*func_x3), mu + 1, t1, t2, t3, down_index1, down_index2, up_index2, object, tensor, tensor_reverse);
		
		for (size_t nu = 0; nu < N; nu++)
		{
			x = createv<double>(3);
			d_x((*func_x1), (*func_x2), (*func_x3), nu + 1, t1, t2, t3, x);
			
			var[0] += tensor_reverse[mu][nu] * div_item * x[0];
			var[1] += tensor_reverse[mu][nu] * div_item * x[1];
			var[2] += tensor_reverse[mu][nu] * div_item * x[2];

			del(x);
		}	
	}
}

template<typename type_S_screen>
inline type_S_screen S_screen(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...) {
	
	short N_i, _k, _l, i1, i2, j1, j2;
	double** tensor, ** tensor_reverse;

	va_list args;
	va_start(args, func_x3);
	N_i = va_arg(args, short);
	_k = va_arg(args, short);
	_l = va_arg(args, short);
	i1 = va_arg(args, short);
	i2 = va_arg(args, short);
	j1 = va_arg(args, short);
	j2 = va_arg(args, short);
	tensor = va_arg(args, double**);
	tensor_reverse = va_arg(args, double**);
	va_end(args);
	
	
	double a = 0, b = 0;
	double c = 0, d = 0;
	
	a = A_screen + i1 * h1_screen; b = _k == 1 ? a + 2.0 * h1_screen : a + h1_screen;
	c = C_screen + i2 * h2_screen; d = _k == 2 ? c + 2.0 * h2_screen : c + h2_screen;
	

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

			grad((*func_x1), (*func_x2), (*func_x3), l1, l2, 0, j1, j2, _l, false, tensor, tensor_reverse, grad_res);
			A_v(k_c, (*func_x1), (*func_x2), (*func_x3), N_i, 0, l1, l2, 0, j1, j2, _l, false, tensor, A_v_res);
			
			
			base_func_rft((*func_x1), (*func_x2), (*func_x3), l1, l2, 0, i1, i2, 0, _k, false, v);
			vector_n((*func_x1), (*func_x2), (*func_x3), l1, l2, 0, normal);

			Sum_ker[0] = grad_res[0] + k0 * k0 * A_v_res[0];
			Sum_ker[1] = grad_res[1] + k0 * k0 * A_v_res[1];
			Sum_ker[2] = grad_res[2] + k0 * k0 * A_v_res[2];
			
			
			complex<double> tem = multv3(Sum_ker, normal);

			Sum_ker[0] -= normal[0] * tem;
			Sum_ker[1] -= normal[1] * tem;
			Sum_ker[2] -= normal[2] * tem;

			
			Sum += multv3<type_S_screen>(Sum_ker, v) * sqrt(det_g((*func_x1), (*func_x2), (*func_x3), l1, l2, 0, tensor));
			del(v); del(normal); del(A_v_res); del(grad_res);
		}
	}
	return h_i * h_j * Sum;
}

template<typename type_S_obj>
inline type_S_obj S_obj(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...) {
	short N_i, _k, _l, i1, i2, j1, j2;
	double** tensor, ** tensor_reverse;
	
	va_list args;
	va_start(args, func_x3);
	N_i = va_arg(args, short);
	_k = va_arg(args, short);
	_l = va_arg(args, short);
	i1 = va_arg(args, short);
	i2 = va_arg(args, short);
	j1 = va_arg(args, short);
	j2 = va_arg(args, short);
	tensor = va_arg(args, double**);
	tensor_reverse = va_arg(args, double**);
	va_end(args);
	
	double a = 0, b = 0;
	double c = 0, d = 0;

	a = A_obj + i1 * h1_obj; b = _k == 1 ? a + 2.0 * h1_obj : a + h1_obj;
	c = C_obj + i2 * h2_obj; d = _k == 2 ? c + 2.0 * h2_obj : c + h2_obj;


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

			
			grad<type_S_obj>((*func_x1), (*func_x2), (*func_x3), l1, l2, 0, j1, j2, _l, true, tensor, tensor_reverse, grad_res);
			
			A_v(k_c, (*func_x1), (*func_x2), (*func_x3), N_i, 0, l1, l2, 0, j1, j2, _l, true, tensor, A_v_res);
			base_func_rft((*func_x1), (*func_x2), (*func_x3), l1, l2, 0, i1, i2, 0, _k, true, v);

			Sum_ker[0] = grad_res[0] + k0 * k0 * A_v_res[0];
			Sum_ker[1] = grad_res[1] + k0 * k0 * A_v_res[1];
			Sum_ker[2] = grad_res[2] + k0 * k0 * A_v_res[2];

			Sum += multv3<type_S_obj>(Sum_ker, v) * sqrt(det_g((*func_x1), (*func_x2), (*func_x3), l1, l2, 0, tensor));


			del(v); del(A_v_res); del(grad_res);
		}
	}
	return h_i * h_j * Sum;
}

template<typename type_f>
inline type_f f_vec(void(*func_ker)(double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...), double(*func_x1)(int, ...), double(*func_x2)(int, ...), double(*func_x3)(int, ...), ...) {
	short N_i, _k, _l, i1, i2;
	double** tensor;
	bool object;

	va_list args;
	va_start(args, func_x3);
	N_i = va_arg(args, short);
	_k = va_arg(args, short);
	i1 = va_arg(args, short);
	i2 = va_arg(args, short);
	object = va_arg(args, bool);
	tensor = va_arg(args, double**);
	va_end(args);
	
	
	double a = 0, b = 0;
	double c = 0, d = 0;

	if (object) {
		a = A_obj + i1 * h1_obj; b = _k == 1 ? a + 2.0 * h1_obj : a + h1_obj;
		c = C_obj + i2 * h2_obj; d = _k == 2 ? c + 2.0 * h2_obj : c + h2_obj;
	}
	else {
		a = A_screen + i1 * h1_screen; b = _k == 1 ? a + 2.0 * h1_screen : a + h1_screen;
		c = C_screen + i2 * h2_screen; d = _k == 2 ? c + 2.0 * h2_screen : c + h2_screen;
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
			
			//base_func((*func_x1), (*func_x2), (*func_x3), l1, l2, 0, i1, i2, 0, _k, v, object);
			//base_func_rft((*func_x1), (*func_x2), (*func_x3), l1, l2, 0, i1, i2, 0, _k, v, object);
			(*func_ker)((*func_x1), (*func_x2), (*func_x3), l1, l2, 0, E0);
			
			Sum += multv3<type_f>(E0, v) * sqrt(det_g((*func_x1), (*func_x2), (*func_x3), l1, l2, 0, tensor));
			
			del(v); del(E0);
		}
	}
	return h_i * h_j * Sum;
}

#endif DIFF_GEOM_H