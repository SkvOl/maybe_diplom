#pragma once


#ifndef MATRIX_H
#define MATRIX_H
#include <Windows.h>
#include <iostream>
#include <malloc.h>
#include <fstream>
#include <complex>
#include <string>
#include <math.h>
#include <mpi.h>

using namespace std;


template<typename type_matrix>
inline type_matrix** createm(size_t M, size_t N, bool mod = false) {
	//выделение пам€ти под указатель
	//mod - нужноли при создании создавать единичную матрицу

	type_matrix** var = (type_matrix**)malloc(M * sizeof(type_matrix*));
	for (size_t i = 0; i < M; i++)
		var[i] = (type_matrix*)malloc(N * sizeof(type_matrix));



	if (mod)
		for (size_t i = 0; i < M; i++)
			for (size_t j = 0; j < N; j++)
				if (i == j) var[i][j] = 1;
				else var[i][j] = 0;

	return var;
}

template<typename type_vector>
inline type_vector* createv(size_t N, bool mod = false) {
	//выделение пам€ти под указатель
	//mod - нужноли при создании создавать единичный вектор
	type_vector* var = (type_vector*)malloc(N * sizeof(type_vector));

	if (mod)
		for (size_t i = 0; i < N; i++)
			var[i] = 1;

	return var;
}

template<typename type_matrix_print>
inline void print(type_matrix_print** var, string c = "", int cM = -1, int cN = -1) {
	//вывод указател€ var в консоль
	//c (color) - цвет главной диагонали при выводе на экран: G, g - зелЄный, B, b - синий, R, r - красный, I,i - интенсивнее серого
	size_t M, N;
	if (cM == -1 && cN == -1) {
		M = _msize(var) / sizeof(var[0]);
		N = _msize(var[0]) / sizeof(var[0][0]);
	}
	else {
		M = cM;
		N = cN;
	}

	const char* type = "";
	if (sizeof(type_matrix_print) == sizeof(int))  type = "%d ";
	if (sizeof(type_matrix_print) == sizeof(double))  type = "%f ";
	if (sizeof(type_matrix_print) == sizeof(complex<double>))  type = "complex";

	HANDLE hConsoleHandle = GetStdHandle(STD_OUTPUT_HANDLE);
	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (i == j) {
				if (c == "g" || c == "G")SetConsoleTextAttribute(hConsoleHandle, FOREGROUND_GREEN);
				if (c == "b" || c == "B")SetConsoleTextAttribute(hConsoleHandle, FOREGROUND_BLUE);
				if (c == "r" || c == "R")SetConsoleTextAttribute(hConsoleHandle, FOREGROUND_RED);
				if (c == "i" || c == "I")SetConsoleTextAttribute(hConsoleHandle, FOREGROUND_INTENSITY);

				if (type != "complex") printf(type, var[i][j]);
				else cout << var[i][j] << " ";

				fflush(stdout);
				SetConsoleTextAttribute(hConsoleHandle, 15);
			}
			else {
				if (type != "complex") printf(type, var[i][j]);
				else cout << var[i][j] << " ";

				fflush(stdout);
			}

		}
		printf("\n");
		fflush(stdout);
	}
	printf("\n");
	fflush(stdout);
}

template<typename type_vector_print>
inline void print(type_vector_print* var) {
	//вывод указател€ var в консоль
	size_t M = _msize(var) / sizeof(var[0]);

	const char* type = "";
	if (sizeof(type_vector_print) == sizeof(int))  type = "%d\n";
	if (sizeof(type_vector_print) == sizeof(double))  type = "%f\n";
	if (sizeof(type_vector_print) == sizeof(complex<double>))  type = "complex";

	for (size_t i = 0; i < M; i++)
	{
		if (type != "complex") printf(type, var[i]);
		else cout << var[i] << "\n";

		fflush(stdout);
	}
	printf("\n");
	fflush(stdout);
}

void space(size_t k = 0) {
	//вывод пробелов в консоль
	for (size_t ind_k = 0; ind_k < k; ind_k++)
	{
		printf("\n");
		fflush(stdout);
	}
}

template<typename type_matrix_size>
int size(type_matrix_size** var) {
	//вычисл€ет объЄм занимаемой пам€ти указателем var
	size_t M = _msize(var) / sizeof(var[0]);
	size_t sum = 0;
	for (size_t i = 0; i < M; i++)
		sum += _msize(var[i]);
	return sum;
}

template<typename type_vector_size>
int size(type_vector_size* var) {
	//вычисл€ет объЄм занимаемой пам€ти указателем var
	return _msize(var);
}

template<class type_gm>
string gm(type_gm**& var, int* count_one_rank = NULL, int _step = 0, int _rank = 0, int _size = 0) {
	//метод √аусса
	MPI_Status status;
	MPI_Request request;
	size_t M = _msize(var) / sizeof(var[0]);
	size_t N = _msize(var[0]) / sizeof(var[0][0]);
	MPI_Datatype type;

	if (sizeof(type_gm) == sizeof(int))  type = MPI_INT;
	if (sizeof(type_gm) == sizeof(double))  type = MPI_DOUBLE;
	if (sizeof(type_gm) == sizeof(complex<double>))  type = MPI_DOUBLE_COMPLEX;


	for (size_t _k = 0; _k < N - 1; _k++) {
		type_gm ed = 1;
		if (_k >= _step && _k < _step + count_one_rank[_rank]) {
			if (var[_k - _step][_k] != ed) {
				type_gm T = var[_k - _step][_k];
				for (size_t j = _k; j < N; j++) {
					var[_k - _step][j] = var[_k - _step][j] / T;
				}
			}

			for (size_t i = _rank + 1; i < _size; i++) {
				MPI_Ssend(var[_k - _step], N, type, i, 1, MPI_COMM_WORLD);
			}

			for (size_t i = _k - _step; i < M; i++) {
				if ((var[i][_k] != ed) && (i != _k - _step)) {
					type_gm T = var[i][_k];
					var[i][_k] = 0;
					for (size_t j = _k + 1; j < N; j++) {
						var[i][j] -= var[_k - _step][j] * T;
					}
				}
			}
		}
		else if (_k < _step) {
			size_t step = 0, head_rank = 0;
			for (head_rank = 0; head_rank < _size; head_rank++) {
				step += count_one_rank[head_rank];
				if (_k < step) break;
			}

			type_gm* _kvar = createv<type_gm>(N);
			MPI_Recv(_kvar, N, type, head_rank, 1, MPI_COMM_WORLD, &status);

			for (size_t i = 0; i < M; i++) {
				type_gm T = var[i][_k];
				var[i][_k] = 0;
				for (size_t j = _k + 1; j < N; j++)
					var[i][j] -= _kvar[j] * T;
			}
		}
	}

	//обратный ход
	int M_in;
	type_gm* f_in = NULL, * f_out = createv<type_gm>(M);
	if (_rank != _size - 1) {
		int step = 0;
		for (size_t _i = 0; _i < _size - 1 - _rank; _i++) {
			MPI_Probe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, type, &M_in);
			f_in = createv<type_gm>(M_in);
			MPI_Recv(f_in, M_in, type, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
			for (int i = M - 1; i >= 0; i--) {
				step = 0;
				for (size_t i = _size - 1; i > status.MPI_SOURCE; i--)
					step += count_one_rank[i];

				for (size_t j = 2 + step; j < M_in + step + 2; j++)
					var[i][N - 1] -= var[i][N - j] * f_in[M_in - j + step + 1];
			}
		}

		step = 0;
		for (size_t i = _size - 1; i > _rank; i--)
			step += count_one_rank[i];

		for (int i = M - 2; i >= 0; i--) {
			for (size_t j = 2 + step; j < N - _step - i; j++)
				var[i][N - 1] -= var[i][N - j] * var[M + 1 - j + step][N - 1];

			f_out[i] = var[i][N - 1];
		}
		f_out[M - 1] = var[M - 1][N - 1];
	}
	else {
		for (int i = M - 2; i >= 0; i--) {
			for (size_t j = 2; j < N - _step - i; j++)
				var[i][N - 1] -= var[i][N - j] * var[M - j + 1][N - 1];

			f_out[i] = var[i][N - 1];
		}
		f_out[M - 1] = var[M - 1][N - 1];
	}

	if (_rank != 0)
		for (int i = _rank - 1; i >= 0; i--)
			MPI_Send(f_out, M, type, i, 2, MPI_COMM_WORLD);


	//обратный ход


	return "rank: " + to_string(_rank) + "  Successfully";
}

template<typename type_fill>
type_fill* fill(int M, type_fill item = 0.0) {
	type_fill* var = createv<type_fill>(M);
	for (size_t i = 0; i < M; i++) {
		var[i] = item;
	}
	return var;
}

template<typename type_simple_iteration>
type_simple_iteration* simple_iteration(type_simple_iteration** var = 0, type_simple_iteration* init_approx = 0, size_t k = 0, int* count_one_rank = NULL, int* arr_step = NULL, int _rank = 0, int _size = 0)
{
	size_t M = _msize(var) / sizeof(var[0]);
	size_t N = _msize(var[0]) / sizeof(var[0][0]);
	MPI_Datatype type;

	if (sizeof(type_simple_iteration) == sizeof(int))  type = MPI_INT;
	if (sizeof(type_simple_iteration) == sizeof(double))  type = MPI_DOUBLE;

	printf("simple iteration M: %d\n", M);

	type_simple_iteration* res = createv<type_simple_iteration>(M);
	size_t k_ind = 0;
	while (true) {
		res = fill<type_simple_iteration>(M);
		k_ind++;
		printf("s_i k: %d", k_ind);
		fflush(stdout);
		for (size_t i = 0; i < M; i++)
		{
			if (var[i][i + arr_step[_rank]] > 0) {
				for (size_t j = 0; j < N - 1; j++)
				{
					if ((i + arr_step[_rank]) == j) res[i] -= (var[i][j] - 1) * init_approx[j];
					if ((i + arr_step[_rank]) != j) res[i] -= var[i][j] * init_approx[j];
				}
				res[i] += var[i][N - 1];
			}
			if (var[i][i + arr_step[_rank]] < 0) {
				for (size_t j = 0; j < N - 1; j++)
				{
					if ((i + arr_step[_rank]) == j) res[i] += (var[i][j] + 1) * init_approx[j];
					if ((i + arr_step[_rank]) != j) res[i] += var[i][j] * init_approx[j];
				}
				res[i] -= var[i][N - 1];
			}
		}
		if (k_ind <= k) {
			for (size_t i = 0; i < _size; i++)
				MPI_Gatherv(res, count_one_rank[_rank], type, init_approx, count_one_rank, arr_step, type, i, MPI_COMM_WORLD);
			
			res = fill<type_simple_iteration>(M);
		}
		else {
			return res;
		}
	}
}

complex<double>* simple_iteration_c(complex<double>** var = 0, complex<double>* init_approx = 0, size_t k = 0, int* count_one_rank = NULL, int* arr_step = NULL, int _rank = 0, int _size = 0)
{
	//¬озможно метод вообще не подходит дл€ комплексных, проверить его правильность не удалось.
	size_t M = _msize(var) / sizeof(var[0]);
	size_t N = _msize(var[0]) / sizeof(var[0][0]);
	complex<double>* res = createv<complex<double>>(M), one = 1.0;
	size_t k_ind = 0;
	while (true) {
		res = fill<complex<double>>(M);
		k_ind++;
		//printf("s_i k: %d\n", k_ind);
		fflush(stdout);
		for (size_t i = 0; i < M; i++)
		{
			if (var[i][i + arr_step[_rank]].real() > 0) {
				for (size_t j = 0; j < N - 1; j++)
				{
					if ((i + arr_step[_rank]) == j) res[i] -= (var[i][j] - one) * init_approx[j];
					if ((i + arr_step[_rank]) != j) res[i] -= var[i][j] * init_approx[j];
				}
				res[i] += var[i][N - 1];
			}
			if (var[i][i + arr_step[_rank]].real() < 0) {
				for (size_t j = 0; j < N - 1; j++)
				{
					if ((i + arr_step[_rank]) == j) res[i] += (var[i][j] + one) * init_approx[j];
					if ((i + arr_step[_rank]) != j) res[i] += var[i][j] * init_approx[j];
				}
				res[i] -= var[i][N - 1];
			}
		}
		if (k_ind <= k) {
			for (size_t i = 0; i < _size; i++)
				MPI_Gatherv(res, count_one_rank[_rank], MPI_DOUBLE_COMPLEX, init_approx, count_one_rank, arr_step, MPI_DOUBLE_COMPLEX, i, MPI_COMM_WORLD);

			res = fill<complex<double>>(M);
		}
		else {
			return res;
		}
	}
}

template<typename type_matrix_del>
void del(type_matrix_del**& var) {
	//очистка пам€ти указател€ var
	size_t M = _msize(var) / sizeof(var[0]);

	for (size_t i = 0; i < M; i++) free(var[i]);
	free(var);
}

template<typename type_vector_del>
void del(double*& var) {
	//очистка пам€ти указател€ var
	free(var);
}

template<typename type_LU>
void LU(type_LU** var, type_LU**& L, type_LU**& U) {
	size_t M = _msize(var) / sizeof(var[0]);

	for (size_t i = 0; i < M; i++) {
		for (size_t j = 0; j < M; j++) {
			if (i == j) L[i][i] = 1;
			else L[i][j] = 0;
			U[i][j] = 0;
		}
	}
	for (short i = 0; i < M; i++) {
		for (short j = 0; j < M; j++) {
			type_LU sumU = 0, sumL = 0;
			if (i <= j) {
				for (short z = 0; z <= i - 1; z++) {
					sumU += L[i][z] * U[z][j];
				}
				U[i][j] = var[i][j] - sumU;
			}

			if (i > j) {
				for (short z = 0; z <= j - 1; z++) {
					sumL += L[i][z] * U[z][j];
				}
				L[i][j] = (var[i][j] - sumL) / U[j][j];
			}
		}
	}
}

template<typename type_mult>
inline void mult(type_mult** var1, type_mult** var2, type_mult**& res) {
	size_t M = _msize(var1) / sizeof(var1[0]);
	size_t N = _msize(var1[0]) / sizeof(var1[0][0]);

	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
			type_mult sum = 0;
			for (size_t k = 0; k < M; k++)
			{
				sum += var1[i][k] * var2[k][j];
			}
			res[i][j] = sum;
		}
	}
}

template<typename type_diag>
type_diag** diag(type_diag** var) {
	size_t M = _msize(var) / sizeof(var[0]);


	type_diag** var_D = createm<type_diag>(M, M), ** var_L = createm<type_diag>(M, M), ** var_U = createm<type_diag>(M, M);

	LU(var, var_L, var_U);
	mult(var_U, var_L, var_D);

	for (size_t k = 1; k <= 10; k++) {
		LU(var_D, var_L, var_U);
		mult(var_U, var_L, var_D);
	}
	del(var_U); del(var_L);

	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < M; j++)
		{
			if (i != j) var_D[i][j] = 0;
		}
	}
	return var_D;
}

template<typename type_eigenvalues>
type_eigenvalues* eigenvalues(type_eigenvalues** var) {
	size_t M = _msize(var) / sizeof(var[0]);

	type_eigenvalues** var_D = diag<type_eigenvalues>(var), * res = createv<type_eigenvalues>(M);
	for (size_t i = 0; i < M; i++)
		res[i] = var_D[i][i];
	
	del(var_D);
	return res;
}

template<typename type_cond>
type_cond cond(type_cond** var) {
	size_t M = _msize(var) / sizeof(var[0]);


	type_cond* var1 = eigenvalues(var);

	type_cond max = 0, min = 0;
	for (size_t i = 0; i < M; i++)
	{
		if (abs(var1[i]) >= max) max = abs(var1[i]);
		if (abs(var1[i]) <= max) min = abs(var1[i]);
	}

	del(var1);
	return max / min;
}

#endif MATRIX_H