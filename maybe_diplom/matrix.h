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

using namespace std;


template<class type_matrix>
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

template<class type_vector>
inline type_vector* createv(size_t N, bool mod = false) {
    //выделение пам€ти под указатель
    //mod - нужноли при создании создавать единичный вектор
    type_vector* var = (type_vector*)malloc(N * sizeof(type_vector));

    if (mod)
        for (size_t i = 0; i < N; i++)
            var[i] = 1;

    return var;
}

template<class type_matrix_print>
inline void print(type_matrix_print** var, string c = "") {
    //вывод указател€ var в консоль
    //c (color) - цвет главной диагонали при выводе на экран: G, g - зелЄный, B, b - синий, R, r - красный, I,i - интенсивнее серого
    size_t M = _msize(var) / sizeof(var[0]);
    size_t N = _msize(var[0]) / sizeof(var[0][0]);

    const char* type = "";
    if (sizeof(type_matrix_print) == sizeof(int))  type = "%d\n";
    if (sizeof(type_matrix_print) == sizeof(double))  type = "%f\n";
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

template<class type_vector_print>
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

template<class type_matrix_size>
int size(type_matrix_size** var) {
    //вычисл€ет объЄм занимаемой пам€ти указателем var
    size_t M = _msize(var) / sizeof(var[0]);
    size_t sum = 0;
    for (size_t i = 0; i < M; i++)
        sum += _msize(var[i]);
    return sum;
}

template<class type_vector_size>
int size(type_vector_size* var) {
    //вычисл€ет объЄм занимаемой пам€ти указателем var
    return _msize(var);
}

template<class type_gm>
string gm(type_gm**& var) {
    //метод √аусса
    size_t M = _msize(var) / sizeof(var[0]);
    size_t N = _msize(var[0]) / sizeof(var[0][0]);

    if (M + 1 != N) return "Error";

    for (size_t k = 0; k < M; k++) {
        if (k == 0) printf("\nk=%i\n", k);
        if (k != 0) printf("k=%i\n", k);
        fflush(stdout);
        type_gm ed = 1;
        if (var[k][k] != ed) {
            type_gm T = var[k][k];
            for (size_t j = k; j < N; j++) {
                var[k][j] = var[k][j] / T;
            }
        }
        for (size_t i = k; i < M; i++) {
            if ((var[i][k] != ed) && (i != k)) {
                type_gm T = var[i][k];
                var[i][k] = 0;
                for (size_t j = k + 1; j < N; j++) {
                    var[i][j] -= var[k][j] * T;
                }
            }
        }
    }
    for (int i = M - 1; i >= 0; i--) {
        type_gm Sum = var[i][M];
        for (size_t j = i + 1.0; j < M; j++) {
            Sum -= var[i][j] * var[j][M];
        }
        var[i][M] = Sum;
    }

    return "Successfully";
}

template<class type_matrix_del>
void del(type_matrix_del**& var) {
    //очистка пам€ти указател€ var
    size_t M = _msize(var) / sizeof(var[0]);

    for (size_t i = 0; i < M; i++) free(var[i]);
    free(var);
}

template<class type_vector_del>
void del(double*& var) {
    //очистка пам€ти указател€ var
    free(var);
}

#endif MATRIX_H