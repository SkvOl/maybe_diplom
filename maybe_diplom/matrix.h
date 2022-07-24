#pragma once


#ifndef MATRIX_H
#define MATRIX_H
#include <Windows.h>
#include <iostream>
#include <malloc.h>
#include <fstream>
#include <math.h>

using namespace std;


double** createm(size_t M, size_t N, bool mod = false) {
    //выделение пам€ти под указатель
    //mod - нужноли при создании создавать единичную матрицу
    double** var = (double**)malloc(M * sizeof(double*));
    for (int i = 0; i < M; i++)
        var[i] = (double*)malloc(N * sizeof(double));

    if(mod)
        for (size_t i = 0; i < M; i++)
            for (size_t j = 0; j < N; j++)
                if (i == j) var[i][j] = 1;
                else var[i][j] = 0;

    return var;
}

double* createv(size_t N, bool mod = false) {
    //выделение пам€ти под указатель
    //mod - нужноли при создании создавать единичный вектор
    double* var = (double*)malloc(N * sizeof(double));

    if (mod)
        for (size_t i = 0; i < N; i++)
            var[i] = 1;

    return var;
}

inline void print(double** var, string c = "") {
    //вывод указател€ var в консоль
    //c (color) - цвет главной диагонали при выводе на экран: G, g - зелЄный, B, b - синий, R, r - красный, I,i - интенсивнее серого
    size_t M = _msize(var) / sizeof(var[0]);
    size_t N = _msize(var[0]) / sizeof(var[0][0]);

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

                printf("%4.*f ", 3, var[i][j]);
                fflush(stdout);
                SetConsoleTextAttribute(hConsoleHandle, 15);
            }
            else {
                printf("%4.*f ", 3, var[i][j]);
                fflush(stdout);
            }

        }
        printf("\n");
        fflush(stdout);
    }
    printf("\n");
    fflush(stdout);
}

inline void print(double* var) {
    //вывод указател€ var в консоль
    size_t M = _msize(var) / sizeof(var[0]);
    
    for (size_t i = 0; i < M; i++)
    {
        printf("%f\n", var[i]);
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

int size(double **var) {
    //вычисл€ет объЄм занимаемой пам€ти указателем var
    size_t M = _msize(var) / sizeof(var[0]);
    size_t sum = 0;
    for (size_t i = 0; i < M; i++)
        sum += _msize(var[i]);
    return sum;
}

int size(double* var) {
    //вычисл€ет объЄм занимаемой пам€ти указателем var
    return _msize(var);
}

string gm(double**& var) {
    //метод √аусса
    size_t M = _msize(var) / sizeof(var[0]);
    size_t N = _msize(var[0]) / sizeof(var[0][0]);

    if (M + 1 != N) return "Error";

    for (size_t k = 0; k < M; k++) {
        if (k == 0) printf("\nk=%i\n", k);
        if (k != 0) printf("k=%i\n", k);
        fflush(stdout);
        double ed = 1;
        if (var[k][k] != ed) {
            double T = var[k][k];
            for (size_t j = k; j < N; j++) {
                var[k][j] = var[k][j] / T;
            }
        }
        for (size_t i = k; i < M; i++) {
            if ((var[i][k] != ed) && (i != k)) {
                double T = var[i][k];
                var[i][k] = 0;
                for (size_t j = k + 1; j < N; j++) {
                    var[i][j] -= var[k][j] * T;
                }
            }
        }
    }
    for (int i = M - 1; i >= 0; i--) {
        double Sum = var[i][M];
        for (size_t j = i + 1.0; j < M; j++) {
            Sum -= var[i][j] * var[j][M];
        }
        var[i][M] = Sum;
    }

    return "Successfully";
}

void del(double**& var) {
    //очистка пам€ти указател€ var
    size_t M = _msize(var) / sizeof(var[0]);
    
    for (size_t i = 0; i < M; i++) free(var[i]);
    free(var);
}

void del(double*& var) {
    //очистка пам€ти указател€ var
    free(var);
}

#endif MATRIX_H