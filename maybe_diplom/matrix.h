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
    double* var = (double*)malloc(N * sizeof(double));

    if (mod)
        for (size_t i = 0; i < N; i++)
            var[i] = 1;

    return var;
}

inline void print(double** var, string c = "") {
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
    for (size_t ind_k = 0; ind_k < k; ind_k++)
    {
        printf("\n");
        fflush(stdout);
    }
}

int size(double **var) {
    size_t M = _msize(var) / sizeof(var[0]);
    size_t sum = 0;
    for (size_t i = 0; i < M; i++)
        sum += _msize(var[i]);
    return sum;
}

int size(double* var) {
    return _msize(var);
}

void del(double**& var) {
    size_t M = _msize(var) / sizeof(var[0]);
    
    for (size_t i = 0; i < M; i++) free(var[i]);
    free(var);
}

void del(double*& var) {
    free(var);
}
#endif MATRIX_H