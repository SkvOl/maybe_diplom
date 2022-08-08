#pragma once


#ifndef MATRIX_H
#define MATRIX_H
#include <Windows.h>
#include <iostream>
#include <malloc.h>
#include <fstream>
#include <string>
#include <math.h>
#include <mpi.h>

using namespace std;


double** createm(size_t M, size_t N, bool mod = false) {
    //выделение памяти под указатель
    //mod - нужноли при создании создавать единичную матрицу

    double** var = (double**)malloc(M * sizeof(double*));
    for (int i = 0; i < M; i++)
        var[i] = (double*)malloc(N * sizeof(double));
    


    if (mod)
        for (size_t i = 0; i < M; i++)
            for (size_t j = 0; j < N; j++)
                if (i == j) var[i][j] = 1;
                else var[i][j] = 0;

    return var;
}

double* createv(size_t N, bool mod = false) {
    //выделение памяти под указатель
    //mod - нужноли при создании создавать единичный вектор
    double* var = (double*)malloc(N * sizeof(double));

    if (mod)
        for (size_t i = 0; i < N; i++)
            var[i] = 1;

    return var;
}

inline void print(double** var, string c = "") {
    //вывод указателя var в консоль
    //c (color) - цвет главной диагонали при выводе на экран: G, g - зелёный, B, b - синий, R, r - красный, I,i - интенсивнее серого
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
    //вывод указателя var в консоль
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
    //вычисляет объём занимаемой памяти указателем var
    size_t M = _msize(var) / sizeof(var[0]);
    size_t sum = 0;
    for (size_t i = 0; i < M; i++)
        sum += _msize(var[i]);
    return sum;
}

int size(double* var) {
    //вычисляет объём занимаемой памяти указателем var
    return _msize(var);
}

string gm(double**& var, int* count_one_rank = NULL, int _step = 0, int _rank = 0, int _size = 0) {
    //метод Гаусса
    MPI_Status status;
    MPI_Request request;
    size_t M = _msize(var) / sizeof(var[0]);
    size_t N = _msize(var[0]) / sizeof(var[0][0]);

    for (size_t _k = 0; _k < N - 1; _k++) {
        double ed = 1;
        if(_k >= _step && _k < _step + count_one_rank[_rank]) {
            if (var[_k - _step][_k] != ed) {
                double T = var[_k - _step][_k];
                for (size_t j = _k; j < N; j++) {
                    var[_k - _step][j] = var[_k - _step][j] / T;
                }
            }
            for (size_t i = 0; i < _size; i++) {
                if(i != _rank) MPI_Send(&_rank, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            for (size_t i = _rank + 1; i < _size; i++) {
                MPI_Send(var[_k - _step], N, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            }

            for (size_t i = _k - _step; i < M; i++) {
                if ((var[i][_k] != ed) && (i != _k - _step)) {
                    double T = var[i][_k];
                    var[i][_k] = 0;
                    for (size_t j = _k + 1; j < N; j++) {
                        var[i][j] -= var[_k - _step][j] * T;
                    }
                }
            }
        }
        else{
            int head_rank;
            MPI_Recv(&head_rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            if (_rank > head_rank) {
                double* _kvar = createv(N);
                MPI_Recv(_kvar, N, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
                for (size_t i = 0; i < M; i++) {
                    double T = var[i][_k];
                    var[i][_k] = 0;
                    for (size_t j = _k + 1; j < N; j++) {
                        var[i][j] -= _kvar[j] * T;
                    }
                }
            }
        }
    }
    
    //обратный ход
    int M_in;
    double* f_in = NULL, *f_out = createv(M);
    if (_rank != _size - 1) {
        int step = 0;
        for (size_t _i = 0; _i < _size - 1 - _rank; _i++) {
            MPI_Probe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &M_in);
            f_in = createv(M_in);
            MPI_Recv(f_in, M_in, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
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
            MPI_Send(f_out, M, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
        

    //обратный ход
    

    return "rank: " + to_string(_rank) + "  Successfully";
}

void del(double**& var) {
    //очистка памяти указателя var
    size_t M = _msize(var) / sizeof(var[0]);
    
    for (size_t i = 0; i < M; i++) free(var[i]);
    free(var);
}

void del(double*& var) {
    //очистка памяти указателя var
    free(var);
}

#endif MATRIX_H