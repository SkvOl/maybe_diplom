#include "matrix.h"
#include "integrals.h"



int _n = 2, _N = _n * _n;


//отрезок для n мерных интегральных уравнений
double A = 0, B = 1;
double C = 0, D = 1;
double E = 0, F = 1;


//шаг  для n мерных интегральных уравнений
double h1 = (B - A) / _n;
double h2 = (D - C) / _n;
double h3 = (F - E) / _n;


inline double u(double y1) {     
    return y1 * y1;
}

void mk(double**& var, size_t type, size_t dim_s = 1) {
    //метод Коллокаций
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    size_t M = _msize(var) / sizeof(var[0]);
    int m = (int) sqrt(M);

    if(dim_s == 1)
        for (size_t i = 0; i < M; i++){
            double ksi = A + (i + 0.5) * h1;

            for (size_t j = 0; j < M; j++){
                double a = A + j * h1, b = A + (j + 1.0) * h1;
                
                var[i][j] = base_func(i, j) * (type - 1.0) - lambda * I_k(100, a, b, ksi);
            }
            var[i][M] = func(ksi);
        }
    else
        for (size_t i = 0; i < M; i++){
            short i1 = i / m, i2 = i % m;
            double ksi1 = A + (i1 + 0.5) * h1, ksi2 = A + (i2 + 0.5) * h2;

            for (size_t j = 0; j < M; j++) {
                short j1 = j / m, j2 = j % m;
                double a = A + j1 * h1, b = A + (j1 + 1.0) * h1, c = C + j2 * h2, d = C + (j2 + 1.0) * h2;

                var[i][j] = base_func(i, j) * (type - 1.0) - lambda * I_k(100, a, b, c, d, ksi1, ksi2);     
            }
            var[i][M] = func(ksi1, ksi2);
        }
}

void mg(double**& var, size_t type, size_t dim_s = 1, int _step = 0) {
    //метод Галёркина
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    size_t M = _msize(var) / sizeof(var[0]);
    size_t N = _msize(var[0]) / sizeof(var[0][0]);

    if (dim_s == 1) {
        for (size_t i = 0; i < M; i++) {
            double a = A + (i + _step) * h1, b = A + ((i + _step) + 1.0) * h1;

            for (size_t j = 0; j < N; j++) {
                double c = A + j * h1, d = A + (j + 1.0) * h1;

                var[i][j] = h1 * base_func((i + _step), j) * (type - 1.0) - lambda * In<double>(3, k, a, b, c, d);
            }
            var[i][N - 1] = In<double>(3, a, b);
        }
    }
    else if (dim_s == 2) {
        int m = (int)sqrt(N - 1);
        for (size_t i = 0; i < M; i++) {
            short i1 = (i + _step) / m, i2 = (i + _step) % m;
            double a = A + i1 * h1, b = A + (i1 + 1.0) * h1;
            double c = C + i2 * h2, d = C + (i2 + 1.0) * h2;

            //printf("i=%d\n", i);
            for (size_t j = 0; j < N; j++) {
                short j1 = j / m, j2 = j % m;
                double e = A + j1 * h1, f = A + (j1 + 1.0) * h1;
                double g = C + j2 * h2, l = C + (j2 + 1.0) * h2;

                var[i][j] = h1 * h2 * base_func(i + _step, j) * (type - 1.0) - lambda * In<double>(3, k, a, b, c, d, e, f, g, l);
            }
            var[i][N - 1] = In<double>(3, func, a, b, c, d);
        }
    }
    else if (dim_s == 3) {
        int m = (int)pow(N - 1, 1.0 / 3.0);

        for (size_t i = 0; i < M; i++) {
            short i1 = (i + _step) % m, i2 = (i + _step) / m - m * ((i + _step) / m / m), i3 = (i + _step) / m / m;
            double a = A + i1 * h1, b = a + h1;
            double c = C + i2 * h2, d = c + h2;
            double e = E + i3 * h3, f = e + h3;
        
            for (size_t j = 0; j < N; j++) {
                short j1 = j % m, j2 = j / m - m * (j / m / m), j3 = j / m / m;

                double g = A + j1 * h1, l = g + h1;
                double o = C + j2 * h2, p = o + h2;
                double q = E + j3 * h3, s = q + h3;

                var[i][j] = h1 * h2 * h3 * base_func(i + _step, j) * (type - 1.0) - lambda * In<double>(5, k, a, b, c, d, e, f, g, l, o, p, q, s);                       
            }
            var[i][N - 1] = In<double>(5, func, a, b, c, d, e, f);
        }
    }
}


int main1() {
    return 0;
    int _rank, _size;
    
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Status status;
    MPI_Request request;

    /// подготовительные работы
    int* count_one_rank = createv<int>(_size);

    for (size_t i = 0; i < _size; i++)
        count_one_rank[i] = _N / _size;

    short mod = _N % _size;
    if(mod != 0)
        for (short i = _size - 1; i >= 0; i--) {
            count_one_rank[i]++;
            mod--;
            if (mod == 0) break;
        }

    int _step = 0;
    for (size_t i = 0; i < _rank; i++)
        _step += count_one_rank[i];
    /// подготовительные работы

    double** a = createm<double>(count_one_rank[_rank], _N + 1);

    double t1 = MPI_Wtime();
    mg(a, 2, 2, _step);
    double t2 = MPI_Wtime() - t1;
    printf("rank: %d  time fill matrix is: %f\n", _rank, t2);
    fflush(stdout);
    //if (_rank == 0)print(a);

    t1 = MPI_Wtime();
    gm(a, count_one_rank, _step, _rank, _size);
    t2 = MPI_Wtime() - t1;
    printf("rank: %d  time Gauss method is: %f\n", _rank, t2);
    fflush(stdout);
    
    
    double* part_res = createv<double>(count_one_rank[_rank]), * res = NULL;
    int* arr_step = NULL;

    
    for (size_t i = 0; i < count_one_rank[_rank]; i++)
        part_res[i] = a[i][_N];
    
    
    if (_rank == 0) {
        del(a);

        res = createv<double>(_N);
        arr_step = createv<int>(_size);

        for (size_t i = 0; i < _size; i++) {
            arr_step[i] = 0;

            for (size_t j = 0; j < i; j++)
                arr_step[i] += count_one_rank[j];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(part_res, count_one_rank[_rank], MPI_DOUBLE, res, count_one_rank, arr_step, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (_rank == 0) 
        for (size_t i = 0; i < _N; i++) {
            if (i % _n == 0 && i != 0) printf("\n");
            printf("%f ", res[i]);
        }  


    MPI_Finalize();
}

int main() {
    //return 0;
    _N *= _n;
    int _rank, _size;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Status status;
    MPI_Request request;

    /// подготовительные работы
    int* count_one_rank = createv<int>(_size);

    for (size_t i = 0; i < _size; i++)
        count_one_rank[i] = _N / _size;

    short mod = _N % _size;
    if (mod != 0)
        for (short i = _size - 1; i >= 0; i--) {
            count_one_rank[i]++;
            mod--;
            if (mod == 0) break;
        }

    int _step = 0;
    for (size_t i = 0; i < _rank; i++)
        _step += count_one_rank[i];
    /// подготовительные работы

    double** a = createm<double>(count_one_rank[_rank], _N + 1);

    double t1 = MPI_Wtime();
    mg(a, 2, 3, _step);
    double t2 = MPI_Wtime() - t1;
    printf("rank: %d  time fill matrix is: %f\n", _rank, t2);
    fflush(stdout);
    //if (_rank == 0) print(a);

    t1 = MPI_Wtime();
    gm(a, count_one_rank, _step, _rank, _size);
    t2 = MPI_Wtime() - t1;
    printf("rank: %d  time Gauss method is: %f\n", _rank, t2);
    fflush(stdout);


    double* part_res = createv<double>(count_one_rank[_rank]), * res = NULL;
    int* arr_step = NULL;


    for (size_t i = 0; i < count_one_rank[_rank]; i++)
        part_res[i] = a[i][_N];


    if (_rank == 0) {
        del(a);

        res = createv<double>(_N);
        arr_step = createv<int>(_size);

        for (size_t i = 0; i < _size; i++) {
            arr_step[i] = 0;

            for (size_t j = 0; j < i; j++)
                arr_step[i] += count_one_rank[j];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(part_res, count_one_rank[_rank], MPI_DOUBLE, res, count_one_rank, arr_step, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (_rank == 0)
        for (size_t i = 0; i < _N; i++) {
            printf("%f\n", res[i]);
        }


    MPI_Finalize();
}