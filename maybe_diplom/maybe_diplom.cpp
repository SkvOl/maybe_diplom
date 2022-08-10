#include "matrix.h"
#include "integrals.h"


const double pi = 3.1415926, Eps = 0.0001;
const double  k0 = 1, k1 = 1.5 * k0;

const int _n = 10, _N = _n * _n;

double R = 5;
double lambda = 1;

//отрезок для двумерных интегральных уравнений
double A = 0, B = 1;
double C = 0, D = 1;

//шаг для одномерных интегральных уравнений
double h = (B - A) / _n;

//отрезок для двумерных интегральных уравнений
double h1 = (B - A) / _n, h2 = (D - C) / _n;


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
            double ksi = A + (i + 0.5) * h;

            for (size_t j = 0; j < M; j++){
                double a = A + j * h, b = A + (j + 1.0) * h;
                
                var[i][j] = base_func(i, j) * (type - 1.0) - lambda * I_k(100, a, b, ksi);
            }
            var[i][M] = f(ksi);
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
            var[i][M] = f(ksi1, ksi2);
        }
}

void mg(double**& var, size_t type, size_t dim_s = 1, int _step = 0) {
    //метод Галёркина
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    size_t M = _msize(var) / sizeof(var[0]);
    size_t N = _msize(var[0]) / sizeof(var[0][0]);

    int m = (int)sqrt(N - 1);

    if (dim_s == 1) {
        double(*function)(double, ...) = NULL;
        function = &k;

        for (size_t i = 0; i < M; i++) {
            double a = A + (i + _step) * h, b = A + ((i + _step) + 1.0) * h;

            for (size_t j = 0; j < N - 1; j++) {
                double c = A + j * h, d = A + (j + 1.0) * h;

                var[i][j] = h * base_func(i + _step, j) * (type - 1.0) - lambda * I(100, function, a, b, c, d);
            }
            var[i][N - 1] = I(100, a, b);
        }
    }
    else {
        double(*function1)(double, ...) = NULL;
        double(*function2)(double, ...) = NULL;
        function1 = &k;
        function2 = &f;

        for (size_t i = 0; i < N - 1; i++) {
            short i1 = i / m, i2 = i % m;
            double a = A + i1 * h1, b = A + (i1 + 1.0) * h1;
            double c = C + i2 * h2, d = C + (i2 + 1.0) * h2;
            /*printf("i=%d\n", i);
            fflush(stdout);*/

            for (size_t j = 0; j < M; j++) {
                short j1 = (j + _step) / m, j2 = (j + _step) % m;
                double e = A + j1 * h1, f = A + (j1 + 1.0) * h1;
                double g = C + j2 * h2, l = C + (j2 + 1.0) * h2;

                var[j][i] = h1 * h2 * base_func(i, j + _step) * (type - 1.0) - lambda * I(10, function1, a, b, c, d, e, f, g, l);
                var[j][N - 1] = I(100, function2, e, f, g, l);
            } 
        }
    }
}


int main() {
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
    /*printf("rank: %d  time fill matrix is: %f\n", _rank, t2);
    fflush(stdout);*/
    //print(a);

    t1 = MPI_Wtime();
    gm(a, count_one_rank, _step, _rank, _size);
    t2 = MPI_Wtime() - t1;
    /*printf("rank: %d  time Gauss method is: %f\n", _rank, t2);
    fflush(stdout);*/
    
    
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