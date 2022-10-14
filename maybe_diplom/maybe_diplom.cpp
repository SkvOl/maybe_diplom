#include "matrix.h"
#include "integrals.h"
#include "diff_geom.h"
#include "constants.h"




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



void mg(complex<double>** var, size_t type, size_t dim_s = 1, int _step = 0) {
    //метод Галёркина
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    size_t M = _msize(var) / sizeof(var[0]);
    size_t N = _msize(var[0]) / sizeof(var[0][0]);
    double** tensor, **tensor_reverse;
    short i1, i2, j1, j2, _k, _l;
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
            i1 = (i + _step) / m, i2 = (i + _step) % m;
            _k = i + _step < (N - 1) / 2 ? 1 : 2;
            //double a = A + i1 * h1, b = A + (i1 + 1.0) * h1;
            //double c = C + i2 * h2, d = C + (i2 + 1.0) * h2;

            if (_step == 0) {
                printf("%d\n", i);
                fflush(stdout);
            }
            tensor = createm<double>(2, 2);
            tensor_reverse = createm<double>(2, 2);
            for (size_t j = 0; j < N; j++) {  
                j1 = j / m, j2 = j % m;
                _l = j + _step < (N - 1) / 2 ? 1 : 2;
                
                //double e = A + j1 * h1, f = A + (j1 + 1.0) * h1;
                //double g = C + j2 * h2, l = C + (j2 + 1.0) * h2;

                //var[i][j] = h1 * h2 * base_func(i + _step, j) * (type - 1.0) - lambda * In<double>(3, k, a, b, c, d, e, f, g, l);
                var[i][j] = lambda * S<complex<double>>(1, tensor, tensor_reverse, x1_screen, x2_screen, x3_screen, _k, _l, i1, i2, j1, j2);
                
            }
            //var[i][N - 1] = In<double>(3, func, a, b, c, d);
            var[i][N - 1] = f_vec<complex<double>>(1, tensor, x1_screen, x2_screen, x3_screen, _k, i1, i2);          
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


int main() {
    //return 0;
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

    complex<double>** a = createm<complex<double>>(count_one_rank[_rank], _N + 1);

    double t1 = MPI_Wtime();
    mg(a, 2, 2, _step);
    double t2 = MPI_Wtime() - t1;
    printf("rank: %d  time fill matrix is: %f\n", _rank, t2);
    fflush(stdout);
    if (_rank == 0)print(a, "g", 10, 10);

    t1 = MPI_Wtime();
    gm(a, count_one_rank, _step, _rank, _size);
    t2 = MPI_Wtime() - t1;
    printf("rank: %d  time Gauss method is: %f\n", _rank, t2);
    fflush(stdout);
    
    
    complex<double>* part_res = createv<complex<double>>(count_one_rank[_rank]), * res = NULL;
    int* arr_step = NULL;

    
    for (size_t i = 0; i < count_one_rank[_rank]; i++)
        part_res[i] = a[i][_N];
    
    
    if (_rank == 0) {
        del(a);

        res = createv<complex<double>>(_N);
        arr_step = createv<int>(_size);

        for (size_t i = 0; i < _size; i++) {
            arr_step[i] = 0;

            for (size_t j = 0; j < i; j++)
                arr_step[i] += count_one_rank[j];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(part_res, count_one_rank[_rank], MPI_DOUBLE_COMPLEX, res, count_one_rank, arr_step, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    if (_rank == 0) {
        ofstream file1("res.txt", ios_base::out);
        file1 << "x1 " << "x2 " << "x3 " << "v\n";
        short i1, i2;
        double t1, t2;
        for (size_t i = 0; i < _N; i++) {
            i1 = i / _n; i2 = i % _n;
            t1 = A + (i1 + 0.5) * h1;
            t2 = C + (i2 + 0.5) * h2;
            if (i % _n == 0 && i != 0) printf("\n");
            cout << abs(res[i]) * pow(10, 7) << " ";
            file1 << x1_screen(t1, t2, NULL, 0) << " " << x2_screen(t1, t2, NULL, 0) << " " << x3_screen(t1, t2, NULL, 0) << " " << abs(res[i]) << "\n";
        }
        file1.close();
    }


    MPI_Finalize();
}

//int main3g() {
//    return 0;
//    _N *= _n;
//    int _rank, _size;
//
//    MPI_Init(NULL, NULL);
//    MPI_Comm_size(MPI_COMM_WORLD, &_size);
//    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
//    MPI_Status status;
//    MPI_Request request;
//
//    /// подготовительные работы
//    int* count_one_rank = createv<int>(_size),* arr_step = createv<int>(_size);
//
//    //массив количества элементов на каждом ранге
//    for (size_t i = 0; i < _size; i++)
//        count_one_rank[i] = _N / _size;
//
//    short mod = _N % _size;
//    if (mod != 0)
//        for (short i = _size - 1; i >= 0; i--) {
//            count_one_rank[i]++;
//            mod--;
//            if (mod == 0) break;
//        }
//    //массив количества элементов на каждом ранге
//    //массив шагов для каждого ранга
//    for (size_t i = 0; i < _size; i++) {
//        arr_step[i] = 0;
//
//        for (size_t j = 0; j < i; j++)
//            arr_step[i] += count_one_rank[j];
//    }
//    //массив шагов для каждого ранга
//    /// подготовительные работы
//
//    double** a = createm<double>(count_one_rank[_rank], _N + 1);
//
//    double t1 = MPI_Wtime();
//    mg(a, 2, 3, arr_step[_rank]);
//    double t2 = MPI_Wtime() - t1;
//    printf("rank: %d  time fill matrix is: %f\n", _rank, t2);
//    fflush(stdout);
//    //if (_rank == 0) print(a);
//
//    t1 = MPI_Wtime();
//    gm(a, count_one_rank, arr_step[_rank], _rank, _size);
//    t2 = MPI_Wtime() - t1;
//    printf("rank: %d  time Gauss method is: %f\n", _rank, t2);
//    fflush(stdout);
//
//
//    double* part_res = createv<double>(count_one_rank[_rank]), * res = NULL;
//
//
//    for (size_t i = 0; i < count_one_rank[_rank]; i++)
//        part_res[i] = a[i][_N];
//
//
//    if (_rank == 0) {
//        del(a);
//        res = createv<double>(_N);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Gatherv(part_res, count_one_rank[_rank], MPI_DOUBLE, res, count_one_rank, arr_step, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//    if (_rank == 0)
//        for (size_t i = 0; i < _N; i++) {
//            printf("%f\n", res[i]);
//        }
//
//
//    MPI_Finalize();
//}
//
//int main3si() {
//    return 0;
//    //проверка правильности работы параллельного метода простой итерации
//    _N *= _n;
//    int _rank, _size;
//
//    MPI_Init(NULL, NULL);
//    MPI_Comm_size(MPI_COMM_WORLD, &_size);
//    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
//    MPI_Status status;
//    MPI_Request request;
//
//    /// подготовительные работы
//    int* count_one_rank = createv<int>(_size), * arr_step = createv<int>(_size);
//    //массив количества элементов на каждом ранге
//    for (size_t i = 0; i < _size; i++)
//        count_one_rank[i] = _N / _size;
//
//    short mod = _N % _size;
//    if (mod != 0)
//        for (short i = _size - 1; i >= 0; i--) {
//            count_one_rank[i]++;
//            mod--;
//            if (mod == 0) break;
//        }
//    //массив количества элементов на каждом ранге
//    //массив шагов для каждого ранга
//    for (size_t i = 0; i < _size; i++) {
//        arr_step[i] = 0;
//
//        for (size_t j = 0; j < i; j++)
//            arr_step[i] += count_one_rank[j];
//    }
//    //массив шагов для каждого ранга
//    /// подготовительные работы
//
//    double** a = createm<double>(count_one_rank[_rank], _N + 1), * init = createv<double>(_N);
//    double* part_res = createv<double>(count_one_rank[_rank]), * res = NULL;
//
//    init[0] = 0.6;
//    init[1] = 1.2;
//    init[2] = 1.2;
//    init[3] = 1.7;
//    init[4] = 1.2;
//    init[5] = 1.7;
//    init[6] = 1.7;
//    init[7] = 1.9;
//
//    double t1 = MPI_Wtime();
//    mg(a, 2, 3, arr_step[_rank]);
//    double t2 = MPI_Wtime() - t1;
//    //printf("rank: %d  time fill matrix is: %f\n", _rank, t2);
//    fflush(stdout);
//    //if (_rank == 0) print(a);
//
//    t1 = MPI_Wtime();
//    //gm(a, count_one_rank, arr_step[_rank], _rank, _size);
//    part_res = simple_iteration(a, init, 2, count_one_rank, arr_step, _rank, _size);
//    t2 = MPI_Wtime() - t1;
//    //printf("rank: %d  time Gauss method is: %f\n", _rank, t2);
//    fflush(stdout);
//
//
//    
//
//
//    /*for (size_t i = 0; i < count_one_rank[_rank]; i++)
//        part_res[i] = a[i][_N];*/
//
//
//    if (_rank == 0) {
//        del(a);
//        res = createv<double>(_N);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Gatherv(part_res, count_one_rank[_rank], MPI_DOUBLE, res, count_one_rank, arr_step, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//    if (_rank == 0)
//        for (size_t i = 0; i < _N; i++) {
//            printf("%f\n", res[i]);
//        }
//
//
//    MPI_Finalize();
//}



int /*main*/Проверка_базисной_функции () {
    return 0;

    double num_t = 1005.0;
    h1 = (B - A) / 4.0;
    h2 = (D - C) / 4.0;

    cout << h1 << " " << h2 << "\n";
    
    ofstream file1("v1.txt", ios_base::out);
    file1 << "x1 " << "x2 " << "x3 " << "v\n";

    ofstream file2("v2.txt", ios_base::out);
    file2 << "x1 " << "x2 " << "x3 " << "v\n";

    ofstream file3("v3.txt", ios_base::out);
    file3 << "x1 " << "x2 " << "x3 " << "v\n";
    for (double t1 = A; t1 < B; t1 += (B - A) / num_t) {
        for (double t2 = C; t2 < D; t2 += (D - C) / num_t) {
            double** res = createm<double>(3, 1);
            base_func(x1_screen, x2_screen, x3_screen, t1, t2, 1, 1, 0, 2, res);
            file1 << x1_screen(t1, t2, NULL, 0) << " " << x2_screen(t1, t2, NULL, 0) << " " << x3_screen(t1, t2, NULL, 0) << " " << res[0][0] << "\n";
            file2 << x1_screen(t1, t2, NULL, 0) << " " << x2_screen(t1, t2, NULL, 0) << " " << x3_screen(t1, t2, NULL, 0) << " " << res[1][0] << "\n";
            file3 << x1_screen(t1, t2, NULL, 0) << " " << x2_screen(t1, t2, NULL, 0) << " " << x3_screen(t1, t2, NULL, 0) << " " << res[2][0] << "\n";
            del(res);
        }
    }
    file1.close();
    file2.close();
    file3.close();
}

//int /*main*/Проверка_интеграла() {
//    return 0;
//    double** tensor = createm<double>(2, 2);
//    cout << "I2: " << I(100, x1_screen, x2_screen, x3_screen, tensor, -pi / 2.0, pi / 2.0, 0.0, 2.0 * pi) << "\n\n";
//
//    double** tensor = createm<double>(2, 2), ** tensor_reverse = createm<double>(2, 2);
//    cout << "INTEGRAL: " << S<complex<double>>(1, tensor, tensor_reverse, x1_screen, x2_screen, x3_screen, 1, 2, 1, 1, 1, 1) << "\n\n";
//    
//    cout << "F INTEGRAL: " << f<complex<double>>(1, tensor, x1_screen, x2_screen, x3_screen, 1, 1, 1) << "\n";
//}

