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
            double ksi = A_obj + (i + 0.5) * h1_obj;

            for (size_t j = 0; j < M; j++){
                double a = A_obj + j * h1_obj, b = A_obj + (j + 1.0) * h1_obj;
                
                var[i][j] = base_func(i, j) * (type - 1.0) - lambda * I_k(100, a, b, ksi);
            }
            var[i][M] = func(ksi);
        }
    else
        for (size_t i = 0; i < M; i++){
            short i1 = i / m, i2 = i % m;
            double ksi1 = A_obj + (i1 + 0.5) * h1_obj, ksi2 = A_obj + (i2 + 0.5) * h2_obj;

            for (size_t j = 0; j < M; j++) {
                short j1 = j / m, j2 = j % m;
                double a = A_obj + j1 * h1_obj, b = A_obj + (j1 + 1.0) * h1_obj, c = C_obj + j2 * h2_obj, d = C_obj + (j2 + 1.0) * h2_obj;

                var[i][j] = base_func(i, j) * (type - 1.0) - lambda * I_k(100, a, b, c, d, ksi1, ksi2);     
            }
            var[i][M] = func(ksi1, ksi2);
        }
}


void mg(complex<double>** var, size_t type, size_t dim_s = 1, int _step = 0, int _size = 0, int* count_one_rank = NULL) {
    //метод Галёркина
    //var - матрица, в которую будут записываться данные 
    //type - тип уравнения Фредгольма, первого или второго рода
    //dim_s - размер измерения в котором считается метод Галёркина
    MPI_Request request;
    MPI_Status status;
    size_t M = _msize(var) / sizeof(var[0]);
    size_t N = _msize(var[0]) / sizeof(var[0][0]);
    double** tensor, ** tensor_reverse;
    short i1, i2, j1, j2, _k, _l;
    if (dim_s == 1) {
        for (size_t i = 0; i < M; i++) {
            double a = A_obj + (i + _step) * h1_obj, b = A_obj + ((i + _step) + 1.0) * h1_obj;

            for (size_t j = 0; j < N; j++) {
                double c = A_obj + j * h1_obj, d = A_obj + (j + 1.0) * h1_obj;

                var[i][j] = h1_obj * base_func((i + _step), j) * (type - 1.0) - lambda * In<double>(3, k, a, b, c, d);
            }
            var[i][N - 1] = In<double>(3, a, b);
        }
    }
    else if (dim_s == 2) {
        int m = (int)(1 + sqrt(N)) / 2;
        //cout << "m: " << m << "\n";
        for (size_t i = 0; i < M; i++) {
            _k = ((i + _step) / ((N - 1) / 4)) % 2 + 1;

            if (_k == 1) { i1 = ((i + _step) / m) % (m - 1); i2 = (i + _step) % m; }
            else { i1 = ((i + _step) / (m - 1)) % m; i2 = (i + _step) % (m - 1); }

            if (_step == 0) {
                printf("%d\n", i);
                fflush(stdout);
            }

            tensor = createm<double>(2, 2);
            tensor_reverse = createm<double>(2, 2);
            for (size_t j = 0; j < (N - 1) / 2; j++) {
                _l = (j / ((N - 1) / 4)) % 2 + 1;

                if (_l == 1) { j1 = (j / m) % (m - 1); j2 = j % m; }
                else { j1 = (j / (m - 1)) % m; j2 = j % (m - 1); }



                //complex<double> i_ksi = f_vec<complex<double>>(2, permittivity_obj, tensor, x1_obj, x2_obj, x3_obj, _k, i1, i2, 1);

                if ((i + _step) < (N - 1) / 2) {
                    var[i][j + (N - 1) / 2] = 1;//-lambda * S_obj_screen<complex<double>>(4, tensor, tensor_reverse, x1_obj, x2_obj, x3_obj, _k, _l, i1, i2, j1, j2);
                    int rank = 0, step = 0;
                    for (rank; rank < _size; rank++) {
                        step += count_one_rank[rank];
                        if (step > (i + _step + (N - 1) / 2)) break;
                    }
                    
                    fflush(stdout);
                    MPI_Isend(&var[i][j + (N - 1) / 2], 1, MPI_DOUBLE_COMPLEX, rank, 3, MPI_COMM_WORLD, &request);
                    
                    var[i][j] = lambda * S_obj<complex<double>>(x1_obj, x2_obj, x3_obj, 4, _k, _l, i1, i2, j1, j2, tensor, tensor_reverse);
                
                }
                else {
                    var[i][j + (N - 1) / 2] = 1; //-lambda * S_screen<complex<double>>(4, tensor, tensor_reverse, x1_screen, x2_screen, x3_screen, _k, _l, i1, i2, j1, j2);
                    

                    MPI_Irecv(&var[i][j], 1, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &request);
                    MPI_Wait(&request, &status);
                }
            }
            if ((i + _step) < (N - 1) / 2) var[i][N - 1] = 1; //-f_vec<complex<double>>(2, func_cv, tensor, x1_obj, x2_obj, x3_obj, _k, i1, i2, 1);
            else var[i][N - 1] = 1;//-f_vec<complex<double>>(2, func_cv_tang, tensor, x1_screen, x2_screen, x3_screen, _k, i1, i2, 0);

            del(tensor); del(tensor_reverse);
        }
    }
    else if (dim_s == 3) {
        int m = (int)pow(N - 1, 1.0 / 3.0);

        for (size_t i = 0; i < M; i++) {
            short i1 = (i + _step) % m, i2 = (i + _step) / m - m * ((i + _step) / m / m), i3 = (i + _step) / m / m;
            double a = A_obj + i1 * h1_obj, b = a + h1_obj;
            double c = C_obj + i2 * h2_obj, d = c + h2_obj;
            double e = E_obj + i3 * h3_obj, f = e + h3_obj;

            for (size_t j = 0; j < N; j++) {
                short j1 = j % m, j2 = j / m - m * (j / m / m), j3 = j / m / m;

                double g = A_obj + j1 * h1_obj, l = g + h1_obj;
                double o = C_obj + j2 * h2_obj, p = o + h2_obj;
                double q = E_obj + j3 * h3_obj, s = q + h3_obj;

                var[i][j] = h1_obj * h2_obj * h3_obj * base_func(i + _step, j) * (type - 1.0) - lambda * In<double>(5, k, a, b, c, d, e, f, g, l, o, p, q, s);
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

    complex<double>** a = createm<complex<double>>(count_one_rank[_rank], _N + 1);

    double t1 = MPI_Wtime();
    mg(a, 2, 2, _step, _size, count_one_rank);
    double t2 = MPI_Wtime() - t1;
    printf("rank: %d  time fill matrix is: %f\n", _rank, t2);
    fflush(stdout);
    //if (_rank == 0) print(absm(a), "g");
    if (_rank == 0) print(absm(a), "g", 8, 8);
    //print(col(a, _N));

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
        ofstream file1("res1.txt", ios_base::out);
        file1 << "x1 x2 x3 v_real v_imag v_abs\n";
        ofstream file2("res2.txt", ios_base::out);
        file2 << "x1 x2 x3 v_real v_imag v_abs\n";
        ofstream file3("res3.txt", ios_base::out);
        file3 << "x1 x2 x3 v_real v_imag v_abs\n";
        ofstream file4("res4.txt", ios_base::out);
        file4 << "x1 x2 x3 v_real v_imag v_abs\n";

        short i1, i2;
        double x1, x2, x3;
        double* v;
        complex<double>* sum; int count1 = 0, count2 = 0;;
        for (double t1 = A_screen; t1 < B_screen; t1 += (B_screen - A_screen) / 100.0) {
            for (double t2 = C_screen; t2 < D_screen; t2 += (D_screen - C_screen) / 100.0) {
                x1 = x1_screen(0, t1, t2, 0, NULL); x2 = x2_screen(0, t1, t2, 0, NULL, 0); x3 = x3_screen(0, t1, t2, 0, NULL);
                sum = createv<complex<double>>(3, true);
                for (size_t i = _N / 2; i < _N; i++) {
                    v = createv<double>(3);
                    if (i < 3 * _N / 4) {
                        i1 = (i / _n) % (_n - 1); i2 = (i) % _n; 
                        base_func_rft(x1_screen, x2_screen, x3_screen, t1, t2, 0, i1, i2, 0, 1, 0, v);
                    }
                    else {
                        i1 = ((i) / (_n - 1)) % _n; i2 = (i) % (_n - 1);
                        base_func_rft(x1_screen, x2_screen, x3_screen, t1, t2, 0, i1, i2, 0, 2, 0, v);
                    }      
                    sum[0] += res[i] * v[0];
                    sum[1] += res[i] * v[1];
                    sum[2] += res[i] * v[2];
                    del(v);
                }

                file1 << x1 << " " << x2 << " " << x3 << " " << sum[0].real() << " " << sum[0].imag() << " " << abs(sum[0]) << "\n";
                file2 << x1 << " " << x2 << " " << x3 << " " << sum[1].real() << " " << sum[1].imag() << " " << abs(sum[1]) << "\n";
                file3 << x1 << " " << x2 << " " << x3 << " " << sum[2].real() << " " << sum[2].imag() << " " << abs(sum[2]) << "\n";
                file4 << x1 << " " << x2 << " " << x3 << " " << sqrt(pow(sum[0].real(), 2) + pow(sum[1].real(), 2) + pow(sum[2].real(), 2)) << " " << sqrt(pow(sum[0].imag(), 2) + pow(sum[1].imag(), 2) + pow(sum[2].imag(), 2)) << " " << sqrt(pow(abs(sum[0]), 2) + pow(abs(sum[1]), 2) + pow(abs(sum[2]), 2)) << "\n";
                del(sum);
            }
        }
        file1.close();
        file2.close();
        file3.close();
        file4.close();
    }


    MPI_Finalize();
}

