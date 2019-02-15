#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef TIMER
    #include <time.h>
#endif

#ifdef OMP
    #include <omp.h>
#endif

#include "solver.h"
#include "algebra.h"


int main(int argc, char** argv) {

    printf("Welcome to Conjugate Gradient================================= \n");

    #ifdef OMP
    printf("OMP Activated\n");
    #pragma omp parallel
    {
        int numberOfThreads = omp_get_num_threads();
        int ID = omp_get_thread_num();
        printf("%d %d with max= \n",ID,numberOfThreads);
    }
    #endif
    
    int size = 1000000;
    int length = size + 2*(size-1) + 2*(size-3);
    double tol = 10E-10;
    int max_iter = 3000;

    double* b = (double*) malloc(size * sizeof(double));
    double* x = (double*) malloc(size * sizeof(double));
    double* m = (double*) malloc(size * sizeof(double));

    sparse_coo matrix;
    init_sparse_coo(&matrix, length);

    for(int i = 0; i < size; i++) {
        b[i] = 1.0;
        x[i] = 0.0;
        #ifdef SIMPLE_CG
            m[i] = 1.0;
        #else
            m[i] = 1.0 / ((double) log(i+1) + 1);
        #endif
    }

    int k = 0;
    for (int i = 0; i < size; i++) {
        matrix.index_1[k] = i;
        matrix.index_2[k] = i;
        matrix.values[k] = ((double) log(i+1) + 1);
        k++;
    }

    for (int i = 0; i < size-1; i++) {
        matrix.index_1[k] = i;
        matrix.index_2[k] = i+1;
        matrix.values[k] = -1.0;
        k++;
        matrix.index_1[k] = i+1;
        matrix.index_2[k] = i;
        matrix.values[k] = -1.0;
        k++;
    }

    for (int i = 0; i < size-3; i++) {
        matrix.index_1[k] = i+2;
        matrix.index_2[k] = i+1;
        matrix.values[k] = -1.0;
        k++;
        matrix.index_1[k] = i+1;
        matrix.index_2[k] = i+2;
        matrix.values[k] = -1.0;
        k++;
    }

    #ifdef TIMER
        clock_t start, end;
        double time_difference;
        start = clock();
    #endif
    solver(&matrix, b, x, m, tol, max_iter, size);
    #ifdef TIMER
        end = clock();
        time_difference = (double) (end - start) / CLOCKS_PER_SEC;
        printf("TIMER: %f seconds\n", time_difference);
    #endif
    printf("End of Conjugate Gradient");

    return 0;
}