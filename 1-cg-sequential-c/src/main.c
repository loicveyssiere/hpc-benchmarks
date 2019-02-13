#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solver.h"
#include "algebra.h"


int main(int argc, char** argv) {

    printf("Welcome to Conjugate Gradient\n");

    int size = 10000;
    int length = size + 2*(size-1) + 2*(size-3);
    double tol = 10E-10;
    int max_iter = 3000;

    double* b = (double*) malloc(size * sizeof(double));
    double* x = (double*) malloc(size * sizeof(double));

    sparse_coo matrix;
    printf("TEST1\n");
    init_sparse_coo(&matrix, length);
    printf("TEST2\n");

    for(int i = 0; i < size; i++) {
        b[i] = 1.0;
        x[i] = 0.0;
    }

    int k = 0;
    for (int i = 0; i < size; i++) {
        matrix.index_1[k] = i;
        matrix.index_2[k] = i;
        matrix.values[k] = 4.0;
        k++;
    }

    for (int i = 0; i < size-1; i++) {
        matrix.index_1[k] = i;
        matrix.index_2[k] = i+1;
        matrix.values[k] = 1.0;
        k++;
        matrix.index_1[k] = i+1;
        matrix.index_2[k] = i;
        matrix.values[k] = 1.0;
        k++;
    }

    for (int i = 0; i < size-3; i++) {
        matrix.index_1[k] = i+2;
        matrix.index_2[k] = i+1;
        matrix.values[k] = 1.0;
        k++;
        matrix.index_1[k] = i+1;
        matrix.index_2[k] = i+2;
        matrix.values[k] = 1.0;
        k++;
    }

    printf("%d : %d\n", k, length);

    solver(&matrix, b, x, tol, max_iter, size);
    printf("End of Conjugate Gradient");

    return 0;
}