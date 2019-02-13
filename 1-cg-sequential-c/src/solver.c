#include <stdio.h>
#include <stdlib.h>
#include "algebra.h"

int solver(sparse_coo* A, double* b, double* x, double tol, int max_iter, int size) {
    
    // INIT
    int count = 0;

    double rho = 0.0;
    double sigma = 0.0;
    double alpha = 0.0;
    double beta = 0.0;

    double* r = (double*) malloc(size * sizeof(double));
    double* p = (double*) malloc(size * sizeof(double));
    double* q = (double*) malloc(size * sizeof(double));

    // r = b - Ax
    matrix_vec_product(A, x, r, size);
    vec_sub(b, r, r, size);

    // p = r
    vec_copy(r, p, size);

    // q = Ap
    matrix_vec_product(A, p, q, size);

    // rho = <r, r>
    rho = dot_product(r, r, size);

    // sigma = <p, q>
    sigma = dot_product(p, q, size);

    printf("%d - Error value: %e\n", count, rho);

    while (rho > tol && count < max_iter) {

        // alpha = rho / sigma
        //printf("%d - Sigma value: %f\n", count, sigma);
        if (sigma == 0.0)
            alpha = 0.0;
        else
            alpha = rho / sigma;

        // x = x + alpha * p
        vec_affine(x, p, alpha, x, size);

        // r = r - alpha * q
        vec_affine(r, q, -1 * alpha, r, size);

        // rho = <r, r>
        alpha = dot_product(r, r, size);

        // beta = rho_new / rho_old
        beta = alpha / rho;
        rho = alpha;

        // p = r + beta * p
        vec_affine(r, p, beta, p, size);

        // q = Ap
        matrix_vec_product(A, p, q, size);

        // sigma = <p, q>
        sigma = dot_product(p, q, size);

        count++;

        printf("%d - Error value: %e\n", count, rho);
    }
    
    return 0;
}