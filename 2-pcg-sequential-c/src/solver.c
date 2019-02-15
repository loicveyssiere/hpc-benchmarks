#include <stdio.h>
#include <stdlib.h>
#include "algebra.h"

int solver(sparse_coo* A, double* b, double* x, double* m, double tol, int max_iter, int size) {
    
    // INIT
    int count = 0;

    double rho = 0.0;
    double sigma = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    double norm2 = 0.0;

    double* r = (double*) malloc(size * sizeof(double));
    double* p = (double*) malloc(size * sizeof(double));
    double* q = (double*) malloc(size * sizeof(double));
    double* w = (double*) malloc(size * sizeof(double));

    // r = b - Ax
    matrix_vec_product(A, x, r, size);
    vec_sub(b, r, r, size);

    // w = M-1 * r
    vec_mul(m, r, w, size);

    // p = w
    vec_copy(w, p, size);

    // q = Ap
    matrix_vec_product(A, p, q, size);

    // rho = <r, w>
    rho = dot_product(r, w, size);

    // norm2 = <r, r>
    norm2 = dot_product(r, r, size);

    // sigma = <p, q>
    sigma = dot_product(p, q, size);

    printf("%d - Error value: %e\n", count, norm2);

    while (norm2 > tol && count < max_iter) {

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

        // w = M-1 * r
        vec_mul(m, r, w, size);

        // rho = <r, w>
        alpha = dot_product(r, w, size);

        // norm2 = <r, r>
        norm2 = dot_product(r, r, size);

        // beta = rho_new / rho_old
        beta = alpha / rho;
        rho = alpha;

        // p = w + beta * p
        vec_affine(w, p, beta, p, size);

        // q = Ap
        matrix_vec_product(A, p, q, size);

        // sigma = <p, q>
        sigma = dot_product(p, q, size);

        count++;

        printf("%d - Error value: %e\n", count, norm2);
    }
    
    return 0;
}