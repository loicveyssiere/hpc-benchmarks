/* =============================================================================
LINEAR ALGEBRA
============================================================================= */

#include <stdio.h>
#include <stdlib.h>
#ifdef OMP
    #include <omp.h>
#endif

#include "algebra.h"

void init_sparse_coo(sparse_coo* matrix, int length) {
    matrix->length = length;
    matrix->index_1 = (int*) malloc(length * sizeof(int));
    matrix->index_2 = (int*) malloc(length * sizeof(int));
    matrix->values = (double*) malloc(length * sizeof(double));
}

inline double dot_product(double* vec1, double* vec2, int size) {
    int i;
    double result = 0.0;
    for (i = 0; i < size; i++) {
        result += vec1[i] * vec2[i];
    }
    return result;
}

inline void matrix_vec_product(sparse_coo* matrix, double* vec, double* result, int size) {
    vec_set(result, 0.0, 0, size);
    int i;
    for(i = 0; i < matrix->length; i++) {
        result[matrix->index_1[i]] += matrix->values[i] * vec[matrix->index_2[i]];
    }
}

inline void vec_copy(double* from, double* to, int size) {
    for(int i = 0; i < size; i++) {
        to[i] = from[i];
    }
}

inline void vec_set(double* vec, double value, int start, int end) {
    for(int i = start; i < end; i++) {
        vec[i] = value;
    }
}

inline void vec_sub(double* vec1, double* vec2, double* result, int size) {
    for (int i = 0; i < size; i++) {
        result[i] = vec1[i] - vec2[i];
    }
}

inline void vec_add(double* vec1, double* vec2, double* result, int size) {
    for (int i = 0; i < size; i++) {
        result[i] = vec1[i] + vec2[i];
    }
}

inline void vec_mul(double* vec1, double* vec2, double* result, int size) {
    for (int i = 0; i < size; i++) {
        result[i] = vec1[i] * vec2[i];
    }
}

inline void scalar_mul(double* vec, double scalar, double* result, int size) {
    for (int i = 0; i < size; i++) {
        result[i] = scalar * vec[i];
    }
}

inline void vec_affine(double* vec1, double* vec2, double scalar, double* result, int size) {
    int i;
    //#pragma omp parallel for shared(vec1, vec2, scalar, result, size) private(i)
    for (i = 0; i < size; i++) {
        result[i] = vec1[i] + scalar * vec2[i];
    }
}