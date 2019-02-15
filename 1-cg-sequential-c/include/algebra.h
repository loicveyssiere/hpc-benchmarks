#ifndef _ALGEBRA_H
#define _ALGEBRA_H

struct sparse_coo {
    int* index_1;
    int* index_2;
    double* values;
    int length;
};

typedef struct sparse_coo sparse_coo;

void init_sparse_coo(sparse_coo* matrix, int length);

/* =============================================================================
LINEAR ALGEBRA
============================================================================= */

inline double dot_product(double* vec1, double* vec2, int size);

inline void matrix_vec_product(sparse_coo* matrix, double* vec, double* result, int size);

inline void vec_copy(double* from, double* to, int size);

inline void vec_set(double* vec, double value, int start, int end);

inline void vec_sub(double* vec1, double* vec2, double* result, int size);

inline void vec_add(double* vec1, double* vec2, double* result, int size);

inline void vec_mul(double* vec1, double* vec2, double* result, int size);

inline void scalar_mul(double* vec, double scalar, double* result, int size);

inline void vec_affine(double* vec1, double* vec2, double scalar, double* result, int size);

#endif