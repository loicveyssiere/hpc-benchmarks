
#ifndef _SOLVER_H
#define _SOLVER_H

#include "algebra.h"

int solver(sparse_coo* A, double* b, double* x, double tol, int max_iter, int size);

#endif
