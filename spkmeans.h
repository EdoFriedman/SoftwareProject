#ifndef SPKMEANS_H
#define SPKMEANS_H

#include "matrix.h"

typedef struct {
    double* eigenvalues;
    Matrix eigenvectors;

} eigen_struct;

eigen_struct jacobi(Matrix A);
/**
 * Computes the Eigenvalues and Eigenvectors of a matrix A.
 */


double sqr_euclidean_dist(double* x1,  double* x2, size_t n_dim);
/**
 * Computes the squared euclidean distance of two n_dim dimensional
 * vectors x1, x2.
 */

Matrix get_w_adj_matrix(double** datapoints, size_t n_datapoints, size_t n_dim);
/**
 * Computes the Weighted Adjacency Matrix.
 */

Matrix get_ddg_matrix(Matrix w_adj_matrix);
/**
 * Computes the Diagonal Degree Matrix.
 */

Matrix get_laplacian_matrix(double **datapoints, size_t n_datapoints, size_t n_dim);
/**
 * Computes the Graph Laplacian Matrix.
 */

#endif
