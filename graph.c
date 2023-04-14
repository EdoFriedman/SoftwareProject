#include "graph.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "matrix.h"

double sqr_euclidean_dist(double *x1, double *x2, size_t n_dim) {
    size_t i;
    double res = 0;
    for (i = 0; i < n_dim; i++) {
        res += pow(x1[i] - x2[i], 2);
    }

    return res;
}

Matrix get_w_adj_matrix(double **datapoints, size_t n_datapoints, size_t n_dim) {
    Matrix w_adj_mat = mat_zeroes(n_datapoints, n_datapoints);

    size_t i, j;
    for (i = 0; i < n_datapoints; i++) {
        for (j = 0; j < n_datapoints; j++) {
            if (i != j) {
                w_adj_mat.data[i][j] = exp(-sqr_euclidean_dist(datapoints[i],
                                                               datapoints[j], n_dim) / 2);
            }
        }
    }

    return w_adj_mat;
}

Matrix get_ddg_matrix(Matrix w_adj_matrix) {
    double w_sum;
    size_t i, j;
    Matrix res = mat_zeroes(w_adj_matrix.shape[0], w_adj_matrix.shape[0]);

    for (i = 0; i < w_adj_matrix.shape[0]; i++) {
        w_sum = 0;
        for (j = 0; j < w_adj_matrix.shape[1]; j++) {
            w_sum += w_adj_matrix.data[i][j];
        }

        res.data[i][i] = w_sum;
    }

    return res;
}

Matrix get_laplacian_matrix(double **datapoints, size_t n_datapoints, size_t n_dim) {
    Matrix w_adj_matrix = get_w_adj_matrix(datapoints, n_datapoints, n_dim);
    Matrix ddg_matrix = get_ddg_matrix(w_adj_matrix);
    Matrix laplacian_matrix = mat_add(ddg_matrix, mat_mul_s(w_adj_matrix, -1));

    return laplacian_matrix;
}