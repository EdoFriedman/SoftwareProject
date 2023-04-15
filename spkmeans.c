#include "spkmeans.h"
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

Matrix wam(double **datapoints, size_t n_datapoints, size_t n_dim) {
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

Matrix ddg(Matrix wam) {
    double w_sum;
    size_t i, j;
    Matrix res = mat_zeroes(wam.shape[0], wam.shape[1]);

    for (i = 0; i < wam.shape[0]; i++) {
        w_sum = 0;
        for (j = 0; j < wam.shape[1]; j++) {
            w_sum += wam.data[i][j];
        }

        res.data[i][i] = w_sum;
    }

    return res;
}

Matrix gl(double **datapoints, size_t n_datapoints, size_t n_dim) {
    Matrix w_adj_matrix = wam(datapoints, n_datapoints, n_dim);
    Matrix ddg_matrix = ddg(w_adj_matrix);
    Matrix laplacian_matrix = mat_add(ddg_matrix, mat_mul_s(w_adj_matrix, -1));

    mat_free(w_adj_matrix);
    mat_free(ddg_matrix);

    return laplacian_matrix;
}

// return a pair of indexes of the pivot according to section #3 in the jacobi algorithm
int* get_pivot(Matrix A) {
    double max_val = 0;
    int* idx = (int*)calloc(2, sizeof(int));
    size_t i, j;

    for (i = 0; i < A.shape[0]; i++) {
        for (j = 0; j < A.shape[1]; j++) {
            if (i != j) {
                if (max_val < A.data[i][j]) {
                    idx[0] = i;
                    idx[1] = j;
                    max_val = A.data[i][j];
                }
            }
        }

    }

    return idx;

}

double off_diagonal_sum(Matrix A){
    double sum = 0;
    size_t i, j;

    for (i = 0; i < A.shape[0]; i++) {
        for (j = 0; j < A.shape[1]; j++) {
            if (i != j) {
                sum += A.data[i][j];
            }
        }

    }

    return sum;

}

Matrix get_rotation_matrix(double c, double s, int i, int j, size_t n) {
    Matrix P = mat_id(n);
    P.data[i][i] = c;
    P.data[j][j] = c;
    P.data[i][j] = s;
    P.data[j][i] = -s;

    return P;
}

eigen_struct jacobi(Matrix A){
    // find pivot DONE
    // get t, c, s, theta DONE
    // get P DONE
    // V_tag = P * V DONE
    // A_tag = P.T * A * P DONE
    // if not converged: A = A_tag, V = V_tag. DONE
    // else return.

    int* pivot_idx;
    size_t i, j;
    size_t ii;
    double c, t, s, theta;
    Matrix A_tag, P, V_tag, temp1, temp2;
    Matrix V = mat_id(A.shape[0]); // A.shape is n, n;
    eigen_struct res;

    A = mat_clone(A);
    for (ii = 0; ii < 100; ii++) {
        pivot_idx = get_pivot(A);
        i = pivot_idx[0];
        j = pivot_idx[1];
        free(pivot_idx);

        theta = (A.data[j][j] - A.data[i][i]) / (2 * A.data[i][j]);
        t = (theta >= 0 ? 1 : -1) / (fabs(theta) + sqrt(theta * theta + 1));
        c = 1 / (sqrt(t * t + 1));
        s = t * c;

        P = get_rotation_matrix(c, s, i, j, A.shape[0]); // A.shape is n, n;
        V_tag = mat_mul_m(V, P);

        temp1 = mat_transpose(P);
        temp2 = mat_mul_m(temp1, A);
        mat_free(temp1);
        A_tag = mat_mul_m(temp2, P);
        mat_free(temp2);


        if (fabs(off_diagonal_sum(A) - off_diagonal_sum(A_tag)) < 0.00001) {
            mat_free(A);
            A = A_tag;
            mat_free(V);
            V = V_tag;
            break;
        }
        else {
            mat_free(A);
            A = A_tag;
            mat_free(V);
            V = V_tag;
        }
    }

    // return eigenvectors and eigenvalues.
    res.eigenvalues = malloc(A.shape[0] * sizeof(double));
    for (ii = 0; i < A.shape[0]; i++) {
        res.eigenvalues[ii] = A.data[ii][ii];
    }
    res.eigenvectors = V;
    mat_free(A);


    return res;
}

int main(){
    return 0;
};