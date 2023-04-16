#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "matrix.h"
#include "spkmeans.h"

#define MAX_LINE_LENGTH 4096

double euclideanDistance(double *point1, double *point2, size_t dim) {
    double sum;
    size_t i;
    sum = 0;
    for (i = 0; i < dim; i++) {
        sum += (point1[i] - point2[i]) * (point1[i] - point2[i]);
    }
    return sqrt(sum);
}

int distanceArgmin(double **centroids, double *point, size_t centroid_count, size_t dim) {
    int min_index;
    double min_distance;
    size_t i;
    double distance;
    min_index = 0;
    min_distance = euclideanDistance(centroids[0], point, dim);
    for (i = 1; i < centroid_count; i++) {
        if (min_distance > (distance = euclideanDistance(centroids[i], point, dim))) {
            min_index = i;
            min_distance = distance;
        }
    }
    return min_index;
}

double **kmeans(kmeans_input input, double **initial_centroids) {
    double **centroids;
    double **new_centroids;
    double **cluster_point_sums;
    size_t *cluster_point_counts;
    int closest_centroid;
    size_t i;
    size_t j;
    size_t k;
    int converged;
    centroids = malloc(input.cluster_count * sizeof(double *));
    for (i = 0; i < input.cluster_count; i++) {
        centroids[i] = malloc(input.dim * sizeof(double));
        for (j = 0; j < input.dim; j++) {
            centroids[i][j] = initial_centroids[i][j];
        }
    }
    cluster_point_sums = malloc(input.cluster_count * sizeof(double *));
    cluster_point_counts = malloc(input.cluster_count * sizeof(size_t));
    for (i = 0; i < input.cluster_count; i++) {
        cluster_point_sums[i] = malloc(input.dim * sizeof(double));
    }
    for (i = 0; i < input.iter_count; i++) {
        memset(cluster_point_counts, 0, input.cluster_count * sizeof(size_t));
        for (j = 0; j < input.cluster_count; j++) {
            memset(cluster_point_sums[j], 0, input.dim * sizeof(double));
        }
        for (j = 0; j < input.datapoint_count; j++) {
            closest_centroid = distanceArgmin(centroids, input.datapoints[j], input.cluster_count,
                                              input.dim);
            for (k = 0; k < input.dim; k++) {
                cluster_point_sums[closest_centroid][k] += input.datapoints[j][k];
            }
            cluster_point_counts[closest_centroid]++;
        }
        converged = !0;
        new_centroids = malloc(input.cluster_count * sizeof(double *));
        for (j = 0; j < input.cluster_count; j++) {
            new_centroids[j] = malloc(input.dim * sizeof(double));
            for (k = 0; k < input.dim; k++) {
                new_centroids[j][k] = cluster_point_sums[j][k] / cluster_point_counts[j];
            }
            converged = converged && euclideanDistance(centroids[j], new_centroids[j], input.dim) < input.epsilon;
            free(centroids[j]);
        }
        free(centroids);
        centroids = new_centroids;
        if (converged) break;
    }
    for (i = 0; i < input.cluster_count; i++) {
        free(cluster_point_sums[i]);
    }
    free(cluster_point_counts);
    free(cluster_point_sums);
    return centroids;
}


void print_matrix(Matrix A) {
    size_t i, j;
    char *delim;
    for (i = 0; i < A.shape[0]; i++) {
        delim = "";
        for (j = 0; j < A.shape[1]; j++) {
            printf("%s%.4f", delim, A.data[i][j]);
            delim = ",";
        }
        printf("\n");
    }
}

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

int *get_pivot(Matrix A) {
    double max_val = 0;
    int *idx = (int *) calloc(2, sizeof(int));
    size_t i, j;

    for (i = 0; i < A.shape[0]; i++) {
        for (j = 0; j < A.shape[1]; j++) {
            if (i != j) {
                if (max_val < fabs(A.data[i][j])) {
                    idx[0] = i;
                    idx[1] = j;
                    max_val = fabs(A.data[i][j]);
                }
            }
        }

    }

    return idx;

}

double off_diagonal_sum(Matrix A) {
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

eigen_struct jacobi(Matrix A) {

    int *pivot_idx;
    size_t i, j;
    size_t ii;
    double c, t, s, theta;
    Matrix A_tag, P, V_tag, temp1, temp2;
    Matrix V = mat_id(A.shape[0]);
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

        P = get_rotation_matrix(c, s, i, j, A.shape[0]);
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
        } else {
            mat_free(A);
            A = A_tag;
            mat_free(V);
            V = V_tag;
        }
    }

    res.eigenvalues = malloc(A.shape[0] * sizeof(double));
    for (ii = 0; ii < A.shape[0]; ii++) {
        res.eigenvalues[ii] = A.data[ii][ii];
    }
    res.eigenvectors = V;
    mat_free(A);


    return res;
}

double **parse_csv(char *filename, size_t *n_datapoints, size_t *n_dim) {
    char *token;
    FILE *file;
    char line[MAX_LINE_LENGTH];
    size_t line_cap = 8;
    size_t dim;
    double **res;
    size_t i, j;

    res = malloc(line_cap * sizeof(double *));


    file = fopen(filename, "r");

    fgets(line, MAX_LINE_LENGTH, file);
    token = strtok(line, ",");

    dim = 0;
    while (token != NULL) {
        dim++;
        token = strtok(NULL, ",");
    }

    rewind(file);

    i = 0;
    while (fgets(line, MAX_LINE_LENGTH, file)) {
        if (i >= line_cap) {
            line_cap *= 2;
            res = realloc(res, line_cap * sizeof(double *));
        }

        res[i] = malloc(dim * sizeof(double));
        j = 0;

        token = strtok(line, ",");

        while (token != NULL) {
            res[i][j] = strtod(token, NULL);

            j++;
            token = strtok(NULL, ",");
        }

        i++;

    }

    fclose(file);

    *n_dim = dim;
    *n_datapoints = i;

    return res;

}

int main(int argc, char **argv) {
    char *goal;
    char *filename;
    double **input;
    size_t n_datapoints, n_dims;
    Matrix wam_output, ddg_output, gl_output, jacobi_input;
    eigen_struct jacobi_output;
    size_t i;
    char *delim;

    (void) argc;

    goal = argv[1];
    filename = argv[2];

    input = parse_csv(filename, &n_datapoints, &n_dims);

    if (strcmp(goal, "wam") == 0) {
        wam_output = get_w_adj_matrix(input, n_datapoints, n_dims);
        print_matrix(wam_output);
        mat_free(wam_output);
    }
    if (strcmp(goal, "ddg") == 0) {
        wam_output = get_w_adj_matrix(input, n_datapoints, n_dims);
        ddg_output = get_ddg_matrix(wam_output);
        print_matrix(ddg_output);
        mat_free(wam_output);
        mat_free(ddg_output);
    }
    if (strcmp(goal, "gl") == 0) {
        gl_output = get_laplacian_matrix(input, n_datapoints, n_dims);
        print_matrix(gl_output);
        mat_free(gl_output);
    }
    if (strcmp(goal, "jacobi") == 0) {
        jacobi_input.data = input;
        jacobi_input.shape[0] = n_datapoints;
        jacobi_input.shape[1] = n_dims;

        jacobi_output = jacobi(jacobi_input);

        delim = "";
        for (i = 0; i < n_datapoints; i++) {
            printf("%s%.4f", delim, jacobi_output.eigenvalues[i]);
            delim = ",";
        }
        printf("\n");
        print_matrix(jacobi_output.eigenvectors);
        mat_free(jacobi_output.eigenvectors);
        free(jacobi_output.eigenvalues);
    }

    for(i = 0; i < n_datapoints; i++) {
        free(input[i]);
    }
    free(input);

    return 0;
}