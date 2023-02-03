#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "matrix.h"

#define _MAT_SET_VALUE(mat, value) {\
    size_t i,j;\
    for(i = 0; i < mat.shape[0]; i++) {\
        for(j = 0; j < mat.shape[1]; j++) {\
            mat.data[i][j] = value;\
        }\
    }\
}

Matrix mat_zeroes(size_t rows, size_t columns) {
    Matrix mat;
    size_t i;

    mat.shape[0] = rows;
    mat.shape[1] = columns;

    mat.data = calloc(rows, sizeof(double*));
    for(i = 0; i < rows; i++) {
        mat.data[i] = calloc(columns, sizeof(double));
    }
    return mat;
}

Matrix mat_id(size_t size) {
    Matrix mat;
    size_t i;

    mat = mat_zeroes(size, size);

    for(i = 0; i < size; i++) {
        mat.data[i][i] = 1;
    }
    return mat;
}

Matrix mat_from_func(size_t rows, size_t columns, double (*func)(double, double)) {
    Matrix mat;
    
    mat = mat_zeroes(rows, columns);

    _MAT_SET_VALUE(mat, func(i, j));

    return mat;
}

Matrix mat_clone(Matrix mat) {
    Matrix new_mat;
    size_t i;

    new_mat = mat_zeroes(mat.shape[0], mat.shape[1]);

    for(i = 0; i < mat.shape[0]; i++) {
        memcpy(new_mat.data[i], mat.data[i], mat.shape[1]);
    }
    
    return new_mat;
}

void mat_free(Matrix mat) {
    size_t i;

    for(i = 0; i < mat.shape[0]; i++) {
        free(mat.data[i]);
    }
    free(mat.data);
}

Matrix mat_mul_s(Matrix mat, double scalar) {
     Matrix res;
     res = mat_clone(mat);
    _MAT_SET_VALUE(res, mat.data[i][j] * scalar);
    return res;
}

Matrix mat_mul_m(Matrix mat1, Matrix mat2) {
    Matrix res;
    size_t i, j, k;
    double sum;

    assert(mat1.shape[1] == mat2.shape[0]);

    res = mat_zeroes(mat1.shape[0], mat2.shape[1]);
    for(i = 0; i < res.shape[0]; i++) {
        for(j = 0; j < res.shape[1]; j++) {
            sum = 0;
            for(k = 0; k < mat1.shape[1]; k++) {
                sum += mat1.data[i][k] * mat2.data[k][j];
            }
            res.data[i][j] = sum;
        }
    }
    return res;
}

Matrix mat_add(Matrix mat1, Matrix mat2) {
    Matrix res;

    assert(mat1.shape[0] == mat2.shape[0] && mat1.shape[1] == mat2.shape[1]);

    res = mat_zeroes(mat1.shape[0], mat2.shape[0]);

    _MAT_SET_VALUE(res, mat1.data[i][j] + mat2.data[i][j]);

    return res;
}

Matrix mat_elem_pow(Matrix mat, double power) {
    Matrix res;

    res = mat_clone(mat);

    _MAT_SET_VALUE(res, pow(res.data[i][j], power));

    return res;
}


Matrix mat_sum_cols(Matrix mat) {
    Matrix res;
    size_t i, j;

    res = mat_zeroes(1, mat.shape[1]);

    for(i = 0; i < mat.shape[0]; i++) {
        for(j = 0; j < mat.shape[1]; j++) {
            res.data[0][j] += mat.data[i][j];
        }
    }
    return res;
}

Matrix mat_sum_rows(Matrix mat) {
    Matrix res;
    size_t i, j;

    res = mat_zeroes(mat.shape[0], 1);

    for(i = 0; i < mat.shape[0]; i++) {
        for(j = 0; j < mat.shape[1]; j++) {
            res.data[i][0] += mat.data[i][j];
        }
    }
    return res;
}

Matrix mat_transpose(Matrix mat) {
    Matrix res;
    size_t i, j;

    res = mat_zeroes(mat.shape[1], mat.shape[0]);

    for(i = 0; i < res.shape[0]; i++) {
        for(j = 0; j < res.shape[1]; j++) {
            res.data[i][j] = mat.data[j][i];
        }
    }
    return res;
}