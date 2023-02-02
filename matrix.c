#include <stdlib.h>

typedef struct {
    size_t shape[2];
    double** data;
} Matrix;

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