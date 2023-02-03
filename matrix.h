#include <stdlib.h>

typedef struct {
    size_t shape[2];
    double** data;
} Matrix;

/**
 * Returns a matrix filled with zeroes  .
 */
Matrix mat_zeroes(size_t rows, size_t columns);

/**
 * Returns an identity matrix.
 */
Matrix mat_id(size_t size);

/**
 * Returns a matrix where the value at position i,j is func(i,j).
 */
Matrix mat_from_func(size_t rows, size_t columns, double (*func)(double, double));

/**
 * Returns a copy of the given matrix.
 */
Matrix mat_clone(Matrix mat);

/**
 * Frees the memory used by the given matrix.
 * Do not use the matrix after calling this function.
 */
void mat_free(Matrix mat);

/**
 * Multiplies all elements in the given matrix by scalar and returns the result.
 * The original matrix remains unchanged.
 */
Matrix mat_mul_s(Matrix mat, double scalar);

/**
 * Multiplies the two matrices together and returns the result.
 */
Matrix mat_mul_m(Matrix mat1, Matrix mat2);

/**
 * Adds the two matrices together and returns the result.
 */
Matrix mat_add(Matrix mat1, Matrix mat2);

/**
 * Raises all of the matrices elements to the given power.
 */
Matrix mat_elem_pow(Matrix mat, double power);

/**
 * Sums all the elements of the matrix's columns and returns a 1 x n sized matrix containing the result
 */
Matrix mat_sum_cols(Matrix mat);

/**
 * Sums all the elements of the matrix's rows and returns an m x 1 sized matrix containing the result
 */
Matrix mat_sum_rows(Matrix mat);

/**
 * Returns the transpose of the given matrix
 */
Matrix mat_transpose(Matrix mat);