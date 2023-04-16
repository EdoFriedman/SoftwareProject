#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "matrix.h"
#include "spkmeans.h"

/**
 * Gets a Python list containing datapoints
 * Returns a C double** containing the datapoints and sets *datapoint_count and *dim accordingly.
 */
double** get_datapoints(PyObject* py_list, size_t* datapoint_count, size_t* dim) {
    PyObject* py_datapoint;
    double** datapoints;
    size_t i;
    size_t j;
    *dim = 0;
    *datapoint_count = PyObject_Length(py_list);
    datapoints = malloc(*datapoint_count * sizeof(double*));
    for(i = 0; i < *datapoint_count; i++) {
        py_datapoint = PyObject_GetItem(py_list, PyLong_FromUnsignedLong(i));
        if(!*dim) *dim = PyObject_Length(py_datapoint);
        datapoints[i] = malloc(*dim * sizeof(double));
        for(j = 0; j < *dim; j++) {
            datapoints[i][j] = PyFloat_AsDouble(PyObject_GetItem(py_datapoint, PyLong_FromUnsignedLong(j)));
        }
    }
    return datapoints;
}

void free_datapoints(double** datapoints, size_t datapoint_count) {
    size_t i;
    for(i = 0; i < datapoint_count; i++) {
        free(datapoints[i]);
    }
    free(datapoints);
}

PyObject* matrix_to_pylist(Matrix* matrix) {
    PyObject* py_matrix;
    PyObject* tmp_list;
    size_t i;
    size_t j;
    py_matrix = PyList_New(matrix->shape[0]);
    for(i = 0; i < matrix->shape[0]; i++) {
        tmp_list = PyList_New(matrix->shape[1]);
        for(j = 0; j < matrix->shape[1]; j++) {
            PyList_SetItem(tmp_list, j, PyFloat_FromDouble(matrix->data[i][j]));
        }
        PyList_SetItem(py_matrix, i, tmp_list);
    }
    return py_matrix;
}

static PyObject* spk_py(PyObject* self, PyObject* args) {
    PyObject* py_initial_centroids;
    PyObject* py_U;
    Matrix U;
    double** initial_centroids;
    size_t cluster_count;
    kmeans_input kmeansInput;
    Matrix kmeansOutput;
    PyObject* kmeansOutputPy;
    if(!PyArg_ParseTuple(args, "OO", &py_initial_centroids, &py_U)) {
        return NULL;
    }
    initial_centroids = get_datapoints(py_initial_centroids, &cluster_count, &U.shape[1]);
    U.data = get_datapoints(py_U, &U.shape[0], &U.shape[1]);

    kmeansInput.epsilon = 0;
    kmeansInput.iter_count = 300;
    kmeansInput.datapoints = U.data;
    kmeansInput.cluster_count = cluster_count;
    kmeansInput.dim = U.shape[1];
    kmeansInput.datapoint_count = U.shape[0];

    kmeansOutput.data = kmeans(kmeansInput, initial_centroids);
    kmeansOutput.shape[0] = cluster_count;
    kmeansOutput.shape[1] = cluster_count;
    free_datapoints(initial_centroids, cluster_count);
    mat_free(U);
    kmeansOutputPy = matrix_to_pylist(&kmeansOutput);
    mat_free(kmeansOutput);
    return kmeansOutputPy;
}

static PyObject* wam_py(PyObject* self, PyObject* args) {
    PyObject* py_datapoints;
    double** datapoints;
    size_t datapoint_count;
    size_t dim;
    Matrix wam;
    PyObject* res;
    if(!PyArg_ParseTuple(args, "O", &py_datapoints)) {
        return NULL;
    }
    datapoints = get_datapoints(py_datapoints, &datapoint_count, &dim);

    wam = get_w_adj_matrix(datapoints, datapoint_count, dim);
    res = matrix_to_pylist(&wam);

    mat_free(wam);
    free_datapoints(datapoints, datapoint_count);
    return res;
}

static PyObject* ddg_py(PyObject* self, PyObject* args) {
    PyObject* py_datapoints;
    double** datapoints;
    size_t datapoint_count;
    size_t dim;
    Matrix wam;
    Matrix ddg;
    PyObject* res;
    if(!PyArg_ParseTuple(args, "O", &py_datapoints)) {
        return NULL;
    }
    datapoints = get_datapoints(py_datapoints, &datapoint_count, &dim);

    wam = get_w_adj_matrix(datapoints, datapoint_count, dim);
    ddg = get_ddg_matrix(wam);
    res = matrix_to_pylist(&ddg);

    free_datapoints(datapoints, datapoint_count);
    mat_free(wam);
    mat_free(ddg);
    return res;
}

static PyObject* gl_py(PyObject* self, PyObject* args) {
    PyObject* py_datapoints;
    double** datapoints;
    size_t datapoint_count;
    size_t dim;
    Matrix gl;
    PyObject* res;
    if(!PyArg_ParseTuple(args, "O", &py_datapoints)) {
        return NULL;
    }
    datapoints = get_datapoints(py_datapoints, &datapoint_count, &dim);

    gl = get_laplacian_matrix(datapoints, datapoint_count, dim);
    res = matrix_to_pylist(&gl);

    free_datapoints(datapoints, datapoint_count);
    mat_free(gl);
    return res;
}

static PyObject* jacobi_py(PyObject* self, PyObject* args) {
    PyObject* py_matrix;
    Matrix matrix;
    eigen_struct jacobi_output;
    size_t i;
    PyObject* eigenvectors;
    PyObject* eigenvalues;
    if(!PyArg_ParseTuple(args, "O", &py_matrix)) {
        return NULL;
    }
    matrix.data = get_datapoints(py_matrix, &matrix.shape[0], &matrix.shape[1]);

    jacobi_output = jacobi(matrix);
    eigenvectors = matrix_to_pylist(&jacobi_output.eigenvectors);
    eigenvalues = PyList_New(jacobi_output.eigenvectors.shape[1]);
    for (i = 0; i < jacobi_output.eigenvectors.shape[1]; i++) {
        PyList_SetItem(eigenvalues, i, PyFloat_FromDouble(jacobi_output.eigenvalues[i]));
    }

    mat_free(jacobi_output.eigenvectors);
    free(jacobi_output.eigenvalues);
    mat_free(matrix);
    return PyTuple_Pack(2, eigenvalues, eigenvectors);
}

static PyMethodDef spkmeansMethods[] = {
    {"spk",
    (PyCFunction) spk_py,
    METH_VARARGS,
    PyDoc_STR("Performs full spectral kmeans.")},

    {"wam",
    (PyCFunction) wam_py,
    METH_VARARGS,
    PyDoc_STR("Calculates the weighted adjacency matrix.")},

    {"ddg",
    (PyCFunction) ddg_py,
    METH_VARARGS,
    PyDoc_STR("Calculates the diagonal degree matrix.")},

    {"gl",
    (PyCFunction) gl_py,
    METH_VARARGS,
    PyDoc_STR("Calculates the graph laplacian.")},

    {"jacobi",
    (PyCFunction) jacobi_py,
    METH_VARARGS,
    PyDoc_STR("Calculates the eigenvalues and eigenvectors of a real symmetric matrix using the Jacobi algorithm.")},

    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef spkmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    spkmeansMethods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&spkmeansmodule);
    if (!m) {
        return NULL;
    }
    return m;
}