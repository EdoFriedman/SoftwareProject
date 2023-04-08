#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "matrix.h"

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
        py_datapoint = PyObject_GetItem(py_list, i);
        if(!*dim) *dim = PyObject_Length(py_datapoint);
        datapoints[i] = malloc(*dim * sizeof(double));
        for(j = 0; j < *dim; j++) {
            datapoints[i][j] = PyFloat_AsDouble(PyObject_GetItem(py_datapoint, j));
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

static PyObject* spk_py(PyObject* self, PyObject* args) {
    PyObject* py_initial_centroids;
    PyObject* py_datapoints;
    double** initial_centroids;
    size_t cluster_count;
    size_t dim;
    double** datapoints;
    size_t datapoint_count;
    if(!PyArg_ParseTuple(args, "OO", &py_initial_centroids, &py_datapoints)) {
        return NULL;
    }
    initial_centroids = get_datapoints(py_initial_centroids, &cluster_count, &dim);
    datapoints = get_datapoints(py_datapoints, &datapoint_count, &dim);


    free_datapoints(initial_centroids, cluster_count);
    free_datapoints(datapoints, datapoint_count);
    return NULL;
}

static PyObject* wam_py(PyObject* self, PyObject* args) {
    PyObject* py_datapoints;
    double** datapoints;
    size_t datapoint_count;
    size_t dim;
    if(!PyArg_ParseTuple(args, "O", &py_datapoints)) {
        return NULL;
    }
    datapoints = get_datapoints(py_datapoints, &datapoint_count, &dim);

    free_datapoints(datapoints, datapoint_count);
    return NULL;
}

static PyObject* ddg_py(PyObject* self, PyObject* args) {
    PyObject* py_datapoints;
    double** datapoints;
    size_t datapoint_count;
    size_t dim;
    if(!PyArg_ParseTuple(args, "O", &py_datapoints)) {
        return NULL;
    }
    datapoints = get_datapoints(py_datapoints, &datapoint_count, &dim);

    free_datapoints(datapoints, datapoint_count);
    return NULL;
}

static PyObject* gl_py(PyObject* self, PyObject* args) {
    PyObject* py_datapoints;
    double** datapoints;
    size_t datapoint_count;
    size_t dim;
    if(!PyArg_ParseTuple(args, "O", &py_datapoints)) {
        return NULL;
    }
    datapoints = get_datapoints(py_datapoints, &datapoint_count, &dim);

    free_datapoints(datapoints, datapoint_count);
    return NULL;
}

static PyObject* jacobi_py(PyObject* self, PyObject* args) {
    PyObject* py_matrix;
    Matrix matrix;
    if(!PyArg_ParseTuple(args, "O", &py_matrix)) {
        return NULL;
    }
    matrix.data = get_datapoints(py_matrix, &matrix.shape[0], &matrix.shape[1]);

    mat_free(matrix);
    return NULL;
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