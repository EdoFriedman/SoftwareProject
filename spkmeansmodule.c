#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject* spk_py(PyObject* self, PyObject* args) {
    return NULL;
}

static PyObject* wam_py(PyObject* self, PyObject* args) {
    return NULL;
}

static PyObject* ddg_py(PyObject* self, PyObject* args) {
    return NULL;
}

static PyObject* gl_py(PyObject* self, PyObject* args) {
    return NULL;
}

static PyObject* jacobi_py(PyObject* self, PyObject* args) {
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