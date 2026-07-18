/*
 * ilreco._core — CPython extension wrapping the ilreco C library.
 *
 * One type is exported: _core.Calorimeter. It owns one ilreco_config and a
 * pool of ilreco_workspace objects. reconstruct_batch() releases the GIL for
 * the whole event loop and takes a workspace from the pool, so any number of
 * Python threads (or a free-threaded interpreter) may call the same
 * Calorimeter concurrently — each in-flight call holds its own workspace.
 *
 * The user-facing API (table dispatch, n_jobs, profile shortcuts) lives in
 * pure Python: python/ilreco/__init__.py.
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <pythread.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <ilreco.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Enough parallel calls for any realistic thread count; workspaces returned
 * to a full pool are destroyed instead of stored. */
#define POOL_CAPACITY 64

/* Upper bound on clusters read out per event. The library itself caps the
 * number of reconstructed objects per event below this. */
#define MAX_CLUSTERS_PER_EVENT 256

typedef struct {
    PyObject_HEAD
    ilreco_config *config;
    ilreco_workspace *pool[POOL_CAPACITY];
    int32_t pool_size;
    PyThread_type_lock pool_lock;
    int32_t started;   /* set on first reconstruct; config setters then fail */
    int32_t n_cols;
    int32_t n_rows;
} CalorimeterObject;

/* ------------------------- workspace pool -------------------------------- */

static ilreco_workspace *pool_acquire(CalorimeterObject *self) {
    ilreco_workspace *workspace = NULL;
    PyThread_acquire_lock(self->pool_lock, WAIT_LOCK);
    if (self->pool_size > 0) {
        self->pool_size -= 1;
        workspace = self->pool[self->pool_size];
    }
    PyThread_release_lock(self->pool_lock);
    if (!workspace) {
        workspace = ilreco_workspace_create(self->config);
    }
    return workspace;
}

static void pool_release(CalorimeterObject *self, ilreco_workspace *workspace) {
    PyThread_acquire_lock(self->pool_lock, WAIT_LOCK);
    if (self->pool_size < POOL_CAPACITY) {
        self->pool[self->pool_size] = workspace;
        self->pool_size += 1;
        workspace = NULL;
    }
    PyThread_release_lock(self->pool_lock);
    if (workspace) {
        ilreco_workspace_destroy(workspace);
    }
}

/* ------------------------- Calorimeter type ------------------------------ */

static int Calorimeter_init(CalorimeterObject *self, PyObject *args,
                            PyObject *kwargs) {
    static char *keyword_names[] = {"n_cols", "n_rows", "profile_path", NULL};
    int n_cols = 0;
    int n_rows = 0;
    const char *profile_path = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iis", keyword_names,
                                     &n_cols, &n_rows, &profile_path)) {
        return -1;
    }

    char error[256] = {0};
    self->config = ilreco_config_create(n_cols, n_rows, profile_path,
                                        error, sizeof error);
    if (!self->config) {
        PyErr_Format(PyExc_ValueError, "ilreco_config_create failed: %s", error);
        return -1;
    }
    self->pool_lock = PyThread_allocate_lock();
    if (!self->pool_lock) {
        ilreco_config_destroy(self->config);
        self->config = NULL;
        PyErr_NoMemory();
        return -1;
    }
    self->pool_size = 0;
    self->started = 0;
    self->n_cols = n_cols;
    self->n_rows = n_rows;
    return 0;
}

static void Calorimeter_dealloc(CalorimeterObject *self) {
    for (int32_t i = 0; i < self->pool_size; ++i) {
        ilreco_workspace_destroy(self->pool[i]);
    }
    if (self->pool_lock) {
        PyThread_free_lock(self->pool_lock);
    }
    ilreco_config_destroy(self->config);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

/* Config setters are only legal before the first reconstruction: after that
 * the config may be shared by concurrent calls and must stay immutable. */
static int require_not_started(CalorimeterObject *self, const char *what) {
    if (self->started) {
        PyErr_Format(PyExc_RuntimeError,
                     "%s must be called before the first reconstruction", what);
        return -1;
    }
    return 0;
}

static PyObject *Calorimeter_set_cell_mask(CalorimeterObject *self,
                                           PyObject *mask_object) {
    if (require_not_started(self, "set_cell_mask") < 0) {
        return NULL;
    }
    PyArrayObject *mask = (PyArrayObject *)PyArray_FROM_OTF(
        mask_object, NPY_UINT8, NPY_ARRAY_IN_ARRAY);
    if (!mask) {
        return NULL;
    }
    const npy_intp expected = (npy_intp)self->n_cols * self->n_rows;
    if (PyArray_SIZE(mask) != expected) {
        Py_DECREF(mask);
        PyErr_Format(PyExc_ValueError,
                     "cell mask must have n_rows*n_cols = %ld entries",
                     (long)expected);
        return NULL;
    }
    const int status = ilreco_config_set_cell_mask(
        self->config, (const unsigned char *)PyArray_DATA(mask));
    Py_DECREF(mask);
    if (status != 0) {
        PyErr_SetString(PyExc_RuntimeError, "ilreco_config_set_cell_mask failed");
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject *Calorimeter_set_hole_classification(CalorimeterObject *self,
                                                     PyObject *enabled_object) {
    if (require_not_started(self, "set_hole_classification") < 0) {
        return NULL;
    }
    const long enabled = PyLong_AsLong(enabled_object);
    if (enabled == -1 && PyErr_Occurred()) {
        return NULL;
    }
    ilreco_config_set_hole_classification(self->config, (int32_t)enabled);
    Py_RETURN_NONE;
}

static PyObject *Calorimeter_set_seed_threshold(CalorimeterObject *self,
                                                PyObject *threshold_object) {
    if (require_not_started(self, "set_seed_threshold") < 0) {
        return NULL;
    }
    const double threshold_gev = PyFloat_AsDouble(threshold_object);
    if (threshold_gev == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    ilreco_config_set_seed_threshold(self->config, threshold_gev);
    Py_RETURN_NONE;
}

static PyObject *Calorimeter_set_zcal(CalorimeterObject *self,
                                      PyObject *zcal_object) {
    if (require_not_started(self, "set_zcal") < 0) {
        return NULL;
    }
    const double zcal_cm = PyFloat_AsDouble(zcal_object);
    if (zcal_cm == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    ilreco_config_set_zcal(self->config, zcal_cm);
    Py_RETURN_NONE;
}

/* Growable per-call output buffers (plain C memory; copied into numpy arrays
 * only after the GIL is re-acquired). */
typedef struct {
    int64_t *event_index;
    double *e;
    double *x;
    double *y;
    double *chi2;
    int32_t *size;
    int32_t *type;
    npy_intp count;
    npy_intp capacity;
} cluster_buffer;

static int cluster_buffer_reserve(cluster_buffer *buffer, npy_intp capacity) {
    if (capacity <= buffer->capacity) {
        return 0;
    }
    int64_t *event_index = realloc(buffer->event_index, capacity * sizeof(int64_t));
    double *e = realloc(buffer->e, capacity * sizeof(double));
    double *x = realloc(buffer->x, capacity * sizeof(double));
    double *y = realloc(buffer->y, capacity * sizeof(double));
    double *chi2 = realloc(buffer->chi2, capacity * sizeof(double));
    int32_t *size = realloc(buffer->size, capacity * sizeof(int32_t));
    int32_t *type = realloc(buffer->type, capacity * sizeof(int32_t));
    if (event_index) buffer->event_index = event_index;
    if (e) buffer->e = e;
    if (x) buffer->x = x;
    if (y) buffer->y = y;
    if (chi2) buffer->chi2 = chi2;
    if (size) buffer->size = size;
    if (type) buffer->type = type;
    if (!event_index || !e || !x || !y || !chi2 || !size || !type) {
        return -1;
    }
    buffer->capacity = capacity;
    return 0;
}

static void cluster_buffer_free(cluster_buffer *buffer) {
    free(buffer->event_index);
    free(buffer->e);
    free(buffer->x);
    free(buffer->y);
    free(buffer->chi2);
    free(buffer->size);
    free(buffer->type);
}

/* Copy one C column into a fresh 1-D numpy array. */
static PyObject *column_to_array(const void *data, npy_intp count,
                                 int numpy_type, size_t element_size) {
    PyObject *array = PyArray_SimpleNew(1, &count, numpy_type);
    if (array) {
        memcpy(PyArray_DATA((PyArrayObject *)array), data, count * element_size);
    }
    return array;
}

/*
 * reconstruct_batch(cols, rows, energies, starts)
 *
 *   cols, rows  int32 arrays of 0-based cell indices, all events concatenated
 *   energies    float64 array [GeV], parallel to cols/rows
 *   starts      int64 offsets, length n_events+1: event i owns rows
 *               starts[i] .. starts[i+1]-1
 *
 * Returns a tuple of per-cluster columns:
 *   (event_index int64, e f8, x f8, y f8, chi2 f8, size i4, type i4)
 * where event_index is the position of the event within THIS batch. The
 * Python layer maps it back to the caller's event ids.
 */
static PyObject *Calorimeter_reconstruct_batch(CalorimeterObject *self,
                                               PyObject *args) {
    PyObject *cols_object = NULL;
    PyObject *rows_object = NULL;
    PyObject *energies_object = NULL;
    PyObject *starts_object = NULL;
    if (!PyArg_ParseTuple(args, "OOOO", &cols_object, &rows_object,
                          &energies_object, &starts_object)) {
        return NULL;
    }

    PyArrayObject *cols = (PyArrayObject *)PyArray_FROM_OTF(
        cols_object, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *rows = (PyArrayObject *)PyArray_FROM_OTF(
        rows_object, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *energies = (PyArrayObject *)PyArray_FROM_OTF(
        energies_object, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *starts = (PyArrayObject *)PyArray_FROM_OTF(
        starts_object, NPY_INT64, NPY_ARRAY_IN_ARRAY);
    if (!cols || !rows || !energies || !starts) {
        goto error_arrays;
    }

    const npy_intp n_hits = PyArray_SIZE(cols);
    const npy_intp n_events = PyArray_SIZE(starts) - 1;
    if (PyArray_SIZE(rows) != n_hits || PyArray_SIZE(energies) != n_hits ||
        n_events < 0) {
        PyErr_SetString(PyExc_ValueError,
                        "cols/rows/energies must have equal length and starts "
                        "must have n_events+1 entries");
        goto error_arrays;
    }

    const int32_t *col_data = (const int32_t *)PyArray_DATA(cols);
    const int32_t *row_data = (const int32_t *)PyArray_DATA(rows);
    const double *energy_data = (const double *)PyArray_DATA(energies);
    const int64_t *start_data = (const int64_t *)PyArray_DATA(starts);

    /* validate offsets and find the largest event (scratch buffer size) */
    npy_intp max_event_hits = 0;
    for (npy_intp event = 0; event < n_events; ++event) {
        const int64_t begin = start_data[event];
        const int64_t end = start_data[event + 1];
        if (begin < 0 || end < begin || end > n_hits) {
            PyErr_Format(PyExc_ValueError,
                         "starts is not a valid offset array at event %ld",
                         (long)event);
            goto error_arrays;
        }
        if (end - begin > max_event_hits) {
            max_event_hits = end - begin;
        }
    }

    ilreco_hit *hit_scratch = NULL;
    if (max_event_hits > 0) {
        hit_scratch = malloc(max_event_hits * sizeof(ilreco_hit));
        if (!hit_scratch) {
            PyErr_NoMemory();
            goto error_arrays;
        }
    }

    cluster_buffer buffer = {0};
    if (cluster_buffer_reserve(&buffer, n_events + 64) < 0) {
        free(hit_scratch);
        PyErr_NoMemory();
        goto error_arrays;
    }

    ilreco_workspace *workspace = pool_acquire(self);
    if (!workspace) {
        cluster_buffer_free(&buffer);
        free(hit_scratch);
        PyErr_NoMemory();
        goto error_arrays;
    }
    self->started = 1;

    npy_intp failed_event = -1;      /* invalid input reported by the library */
    int out_of_memory = 0;
    ilreco_cluster event_clusters[MAX_CLUSTERS_PER_EVENT];

    Py_BEGIN_ALLOW_THREADS
    for (npy_intp event = 0; event < n_events; ++event) {
        const int64_t begin = start_data[event];
        const int32_t event_hits = (int32_t)(start_data[event + 1] - begin);
        for (int32_t i = 0; i < event_hits; ++i) {
            hit_scratch[i].col = col_data[begin + i];
            hit_scratch[i].row = row_data[begin + i];
            hit_scratch[i].e = energy_data[begin + i];
        }
        const int32_t n_found = ilreco_reconstruct(workspace, hit_scratch,
                                                   event_hits, event_clusters,
                                                   MAX_CLUSTERS_PER_EVENT);
        if (n_found < 0) {
            failed_event = event;
            break;
        }
        const int32_t n_keep = n_found < MAX_CLUSTERS_PER_EVENT
                                   ? n_found
                                   : MAX_CLUSTERS_PER_EVENT;
        if (buffer.count + n_keep > buffer.capacity &&
            cluster_buffer_reserve(&buffer, 2 * buffer.capacity + n_keep) < 0) {
            out_of_memory = 1;
            break;
        }
        for (int32_t k = 0; k < n_keep; ++k) {
            buffer.event_index[buffer.count] = event;
            buffer.e[buffer.count] = event_clusters[k].e;
            buffer.x[buffer.count] = event_clusters[k].x;
            buffer.y[buffer.count] = event_clusters[k].y;
            buffer.chi2[buffer.count] = event_clusters[k].chi2;
            buffer.size[buffer.count] = event_clusters[k].size;
            buffer.type[buffer.count] = event_clusters[k].type;
            buffer.count += 1;
        }
    }
    Py_END_ALLOW_THREADS

    pool_release(self, workspace);
    free(hit_scratch);

    if (failed_event >= 0) {
        PyErr_Format(PyExc_ValueError,
                     "invalid hit (outside the %dx%d grid or on a masked-out "
                     "cell) in batch event %ld",
                     self->n_cols, self->n_rows, (long)failed_event);
        cluster_buffer_free(&buffer);
        goto error_arrays;
    }
    if (out_of_memory) {
        PyErr_NoMemory();
        cluster_buffer_free(&buffer);
        goto error_arrays;
    }

    PyObject *result = PyTuple_New(7);
    PyObject *event_column = column_to_array(buffer.event_index, buffer.count,
                                             NPY_INT64, sizeof(int64_t));
    PyObject *e_column = column_to_array(buffer.e, buffer.count,
                                         NPY_FLOAT64, sizeof(double));
    PyObject *x_column = column_to_array(buffer.x, buffer.count,
                                         NPY_FLOAT64, sizeof(double));
    PyObject *y_column = column_to_array(buffer.y, buffer.count,
                                         NPY_FLOAT64, sizeof(double));
    PyObject *chi2_column = column_to_array(buffer.chi2, buffer.count,
                                            NPY_FLOAT64, sizeof(double));
    PyObject *size_column = column_to_array(buffer.size, buffer.count,
                                            NPY_INT32, sizeof(int32_t));
    PyObject *type_column = column_to_array(buffer.type, buffer.count,
                                            NPY_INT32, sizeof(int32_t));
    cluster_buffer_free(&buffer);
    if (!result || !event_column || !e_column || !x_column || !y_column ||
        !chi2_column || !size_column || !type_column) {
        Py_XDECREF(result);
        Py_XDECREF(event_column);
        Py_XDECREF(e_column);
        Py_XDECREF(x_column);
        Py_XDECREF(y_column);
        Py_XDECREF(chi2_column);
        Py_XDECREF(size_column);
        Py_XDECREF(type_column);
        goto error_arrays;
    }
    PyTuple_SET_ITEM(result, 0, event_column);
    PyTuple_SET_ITEM(result, 1, e_column);
    PyTuple_SET_ITEM(result, 2, x_column);
    PyTuple_SET_ITEM(result, 3, y_column);
    PyTuple_SET_ITEM(result, 4, chi2_column);
    PyTuple_SET_ITEM(result, 5, size_column);
    PyTuple_SET_ITEM(result, 6, type_column);

    Py_DECREF(cols);
    Py_DECREF(rows);
    Py_DECREF(energies);
    Py_DECREF(starts);
    return result;

error_arrays:
    Py_XDECREF(cols);
    Py_XDECREF(rows);
    Py_XDECREF(energies);
    Py_XDECREF(starts);
    return NULL;
}

static PyObject *Calorimeter_get_n_cols(CalorimeterObject *self, void *closure) {
    (void)closure;
    return PyLong_FromLong(self->n_cols);
}

static PyObject *Calorimeter_get_n_rows(CalorimeterObject *self, void *closure) {
    (void)closure;
    return PyLong_FromLong(self->n_rows);
}

static PyMethodDef Calorimeter_methods[] = {
    {"reconstruct_batch", (PyCFunction)Calorimeter_reconstruct_batch,
     METH_VARARGS,
     "reconstruct_batch(cols, rows, energies, starts) -> per-cluster column "
     "tuple; see the module source for the exact contract"},
    {"set_cell_mask", (PyCFunction)Calorimeter_set_cell_mask, METH_O,
     "set_cell_mask(mask_uint8_row_major) — nonzero = cell exists; must be "
     "called before the first reconstruction"},
    {"set_hole_classification", (PyCFunction)Calorimeter_set_hole_classification,
     METH_O, "enable/disable built-in central-2x2 hole labeling"},
    {"set_seed_threshold", (PyCFunction)Calorimeter_set_seed_threshold, METH_O,
     "cluster seed threshold in GeV"},
    {"set_zcal", (PyCFunction)Calorimeter_set_zcal, METH_O,
     "target-to-calorimeter distance in cm (two-gamma separation cut)"},
    {NULL, NULL, 0, NULL}};

static PyGetSetDef Calorimeter_getset[] = {
    {"n_cols", (getter)Calorimeter_get_n_cols, NULL, "grid columns", NULL},
    {"n_rows", (getter)Calorimeter_get_n_rows, NULL, "grid rows", NULL},
    {NULL, NULL, NULL, NULL, NULL}};

static PyTypeObject CalorimeterType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "ilreco._core.Calorimeter",
    .tp_basicsize = sizeof(CalorimeterObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = "Low-level calorimeter context: config + workspace pool. "
              "Use ilreco.Calorimeter (the Python wrapper) instead.",
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc)Calorimeter_init,
    .tp_dealloc = (destructor)Calorimeter_dealloc,
    .tp_methods = Calorimeter_methods,
    .tp_getset = Calorimeter_getset,
};

static PyModuleDef core_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "ilreco._core",
    .m_doc = "C extension core of the ilreco Python binding",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit__core(void) {
    import_array();
    if (PyType_Ready(&CalorimeterType) < 0) {
        return NULL;
    }
    PyObject *module = PyModule_Create(&core_module);
    if (!module) {
        return NULL;
    }
    Py_INCREF(&CalorimeterType);
    if (PyModule_AddObject(module, "Calorimeter",
                           (PyObject *)&CalorimeterType) < 0) {
        Py_DECREF(&CalorimeterType);
        Py_DECREF(module);
        return NULL;
    }
    return module;
}
