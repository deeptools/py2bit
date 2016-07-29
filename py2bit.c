#include <Python.h>
#include <inttypes.h>
#include "py2bit.h"

static PyObject *py2bitOpen(PyObject *self, PyObject *args, PyObject *kwds) {
    char *fname = NULL;
    PyObject *storeMaskedO = Py_False;
    pyTwoBit_t *pytb;
    int storeMasked = 0;
    TwoBit *tb = NULL;
    static char *kwd_list[] = {"fname", "storeMasked", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|O", kwd_list, &fname, &storeMaskedO)) goto error;

    if(storeMaskedO == Py_True) storeMasked = 1;

    //Open the file
    tb = twobitOpen(fname, storeMasked);
    if(!tb) goto error;

    pytb = PyObject_New(pyTwoBit_t, &pyTwoBit);
    if(!pytb) goto error;
    pytb->storeMasked = storeMasked;
    pytb->tb = tb;

    return (PyObject*) pytb;

error:
    if(tb) twobitClose(tb);
    PyErr_SetString(PyExc_RuntimeError, "Received an error during file opening!");
    return NULL;
}

static void py2bitDealloc(pyTwoBit_t *self) {
    if(self->tb) twobitClose(self->tb);
    PyObject_DEL(self);
}

static PyObject *py2bitClose(pyTwoBit_t *self, PyObject *args) {
    if(self->tb) twobitClose(self->tb);
    self->tb = NULL;
    Py_INCREF(Py_None);
    return Py_None;
}

//Returns the file size, number of chromosomes/contigs, total sequence length and total masked length
static PyObject *py2bitInfo(pyTwoBit_t *self, PyObject *args) {
    TwoBit *tb = self->tb;
    PyObject *ret = NULL, *val = NULL;
    uint32_t i, j, foo;

    ret = PyDict_New();

    //file size
    val = PyLong_FromUnsignedLongLong(tb->sz);
    if(!val) goto error;
    if(PyDict_SetItemString(ret, "file size", val) == -1) goto error;
    Py_DECREF(val);

    //nContigs
    val = PyLong_FromUnsignedLong(tb->hdr->nChroms);
    if(!val) goto error;
    if(PyDict_SetItemString(ret, "nChroms", val) == -1) goto error;
    Py_DECREF(val);

    //sequence length
    foo = 0;
    for(i=0; i<tb->hdr->nChroms; i++) foo += tb->idx->size[i];
    val = PyLong_FromUnsignedLong(foo);
    if(!val) goto error;
    if(PyDict_SetItemString(ret, "sequence length", val) == -1) goto error;
    Py_DECREF(val);

    //hard-masked length
    foo = 0;
    for(i=0; i<tb->hdr->nChroms; i++) {
        for(j=0; j<tb->idx->nBlockCount[i]; j++) {
            foo += tb->idx->nBlockSizes[i][j];
        }
    }
    val = PyLong_FromUnsignedLong(foo);
    if(!val) goto error;
    if(PyDict_SetItemString(ret, "hard-masked length", val) == -1) goto error;
    Py_DECREF(val);

    //soft-masked length
    if(tb->idx->maskBlockStart) {
        foo = 0;
        for(i=0; i<tb->hdr->nChroms; i++) {
            for(j=0; j<tb->idx->maskBlockCount[i]; j++) {
                foo += tb->idx->maskBlockSizes[i][j];
            }
        }

        val = PyLong_FromUnsignedLong(foo);
        if(!val) goto error;
        if(PyDict_SetItemString(ret, "soft-masked length", val) == -1) goto error;
        Py_DECREF(val);
    }

    return ret;

error:
    Py_XDECREF(val);
    Py_XDECREF(ret);
    PyErr_SetString(PyExc_RuntimeError, "Received an error while gathering information on the 2bit file!");
    return NULL;
}

static PyObject *py2bitChroms(pyTwoBit_t *self, PyObject *args) {
    PyObject *ret = NULL, *val = NULL;
    TwoBit *tb = self->tb;
    char *chrom = NULL;
    uint32_t i;

    if(!(PyArg_ParseTuple(args, "|s", &chrom)) || !chrom) {
        ret = PyDict_New();
        if(!ret) goto error;
        for(i=0; i<tb->hdr->nChroms; i++) {
            val = PyLong_FromUnsignedLong(tb->idx->size[i]);
            if(!val) goto error;
            if(PyDict_SetItemString(ret, tb->cl->chrom[i], val) == -1) goto error;
            Py_DECREF(val);
        }
    } else {
        for(i=0; i<tb->hdr->nChroms; i++) {
            if(strcmp(tb->cl->chrom[i], chrom) == 0) {
                ret = PyLong_FromUnsignedLong(tb->idx->size[i]);
                if(!ret) goto error;
                break;
            }
        }
    }

    if(!ret) {
        Py_INCREF(Py_None);
        ret = Py_None;
    }

    return ret;

error :
    Py_XDECREF(val);
    Py_XDECREF(ret);
    PyErr_SetString(PyExc_RuntimeError, "Received an error while adding an item to the output dictionary!");
    return NULL;
}

#if PY_MAJOR_VERSION >= 3
PyObject *PyString_FromString(char *seq) {
    return PyUnicode_FromString(seq);
}
#endif

static PyObject *py2bitSequence(pyTwoBit_t *self, PyObject *args, PyObject *kwds) {
    PyObject *ret = NULL;
    TwoBit *tb = self->tb;
    char *seq, *chrom;
    unsigned long startl = 0, endl = 0;
    uint32_t start, end, len;
    static char *kwd_list[] = {"chrom", "start", "end", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|kk", kwd_list, &chrom, &startl, &endl)) {
        PyErr_SetString(PyExc_RuntimeError, "You must supply at least a chromosome!");
        return NULL;
    }

    len = twobitChromLen(tb, chrom);
    if(len == 0) {
        PyErr_SetString(PyExc_RuntimeError, "The specified chromosome doesn't exist in the 2bit file!");
        return NULL;
    }
    if(endl > len) endl = len;
    end = (uint32_t) endl;
    if(startl >= endl && startl > 0) {
        PyErr_SetString(PyExc_RuntimeError, "The start value must be less then the end value (and the end of the chromosome");
        return NULL;
    }
    start = (uint32_t) startl;
    seq = twobitSequence(tb, chrom, start, end);
    if(!seq) {
        PyErr_SetString(PyExc_RuntimeError, "There was an error while fetching the sequence!");
        return NULL;
    }

    ret = PyString_FromString(seq);
    free(seq);
    if(!ret) {
        PyErr_SetString(PyExc_RuntimeError, "Received an error while converting the C-level char array to a python string!");
        return NULL;
    }

    return ret;
}

static PyObject *py2bitBases(pyTwoBit_t *self, PyObject *args, PyObject *kwds) {
    PyObject *ret = NULL, *val = NULL;
    PyObject *fractionO = Py_True;
    TwoBit *tb = self->tb;
    char *chrom;
    void *o = NULL;
    unsigned long startl = 0, endl = 0;
    uint32_t start, end, len;
    static char *kwd_list[] = {"chrom", "start", "end", "fraction", NULL};
    int fraction = 1;

    if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|kkO", kwd_list, &chrom, &startl, &endl, &fractionO)) {
        PyErr_SetString(PyExc_RuntimeError, "You must supply at least a chromosome!");
        return NULL;
    }

    len = twobitChromLen(tb, chrom);
    if(len == 0) {
        PyErr_SetString(PyExc_RuntimeError, "The specified chromosome doesn't exist in the 2bit file!");
        return NULL;
    }
    if(endl > len) endl = len;
    end = (uint32_t) endl;
    if(startl >= endl && startl > 0) {
        PyErr_SetString(PyExc_RuntimeError, "The start value must be less then the end value (and the end of the chromosome");
        return NULL;
    }
    start = (uint32_t) startl;

    if(fractionO == Py_False) fraction = 0;

    o = twobitBases(tb, chrom, start, end, fraction);
    if(!o) {
        PyErr_SetString(PyExc_RuntimeError, "Received an error while determining the per-base metrics.");
        return NULL;
    }

    ret = PyDict_New();
    if(!ret) goto error;

    //A
    if(fraction) val = PyFloat_FromDouble(((double*)o)[0]);
    else val = PyLong_FromUnsignedLong(((uint32_t*)o)[0]);
    if(!val) goto error;
    if(PyDict_SetItemString(ret, "A", val) == -1) goto error;
    Py_DECREF(val);

    //C
    if(fraction) val = PyFloat_FromDouble(((double*)o)[1]);
    else val = PyLong_FromUnsignedLong(((uint32_t*)o)[1]);
    if(!val) goto error;
    if(PyDict_SetItemString(ret, "C", val) == -1) goto error;
    Py_DECREF(val);

    //T
    if(fraction) val = PyFloat_FromDouble(((double*)o)[2]);
    else val = PyLong_FromUnsignedLong(((uint32_t*)o)[2]);
    if(!val) goto error;
    if(PyDict_SetItemString(ret, "T", val) == -1) goto error;
    Py_DECREF(val);

    //G
    if(fraction) val = PyFloat_FromDouble(((double*)o)[3]);
    else val = PyLong_FromUnsignedLong(((uint32_t*)o)[3]);
    if(!val) goto error;
    if(PyDict_SetItemString(ret, "G", val) == -1) goto error;
    Py_DECREF(val);

    free(o);

    return ret;

error:
    if(o) free(o);
    Py_XDECREF(ret);
    Py_XDECREF(val);
    PyErr_SetString(PyExc_RuntimeError, "Received an error while constructing the output dictionary!");
    return NULL;
}

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_py2bit(void) {
    PyObject *res;

    if(PyType_Ready(&pyTwoBit) < 0) return NULL;
    res = PyModule_Create(&py2bitmodule);
    if(!res) return NULL;

    Py_INCREF(&pyTwoBit);
    PyModule_AddObject(res, "py2bit", (PyObject *) &pyTwoBit);

    return res;
}
#else
//Python2 initialization
PyMODINIT_FUNC initpy2bit(void) {
    if(PyType_Ready(&pyTwoBit) < 0) return;
    Py_InitModule3("py2bit", tbMethods, "A module for handling 2bit files");
}
#endif
