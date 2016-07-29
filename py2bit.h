#include <Python.h>
#include "2bit.h"

typedef struct {
    PyObject_HEAD
    TwoBit *tb;
    int storeMasked; //Whether storeMasked was set. 0 = False, 1 = True
} pyTwoBit_t;

static PyObject* py2bitOpen(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject *py2bitInfo(pyTwoBit_t *pybw, PyObject *args);
static PyObject* py2bitClose(pyTwoBit_t *pybw, PyObject *args);
static PyObject* py2bitChroms(pyTwoBit_t *pybw, PyObject *args);
static PyObject *py2bitSequence(pyTwoBit_t *pybw, PyObject *args, PyObject *kwds);
static PyObject *py2bitBases(pyTwoBit_t *pybw, PyObject *args, PyObject *kwds);
static void py2bitDealloc(pyTwoBit_t *pybw);

static PyMethodDef tbMethods[] = {
    {"open", (PyCFunction)py2bitOpen, METH_VARARGS|METH_KEYWORDS,
"Open a 2bit file.\n\
\n\
Returns:\n\
   A TwoBit object on success, otherwise None.\n\
\n\
Arguments:\n\
    file: The name of a 2bit file.\n\
\n\
Optional arguments:\n\
    storeMasked: Whether to store information about soft-masking (default False).\n\
\n\
Note that storing soft-masking information can be memory intensive and doing so\n\
will result in soft-masked bases being lower case if the sequence is fetched\n\
(see the sequence() function)\n\
\n\
>>> import py2bit\n\
>>> tb = py2bit.open(\"some_file.2bit\")\n\
\n\
To store soft-masking information:\n\
>>> tb = py2bit.open(\"some_file.2bit\", True)"},
    {"info", (PyCFunction)py2bitInfo, METH_VARARGS,
"Returns a dictionary containing the following key:value pairs: \n\
\n\
  * The file size, in bytes ('file size').\n\
  * The number of chromosomes/contigs ('nChroms').\n\
  * The total sequence length ('sequence length').\n\
  * The total hard-masked length ('hard-masked length').\n\
  * The total soft-masked length, if available ('soft-masked length').\n\
\n\
A base is hard-masked if it is an N and soft-masked if it's lower case. Note that soft-masking is ignored by default (you must specify 'storeMasked=True' when you open the file.\n\
\n\
>>> import py2bit\n\
>>> tb = py2bit.open(\"some_file.2bit\")\n\
>>> tb.info()\n\
{'file size': 160L, 'nChroms': 2L, 'sequence length': 250L, 'hard-masked length': 150L, 'soft-masked length': 8L}\n\
>>> tb.close()\n"},
    {"close", (PyCFunction)py2bitClose, METH_VARARGS,
"Close a 2bit file.\n\
\n\
>>> import py2bit\n\
>>> tb = py2bit.open(\"some_file.2bit\")\n\
>>> tb.close()\n"},
    {"chroms", (PyCFunction)py2bitChroms, METH_VARARGS,
"Return a chromosome: length dictionary. The order is typically not\n\
alphabetical and the lengths are long (thus the 'L' suffix).\n\
\n\
Optional arguments:\n\
    chrom: An optional chromosome name\n\
\n\
Returns:\n\
    A list of chromosome lengths or a dictionary of them.\n\
\n\
>>> import py2bit\n\
>>> tb = py2bit.open(\"test/test.2bit\")\n\
>>> tb.chroms()\n\
{'chr1': 150L, 'chr2': 100L}\n\
\n\
Note that you may optionally supply a specific chromosome:\n\
\n\
>>> tb.chroms(\"chr1\")\n\
150L\n\
\n\
If you specify a non-existant chromosome then no output is produced:\n\
\n\
>>> tb.chroms(\"foo\")\n\
>>>\n"},
    {"sequence", (PyCFunction)py2bitSequence, METH_VARARGS|METH_KEYWORDS,
"Retrieve the sequence of a chromosome, or subset of it. On error, a runtime\n\
exception is thrown.\n\
\n\
Positional arguments:\n\
    chr:   Chromosome name\n\
\n\
Keyword arguments:\n\
    start: Starting position (0-based)\n\
    end:   Ending position (1-based)\n\
\n\
Returns:\n\
    A string containing the sequence.\n\
\n\
If start and end aren't specified, the entire chromosome is returned. If the\n\
end value is beyond the end of the chromosome then it is adjusted accordingly.\n\
\n\
>>> import py2bit\n\
>>> tb = py2bit.open(\"test/test.2bit\")\n\
>>> tb.sequence(\"chr1\")\n\
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATCGATCGTAGCTAGCTAGCTAGCTGATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n\
>>> tb.sequence(\"chr1\", 24, 74)\n\
NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATC\n\
>>> tb.close()"},
    {"bases", (PyCFunction)py2bitBases, METH_VARARGS|METH_KEYWORDS,
"Retrieve the percentage or number of A, C, T, and Gs in a chromosome or subset\n\
thereof. On error, a runtime exception is thrown.\n\
\n\
Positional arguments:\n\
    chr:   Chromosome name\n\
\n\
Optional keyword arguments:\n\
    start: Starting position (0-based)\n\
    end:   Ending position (1-based)\n\
    fraction: Whether to return fractional or integer values (default 'True',\n\
              so fractional values are returned)\n\
\n\
Returns:\n\
    A dictionary with nucleotide as the key and fraction (or count) as the\n\
    value.\n\
\n\
If start and end aren't specified, the entire chromosome is returned. If the\n\
end value is beyond the end of the chromosome then it is adjusted accordingly.\n\
\n\
Note that the fractions will sum to much less than 1 if there are hard-masked\n\
bases. Counts may sum to less than the length of the region for the same reason.\n\
\n\
>>> import py2bit\n\
>>> tb = py2bit.open(\"test/test.2bit\")\n\
>>> tb.bases(tb, \"chr1\")\n\
{'A': 0.08, 'C': 0.08, 'T': 0.08666666666666667, 'G': 0.08666666666666667}\n\
>>> tb.bases(tb, \"chr1\", 24, 74)\n\
{'A': 0.12, 'C': 0.12, 'T': 0.12, 'G': 0.12}\n\
>>> tb.bases(tb, \"chr1\", 24, 74, True)\n\
{'A': 6, 'C': 6, 'T': 6, 'G': 6}\n\
>>> tb.close()"},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
struct py2bitmodule_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct py2bitmodule_state*)PyModule_GetState(m))

static PyModuleDef py2bitmodule = {
    PyModuleDef_HEAD_INIT,
    "py2bit",
    "A python module for accessing 2bit files",
    -1,
    tbMethods,
    NULL, NULL, NULL, NULL
};
#endif

static PyTypeObject pyTwoBit = {
#if PY_MAJOR_VERSION >= 3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,              /*ob_size*/
#endif
    "py2bit.pyTwoBit",         /*tp_name*/
    sizeof(pyTwoBit),          /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)py2bitDealloc,     /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash*/
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    PyObject_GenericGetAttr, /*tp_getattro*/
    PyObject_GenericSetAttr, /*tp_setattro*/
    0,                         /*tp_as_buffer*/
#if PY_MAJOR_VERSION >= 3
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
#else
    Py_TPFLAGS_HAVE_CLASS,     /*tp_flags*/
#endif
    "bigWig File",             /*tp_doc*/
    0,                         /*tp_traverse*/
    0,                         /*tp_clear*/
    0,                         /*tp_richcompare*/
    0,                         /*tp_weaklistoffset*/
    0,                         /*tp_iter*/
    0,                         /*tp_iternext*/
    tbMethods,                 /*tp_methods*/
    0,                         /*tp_members*/
    0,                         /*tp_getset*/
    0,                         /*tp_base*/
    0,                         /*tp_dict*/
    0,                         /*tp_descr_get*/
    0,                         /*tp_descr_set*/
    0,                         /*tp_dictoffset*/
    0,                         /*tp_init*/
    0,                         /*tp_alloc*/
    0,                         /*tp_new*/
    0,0,0,0,0,0
};
