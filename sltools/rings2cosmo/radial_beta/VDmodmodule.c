#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "VDIntegral.h"

static PyObject* VDmod_integral(PyObject* self, PyObject *args){
    double b, alpha, delta, ra; 
    double VD = 0;

    if(!PyArg_ParseTuple(args,"dddd", &b, &alpha, &delta, &ra)){
        return NULL; //return error if none found
    }

    VD = Integral_func (b,  alpha, delta, ra); //Integral_func (double b, double alpha, double delta, double ra)


    return Py_BuildValue("d",VD);
}

static PyMethodDef VDmod_methods[] = {
    {"integralVD",  VDmod_integral, METH_VARARGS, "Returns the result for the integral part of the vd"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef VDmod = {
    PyModuleDef_HEAD_INIT,
    "VDmod",   /* name of module */
    "module for sove the integral for VD of ETGs", /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    VDmod_methods
};

PyMODINIT_FUNC PyInit_VDmod(void)
{
    return PyModule_Create(&VDmod);
}