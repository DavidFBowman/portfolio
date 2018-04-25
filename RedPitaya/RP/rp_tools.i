/*rp_tools.i */
%module rp_tools
%{
#define SWIG_FILE_WITH_INIT
#include <python2.7/Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#include "./rp.h"
#include "./rp_tools.h"

extern float *acquire_buffer();
extern void Initialize ();
extern void Finalize ();
extern void blinker();

%}

extern void Initialize ();
extern void Finalize ();

extern void sine_wave();
extern void stop_signal();
extern void blinker();

//extern float *acquire_buffer();

%typemap(out) float* acquire_buffer{
  int i;
  //$1, $1_dim0, $1_dim1
  $result = PyList_New(16384);
  for (i = 0; i < 16384; i++) {
    PyObject *o = PyFloat_FromDouble((double) $1[i]);
    PyList_SetItem($result,i,o);
  }
}

%include "rp_tools.c"