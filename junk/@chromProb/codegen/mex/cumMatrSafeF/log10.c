/*
 * log10.c
 *
 * Code generation for function 'log10'
 *
 * C source code generated on: Wed Jun 11 19:52:42 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cumMatrSafeF.h"
#include "log10.h"

/* Function Definitions */

/*
 *
 */
void b_log10(emxArray_real_T *x)
{
  int32_T i4;
  int32_T k;
  i4 = x->size[1];
  for (k = 0; k < i4; k++) {
    x->data[k] = muDoubleScalarLog10(x->data[k]);
  }
}

/* End of code generation (log10.c) */
