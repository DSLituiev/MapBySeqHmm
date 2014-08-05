/*
 * cumMatrSafeF_mexutil.c
 *
 * Code generation for function 'cumMatrSafeF_mexutil'
 *
 * C source code generated on: Wed Jun 11 19:52:42 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cumMatrSafeF.h"
#include "cumMatrSafeF_mexutil.h"

/* Function Definitions */
void error(const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 0, NULL, 1, &pArray, "error", TRUE,
                        location);
}

const mxArray *message(const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  const mxArray *m4;
  pArray = b;
  return emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 1, &m4, 1, &pArray, "message",
    TRUE, location);
}

/* End of code generation (cumMatrSafeF_mexutil.c) */
