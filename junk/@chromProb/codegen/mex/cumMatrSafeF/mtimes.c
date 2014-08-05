/*
 * mtimes.c
 *
 * Code generation for function 'mtimes'
 *
 * C source code generated on: Wed Jun 11 19:52:42 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cumMatrSafeF.h"
#include "mtimes.h"
#include "cumMatrSafeF_mexutil.h"

/* Variable Definitions */
static emlrtMCInfo emlrtMCI = { 94, 13, "mtimes",
  "/usr/local/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo b_emlrtMCI = { 93, 23, "mtimes",
  "/usr/local/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo c_emlrtMCI = { 99, 13, "mtimes",
  "/usr/local/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo d_emlrtMCI = { 98, 23, "mtimes",
  "/usr/local/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

/* Function Definitions */

/*
 *
 */
void b_dynamic_size_checks(const emxArray_real_T *a, const emxArray_real_T *b)
{
  const mxArray *y;
  static const int32_T iv6[2] = { 1, 45 };

  const mxArray *m2;
  char_T cv8[45];
  int32_T i;
  static const char_T cv9[45] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'm', 't', 'i', 'm', 'e', 's', '_', 'n', 'o', 'D',
    'y', 'n', 'a', 'm', 'i', 'c', 'S', 'c', 'a', 'l', 'a', 'r', 'E', 'x', 'p',
    'a', 'n', 's', 'i', 'o', 'n' };

  const mxArray *b_y;
  static const int32_T iv7[2] = { 1, 21 };

  char_T cv10[21];
  static const char_T cv11[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'i', 'n', 'n', 'e', 'r', 'd', 'i', 'm' };

  if (!(a->size[1] == b->size[0])) {
    if ((a->size[1] == 1) || ((b->size[0] == 1) && (b->size[1] == 1))) {
      y = NULL;
      m2 = mxCreateCharArray(2, iv6);
      for (i = 0; i < 45; i++) {
        cv8[i] = cv9[i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 45, m2, cv8);
      emlrtAssign(&y, m2);
      error(message(y, &emlrtMCI), &b_emlrtMCI);
    } else {
      b_y = NULL;
      m2 = mxCreateCharArray(2, iv7);
      for (i = 0; i < 21; i++) {
        cv10[i] = cv11[i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 21, m2, cv10);
      emlrtAssign(&b_y, m2);
      error(message(b_y, &c_emlrtMCI), &d_emlrtMCI);
    }
  }
}

/*
 *
 */
void dynamic_size_checks(const emxArray_real_T *a, const emxArray_real_T *b)
{
  const mxArray *y;
  static const int32_T iv4[2] = { 1, 45 };

  const mxArray *m1;
  char_T cv4[45];
  int32_T i;
  static const char_T cv5[45] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'm', 't', 'i', 'm', 'e', 's', '_', 'n', 'o', 'D',
    'y', 'n', 'a', 'm', 'i', 'c', 'S', 'c', 'a', 'l', 'a', 'r', 'E', 'x', 'p',
    'a', 'n', 's', 'i', 'o', 'n' };

  const mxArray *b_y;
  static const int32_T iv5[2] = { 1, 21 };

  char_T cv6[21];
  static const char_T cv7[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'i', 'n', 'n', 'e', 'r', 'd', 'i', 'm' };

  if (!(a->size[1] == b->size[0])) {
    if (((a->size[0] == 1) && (a->size[1] == 1)) || (b->size[0] == 1)) {
      y = NULL;
      m1 = mxCreateCharArray(2, iv4);
      for (i = 0; i < 45; i++) {
        cv4[i] = cv5[i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 45, m1, cv4);
      emlrtAssign(&y, m1);
      error(message(y, &emlrtMCI), &b_emlrtMCI);
    } else {
      b_y = NULL;
      m1 = mxCreateCharArray(2, iv5);
      for (i = 0; i < 21; i++) {
        cv6[i] = cv7[i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 21, m1, cv6);
      emlrtAssign(&b_y, m1);
      error(message(b_y, &c_emlrtMCI), &d_emlrtMCI);
    }
  }
}

/* End of code generation (mtimes.c) */
