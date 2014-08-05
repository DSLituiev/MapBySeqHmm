/*
 * cumMatrSafeF_api.c
 *
 * Code generation for function 'cumMatrSafeF_api'
 *
 * C source code generated on: Wed Jun 11 19:52:42 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cumMatrSafeF.h"
#include "cumMatrSafeF_api.h"
#include "cumMatrSafeF_emxutil.h"

/* Function Declarations */
static uint32_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId);
static int32_T c_emlrt_marshallIn(const mxArray *M, const char_T *identifier);
static int32_T d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId);
static void e_emlrt_marshallIn(const mxArray *E, const char_T *identifier,
  emxArray_real_T *y);
static uint32_T emlrt_marshallIn(const mxArray *Np, const char_T *identifier);
static const mxArray *emlrt_marshallOut(emxArray_real_T *u);
static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static void g_emlrt_marshallIn(const mxArray *A, const char_T *identifier,
  emxArray_real_T *y);
static void h_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static uint32_T i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *
  msgId);
static int32_T j_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId);
static void k_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);
static void l_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);

/* Function Definitions */
static uint32_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId)
{
  uint32_T y;
  y = i_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static int32_T c_emlrt_marshallIn(const mxArray *M, const char_T *identifier)
{
  int32_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = d_emlrt_marshallIn(emlrtAlias(M), &thisId);
  emlrtDestroyArray(&M);
  return y;
}

static int32_T d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId)
{
  int32_T y;
  y = j_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void e_emlrt_marshallIn(const mxArray *E, const char_T *identifier,
  emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  f_emlrt_marshallIn(emlrtAlias(E), &thisId, y);
  emlrtDestroyArray(&E);
}

static uint32_T emlrt_marshallIn(const mxArray *Np, const char_T *identifier)
{
  uint32_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = b_emlrt_marshallIn(emlrtAlias(Np), &thisId);
  emlrtDestroyArray(&Np);
  return y;
}

static const mxArray *emlrt_marshallOut(emxArray_real_T *u)
{
  const mxArray *y;
  static const int32_T iv8[2] = { 0, 0 };

  const mxArray *m3;
  y = NULL;
  m3 = mxCreateNumericArray(2, (int32_T *)&iv8, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m3, (void *)u->data);
  mxSetDimensions((mxArray *)m3, u->size, 2);
  emlrtAssign(&y, m3);
  return y;
}

static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  k_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const mxArray *A, const char_T *identifier,
  emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  h_emlrt_marshallIn(emlrtAlias(A), &thisId, y);
  emlrtDestroyArray(&A);
}

static void h_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  l_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static uint32_T i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *
  msgId)
{
  uint32_T ret;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "uint32", FALSE, 0U, 0);
  ret = *(uint32_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static int32_T j_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId)
{
  int32_T ret;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "int32", FALSE, 0U, 0);
  ret = *(int32_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void k_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv9[2];
  boolean_T bv0[2];
  int32_T i3;
  int32_T iv10[2];
  for (i3 = 0; i3 < 2; i3++) {
    iv9[i3] = 10000 + -9500 * i3;
    bv0[i3] = TRUE;
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv9, bv0, iv10);
  ret->size[0] = iv10[0];
  ret->size[1] = iv10[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static void l_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv11[3];
  boolean_T bv1[3];
  int32_T i;
  int32_T iv12[3];
  for (i = 0; i < 3; i++) {
    iv11[i] = -1;
    bv1[i] = TRUE;
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 3U,
    iv11, bv1, iv12);
  ret->size[0] = iv12[0];
  ret->size[1] = iv12[1];
  ret->size[2] = iv12[2];
  ret->allocatedSize = ret->size[0] * ret->size[1] * ret->size[2];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

void cumMatrSafeF_api(const mxArray * const prhs[4], const mxArray *plhs[2])
{
  emxArray_real_T *E;
  emxArray_real_T *A;
  emxArray_real_T *logAlpha;
  emxArray_real_T *logBeta;
  uint32_T Np;
  int32_T M;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  b_emxInit_real_T(&E, 2, TRUE);
  c_emxInit_real_T(&A, 3, TRUE);
  b_emxInit_real_T(&logAlpha, 2, TRUE);
  b_emxInit_real_T(&logBeta, 2, TRUE);

  /* Marshall function inputs */
  Np = emlrt_marshallIn(emlrtAliasP(prhs[0]), "Np");
  M = c_emlrt_marshallIn(emlrtAliasP(prhs[1]), "M");
  e_emlrt_marshallIn(emlrtAlias(prhs[2]), "E", E);
  g_emlrt_marshallIn(emlrtAlias(prhs[3]), "A", A);

  /* Invoke the target function */
  cumMatrSafeF(Np, M, E, A, logAlpha, logBeta);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(logAlpha);
  plhs[1] = emlrt_marshallOut(logBeta);
  logBeta->canFreeData = FALSE;
  emxFree_real_T(&logBeta);
  logAlpha->canFreeData = FALSE;
  emxFree_real_T(&logAlpha);
  A->canFreeData = FALSE;
  emxFree_real_T(&A);
  E->canFreeData = FALSE;
  emxFree_real_T(&E);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (cumMatrSafeF_api.c) */
