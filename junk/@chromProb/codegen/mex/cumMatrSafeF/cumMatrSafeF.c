/*
 * cumMatrSafeF.c
 *
 * Code generation for function 'cumMatrSafeF'
 *
 * C source code generated on: Wed Jun 11 19:52:42 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cumMatrSafeF.h"
#include "cumMatrSafeF_emxutil.h"
#include "mtimes.h"
#include "log10.h"
#include "cumMatrSafeF_mexutil.h"
#include "cumMatrSafeF_data.h"

/* Variable Definitions */
static emlrtMCInfo e_emlrtMCI = { 41, 9, "eml_min_or_max",
  "/usr/local/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtMCInfo f_emlrtMCI = { 38, 19, "eml_min_or_max",
  "/usr/local/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtMCInfo g_emlrtMCI = { 74, 9, "eml_min_or_max",
  "/usr/local/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtMCInfo h_emlrtMCI = { 73, 19, "eml_min_or_max",
  "/usr/local/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

/* Function Declarations */
static void b_eml_xgemm(int32_T n, int32_T k, const emxArray_real_T *A, const
  emxArray_real_T *B, int32_T ldb, emxArray_real_T *C);
static void eml_xgemm(int32_T m, int32_T k, const emxArray_real_T *A, int32_T
                      lda, const emxArray_real_T *B, int32_T ldb,
                      emxArray_real_T *C, int32_T ldc);

/* Function Definitions */

/*
 *
 */
static void b_eml_xgemm(int32_T n, int32_T k, const emxArray_real_T *A, const
  emxArray_real_T *B, int32_T ldb, emxArray_real_T *C)
{
  real_T alpha1;
  real_T beta1;
  char_T TRANSB;
  char_T TRANSA;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  double * alpha1_t;
  double * Aia0_t;
  double * Bib0_t;
  double * beta1_t;
  double * Cic0_t;
  if ((n < 1) || (k < 1)) {
  } else {
    alpha1 = 1.0;
    beta1 = 0.0;
    TRANSB = 'N';
    TRANSA = 'N';
    m_t = (ptrdiff_t)(1);
    n_t = (ptrdiff_t)(n);
    k_t = (ptrdiff_t)(k);
    lda_t = (ptrdiff_t)(1);
    ldb_t = (ptrdiff_t)(ldb);
    ldc_t = (ptrdiff_t)(1);
    alpha1_t = (double *)(&alpha1);
    Aia0_t = (double *)(&A->data[0]);
    Bib0_t = (double *)(&B->data[0]);
    beta1_t = (double *)(&beta1);
    Cic0_t = (double *)(&C->data[0]);
    dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, alpha1_t, Aia0_t, &lda_t, Bib0_t,
          &ldb_t, beta1_t, Cic0_t, &ldc_t);
  }
}

/*
 *
 */
static void eml_xgemm(int32_T m, int32_T k, const emxArray_real_T *A, int32_T
                      lda, const emxArray_real_T *B, int32_T ldb,
                      emxArray_real_T *C, int32_T ldc)
{
  real_T alpha1;
  real_T beta1;
  char_T TRANSB;
  char_T TRANSA;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  double * alpha1_t;
  double * Aia0_t;
  double * Bib0_t;
  double * beta1_t;
  double * Cic0_t;
  if ((m < 1) || (k < 1)) {
  } else {
    alpha1 = 1.0;
    beta1 = 0.0;
    TRANSB = 'N';
    TRANSA = 'N';
    m_t = (ptrdiff_t)(m);
    n_t = (ptrdiff_t)(1);
    k_t = (ptrdiff_t)(k);
    lda_t = (ptrdiff_t)(lda);
    ldb_t = (ptrdiff_t)(ldb);
    ldc_t = (ptrdiff_t)(ldc);
    alpha1_t = (double *)(&alpha1);
    Aia0_t = (double *)(&A->data[0]);
    Bib0_t = (double *)(&B->data[0]);
    beta1_t = (double *)(&beta1);
    Cic0_t = (double *)(&C->data[0]);
    dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, alpha1_t, Aia0_t, &lda_t, Bib0_t,
          &ldb_t, beta1_t, Cic0_t, &ldc_t);
  }
}

/*
 * function [logAlpha, logBeta] = cumMatrSafeF(Np, M, E, A)
 */
void cumMatrSafeF(uint32_T Np, int32_T M, const emxArray_real_T *E, const
                  emxArray_real_T *A, emxArray_real_T *logAlpha, emxArray_real_T
                  *logBeta)
{
  emxArray_real_T *Ac;
  emxArray_real_T *y;
  int32_T i0;
  int32_T loop_ub;
  int32_T i;
  emxArray_real_T *b_Ac;
  real_T tmp_data[500];
  emxArray_real_T *scaleA;
  real_T scalePrev;
  int32_T m;
  emxArray_int32_T *r0;
  emxArray_real_T *b_y;
  emxArray_real_T *a;
  emxArray_real_T *c_y;
  int32_T i1;
  int32_T i2;
  boolean_T guard2 = FALSE;
  boolean_T b0;
  const mxArray *d_y;
  static const int32_T iv0[2] = { 1, 36 };

  const mxArray *m0;
  char_T cv0[36];
  static const char_T cv1[36] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'a', 'u', 't', 'o', 'D', 'i', 'm', 'I', 'n', 'c',
    'o', 'm', 'p', 'a', 't', 'i', 'b', 'i', 'l', 'i', 't', 'y' };

  const mxArray *e_y;
  static const int32_T iv1[2] = { 1, 39 };

  char_T cv2[39];
  static const char_T cv3[39] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'e', 'm', 'l', '_', 'm', 'i', 'n', '_', 'o', 'r',
    '_', 'm', 'a', 'x', '_', 'v', 'a', 'r', 'D', 'i', 'm', 'Z', 'e', 'r', 'o' };

  real_T mtmp;
  int32_T exitg4;
  int32_T exitg3;
  emxArray_real_T *Bc;
  emxArray_real_T *b;
  int16_T unnamed_idx_1;
  boolean_T guard1 = FALSE;
  const mxArray *f_y;
  static const int32_T iv2[2] = { 1, 36 };

  const mxArray *g_y;
  static const int32_T iv3[2] = { 1, 39 };

  int32_T exitg2;
  int32_T exitg1;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  emxInit_real_T(&Ac, 1, TRUE);
  b_emxInit_real_T(&y, 2, TRUE);

  /*  calculate the cumulative matrices in a safe way */
  /*  assert(isa(Np,'double')); */
  /* 'cumMatrSafeF:6' assert(isscalar(Np)); */
  /*  assert(isa(M,'double')); */
  /* 'cumMatrSafeF:8' assert(isscalar(M)); */
  /* 'cumMatrSafeF:9' assert(isa(E,'double')); */
  /* 'cumMatrSafeF:10' assert(isa(A,'double')); */
  /* 'cumMatrSafeF:12' if ~isempty(Np) */
  /*      if ~isempty(T) */
  /*          A =  bsxfun(@times, ... */
  /*              permute( E(1:end-1, :), [2, 3, 1] ), T ) ; */
  /*          clear T */
  /*      else */
  /*          warning('crossMatr:emptyT', 'define transition matrix T first!'); */
  /*      end */
  /* 'cumMatrSafeF:20' logAlpha = zeros(M, Np); */
  /* 'cumMatrSafeF:21' logBeta  = zeros(M, Np); */
  /* 'cumMatrSafeF:23' logAlpha = -inf(M, Np); */
  i0 = logAlpha->size[0] * logAlpha->size[1];
  logAlpha->size[0] = M;
  logAlpha->size[1] = (int32_T)Np;
  emxEnsureCapacity((emxArray__common *)logAlpha, i0, (int32_T)sizeof(real_T));
  loop_ub = M * (int32_T)Np;
  for (i0 = 0; i0 < loop_ub; i0++) {
    logAlpha->data[i0] = rtMinusInf;
  }

  /* 'cumMatrSafeF:24' logBeta  = -inf(M, Np); */
  i0 = logBeta->size[0] * logBeta->size[1];
  logBeta->size[0] = M;
  logBeta->size[1] = (int32_T)Np;
  emxEnsureCapacity((emxArray__common *)logBeta, i0, (int32_T)sizeof(real_T));
  loop_ub = M * (int32_T)Np;
  for (i0 = 0; i0 < loop_ub; i0++) {
    logBeta->data[i0] = rtMinusInf;
  }

  /*     %% forward */
  /* 'cumMatrSafeF:27' Ac = E( end, :)'; */
  loop_ub = E->size[1];
  i = E->size[0];
  i0 = Ac->size[0];
  Ac->size[0] = loop_ub;
  emxEnsureCapacity((emxArray__common *)Ac, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < loop_ub; i0++) {
    Ac->data[i0] = E->data[(i + E->size[0] * i0) - 1];
  }

  b_emxInit_real_T(&b_Ac, 2, TRUE);

  /* 'cumMatrSafeF:28' logAlpha(M, :) = log10(Ac'); */
  i0 = b_Ac->size[0] * b_Ac->size[1];
  b_Ac->size[0] = 1;
  b_Ac->size[1] = Ac->size[0];
  emxEnsureCapacity((emxArray__common *)b_Ac, i0, (int32_T)sizeof(real_T));
  loop_ub = Ac->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_Ac->data[b_Ac->size[0] * i0] = Ac->data[i0];
  }

  i = Ac->size[0];
  i0 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = i;
  emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < i; i0++) {
    y->data[y->size[0] * i0] = b_Ac->data[i0];
  }

  emxFree_real_T(&b_Ac);
  b_log10(y);
  i = y->size[1];
  loop_ub = y->size[0] * y->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_data[i0] = y->data[i0];
  }

  for (i0 = 0; i0 < i; i0++) {
    logAlpha->data[(M + logAlpha->size[0] * i0) - 1] = tmp_data[i0];
  }

  emxInit_real_T(&scaleA, 1, TRUE);

  /* 'cumMatrSafeF:30' scalePrev = 0; */
  scalePrev = 0.0;

  /* 'cumMatrSafeF:31' scaleA = zeros(M, 1); */
  i0 = scaleA->size[0];
  scaleA->size[0] = M;
  emxEnsureCapacity((emxArray__common *)scaleA, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < M; i0++) {
    scaleA->data[i0] = 0.0;
  }

  /*  for-loop frwd */
  /* 'cumMatrSafeF:33' for m = (M-1):-1: 1 */
  m = M - 2;
  emxInit_int32_T(&r0, 1, TRUE);
  emxInit_real_T(&b_y, 1, TRUE);
  b_emxInit_real_T(&a, 2, TRUE);
  emxInit_real_T(&c_y, 1, TRUE);
  while (m + 1 > 0) {
    /* 'cumMatrSafeF:34' Ac = A(:,:, m) * (Ac * 10.^(scaleA(m+1) - scalePrev) ); */
    scalePrev = scaleA->data[m + 1] - scalePrev;
    scalePrev = muDoubleScalarPower(10.0, scalePrev);
    i0 = b_y->size[0];
    b_y->size[0] = Ac->size[0];
    emxEnsureCapacity((emxArray__common *)b_y, i0, (int32_T)sizeof(real_T));
    loop_ub = Ac->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_y->data[i0] = Ac->data[i0] * scalePrev;
    }

    loop_ub = A->size[0];
    i = A->size[1];
    i0 = a->size[0] * a->size[1];
    a->size[0] = loop_ub;
    a->size[1] = i;
    emxEnsureCapacity((emxArray__common *)a, i0, (int32_T)sizeof(real_T));
    for (i0 = 0; i0 < i; i0++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1 + a->size[0] * i0] = A->data[(i1 + A->size[0] * i0) + A->
          size[0] * A->size[1] * m];
      }
    }

    dynamic_size_checks(a, b_y);
    i0 = A->size[1];
    if ((i0 == 1) || (b_y->size[0] == 1)) {
      i0 = Ac->size[0];
      Ac->size[0] = a->size[0];
      emxEnsureCapacity((emxArray__common *)Ac, i0, (int32_T)sizeof(real_T));
      loop_ub = a->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        Ac->data[i0] = 0.0;
        i = a->size[1];
        for (i1 = 0; i1 < i; i1++) {
          Ac->data[i0] += a->data[i0 + a->size[0] * i1] * b_y->data[i1];
        }
      }
    } else {
      i0 = A->size[0];
      i1 = Ac->size[0];
      Ac->size[0] = (int16_T)i0;
      emxEnsureCapacity((emxArray__common *)Ac, i1, (int32_T)sizeof(real_T));
      loop_ub = (int16_T)i0;
      for (i0 = 0; i0 < loop_ub; i0++) {
        Ac->data[i0] = 0.0;
      }

      i0 = A->size[0];
      i1 = A->size[1];
      i = A->size[0];
      loop_ub = A->size[1];
      i2 = A->size[0];
      eml_xgemm(i0, i1, a, i, b_y, loop_ub, Ac, i2);
    }

    /* 'cumMatrSafeF:35' logAlpha(m, :)= log10(Ac) - scaleA(m+1); */
    loop_ub = logAlpha->size[1];
    i0 = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity((emxArray__common *)r0, i0, (int32_T)sizeof(int32_T));
    for (i0 = 0; i0 < loop_ub; i0++) {
      r0->data[i0] = i0;
    }

    i0 = b_y->size[0];
    b_y->size[0] = Ac->size[0];
    emxEnsureCapacity((emxArray__common *)b_y, i0, (int32_T)sizeof(real_T));
    loop_ub = Ac->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_y->data[i0] = Ac->data[i0];
    }

    for (i = 0; i < Ac->size[0]; i++) {
      b_y->data[i] = muDoubleScalarLog10(b_y->data[i]);
    }

    scalePrev = scaleA->data[m + 1];
    i0 = c_y->size[0];
    c_y->size[0] = b_y->size[0];
    emxEnsureCapacity((emxArray__common *)c_y, i0, (int32_T)sizeof(real_T));
    loop_ub = b_y->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_y->data[i0] = b_y->data[i0] - scalePrev;
    }

    i = r0->size[0];
    for (i0 = 0; i0 < i; i0++) {
      logAlpha->data[m + logAlpha->size[0] * r0->data[i0]] = c_y->data[i0];
    }

    /* 'cumMatrSafeF:36' scalePrev = scaleA(m+1); */
    scalePrev = scaleA->data[m + 1];

    /* 'cumMatrSafeF:37' scaleA(m) = - max( logAlpha( m, :) ); */
    i0 = logAlpha->size[1];
    guard2 = FALSE;
    if (i0 == 1) {
      guard2 = TRUE;
    } else {
      i0 = logAlpha->size[1];
      if (i0 != 1) {
        guard2 = TRUE;
      } else {
        b0 = FALSE;
      }
    }

    if (guard2 == TRUE) {
      b0 = TRUE;
    }

    if (b0) {
    } else {
      d_y = NULL;
      m0 = mxCreateCharArray(2, iv0);
      for (i = 0; i < 36; i++) {
        cv0[i] = cv1[i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 36, m0, cv0);
      emlrtAssign(&d_y, m0);
      error(message(d_y, &e_emlrtMCI), &f_emlrtMCI);
    }

    i0 = logAlpha->size[1];
    if (i0 > 0) {
    } else {
      e_y = NULL;
      m0 = mxCreateCharArray(2, iv1);
      for (i = 0; i < 39; i++) {
        cv2[i] = cv3[i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 39, m0, cv2);
      emlrtAssign(&e_y, m0);
      error(message(e_y, &g_emlrtMCI), &h_emlrtMCI);
    }

    i = 1;
    mtmp = logAlpha->data[m];
    i0 = logAlpha->size[1];
    if (i0 > 1) {
      if (muDoubleScalarIsNaN(mtmp)) {
        loop_ub = 2;
        do {
          exitg4 = 0;
          i0 = logAlpha->size[1];
          if (loop_ub <= i0) {
            i = loop_ub;
            if (!muDoubleScalarIsNaN(logAlpha->data[m + logAlpha->size[0] *
                 (loop_ub - 1)])) {
              mtmp = logAlpha->data[m + logAlpha->size[0] * (loop_ub - 1)];
              exitg4 = 1;
            } else {
              loop_ub++;
            }
          } else {
            exitg4 = 1;
          }
        } while (exitg4 == 0);
      }

      i0 = logAlpha->size[1];
      if (i < i0) {
        do {
          exitg3 = 0;
          i0 = logAlpha->size[1];
          if (i + 1 <= i0) {
            if (logAlpha->data[m + logAlpha->size[0] * i] > mtmp) {
              mtmp = logAlpha->data[m + logAlpha->size[0] * i];
            }

            i++;
          } else {
            exitg3 = 1;
          }
        } while (exitg3 == 0);
      }
    }

    scaleA->data[m] = -mtmp;
    m--;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, emlrtRootTLSGlobal);
  }

  emxFree_real_T(&c_y);
  emxFree_real_T(&a);
  emxFree_real_T(&b_y);
  emxFree_int32_T(&r0);
  b_emxInit_real_T(&Bc, 2, TRUE);

  /*     %% backward */
  /* 'cumMatrSafeF:41' Bc = ones(1, Np); */
  i0 = Bc->size[0] * Bc->size[1];
  Bc->size[0] = 1;
  Bc->size[1] = (int32_T)Np;
  emxEnsureCapacity((emxArray__common *)Bc, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)Np;
  for (i0 = 0; i0 < loop_ub; i0++) {
    Bc->data[i0] = 1.0;
  }

  /* 'cumMatrSafeF:42' logBeta(1, 1: Np) = 0; */
  if (1U > Np) {
    loop_ub = 0;
  } else {
    loop_ub = (int32_T)Np;
  }

  for (i0 = 0; i0 < loop_ub; i0++) {
    logBeta->data[logBeta->size[0] * i0] = 0.0;
  }

  /* 'cumMatrSafeF:44' scalePrev = 0; */
  scalePrev = 0.0;

  /* 'cumMatrSafeF:45' scaleB = zeros(M-1, 1); */
  i0 = scaleA->size[0];
  scaleA->size[0] = M - 1;
  emxEnsureCapacity((emxArray__common *)scaleA, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 <= M - 2; i0++) {
    scaleA->data[i0] = 0.0;
  }

  /* 'cumMatrSafeF:47' for m = 2:1:M */
  m = 1;
  b_emxInit_real_T(&b, 2, TRUE);
  while (m + 1 <= M) {
    /* 'cumMatrSafeF:48' Bc = Bc * A(:,:, m-1)' * (10.^(scaleB(m-1) - scalePrev)); */
    loop_ub = A->size[0];
    i = A->size[1];
    i0 = b->size[0] * b->size[1];
    b->size[0] = i;
    b->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)b, i0, (int32_T)sizeof(real_T));
    for (i0 = 0; i0 < loop_ub; i0++) {
      for (i1 = 0; i1 < i; i1++) {
        b->data[i1 + b->size[0] * i0] = A->data[(i0 + A->size[0] * i1) + A->
          size[0] * A->size[1] * (m - 1)];
      }
    }

    b_dynamic_size_checks(Bc, b);
    if ((Bc->size[1] == 1) || (b->size[0] == 1)) {
      i0 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = b->size[1];
      emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
      loop_ub = b->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        y->data[y->size[0] * i0] = 0.0;
        i = Bc->size[1];
        for (i1 = 0; i1 < i; i1++) {
          y->data[y->size[0] * i0] += Bc->data[Bc->size[0] * i1] * b->data[i1 +
            b->size[0] * i0];
        }
      }
    } else {
      unnamed_idx_1 = (int16_T)b->size[1];
      i0 = y->size[0] * y->size[1];
      y->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
      i0 = y->size[0] * y->size[1];
      y->size[1] = unnamed_idx_1;
      emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
      loop_ub = unnamed_idx_1;
      for (i0 = 0; i0 < loop_ub; i0++) {
        y->data[i0] = 0.0;
      }

      b_eml_xgemm(b->size[1], Bc->size[1], Bc, b, Bc->size[1], y);
    }

    scalePrev = scaleA->data[m - 1] - scalePrev;
    scalePrev = muDoubleScalarPower(10.0, scalePrev);
    i0 = Bc->size[0] * Bc->size[1];
    Bc->size[0] = 1;
    Bc->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)Bc, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[0] * y->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      Bc->data[i0] = y->data[i0] * scalePrev;
    }

    /* 'cumMatrSafeF:49' logBeta(m, :) = log10(Bc) - scaleB(m-1); */
    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    i = y->size[0];
    loop_ub = y->size[1];
    loop_ub *= i;
    for (i0 = 0; i0 < loop_ub; i0++) {
      y->data[i0] *= scalePrev;
    }

    b_log10(y);
    scalePrev = scaleA->data[m - 1];
    loop_ub = y->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      logBeta->data[m + logBeta->size[0] * i0] = y->data[y->size[0] * i0] -
        scalePrev;
    }

    /* 'cumMatrSafeF:50' scalePrev = scaleB(m-1); */
    scalePrev = scaleA->data[m - 1];

    /* 'cumMatrSafeF:51' scaleB(m) = - max(logBeta(m, :)); */
    i0 = logBeta->size[1];
    guard1 = FALSE;
    if (i0 == 1) {
      guard1 = TRUE;
    } else {
      i0 = logBeta->size[1];
      if (i0 != 1) {
        guard1 = TRUE;
      } else {
        b0 = FALSE;
      }
    }

    if (guard1 == TRUE) {
      b0 = TRUE;
    }

    if (b0) {
    } else {
      f_y = NULL;
      m0 = mxCreateCharArray(2, iv2);
      for (i = 0; i < 36; i++) {
        cv0[i] = cv1[i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 36, m0, cv0);
      emlrtAssign(&f_y, m0);
      error(message(f_y, &e_emlrtMCI), &f_emlrtMCI);
    }

    i0 = logBeta->size[1];
    if (i0 > 0) {
    } else {
      g_y = NULL;
      m0 = mxCreateCharArray(2, iv3);
      for (i = 0; i < 39; i++) {
        cv2[i] = cv3[i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 39, m0, cv2);
      emlrtAssign(&g_y, m0);
      error(message(g_y, &g_emlrtMCI), &h_emlrtMCI);
    }

    i = 1;
    mtmp = logBeta->data[m];
    i0 = logBeta->size[1];
    if (i0 > 1) {
      if (muDoubleScalarIsNaN(mtmp)) {
        loop_ub = 2;
        do {
          exitg2 = 0;
          i0 = logBeta->size[1];
          if (loop_ub <= i0) {
            i = loop_ub;
            if (!muDoubleScalarIsNaN(logBeta->data[m + logBeta->size[0] *
                 (loop_ub - 1)])) {
              mtmp = logBeta->data[m + logBeta->size[0] * (loop_ub - 1)];
              exitg2 = 1;
            } else {
              loop_ub++;
            }
          } else {
            exitg2 = 1;
          }
        } while (exitg2 == 0);
      }

      i0 = logBeta->size[1];
      if (i < i0) {
        do {
          exitg1 = 0;
          i0 = logBeta->size[1];
          if (i + 1 <= i0) {
            if (logBeta->data[m + logBeta->size[0] * i] > mtmp) {
              mtmp = logBeta->data[m + logBeta->size[0] * i];
            }

            i++;
          } else {
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      }
    }

    scaleA->data[m] = -mtmp;
    m++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, emlrtRootTLSGlobal);
  }

  emxFree_real_T(&b);
  emxFree_real_T(&Bc);
  emxFree_real_T(&scaleA);
  emxFree_real_T(&y);
  emxFree_real_T(&Ac);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (cumMatrSafeF.c) */
