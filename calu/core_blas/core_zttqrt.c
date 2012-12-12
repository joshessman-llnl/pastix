/**
 *
 * @file core_zttqrt.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"
#undef REAL
#define COMPLEX

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zttqrt computes a QR factorization of a rectangular matrix
 *  formed by coupling a complex N-by-N upper triangular tile A1
 *  on top of a complex M-by-N upper trapezoidal tile A2:
 *
 *    | A1 | = Q * R
 *    | A2 |
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k), where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A2(1:m,i),
 *  and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A2.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A1 and A2.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the N-by-N tile A1.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the N-by-N upper trapezoidal tile R;
 *         the elements below the diagonal are not referenced.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1.  LDA1 >= max(1,N).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N upper triangular tile A2.
 *         On exit, the elements on and above the diagonal of the array
 *         with the array TAU, represent
 *         the unitary tile Q as a product of elementary reflectors
 *         (see Further Details).
 *
 * @param[in] LDA2
 *         The leading dimension of the array A2.  LDA2 >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] TAU
 *         The scalar factors of the elementary reflectors (see Further
 *         Details).
 *
 * @param[in,out] WORK
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zttqrt = PCORE_zttqrt
#define CORE_zttqrt PCORE_zttqrt
/* Trick to call the version without tracing */
#define CORE_zlaset PCORE_zlaset
void
CORE_zlaset(PLASMA_enum uplo, int n1, int n2,
            PLASMA_Complex64_t alpha, PLASMA_Complex64_t beta,
            PLASMA_Complex64_t *tileA, int ldtilea);
#define CORE_zpemv PCORE_zpemv
int
CORE_zpemv(PLASMA_enum trans, int storev,
           int M, int N, int L,
           PLASMA_Complex64_t ALPHA,
           const PLASMA_Complex64_t *A, int LDA,
           const PLASMA_Complex64_t *X, int INCX,
           PLASMA_Complex64_t BETA,
           PLASMA_Complex64_t *Y, int INCY,
           PLASMA_Complex64_t *WORK);
#endif
int CORE_zttqrt(int M, int N, int IB,
                PLASMA_Complex64_t *A1, int LDA1,
                PLASMA_Complex64_t *A2, int LDA2,
                PLASMA_Complex64_t *T, int LDT,
                PLASMA_Complex64_t *TAU, PLASMA_Complex64_t *WORK)
{
    static int                ione  = 1;
    static PLASMA_Complex64_t zone  = 1.0;
    static PLASMA_Complex64_t zzero = 0.0;

    PLASMA_Complex64_t alpha;
    int i, j, l, ii, sb, mi, ni;

    /* Check input arguments */
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if (IB < 0) {
        coreblas_error(3, "Illegal value of IB");
        return -3;
    }
    if ((LDA2 < max(1,M)) && (M > 0)) {
        coreblas_error(7, "Illegal value of LDA2");
        return -7;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    /* TODO: Need to check why some cases require
     *  this to not have uninitialized values */
    CORE_zlaset( PlasmaUpperLower, IB, N,
                 0., 0., T, LDT);

    for (ii = 0; ii < N; ii += IB) {
        sb = min(N-ii, IB);
        for (i = 0; i < sb; i++) {
            j  = ii + i;
            mi = min( j + 1, M );
            ni = sb-i-1;

            /*
             * Generate elementary reflector H( II*IB+I ) to annihilate
             * A( II*IB+I:mi, II*IB+I ).
             */
            LAPACKE_zlarfg_work(
                    mi+1, &A1[LDA1*j+j], &A2[LDA2*j], ione, &TAU[j]);

            if (ni > 0) {
                /*
                 * Apply H( II*IB+I ) to A( II*IB+I:M, II*IB+I+1:II*IB+IB ) from the left.
                 */
                cblas_zcopy(
                    ni,
                    &A1[LDA1*(j+1)+j], LDA1,
                    WORK, 1);

#ifdef COMPLEX
                LAPACKE_zlacgv_work(ni, WORK, 1);
#endif
                cblas_zgemv(
                    CblasColMajor, (CBLAS_TRANSPOSE)PlasmaConjTrans,
                    mi, ni,
                    CBLAS_SADDR(zone), &A2[LDA2*(j+1)], LDA2,
                                       &A2[LDA2*j],     1,
                    CBLAS_SADDR(zone), WORK,            1);
#ifdef COMPLEX
                LAPACKE_zlacgv_work(ni, WORK, 1);
#endif
                alpha = -conj(TAU[j]);
                cblas_zaxpy(
                    ni, CBLAS_SADDR(alpha),
                    WORK, 1,
                    &A1[LDA1*(j+1)+j], LDA1);
#ifdef COMPLEX
                LAPACKE_zlacgv_work(ni, WORK, 1);
#endif
                cblas_zgerc(
                    CblasColMajor, mi, ni,
                    CBLAS_SADDR(alpha), &A2[LDA2*j], 1,
                    WORK, 1,
                    &A2[LDA2*(j+1)], LDA2);
            }

            /*
             * Calculate T
             *
             * T(0:i-1, j) = alpha * A2(0:M-1, ii:j-1)' * A2(0:M-1, j)
             */

            if ( i > 0 ) {

                l = min(i, max(0, M-ii));
                alpha = -(TAU[j]);

                CORE_zpemv(
                        PlasmaConjTrans, PlasmaColumnwise,
                        min(j, M), i, l,
                        alpha, &A2[LDA2*ii], LDA2,
                               &A2[LDA2*j],  1,
                        zzero, &T[LDT*j],    1,
                        WORK);

                /* T(0:i-1, j) = T(0:i-1, ii:j-1) * T(0:i-1, j) */
                cblas_ztrmv(
                        CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                        (CBLAS_TRANSPOSE)PlasmaNoTrans,
                        (CBLAS_DIAG)PlasmaNonUnit,
                        i, &T[LDT*ii], LDT,
                           &T[LDT*j], 1);

            }

            T[LDT*j+i] = TAU[j];
        }

        /* Apply Q' to the rest of the matrix to the left  */
        if (N > ii+sb) {
            mi = min(ii+sb, M);
            ni = N-(ii+sb);
            l  = min(sb, max(0, mi-ii));
            CORE_zparfb(
                PlasmaLeft, PlasmaConjTrans,
                PlasmaForward, PlasmaColumnwise,
                IB, ni, mi, ni, sb, l,             //replaced sb by IB
                &A1[LDA1*(ii+sb)+ii], LDA1,
                &A2[LDA2*(ii+sb)], LDA2,
                &A2[LDA2*ii], LDA2,
                &T[LDT*ii], LDT,
                WORK, sb);
        }
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zttqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex64_t *A1, int lda1,
                       PLASMA_Complex64_t *A2, int lda2,
                       PLASMA_Complex64_t *T, int ldt)
{
    DAG_CORE_TTQRT;
    QUARK_Insert_Task(quark, CORE_zttqrt_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A1,            INOUT|QUARK_REGION_D|QUARK_REGION_U,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A2,            INOUT|QUARK_REGION_D|QUARK_REGION_U|LOCALITY,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,    T,             OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb,       NULL,          SCRATCH,
        sizeof(PLASMA_Complex64_t)*ib*nb,    NULL,          SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zttqrt_quark = PCORE_zttqrt_quark
#define CORE_zttqrt_quark PCORE_zttqrt_quark
#endif
void CORE_zttqrt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    PLASMA_Complex64_t *A1;
    int lda1;
    PLASMA_Complex64_t *A2;
    int lda2;
    PLASMA_Complex64_t *T;
    int ldt;
    PLASMA_Complex64_t *TAU;
    PLASMA_Complex64_t *WORK;

    quark_unpack_args_11(quark, m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
    CORE_zttqrt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
}
