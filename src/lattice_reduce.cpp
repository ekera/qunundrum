/*!
 * \file    lattice_reduce.cpp
 * \ingroup lattice
 *
 * \brief   The definition of functions for reducing lattice bases.
 */

#include "lattice_reduce.h"

#include "lattice.h"
#include "lattice_gso.h"
#include "parameters.h"
#include "diagonal_parameters.h"
#include "errors.h"
#include "math.h"

#include <fplll/fplll.h>

#include <gmp.h>

#include <stdint.h>

using namespace std;
using namespace fplll;

/*!
 * \brief   A more conservative alternative to the bkz_reduce() function in
 *          fpLLL that seemingly manages to reduce bases that causes the
 *          bkz_reduce() to freeze up and never return.
 *
 * This function will be removed in the future when the reason for which fpLLL
 * freezes up has been identified and corrected. There is not very much
 * documentation available for fpLLL and it may well be that we do not properly
 * call the function in fpLLL, hence causing it to freeze up.
 *
 * \internal
 *
 * \param[in, out] b    The basis to reduce.
 * \param[in] param     The BKZ parameters for fpLLL.
 * \param[in] lll_delta The delta parameter in fpLLL.
 * \param[in] u         An empty transformation matrix.
 * \param[in] u_inv     An empty inverse transformation matrix.
 *
 * \return If reduction succeeds #RED_SUCCESS is returned, otherwise some
 *         other value is returned.
 */
static int conservative_bkz_reduction(
  ZZ_mat<mpz_t> &b,
  const BKZParam &param,
  double lll_delta,
  ZZ_mat<mpz_t> &u,
  ZZ_mat<mpz_t> &u_inv)
{
  const int gso_flags = GSO_INT_GRAM;

  if ((b.get_rows() == 0) || (b.get_cols() == 0)) {
    return RED_SUCCESS;
  }

  /* LLL-reduce the basis. */
  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>>
    m_gso(b, u, u_inv, gso_flags);

  LLLReduction<Z_NR<mpz_t>, FP_NR<mpfr_t>>
    lll_obj(m_gso, lll_delta, LLL_DEF_ETA, LLL_DEFAULT);

  lll_obj.lll();

  if (RED_SUCCESS != lll_obj.status) {
    return lll_obj.status;
  }

  /* BKZ-reduce the LLL-reduced basis. */
  LLLReduction<Z_NR<mpz_t>, FP_NR<mpfr_t>>
    lll_obj2(m_gso, lll_delta, LLL_DEF_ETA, LLL_DEFAULT);

  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>>
    m_gso2(b, u, u_inv, gso_flags);

  BKZReduction<Z_NR<mpz_t>, FP_NR<mpfr_t>>
    bkz_obj(m_gso2, lll_obj2, param);

  bkz_obj.bkz();

  return bkz_obj.status;
}

void lattice_compute_reduced_basis(
  ZZ_mat<mpz_t> &A,
  const mpz_t * const js,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t red_precision)
{
  /* Setup constants. */
  mpz_t pow2lm;
  mpz_init(pow2lm);
  mpz_set_ui(pow2lm, 0);
  mpz_setbit(pow2lm, parameters->l + parameters->m);

  /* Setup the A matrix. */
  A.resize(n + 1, n + 1);

  for (uint32_t j = 0; j < n; j++) {
    mpz_set(A[0][j].get_data(), js[j]);
  }

  mpz_set_ui(A[0][n].get_data(), 1);

  for (uint32_t i = 1; i <= n; i++) {
    for (uint32_t j = 0; j <= n; j++) {
      mpz_set_ui(A[i][j].get_data(), 0);
    }

    mpz_set(A[i][i - 1].get_data(), pow2lm);
  }

  /* Reduce the matrix. */
  int status;

  if (REDUCTION_ALGORITHM_LLL == algorithm) {
    /* Reduce the matrix using Lenstra-Lenstra-Lovász reduction. */
    status = lll_reduction(
                A,
                LLL_DEF_DELTA, /* delta */
                LLL_DEF_ETA, /* eta */
                LM_WRAPPER, /* method */
                FT_DEFAULT, /* floating point */
                0, /* precision */
                LLL_DEFAULT); /* flags */
    if (RED_SUCCESS != status) {
      critical("lattice_compute_reduced_basis(): "
        "Failed to reduce matrix A using LLL.");
    }
  } else if (REDUCTION_ALGORITHM_BKZ == algorithm) {
    /* Reduce the matrix using block Korkin-Zolotarev reduction. */

    int old_precision = FP_NR<mpfr_t>::set_prec(red_precision);

    #ifndef BKZ_FPLLL
    ZZ_mat<mpz_t> empty_mat;

    ZZ_mat<mpz_t> &u = empty_mat;
    ZZ_mat<mpz_t> &u_inv = empty_mat;

    /* Setup a simple BKZ reduction strategy. */
    vector<Strategy> strategies;
    BKZParam param(((n + 1) > LATTICE_MAX_BLOCK_SIZE_BKZ) ?
      LATTICE_MAX_BLOCK_SIZE_BKZ : n + 1, strategies);

    status = conservative_bkz_reduction(
              A,
              param,
              LLL_DEF_DELTA,
              u,
              u_inv);

    empty_mat.clear();
    #else
    status = bkz_reduction(
                A,
                ((n + 1) > LATTICE_MAX_BLOCK_SIZE_BKZ) ?
                  LATTICE_MAX_BLOCK_SIZE_BKZ : n + 1, /* block size */
                BKZ_DEFAULT, /* flags */
                FT_MPFR, /* float type */
                red_precision); /* float precision */
    #endif

    FP_NR<mpfr_t>::set_prec(old_precision);

    if (RED_SUCCESS != status) {
      critical("lattice_compute_reduced_basis(): "
        "Failed to reduce matrix A using BKZ.");
    }
  } else if (REDUCTION_ALGORITHM_HKZ == algorithm) {
    /* Reduce the matrix using Hermite Korkin-Zolotarev reduction. */
    status = hkz_reduction(A, HKZ_DEFAULT, FT_MPFR, red_precision);
    if (RED_SUCCESS != status) {
      critical("lattice_compute_reduced_basis(): "
        "Failed to reduce matrix A using HKZ.");
    }
  } else {
    critical("lattice_compute_reduced_basis(): Unknown reduction algorithm.");
  }

  /* Clear memory. */
  mpz_clear(pow2lm);
}

void lattice_compute_reduced_diagonal_basis(
  ZZ_mat<mpz_t> &A,
  const mpz_t * const js,
  const uint32_t n,
  const Diagonal_Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t red_precision)
{
  /* Setup variables. */
  mpz_t alpha;
  mpz_init(alpha);

  /* Setup constants. */
  mpz_t pow2msigma;
  mpz_init(pow2msigma);
  mpz_set_ui(pow2msigma, 0);
  mpz_setbit(pow2msigma, parameters->m + parameters->sigma);

  mpz_t pow2lsigma;
  mpz_init(pow2lsigma);
  mpz_set_ui(pow2lsigma, 0);
  mpz_setbit(pow2lsigma, parameters->l - parameters->sigma);

  /* Setup the A matrix. */
  A.resize(n + 1, n + 1);

  for (uint32_t j = 0; j < n; j++) {
    mpz_mul(alpha, parameters->r, js[j]); /* alpha = rj */

    mpz_set(A[0][j].get_data(), alpha); /* alpha = rj */

    mod_reduce(alpha, pow2msigma);

    mpz_sub(A[0][j].get_data(), A[0][j].get_data(), alpha);

    mpz_mul(A[0][j].get_data(), A[0][j].get_data(), pow2lsigma);
  }

  mpz_set(A[0][n].get_data(), parameters->r);

  for (uint32_t i = 1; i <= n; i++) {
    for (uint32_t j = 0; j <= n; j++) {
      mpz_set_ui(A[i][j].get_data(), 0);
    }

    mpz_set(A[i][i - 1].get_data(), pow2msigma);
    mpz_mul(A[i][i - 1].get_data(), A[i][i - 1].get_data(), parameters->r);
    mpz_mul(A[i][i - 1].get_data(), A[i][i - 1].get_data(), pow2lsigma);
  }

  /* Reduce the matrix. */
  int status;

  if (REDUCTION_ALGORITHM_LLL == algorithm) {
    /* Reduce the matrix using Lenstra-Lenstra-Lovász reduction. */
    status = lll_reduction(
                A,
                LLL_DEF_DELTA, /* delta */
                LLL_DEF_ETA, /* eta */
                LM_WRAPPER, /* method */
                FT_DEFAULT, /* floating point */
                0, /* precision */
                LLL_DEFAULT); /* flags */
    if (RED_SUCCESS != status) {
      critical("lattice_compute_reduced_diagonal_basis(): "
        "Failed to reduce matrix A using LLL.");
    }
  } else if (REDUCTION_ALGORITHM_BKZ == algorithm) {
    /* Reduce the matrix using block Korkin-Zolotarev reduction. */

    int old_precision = FP_NR<mpfr_t>::set_prec(red_precision);

    #ifndef BKZ_FPLLL
    ZZ_mat<mpz_t> empty_mat;

    ZZ_mat<mpz_t> &u = empty_mat;
    ZZ_mat<mpz_t> &u_inv = empty_mat;

    /* Setup a simple BKZ reduction strategy. */
    vector<Strategy> strategies;
    BKZParam param(((n + 1) > LATTICE_MAX_BLOCK_SIZE_BKZ) ?
      LATTICE_MAX_BLOCK_SIZE_BKZ : n + 1, strategies);

    status = conservative_bkz_reduction(
              A,
              param,
              LLL_DEF_DELTA,
              u,
              u_inv);

    empty_mat.clear();
    #else
    status = bkz_reduction(
                A,
                ((n + 1) > LATTICE_MAX_BLOCK_SIZE_BKZ) ?
                  LATTICE_MAX_BLOCK_SIZE_BKZ : n + 1, /* block size */
                BKZ_DEFAULT, /* flags */
                FT_MPFR, /* float type */
                red_precision); /* float precision */
    #endif

    FP_NR<mpfr_t>::set_prec(old_precision);

    if (RED_SUCCESS != status) {
      critical("lattice_compute_reduced_diagonal_basis(): "
        "Failed to reduce matrix A using BKZ.");
    }
  } else if (REDUCTION_ALGORITHM_HKZ == algorithm) {
    /* Reduce the matrix using Hermite Korkin-Zolotarev reduction. */
    status = hkz_reduction(A, HKZ_DEFAULT, FT_MPFR, red_precision);
    if (RED_SUCCESS != status) {
      critical("lattice_compute_reduced_diagonal_basis(): "
        "Failed to reduce matrix A using HKZ.");
    }
  } else {
    critical("lattice_compute_reduced_diagonal_basis(): "
      "Unknown reduction algorithm.");
  }

  /* Clear memory. */
  mpz_clear(alpha);
  mpz_clear(pow2msigma);
  mpz_clear(pow2lsigma);
}
