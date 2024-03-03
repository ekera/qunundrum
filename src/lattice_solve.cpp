/*!
 * \file    lattice_solve.cpp
 * \ingroup lattice_solve
 *
 * \brief   The definition of functions for recovering d and r by using nearest
 *          plane solvers, and similar techniques, that do not require
 *          enumerating vectors in the lattice L.
 */

#include "lattice_solve.h"

#include "common.h"
#include "diagonal_parameters.h"
#include "errors.h"
#include "lattice.h"
#include "lattice_algebra.h"
#include "lattice_babai.h"
#include "lattice_gso.h"
#include "lattice_reduce.h"
#include "lattice_smoothness.h"
#include "math.h"
#include "parameters.h"
#include "sample.h"
#include "timer.h"

#include <gmp.h>
#include <mpfr.h>

#include <fplll/fplll.h>

#include <stdint.h>

#include <vector>

using namespace std;
using namespace fplll;

/*!
 * \brief   A bound on the largest multiple of the last component of the
 *          shortest non-zero vector in the reduced basis to be added to or
 *          subtracted from the candidate for d.
 *
 * Note that the last component of the shortest non-zero vector in the reduced
 * lattice basis is expected to be a divisor of r when solving for an order or
 * for a general discrete logarithm.
 */
#define SEARCH_BOUND_SHORTEST_VECTOR_MULTIPLE 256

void lattice_solve_reduced_basis_for_d(
  Lattice_Status_Recovery * const status_d,
  const ZZ_mat<mpz_t> &A,
  const FP_mat<mpfr_t> &G,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  const uint32_t precision,
  const bool detect_smooth_r)
{
  /* Check dimensions. */
  if ((A.get_rows() != (int)(n + 1)) || (A.get_cols() != (int)(n + 1))) {
    critical("lattice_solve_reduced_basis_for_d(): "
      "Incorrect dimensions for the matrix A.");
  }

  if ((G.get_rows() != (int)(n + 1)) || (G.get_cols() != (int)(n + 1))) {
    critical("lattice_solve_reduced_basis_for_d(): "
      "Incorrect dimensions for the matrix G.");
  }

  /* Setup constants. */
  mpz_t pow2m;
  mpz_init(pow2m);
  mpz_set_ui(pow2m, 0);
  mpz_setbit(pow2m, parameters->m);

  mpz_t pow2lm;
  mpz_init(pow2lm);
  mpz_set_ui(pow2lm, 0);
  mpz_setbit(pow2lm, parameters->l + parameters->m);

  mpz_t pow2lm_half;
  mpz_init(pow2lm_half);
  mpz_set_ui(pow2lm_half, 0);
  mpz_setbit(pow2lm_half, parameters->l + parameters->m - 1);

  /* Setup variables. */
  mpz_t candidate_d;
  mpz_init(candidate_d);

  /* Setup the target vector. */
  vector<Z_NR<mpz_t>> target;
  target.resize(n + 1);

  for (uint32_t j = 0; j < n; j++) {
    mpz_mul(target[j].get_data(), ks[j], pow2m);
    mpz_neg(target[j].get_data(), target[j].get_data());
    mpz_mod(target[j].get_data(), target[j].get_data(), pow2lm);

    if (mpz_cmp(target[j].get_data(), pow2lm_half) >= 0) {
      mpz_sub(target[j].get_data(), target[j].get_data(), pow2lm);
    }
  }

  mpz_set_ui(target[n].get_data(), 0);

  /* Solve for the closest vector. */
  vector<Z_NR<mpz_t>> solution;

  /* Compute the closest vector. */
  babai_closest_vector(solution, target, G, A, n, precision);

  /* Extract the candidate d. */
  mpz_abs(candidate_d, solution[n].get_data());

  if (mpz_cmp(candidate_d, parameters->d) == 0) {
    (*status_d) = LATTICE_STATUS_RECOVERED_IMMEDIATE;
  } else {
    /* Setup variables. */
    mpz_t multiple;
    mpz_init(multiple);

    mpz_t tmp_d;
    mpz_init(tmp_d);

    (*status_d) = LATTICE_STATUS_NOT_RECOVERED;

    mpz_set(candidate_d, solution[n].get_data());

    for (uint32_t i = 1; i <= SEARCH_BOUND_SHORTEST_VECTOR_MULTIPLE; i++) {
      mpz_mul_ui(multiple, A[0][n].get_data(), i);

      mpz_add(tmp_d, candidate_d, multiple);
      mpz_abs(tmp_d, tmp_d);
      if (mpz_cmp(tmp_d, parameters->d) == 0) {
        (*status_d) = LATTICE_STATUS_RECOVERED_SEARCH;
        break;
      }

      mpz_sub(tmp_d, candidate_d, multiple);
      mpz_abs(tmp_d, tmp_d);
      if (mpz_cmp(tmp_d, parameters->d) == 0) {
        (*status_d) = LATTICE_STATUS_RECOVERED_SEARCH;
        break;
      }
    }

    /* Clear memory. */
    mpz_clear(multiple);
    mpz_clear(tmp_d);
  }

  if (detect_smooth_r && (LATTICE_STATUS_NOT_RECOVERED == (*status_d))) {
    /* Setup variables. */
    mpz_t candidate_r;
    mpz_init(candidate_r);

    mpz_t test_r;
    mpz_init(test_r);

    /* Extract the candidate r. */
    mpz_abs(candidate_r, A[0][n].get_data());

    /* Test if r is a multiple of the candidate r. */
    mpz_mod(test_r, parameters->r, candidate_r);

    if (mpz_cmp_ui(test_r, 0) == 0) { /* candidate_r divides r */
      /* Let test_r be the multiple of the candidate r required to form r. */
      mpz_div(test_r, parameters->r, candidate_r);
        /* test_r = r / candidate_r */

      /* If the multiple is greater than one but smooth, check if we have
       * recovered d modulo the non-smooth part of r. */
      if (mpz_cmp_ui(test_r, 1) > 0) {
        if (lattice_smoothness_is_smooth(
              test_r, LATTICE_SMOOTHNESS_CONSTANT_C, parameters->m))
        {
          mpz_t reduced_d;
          mpz_init(reduced_d);

          mpz_t reduced_candidate_d;
          mpz_init(reduced_candidate_d);

          mpz_mod(reduced_candidate_d, candidate_d, candidate_r);
          mpz_mod(reduced_d, parameters->d, candidate_r);

          if (mpz_cmp(reduced_d, reduced_candidate_d) == 0) {
            (*status_d) = LATTICE_STATUS_RECOVERED_IMMEDIATE_SMOOTH;
          }

          /* Clear memory. */
          mpz_clear(reduced_d);
          mpz_clear(reduced_candidate_d);
        }
      }
    }

    /* Clear memory. */
    mpz_clear(candidate_r);
    mpz_clear(test_r);
  }

  /* Clear memory. */
  target.clear();
  solution.clear();

  mpz_clear(pow2m);
  mpz_clear(pow2lm);
  mpz_clear(pow2lm_half);
  mpz_clear(candidate_d);
}

void lattice_solve_reduced_basis_for_r(
  Lattice_Status_Recovery * const status_r,
  const ZZ_mat<mpz_t> &A,
  const uint32_t n,
  const Parameters * const parameters,
  const bool detect_smooth_r)
{
  /* Check dimensions. */
  if ((A.get_rows() != (int)(n + 1)) || (A.get_cols() != (int)(n + 1))) {
    critical("lattice_solve_reduced_basis_for_r(): "
      "Incorrect dimensions for the matrix A.");
  }

  /* Setup variables. */
  mpz_t candidate_r;
  mpz_init(candidate_r);

  mpz_t test_r;
  mpz_init(test_r);

  /* Extract the candidate r. */
  mpz_abs(candidate_r, A[0][n].get_data());

  /* Test if r is a multiple of the candidate r. */
  mpz_mod(test_r, parameters->r, candidate_r);

  if (mpz_cmp_ui(test_r, 0) != 0) {
    /* If not, we fail to recover r. */
    (*status_r) = LATTICE_STATUS_NOT_RECOVERED;
  } else {
    /* Let test_r be the multiple of the candidate r required to form r. */
    mpz_div(test_r, parameters->r, candidate_r); /* test_r = r / candidate_r */

    if (mpz_cmp_ui(test_r, 1) == 0) {
      /* If the multiple is one, we recover r immediately. */
      (*status_r) = LATTICE_STATUS_RECOVERED_IMMEDIATE;
    } else if (mpz_cmp_ui(test_r, SEARCH_BOUND_SHORTEST_VECTOR_MULTIPLE) <= 0) {
      /* If the multiple is small, we recover r by performing a small search. */
      (*status_r) = LATTICE_STATUS_RECOVERED_SEARCH;
    } else if (detect_smooth_r &&
                lattice_smoothness_is_smooth(
                  test_r, LATTICE_SMOOTHNESS_CONSTANT_C, parameters->m))
    {
      /* If the multiple is smooth, we recover by exponentiating to powers. */
      (*status_r) = LATTICE_STATUS_RECOVERED_IMMEDIATE_SMOOTH;
    } else {
      /* Otherwise, we fail to recover r. */
      (*status_r) = LATTICE_STATUS_NOT_RECOVERED;
    }
  }

  /* Clear memory. */
  mpz_clear(candidate_r);
  mpz_clear(test_r);
}

void lattice_solve_for_d(
  Lattice_Status_Recovery * const status_d,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  uint32_t precision,
  const bool detect_smooth_r)
{
  /* Setup the A matrix. */
  ZZ_mat<mpz_t> A(n + 1, n + 1);

  /* Setup the Gram-Schmidt matrices. */
  FP_mat<mpfr_t> G(n + 1, n + 1);
  FP_mat<mpfr_t> M(n + 1, n + 1);

  while (TRUE) {
    /* Reduce the A matrix. */
    lattice_compute_reduced_basis(
      A,
      js,
      n,
      parameters,
      (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) ?
        REDUCTION_ALGORITHM_LLL : algorithm,
      parameters->m); /* red_precision */

    /* Compute the Gram-Schmidt matrices. */
    gram_schmidt_orthogonalization(M, G, A, n, precision);

    /* Solve for d. */
    lattice_solve_reduced_basis_for_d(
      status_d,
      A,
      G,
      ks,
      n,
      parameters,
      precision,
      detect_smooth_r);

    if (LATTICE_STATUS_NOT_RECOVERED != (*status_d)) {
      break;
    }

    if (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) {
      algorithm = REDUCTION_ALGORITHM_BKZ;
    } else {
      break;
    }
  }

  /* Clear memory. */
  A.clear();
  G.clear();
  M.clear();
}

void lattice_solve_for_r(
  Lattice_Status_Recovery * const status_r,
  const mpz_t * const js,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const bool detect_smooth_r)
{
  /* Setup the A matrix. */
  ZZ_mat<mpz_t> A(n + 1, n + 1);

  while (TRUE) {
    /* Reduce the A matrix. */
    lattice_compute_reduced_basis(
      A,
      js,
      n,
      parameters,
      (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) ?
        REDUCTION_ALGORITHM_LLL : algorithm,
      parameters->m); /* red_precision */

    /* Solve for r. */
    lattice_solve_reduced_basis_for_r(
      status_r,
      A,
      n,
      parameters,
      detect_smooth_r);

    if (LATTICE_STATUS_NOT_RECOVERED != (*status_r)) {
      break;
    }

    if (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) {
      algorithm = REDUCTION_ALGORITHM_BKZ;
    } else {
      break;
    }
  }

  /* Clear memory. */
  A.clear();
}

void lattice_solve_for_d_r(
  Lattice_Status_Recovery * const status_d,
  Lattice_Status_Recovery * const status_r,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const bool detect_smooth_r)
{
  /* Setup the A matrix. */
  ZZ_mat<mpz_t> A(n + 1, n + 1);

  /* Setup the Gram-Schmidt matrices. */
  FP_mat<mpfr_t> G(n + 1, n + 1);
  FP_mat<mpfr_t> M(n + 1, n + 1);

  /* Initially set the statuses to not recovered. */
  (*status_d) = LATTICE_STATUS_NOT_RECOVERED;
  (*status_r) = LATTICE_STATUS_NOT_RECOVERED;

  while (TRUE) {
    /* Reduce the A matrix. */
    lattice_compute_reduced_basis(
      A,
      js,
      n,
      parameters,
      (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) ?
        REDUCTION_ALGORITHM_LLL : algorithm,
      parameters->m); /* red_precision */

    if (LATTICE_STATUS_NOT_RECOVERED == (*status_r)) {
      /* Solve for r. */
      lattice_solve_reduced_basis_for_r(
        status_r,
        A,
        n,
        parameters,
        detect_smooth_r);
    }

    if (LATTICE_STATUS_NOT_RECOVERED == (*status_d)) {
      /* Compute the Gram-Schmidt matrices. */
      gram_schmidt_orthogonalization(M, G, A, n, precision);

      /* Solve for d. */
      lattice_solve_reduced_basis_for_d(
        status_d,
        A,
        G,
        ks,
        n,
        parameters,
        precision,
        detect_smooth_r);
    }

    /* Check if both d and r have been recovered. Note that timeouts cannot
     * occur when solving without enumerating, so the below check suffices. */
    if ((LATTICE_STATUS_NOT_RECOVERED != (*status_d)) &&
        (LATTICE_STATUS_NOT_RECOVERED != (*status_r)))
    {
      break;
    }

    if (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) {
      algorithm = REDUCTION_ALGORITHM_BKZ;
    } else {
      break;
    }
  }

  /* Clear memory. */
  A.clear();
  G.clear();
  M.clear();
}

void lattice_solve_reduced_basis_for_d_given_r(
  Lattice_Status_Recovery * const status_d,
  const ZZ_mat<mpz_t> &A,
  const FP_mat<mpfr_t> &G,
  const mpz_t * const ks,
  const uint32_t n,
  const Diagonal_Parameters * const parameters,
  const uint32_t precision)
{
  /* Check dimensions. */
  if ((A.get_rows() != (int)(n + 1)) || (A.get_cols() != (int)(n + 1))) {
    critical("lattice_solve_reduced_basis_for_d_given_r(): "
      "Incorrect dimensions for the matrix A.");
  }

  if ((G.get_rows() != (int)(n + 1)) || (G.get_cols() != (int)(n + 1))) {
    critical("lattice_solve_reduced_basis_for_d_given_r(): "
      "Incorrect dimensions for the matrix G.");
  }

  /* Setup constants. */
  mpz_t pow2msigma;
  mpz_init(pow2msigma);
  mpz_set_ui(pow2msigma, 0);
  mpz_setbit(pow2msigma, parameters->m + parameters->sigma);

  /* Setup variables. */
  mpz_t candidate_d;
  mpz_init(candidate_d);

  mpz_t target_d;
  mpz_init(target_d);

  /* Setup the target vector. */
  vector<Z_NR<mpz_t>> target;
  target.resize(n + 1);

  for (uint32_t j = 0; j < n; j++) {
    mpz_mul(target[j].get_data(), ks[j], pow2msigma);
      /* target[j] = 2^(m + sigma) k */
    mpz_mul(target[j].get_data(), target[j].get_data(), parameters->r);
      /* target[j] = 2^(m + sigma) k r */
    mpz_neg(target[j].get_data(), target[j].get_data());
      /* target[j] = -2^(m + sigma) k r */
  }

  mpz_set_ui(target[n].get_data(), 0);

  /* Solve for the closest vector. */
  vector<Z_NR<mpz_t>> solution;

  /* Compute the closest vector. */
  babai_closest_vector(solution, target, G, A, n, precision);

  /* Compute the target. */
  mpz_set_ui(target_d, 0);
    /* target_d = 0 */
  mpz_setbit(target_d, parameters->sigma);
    /* target_d = 2^sigma */
  mpz_mul(target_d, target_d, parameters->d);
    /* target_d = 2^sigma d */
  mpz_mul(target_d, target_d, parameters->r);
    /* target_d = 2^sigma d r */

  /* Extract the candidate d. */
  mpz_abs(candidate_d, solution[n].get_data());

  if (mpz_cmp(candidate_d, target_d) == 0) {
    (*status_d) = LATTICE_STATUS_RECOVERED_IMMEDIATE;
  } else {
    /* Setup variables. */
    mpz_t multiple;
    mpz_init(multiple);

    mpz_t tmp_d;
    mpz_init(tmp_d);

    (*status_d) = LATTICE_STATUS_NOT_RECOVERED;

    mpz_set(candidate_d, solution[n].get_data());

    for (uint32_t i = 1; i <= SEARCH_BOUND_SHORTEST_VECTOR_MULTIPLE; i++) {
      mpz_mul_ui(multiple, A[0][n].get_data(), i);

      mpz_add(tmp_d, candidate_d, multiple);
      mpz_abs(tmp_d, tmp_d);
      if (mpz_cmp(tmp_d, target_d) == 0) {
        (*status_d) = LATTICE_STATUS_RECOVERED_SEARCH;
        break;
      }

      mpz_sub(tmp_d, candidate_d, multiple);
      mpz_abs(tmp_d, tmp_d);
      if (mpz_cmp(tmp_d, target_d) == 0) {
        (*status_d) = LATTICE_STATUS_RECOVERED_SEARCH;
        break;
      }
    }

    /* Clear memory. */
    mpz_clear(multiple);
    mpz_clear(tmp_d);
  }

  /* Clear memory. */
  target.clear();
  solution.clear();

  mpz_clear(pow2msigma);

  mpz_clear(target_d);
  mpz_clear(candidate_d);
}

void lattice_solve_for_d_given_r(
  Lattice_Status_Recovery * const status_d,
  const mpz_t * const js,
  const mpz_t * const ks,
  const int32_t * const etas,
  const uint32_t n,
  const Diagonal_Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision)
{
  /* Setup the A matrix. */
  ZZ_mat<mpz_t> A(n + 1, n + 1);

  /* Setup the Gram-Schmidt matrices. */
  FP_mat<mpfr_t> G(n + 1, n + 1);
  FP_mat<mpfr_t> M(n + 1, n + 1);

  /* Initially set the status to not recovered. */
  (*status_d) = LATTICE_STATUS_NOT_RECOVERED;

  while (TRUE) {
    /* Reduce the A matrix. */
    lattice_compute_reduced_diagonal_basis(
      A,
      js,
      etas,
      n,
      parameters,
      (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) ?
        REDUCTION_ALGORITHM_LLL : algorithm,
      2 * parameters->m); /* red_precision */

    if (LATTICE_STATUS_NOT_RECOVERED == (*status_d)) {
      /* Compute the Gram-Schmidt matrices. */
      gram_schmidt_orthogonalization(M, G, A, n, precision);

      /* Solve for d. */
      lattice_solve_reduced_basis_for_d_given_r(
        status_d,
        A,
        G,
        ks,
        n,
        parameters,
        precision);
    }

    if (LATTICE_STATUS_NOT_RECOVERED != (*status_d))
    {
      break;
    }

    if (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) {
      algorithm = REDUCTION_ALGORITHM_BKZ;
    } else {
      break;
    }
  }

  /* Clear memory. */
  A.clear();
  G.clear();
  M.clear();
}
