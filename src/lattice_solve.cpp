/*!
 * \file    lattice_solve.cpp
 * \ingroup lattice_solve
 *
 * \brief   The definition of functions for solving for d and r using Babai's
 *          nearest plane algorithm and lattice basis reduction techniques.
 */

#include "lattice_solve.h"

#include "lattice.h"
#include "lattice_gso.h"
#include "lattice_babai.h"
#include "lattice_algebra.h"
#include "lattice_reduce.h"

#include "parameters.h"
#include "sample.h"
#include "errors.h"
#include "timer.h"
#include "math.h"

#include <fplll/fplll.h>

#include <gmp.h>
#include <mpfr.h>

#include <vector>

#include <stdint.h>

using namespace std;
using namespace fplll;

/*!
 * \brief A bound on the largest multiple of the last component of the 
 *        shortest vector in the reduced basis matrix of the lattice to be 
 *        added to or subtracted from the candidate for d.
 * 
 * The last component of the shortest vector in the reduced basis matrix of 
 * the lattice is expected to be a divisor of r.
 */
#define SEARCH_BOUND_SHORTEST_VECTOR_MULTIPLE 256

void lattice_solve_reduced_basis_for_d(
  Lattice_Status_Recovery * const status_d,
  const ZZ_mat<mpz_t> &A,
  const FP_mat<mpfr_t> &G,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  const uint32_t precision)
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

  #ifdef CLOSEST_VECTOR_FPLLL
    vector<Z_NR<mpz_t>> solution_coordinates;

    status = closest_vector(A, target, solution_coordinates, flags);
    if (0 != status) {
      critical("lattice_solve_reduced_basis_for_d(): "
        "Failed to find the closest vector.");
    }

    vector_matrix_product(solution, solution_coordinates, A);

    /* Clear memory. */
    solution_coordinates.clear();
  #else
    /* Compute the closest vector. */
    babai_closest_vector(solution, target, G, A, n, precision);
  #endif

  /* Extract the candidate d. */
  mpz_abs(candidate_d, solution[n].get_data());

  if (mpz_cmp(candidate_d, parameters->d) == 0) {
    (*status_d) = LATTICE_STATUS_RECOVERED_IMMEDIATE;
  } else {
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

  /* Clear memory. */
  target.clear();
  solution.clear();

  mpz_clear(pow2m);
  mpz_clear(pow2lm);
  mpz_clear(pow2lm_half);
  mpz_clear(candidate_d);
}

void lattice_solve_for_d(
  Lattice_Status_Recovery * const status_d,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  uint32_t precision)
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
      precision);

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

void lattice_solve_reduced_basis_for_r(
  Lattice_Status_Recovery * const status_r,
  const ZZ_mat<mpz_t> &A,
  const uint32_t n,
  const Parameters * const parameters)
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
  mpz_mod(test_r, parameters->r, candidate_r);

  if (mpz_cmp_ui(test_r, 0) != 0) {
    (*status_r) = LATTICE_STATUS_NOT_RECOVERED;
  } else {
    mpz_div(test_r, parameters->r, candidate_r); /* r / candidate_r */

    if (mpz_cmp_ui(test_r, 1) == 0) {
      (*status_r) = LATTICE_STATUS_RECOVERED_IMMEDIATE;
    } else if (mpz_cmp_ui(test_r, SEARCH_BOUND_SHORTEST_VECTOR_MULTIPLE) <= 0) {
      (*status_r) = LATTICE_STATUS_RECOVERED_SEARCH;
    } else {
      (*status_r) = LATTICE_STATUS_NOT_RECOVERED;
    }
  }

  /* Clean up memory. */
  mpz_clear(candidate_r);
  mpz_clear(test_r);
}

void lattice_solve_for_r(
  Lattice_Status_Recovery * const status_r,
  const mpz_t * const js,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm)
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
      parameters);

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
  const uint32_t precision)
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
        parameters);
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
        precision);
    }

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
