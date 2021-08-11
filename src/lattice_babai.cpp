/*!
 * \file    lattice_babai.cpp
 * \ingroup lattice
 *
 * \brief   The definition of functions for implementing Babai's nearest plane
 *          algorithm that may be used to find the closest vector to a target
 *          vector in an integer lattice.
 */

#include "lattice_babai.h"

#include "errors.h"
#include "lattice.h"

#include <gmp.h>
#include <mpfr.h>

#include <fplll/fplll.h>

#include <stdint.h>

#include <vector>

using namespace std;
using namespace fplll;

void babai_closest_vector(
  vector<Z_NR<mpz_t>> &solution,
  const vector<Z_NR<mpz_t>> &target,
  const FP_mat<mpfr_t> &G,
  const ZZ_mat<mpz_t> &A,
  const uint32_t n,
  const uint32_t precision)
{
  /* Sanity checks of the dimensions of the A and G matrices. */
  if ((A.get_rows() != (int)(n + 1)) || (A.get_cols() != (int)(n + 1))) {
    critical("babai_closest_vector(): Matrix A has incorrect dimensions.");
  }

  if ((G.get_rows() != (int)(n + 1)) || (G.get_cols() != (int)(n + 1))) {
    critical("babai_closest_vector(): Matrix G has incorrect dimensions.");
  }

  /* Sanity checks of the size of the target vector. */
  if (target.size() != (n + 1)) {
    critical("babai_closest_vector(): "
      "The target vector has incorrect dimensions.");
  }

  /* Resize the solution vector if necessary. */
  if (solution.size() != (n + 1)) {
    solution.resize(n + 1);
  }

  /* Initialize temporary variables. */
  mpfr_t x;
  mpfr_init2(x, precision);

  mpfr_t sum;
  mpfr_init2(sum, precision);

  mpfr_t factor_f;
  mpfr_init2(factor_f, precision);
  
  mpz_t tmp_z;
  mpz_init(tmp_z);

  mpz_t factor_z;
  mpz_init(factor_z);

  /* Copy the target to the solution. */
  for (uint32_t j = 0; j <= n; j++) {
    mpz_set(solution[j].get_data(), target[j].get_data());
  }

  for (int32_t i = (int32_t)n; i >= 0; i--) {
    /* Compute the projection factor. */
    mpfr_set_ui(sum, 0, MPFR_RNDN);

    for (uint32_t j = 0; j <= n; j++) {
      mpfr_mul_z(x, G[i][j].get_data(), solution[j].get_data(), MPFR_RNDN);
      mpfr_add(sum, sum, x, MPFR_RNDN);
    }

    mpfr_set(factor_f, sum, MPFR_RNDN);

    mpfr_set_ui(sum, 0, MPFR_RNDN);

    for (uint32_t j = 0; j <= n; j++) {
      mpfr_mul(x, G[i][j].get_data(), G[i][j].get_data(), MPFR_RNDN);
      mpfr_add(sum, sum, x, MPFR_RNDN);
    }

    mpfr_div(factor_f, factor_f, sum, MPFR_RNDN);
    mpfr_round(factor_f, factor_f);
    mpfr_get_z(factor_z, factor_f, MPFR_RNDN);

    /* Subtract the projection factor times the basis vector. */
    for (uint32_t j = 0; j <= n; j++) {
      mpz_mul(tmp_z, factor_z, A[i][j].get_data());
      mpz_sub(solution[j].get_data(), solution[j].get_data(), tmp_z);
    }
  }

  for (uint32_t j = 0; j <= n; j++) {
    mpz_sub(solution[j].get_data(),
            target[j].get_data(),
            solution[j].get_data());
  }

  /* Clear memory. */
  mpfr_clear(x);
  mpfr_clear(sum);
  mpfr_clear(factor_f);

  mpz_clear(tmp_z);
  mpz_clear(factor_z);
}
