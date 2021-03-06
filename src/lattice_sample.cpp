/*!
 * \file    lattice_sample.cpp
 * \ingroup lattice_sample
 *
 * \brief   The definition of functions for mapping any argument pair to the
 *          closest admissible argument pair in the argument plane.
 */

#include "lattice_sample.h"

#include "lattice_gso.h"
#include "lattice_babai.h"

#include "parameters.h"

#include "math.h"
#include "errors.h"

#include <fplll/fplll.h>

#include <mpfr.h>
#include <gmp.h>

#include <vector>

#include <stdint.h>

using namespace fplll;
using namespace std;

void lattice_alpha_init(
  Lattice_Alpha * const lattice,
  const Parameters * const parameters)
{
  /* Compute kappa_d and kappa_r. */
  const uint32_t kappa_d = kappa(parameters->d);
  const uint32_t kappa_r = kappa(parameters->r);

  /* Compute gamma. (Should be zero almost all the time.) */
  const uint32_t gamma =
    (uint32_t)max_i(0, (int32_t)kappa_r - (int32_t)(parameters->l + kappa_d));

  /* Compute 2^(m - gamma). */
  mpz_t pow_m_gamma;
  mpz_init_set_ui(pow_m_gamma, 0);
  mpz_setbit(pow_m_gamma, parameters->m - gamma);

  /* Compute 2^kappa_r. */
  mpz_t pow_kappa_r;
  mpz_init_set_ui(pow_kappa_r, 0);
  mpz_setbit(pow_kappa_r, kappa_r);

  /* Compute delta_r. */
  mpz_t delta_r;
  mpz_init(delta_r);
  mpz_div(delta_r, parameters->r, pow_kappa_r);
  mpz_invert(delta_r, delta_r, pow_m_gamma);
  mpz_mul(delta_r, parameters->d, delta_r);
  mpz_mod(delta_r, delta_r, pow_m_gamma);

  /* Setup the basis matrix. */
  lattice->A.resize(2, 2);
  mpz_set(lattice->A[0][0].get_data(), delta_r);
  mpz_set(lattice->A[0][1].get_data(), pow_kappa_r);
  mpz_set(lattice->A[1][0].get_data(), pow_m_gamma);
  mpz_set_ui(lattice->A[1][1].get_data(), 0);

  /* Reduce the basis matrix. */
  int status = hkz_reduction(lattice->A, HKZ_DEFAULT, FT_MPFR, parameters->m);
  if (RED_SUCCESS != status) {
    critical("lattice_alpha_init(): Failed to reduce matrix A using HKZ.");
  }

  /* Ideally, we would use Gauss-Lagrange for this two-dimensional matrix, but
   * HKZ should be equivalent for such matrices, so this does not matter. */

  /* Compute a Gram-Schmidt orthogonalized basis. */
  FP_mat<mpfr_t> M(2, 2);

  lattice->G.resize(2, 2);

  gram_schmidt_orthogonalization(M, lattice->G, lattice->A, 1, parameters->m);

  /* Clean up. */
  M.clear();

  mpz_clear(pow_m_gamma);
  mpz_clear(pow_kappa_r);
  mpz_clear(delta_r);
}

void lattice_alpha_clear(
  Lattice_Alpha * const lattice)
{
  lattice->A.clear();
  lattice->G.clear();
}

void lattice_alpha_map(
  mpz_t alpha_d,
  mpz_t alpha_r,
  const Lattice_Alpha * const lattice,
  const Parameters * const parameters)
{
  /* Setup the target vector. */
  vector<Z_NR<mpz_t>> target;
  target.resize(2);

  mpz_set(target[0].get_data(), alpha_d);
  mpz_set(target[1].get_data(), alpha_r);

  /* Execute Babai's algorithm to find the closest vector to the target. */
  vector<Z_NR<mpz_t>> solution;
  solution.resize(2);

  babai_closest_vector(
    solution,
    target,
    lattice->G,
    lattice->A,
    1, /* n */
    parameters->m);

  /* Export the solution. */
  mpz_set(alpha_d, solution[0].get_data());
  mpz_set(alpha_r, solution[1].get_data());

  /* Clean up. */
  target.clear();
  solution.clear();
}
