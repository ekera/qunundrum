/*!
 * \file    lattice_smoothness.cpp
 * \ingroup lattice_solve
 *
 * \brief   The definition of functions for checking if an integer is smooth 
 *          and for removing smooth factors from partially smooth integers.
 */

#include "lattice_smoothness.h"

#include "common.h"

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>

void lattice_smoothness_remove_smooth_factors(
  mpz_t reduced_z,
  const mpz_t z,
  const double c,
  const uint32_t m)
{
  /* Setup variables. */
  mpz_t q;
  mpz_init(q);

  mpz_t pow_q;
  mpz_init(pow_q);

  mpz_t tmp_z;
  mpz_init(tmp_z);

  mpfr_t B;
  mpfr_init(B);

  mpfr_t tmp_f;
  mpfr_init(tmp_f);

  /* Set B = 2^(cm). */
  mpfr_set_d(B, c * m, MPFR_RNDN);
  mpfr_set_ui(tmp_f, 2, MPFR_RNDN);
  mpfr_pow(B, tmp_f, B, MPFR_RNDN);

  /* Reduce r. */
  mpz_set(reduced_z, z);
  mpz_set_ui(q, 2);

  while (TRUE) {
    mpz_mod(tmp_z, reduced_z, q);
    if (mpz_cmp_ui(tmp_z, 0) == 0) {
      mpz_set(pow_q, q);

      while (TRUE) {
        mpz_mul(tmp_z, pow_q, q);

        if (mpz_sizeinbase(tmp_z, 2) >= m) {
          break;
        }

        mpz_mod(tmp_z, reduced_z, tmp_z);
        if (mpz_cmp_ui(tmp_z, 0) != 0) {
          break;
        }

        mpz_mul(pow_q, pow_q, q);
      }

      mpz_div(reduced_z, reduced_z, pow_q);
    }

    /* Process the next prime. */
    mpz_nextprime(q, q);

    if (mpz_cmp_d(q, c * m) > 0) {
      break;
    }
  }

  /* Clear memory. */
  mpz_clear(q);
  mpz_clear(pow_q);
  mpz_clear(tmp_z);

  mpfr_clear(B);
  mpfr_clear(tmp_f);
}

bool lattice_smoothness_is_smooth(
  const mpz_t z,
  const double c,
  const uint32_t m)
{
  /* Setup variables. */
  mpz_t reduced_z;
  mpz_init(reduced_z);

  /* Remove smooth cofactors from z. */
  lattice_smoothness_remove_smooth_factors(reduced_z, z, c, m);

  /* Check if reduced z is one. */
  const bool result = (mpz_cmp_ui(reduced_z, 1) == 0) ? TRUE : FALSE;

  /* Clear memory. */
  mpz_clear(reduced_z);

  /* Signal result. */
  return result;
}
