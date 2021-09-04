/*!
 * \file    parameters_selection.cpp
 * \ingroup parameters
 *
 * \brief   The definition of functions for selecting the discrete logarithm d
 *          and order r.
 */

#include "parameters_selection.h"

#include "common.h"
#include "errors.h"
#include "random.h"
#include "rsa.h"

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>

void parameters_selection_deterministic_d_r(
  mpz_t d,
  mpz_t r,
  const uint32_t m)
{
  if (m > 8192) {
    critical("parameters_selection_deterministic_d_r(): "
      "The parameter m must be less than or equal to 8192.");
  }

  mpz_t data_z;
  mpz_init(data_z);

  {
    const uint32_t precision = 3 * 8192;

    mpfr_t data;
    mpfr_init2(data, precision);
    mpfr_const_catalan(data, MPFR_RNDN);

    mpfr_t pow2;
    mpfr_init2(pow2, precision);
    mpfr_set_ui_2exp(pow2, 1, (mpfr_exp_t)precision, MPFR_RNDN);
    mpfr_mul(data, data, pow2, MPFR_RNDN);
    mpfr_floor(data, data);

    mpfr_get_z(data_z, data, MPFR_RNDN);

    /* Clear memory. */
    mpfr_clear(data);
    mpfr_clear(pow2);
  }

  uint32_t length = mpz_sizeinbase(data_z, 2);

  mpz_set_ui(r, 0);
  for (uint32_t i = 0; i < m - 1; i++) {
    if (0 != mpz_tstbit(data_z, (length - 1) - i)) {
      mpz_setbit(r, m - 2 - i);
    }
  }

  mpz_set_ui(d, 0);
  for (uint32_t i = 0; i < m - 1; i++) {
    if (0 != mpz_tstbit(data_z, (length - 1) - (8192 - 1) - i)) {
      mpz_setbit(d, m - 2 - i);
    }
  }

  mpz_mod(d, d, r);

  mpz_setbit(r, m - 1);
  mpz_setbit(d, m - 1);

  /* Clear memory. */
  mpz_clear(data_z);
}

void parameters_selection_random_d_or_r(
  mpz_t value,
  const uint32_t m)
{
  /* Setup constants. */
  mpz_t half_modulus_minus_one;
  mpz_init(half_modulus_minus_one);
  mpz_set_ui(half_modulus_minus_one, 0);
  mpz_setbit(half_modulus_minus_one, m - 1);
  mpz_sub_ui(half_modulus_minus_one, half_modulus_minus_one, 1);

  /* Setup a random state. */
  Random_State random_state;
  random_init(&random_state);

  /* Select value uniformly at random from [0, 2^(m-1) - 2]. */
  random_generate_mpz(value, half_modulus_minus_one, &random_state);

  /* Add 2^(m-1) + 1 to value, so that it is on
   *
   *    [2^(m-1) + 1, 2^(m-1) - 2 + 2^(m-1) + 1] =
   *      [2^(m-1) + 1, 2^m - 1] = (2^(m-1), 2^m).
   */
  mpz_add(value, value, half_modulus_minus_one);
  mpz_add_ui(value, value, 2);

  /* Clear memory. */
  random_close(&random_state);

  mpz_clear(half_modulus_minus_one);
}

void parameters_selection_random_d_and_r(
  mpz_t d,
  mpz_t r,
  const uint32_t m)
{
  /* Select r uniformly at random from the interval (2^(m-1), 2^m). */
  parameters_selection_random_d_or_r(r, m);

  /* Select d uniformly at random from the interval [r / 2, r). */
  mpz_t half_r;
  mpz_init(half_r);
  mpz_set(half_r, r);
  mpz_div_ui(half_r, half_r, 2);

  /* Setup a random state. */
  Random_State random_state;
  random_init(&random_state);

  /* Randomize d. */
  random_generate_mpz(d, half_r, &random_state);
  mpz_add(d, d, half_r);

  /* Clear memory. */
  random_close(&random_state);

  mpz_clear(half_r);
}

void parameters_selection_random_rsa(
  mpz_t d,
  mpz_t p,
  mpz_t q,
  const uint32_t n)
{
  if ((n < 16) || ((n % 2) != 0)) {
    critical("parameters_selection_random_rsa(): "
      "The modulus length must be even and greater than or equal to 16.");
  }

  /* Setup a random state. */
  Random_State random_state;
  random_init(&random_state);

  rsa_generate_modulus(p, q, n, TRUE, &random_state);

  parameters_selection_explicit_rsa(d, p, q, n);

  /* Clear memory. */
  random_close(&random_state);
}

void parameters_selection_explicit_rsa(
  mpz_t d,
  const mpz_t p,
  const mpz_t q,
  const uint32_t n)
{
  mpz_t pow;
  mpz_init(pow);

  if ((n < 16) || ((n % 2) != 0)) {
    critical("parameters_selection_explicit_rsa(): "
      "The modulus length must be even and greater than or equal to 16.");
  }

  if ((n != 2 * mpz_sizeinbase(p, 2)) || (n != 2 * mpz_sizeinbase(q, 2))) {
    critical("parameters_selection_explicit_rsa(): "
      "Incorrect size of p and/or q given n.");
  }

  /* Compute d. */
  mpz_set_ui(pow, 0);
  mpz_setbit(pow, n / 2 - 1);

  mpz_add(d, p, q);
  mpz_sub_ui(d, d, 2);
  mpz_div_ui(d, d, 2);
  mpz_sub(d, d, pow);

  /* Sanity check. */
  if (mpz_sizeinbase(d, 2) >= n / 2) {
    /* Should never occur, given the above tests. */
    critical("parameters_selection_explicit_rsa(): "
      "Unexpected result. Incorrect p and/or q or internal error.");
  }

  /* Clear memory. */
  mpz_clear(pow);
}
