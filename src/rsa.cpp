/*!
 * \file    rsa.cpp
 * \ingroup rsa
 *
 * \brief   The definition of functions for generating RSA moduli.
 */

#include "rsa.h"

#include "common.h"
#include "errors.h"
#include "random.h"

#include <gmp.h>

#include <stdint.h>

void rsa_generate_modulus(
  mpz_t p,
  mpz_t q,
  const uint32_t n,
  const bool check_modulus_size,
  Random_State * const random_state)
{
  if ((n < 16) || ((n % 2) != 0)) {
    critical("rsa_generate_modulus(): Incorrect modulus length.");
  }

  mpz_t N;
  mpz_init(N);

  mpz_t modulus;
  mpz_init(modulus);
  mpz_set_ui(modulus, 0);
  mpz_setbit(modulus, n / 2);

  int result;

  while (TRUE) {
    mpz_set_ui(p, 0);
    random_generate_mpz(p, modulus, random_state);
    mpz_setbit(p, n / 2 - 1);
    mpz_setbit(p, 0);

    result = mpz_probab_prime_p(p, RSA_MILLER_RABIN_ITERATIONS);
    /* Note: Returns 2 if n is prime, return 1 if n is probably prime (without
     * being certain), or return 0 if n is definitely non-prime. */
    if ((1 != result) && (2 != result)) {
      continue;
    }

    while (TRUE) {
      mpz_set_ui(q, 0);
      random_generate_mpz(q, modulus, random_state);
      mpz_setbit(q, n / 2 - 1);
      mpz_setbit(q, 0);

      result = mpz_probab_prime_p(q, RSA_MILLER_RABIN_ITERATIONS);
      /* Note: Returns 2 if n is prime, return 1 if n is probably prime (without
       * being certain), or return 0 if n is definitely non-prime. */
      if ((1 != result) && (2 != result)) {
        continue;
      }

      /* Check that p is not equal to q. */
      if (0 == mpz_cmp(p, q)) {
        continue;
      }

      break;
    }

    if (TRUE == check_modulus_size) {
      mpz_mul(N, p, q);
      if (n != mpz_sizeinbase(N, 2)) {
        continue;
      }
    }

    break;
  }

  /* Clear memory. */
  mpz_clear(modulus);
  mpz_clear(N);
}
