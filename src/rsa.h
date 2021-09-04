/*!
 * \file    rsa.h
 * \ingroup rsa
 *
 * \brief   The declaration of functions for generating RSA moduli.
 */

/*!
 * \defgroup rsa RSA
 * \ingroup  math
 *
 * \brief    A module for functions for generating RSA moduli.
 */

#ifndef RSA_H
#define RSA_H

#include "random.h"

#include <gmp.h>

#include <stdint.h>

/*!
 * \brief   The number of iterations to run of the Miller–Rabin probabilistic
 *          primality test when generating and verifying RSA moduli.
 */
#define RSA_MILLER_RABIN_ITERATIONS       200

/*!
 * \brief   Generates a random RSA modulus N = pq bits.
 *
 * This function requires that n is even and that n >= 16.
 *
 * This function guarantees that 2^(n / 2 - 1) <= p, q < 2^(n/2).
 *
 * If the check_modulus_size flag is set to #TRUE this function also guarantees
 * that the modulus is of length n bit.
 *
 * \param[in, out] p              The prime factor p.
 * \param[in, out] q              The prime factor q.
 * \param[in] n                   The length parameter n.
 * \param[in] check_modulus_size  A flag that should be set to #TRUE if the
 *                                modulus must be of length exactly n bit, and
 *                                to #FALSE otherwise.
 * \param[in, out] random_state   The random state to use to generate p and q.
 */
void rsa_generate_modulus(
  mpz_t p,
  mpz_t q,
  const uint32_t n,
  const bool check_modulus_size,
  Random_State * const random_state);

#endif /* RSA_H */
