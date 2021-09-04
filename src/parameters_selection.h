/*!
 * \file    parameters_selection.h
 * \ingroup parameters
 *
 * \brief   The definition of functions for selecting the discrete logarithm d
 *          and order r.
 */

#ifndef PARAMETER_SELECTION_H
#define PARAMETER_SELECTION_H

#include "random.h"

#include <gmp.h>

#include <stdint.h>

/*!
 * \brief   Deterministically selects the logarithm d and the order r as m bit 
 *          integers by reading digits from Catalan's constant
 * 
 * \param[in, out] d    The logarithm d.
 * \param[in, out] r    The order r.
 * \param[in] m         The bit length of d and r.
 */
void parameters_selection_deterministic_d_r(
  mpz_t d,
  mpz_t r,
  const uint32_t m);

/*!
 * \brief   Selects the logarithm d or the order r uniformly at random from
 *          the interval (2^(m-1), 2^m).
 *
 * \param[in, out] value    The logarithm d or order r.
 * \param[in] m             The bit length of the value.
 */
void parameters_selection_random_d_or_r(
  mpz_t value,
  const uint32_t m);

/*!
 * \brief   Selects the order r uniformly at random from the interval
 *          (2^(m-1), 2^m) and the logarithm d uniformly at random from the
 *          interval [r/2, r).
 *
 * \param[in, out] d        The logarithm d.
 * \param[in, out] r        The order r.
 * \param[in] m             The bit length of the order r.
 */
void parameters_selection_random_d_and_r(
  mpz_t d,
  mpz_t r,
  const uint32_t m);

/*!
 * \brief   Selects the logarithm d from a given RSA modulus N = pq of length n
 *          bits, by setting d = (p + q - 2) / 2 - 2^(n / 2 - 1).
 *
 * This function requires that 2^(n-1) <= p, q < 2^n.
 *
 * It furthermore requires that n is even and that n >= 16.
 *
 * \param[in, out] d        The logarithm d.
 * \param[in] p             The prime p.
 * \param[in] q             The prime q.
 * \param[in] n             The bit length of the modulus N = pq.
 */
void parameters_selection_explicit_rsa(
  mpz_t d,
  const mpz_t p,
  const mpz_t q,
  const uint32_t n);

/*!
 * \brief   Selects the logarithm d by first selecting an RSA modulus N = pq
 *          of length n bits, and setting d = (p + q - 2) / 2 - 2^(n / 2 - 1).
 *
 * This function requires that n is even and that n >= 16.
 *
 * This conveniency function calls rsa_generate_modulus() to generate p and q, 
 * and parameters_selection_explicit_rsa() to select the logarithm from p and q.
 *
 * \param[in, out] d        The logarithm d.
 * \param[in, out] p        The prime p.
 * \param[in, out] q        The prime q.
 * \param[in] n             The bit length of the modulus N = pq.
 */
void parameters_selection_random_rsa(
  mpz_t d,
  mpz_t p,
  mpz_t q,
  const uint32_t n);

#endif /* PARAMETER_SELECTION_H */
