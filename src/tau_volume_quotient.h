/*!
 * \file    tau_volume_quotient.h
 * \ingroup estimating_volume_quotients
 *
 * \brief   The declaration of functions for converting tau estimates into
 *          volume quotient estimates.
 */

/*!
 * \defgroup estimating_volume_quotients Estimating volume quotients
 * \ingroup  math
 *
 * \brief    A module for functions for estimating volume quotients.
 */

#ifndef TAU_VOLUME_QUOTIENT_H
#define TAU_VOLUME_QUOTIENT_H

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>

/*!
 * \brief  Converts a tau estimate into a volume quotient estimate.
 *
 * Given
 *
 *    tau = log_( (1 / n) sum_{i = 1}^{n} alpha^2_i ) / 2 - m
 *
 * this function first computes the radius
 *
 *    R = sqrt( sum_{i = 1}^{n} alpha^2_i + ( d^2 or r^2) )
 *
 * and from it the volume quotient
 *
 *    v = pi^(D / 2) R^D / Gamma(D / 2 + 1) / 2^((l + m) n),
 *
 * where D = n + 1 is the dimension of the lattice, as described in the paper.
 *
 * \param[in] m       The parameter m.
 * \param[in] l       The parameter l.
 * \param[in] n       The number of runs n.
 * \param[in] tau     The tau estimate.
 * \param[in] d_or_r  The logarithm d or order r, depending on the context in
 *                    which the volume quotient is computed. May be set to NULL
 *                    to slightly over-estimate the volume quotient by forcing
 *                    the last component of the sum to 2^m.
 * \param[in, out] v  The volume quotient estimate.
 *
 * \return  #TRUE if the volume quotient is less than two, #FALSE otherwise.
 */
bool tau_volume_quotient(
  const uint32_t m,
  const uint32_t l,
  const uint32_t n,
  const mpz_t d_or_r,
  const long double tau,
  mpfr_t v);

/*!
 * \brief  Converts a tau estimate into a volume quotient estimate.
 *
 * Given
 *
 *    tau = log_( (1 / n) sum_{i = 1}^{n} alpha^2_i ) / 2 - (m + sigma - l)
 *
 * this function first computes the radius
 *
 *    R = sqrt( sum_{i = 1}^{n} alpha^2_i + 2^(sigma - l) d^2 )
 *
 * and from it the volume quotient
 *
 *    v = pi^(D / 2) R^D / Gamma(D / 2 + 1) / 2^((m + sigma) n + sigma - l),
 *
 * where D = n + 1 is the dimension of the lattice, as described in the paper.
 *
 * \param[in] m       The parameter m.
 * \param[in] sigma   The parameter sigma.
 * \param[in] l       The parameter l.
 * \param[in] n       The number of runs n.
 * \param[in] tau     The tau estimate.
 * \param[in] d       The logarithm d. May be set to NULL to slightly
 *                    over-estimate the volume quotient by forcing the last
 *                    component of the sum to 2^m.
 * \param[in, out] v  The volume quotient estimate.
 *
 * \return  #TRUE if the volume quotient is less than two, #FALSE otherwise.
 */
bool tau_volume_quotient_diagonal(
  const uint32_t m,
  const uint32_t sigma,
  const uint32_t l,
  const uint32_t n,
  const mpz_t d,
  const long double tau,
  mpfr_t v);

#endif /* TAU_VOLUME_QUOTIENT_H */
