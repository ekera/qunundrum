/*!
 * \file    probability.h
 * \ingroup two_dimensional_distribution_probability
 *
 * \brief   The declaration of functions for computing the probability of
 *          observing a pair (j, k) with angle pair (theta_d, theta_r).
 */

/*!
 * \defgroup two_dimensional_distribution_probability \
 *           Two-dimensional probability functions
 * \ingroup  two_dimensional_distribution
 *
 * \brief   A module for functions for computing the probability of observing a
 *          pair (j, k) with angle pair (theta_d, theta_r).
 */

#ifndef PROBABILITY_H
#define PROBABILITY_H

#include "parameters.h"

#include <mpfr.h>

#include <stdint.h>

/*!
 * \brief   Computes the approximate probability of observing a pair (j, k)
 *          with angle pair (theta_d, theta_r) for fixed sigma, and an upper
 *          bound on the error in the approximation.
 *
 * This function uses the approximation in the paper [1] by Ekerå on computing
 * general discrete logarithms and orders with tradeoffs. The parameter sigma
 * that enters into the approximation must be explicitly specified as an
 * argument to this function. It is not modified by this function.
 *
 * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms
 * and orders with tradeoffs. J. Math. Cryptol. 15, pp. 359–407 (2021).
 *
 * This function uses the default precision #PRECISION internally. There is no
 * need to go beyond the default precision in theta_d or theta_r when calling
 * the function.
 *
 * Note that this function is not suitable for evaluating the probability far
 * out on the diagonal where theta_d is very close to theta_r d/r. This would
 * require increasing the internal precision, and the precision in theta_d and
 * theta_r when calling the function. As a very small fraction of the
 * probability mass is far out on the diagonal, we can afford not to capture it.
 *
 * \param[in, out] norm     The approximate probability for the pair.
 * \param[in, out] error    An upper bound for the error in the approximation.
 * \param[in] sigma         The integer parameter sigma.
 * \param[in] theta_d       The angle theta_d.
 * \param[in] theta_r       The angle theta_r.
 * \param[in] parameters    The parameters for the probability distribution.
 *
 * \return  Returns #TRUE if the bound on the error in the approximation is at
 *          most one percent of the approximate probability, #FALSE otherwise.
 */
bool probability_approx(
  mpfr_t norm,
  mpfr_t error,
  const uint32_t sigma,
  const mpfr_t theta_d,
  const mpfr_t theta_r,
  const Parameters * const parameters);

/*!
 * \brief   Computes the approximate probability of observing a pair (j, k)
 *          with angle pair (theta_d, theta_r) for an initial guess of sigma
 *          that is adjusted by the function to locally minimize the error.
 *
 * This function uses the approximation in the paper by Ekerå [1] on computing
 * general discrete logarithms and orders with tradeoffs. An initial guess of
 * the parameter sigma that enters into the approximation must be explicitly
 * specified as an argument to this function. This function computes the
 * probability and an associated error bound by calling probability_approx().
 * It then checks, recursively, if increasing or decreasing sigma would yields a
 * better error bound. If this is the case, sigma is updated in place and the
 * probability and error bound for this sigma returned.
 *
 * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms
 * and orders with tradeoffs. J. Math. Cryptol. 15, pp. 359–407 (2021).
 *
 * This function calls probability_approx(). See the documentation of this
 * function for information on precision requirements on theta_d and theta_r.
 *
 * \param[in, out] norm   The approximate probability for the pair.
 * \param[in, out] error  An upper bound for the error in the approximation.
 * \param[in, out] sigma  An initial guess for the integer parameter sigma that
 *                        is adjusted by this function if necessary.
 * \param[in] theta_d     The angle theta_d.
 * \param[in] theta_r     The angle theta_r.
 * \param[in] parameters  The parameters for the probability distribution.
 *
 * \return  Returns #TRUE if the bound on the error in the approximation is at
 *          most one percent of the approximate probability, #FALSE otherwise.
 */
bool probability_approx_adjust_sigma(
  mpfr_t norm,
  mpfr_t error,
  uint32_t &sigma,
  const mpfr_t theta_d,
  const mpfr_t theta_r,
  const Parameters * const parameters);

/*!
 * \brief   Computes the approximate probability of observing a pair (j, k)
 *          with angle pair (theta_d, theta_r) for optimal sigma in a range.
 *
 * This function uses the approximation [1] in the paper by Ekerå on computing
 * general discrete logarithms and orders with tradeoffs. This function computes
 * the probability and an associated error bound by calling probability_approx()
 * for all permissible values of the parameter sigma. It returns the value of
 * sigma that yields the best error bound, and the probability associated with
 * this value of sigma.
 *
 * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms
 * and orders with tradeoffs. J. Math. Cryptol. 15, pp. 359–407 (2021).
 *
 * This function calls probability_approx(). See the documentation of this
 * function for information on precision requirements on theta_d and theta_r.
 *
 * \param[in, out] norm   The approximate probability for the pair.
 * \param[in, out] error  An upper bound for the error in the approximation.
 * \param[out] sigma      The optimal sigma.
 * \param[in] theta_d     The angle theta_d.
 * \param[in] theta_r     The angle theta_r.
 * \param[in] parameters  The parameters for the probability distribution.
 *
 * \return  Returns #TRUE if the bound on the error in the approximation is at
 *          most one percent of the approximate probability, #FALSE otherwise.
 */
bool probability_approx_optimal_sigma(
  mpfr_t norm,
  mpfr_t error,
  uint32_t &sigma,
  const mpfr_t theta_d,
  const mpfr_t theta_r,
  const Parameters * const parameters);

/*!
 * \brief   Computes the approximate probability of observing a pair (j, k)
 *          with angle pair (theta_d, theta_r).
 *
 * This function uses the quick and dirty approximation in the introduction to
 * the paper [1] by Ekerå on computing general discrete logarithms and orders
 * with tradeoffs. It produces a good approximation in most cases. However, it
 * is seemingly non-trivial to bound the error in the approximation.
 *
 * This function is present mainly for completeness and verification purposes.
 *
 * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms
 * and orders with tradeoffs. J. Math. Cryptol. 15, pp. 359–407 (2021).
 *
 * This function uses the default precision #PRECISION internally. There is no
 * need to go beyond the default precision in theta_d or theta_r when calling
 * the function.
 *
 * Note that this function is not suitable for evaluating the probability far
 * out on the diagonal where theta_d is very close to theta_r d/r. This would
 * require increasing the internal precision, and the precision in theta_d and
 * theta_r when calling the function. As a very small fraction of the
 * probability mass is far out on the diagonal, we can afford not to capture it.
 *
 * \param[in, out] norm   The approximate probability for the pair.
 * \param[in] theta_d     The angle theta_d.
 * \param[in] theta_r     The angle theta_r.
 * \param[in] parameters  The parameters for the probability distribution.
 */
void probability_approx_quick(
  mpfr_t norm,
  const mpfr_t theta_d,
  const mpfr_t theta_r,
  const Parameters * const parameters);

#endif /* PROBABILITY_H */
