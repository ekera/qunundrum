/*!
 * \file    diagonal_probability.h
 * \ingroup diagonal_distribution_probability
 *
 * \brief   The declaration of functions for computing the probability of
 *          observing a pair (j, k) with angle pair (theta_d, theta_r) for
 *          a general discrete logarithm.
 */

/*!
 * \defgroup diagonal_distribution_probability \
 *           Diagonal probability functions
 * \ingroup  diagonal_distribution
 *
 * \brief   A module for functions for computing the probability of observing a
 *          pair (j, k) with angle pair (theta_d, theta_r).
 */

#ifndef DIAGONAL_PROBABILITY_H
#define DIAGONAL_PROBABILITY_H

#include "diagonal_parameters.h"

#include <mpfr.h>

/*!
 * \brief   Computes the approximate probability f_eta(theta_r) of observing a
 *          pair (j, k) with angle theta_r summed over all angles theta_d, where
 *          eta is the peak index.
 *
 * This function uses the heuristic approximation f in the paper [1] by Ekerå on
 * revisiting Shor's algorithm for general discrete logarithms.
 *
 * [1] Ekerå, M.: Revisiting Shor's quantum algorithm for computing general
 * discrete logarithms. In: ArXiv Pre-Print 1905.09084v3.
 *
 * To ensure numeric stability, this function uses higher precision than the
 * default precision of #PRECISION internally. However, there is no need to go
 * beyond the default precision in theta_r when calling the function.
 *
 * \param[in, out] norm     The approximate probability for the pair.
 * \param[in] theta_r       The angle theta_r.
 * \param[in] eta           The peak index eta.
 * \param[in] parameters    The parameters for the probability distribution.
 */
void diagonal_probability_approx_f_eta(
  mpfr_t norm,
  const mpfr_t theta_r,
  const int32_t eta,
  const Diagonal_Parameters * const parameters);

/*!
 * \brief   Computes h(phi).
 *
 * This function uses the heuristic approximation h in the paper [1] by Ekerå on
 * revisiting Shor's algorithm for general discrete logarithms.
 *
 * [1] Ekerå, M.: Revisiting Shor's quantum algorithm for computing general
 * discrete logarithms. In: ArXiv Pre-Print 1905.09084v3.
 *
 * Since the probability is high when theta_r d/r is very close to theta_d, this
 * function uses higher precision than the default precision of #PRECISION
 * internally to ensure numeric stability.
 *
 * You typically need to set a high precision, on the order of l bits, in phi
 * when calling this function. This function has been tested with a precision
 * of 3l bits in phi. This is the recommended precision setting.
 *
 * \param[in, out] norm     The approximate probability for the pair.
 * \param[in] phi           The angle phi.
 * \param[in] parameters    The parameters for the probability distribution.
 */
void diagonal_probability_approx_h(
  mpfr_t norm,
  const mpfr_t phi,
  const Diagonal_Parameters * const parameters);

#endif /* DIAGONAL_PROBABILITY_H */
