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

#include "parameters.h"

#include <mpfr.h>

/*!
 * \brief   Computes the approximate probability of observing a pair (j, k)
 *          with angle pair (theta_d, theta_r).
 *
 * This function uses the heuristic approximation in the paper [1] by Ekerå on
 * revisiting Shor's algorithm for general discrete logarithms. Note that there 
 * is currently no error bound available for this approximation.
 * 
 * [1] Ekerå, M.: Revisiting Shor's quantum algorithm for computing general 
 * discrete logarithms. In: ArXiv Pre-Print 1905.09084.
 * 
 * Since the probability is high when theta_r d/r is very close to theta_d, this
 * function uses higher precision than the default precision of #PRECISION 
 * internally to ensure numeric stability.
 * 
 * You typically need to set a high precision, on the order of m + l bits, in 
 * theta_d and theta_r when calling this function. This function has been tested
 * when a precision of 2 (m + l) bits in theta_d and theta_r respectively. This
 * is the recommended precision setting.
 *
 * \param[in, out] norm     The approximate probability for the pair.
 * \param[in] theta_d       The angle theta_d.
 * \param[in] theta_r       The angle theta_r.
 * \param[in] parameters    The parameters for the probability distribution.
 */
void diagonal_probability_approx(
  mpfr_t norm,
  const mpfr_t theta_d,
  const mpfr_t theta_r,
  const Parameters * const parameters);

#endif /* DIAGONAL_PROBABILITY_H */
