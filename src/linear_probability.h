/*!
 * \file    linear_probability.h
 * \ingroup linear_distribution_probability
 *
 * \brief   The declaration of functions for computing the probability of
 *          observing a pair (j, k) with angle theta_d, in the short discrete
 *          logarithm case, or an integer j with angle theta_r in the order
 *          case.
 */

/*!
 * \defgroup linear_distribution_probability Linear probability functions
 * \ingroup  linear_distribution
 *
 * \brief    A module for functions for computing the probability of observing a
 *           pair (j, k) with angle theta_d, in the short discrete logarithm
 *           case, or an integer j with angle theta_r in the order case.
 */

#ifndef LINEAR_PROBABILITY_H
#define LINEAR_PROBABILITY_H

#include "parameters.h"

#include <mpfr.h>

/*!
 * \brief   Computes the exact probability of observing a pair (j, k) with
 *          angle theta_d using a closed-form expression.
 * 
 * This function uses the closed-form expression in the paper [1] by Ekerå on 
 * post-processing in Ekerå's and Håstad's quantum algorithm [2] for computing 
 * short discrete logarithms and factoring RSA integers.
 * 
 * [1] Ekerå, M.: On post-processing in the quantum algorithm for computing 
 * short discrete logarithms. In: IACR ePrint Archive, 2017/1122.
 * 
 * [2] Ekerå, M. and Håstad, J.: Quantum algorithms for computing short discrete
 * logarithms and factor RSA integers. In: PQCrypto 2017, Springer LNCS 10346,
 * pp. 347-363 (2017). DOI: https://doi.org/10.1007/978-3-319-59879-6_20.
 * 
 * To ensure numeric stability, this function uses higher precision than the
 * default precision of #PRECISION internally. However, there is no need to go 
 * beyond the default precision in theta_d when calling the function.
 * 
 * \param[in, out] norm   The approximate probability for the pair.
 * \param[in] theta_d     The angle theta_d.
 * \param[in] parameters  The parameters for the probability distribution.
 */
void linear_probability_d(
  mpfr_t norm,
  const mpfr_t theta_d,
  const Parameters * const parameters);

/*!
 * \brief   Computes the exact probability of observing an integer j with angle
 *          theta_r using a closed-form expression.
 *
 * This function uses the closed-form expression in appendix A to the paper [1] 
 * by Ekerå on computing general discrete logarithms and orders with tradeoffs.
 * 
 * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms 
 * and orders with tradeoffs. In: IACR ePrint Archive, 2018/797.
 * 
 * This function uses the default precision #PRECISION internally. There is no 
 * need to go beyond the default precision in theta_r when calling the function.
 * 
 * \param[in, out] norm   The approximate probability for the pair.
 * \param[in] theta_r     The angle theta_r.
 * \param[in] parameters  The parameters for the probability distribution.
 */
void linear_probability_r(
  mpfr_t norm,
  const mpfr_t theta_r,
  const Parameters * const parameters);

#endif /* LINEAR_PROBABILITY_H */
