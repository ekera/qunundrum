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
 * This function uses the closed-form expression in [1] for the probability of 
 * observing a pair (j, k) with angle theta_d in Ekerå's and Ekerå–Håstad's 
 * algorithms [2] for computing short discrete logarithms.
 * 
 * [1] Ekerå, M.: On post-processing in the quantum algorithm for computing 
 * short discrete logarithms. Des. Codes, Cryptogr. 88, pp. 2313–2335 (2020).
 * 
 * [2] Ekerå, M. and Håstad, J.: Quantum algorithms for computing short discrete
 * logarithms and factor RSA integers. In: PQCrypto 2017, Springer LNCS 10346,
 * pp. 347-363 (2017).
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
 * This function uses the closed-form expression in appendix A to [1] for the 
 * probability of observing an integer j with angle theta_r in Shor's [2] and 
 * Seifert's [3] order-finding algorithms.
 * 
 * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms
 * and orders with tradeoffs. J. Math. Cryptol. 15, pp. 359–407 (2021).
 * 
 * [2] Shor, P.W.: Polynomial-time algorithms for prime factorization and 
 * discrete logarithms on a quantum computer. In: SIAM Journal on Scientific 
 * Computing (SISC), volume 26(5), pp. 1484 (1997). 
 * 
 * [3] Seifert, J.-P.: Using fewer Qubits in Shor's factorization algorithm via 
 * simultaneous Diophantine approximation. In: CT-RSA 2001, Spring LNCS 2020, 
 * pp. 319-327 (2001).
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
