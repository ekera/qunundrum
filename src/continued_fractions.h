/*!
 * \file    continued_fractions.h
 * \ingroup continued_fractions
 *
 * \brief   The declaration of functions for computing continued fractions.
 *
 * These functions are used to implement Shor's solver for order finding.
 */

/*!
 * \defgroup continued_fractions Continued fractions
 * \ingroup math
 *
 * \brief A module for functions for computing continued fractions.
 */

#ifndef CONTINUED_FRACTIONS_H
#define CONTINUED_FRACTIONS_H

#include "parameters.h"

#include <mpfr.h>
#include <gmp.h>

#include <stdint.h>

/*!
 * \brief   Solves for an order using Shor's original post-processing algorithm
 *          based on computing a continued fractions expansion.
 *
 * The post-processing algorithm is described in the paper [1] by Shor. This 
 * function uses the optimization proposed by Andrew Odlyzko to Shor and 
 * described in [1]. It furthermore features an optimized implementation of the
 * continued fractions expansion algorithm that only computes denominators and 
 * does so efficiently.
 * 
 * [1] Shor, P.W.: Polynomial-time algorithms for prime factorization and 
 * discrete logarithms on a quantum computer. In: SIAM Journal on Scientific 
 * Computing (SISC), volume 26(5), pp. 1484 (1997). 
 * 
 * \param[in] j                      The integer j.
 * \param[in] search_bound_cofactor  A search bound on the cofactor. Recall that
 *                                   it may be the case that the denominator  
 *                                   obtained and r share factors, in which case
 *                                   we need to exhaustively search cofactors.
 * \param[in] parameters             The parameters for the probability 
 *                                   distribution from which the integer j was
 *                                   sampled. Contains r and m.
 * 
 * \return  Returns #TRUE if the order could be recovered, #FALSE otherwise.
 */
bool continued_fractions_solve(
  const mpz_t j,
  const uint32_t search_bound_cofactor,
  const Parameters * const parameters);

#endif /* CONTINUED_FRACTIONS_H */
