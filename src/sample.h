/*!
 * \file    sample.h
 * \ingroup sample_distribution
 *
 * \brief   The declaration of functions for sampling probability distributions.
 */

/*!
 * \defgroup sample_distribution Sampling distributions
 * \ingroup  distribution
 *
 * \brief    A module for functions for sampling probability distributions.
 */

#ifndef SAMPLE_H
#define SAMPLE_H

#include "parameters.h"
#include "diagonal_parameters.h"
#include "random.h"

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>

/*!
 * \name Sampling alpha
 * \{
 */

/*!
 * \brief   Coarsely samples alpha uniformly at random from a given interval.
 *
 * Alpha is coarsely sampled, and is not guaranteed to be admissible.
 *
 * \param[in, out] alpha          The argument alpha sampled.
 * \param[in] min_log_alpha       The minimum signed logarithmic alpha.
 * \param[in] max_log_alpha       The maximum signed logarithmic alpha.
 * \param[in, out] random_state   The random state to use when sampling.
 */
void sample_approximate_alpha_from_region(
  mpfr_t alpha,
  const double min_log_alpha,
  const double max_log_alpha,
  Random_State * const random_state);

/*!
 * \brief   Samples alpha uniformly at random from an interval, with the
 *          restriction that 2^kappa must divide alpha.
 *
 * Alpha is sampled with high precision, uniformly at random from the interval
 *
 *  [2^abs(min_log_alpha), 2^abs(max_log_alpha)),
 * 
 * after which alpha is reduced by (alpha % 2^kappa) to ensure that that 
 * alpha is admissible. The sign is then flipped if min_log_alpha is negative.
 *
 * \param[in, out] alpha          The argument alpha sampled.
 * \param[in] min_log_alpha       The minimum signed logarithmic alpha.
 * \param[in] max_log_alpha       The maximum signed logarithmic alpha.
 * \param[in] kappa               The integer kappa.
 * \param[in, out] random_state   The random state to use when sampling.
 */
void sample_alpha_from_region(
  mpz_t alpha,
  const double min_log_alpha,
  const double max_log_alpha,
  const uint32_t kappa,
  Random_State * const random_state);

/*!
 * \}
 */

/*!
 * \name Sampling outputs
 * \{
 */

/*!
 * \brief   Samples an integer j uniformly at random from the set of all
 *          integers that yield the admissible argument alpha_r.
 *
 * This function assumes alpha_r to be from a linear distribution computed
 * by collapsing a two-dimensional distribution to a marginal distribution, or
 * computed for an order.
 *
 * \param[in, out] j              The integer j.
 * \param[in] alpha_r             The argument alpha_r.
 * \param[in] parameters          The probability distribution parameters.
 * \param[in, out] random_state   The random state to use when sampling.
 */
void sample_j_from_alpha_r(
  mpz_t j,
  const mpz_t alpha_r,
  const Parameters * const parameters,
  Random_State * const random_state);

/*!
 * \brief   Samples an integer pair (j, k) uniformly at random from the set of
 *          all pairs that yield the admissible argument alpha_d.
 *
 * This function assumes alpha_d to be from a linear distribution computed
 * by collapsing a two-dimensional distribution to a marginal distribution, or
 * computed for a short discrete logarithm.
 *
 * \param[in, out] j              The integer j.
 * \param[in, out] k              The integer k.
 * \param[in] alpha_d             The argument alpha_d.
 * \param[in] parameters          The probability distribution parameters.
 * \param[in, out] random_state   The random state to use when sampling.
 */
void sample_j_k_from_alpha_d(
  mpz_t j,
  mpz_t k,
  const mpz_t alpha_d,
  const Parameters * const parameters,
  Random_State * const random_state);

/*!
 * \brief   Samples an integer pair (j, k) uniformly at random from the set of
 *          all pairs that yield the admissible argument pair 
 *          (alpha_d, alpha_r).
 *
 * This function assumes (alpha_d, alpha_r) to be from a two-dimensional
 * distribution. The argument pair must be admissible. You may use lattice-based
 * techniques to map high resolution independent estimate of alpha_d and 
 * alpha_r, respectively, to an admissible pair (alpha_d, alpha_r), from which 
 * (j, k) may then be sampled using this function. Calling this function with 
 * inadmissible argument pairs may produce unexpected behavior.
 *
 * \param[in, out] j              The integer j.
 * \param[in, out] k              The integer k.
 * \param[in] alpha_d             The argument alpha_d.
 * \param[in] alpha_r             The argument alpha_r.
 * \param[in] parameters          The probability distribution parameters.
 * \param[in, out] random_state   The random state to use when sampling.
 */
void sample_j_k_from_alpha_d_r(
  mpz_t j,
  mpz_t k,
  const mpz_t alpha_d,
  const mpz_t alpha_r,
  const Parameters * const parameters,
  Random_State * const random_state);

/*!
 * \brief   Samples an integer j uniformly at random from the set of all 
 *          integers j that yield the admissible argument alpha_r.
 *
 * This function assumes alpha_r to be from a diagonal distribution.
 *
 * \param[in, out] j              The integer j.
 * \param[in] alpha_r             The argument alpha_r.
 * \param[in] parameters          The probability distribution parameters.
 * \param[in, out] random_state   The random state to use when sampling.
 */
void sample_j_from_diagonal_alpha_r(
  mpz_t j,
  const mpz_t alpha_r,
  const Diagonal_Parameters * const parameters,
  Random_State * const random_state);

/*!
 * \}
 */

#endif /* SAMPLE_H */
