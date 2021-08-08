/*!
 * \file    tau_estimate.h
 * \ingroup estimating_tau
 *
 * \brief   The declaration of functions for computing tau estimates that may
 *          be converted into volume quotient estimates.
 */

/*!
 * \defgroup estimating_tau Estimating tau
 * \ingroup  estimating_volume_quotients
 *
 * \brief    A module for functions for estimating tau.
 */

#ifndef TAU_ESTIMATE_H
#define TAU_ESTIMATE_H

#include "distribution.h"
#include "linear_distribution.h"
#include "parameters.h"
#include "random.h"

#include <stdint.h>

/*!
 * \brief   Computes tau_d and tau_r from n argument pairs (alpha_d, alpha_r)
 *          sampled from the distribution.
 * 
 * If at least one of the n argument pairs sampled is out of bounds, tau_d
 * and tau_r are both set to MAX_DBL and #FALSE is returned. Otherwise,
 * tau_d and tau_r are valid estimates and #TRUE is returned.
 *
 * For x in {d, r}, we define tau_x as
 * 
 *  tau_x = log_( (1 / n) sum_{i = 1}^{n} alpha^2_{x, i} ) / 2 - m.
 * 
 * To compute the volume quotients v_d and v_r, you may convert the output from
 * this function into volume quotients using tau_volume_quotient().
 *
 * \param[in] distribution      The distribution from which to sample pairs.
 * \param[in, out] random_state The random state to use when sampling.
 * \param[in] n                 The integer n.
 * \param[out] tau_d            The computed value of tau_d.
 * \param[out] tau_r            The computed value of tau_r.
 *
 * \return  Returns #FALSE if at least one of the n argument pairs sampled are
 *          out of bounds, #TRUE otherwise.
 */
bool tau_estimate(
  const Distribution * const distribution,
  Random_State * const random_state,
  const uint32_t n,
  long double &tau_d,
  long double &tau_r);

/*!
 * \brief   Computes tau from n arguments alpha sampled from the distribution.
 *
 * If at least one of the n arguments alpha sampled is out of bounds, tau
 * is set to MAX_DBL and #FALSE is returned. Otherwise, tau is a valid estimate
 * and #TRUE is returned.
 *
 * We define tau = log_( (1 / n) sum_{i = 1}^{n} alpha^2_i ) / 2 - m.
 * 
 * To compute the volume quotient v, you may convert the output from this
 * function into a volume quotient using tau_volume_quotient().
 *
 * \param[in] distribution      The distribution from which to sample alpha.
 * \param[in, out] random_state The random state to use when sampling.
 * \param[in] n                 The integer n.
 * \param[out] tau              The computed value of tau.
 *
 * \return  Returns #FALSE if at least one of the n arguments alpha sampled is
 *          out of bounds, #TRUE otherwise.
 */
bool tau_estimate_linear(
  const Linear_Distribution * const distribution,
  Random_State * const random_state,
  const uint32_t n,
  long double &tau);

#endif /* TAU_ESTIMATE_H */
