/*!
 * \file    lattice_smoothness.h
 * \ingroup lattice_solve
 *
 * \brief   The declaration of functions for checking if an integer is smooth 
 *          and for removing smooth factors from partially smooth integers.
 */

#ifndef LATTICE_SMOOTHNESS_H
#define LATTICE_SMOOTHNESS_H

#include <gmp.h>

#include <stdint.h>

/*!
 * \brief   Removes all prime factor powers q^e < 2^m from z, where it is 
 *          furthermore required that q <= cm, for c and m constants.
 * 
 * \param[in, out] reduced_z  The integer z with smooth factors removed.
 * \param[in] z               The integer z.
 * \param[in] c               The constant c.
 * \param[in] m               The constant m.
 */
void lattice_smoothness_remove_smooth_factors(
  mpz_t reduced_z,
  const mpz_t z,
  const double c,
  const uint32_t m);

/*!
 * \brief   Checks if z is smooth in the sense that it is a product of prime 
 *          factor powers q^e < 2^m from z, where it is furthermore required 
 *          that q <= cm, for c and m constants.
 * 
 * \param[in] z   The integer z.
 * \param[in] c   The constant c.
 * \param[in] m   The constant m.
 * 
 * \return  Returns #TRUE if z is smooth, #FALSE otherwise.
 */
bool lattice_smoothness_is_smooth(
  const mpz_t z,
  const double c,
  const uint32_t m);

#endif /* LATTICE_SMOOTHNESS_H */
