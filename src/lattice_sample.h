/*!
 * \file    lattice_sample.h
 * \ingroup lattice_sample
 *
 * \brief   The declaration of functions and data structures for mapping any
 *          argument pair to the closest admissible argument pair in the
 *          argument plane.
 */

/*!
 * \defgroup lattice_sample Lattice-based sampling
 * \ingroup  lattice
 * \ingroup  sample_distribution
 *
 * \brief    A module for using lattice-based techniques to sample 
 *           two-dimensional probability distributions.
 */

#ifndef LATTICE_SAMPLE_H
#define LATTICE_SAMPLE_H

#include "parameters.h"

#include <fplll/fplll.h>

#include <mpfr.h>
#include <gmp.h>

/*!
 * \brief   A data structure representing the alpha lattice.
 *
 * This lattice is used to map any argument pair to the closest admissible
 * argument pair in the argument plane. This is required when sampling argument
 * pairs and integer pairs (j, k) from two-dimensional probability
 * distributions. As the lattice depends only on the parameter used to compute
 * the probability distribution, and not on the actual argument pair that is to
 * be mapped, the lattice may be pre-computed. This is beneficial as we may
 * avoid computing it over and over again in each sampling operation.
 */
typedef struct {
  /*!
   * \brief   A reduced basis for the lattice.
   */
  fplll::ZZ_mat<mpz_t> A;

  /*!
   * \brief   The Gram-Schmidt-orthogonalized matrix for the basis A.
   */
  fplll::FP_mat<mpfr_t> G;
} Lattice_Alpha;

/*!
 * \name Initialization
 * \{
 */

/*!
 * \brief   Initializes an alpha lattice.
 *
 * \param[in, out] lattice    The alpha lattice to initialize.
 * \param[in] parameters      The distribution parameters to use to initialized
 *                            the alpha lattice.
 */
void lattice_alpha_init(
  Lattice_Alpha * const lattice,
  const Parameters * const parameters);

/*!
 * \brief   Clears an alpha lattice.
 *
 * \param[in, out] lattice    The alpha lattice to clear.
 */
void lattice_alpha_clear(
  Lattice_Alpha * const lattice);

/*!
 * \}
 */

/*!
 * \name Mapping
 * \{
 */

/*!
 * \brief   Maps any argument pair to the closest admissible argument pair in
 *          the argument plane.
 *
 * \param[in, out] alpha_d    The argument alpha_d.
 * \param[in, out] alpha_r    The argument alpha_r.
 * \param[in] lattice         The alpha lattice for the distribution parameters.
 * \param[in] parameters      The distribution parameters.
 */
void lattice_alpha_map(
  mpz_t alpha_d,
  mpz_t alpha_r,
  const Lattice_Alpha * const lattice,
  const Parameters * const parameters);

/*!
 * \}
 */

#endif /* LATTICE_SAMPLE_H */
