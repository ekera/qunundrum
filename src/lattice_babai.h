/*!
 * \file    lattice_babai.h
 * \ingroup lattice
 *
 * \brief   The declaration of functions for implementing Babai's nearest plane
 *          algorithm that may be used to find the closest vector to a target
 *          vector in an integer lattice.
 */

#ifndef LATTICE_BABAI_H
#define LATTICE_BABAI_H

#include <fplll/fplll.h>

#include <gmp.h>
#include <mpfr.h>

#include <vector>

#include <stdint.h>

/*!
 * \brief Applies Babai's nearest plane algorithm to find the closest vector to
 *        a target vector in an integer lattice defined by a basis matrix A.
 *
 * The basis matrix must be LLL- or BKZ-reduced to form a nearly orthogonal
 * basis. The caller must provide an orthogonal basis matrix G such that
 * A = G * M where M is triangular matrix of projection factors as computed
 * by calling gram_schmidt().
 *
 * \param[in, out] solution   A vector in L that is close to the target vector.
 * \param[in] target          The target vector in (n + 1)-dimensional Z-space.
 * \param[in] G               An Gram-Schmidt orthogonalized basis matrix G
 *                            such that A = G * M, where M is a triangular
 *                            matrix of Gram-Schmidt projection factors, see
 *                            the gram_schmidt_orthogonalization() function.
 * \param[in] A               An  LLL- or BKZ-reduced full rank basis matrix
 *                            of dimension (n + 1) x (n + 1) defining the
 *                            (n + 1) x (n + 1) dimensional integer lattice L.
 * \param[in] n               The integer n.
 * \param[in] precision       The required floating point precision.
 */
void babai_closest_vector(
  std::vector< fplll::Z_NR<mpz_t> > &solution,
  const std::vector< fplll::Z_NR<mpz_t> > &target,
  const fplll::FP_mat<mpfr_t> &G,
  const fplll::ZZ_mat<mpz_t> &A,
  const uint32_t n,
  const uint32_t precision);

#endif /* LATTICE_BABAI_H */
