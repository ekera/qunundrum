/*!
 * \file    lattice_algebra.h
 * \ingroup lattice
 *
 * \brief   The declaration of functions for basic lattice basis matrix algebra.
 */

#ifndef LATTICE_ALGEBRA_H
#define LATTICE_ALGEBRA_H

#include <fplll/fplll.h>

#include <gmp.h>
#include <mpfr.h>

#include <vector>

#include <stdint.h>

using namespace fplll;

/*!
 * \brief   Solves a lattice basis matrix to the left for a coordinate vector.
 *
 * For A a (reduced) basis matrix for a lattice L and a target vector t this
 * function solves for a coordinate vector c so that c A = t. If t is in L,
 * the entries in the coordinate vector will be integers, otherwise they will
 * be quotients. An arbitrary precision floating point data type is used to
 * represent the entries to support both cases.
 *
 * \param[in, out] coordinates    The coordinate vector.
 * \param[in] target              Target vector in the lattice.
 * \param[in] A                   The basis for the lattice.
 * \param[in] n                   The number of samples used to construct the
 *                                basis matrix A. Otherwise put n + 1 is the 
 *                                dimension of A.
 * \param[in] precision           The required floating point precision.
 */
void solve_left_coordinates(
  std::vector< fplll::FP_NR<mpfr_t> > &coordinates,
  const std::vector< fplll::Z_NR<mpz_t> > &target,
  const fplll::ZZ_mat<mpz_t> &A,
  const uint32_t n,
  const uint32_t precision);

#endif /* LATTICE_ALGEBRA_H */
