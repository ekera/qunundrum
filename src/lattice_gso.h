/*!
 * \file    lattice_gso.h
 * \ingroup lattice
 *
 * \brief   The declaration of functions for computing the Gram-Schmidt
 *          orthogonalized basis of a lattice basis matrix.
 */

#ifndef LATTICE_GSO_H
#define LATTICE_GSO_H

#include <gmp.h>
#include <mpfr.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wtemplate-id-cdtor"
#endif
#include <fplll/fplll.h>
#pragma GCC diagnostic pop

#include <stdint.h>

/*!
 * \brief Computes a Gram-Schmidt orthogonalized basis matrix G for a full rank
 *        (n + 1) x (n + 1) basis matrix A and a triangular matrix M of
 *        Gram-Schmidt projection factors such that A = M * G.
 *
 * \param[in, out] M          The (n + 1) x (n + 1) triangular matrix M of
 *                            Gram-Schmidt projection factors.
 * \param[in, out] G          The Gram-Schmidt (n + 1) x (n + 1) orthogonalized
 *                            basis matrix G.
 * \param[in] A               A full rank (n + 1) x (n + 1) basis matrix.
 * \param[in] n               The integer n.
 * \param[in] precision       The required floating point precision.
 */
void gram_schmidt_orthogonalization(
  fplll::FP_mat<mpfr_t> &M,
  fplll::FP_mat<mpfr_t> &G,
  const fplll::ZZ_mat<mpz_t> &A,
  const uint32_t n,
  const uint32_t precision);

#endif /* LATTICE_GSO_H */
