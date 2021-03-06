/*!
 * \file    lattice_enumerate.h
 * \ingroup lattice_enumerate
 *
 * \brief   The declaration of functions for solving for d and r using Kannan's
 *          enumeration algorithm and lattice basis reduction techniques.
 */

/*!
 * \defgroup lattice_enumerate Enumeration solvers
 * \ingroup  lattice
 *
 * \brief    A module for functions for solving for d and r using Kannan's
 *           enumeration algorithm and lattice basis reduction techniques.
 */

#ifndef LATTICE_ENUMERATE_H
#define LATTICE_ENUMERATE_H

#include "lattice.h"

#include "parameters.h"
#include "random.h"

#include <fplll/fplll.h>

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>

/*!
 * \name Enumerating pre-reduced bases
 * \{
 */

/*!
 * \brief   Given n integers k, and a reduced basis A for a lattice L, this
 *          function attempts to recover d as the last component of an unknown
 *          vector u in L, by enumerating all vectors in L within a ball
 *          centered on a vector v that may be constructed from k and that is
 *          known to be close to u.
 *
 * For further details, see the paper.
 *
 * This function enumerates the lattice using an algorithm derived from Kannan,
 * returning status information. For each vector enumerated, it checks if the
 * vector yields d by testing against d stored in the parameters data structure.
 *
 * \param[out] status_d     A pointer to an enumeration entry in which to store
 *                          status information on the recovery of d.
 * \param[in] A             The (n + 1) x (n + 1) reduced basis matrix A for
 *                          the lattice L.
 * \param[in] G             The Gram-Schmidt (n + 1) x (n + 1) orthogonalized
 *                          basis matrix G for the matrix A.
 * \param[in] M             The (n + 1) x (n + 1) triangular matrix M of
 *                          Gram-Schmidt projection factors for the matrix A.
 * \param[in] ks            The n samples of integers k.
 * \param[in] n             The integer n.
 * \param[in] parameters    The parameters of the distribution. These
 *                          parameters in particular contain d.
 * \param[in] precision     The precision to use when enumerating.
 * \param[in] timeout       A timeout in seconds after which the enumeration
 *                          will be aborted if d has not been recovered. May be
 *                          set to zero to disable the timeout.
 */
void lattice_enumerate_reduced_basis_for_d(
  Lattice_Status_Recovery * const status_d,
  const fplll::ZZ_mat<mpz_t> &A,
  const fplll::FP_mat<mpfr_t> &G,
  const fplll::FP_mat<mpfr_t> &M,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  const uint32_t precision,
  const uint64_t timeout = 0);

/*!
 * \brief   Given n integers k, and a reduced basis A for a lattice L, this
 *          function attempts to recover r as the last component of an unknown
 *          short vector u in L, by enumerating short vectors in L.
 *
 * For further details, see the paper for details.
 *
 * This function enumerates the lattice using an algorithm derived from Kannan,
 * returning status information. For each vector enumerated, it checks if the
 * vector yields r by testing against r stored in the parameters data structure.
 *
 * \param[out] status_r     A pointer to an enumeration entry in which to store
 *                          status information on the recovery of r.
 * \param[in] A             The (n + 1) x (n + 1) reduced basis matrix A for
 *                          the lattice L.
 * \param[in] G             The Gram-Schmidt (n + 1) x (n + 1) orthogonalized
 *                          basis matrix G for the matrix A.
 * \param[in] M             The (n + 1) x (n + 1) triangular matrix M of
 *                          Gram-Schmidt projection factors for the matrix A.
 * \param[in] n             The integer n.
 * \param[in] parameters    The parameters of the distribution from which the
 *                          (j, k) pairs were sampled. These parameters in
 *                          particular contain r.
 * \param[in] precision     The precision to use when enumerating.
 * \param[in] timeout       A timeout in seconds after which the enumeration
 *                          will be aborted if r has not been recovered. May be
 *                          set to zero to disable the timeout.
 */
void lattice_enumerate_reduced_basis_for_r(
  Lattice_Status_Recovery * const status_r,
  const fplll::ZZ_mat<mpz_t> &A,
  const fplll::FP_mat<mpfr_t> &G,
  const fplll::FP_mat<mpfr_t> &M,
  const uint32_t n,
  const Parameters * const parameters,
  const uint32_t precision,
  const uint64_t timeout = 0);

/*!
 * \}
 */

/*!
 * \name Reducing and enumerating
 * \{
 */

/*!
 * \brief   Given n pairs (j, k) this function attempts to recover d by
 *          constructing, reducing and enumerating a lattice.
 *
 * This function enumerates the lattice using an algorithm derived from Kannan,
 * returning status information. For each vector enumerated, it checks if the
 * vector yields d by testing against d stored in the parameters data structure.
 *
 * \param[out] status_d   A pointer to an enumeration entry in which to store
 *                        status information on the recovery of d.
 * \param[in] js          The n samples of the j entry in the (j, k) pairs.
 * \param[in] ks          The n samples of the k entry in the (j, k) pairs.
 * \param[in] n           The integer n.
 * \param[in] parameters  The parameters of the distribution from which the
 *                        (j, k) pairs were sampled. These parameters in
 *                        particular contain d.
 * \param[in] algorithm   An enumeration entry that specifies the lattice basis
 *                        reduction algorithm, or combination of such
 *                        algorithms, to use when attempting recovery.
 * \param[in] precision   The precision to use when performing Gram-Schmidt
 *                        orthogonalization and executing Babai's algorithm.
 * \param[in] timeout     A timeout in seconds after which the enumeration will
 *                        aborted if d has not been recovered. May be set to
 *                        zero to disable the timeout.
 */
void lattice_enumerate_for_d(
  Lattice_Status_Recovery * const status_d,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const uint64_t timeout = 0);

/*!
 * \brief   Given n integers j, this function attempts to recover r by
 *          constructing, reducing and enumerating a lattice.
 *
 * This function enumerates the lattice using an algorithm derived from Kannan,
 * returning status information. For each vector enumerated, it checks if the
 * vector yields r by testing against r stored in the parameters data structure.
 *
 * \param[out] status_r   A pointer to an enumeration entry in which to store
 *                        status information on the recovery of r.
 * \param[in] js          The n samples of integers j.
 * \param[in] n           The integer n.
 * \param[in] parameters  The parameters of the distribution from which the
 *                        integers j were sampled. These parameters in
 *                        particular contain r.
 * \param[in] algorithm   An enumeration entry that specifies the lattice basis
 *                        reduction algorithm, or combination of such
 *                        algorithms, to use when attempting recovery.
 * \param[in] precision   The precision to use when performing Gram-Schmidt
 *                        orthogonalization and executing Babai's algorithm.
 * \param[in] timeout     A timeout in seconds after which the enumeration will
 *                        aborted if r has not been recovered. May be set to
 *                        zero to disable the timeout.
 */
void lattice_enumerate_for_r(
  Lattice_Status_Recovery * const status_r,
  const mpz_t * const js,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const uint64_t timeout = 0);

/*!
 * \brief   Given n pairs (j, k) this function attempts to recover d and r by
 *          constructing, reducing and enumerating a lattice.
 *
 * This function enumerates the lattice using an algorithm derived from Kannan,
 * returning status information. For each vector enumerated, it checks if the
 * vector yields d or r, respectively, by testing against d or r stored in the
 * parameters data structure.
 *
 * \param[out] status_d   A pointer to an enumeration entry in which to store
 *                        status information on the recovery of d.
 * \param[out] status_r   A pointer to an enumeration entry in which to store
 *                        status information on the recovery of r.
 * \param[in] js          The n samples of the j entry in the (j, k) pairs.
 * \param[in] ks          The n samples of the k entry in the (j, k) pairs.
 * \param[in] n           The integer n.
 * \param[in] parameters  The parameters of the distribution from which the
 *                        (j, k) pairs were sampled. These parameters in
 *                        particular contain d and r.
 * \param[in] algorithm   An enumeration entry that specifies the lattice basis
 *                        reduction algorithm, or combination of such
 *                        algorithms, to use when attempting recovery.
 * \param[in] precision   The precision to use when performing Gram-Schmidt
 *                        orthogonalization and executing Babai's algorithm.
 * \param[in] timeout     A timeout in seconds after which the enumeration will
 *                        aborted if d and r, respectively, has not been
 *                        recovered. May be set to zero to disable the timeout.
 */
void lattice_enumerate_for_d_r(
  Lattice_Status_Recovery * const status_d,
  Lattice_Status_Recovery * const status_r,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const uint64_t timeout = 0);

/*!
 * \}
 */

#endif /* LATTICE_ENUMERATE_H */
