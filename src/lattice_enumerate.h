/*!
 * \file    lattice_enumerate.h
 * \ingroup lattice_enumerate
 *
 * \brief   The declaration of functions for recovering d and r by attempting to
 *          enumerate all vectors within a ball in the lattice L.
 */

/*!
 * \defgroup lattice_enumerate Enumerating solvers
 * \ingroup  lattice
 *
 * \brief    A module for functions for recovering d and r by attempting to
 *           enumerate all vectors within a ball in the lattice L.
 */

#ifndef LATTICE_ENUMERATE_H
#define LATTICE_ENUMERATE_H

#include "lattice.h"

#include "common.h"
#include "diagonal_parameters.h"
#include "parameters.h"

#include <gmp.h>
#include <mpfr.h>

#include <fplll/fplll.h>

#include <stdint.h>

/*!
 * \name Enumerating pre-reduced bases
 * \{
 */

/*!
 * \brief   Given n integers k, and a reduced basis A for the lattice L, this
 *          function attempts to recover d by enumerating vectors in L.
 *
 * More specifically, d is the last component of an unknown vector u in L. The
 * unknown vector u is close to a known vector v that may be constructed from k.
 *
 * This function uses the reduced basis A to attempt to enumerate all vectors in
 * L in a ball centered on v, with the aim of recovering u and by extension d.
 * For further details, see [1, 2, 3].
 *
 * [1] Ekerå, M. and Håstad, J.: Quantum algorithms for computing short discrete
 * logarithms and factor RSA integers. In: PQCrypto 2017, Springer LNCS 10346,
 * pp. 347-363 (2017).
 *
 * [2] Ekerå, M.: On post-processing in the quantum algorithm for computing
 * short discrete logarithms. Des. Codes, Cryptogr. 88, pp. 2313–2335 (2020).
 *
 * [3] Ekerå, M.: Quantum algorithms for computing general discrete logarithms
 * and orders with tradeoffs. J. Math. Cryptol. 15, pp. 359–407 (2021).
 *
 * \param[out] status_d       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of d.
 * \param[in] A               The (n + 1) x (n + 1) reduced basis matrix A for
 *                            the lattice L.
 * \param[in] G               The Gram-Schmidt (n + 1) x (n + 1) orthogonalized
 *                            basis matrix G for the matrix A.
 * \param[in] M               The (n + 1) x (n + 1) triangular matrix M of
 *                            Gram-Schmidt projection factors for the matrix A.
 * \param[in] ks              The n entries for k in the (j, k) pairs.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution. These
 *                            parameters in particular contain d.
 * \param[in] precision       The precision to use when enumerating.
 * \param[in] detect_smooth_r A flag that may be set to #TRUE to detect and
 *                            handle cases where the shortest non-zero vector in
 *                            the reduced lattice basis is on the form u_r / z,
 *                            for smooth z and u_r the vector yielding r, in
 *                            which case r / z and d mod (r / z) may be
 *                            recovered instead of d and r when solving for a
 *                            general discrete logarithm d. Setting this flag to
 *                            #TRUE hence increases the probability of solving
 *                            for a general discrete logarithm d, when the group
 *                            order r is very smooth and n is close to one. This
 *                            flag has no effect for short discrete logarithms.
 *                            There is therefore no default value for this flag.
 * \param[in] timeout         A timeout in seconds after which the enumeration
 *                            will be aborted if d has not been recovered. May
 *                            be set to zero to disable the timeout.
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
  const bool detect_smooth_r,
  const uint64_t timeout = 0);

/*!
 * \brief   Given a reduced basis A for the lattice L, this function attempts to
 *          recover r by enumerating vectors in L.
 *
 * More specifically, r is the last component of an unknown short vector u in L.
 *
 * This function uses the reduced basis A to attempt to enumerate all vectors in
 * L in a ball centered on the origin, with the aim of recovering u and by
 * extension r. For further details, see [1].
 *
 * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms
 * and orders with tradeoffs. J. Math. Cryptol. 15, pp. 359–407 (2021).
 *
 * \param[out] status_r       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of r.
 * \param[in] A               The (n + 1) x (n + 1) reduced basis matrix A for
 *                            the lattice L.
 * \param[in] G               The Gram-Schmidt (n + 1) x (n + 1) orthogonalized
 *                            basis matrix G for the matrix A.
 * \param[in] M               The (n + 1) x (n + 1) triangular matrix M of
 *                            Gram-Schmidt projection factors for the matrix A.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution. These
 *                            parameters in particular contain r.
 * \param[in] precision       The precision to use when enumerating.
 * \param[in] detect_smooth_r A flag that may be set to #TRUE to detect and
 *                            handle cases where the shortest non-zero vector in
 *                            the reduced lattice basis is on the form u_r / z,
 *                            for smooth z and u_r the vector yielding r, in
 *                            which case r / z and d mod (r / z) may be
 *                            recovered instead of d and r. Setting this flag to
 *                            #TRUE hence increases the probability of solving
 *                            for r, when r is very smooth and n is close to
 *                            one. Defaults to #TRUE.
 * \param[in] timeout         A timeout in seconds after which the enumeration
 *                            will be aborted if r has not been recovered. May
 *                            be set to zero to disable the timeout.
 */
void lattice_enumerate_reduced_basis_for_r(
  Lattice_Status_Recovery * const status_r,
  const fplll::ZZ_mat<mpz_t> &A,
  const fplll::FP_mat<mpfr_t> &G,
  const fplll::FP_mat<mpfr_t> &M,
  const uint32_t n,
  const Parameters * const parameters,
  const uint32_t precision,
  const bool detect_smooth_r = TRUE,
  const uint64_t timeout = 0);

/*!
 * \brief   Given n integers k, and a reduced basis A for the lattice L, this
 *          function attempts to recover d given r by enumerating vectors in L.
 *
 * More specifically, dr is the last component of an unknown vector u in L. The
 * unknown vector u is close to a known vector v that may be constructed from k.
 *
 * This function uses the reduced basis A to attempt to enumerate all vectors in
 * L in a ball centered on v, with the aim of recovering u and by extension dr,
 * which in turn yields d since r is given. For further details, see [1].
 *
 * [1] Ekerå, M.: Revisiting Shor's quantum algorithm for computing general
 * discrete logarithms. In: ArXiv Pre-Print 1905.09084v2.
 *
 * \param[out] status_d       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of d.
 * \param[in] A               The (n + 1) x (n + 1) reduced basis matrix A for
 *                            the lattice L.
 * \param[in] G               The Gram-Schmidt (n + 1) x (n + 1) orthogonalized
 *                            basis matrix G for the matrix A.
 * \param[in] M               The (n + 1) x (n + 1) triangular matrix M of
 *                            Gram-Schmidt projection factors for the matrix A.
 * \param[in] ks              The n entries for k in the (j, k) pairs.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution. These
 *                            parameters in particular contain d and r.
 * \param[in] precision       The precision to use when enumerating.
 * \param[in] timeout         A timeout in seconds after which the enumeration
 *                            will be aborted if d has not been recovered. May
 *                            be set to zero to disable the timeout.
 */
void lattice_enumerate_reduced_basis_for_d_given_r(
  Lattice_Status_Recovery * const status_d,
  const fplll::ZZ_mat<mpz_t> &A,
  const fplll::FP_mat<mpfr_t> &G,
  const fplll::FP_mat<mpfr_t> &M,
  const mpz_t * const ks,
  const uint32_t n,
  const Diagonal_Parameters * const parameters,
  const uint32_t precision,
  const uint64_t timeout);

/*!
 * \}
 */

/*!
 * \name Reducing and enumerating
 * \{
 */

/*!
 * \brief   Given n pairs (j, k), this function attempts to recover d by
 *          constructing a basis A for the lattice L given j, reducing A, and
 *          enumerating L given A and k.
 *
 * This function calls lattice_compute_reduced_basis() to setup and reduce the
 * basis matrix A for the lattice L.
 *
 * It then calls lattice_enumerate_reduced_basis_for_d() to solve for d given A.
 *
 * For further details, see the documentation for said functions.
 *
 * \param[out] status_d       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of d.
 * \param[in] js              The n entries for j in the (j, k) pairs.
 * \param[in] ks              The n entries for k in the (j, k) pairs.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution. These
 *                            parameters in particular contain d.
 * \param[in] algorithm       An enumeration entry that specifies the lattice
 *                            basis reduction algorithm, or combination of such
 *                            algorithms, to use when attempting recovery.
 * \param[in] precision       The precision to use when performing Gram-Schmidt
 *                            orthogonalization and executing Babai's algorithm.
 * \param[in] detect_smooth_r A flag that may be set to #TRUE to detect and
 *                            handle cases where the shortest non-zero vector in
 *                            the reduced lattice basis is on the form u_r / z,
 *                            for smooth z and u_r the vector yielding r, in
 *                            which case r / z and d mod (r / z) may be
 *                            recovered instead of d and r when solving for a
 *                            general discrete logarithm d. Setting this flag to
 *                            #TRUE hence increases the probability of solving
 *                            for a general discrete logarithm d, when the group
 *                            order r is very smooth and n is close to one. The
 *                            flag has no effect for short discrete logarithms.
 *                            There is therefore no default value for the flag.
 * \param[in] timeout         A timeout in seconds after which the enumeration
 *                            will be aborted if d has not been recovered. May
 *                            be set to zero to disable the timeout.
 */
void lattice_enumerate_for_d(
  Lattice_Status_Recovery * const status_d,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const bool detect_smooth_r,
  const uint64_t timeout = 0);

/*!
 * \brief   Given n integers j, this function attempts to recover r by
 *          constructing a basis A for the lattice L given j, reducing A, and
 *          enumerating L given A.
 *
 * This function calls lattice_compute_reduced_basis() to setup and reduce the
 * basis matrix A for the lattice L.
 *
 * It then calls lattice_enumerate_reduced_basis_for_r() to solve for r given A.
 *
 * For further details, see the documentation for said functions.
 *
 * \param[out] status_r       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of r.
 * \param[in] js              The n samples of integers j.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution. These
 *                            parameters in particular contain r.
 * \param[in] algorithm       An enumeration entry that specifies the lattice
 *                            basis reduction algorithm, or combination of such
 *                            algorithms, to use when attempting recovery.
 * \param[in] precision       The precision to use when performing Gram-Schmidt
 *                            orthogonalization and executing Babai's algorithm.
 * \param[in] detect_smooth_r A flag that may be set to #TRUE to detect and
 *                            handle cases where the shortest non-zero vector in
 *                            the reduced lattice basis is on the form u_r / z,
 *                            for smooth z and u_r the vector yielding r, in
 *                            which case r / z and d mod (r / z) may be
 *                            recovered instead of d and r. Setting this flag to
 *                            #TRUE hence increases the probability of solving
 *                            for r, when r is very smooth and n is close to
 *                            one. Defaults to #TRUE.
 * \param[in] timeout         A timeout in seconds after which the enumeration
 *                            will be aborted if r has not been recovered. May
 *                            be set to zero to disable the timeout.
 */
void lattice_enumerate_for_r(
  Lattice_Status_Recovery * const status_r,
  const mpz_t * const js,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const bool detect_smooth_r = TRUE,
  const uint64_t timeout = 0);

/*!
 * \brief   Given n pairs (j, k), this function attempts to recover d and r by
 *          constructing a basis A for the lattice L given j, reducing A, and
 *          enumerating L given A given k.
 *
 * This function calls lattice_compute_reduced_basis() to setup and reduce the
 * basis matrix A for the lattice L.
 *
 * It then calls the functions lattice_enumerate_reduced_basis_for_d() and
 * lattice_enumerate_reduced_basis_for_r() to solve for d and r given A.
 *
 * For further details, see the documentation for said functions.
 *
 * \param[out] status_d       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of d.
 * \param[out] status_r       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of r.
 * \param[in] js              The n entries for j in the (j, k) pairs.
 * \param[in] ks              The n entries for k in the (j, k) pairs.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution. These
 *                            parameters in particular contain d and r.
 * \param[in] algorithm       An enumeration entry that specifies the lattice
 *                            basis reduction algorithm, or combination of such
 *                            algorithms, to use when attempting recovery.
 * \param[in] precision       The precision to use when performing Gram-Schmidt
 *                            orthogonalization and executing Babai's algorithm.
 * \param[in] detect_smooth_r A flag that may be set to #TRUE to detect and
 *                            handle cases where the shortest non-zero vector in
 *                            the reduced lattice basis is on the form u_r / z,
 *                            for smooth z and u_r the vector yielding r, in
 *                            which case r / z and d mod (r / z) may be
 *                            recovered instead of d and r. Setting this flag to
 *                            #TRUE hence increases the probability of solving
 *                            for a general discrete logarithm d and order r,
 *                            when r is very smooth and n is close to one.
 *                            Defaults to #TRUE.
 * \param[in] timeout         A timeout in seconds after which the enumeration
 *                            will be aborted if d and r, respectively, has not
 *                            been recovered. May be set to zero to disable the
 *                            timeout.
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
  const bool detect_smooth_r = TRUE,
  const uint64_t timeout = 0);

/*!
 * \brief   Given n pairs (j, k), this function attempts to recover d by
 *          constructing a basis A for the lattice L given j and r, reducing A,
 *          and enumerating L given A and k.
 *
 * This function calls lattice_compute_reduced_diagonal_basis() to setup and
 * reduce the basis matrix A for the lattice L.
 *
 * It then calls lattice_enumerate_reduced_basis_for_d_given_r() to solve for d
 * given r and A.
 *
 * For further details, see the documentation for said functions.
 *
 * \param[out] status_d       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of d.
 * \param[in] js              The n entries for j in the (j, k) pairs.
 * \param[in] ks              The n entries for k in the (j, k) pairs.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution. These
 *                            parameters in particular contain d and r.
 * \param[in] algorithm       An enumeration entry that specifies the lattice
 *                            basis reduction algorithm, or combination of such
 *                            algorithms, to use when attempting recovery.
 * \param[in] precision       The precision to use when performing Gram-Schmidt
 *                            orthogonalization and executing Babai's algorithm.
 * \param[in] timeout         A timeout in seconds after which the enumeration
 *                            will be aborted if d has not been recovered. May
 *                            be set to zero to disable the timeout.
 */
void lattice_enumerate_for_d_given_r(
  Lattice_Status_Recovery * const status_d,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Diagonal_Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const uint64_t timeout);

/*!
 * \}
 */

#endif /* LATTICE_ENUMERATE_H */
