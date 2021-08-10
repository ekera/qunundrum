/*!
 * \file    lattice_solve.h
 * \ingroup lattice_solve
 *
 * \brief   The declaration of functions for solving for d and r using Babai's
 *          nearest plane algorithm and lattice basis reduction techniques.
 */

/*!
 * \defgroup lattice_solve Nearest plane solvers
 * \ingroup  lattice
 *
 * \brief    A module for functions for solving for d and r using Babai's
 *           nearest plane algorithm and lattice basis reduction techniques.
 */

#ifndef LATTICE_SOLVE_H
#define LATTICE_SOLVE_H

#include "lattice.h"

#include "parameters.h"
#include "diagonal_parameters.h"
#include "random.h"

#include <fplll/fplll.h>

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>

/*!
 * \brief The constant c used to setup the smoothness bound for the order r.
 */
#define LATTICE_SMOOTHNESS_CONSTANT_C 1

/*!
 * \name Solving pre-reduced bases
 * \{
 */

/*!
 * \brief   Given n integers k, and a reduced basis A for a lattice L, this
 *          function attempts to recover d by solving a closest vector problem 
 *          in L.
 * 
 * More specifically, d is the last component of an unknown vector u in L. The
 * unknown vector u is close to a known vector v that may be constructed from k.
 * This function finds the closest vector to v in L, with the aim of recovering 
 * u and by extension d. For further details, see the paper.
 *
 * \param[out] status_d       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of d.
 * \param[in] A               The (n + 1) x (n + 1) reduced basis matrix A for
 *                            the lattice L.
 * \param[in] G               The Gram-Schmidt (n + 1) x (n + 1) orthogonalized
 *                            basis matrix G for the matrix A.
 * \param[in] ks              The n entries for k in the (j, k) pairs.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution. These
 *                            parameters in particular contain d.
 * \param[in] precision       The precision to use when performing Gram-Schmidt
 *                            orthogonalization and executing Babai's algorithm.
 * \param[in] detect_smooth_r A flag that may be set to #TRUE to detect if the 
 *                            shortest vector in the lattice is on the form
 *                            u_r / z for smooth z, in which case d mod r / z 
 *                            may be recovered instead of d. This increases the
 *                            probability of solving for d when computing
 *                            general discrete logarithms, when r is very smooth
 *                            and n is close to one.
 */
void lattice_solve_reduced_basis_for_d(
  Lattice_Status_Recovery * const status_d,
  const fplll::ZZ_mat<mpz_t> &A,
  const fplll::FP_mat<mpfr_t> &G,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  const uint32_t precision,
  const bool detect_smooth_r);

/*!
 * \brief   Given a reduced basis A for a lattice L, this function attempts to 
 *          recover r by solving a shortest non-zero vector problem in L.
 *
 * More specifically, r is the last component of an unknown short vector u in L.
 * This function finds the shortest non-zero vector in L,  with the aim of 
 * recovering u and by extension r. For further details, see the paper.
 *
 * \param[out] status_r       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of r.
 * \param[in] A               The (n + 1) x (n + 1) reduced basis matrix A for
 *                            the lattice L.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution. These
 *                            parameters in particular contain r.
 * \param[in] detect_smooth_r A flag that may be set to #TRUE to detect if the
 *                            shortest vector in the lattice is on the form
 *                            u_r / z for smooth z, in which case r / z may be
 *                            recovered instead of r. This increases the
 *                            probability of solving for r, when r is very
 *                            smooth and n is close to one. Defaults to #TRUE.
 */
void lattice_solve_reduced_basis_for_r(
  Lattice_Status_Recovery * const status_r,
  const fplll::ZZ_mat<mpz_t> &A,
  const uint32_t n,
  const Parameters * const parameters,
  const bool detect_smooth_r = TRUE);

/*!
 * \brief   Given r, n integers k, and a reduced basis A for the lattice L, this 
 *          function attempts to recover d by solving a closest vector problem 
 *          in L.
 * 
 * More specifically, dr is the last component of an unknown vector u in L. The
 * unknown vector u is close to a known vector v that may be constructed from k.
 * This function finds the closest vector to v in L, with the aim of recovering 
 * u and by extension dr, where r is given. For further details, see the paper.
 */
void lattice_solve_reduced_basis_for_d_given_r(
  Lattice_Status_Recovery * const status_d,
  const fplll::ZZ_mat<mpz_t> &A,
  const fplll::FP_mat<mpfr_t> &G,
  const mpz_t * const ks,
  const uint32_t n,
  const Diagonal_Parameters * const parameters,
  const uint32_t precision);

/*!
 * \}
 */

/*!
 * \name Reducing and solving
 * \{
 */

/*!
 * \brief   Given n pairs (j, k), this function attempts to recover d by
 *          using j to construct a basis for a lattice L, reducing the basis, 
 *          and solving a closest vector problem in L.
 *
 * This function calls lattice_compute_reduced_basis() to setup and reduce the 
 * basis. It then calls lattice_solve_reduced_basis_for_d() to solve for d. For
 * further details, see the documentation for said functions, and the paper.
 *
 * \param[out] status_d       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of d.
 * \param[in] js              The n entries for j in the (j, k) pairs.
 * \param[in] ks              The n entries for k in the (j, k) pairs.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution from which the
 *                            (j, k) pairs were sampled. These parameters in
 *                            particular contain d.
 * \param[in] algorithm       An enumeration entry that specifies the lattice
 *                            basis reduction algorithm, or combination of such
 *                            algorithms, to use when attempting recovery.
 * \param[in] precision       The precision to use when performing Gram-Schmidt
 *                            orthogonalization and executing Babai's algorithm.
 * \param[in] detect_smooth_r A flag that may be set to #TRUE to detect if the
 *                            shortest vector in the lattice is on the form
 *                            u_r / z for smooth z, in which case d mod r / z
 *                            may be recovered instead of d. This increases the
 *                            probability of solving for d when computing
 *                            general discrete logarithms, when r is very smooth
 *                            and n is close to one. Set to #FALSE for short 
 *                            discrete logarithms. Defaults to #FALSE.
 */
void lattice_solve_for_d(
  Lattice_Status_Recovery * const status_d,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const bool detect_smooth_r = FALSE);

/*!
 * \brief   Given n integers j, this function attempts to recover d by using j 
 *          to construct a basis for a lattice L, reducing the basis, and 
 *          solving a shortest non-zero vector problem in L.
 * 
 * This function calls lattice_compute_reduced_basis() to setup and reduce the 
 * basis. It then calls lattice_solve_reduced_basis_for_r() to solve for r. For
 * further details, see the documentation for said functions, and the paper.
 * 
 * \param[out] status_r       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of r.
 * \param[in] js              The n samples of integers j.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution from which the
 *                            integers j were sampled. These parameters in
 *                            particular contain r.
 * \param[in] algorithm       An enumeration entry that specifies the lattice
 *                            basis reduction algorithm, or combination of such
 *                            algorithms, to use when attempting recovery.
 * \param[in] detect_smooth_r A flag that may be set to #TRUE to detect if the
 *                            shortest vector in the lattice is on the form
 *                            u_r / z for smooth z, in which case r / z may be
 *                            recovered instead of r. This increases the
 *                            probability of solving for r, when r is very
 *                            smooth and n is close to one. Defaults to #TRUE.
 */
void lattice_solve_for_r(
  Lattice_Status_Recovery * const status_r,
  const mpz_t * const js,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const bool detect_smooth_r = TRUE);

/*!
 * \brief   Given n pairs (j, k), this function attempts to recover d and r by
 *          using j to construct a basis for a lattice L, reducing the basis, 
 *          and solving both a shortest non-zero vector and a closest vector 
 *          problem in L.
 *
 * This function calls lattice_compute_reduced_basis() to setup and reduce the 
 * basis. It then calls the functions lattice_solve_reduced_basis_for_d() and
 * lattice_solve_reduced_basis_for_r() to solve for d and r, respectively. For
 * further details, see the documentation for said functions, and the paper.
 *
 * \param[out] status_d       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of d.
 * \param[out] status_r       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of r.
 * \param[in] js              The n entries for j in the (j, k) pairs.
 * \param[in] ks              The n entries for k in the (j, k) pairs.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution from which the
 *                            (j, k) pairs were sampled. These parameters in
 *                            particular contain d and r.
 * \param[in] algorithm       An enumeration entry that specifies the lattice
 *                            basis reduction algorithm, or combination of such
 *                            algorithms, to use when attempting recovery.
 * \param[in] precision       The precision to use when performing Gram-Schmidt
 *                            orthogonalization and executing Babai's algorithm.
 * \param[in] detect_smooth_r A flag that may be set to #TRUE to detect if the
 *                            shortest vector in the lattice is on the form
 *                            u_r / z for smooth z, in which case d mod r / z
 *                            and r / z may be recovered instead of d and r.
 *                            This increases the probability of solving for d
 *                            and r when computing general discrete logarithms,
 *                            when r is very smooth and n is close to one.
 *                            Defaults to #TRUE.
 */
void lattice_solve_for_d_r(
  Lattice_Status_Recovery * const status_d,
  Lattice_Status_Recovery * const status_r,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const bool detect_smooth_r = TRUE);

/*!
 * \brief   Given r, and n pairs (j, k), this function attempts to recover d by
 *          using r and j to construct a basis for a lattice L, reducing the 
 *          basis, and solving a closest vector problem in L.
 *
 * This function calls lattice_compute_reduced_diagonal_basis() to setup and 
 * reduce the basis. It then calls lattice_solve_reduced_basis_for_d_given_r() 
 * to solve for d given r. For further details, see the documentation for said 
 * functions, and the paper.
 *
 * \param[out] status_d       A pointer to an enumeration entry in which to
 *                            store status information on the recovery of d.
 * \param[in] js              The n entries for j in the (j, k) pairs.
 * \param[in] ks              The n entries for k in the (j, k) pairs.
 * \param[in] n               The integer n.
 * \param[in] parameters      The parameters of the distribution from which the
 *                            (j, k) pairs were sampled. These parameters in
 *                            particular contain d and r.
 * \param[in] algorithm       An enumeration entry that specifies the lattice
 *                            basis reduction algorithm, or combination of such
 *                            algorithms, to use when attempting recovery.
 * \param[in] precision       The precision to use when performing Gram-Schmidt
 *                            orthogonalization and executing Babai's algorithm.
 */
void lattice_solve_for_d_given_r(
  Lattice_Status_Recovery * const status_d,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Diagonal_Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision);

/*!
 * \}
 */

#endif /* LATTICE_SOLVE_H */
