/*!
 * \file    lattice_gso.cpp
 * \ingroup lattice
 *
 * \brief   The definition of functions for computing the Gram-Schmidt
 *          orthogonalized basis of a lattice basis matrix.
 */

#include "lattice_gso.h"

#include "lattice.h"
#include "errors.h"

#include <fplll/fplll.h>

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>

using namespace fplll;

/*!
 * \brief Computes a Gram-Schmidt projection factor by taking dot products
 *        between rows in the (n + 1) x (n + 1) basis matrix A and the
 *        (n + 1) x (n + 1) Gram-Schmidt orthogonalized basis matrix G.
 *
 * \param[in, out] factor     The Gram-Schmidt projection factor.
 * \param[in] A               The full rank (n + 1) x (n + 1) basis matrix.
 * \param[in] row_A           The index of the row in the A matrix.
 * \param[in] G               The Gram-Schmidt (n + 1) x (n + 1) orthogonalized
 *                            basis matrix G.
 * \param[in] row_G           The index of the row in the G matrix.
 * \param[in] n               The integer n.
 * \param[in] precision       The required floating point precision.
 */
static void gram_schmidt_projection(
  mpfr_t factor,
  const ZZ_mat<mpz_t> &A,
  const uint32_t row_A,
  const FP_mat<mpfr_t> &G,
  const uint32_t row_G,
  const uint32_t n,
  const uint32_t precision)
{
  /* Sanity checks of the dimensions of the A and G matrices. */
  if ((A.get_rows() != (int)(n + 1)) || (A.get_cols() != (int)(n + 1))) {
    critical("gram_schmidt_projection(): Matrix A has incorrect dimensions.");
  }

  if ((G.get_rows() != (int)(n + 1)) || (G.get_cols() != (int)(n + 1))) {
    critical("gram_schmidt_projection(): Matrix G has incorrect dimensions.");
  }

  /* Set the precision of the output variable. */
  mpfr_set_prec(factor, precision);

  /* Declare temporary variables. */
  mpfr_t x;
  mpfr_init2(x, precision);

  mpfr_t sum;
  mpfr_init2(sum, precision);

  /* Compute the projection factors. */
  mpfr_set_ui(sum, 0, MPFR_RNDN);

  for (uint32_t i = 0; i <= n; i++) {
    mpfr_mul_z(x, G[row_G][i].get_data(), A[row_A][i].get_data(), MPFR_RNDN);
    mpfr_add(sum, sum, x, MPFR_RNDN);
  }

  mpfr_set(factor, sum, MPFR_RNDN); /* factor = <u, v> */

  mpfr_set_ui(sum, 0, MPFR_RNDN);

  for (uint32_t i = 0; i <= n; i++) {
    mpfr_mul(x, G[row_G][i].get_data(), G[row_G][i].get_data(), MPFR_RNDN);
    mpfr_add(sum, sum, x, MPFR_RNDN);
  }

  mpfr_div(factor, factor, sum, MPFR_RNDN);  /* factor = <u, v> / <u, u> */

  /* Clear memory. */
  mpfr_clear(x);
  mpfr_clear(sum);
}

void gram_schmidt_orthogonalization(
  FP_mat<mpfr_t> &M,
  FP_mat<mpfr_t> &G,
  const ZZ_mat<mpz_t> &A,
  const uint32_t n,
  const uint32_t precision)
{
  /* Resize the M and G matrices if they have incorrect dimensions. */
  if ((M.get_rows() != (int)(n + 1)) || (M.get_cols() != (int)(n + 1))) {
    M.resize(n + 1, n + 1);
  }

  if ((G.get_rows() != (int)(n + 1)) || (G.get_cols() != (int)(n + 1))) {
    G.resize(n + 1, n + 1);
  }

  /* Set the precision of the elements of M and G. */
  for (uint32_t i = 0; i <= n; i++) {
    for (uint32_t j = 0; j <= n; j++) {
      mpfr_set_prec(G[i][j].get_data(), precision);
      mpfr_set_prec(M[i][j].get_data(), precision);
    }
  }

  /* Setup the projection factors. */
  mpz_t tmp_z;
  mpz_init(tmp_z);

  mpfr_t tmp_f;
  mpfr_init2(tmp_f, precision);

  /* Ones on the top left to bottom right diagonal. Zeroes above diagonal. */
  for (uint32_t i = 0; i <= n; i++) {
    mpfr_set_ui(M[i][i].get_data(), 1, MPFR_RNDN);
    for (uint32_t j = i + 1; j <= n; j++) {
      mpfr_set_ui(M[i][j].get_data(), 0, MPFR_RNDN);
    }
  }

  /* Compute the projection factors in M and the orthogonal matrix G. */
  for (uint32_t i = 0; i <= n; i++) {
    for (uint32_t j = 0; j <= n; j++) {
        /* Copy row i from the A matrix into the G matrix. */
        mpfr_set_z(G[i][j].get_data(), A[i][j].get_data(), MPFR_RNDN);
    }

    /* Subtract previous rows from the A matrix, weighted by projection
     * factors as computed and inserted into the M matrix. */

    /* Subtract projections with previous elements. */
    for (uint32_t k = 0; k < i; k++) {
      /* Compute the projection factor and insert it into the M matrix. */
      gram_schmidt_projection(M[i][k].get_data(), A, i, G, k, n, precision);

      for (uint32_t j = 0; j <= n; j++) {
          /* Copy row k from the A matrix multiplied by M[i, k]. */
          mpfr_mul(tmp_f, M[i][k].get_data(), G[k][j].get_data(), MPFR_RNDN);
          mpfr_sub(G[i][j].get_data(), G[i][j].get_data(), tmp_f, MPFR_RNDN);
      }
    }
  }

  /* Clear memory. */
  mpz_clear(tmp_z);
  mpfr_clear(tmp_f);
}
