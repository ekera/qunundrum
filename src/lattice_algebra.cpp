/*!
 * \file    lattice_algebra.cpp
 * \ingroup lattice
 *
 * \brief   The definition of functions for basic lattice basis matrix algebra.
 */

#include "lattice_algebra.h"

#include "lattice.h"
#include "errors.h"

#include <fplll/fplll.h>

#include <gmp.h>
#include <mpfr.h>

#include <vector>

#include <stdint.h>

using namespace std;
using namespace fplll;

void solve_left_coordinates(
  vector<FP_NR<mpfr_t>> &coordinates,
  const vector<Z_NR<mpz_t>> &target,
  const ZZ_mat<mpz_t> &A,
  const uint32_t n,
  const uint32_t precision)
{
  /* Sanity checks of the dimensions of the A and G matrices. */
  if ((A.get_rows() != (int)(n + 1)) || (A.get_cols() != (int)(n + 1))) {
    critical("solve_left_coordinates(): Matrix A has incorrect dimensions.");
  }

  /* Sanity checks of the size of the target vector. */
  if (target.size() != (n + 1)) {
    critical("solve_left_coordinates(): "
      "The target vector has incorrect dimensions.");
  }

  /* Sanity checks of the size of the coordinate vector. */
  if (coordinates.size() != (n + 1)) {
    coordinates.resize(n + 1);
  }

  /* Set the precision of the coordinate vector. */
  for (uint32_t j = 0; j <= n; j++) {
    mpfr_set_prec(coordinates[j].get_data(), precision);
  }

  /* Copy the A matrix to a temporary matrix A_tmp. */
  FP_mat<mpfr_t> A_tmp(n + 1, n + 1);

  for (uint32_t i = 0; i <= n; i++) {
    for (uint32_t j = 0; j <= n; j++) {
      mpfr_set_prec(A_tmp[i][j].get_data(), precision);

      mpfr_set_z(A_tmp[i][j].get_data(), A[i][j].get_data(), MPFR_RNDN);
    }
  }

  /* Initialize A_inv to the identity. */
  FP_mat<mpfr_t> A_inv(n + 1, n + 1);

  for (uint32_t i = 0; i <= n; i++) {
    for (uint32_t j = 0; j <= n; j++) {
      mpfr_set_prec(A_inv[i][j].get_data(), precision);

      mpfr_set_ui(A_inv[i][j].get_data(), (i == j) ? 1 : 0, MPFR_RNDN);
    }
  }

  /* Perform Gaussian elimination to compute the inverse of A. */
  mpfr_t c;
  mpfr_init2(c, precision);

  mpfr_t x;
  mpfr_init2(x, precision);

  for (uint32_t i = 0; i <= n; i++) {
    /* Let tmp = 1/x for x the element at index (i, i) in the matrix A_tmp. */
    mpfr_set_ui(c, 1, MPFR_RNDN);

    /* Test if x is zero. */
    if (0 == mpfr_cmp_ui(A_tmp[i][i].get_data(), 0)) {
      /* Add a later row to this row to make x non-zero. */
      for (uint32_t j = i + 1; j <= n; j++) {
        if (0 != mpfr_cmp_ui(A_tmp[j][i].get_data(), 0)) {
          /* Add row j to row i to make the coefficient non-zero.. */
          for (uint32_t k = i; k <= n; k++) {
            mpfr_add(A_tmp[i][k].get_data(), A_tmp[i][k].get_data(), 
              A_tmp[j][k].get_data(), MPFR_RNDN);
            mpfr_add(A_inv[i][k].get_data(), A_inv[i][k].get_data(), 
              A_inv[j][k].get_data(), MPFR_RNDN);
          }

          break;
        }
      }

      /* Not full rank. */
      if (0 == mpfr_cmp_ui(A_tmp[i][i].get_data(), 0)) {
        critical("solve_left_coordinates(): Division by zero.");
      }
    }

    mpfr_div(c, c, A_tmp[i][i].get_data(), MPFR_RNDN);

    /* Multiply row i of the matrix A_tmp with c. */
    for (uint32_t j = 0; j <= n; j++) {
      mpfr_mul(A_tmp[i][j].get_data(), A_tmp[i][j].get_data(), c, MPFR_RNDN);
      mpfr_mul(A_inv[i][j].get_data(), A_inv[i][j].get_data(), c, MPFR_RNDN);
    }

    /* Subtract suitable multiples of row i from rows j = 0, 1, .., u - 1. */
    for (uint32_t j = 0; j < i; j++) {
      /* Take minus the element in column i on row j as the multiple c. */
      mpfr_neg(c, A_tmp[j][i].get_data(), MPFR_RNDN);

      /* Subtract c times row i from row j. */
      for (uint32_t k = 0; k <= n; k++) {
        mpfr_mul(x, c, A_tmp[i][k].get_data(), MPFR_RNDN);
        mpfr_add(A_tmp[j][k].get_data(), A_tmp[j][k].get_data(), x, MPFR_RNDN);

        mpfr_mul(x, c, A_inv[i][k].get_data(), MPFR_RNDN);
        mpfr_add(A_inv[j][k].get_data(), A_inv[j][k].get_data(), x, MPFR_RNDN);
      }
    }

    /* Subtract suitable multiples of row i from rows j = i + 1, ..., n. */
    for (uint32_t j = i + 1; j <= n; j++) {
      /* Take minus the element in column i on row j as the multiple c. */
      mpfr_neg(c, A_tmp[j][i].get_data(), MPFR_RNDN);

      /* Subtract c times row i from row j. */
      for (uint32_t k = 0; k <= n; k++) {
        mpfr_mul(x, c, A_tmp[i][k].get_data(), MPFR_RNDN);
        mpfr_add(A_tmp[j][k].get_data(), A_tmp[j][k].get_data(), x, MPFR_RNDN);

        mpfr_mul(x, c, A_inv[i][k].get_data(), MPFR_RNDN);
        mpfr_add(A_inv[j][k].get_data(), A_inv[j][k].get_data(), x, MPFR_RNDN);
      }
    }
  }

  /* Compute target * A^-1 to find the solution. */
  for (uint32_t j = 0; j <= n; j++) {
    mpfr_set_ui(coordinates[j].get_data(), 0, MPFR_RNDN);

    for (uint32_t i = 0; i <= n; i++) {
      mpfr_mul_z(x, A_inv[i][j].get_data(), target[i].get_data(), MPFR_RNDN);

      mpfr_add(
        coordinates[j].get_data(),
        coordinates[j].get_data(),
        x,
        MPFR_RNDN);
    }
  }

  /* Clear memory. */
  mpfr_clear(x);
  mpfr_clear(c);

  A_tmp.clear();
  A_inv.clear();
}
