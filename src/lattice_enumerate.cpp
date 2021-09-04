/*!
 * \file    lattice_enumerate.cpp
 * \ingroup lattice_enumerate
 *
 * \brief   The definition of functions for recovering d and r by attempting to
 *          enumerate all vectors within a ball in the lattice L.
 */

#include "lattice_enumerate.h"

#include "common.h"
#include "errors.h"
#include "lattice.h"
#include "lattice_algebra.h"
#include "lattice_babai.h"
#include "lattice_gso.h"
#include "lattice_reduce.h"
#include "lattice_smoothness.h"
#include "lattice_solve.h"
#include "math.h"
#include "parameters.h"
#include "sample.h"
#include "timer.h"

#include <gmp.h>
#include <mpfr.h>

#include <fplll/fplll.h>

#include <stdint.h>

#include <vector>

using namespace std;
using namespace fplll;

/*!
 * \brief   The factor by which the estimated square norm of the shortest
 *          vector in the lattice should be reduced.
 *
 * This reduction is performed so as to attempt to ensure that we do not miss
 * the shortest non-zero vector in the lattice.
 */
#define INITIAL_RADIUS_REDUCTION_FACTOR 1024

/*!
 * \brief   The maximum number of radius doublings.
 *
 * Enumeration is performed in balls of increasingly greater radius. The initial
 * radius is given by the reduced estimated norm of the shortest non-zero vector
 * in the reduced lattice basis. This norm is then doubled repeatedly.
 */
#define MAX_RADIUS_DOUBLINGS 128

/*!
 * \brief   Computes uhat, the optimum for the k:th entry in the cu_coordinates
 *          vector, and the minimum and maximum offsets in uhat.
 *
 * Let A be the (n + 1) x (n + 1) basis matrix of the lattice L.
 *
 * Let G be the Gram-Schmidt orthogonalized (GSO) basis matrix for A, and M the
 * triangular matrix of Gram-Schmidt projection factors, so that A = M * G. Note
 * that gram_schmidt_orthogonalization() computes these matrices given A.
 *
 * Let u = cu_coordinates * A and v = cv_coordinates * A, where cu_coordinates
 * and cv_coordinates are row vectors.
 *
 * Let k be an integer index on the interval between 1 and n + 1.
 *
 * Assume that cv_coordinates is completely populated, and that the entries at
 * indices k, .., n + 1 of cu_coordinates are populated.
 *
 * This function then computes the minimum and maximum offsets in the k:th entry
 * in the cu_coordinates vector, such that
 *
 *   | \\pi_k(u - v) |^2 =
 *     | \\pi_k((cu_coordinates - cv_coordinates) * A) |^2 =
 *       | \\pi_k((cu_coordinates - cv_coordinates) * M * G) |^2 < R^2.
 *
 * where R^2 is square_radius. This is accomplished using standard techniques
 * of looking at the Gram-Schmidt projection of u - v onto planes k, .., n + 1.
 *
 * \internal
 *
 * \param[out]     uhat             An integer to set to the optimum for the
 *                                  k:th entry in the cu_coordinates vector.
 * \param[out]     min_uhat_offset  An integer to set to the minimum offset in
 *                                  the k:th entry in the cu_coordinates vector.
 * \param[out]     max_uhat_offset  An integer to set to the maximum offset in
 *                                  the k:th entry in the cu_coordinates vector.
 * \param[in]      square_radius    The square radius R^2 of the ball to
 *                                  enumerate.
 * \param[in, out] cu_coordinates   The coordinates of u in the A basis. The
 *                                  k:th entry will be set to uhat.
 * \param[in]      cv_coordinates   The coordinates of v in the A basis.
 * \param[in]      M                The triangular matrix M of Gram-Schmidt
 *                                  projection factors, such that A = M * G for
 *                                  G the Gram-Schmidt orthogonalized (GSO)
 *                                  basis matrix for the A basis.
 * \param[in]      G_square_norms   The square norms of the rows of G.
 * \param[in]      k                The index into the cu_coordinates vector.
 *                                  Must  be on the interval [1, n + 1]. The
 *                                  entries at indices k + 1, .., n + 1 of the
 *                                  cu_oordinates vector must be populated.
 * \param[in]      n                The number of runs used to setup the A
 *                                  matrix from which G and M are computed.
 * \param[in]      precision        The floating point precision.
 *
 * \return Returns #TRUE if | \\pi_k(u - v) |^2 > R^2, and #FALSE otherwise.
 */
static bool compute_uhat_offsets(
  mpz_t uhat,
  mpz_t min_uhat_offset,
  mpz_t max_uhat_offset,
  const mpfr_t square_radius,
  vector<Z_NR<mpz_t>> &cu_coordinates,
  const vector<FP_NR<mpfr_t>> &cv_coordinates,
  const FP_mat<mpfr_t> &M,
  const vector<FP_NR<mpfr_t>> &G_square_norms,
  const uint32_t k,
  const uint32_t n,
  const uint32_t precision)
{
  bool result = TRUE;

  mpfr_t square_norm;
  mpfr_init2(square_norm, precision);

  mpfr_t tmp;
  mpfr_init2(tmp, precision);

  mpfr_t tmp2;
  mpfr_init2(tmp2, precision);

  /* 1. Compute uhat (using tmp as temporary storage for uhat). */
  mpfr_set(tmp, cv_coordinates[k - 1].get_data(), MPFR_RNDN);

  for (uint32_t j = k + 1; j <= n + 1; j++) {
    mpfr_sub_z(
      tmp2,
      cv_coordinates[j - 1].get_data(),
      cu_coordinates[j - 1].get_data(),
      MPFR_RNDN);

    mpfr_mul(tmp2, tmp2, M[k - 1][j - 1].get_data(), MPFR_RNDN);
    mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
  }

  mpfr_round(tmp, tmp);
  mpfr_get_z(uhat, tmp, MPFR_RNDN);

  /* 1.1 Store uhat in the k:th position in the cu_coordinates vector. */
  mpz_set(cu_coordinates[k - 1].get_data(), uhat);
    /* Note: k is indexed from one. */

  /* Note: The parameter k is indexed from one, whilst the indices i and j below
   * are indexed from zero. Therefore, we subtract one from k when used. */

  /* 2. Let c = c_u - c_v. */
  vector<FP_NR<mpfr_t>> c(n + 1);

  for (uint32_t j = k - 1; j <= n; j++) { /* Note: k is indexed from one. */
    mpfr_set_prec(c[j].get_data(), precision);

    mpfr_sub_z(
      c[j].get_data(),
      cv_coordinates[j].get_data(),
      cu_coordinates[j].get_data(),
      MPFR_RNDN);
    mpfr_neg(c[j].get_data(), c[j].get_data(), MPFR_RNDN);
  }

  /* 3. Compute the square norm for the projection onto entries k + 1, .., n. */
  mpfr_set_ui(square_norm, 0, MPFR_RNDN);

  for (uint32_t j = k; j <= n; j++) { /* Note: k is indexed from one. */
    mpfr_set(tmp, c[j].get_data(), MPFR_RNDN);

    for (uint32_t i = j + 1; i <= n; i++) {
      mpfr_mul(tmp2, M[i][j].get_data(), c[i].get_data(), MPFR_RNDN);
      mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
    }

    mpfr_sqr(tmp, tmp, MPFR_RNDN);
    mpfr_mul(tmp, tmp, G_square_norms[j].get_data(), MPFR_RNDN);

    mpfr_add(square_norm, square_norm, tmp, MPFR_RNDN);
  }

  /* 4. Compute the contribution for the projection onto entry k to tmp. */
  mpfr_set(tmp, c[k - 1].get_data(), MPFR_RNDN);
  for (uint32_t i = k; i <= n; i++) { /* Note: k is indexed from one. */
    mpfr_mul(tmp2, M[i][k - 1].get_data(), c[i].get_data(), MPFR_RNDN);
    mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
  }

  mpfr_sqr(tmp2, tmp, MPFR_RNDN);
  mpfr_mul(tmp2, tmp2, G_square_norms[k - 1].get_data(), MPFR_RNDN);

  /* Sum up the square norm for the projection onto entries k, .., n to tmp2. */
  mpfr_add(tmp2, tmp2, square_norm, MPFR_RNDN);

  /* 5. Check if the square norm exceeds the square radius R^2. */
  if (mpfr_cmp(tmp2, square_radius) > 0) {
    /* 5.1 If so, set the minimum and maximum offset to zero. */
    mpz_set_ui(min_uhat_offset, 0);
    mpz_set_ui(max_uhat_offset, 0);

    result = FALSE; /* Empty interval. */
  } else {
    /* 5.2 If not, compute the minimum and maximum offsets. */

    mpfr_t residual;
    mpfr_init2(residual, precision);

    mpfr_sub(residual, square_radius, square_norm, MPFR_RNDN);
    mpfr_div(residual, residual, G_square_norms[k - 1].get_data(), MPFR_RNDN);
    mpfr_sqrt(residual, residual, MPFR_RNDN);

    /* Note: We need | offset + tmp | < residual, where tmp may be negative. */

    /* 5.2.1. Compute min_uhat_offset, the minimum offset. */
    mpfr_neg(tmp2, residual, MPFR_RNDN); /* tmp2 = -residual */
    mpfr_sub(tmp2, tmp2, tmp, MPFR_RNDN); /* tmp2 = -residual - tmp */

    if (mpfr_sgn(tmp2) < 0) {
      mpfr_ceil(tmp2, tmp2);
    } else {
      mpfr_floor(tmp2, tmp2);
    }
    mpfr_get_z(min_uhat_offset, tmp2, MPFR_RNDN);

    /* 5.2.2. Compute max_uhat_offset, the maximum offset. */
    mpfr_sub(tmp2, residual, tmp, MPFR_RNDN); /* tmp2 = residual - tmp */

    if (mpfr_sgn(tmp2) < 0) {
      mpfr_ceil(tmp2, tmp2);
    } else {
      mpfr_floor(tmp2, tmp2);
    }
    mpfr_get_z(max_uhat_offset, tmp2, MPFR_RNDN);

    result = TRUE; /* Non-empty interval. */

    /* Clear temporary memory. */
    mpfr_clear(residual);
  }

  /* Clear memory. */
  mpfr_clear(square_norm);
  mpfr_clear(tmp);
  mpfr_clear(tmp2);

  c.clear();

  /* Signal result. */
  return result;
}

/*!
 * \brief   An internal function that is called recursively to enumerate L.
 *
 * This function implements lattice enumeration, by building on ideas originally
 * described by Kannan [1], with a few optimizations:
 *
 *   1. We use standard techniques from the literature, such as incremental
 *      Gram-Schmidt projections, to enumerate all vectors in L within a ball of
 *      a given radius, see the internal function offset_uhat(). We test each
 *      plausible candidate against the known solution d_or_r. (This explains
 *      why d_or_r must be passed to this function.)
 *
 *   2. We use that the last component is d or r (or dr in the scaled basis for
 *      diagonal distributions), and hence on a restricted known interval
 *      min_d_or_r <= d, r <= max_d_or_r.
 *
 *   3. We support detecting partially smooth r.
 *
 * [1] Kannan, R.: Improved algorithms for integer programming and related
 * lattice problems. In: Proceedings of the 15th Symposium on the Theory of
 * Computing. STOC 1983. ACM Press, pp. 99—108 (1983).
 *
 * Note that this a fairly straightforward implementation of lattice enumeration
 * that has proven sufficient for these purposes. There are more efficient
 * methods for enumerating lattices, and it is possible to implement this
 * function more efficiently. However, doing so has not proved necessary.
 *
 * \internal
 *
 * \param[out]     status          An entry in the #Lattice_Status_Recovery
 *                                 enumeration that is used to signal
 *                                 information on whether d_or_r was
 *                                 successfully recovered.
 * \param[in]      square_radius   The square radius R^2 of the ball to
 *                                 enumerate.
 * \param[in, out] cu_coordinates  The coordinates of u in the A basis.
 * \param[in]      cv_coordinates  The coordinates of v in the A basis.
 * \param[in]      M               The triangular matrix M of Gram-Schmidt
 *                                 projection factors, such that A = M * G for
 *                                 G the Gram-Schmidt orthogonalized (GSO) basis
 *                                 matrix for the A basis.
 * \param[in]      G_square_norms  The square norms of the rows of G.
 * \param[in]      k               The index into the cu_coordinates vector.
 *                                 Specifies the depth in the enumeration tree.
 * \param[in]      n               The number of runs used to setup the A
 *                                 matrix from which G and M are computed.
 * \param[in]      min_d_or_r      The minimum value of the correct solution.
 *                                 Used to prune the enumeration at the last
 *                                 level of the enumeration tree.
 * \param[in]      max_d_or_r      The maximum value of the correct solution.
 *                                 Used to prune the enumeration at the last
 *                                 level of the enumeration tree.
 * \param[in]      d_or_r          The correct solution. Used only to test if
 *                                 the candidate solutions found are correct.
 * \param[in]      m               The length m in bits of r. Used only when
 *                                 solving for r, and detect_smooth is set to
 *                                 #TRUE.
 * \param[in]      precision       The floating point precision.
 * \param[in]      detect_smooth   A flag that should be set to #TRUE if the
 *                                 function should return d_or_r' for
 *                                 d_or_r = d_or_r' * z, where z is smooth.
 * \param[in]      timeout         The timeout after which the enumeration will
 *                                 be aborted. Set to 0 to disable the timeout.
 * \param[in]      timer           A timer used to measure how much time has
 *                                 elapsed since the first call to this
 *                                 function. Must be initialized and started by
 *                                 the caller.
 */
static void lattice_enumerate_inner(
  Lattice_Status_Recovery * const status,
  const mpfr_t square_radius,
  vector<Z_NR<mpz_t>> &cu_coordinates,
  const vector<FP_NR<mpfr_t>> &cv_coordinates,
  const ZZ_mat<mpz_t> &A,
  const vector<FP_NR<mpfr_t>> &G_square_norms,
  const FP_mat<mpfr_t> &M,
  const uint32_t k,
  const uint32_t n,
  const mpz_t min_d_or_r,
  const mpz_t max_d_or_r,
  const mpz_t d_or_r,
  const uint32_t m,
  const uint32_t precision,
  const bool detect_smooth,
  const uint64_t timeout,
  Timer * const timer)
{
  /* Sanity checks. */
  if ((mpz_cmp(min_d_or_r, d_or_r) > 0) || (mpz_cmp(max_d_or_r, d_or_r) < 0)) {
    critical("lattice_enumerate_inner(): Sanity checks not respected.");
  }

  /* Timeout check. */
  if ((0 != timeout) && (timer_stop(timer) / (1000 * 1000) > timeout)) {
    (*status) = LATTICE_STATUS_TIMEOUT;

    return; /* Stop. */
  }

  /* Initially set the status to not recovered. */
  (*status) = LATTICE_STATUS_NOT_RECOVERED;

  /* Declare variables and constants. */
  mpz_t uhat;
  mpz_init(uhat);

  mpz_t min_uhat_offset;
  mpz_init(min_uhat_offset);

  mpz_t abs_min_uhat_offset;
  mpz_init(abs_min_uhat_offset);

  mpz_t max_uhat_offset;
  mpz_init(max_uhat_offset);

  mpz_t bound;
  mpz_init(bound);

  mpz_t offset;
  mpz_init(offset);

  mpz_t candidate_d_or_r;
  mpz_init(candidate_d_or_r);

  mpz_t increment_d_or_r;
  mpz_init(increment_d_or_r);

  mpz_t last_component_u;
  mpz_init(last_component_u);

  mpz_t tmp_z;
  mpz_init(tmp_z);

  mpfr_t tmp_f;
  mpfr_init2(tmp_f, precision);

  /* Note: The above variables are cleared at the end of this function. */

  /* Compute uhat, and the minimum and maximum offsets in uhat. */
  const bool result = compute_uhat_offsets(
                        uhat,
                        min_uhat_offset,
                        max_uhat_offset,
                        square_radius,
                        cu_coordinates,
                        cv_coordinates,
                        M,
                        G_square_norms,
                        k,
                        n,
                        precision);
  if (TRUE != result) {
    /* Stop. */
    goto lattice_enumerate_inner_clear;
  }

  if (k > 1) {
    /* Recursively call lattice_enumerate_inner() for the interval. */

    /* Compute bound = max(abs(min_uhat_offset), max_uhat_offset). */
    if (mpz_cmp_ui(max_uhat_offset, 0) < 0) {
      critical("lattice_enumerate_inner(): Expected max_uhat_offset >= 0.");
    }

    mpz_abs(abs_min_uhat_offset, min_uhat_offset);
    if (mpz_cmp(max_uhat_offset, abs_min_uhat_offset) > 0) {
      mpz_set(bound, max_uhat_offset);
    } else {
      mpz_set(bound, abs_min_uhat_offset);
    }

    /* Iterate with offset 0, .., bound with positive and negative sign. */
    mpz_set_si(offset, -1);

    while (mpz_cmp(offset, bound) <= 0) {
      mpz_add_ui(offset, offset, 1);

      if (mpz_cmp(offset, max_uhat_offset) <= 0) {
        mpz_add(cu_coordinates[k - 1].get_data(), uhat, offset);

        lattice_enumerate_inner(
          status,
          square_radius,
          cu_coordinates,
          cv_coordinates,
          A,
          G_square_norms,
          M,
          k - 1,
          n,
          min_d_or_r,
          max_d_or_r,
          d_or_r,
          m,
          precision,
          detect_smooth,
          timeout,
          timer);
        if (LATTICE_STATUS_NOT_RECOVERED != (*status)) {
          break; /* Recovery was successful or an error occurred. */
        }
      }

      if (mpz_cmp_ui(offset, 0) == 0) {
        continue;
      }

      if (mpz_cmp(offset, abs_min_uhat_offset) <= 0) {
        mpz_sub(cu_coordinates[k - 1].get_data(), uhat, offset);

        lattice_enumerate_inner(
          status,
          square_radius,
          cu_coordinates,
          cv_coordinates,
          A,
          G_square_norms,
          M,
          k - 1,
          n,
          min_d_or_r,
          max_d_or_r,
          d_or_r,
          m,
          precision,
          detect_smooth,
          timeout,
          timer);
        if (LATTICE_STATUS_NOT_RECOVERED != (*status)) {
          break; /* Recovery was successful or an error occurred. */
        }
      }
    }

    /* Done. */

  } else if (k == 1) {
    /* 1. Compute the last component of the vector u when there is no offset. */
    mpz_set(cu_coordinates[0].get_data(), uhat);

    mpz_set_ui(last_component_u, 0);

    for (uint32_t i = 0; i <= n; i++) {
      mpz_mul(tmp_z, A[i][n].get_data(), cu_coordinates[i].get_data());
      mpz_add(last_component_u, last_component_u, tmp_z);
    }

    mpz_set(candidate_d_or_r, last_component_u);

    /* 2. Find the increment in the last component in the offset. */
    mpz_abs(increment_d_or_r, A[0][n].get_data());

    /* 3. Check if there are offsets on [min_uhat_offset, max_uhat_offset] that
     * respect the requirement that the last component of u is on the interval
     * [min_d_or_r, max_d_or_r], and check if any of these yield d_or_r. */

    /* 3.1 Adjust the minimum offset to the requirement in 3. */
    mpz_mul(tmp_z, min_uhat_offset, increment_d_or_r);
    mpz_add(tmp_z, last_component_u, tmp_z);

    if (mpz_cmp(tmp_z, min_d_or_r) < 0) {
      /* Check how much min_uhat_offset must be incremented. */
      mpz_sub(tmp_z, min_d_or_r, tmp_z);
      mpfr_set_z(tmp_f, tmp_z, MPFR_RNDN);
      mpfr_div_z(tmp_f, tmp_f, increment_d_or_r, MPFR_RNDN);
      mpfr_ceil(tmp_f, tmp_f);
      mpfr_get_z(tmp_z, tmp_f, MPFR_RNDN);

      mpz_add(min_uhat_offset, min_uhat_offset, tmp_z);
    }

    if (mpz_cmp(min_uhat_offset, max_uhat_offset) > 0) {
      /* Stop. Empty offset interval. */
      goto lattice_enumerate_inner_clear;
    }

    /* 3.2. Adjust the maximum offset to the requirement in 3. */
    mpz_mul(tmp_z, max_uhat_offset, increment_d_or_r);
    mpz_add(tmp_z, last_component_u, tmp_z);

    if (mpz_cmp(tmp_z, max_d_or_r) > 0) {
      /* Check how much max_uhat_offset must be decremented. */
      mpz_sub(tmp_z, tmp_z, max_d_or_r);
      mpfr_set_z(tmp_f, tmp_z, MPFR_RNDN);
      mpfr_div_z(tmp_f, tmp_f, increment_d_or_r, MPFR_RNDN);
      mpfr_ceil(tmp_f, tmp_f);
      mpfr_get_z(tmp_z, tmp_f, MPFR_RNDN);

      mpz_sub(max_uhat_offset, max_uhat_offset, tmp_z);
    }

    if (mpz_cmp(min_uhat_offset, max_uhat_offset) > 0) {
      /* Stop. Empty offset interval. */
      goto lattice_enumerate_inner_clear;
    }

    /* Iterate with offset 0, .., bound with positive and negative sign. */
    mpz_set(offset, min_uhat_offset);

    mpz_mul(tmp_z, min_uhat_offset, increment_d_or_r);
    mpz_add(candidate_d_or_r, last_component_u, tmp_z);

    uint32_t timeout_counter = 0;

    while (mpz_cmp(offset, max_uhat_offset) <= 0) {
      /* Candidate check. */
      if (mpz_cmp(candidate_d_or_r, d_or_r) == 0) {
        /* Set the status to recovered. */
        (*status) = LATTICE_STATUS_RECOVERED_ENUMERATION;

        /* Stop. */
        goto lattice_enumerate_inner_clear;
      } else if (detect_smooth && (0 != mpz_cmp_ui(candidate_d_or_r, 0))) {
        mpz_mod(tmp_z, d_or_r, candidate_d_or_r);

        if (mpz_cmp_ui(tmp_z, 0) == 0) { /* candidate_d_or_r divides d_or_r */
          mpz_div(tmp_z, d_or_r, candidate_d_or_r);

          const bool smooth = lattice_smoothness_is_smooth(
                                tmp_z,
                                LATTICE_SMOOTHNESS_CONSTANT_C,
                                m);

          if (smooth) {
            /* Set the status to recovered. */
            (*status) = LATTICE_STATUS_RECOVERED_ENUMERATION_SMOOTH;

            /* Stop. */
            goto lattice_enumerate_inner_clear;
          }
        }
      }

      /* Timeout check. */
      timeout_counter = (timeout_counter + 1) & 0xff;
      if (0 == timeout_counter) {
        if ((0 != timeout) && (timer_stop(timer) / (1000 * 1000) > timeout)) {
          /* Set the status to timeout. */
          (*status) = LATTICE_STATUS_TIMEOUT;

          /* Stop. */
          goto lattice_enumerate_inner_clear;
        }
      }

      /* Increment the offset and the candidate. */
      mpz_add_ui(offset, offset, 1);
      mpz_add(candidate_d_or_r, candidate_d_or_r, increment_d_or_r);
    }

    /* Done. */

  } else {
    critical("lattice_enumerate_inner(): Expected k to be on [1, n + 1].");
  }

lattice_enumerate_inner_clear:

  switch (*status) {
    case LATTICE_STATUS_NOT_RECOVERED:
    case LATTICE_STATUS_TIMEOUT:
      /* Zeroize the k:th entry in the cu_coordinates vectors. This is not
       * necessary, but avoids leaving data in the vector when backtracking. */
      mpz_set_ui(cu_coordinates[k - 1].get_data(), 0);
      break;

    default:
      break;
  }

  /* Clear memory. */
  mpz_clear(uhat);
  mpz_clear(min_uhat_offset);
  mpz_clear(abs_min_uhat_offset);
  mpz_clear(max_uhat_offset);
  mpz_clear(bound);
  mpz_clear(offset);

  mpz_clear(candidate_d_or_r);
  mpz_clear(increment_d_or_r);
  mpz_clear(last_component_u);

  mpz_clear(tmp_z);
  mpfr_clear(tmp_f);
}

void lattice_enumerate_reduced_basis_for_d(
  Lattice_Status_Recovery * const status_d,
  const ZZ_mat<mpz_t> &A,
  const FP_mat<mpfr_t> &G,
  const FP_mat<mpfr_t> &M,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  const uint32_t precision,
  const bool detect_smooth_r,
  const uint64_t timeout)
{
  /* Check dimensions. */
  if ((A.get_rows() != (int)(n + 1)) || (A.get_cols() != (int)(n + 1))) {
    critical("lattice_enumerate_reduced_basis_for_d(): "
      "Incorrect dimensions for the matrix A.");
  }

  if ((G.get_rows() != (int)(n + 1)) || (G.get_cols() != (int)(n + 1))) {
    critical("lattice_enumerate_reduced_basis_for_d(): "
      "Incorrect dimensions for the matrix G.");
  }

  if ((M.get_rows() != (int)(n + 1)) || (M.get_cols() != (int)(n + 1))) {
    critical("lattice_enumerate_reduced_basis_for_d(): "
      "Incorrect dimensions for the matrix M.");
  }

  /* Setup constants. */
  mpz_t pow2m;
  mpz_init(pow2m);
  mpz_set_ui(pow2m, 0);
  mpz_setbit(pow2m, parameters->m);

  mpz_t pow2lm;
  mpz_init(pow2lm);
  mpz_set_ui(pow2lm, 0);
  mpz_setbit(pow2lm, parameters->l + parameters->m);

  mpz_t pow2lm_half;
  mpz_init(pow2lm_half);
  mpz_set_ui(pow2lm_half, 0);
  mpz_setbit(pow2lm_half, parameters->l + parameters->m - 1);

  /* Setup variables. */
  mpz_t candidate_d;
  mpz_init(candidate_d);

  mpfr_t tmp;
  mpfr_init2(tmp, precision);

  /* Compute the square norms of the rows of the G matrix. */
  vector<FP_NR<mpfr_t>> G_square_norms(n + 1);

  for (uint32_t i = 0; i <= n; i++) {
    mpfr_set_prec(G_square_norms[i].get_data(), precision);
    mpfr_set_ui(G_square_norms[i].get_data(), 0, MPFR_RNDN);

    for (uint32_t j = 0; j <= n; j++) {
      mpfr_mul(tmp, G[i][j].get_data(), G[i][j].get_data(), MPFR_RNDN);

      mpfr_add(
          G_square_norms[i].get_data(),
          G_square_norms[i].get_data(),
          tmp,
          MPFR_RNDN);
    }
  }

  /* Setup the vector
   *
   *   v = ({-2^m k_1}_{2^{m+l}}, .., {-2^m k_n}_{2^{m+l}}, 0)
   * 
   * that defines the origin of the ball to enumerate. */
  vector<Z_NR<mpz_t>> v(n + 1);

  for (uint32_t j = 0; j < n; j++) {
    mpz_mul(v[j].get_data(), ks[j], pow2m);
    mpz_neg(v[j].get_data(), v[j].get_data());
    mpz_mod(v[j].get_data(), v[j].get_data(), pow2lm);

    if (mpz_cmp(v[j].get_data(), pow2lm_half) >= 0) {
      mpz_sub(v[j].get_data(), v[j].get_data(), pow2lm);
    }
  }

  mpz_set_ui(v[n].get_data(), 0);

  /* Solve v = c_v A for c_v. */
  vector<FP_NR<mpfr_t>> cv_coordinates(n + 1);
  solve_left_coordinates(cv_coordinates, v, A, n, precision);

  /* Compute the closest vector. */
  vector<Z_NR<mpz_t>> solution(n + 1);
  babai_closest_vector(solution, v, G, A, n, precision);

  /* Extract the candidate d. */
  mpz_abs(candidate_d, solution[n].get_data());

  /* Set the status to indicate failed recovery. */
  (*status_d) = LATTICE_STATUS_NOT_RECOVERED;

  if (mpz_cmp(candidate_d, parameters->d) == 0) {
    (*status_d) = LATTICE_STATUS_RECOVERED_IMMEDIATE;
  } else {
    mpz_t multiple;
    mpz_init(multiple);

    mpz_t tmp_d;
    mpz_init(tmp_d);

    mpz_set(candidate_d, solution[n].get_data());

    for (uint32_t i = 1; i <= 256; i++) {
      mpz_mul_ui(multiple, A[0][n].get_data(), i);

      mpz_add(tmp_d, candidate_d, multiple);
      mpz_abs(tmp_d, tmp_d);
      if (mpz_cmp(tmp_d, parameters->d) == 0) {
        (*status_d) = LATTICE_STATUS_RECOVERED_SEARCH;
        break;
      }

      mpz_sub(tmp_d, candidate_d, multiple);
      mpz_abs(tmp_d, tmp_d);
      if (mpz_cmp(tmp_d, parameters->d) == 0) {
        (*status_d) = LATTICE_STATUS_RECOVERED_SEARCH;
        break;
      }
    }

    mpz_clear(multiple);
    mpz_clear(tmp_d);
  }

  if (LATTICE_STATUS_NOT_RECOVERED == (*status_d)) {
    /* Step 2: Let c_u be the zero coordinate vector. */
    vector<Z_NR<mpz_t>> cu_coordinates(n + 1);
    for (uint32_t j = 0; j <= n; j++) {
      mpz_set_ui(cu_coordinates[j].get_data(), 0);
    }

    /* Compute the square norm of the shortest non-zero vector in the reduced
     * lattice basis A. */
    mpfr_t square_radius;
    mpfr_init2(square_radius, precision);
    mpfr_set_ui(square_radius, 0, MPFR_RNDN);

    for (uint32_t j = 0; j <= n; j++) {
      mpfr_set_z(tmp, A[0][j].get_data(), MPFR_RNDN);
      mpfr_mul(tmp, tmp, tmp, MPFR_RNDN);
      mpfr_add(square_radius, square_radius, tmp, MPFR_RNDN);
    }

    /* To attempt to ensure that we do not miss the shortest vector in the
     * lattice, reduce the square norm by #INITIAL_RADIUS_REDUCTION_FACTOR. */
    mpfr_div_ui(
      square_radius,
      square_radius,
      INITIAL_RADIUS_REDUCTION_FACTOR,
      MPFR_RNDN);

    /* Start a timer. */
    Timer timer;
    timer_start(&timer);

    /* Execute the enumeration. */
    (*status_d) = LATTICE_STATUS_NOT_RECOVERED;

    mpz_t min_d;
    mpz_init(min_d);
    mpz_set_ui(min_d, 0); /* min_d = 0 */

    mpz_t max_d;
    mpz_init(max_d);
    mpz_set_ui(max_d, 0);
    mpz_setbit(max_d, parameters->m); /* max_d = 2^m */

    mpz_t d;
    mpz_init(d);
    mpz_set(d, parameters->d);

    if (detect_smooth_r) {
      /* Check if there is a small multiple of r in the lattice. This may occur
       * when solving for a general discrete logarithm if r is smooth. */
      mpz_t tmp_r;
      mpz_init(tmp_r);

      mpz_abs(tmp_r, A[0][n].get_data());
      mpz_mod(tmp_r, parameters->r, tmp_r);

      if (0 == mpz_cmp_ui(tmp_r, 0)) {
        mpz_abs(tmp_r, A[0][n].get_data());
        mpz_div(tmp_r, parameters->r, tmp_r); /* tmp_r = r / r' */

        if (mpz_cmp_ui(tmp_r, 1) > 0) {
          const bool smooth = lattice_smoothness_is_smooth(
                                tmp_r,
                                LATTICE_SMOOTHNESS_CONSTANT_C,
                                parameters->m);
          if (smooth) {
            mpz_div(max_d, max_d, tmp_r); /* max_d = 2^m / (r / r') */
          }
        }
      }

      mpz_clear(tmp_r);
    }

    for (uint32_t t = 0; t < MAX_RADIUS_DOUBLINGS; t++) {
      lattice_enumerate_inner(
        status_d,
        square_radius,
        cu_coordinates,
        cv_coordinates,
        A,
        G_square_norms,
        M,
        n + 1, /* = k */
        n,
        min_d,
        max_d,
        d,
        parameters->m, /* = m */
        precision,
        FALSE, /* = detect_smooth */
        timeout,
        &timer);

      /* Note that we set detect_smooth to FALSE above. If r is partially very
       * smooth, there may be an artificially short vector in the lattice,
       * leading us to enumerate very many vectors. To handle this, however, we
       * reduce the upper interval for d, see the above code. The detect_smooth
       * flag is only relevant when enumerating for r. */

      if (LATTICE_STATUS_NOT_RECOVERED != (*status_d)) {
        /* If the fact that r has smooth factors was used, report it. */
        if (detect_smooth_r) {
          if (LATTICE_STATUS_RECOVERED_ENUMERATION == (*status_d)) {
            if (mpz_cmp(d, parameters->d) < 0) {
              (*status_d) = LATTICE_STATUS_RECOVERED_ENUMERATION_SMOOTH;
            }
          }
        }

        /* A solution was found or there was a timeout. */
        break;
      }

      /* Double the radius explored in the next iteration. */
      mpfr_mul_ui(square_radius, square_radius, 4, MPFR_RNDN);
    }

    /* Clear memory. */
    mpz_clear(d);

    mpz_clear(min_d);
    mpz_clear(max_d);

    cu_coordinates.clear();

    mpfr_clear(square_radius);
  }

  /* Clear memory. */
  G_square_norms.clear();
  cv_coordinates.clear();
  v.clear();

  mpfr_clear(tmp);

  mpz_clear(pow2m);
  mpz_clear(pow2lm);
  mpz_clear(pow2lm_half);
  mpz_clear(candidate_d);
}

void lattice_enumerate_reduced_basis_for_d_given_r(
  Lattice_Status_Recovery * const status_d,
  const ZZ_mat<mpz_t> &A,
  const FP_mat<mpfr_t> &G,
  const FP_mat<mpfr_t> &M,
  const mpz_t * const ks,
  const uint32_t n,
  const Diagonal_Parameters * const parameters,
  const uint32_t precision,
  const uint64_t timeout)
{
  /* Check dimensions. */
  if ((A.get_rows() != (int)(n + 1)) || (A.get_cols() != (int)(n + 1))) {
    critical("lattice_enumerate_reduced_basis_for_d_given_r(): "
      "Incorrect dimensions for the matrix A.");
  }

  if ((G.get_rows() != (int)(n + 1)) || (G.get_cols() != (int)(n + 1))) {
    critical("lattice_enumerate_reduced_basis_for_d_given_r(): "
      "Incorrect dimensions for the matrix G.");
  }

  if ((M.get_rows() != (int)(n + 1)) || (M.get_cols() != (int)(n + 1))) {
    critical("lattice_enumerate_reduced_basis_for_d_given_r(): "
      "Incorrect dimensions for the matrix M.");
  }

  /* Setup constants. */
  mpz_t pow2m;
  mpz_init(pow2m);
  mpz_set_ui(pow2m, 0);
  mpz_setbit(pow2m, parameters->m);

  mpz_t rpow2lm;
  mpz_init(rpow2lm);
  mpz_set_ui(rpow2lm, 0);
  mpz_setbit(rpow2lm, parameters->l + parameters->m);
  mpz_mul(rpow2lm, rpow2lm, parameters->r);

  mpz_t rpow2lsigma;
  mpz_init(rpow2lsigma);
  mpz_set_ui(rpow2lsigma, 0);
  mpz_setbit(rpow2lsigma, parameters->l - parameters->sigma);
  mpz_mul(rpow2lsigma, rpow2lsigma, parameters->r);

  /* Setup variables. */
  mpz_t candidate_d;
  mpz_init(candidate_d);

  mpz_t target_d;
  mpz_init(target_d);

  mpfr_t tmp;
  mpfr_init2(tmp, precision);

  /* Compute the square norms of the rows of the G matrix. */
  vector<FP_NR<mpfr_t>> G_square_norms(n + 1);

  for (uint32_t i = 0; i <= n; i++) {
    mpfr_set_prec(G_square_norms[i].get_data(), precision);
    mpfr_set_ui(G_square_norms[i].get_data(), 0, MPFR_RNDN);

    for (uint32_t j = 0; j <= n; j++) {
      mpfr_mul(tmp, G[i][j].get_data(), G[i][j].get_data(), MPFR_RNDN);

      mpfr_add(
          G_square_norms[i].get_data(),
          G_square_norms[i].get_data(),
          tmp,
          MPFR_RNDN);
    }
  }

  /* Setup the vector v that defines the origin of the ball to enumerate. */
  vector<Z_NR<mpz_t>> v(n + 1);

  for (uint32_t j = 0; j < n; j++) {
    mpz_mul(v[j].get_data(), ks[j], pow2m);
    mpz_mul(v[j].get_data(), v[j].get_data(), parameters->r);
    mpz_neg(v[j].get_data(), v[j].get_data());

    mod_reduce(v[j].get_data(), rpow2lm);
  }

  mpz_set_ui(v[n].get_data(), 0);

  /* Solve v = c_v A for c_v. */
  vector<FP_NR<mpfr_t>> cv_coordinates(n + 1);
  solve_left_coordinates(cv_coordinates, v, A, n, precision);

  /* Compute the closest vector. */
  vector<Z_NR<mpz_t>> solution(n + 1);
  babai_closest_vector(solution, v, G, A, n, precision);

  /* Compute the target. */
  mpz_mul(target_d, parameters->d, parameters->r);

  /* Extract the candidate d. */
  mpz_abs(candidate_d, solution[n].get_data());

  /* Set the status to indicate failed recovery. */
  (*status_d) = LATTICE_STATUS_NOT_RECOVERED;

  if (mpz_cmp(candidate_d, target_d) == 0) {
    (*status_d) = LATTICE_STATUS_RECOVERED_IMMEDIATE;
  } else {
    mpz_t multiple;
    mpz_init(multiple);

    mpz_t tmp_d;
    mpz_init(tmp_d);

    mpz_set(candidate_d, solution[n].get_data());

    for (uint32_t i = 1; i <= 256; i++) {
      mpz_mul_ui(multiple, A[0][n].get_data(), i);

      mpz_add(tmp_d, candidate_d, multiple);
      mpz_abs(tmp_d, tmp_d);
      if (mpz_cmp(tmp_d, target_d) == 0) {
        (*status_d) = LATTICE_STATUS_RECOVERED_SEARCH;
        break;
      }

      mpz_sub(tmp_d, candidate_d, multiple);
      mpz_abs(tmp_d, tmp_d);
      if (mpz_cmp(tmp_d, target_d) == 0) {
        (*status_d) = LATTICE_STATUS_RECOVERED_SEARCH;
        break;
      }
    }

    mpz_clear(multiple);
    mpz_clear(tmp_d);
  }

  if (LATTICE_STATUS_NOT_RECOVERED == (*status_d)) {
    /* Step 2: Let c_u be the zero coordinate vector. */
    vector<Z_NR<mpz_t>> cu_coordinates(n + 1);
    for (uint32_t j = 0; j <= n; j++) {
      mpz_set_ui(cu_coordinates[j].get_data(), 0);
    }

    /* Compute the square norm of the shortest non-zero vector in the reduced
     * lattice basis A. */
    mpfr_t square_radius;
    mpfr_init2(square_radius, precision);
    mpfr_set_ui(square_radius, 0, MPFR_RNDN);

    for (uint32_t j = 0; j <= n; j++) {
      mpfr_set_z(tmp, A[0][j].get_data(), MPFR_RNDN);
      mpfr_mul(tmp, tmp, tmp, MPFR_RNDN);
      mpfr_add(square_radius, square_radius, tmp, MPFR_RNDN);
    }

    /* To attempt to ensure that we do not miss the shortest vector in the
     * lattice, reduce the square norm by #INITIAL_RADIUS_REDUCTION_FACTOR. */
    mpfr_div_ui(
      square_radius,
      square_radius,
      INITIAL_RADIUS_REDUCTION_FACTOR,
      MPFR_RNDN);

    /* Start a timer. */
    Timer timer;
    timer_start(&timer);

    /* Execute the enumeration. */
    (*status_d) = LATTICE_STATUS_NOT_RECOVERED;

    mpz_t min_dr;
    mpz_init(min_dr);
    mpz_set_ui(min_dr, 0); /* min_dr = 0 */

    mpz_t max_dr;
    mpz_init(max_dr);
    mpz_set_ui(max_dr, 0);
    mpz_setbit(max_dr, 2 * parameters->m); /* max_dr = 2^(2m) > d * r */

    mpz_t dr;
    mpz_init(dr);
    mpz_mul(dr, parameters->d, parameters->r);

    for (uint32_t t = 0; t < MAX_RADIUS_DOUBLINGS; t++) {
      lattice_enumerate_inner(
        status_d,
        square_radius,
        cu_coordinates,
        cv_coordinates,
        A,
        G_square_norms,
        M,
        n + 1, /* = k */
        n,
        min_dr,
        max_dr,
        dr,
        parameters->m, /* = m */
        precision,
        FALSE, /* = detect_smooth */
        timeout,
        &timer);

      /* Note that we set detect_smooth to FALSE above. If r is partially very
       * smooth, there may be an artificially short vector in the lattice,
       * leading us to enumerate very many vectors. To handle this, however, we
       * reduce the upper interval for d, see the above code. The detect_smooth
       * flag is only relevant when enumerating for r. */

      if (LATTICE_STATUS_NOT_RECOVERED != (*status_d)) {
        /* A solution was found or there was a timeout. */
        break;
      }

      /* Double the radius explored in the next iteration. */
      mpfr_mul_ui(square_radius, square_radius, 4, MPFR_RNDN);
    }

    /* Clear memory. */
    mpz_clear(dr);
    mpz_clear(min_dr);
    mpz_clear(max_dr);

    cu_coordinates.clear();

    mpfr_clear(square_radius);
  }

  /* Clear memory. */
  G_square_norms.clear();
  cv_coordinates.clear();
  v.clear();

  mpz_clear(candidate_d);
  mpz_clear(target_d);
  mpfr_clear(tmp);

  mpz_clear(pow2m);
  mpz_clear(rpow2lm);
  mpz_clear(rpow2lsigma);
}

void lattice_enumerate_reduced_basis_for_r(
  Lattice_Status_Recovery * const status_r,
  const ZZ_mat<mpz_t> &A,
  const FP_mat<mpfr_t> &G,
  const FP_mat<mpfr_t> &M,
  const uint32_t n,
  const Parameters * const parameters,
  const uint32_t precision,
  const bool detect_smooth_r,
  const uint64_t timeout)
{
  /* Check dimensions. */
  if ((A.get_rows() != (int)(n + 1)) || (A.get_cols() != (int)(n + 1))) {
    critical("lattice_enumerate_reduced_basis_for_r(): "
      "Incorrect dimensions for the matrix A.");
  }

  if ((G.get_rows() != (int)(n + 1)) || (G.get_cols() != (int)(n + 1))) {
    critical("lattice_enumerate_reduced_basis_for_r(): "
      "Incorrect dimensions for the matrix G.");
  }

  if ((M.get_rows() != (int)(n + 1)) || (M.get_cols() != (int)(n + 1))) {
    critical("lattice_enumerate_reduced_basis_for_r(): "
      "Incorrect dimensions for the matrix M.");
  }

  /* Setup variables. */
  mpfr_t tmp;
  mpfr_init2(tmp, precision);

  /* Compute the square norms of the rows of the G matrix. */
  vector<FP_NR<mpfr_t>> G_square_norms(n + 1);

  for (uint32_t i = 0; i <= n; i++) {
    mpfr_set_prec(G_square_norms[i].get_data(), precision);
    mpfr_set_ui(G_square_norms[i].get_data(), 0, MPFR_RNDN);

    for (uint32_t j = 0; j <= n; j++) {
      mpfr_mul(tmp, G[i][j].get_data(), G[i][j].get_data(), MPFR_RNDN);

      mpfr_add(
        G_square_norms[i].get_data(),
        G_square_norms[i].get_data(),
        tmp,
        MPFR_RNDN);
    }
  }

  /* Setup c_v for the target vector v = 0. */
  vector<FP_NR<mpfr_t>> cv_coordinates(n + 1);
  for (uint32_t j = 0; j <= n; j++) {
    mpfr_set_ui(cv_coordinates[j].get_data(), 0, MPFR_RNDN);
  }

  /* Setup c_u for the zero coordinate vector. */
  vector<Z_NR<mpz_t>> cu_coordinates(n + 1);
  for (uint32_t j = 0; j <= n; j++) {
    mpz_set_ui(cu_coordinates[j].get_data(), 0);
  }

  /* Compute the square norm of the shortest non-zero vector in the reduced
   * lattice basis A. */
  mpfr_t square_radius;
  mpfr_init2(square_radius, precision);
  mpfr_set_ui(square_radius, 0, MPFR_RNDN);

  for (uint32_t j = 0; j <= n; j++) {
    mpfr_set_z(tmp, A[0][j].get_data(), MPFR_RNDN);
    mpfr_mul(tmp, tmp, tmp, MPFR_RNDN);
    mpfr_add(square_radius, square_radius, tmp, MPFR_RNDN);
  }

  /* To attempt to ensure that we do not miss the shortest vector in the
   * lattice, reduce the square norm by #INITIAL_RADIUS_REDUCTION_FACTOR. */
    mpfr_div_ui(
      square_radius,
      square_radius,
      INITIAL_RADIUS_REDUCTION_FACTOR,
      MPFR_RNDN);

  /* Start a timer. */
  Timer timer;
  timer_start(&timer);

  /* Execute the enumeration. */
  (*status_r) = LATTICE_STATUS_NOT_RECOVERED;

  mpz_t min_r;
  mpz_init(min_r);
  mpz_set_ui(min_r, 0); /* min_r = 0 */

  mpz_t max_r;
  mpz_init(max_r);
  mpz_set_ui(max_r, 0);
  mpz_setbit(max_r, parameters->m); /* max_r = 2^m */

  for (uint32_t t = 0; t < MAX_RADIUS_DOUBLINGS; t++) {
    lattice_enumerate_inner(
      status_r,
      square_radius,
      cu_coordinates,
      cv_coordinates,
      A,
      G_square_norms,
      M,
      n + 1, /* = k */
      n,
      min_r,
      max_r,
      parameters->r, /* = r */
      parameters->m, /* = m */
      precision,
      detect_smooth_r,
      timeout,
      &timer);

    if (LATTICE_STATUS_NOT_RECOVERED != (*status_r)) {
      break;
    }

    /* Double the radius explored in the next iteration. */
    mpfr_mul_ui(square_radius, square_radius, 4, MPFR_RNDN);
  }

  /* Clear memory. */
  mpz_clear(min_r);
  mpz_clear(max_r);

  cv_coordinates.clear();
  cu_coordinates.clear();

  G_square_norms.clear();

  mpfr_clear(square_radius);
  mpfr_clear(tmp);
}

void lattice_enumerate_for_d(
  Lattice_Status_Recovery * const status_d,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const bool detect_smooth_r,
  const uint64_t timeout)
{
  /* Setup the A matrix. */
  ZZ_mat<mpz_t> A(n + 1, n + 1);

  /* Setup the Gram-Schmidt matrices. */
  FP_mat<mpfr_t> G(n + 1, n + 1);
  FP_mat<mpfr_t> M(n + 1, n + 1);

  while (TRUE) {
    /* Reduce the A matrix. */
    lattice_compute_reduced_basis(
      A,
      js,
      n,
      parameters,
      (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) ?
        REDUCTION_ALGORITHM_LLL : algorithm,
      parameters->m); /* red_precision */

    /* Compute the Gram-Schmidt matrices. */
    gram_schmidt_orthogonalization(M, G, A, n, precision);

    /* Solve for d. */
    lattice_solve_reduced_basis_for_d(
      status_d,
      A,
      G,
      ks,
      n,
      parameters,
      precision,
      detect_smooth_r);

    if (LATTICE_STATUS_NOT_RECOVERED != (*status_d)) {
      A.clear();
      G.clear();
      M.clear();

      return;
    }

    if (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) {
      algorithm = REDUCTION_ALGORITHM_BKZ;
    } else {
      break;
    }
  }

  /* Enumerate for d. */
  lattice_enumerate_reduced_basis_for_d(
    status_d,
    A,
    G,
    M,
    ks,
    n,
    parameters,
    precision,
    detect_smooth_r,
    timeout);

  /* Clear memory. */
  A.clear();
  G.clear();
  M.clear();
}

void lattice_enumerate_for_d_given_r(
  Lattice_Status_Recovery * const status_d,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Diagonal_Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const uint64_t timeout)
{
  /* Setup the A matrix. */
  ZZ_mat<mpz_t> A(n + 1, n + 1);

  /* Setup the Gram-Schmidt matrices. */
  FP_mat<mpfr_t> G(n + 1, n + 1);
  FP_mat<mpfr_t> M(n + 1, n + 1);

  while (TRUE) {
    /* Reduce the A matrix. */
    lattice_compute_reduced_diagonal_basis(
      A,
      js,
      n,
      parameters,
      (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) ?
        REDUCTION_ALGORITHM_LLL : algorithm,
      2 * parameters->m); /* red_precision */

    /* Compute the Gram-Schmidt matrices. */
    gram_schmidt_orthogonalization(M, G, A, n, precision);

    /* Solve for d. */
    lattice_solve_reduced_basis_for_d_given_r(
      status_d,
      A,
      G,
      ks,
      n,
      parameters,
      precision);

    if (LATTICE_STATUS_NOT_RECOVERED != (*status_d)) {
      A.clear();
      G.clear();
      M.clear();

      return;
    }

    if (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) {
      algorithm = REDUCTION_ALGORITHM_BKZ;
    } else {
      break;
    }
  }

  /* Enumerate for d. */
  lattice_enumerate_reduced_basis_for_d_given_r(
    status_d,
    A,
    G,
    M,
    ks,
    n,
    parameters,
    precision,
    timeout);

  /* Clear memory. */
  A.clear();
  G.clear();
  M.clear();
}

void lattice_enumerate_for_r(
  Lattice_Status_Recovery * const status_r,
  const mpz_t * const js,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const bool detect_smooth_r,
  const uint64_t timeout)
{
  /* Setup the A matrix. */
  ZZ_mat<mpz_t> A(n + 1, n + 1);

  while (TRUE) {
    /* Reduce the A matrix. */
    lattice_compute_reduced_basis(
      A,
      js,
      n,
      parameters,
      (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) ?
        REDUCTION_ALGORITHM_LLL : algorithm,
      parameters->m); /* red_precision */

    /* Solve for r. */
    lattice_solve_reduced_basis_for_r(
      status_r,
      A,
      n,
      parameters,
      detect_smooth_r);

    if (LATTICE_STATUS_NOT_RECOVERED != (*status_r)) {
      A.clear();

      return;
    }

    if (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) {
      algorithm = REDUCTION_ALGORITHM_BKZ;
    } else {
      break;
    }
  }

  /* Compute the Gram-Schmidt matrices. */
  FP_mat<mpfr_t> G(n + 1, n + 1);
  FP_mat<mpfr_t> M(n + 1, n + 1);

  gram_schmidt_orthogonalization(M, G, A, n, precision);

  /* Enumerate for r. */
  lattice_enumerate_reduced_basis_for_r(
    status_r,
    A,
    G,
    M,
    n,
    parameters,
    precision,
    detect_smooth_r,
    timeout);

  /* Clear memory. */
  A.clear();
  G.clear();
  M.clear();
}

void lattice_enumerate_for_d_r(
  Lattice_Status_Recovery * const status_d,
  Lattice_Status_Recovery * const status_r,
  const mpz_t * const js,
  const mpz_t * const ks,
  const uint32_t n,
  const Parameters * const parameters,
  Lattice_Reduction_Algorithm algorithm,
  const uint32_t precision,
  const bool detect_smooth_r,
  const uint64_t timeout)
{
  /* Setup the A matrix. */
  ZZ_mat<mpz_t> A(n + 1, n + 1);

  /* Setup the Gram-Schmidt matrices. */
  FP_mat<mpfr_t> G(n + 1, n + 1);
  FP_mat<mpfr_t> M(n + 1, n + 1);

  /* Initially set the statuses to not recovered. */
  (*status_d) = LATTICE_STATUS_NOT_RECOVERED;
  (*status_r) = LATTICE_STATUS_NOT_RECOVERED;

  bool gso_updated;

  while (TRUE) {
    /* Reduce the A matrix. */
    lattice_compute_reduced_basis(
      A,
      js,
      n,
      parameters,
      (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) ?
        REDUCTION_ALGORITHM_LLL : algorithm,
      parameters->m); /* red_precision */

    gso_updated = FALSE;

    if (LATTICE_STATUS_NOT_RECOVERED == (*status_r)) {
      /* Solve for r. */
      lattice_solve_reduced_basis_for_r(
        status_r,
        A,
        n,
        parameters,
        detect_smooth_r);
    }

    if (LATTICE_STATUS_NOT_RECOVERED == (*status_d)) {
      /* Compute the Gram-Schmidt matrices. */
      gram_schmidt_orthogonalization(M, G, A, n, precision);

      gso_updated = TRUE;

      /* Solve for d. */
      lattice_solve_reduced_basis_for_d(
        status_d,
        A,
        G,
        ks,
        n,
        parameters,
        precision,
        detect_smooth_r);
    }

    if ((LATTICE_STATUS_NOT_RECOVERED != (*status_d)) &&
        (LATTICE_STATUS_NOT_RECOVERED != (*status_r)))
    {
      A.clear();
      G.clear();
      M.clear();

      return;
    }

    if (REDUCTION_ALGORITHM_LLL_BKZ == algorithm) {
      algorithm = REDUCTION_ALGORITHM_BKZ;
    } else {
      break;
    }
  }

  if (TRUE != gso_updated) {
    /* Compute the Gram-Schmidt matrices. */
    gram_schmidt_orthogonalization(M, G, A, n, precision);

    gso_updated = TRUE;
  }

  if (LATTICE_STATUS_NOT_RECOVERED == (*status_d)) {
    /* Enumerate for d. */
    lattice_enumerate_reduced_basis_for_d(
      status_d,
      A,
      G,
      M,
      ks,
      n,
      parameters,
      precision,
      detect_smooth_r,
      timeout);
  }

  if (LATTICE_STATUS_NOT_RECOVERED == (*status_r)) {
    /* Enumerate for r. */
    lattice_enumerate_reduced_basis_for_r(
      status_r,
      A,
      G,
      M,
      n,
      parameters,
      precision,
      detect_smooth_r,
      timeout);
  }

  /* Clear memory. */
  A.clear();
  G.clear();
  M.clear();
}
