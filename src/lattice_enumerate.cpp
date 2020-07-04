/*!
 * \file    lattice_enumerate.cpp
 * \ingroup lattice_enumerate
 *
 * \brief   The definition of functions for solving for d and r using Kannan's
 *          enumeration algorithm and lattice basis reduction techniques.
 */

#include "lattice_enumerate.h"

#include "lattice.h"
#include "lattice_gso.h"
#include "lattice_babai.h"
#include "lattice_algebra.h"
#include "lattice_reduce.h"
#include "lattice_solve.h"

#include "parameters.h"
#include "sample.h"
#include "errors.h"
#include "timer.h"
#include "math.h"

#include <fplll/fplll.h>

#include <gmp.h>
#include <mpfr.h>

#include <vector>

#include <stdint.h>

using namespace std;
using namespace fplll;

/*!
 * \brief   The factor by which the estimated square norm of the shortest 
 *          vector in the lattice should be reduced.
 *
 * This reduction is performed so as to ensure that we do not accidently miss
 * the shortest vector in the lattice.
 */
#define INITIAL_RADIUS_REDUCTION_FACTOR 1024

/*!
 * \brief   The maximum number of radius doublings.
 * 
 * Enumeration is performed in balls of increasingly greater radius. The initial
 * radius is given by the reduced estimated norm of the shortest vector in the 
 * lattice. This norm is then doubled repeatedly.
 */
#define MAX_RADIUS_DOUBLINGS 128

/*!
 * \brief   Executes the pi projection test proposed by Kannan.
 *
 * For further details, please see [1] and Appendix C in [2].
 * 
 * [1] Kannan, R.: Improved algorithms for integer programming and related 
 * lattice problems. In: Proceedings of the 15th Symposium on the Theory of 
 * Computing. STOC 1983. ACM Press, pp. 99—108 (1983).
 * 
 * [2] Ekerå, M.: On post-processing in the quantum algorithm for computing 
 * short discrete logarithms. In: IACR ePrint Archive, 2017/1122, 1st revision.
 * 
 * \internal
 *
 * \param[in] square_radius     The square radius.
 * \param[in] cu_coordinates    The coordinates of the u vector.
 * \param[in] cv_coordinates    The coordinates of the v vector.
 * \param[in] M                 The M matrix.
 * \param[in] G_square_norms    The square norms of the rows of the G matrix.
 * \param[in]  k                The index of the component currently processed.
 *                              Specifies the depth in the enumeration tree.
 * \param[in] n                 The number of runs used to setup the A matrix
 *                              from which the G and M matrices are computed.
 * \param[in] precision         The floating point precision.
 */
static bool pi_projection_test(
  const mpfr_t square_radius,
  const vector<Z_NR<mpz_t>> &cu_coordinates,
  const vector<FP_NR<mpfr_t>> &cv_coordinates,
  const FP_mat<mpfr_t> &M,
  const vector<FP_NR<mpfr_t>> &G_square_norms,
  const uint32_t k,
  const uint32_t n,
  const uint32_t precision)
{
  bool result = TRUE;

  vector<FP_NR<mpfr_t>> u(n + 1);

  mpfr_t square_norm;
  mpfr_init2(square_norm, precision);

  mpfr_t tmp;
  mpfr_init2(tmp, precision);

  mpfr_t tmp2;
  mpfr_init2(tmp2, precision);
  
  /* Grow the radius by 1 percent to be on the safe side. */
  mpfr_t grown_square_radius;
  mpfr_init2(grown_square_radius, precision);
  mpfr_mul_d(grown_square_radius, square_radius, (double)1.01f, MPFR_RNDN);

  for (uint32_t j = 0; j <= n; j++) {
    mpfr_set_prec(u[j].get_data(), precision);

    /* Let u = c_u - c_v. */
    mpfr_sub_z(
      u[j].get_data(),
      cv_coordinates[j].get_data(),
      cu_coordinates[j].get_data(),
      MPFR_RNDN);
    mpfr_neg(u[j].get_data(), u[j].get_data(), MPFR_RNDN);
  }

  /* Compute the norms. */

  /* Note that the parameter k is indexed from one, whilst the indices i and j
   * are indexed from zero; hence we subtract one from k when used. */
  mpfr_set_ui(square_norm, 0, MPFR_RNDN);

  for (uint32_t j = k - 1; j <= n; j++) { /* k is indexed from one */
    mpfr_set(tmp, u[j].get_data(), MPFR_RNDN);

    for (uint32_t i = j + 1; i <= n; i++) {
      mpfr_mul(tmp2, M[i][j].get_data(), u[i].get_data(), MPFR_RNDN);
      mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
    }

    mpfr_mul(tmp, tmp, tmp, MPFR_RNDN);
    mpfr_mul(tmp, tmp, G_square_norms[j].get_data(), MPFR_RNDN);

    mpfr_add(square_norm, square_norm, tmp, MPFR_RNDN);

    if (mpfr_cmp(square_norm, grown_square_radius) > 0) {
      result = FALSE;
      break;
    }
  }

  /* Clean up memory. */
  mpfr_clear(grown_square_radius);
  mpfr_clear(square_norm);
  mpfr_clear(tmp);
  mpfr_clear(tmp2);

  u.clear();

  /* Signal result. */
  return result;
}

/*!
 * \brief   An internal function called recursively to perform lattice
 *          enumeration using ideas from Kannan.
 *
 * This function implements lattice enumeration. It is implemented naïvely 
 * following Kannan [1] as described in appendix C to [2].
 * 
 * [1] Kannan, R.: Improved algorithms for integer programming and related 
 * lattice problems. In: Proceedings of the 15th Symposium on the Theory of 
 * Computing. STOC 1983. ACM Press, pp. 99—108 (1983).
 * 
 * [2] Ekerå, M.: On post-processing in the quantum algorithm for computing 
 * short discrete logarithms. In: IACR ePrint Archive, 2017/1122, 1st revision.
 *
 * Note however that vectors and matrices are indexed from one in the paper,
 * but indexed from zero in C and C++. Note furthermore that we do not collect
 * a set of potential solutions; instead we test all potential solutions as
 * they are found until a correct solution is found.
 *
 * This is why the solution must be explicitly passed to this function.
 *
 * To speed up enumeration, we make some further optimizations:
 *
 * (i)  We reduce the  number of pi projection tests that we perform in each
 *      component, primarily to speed up enumeration of low-dimensional
 *      lattices. We seek for the bounds on the index in each component first
 *      and then explore the whole range in this component without testing.
 *
 * (ii) We use that the last component is d or r and hence on 0 <= d, r < 2^m.
 *
 * Note that this a fairly straightforward implementation of lattice enumeration
 * that has proven sufficient for these purposes. There are more efficient
 * methods for enumerating lattices, and it is possible to implement this
 * function more efficiently. However, it has not proved necessary.
 *
 * \internal
 *
 * \param[out] status           Used to report status of the enumeration.
 * \param[in]  square_radius    The square radius of the hypersphere to
 *                              enumerate.
 * \param[in]  cu_coordinates   The coordinates of the u vector.
 * \param[in]  cv_coordinates   The coordinates of the v vector.
 * \param[in]  A                The A matrix.
 * \param[in]  G_square_norms   The square norms of the rows of the G matrix.
 * \param[in]  M                The M matrix.
 * \param[in]  k                The index of the component currently processed.
 *                              Specifies the depth in the enumeration tree.
 * \param[in]  m                An integer m such that 0 <= d, r < 2^m.
 * \param[in]  n                The number of runs n used to setup A so that A,
 *                              G and M are all (n+1) x (n+1) matrices.
 * \param[in]  d                The correct solution. Used only to test if the
 *                              candidate solutions found are correct.
 * \param[in]  precision        The floating point precision.
 * \param[in]  timer            A timer used to measure how much time has
 *                              elapsed since the first call this this function.
 *                              Must be initialized and started by the caller.
 * \param[in]  timeout          The timeout after which the enumeration will be
 *                              aborted. Set to zero to disable the timeout.
 */
static void lattice_enumerate_inner(
  Lattice_Status_Recovery * const status,
  const mpfr_t square_radius,
  const vector<Z_NR<mpz_t>> &cu_coordinates,
  const vector<FP_NR<mpfr_t>> &cv_coordinates,
  const ZZ_mat<mpz_t> &A,
  const vector<FP_NR<mpfr_t>> &G_square_norms,
  const FP_mat<mpfr_t> &M,
  const uint32_t k,
  const uint32_t m,
  const uint32_t n,
  const mpz_t d,
  const uint32_t precision,
  Timer * const timer,
  const uint64_t timeout)
{
  /* Check the timeout. */
  if ((0 != timeout) && (timer_stop(timer) / (1000 * 1000) > timeout)) {
    (*status) = LATTICE_STATUS_TIMEOUT;

    return;
  }

  /* Initially set the status to not recovered. */
  (*status) = LATTICE_STATUS_NOT_RECOVERED;

  /* If k is equal to zero we have found a potential solution. */
  if (0 == k) {
    /* Note: This code may in principle be removed. It is no longer needed. */

    /* Initialize variables. */
    mpz_t tmp_z;
    mpz_init(tmp_z);

    mpz_t last_component_u;
    mpz_init(last_component_u);

    /* Test the last component. */
    for (uint32_t i = 0; i <= n; i++) {
      mpz_mul(tmp_z, A[i][n].get_data(), cu_coordinates[i].get_data());
      mpz_add(last_component_u, last_component_u, tmp_z);
    }

    mpz_abs(last_component_u, last_component_u);

    const bool result = (mpz_cmp(last_component_u, d) == 0);

    /* Clean up memory. */
    mpz_clear(tmp_z);
    mpz_clear(last_component_u);

    if (result) {
      /* Set the status to recovered. */
      (*status) = LATTICE_STATUS_RECOVERED_ENUMERATION;
    }

    /* Stop. */
    return;
  }

  /* Initialize variables. */
  mpz_t uhat;
  mpz_init(uhat);

  {
    mpfr_t tmp;
    mpfr_init2(tmp, precision);

    mpfr_t tmp2;
    mpfr_init2(tmp2, precision);

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

    /* Clean up memory. */
    mpfr_clear(tmp);
    mpfr_clear(tmp2);
  }

  /* Setup the new coordinate vector. */
  vector<Z_NR<mpz_t>> new_cu_coordinates(n + 1);

  for (uint32_t j = 0; j <= n; j++) {
    mpz_set(new_cu_coordinates[j].get_data(), cu_coordinates[j].get_data());
  }

  mpz_set(new_cu_coordinates[k - 1].get_data(), uhat);

  /* Find the maximum and minimum value in the k:th position. */
  mpz_t uhat_max;
  mpz_init(uhat_max);

  mpz_t uhat_min;
  mpz_init(uhat_min);

  {
    /* The maximum in the k:th component is = uhat + max_cu. Set max_cu to
     * zero initially, then to one, and then double max_cu in every iteration
     * thereafter until we no longer find ourselves inside the radius of the
     * hypersphere. This yields a coarse estimate that we then refine below. */
    mpz_t max_cu;
    mpz_init_set_ui(max_cu, 0);

    while (TRUE) {
      mpz_add(new_cu_coordinates[k - 1].get_data(), uhat, max_cu);

      bool result = pi_projection_test(
          square_radius,
          new_cu_coordinates,
          cv_coordinates,
          M,
          G_square_norms,
          k,
          n,
          precision);
      if (!result) {
        break;
      }

      if (mpz_cmp_ui(max_cu, 0) == 0) {
        mpz_set_ui(max_cu, 1);
      } else {
        mpz_mul_ui(max_cu, max_cu, 2);
      }
    }

    if (mpz_cmp_ui(max_cu, 0) == 0) {
      /* Abort the run. */

      /* Clean up memory. */
      new_cu_coordinates.clear();

      mpz_clear(uhat);
      mpz_clear(uhat_min);
      mpz_clear(uhat_max);
      mpz_clear(max_cu);

      return;
    } else if (mpz_cmp_ui(max_cu, 2) <= 0) {
      /* We need not do anything; as we explicitly test max_cu in {0, 1, 2} we
       * are guaranteed that max_cu is a tight bound. */
    } else {
      /* Refine the search by using that we know that the minimum is acceptable,
       * whilst the maximum is unacceptable. Select a pivot element inbetween
       * the two values and test it to find a new minimum and maximum. Repeat
       * this process recursively until a tight bound is established. */

      /* Initially set min_cu = max_cu / 2. */
      mpz_t min_cu;
      mpz_init(min_cu);
      mpz_div_ui(min_cu, max_cu, 2);

      mpz_t pivot_cu;
      mpz_init(pivot_cu);

      while (TRUE) {
        mpz_sub(pivot_cu, max_cu, min_cu);
        if (mpz_cmp_ui(pivot_cu, 2) < 0) {
          break;
        }

        mpz_add(pivot_cu, max_cu, min_cu);
        mpz_div_ui(pivot_cu, pivot_cu, 2);

        mpz_add(new_cu_coordinates[k - 1].get_data(), uhat, pivot_cu);

        bool result = pi_projection_test(
            square_radius,
            new_cu_coordinates,
            cv_coordinates,
            M,
            G_square_norms,
            k,
            n,
            precision);
        if (result) {
          mpz_set(min_cu, pivot_cu);
        } else {
          mpz_set(max_cu, pivot_cu);
        }
      }

      /* Clean up memory. */
      mpz_clear(pivot_cu);
      mpz_clear(min_cu);
    }

    /* Export the result. */
    mpz_add(uhat_max, uhat, max_cu);

    /* Clean up memory. */
    mpz_clear(max_cu);
  }

  {
   /* The minimum in the k:th component is min = uhat + min_cu. Set min_cu to
    * -1 initially, and then double max_cu in every iteration thereafter until
    * we no longer find ourselves inside the radius of the hypersphere. This
    * yields a coarse estimate that we then refine below. */

    mpz_t min_cu;
    mpz_init_set_si(min_cu, -1);

    while (TRUE) {
      mpz_add(new_cu_coordinates[k - 1].get_data(), uhat, min_cu);

      bool result = pi_projection_test(
          square_radius,
          new_cu_coordinates,
          cv_coordinates,
          M,
          G_square_norms,
          k,
          n,
          precision);
      if (!result) {
        break;
      }

      mpz_mul_ui(min_cu, min_cu, 2);
    }

    if (mpz_cmp_si(min_cu, -2) >= 0) {
      /* We need not do anything; as we explicitly test min_cu in {-1, -2}, and
       * previously min_cu = 0, we know that min_cu is a tight bound. */
    } else {
      /* Refine the search by using that we know that the maximum is acceptable,
       * whilst the minimum is unacceptable. Select a pivot element inbetween
       * the two values and test it to find a new minimum and maximum. Repeat
       * this process recursively until a tight bound is established. */

      /* Initially set max_cu = min_cu / 2. */
      mpz_t max_cu;
      mpz_init(max_cu);
      mpz_div_ui(max_cu, min_cu, 2);

      mpz_t pivot_cu;
      mpz_init(pivot_cu);

      while (TRUE) {
        mpz_sub(pivot_cu, min_cu, max_cu);
        if (mpz_cmp_si(pivot_cu, -2) > 0) {
          break;
        }

        mpz_add(pivot_cu, max_cu, min_cu);
        mpz_div_ui(pivot_cu, pivot_cu, 2);

        mpz_add(new_cu_coordinates[k - 1].get_data(), uhat, pivot_cu);

        bool result = pi_projection_test(
            square_radius,
            new_cu_coordinates,
            cv_coordinates,
            M,
            G_square_norms,
            k,
            n,
            precision);
        if (result) {
          mpz_set(max_cu, pivot_cu);
        } else {
          mpz_set(min_cu, pivot_cu);
        }
      }

      /* Clean up memory. */
      mpz_clear(pivot_cu);
      mpz_clear(max_cu);
    }

    /* Export the result. */
    mpz_add(uhat_min, uhat, min_cu);

    /* Clean up memory. */
    mpz_clear(min_cu);
  }

  /* Begin: Optimization for d and r on 0 <= d, r < 2^m --------------------- */
  /* We use _d suffixes below, however the situation is analogous in d and r. */

  if ((1 == k) && TRUE) {
    /* Initialize variables. */
    mpz_t last_component_u;
    mpz_init(last_component_u);

    /* The interval is now (uhat_min, uhat_max) around uhat. */
    mpz_set(new_cu_coordinates[0].get_data(), uhat);

    {
      mpz_t tmp;
      mpz_init(tmp);

      for (uint32_t i = 0; i <= n; i++) {
        mpz_mul(
          tmp,
          A[i][n].get_data(),
          new_cu_coordinates[i].get_data());
        mpz_add(
          last_component_u,
          last_component_u,
          tmp);
      }

      /* Clean memory. */
      mpz_clear(tmp);
    }

    mpz_t min_d;
    mpz_init(min_d);
    mpz_set_ui(min_d, 0);

    mpz_t max_d;
    mpz_init(max_d);
    mpz_setbit(max_d, m);

    const int sign = (mpz_sgn(A[0][n].get_data()) < 0) ? -1 : 1;

    mpz_t candidate_d;
    mpz_init_set(candidate_d, last_component_u);

    mpz_t increment_d;
    mpz_init(increment_d);
    mpz_abs(increment_d, A[0][n].get_data());

    /* Move ahead if necessary. */
    if (mpz_cmp(candidate_d, min_d) < 0) {
      mpz_t tmp;
      mpz_init(tmp);

      mpz_sub(tmp, min_d, candidate_d);
      mpz_add(tmp, tmp, increment_d);
      mpz_sub_ui(tmp, tmp, 1);
      mpz_div(tmp, tmp, increment_d);

      if (sign < 1) {
        mpz_sub(
          new_cu_coordinates[k - 1].get_data(),
          new_cu_coordinates[k - 1].get_data(),
          tmp);
      } else {
        mpz_add(
          new_cu_coordinates[k - 1].get_data(),
          new_cu_coordinates[k - 1].get_data(),
          tmp);
      }

      mpz_mul(tmp, tmp, increment_d);
      mpz_add(candidate_d, candidate_d, tmp);

      mpz_clear(tmp);
    }

    /* Iterate. */
    while (TRUE) {
      const int result_min_d = mpz_cmp(candidate_d, min_d);
      const int result_max_d = mpz_cmp(candidate_d, max_d);

      if (result_max_d > 0) {
        break;
      }

      if ((result_min_d >= 0) && (result_max_d <= 0)) {
        if (mpz_cmp(candidate_d, d) == 0) {
          /* Set the status to recovered. */
          (*status) = LATTICE_STATUS_RECOVERED_ENUMERATION;

          /* Clean up memory. */
          mpz_clear(max_d);
          mpz_clear(min_d);
          mpz_clear(candidate_d);
          mpz_clear(increment_d);

          new_cu_coordinates.clear();

          mpz_clear(uhat);
          mpz_clear(uhat_min);
          mpz_clear(uhat_max);

          /* Stop. */
          return;
        }
      }

      /* Increment the coordinate. */
      if (sign < 0) {
        mpz_sub_ui(
          new_cu_coordinates[k - 1].get_data(),
          new_cu_coordinates[k - 1].get_data(),
          1);
      } else {
        mpz_add_ui(
          new_cu_coordinates[k - 1].get_data(),
          new_cu_coordinates[k - 1].get_data(),
          1);
      }

      if (mpz_cmp(new_cu_coordinates[k - 1].get_data(), uhat_max) >= 0) {
        break;
      }

      /* Increment the candidate d. */
      mpz_add(candidate_d, candidate_d, increment_d);
    }

    /* Reset the counters. */
    mpz_set(candidate_d, last_component_u);
    mpz_set(new_cu_coordinates[k - 1].get_data(), uhat);

    /* Move ahead if necessary. */
    if (mpz_cmp(candidate_d, max_d) > 0) {
      mpz_t tmp;
      mpz_init(tmp);

      mpz_sub(tmp, candidate_d, max_d);
      mpz_sub(tmp, tmp, increment_d);
      mpz_add_ui(tmp, tmp, 1);
      mpz_div(tmp, tmp, increment_d);

      if (sign < 1) {
        mpz_add(
          new_cu_coordinates[k - 1].get_data(),
          new_cu_coordinates[k - 1].get_data(),
          tmp);
      } else {
        mpz_sub(
          new_cu_coordinates[k - 1].get_data(),
          new_cu_coordinates[k - 1].get_data(),
          tmp);
      }

      mpz_mul(tmp, tmp, increment_d);
      mpz_sub(candidate_d, candidate_d, tmp);

      mpz_clear(tmp);
    }

    /* Iterate. */
    while (TRUE) {
      /* Decrement the candidate d. */
      mpz_sub(candidate_d, candidate_d, increment_d);

      /* Decrement the coordinate. */
      if (sign < 0) {
        mpz_add_ui(
          new_cu_coordinates[k - 1].get_data(),
          new_cu_coordinates[k - 1].get_data(),
          1);
      } else {
        mpz_sub_ui(
          new_cu_coordinates[k - 1].get_data(),
          new_cu_coordinates[k - 1].get_data(),
          1);
      }

      if (mpz_cmp(new_cu_coordinates[k - 1].get_data(), uhat_min) <= 0) {
        break;
      }

      const int result_min_d = mpz_cmp(candidate_d, min_d);
      const int result_max_d = mpz_cmp(candidate_d, max_d);

      if (result_min_d < 0) {
        break;
      }

      if ((result_min_d >= 0) && (result_max_d <= 0)) {
        if (mpz_cmp(candidate_d, d) == 0) {
          /* Set the status to recovered. */
          (*status) = LATTICE_STATUS_RECOVERED_ENUMERATION;

          /* Clean up memory. */
          mpz_clear(max_d);
          mpz_clear(min_d);
          mpz_clear(candidate_d);
          mpz_clear(increment_d);

          new_cu_coordinates.clear();

          mpz_clear(uhat);
          mpz_clear(uhat_min);
          mpz_clear(uhat_max);

          /* Stop. */
          return;
        }
      }
    }

    /* Clean up memory. */
    mpz_clear(max_d);
    mpz_clear(min_d);
    mpz_clear(candidate_d);
    mpz_clear(increment_d);

    new_cu_coordinates.clear();

    mpz_clear(uhat);
    mpz_clear(uhat_min);
    mpz_clear(uhat_max);

    return; /* No solution found. */
  }

  /* End: Optimization for d and r. ----------------------------------------- */

  /* The interval is now (uhat_min, .., uhat_max). */
  mpz_t counter;
  mpz_init_set_ui(counter, 0);

  while (TRUE) {
    mpz_add(new_cu_coordinates[k - 1].get_data(), uhat, counter);

    bool exceeds_max =
      (mpz_cmp(new_cu_coordinates[k - 1].get_data(), uhat_max) >= 0);

    if (!exceeds_max) {
      /* Recursively enumerate the lattice. */
      lattice_enumerate_inner(
        status,
        square_radius,
        new_cu_coordinates,
        cv_coordinates,
        A,
        G_square_norms,
        M,
        k - 1,
        m,
        n,
        d,
        precision,
        timer,
        timeout);
      if (LATTICE_STATUS_NOT_RECOVERED != (*status)) {
        break; /* Recovery was successful or an error occurred. */
      }
    }

    mpz_sub(new_cu_coordinates[k - 1].get_data(), uhat, counter);

    bool exceeds_min =
      (mpz_cmp(new_cu_coordinates[k - 1].get_data(), uhat_min) <= 0);

    if ((!exceeds_min) && (mpz_cmp_ui(counter, 0) != 0)) {
      /* Recursively enumerate the lattice. */
      lattice_enumerate_inner(
        status,
        square_radius,
        new_cu_coordinates,
        cv_coordinates,
        A,
        G_square_norms,
        M,
        k - 1,
        m,
        n,
        d,
        precision,
        timer,
        timeout);
      if (LATTICE_STATUS_NOT_RECOVERED != (*status)) {
        break; /* Recovery was successful or an error occurred. */
      }
    }

    if (exceeds_min && exceeds_max) {
      break;
    }

    mpz_add_ui(counter, counter, 1);
  }

  mpz_clear(counter);

  /* Clean up memory. */
  new_cu_coordinates.clear();

  mpz_clear(uhat);
  mpz_clear(uhat_min);
  mpz_clear(uhat_max);
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

  /* Setup the starting vector:
   * v = ({-2^m k_1}_{2^{m+l}}, .., {-2^m k_n}_{2^{m+l}}, 0). */
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

    /* Initially set the radius to the square norm of the shortest vector. */
    mpfr_t square_radius;
    mpfr_init2(square_radius, precision);
    mpfr_set_ui(square_radius, 0, MPFR_RNDN);

    for (uint32_t j = 0; j <= n; j++) {
      mpfr_set_z(tmp, A[0][j].get_data(), MPFR_RNDN);
      mpfr_mul(tmp, tmp, tmp, MPFR_RNDN);
      mpfr_add(square_radius, square_radius, tmp, MPFR_RNDN);
    }

    /* To make sure we do not miss the vector, reduce radius by 2^10. */
    mpfr_div_ui(square_radius, square_radius, 1024, MPFR_RNDN);

    /* Start a timer. */
    Timer timer;
    timer_start(&timer);

    /* Execute the enumeration. */
    (*status_d) = LATTICE_STATUS_NOT_RECOVERED;

    for (uint32_t t = 0; t < 128; t++) {
      lattice_enumerate_inner(
        status_d,
        square_radius,
        cu_coordinates,
        cv_coordinates,
        A,
        G_square_norms,
        M,
        n + 1, /* = k */
        parameters->m,
        n,
        parameters->d,
        precision,
        &timer,
        timeout);

      if (LATTICE_STATUS_NOT_RECOVERED != (*status_d)) {
        break;
      }

      /* Double the radius explored in the next iteration. */
      mpfr_mul_ui(square_radius, square_radius, 4, MPFR_RNDN);
    }

    /* Clear memory. */
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

void lattice_enumerate_reduced_basis_for_r(
  Lattice_Status_Recovery * const status_r,
  const ZZ_mat<mpz_t> &A,
  const FP_mat<mpfr_t> &G,
  const FP_mat<mpfr_t> &M,
  const uint32_t n,
  const Parameters * const parameters,
  const uint32_t precision,
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

  /* Initially set the radius to the square norm of the shortest vector. */
  mpfr_t square_radius;
  mpfr_init2(square_radius, precision);
  mpfr_set_ui(square_radius, 0, MPFR_RNDN);

  for (uint32_t j = 0; j <= n; j++) {
    mpfr_set_z(tmp, A[0][j].get_data(), MPFR_RNDN);
    mpfr_mul(tmp, tmp, tmp, MPFR_RNDN);
    mpfr_add(square_radius, square_radius, tmp, MPFR_RNDN);
  }

  /* To make sure we do not miss the vector, reduce radius by 2^10. */
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
      parameters->m,
      n,
      parameters->r,
      precision,
      &timer,
      timeout);

    if (LATTICE_STATUS_NOT_RECOVERED != (*status_r)) {
      break;
    }

    /* Double the radius explored in the next iteration. */
    mpfr_mul_ui(square_radius, square_radius, 4, MPFR_RNDN);
  }

  /* Clear memory. */
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
  lattice_enumerate_reduced_basis_for_d(
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
      parameters);

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
        parameters);
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
        precision);
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
      timeout);
  }

  /* Clear memory. */
  A.clear();
  G.clear();
  M.clear();
}
