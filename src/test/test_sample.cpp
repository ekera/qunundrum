/*!
 * \file    test/test_sample.cpp
 * \ingroup unit_tests_sample_distribution
 *
 * \brief   The definition of unit tests for the \ref sample_distribution
 *          module.
 */

#include "test_sample.h"

#include "test_common.h"

#include "../common.h"
#include "../diagonal_parameters.h"
#include "../distribution.h"
#include "../errors.h"
#include "../lattice.h"
#include "../math.h"
#include "../parameters.h"
#include "../parameters_selection.h"
#include "../random.h"
#include "../sample.h"

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>
#include <stdio.h>

/*!
 * \brief   The maximum buffer size in bytes.
 */
#define MAX_BUFFER_SIZE 8192

void test_sample_approximate_alpha_from_region() {
  printf("Testing sample_approximate_alpha_from_region()...\n");

  mpfr_t alpha;
  mpfr_init2(alpha, PRECISION);

  Random_State random_state;
  random_init(&random_state);

  for (uint32_t i = 2048 - 30; i < 2048 + 30; i++) {
    for (uint32_t j = 0; j < 10; j++) {
      double min_log_alpha = 10 * i + j;     min_log_alpha /= 10;
      double max_log_alpha = 10 * i + j + 1; max_log_alpha /= 10;

      /* Positive. */
      sample_approximate_alpha_from_region(
        alpha, min_log_alpha, max_log_alpha, &random_state);

      if (1 != mpfr_sgn(alpha)) {
        critical("test_sample_approximate_alpha_from_region(): "
          "Failed to sample. Incorrect sign.");
      }

      mpfr_abs(alpha, alpha, MPFR_RNDN);
      mpfr_log2(alpha, alpha, MPFR_RNDN);

      double log_alpha = mpfr_get_d(alpha, MPFR_RNDN);

      if ((log_alpha < min_log_alpha) || (log_alpha > max_log_alpha)) {
        critical("test_sample_approximate_alpha_from_region(): "
          "Failed to sample. Incorrect magnitude.");
      }

      /* Negative. */
      sample_approximate_alpha_from_region(
        alpha, -min_log_alpha, -max_log_alpha, &random_state);

      if (-1 != mpfr_sgn(alpha)) {
        critical("test_sample_approximate_alpha_from_region(): "
          "Failed to sample. Incorrect sign.");
      }

      mpfr_abs(alpha, alpha, MPFR_RNDN);
      mpfr_log2(alpha, alpha, MPFR_RNDN);

      log_alpha = mpfr_get_d(alpha, MPFR_RNDN);

      if ((log_alpha < min_log_alpha) || (log_alpha > max_log_alpha)) {
        critical("test_sample_approximate_alpha_from_region(): "
          "Failed to sample. Incorrect magnitude.");
      }

      /* Note: Calling with incorrect order would trigger a critical error. */
    }
  }

  /* Clear memory. */
  mpfr_clear(alpha);

  random_close(&random_state);
}

void test_sample_alpha_from_region() {
  printf("Testing sample_alpha_from_region()...\n");

  mpz_t alpha;
  mpz_init(alpha);

  mpfr_t alpha_f;
  mpfr_init2(alpha_f, PRECISION);

  Random_State random_state;
  random_init(&random_state);

  for (uint32_t kappa = 0; kappa < 10; kappa++) {
    for (uint32_t i = 2048 - 30; i < 2048 + 30; i++) {
      for (uint32_t j = 0; j < 10; j++) {
        double min_log_alpha = 10 * i + j;     min_log_alpha /= 10;
        double max_log_alpha = 10 * i + j + 1; max_log_alpha /= 10;

        /* Positive. */
        sample_alpha_from_region(
          alpha, min_log_alpha, max_log_alpha, kappa, &random_state);

        if (mpz_sgn(alpha) < 1) {
          critical("test_sample_alpha_from_region(): "
            "Failed to sample. Incorrect sign.");
        }

        mpfr_set_z(alpha_f, alpha, MPFR_RNDN);

        mpfr_abs(alpha_f, alpha_f, MPFR_RNDN);
        mpfr_log2(alpha_f, alpha_f, MPFR_RNDN);

        double log_alpha = mpfr_get_d(alpha_f, MPFR_RNDN);

        if ((log_alpha < min_log_alpha) || (log_alpha > max_log_alpha)) {
          critical("test_sample_alpha_from_region(): "
            "Failed to sample. Incorrect magnitude.");
        }

        mpz_mod_ui(alpha, alpha, 1UL << kappa);
        if (0 != mpz_cmp_ui(alpha, 0)) {
          critical("test_sample_alpha_from_region(): "
            "Failed to sample. Not divisible by 2^kappa.");
        }

        /* Negative. */
        sample_alpha_from_region(
          alpha, -min_log_alpha, -max_log_alpha, kappa, &random_state);

        if (mpz_sgn(alpha) >= 0) {
          critical("test_sample_alpha_from_region(): "
            "Failed to sample. Incorrect sign.");
        }

        mpfr_set_z(alpha_f, alpha, MPFR_RNDN);

        mpfr_abs(alpha_f, alpha_f, MPFR_RNDN);
        mpfr_log2(alpha_f, alpha_f, MPFR_RNDN);

        log_alpha = mpfr_get_d(alpha_f, MPFR_RNDN);

        if ((log_alpha < min_log_alpha) || (log_alpha > max_log_alpha)) {
          critical("test_sample_alpha_from_region(): "
            "Failed to sample. Incorrect magnitude.");
        }

        mpz_mod_ui(alpha, alpha, 1UL << kappa);
        if (0 != mpz_cmp_ui(alpha, 0)) {
          critical("test_sample_alpha_from_region(): "
            "Failed to sample. Not divisible by 2^kappa.");
        }

        /* Note: Calling with incorrect order would trigger a critical error. */
      }
    }
  }

  /* Clear memory. */
  mpz_clear(alpha);
  mpfr_clear(alpha_f);

  random_close(&random_state);
}

void test_sample_j_from_alpha_r() {
  printf("Testing sample_j_from_alpha_r()...\n");

  const uint32_t m = 2048;
  const uint32_t l = m / 8;
  const uint32_t t = 30;

  Random_State random_state;
  random_init(&random_state);

  mpz_t d;
  mpz_init(d);

  mpz_t r;
  mpz_init(r);

  mpz_t alpha;
  mpz_init(alpha);

  mpz_t j;
  mpz_init(j);

  mpz_t tmp;
  mpz_init(tmp);

  mpz_t pow2ml;
  mpz_init_set_ui(pow2ml, 0);
  mpz_setbit(pow2ml, m + l);

  parameters_selection_random_d_or_r(r, m);
  mpz_set(d, r); /* by convention */

  Parameters parameters;
  parameters_init(&parameters);
  parameters_explicit_m_l(&parameters, d, r, m, l, t);

  for (uint32_t i = m - 30; i < m + 30; i++) {
    for (uint32_t i2 = 0; i2 < 10; i2++) {
      double min_log_alpha = 10 * i + i2;     min_log_alpha /= 10;
      double max_log_alpha = 10 * i + i2 + 1; max_log_alpha /= 10;

      /* Standard order. Positive. */
      sample_alpha_from_region(
        alpha, min_log_alpha, max_log_alpha, kappa(r), &random_state);
      sample_j_from_alpha_r(j, alpha, &parameters, &random_state);

      mpz_mul(tmp, j, r);
      mod_reduce(tmp, pow2ml);

      if (0 != mpz_cmp(tmp, alpha)) {
        critical("test_sample_j_from_alpha_r(): "
          "Failed to correctly sample j.");
      }

      /* Standard order. Negative. */
      sample_alpha_from_region(
        alpha, -min_log_alpha, -max_log_alpha, kappa(r), &random_state);
      sample_j_from_alpha_r(j, alpha, &parameters, &random_state);

      mpz_mul(tmp, j, r);
      mod_reduce(tmp, pow2ml);

      if (0 != mpz_cmp(tmp, alpha)) {
        critical("test_sample_j_from_alpha_r(): "
          "Failed to correctly sample j.");
      }

      /* Note: Calling with incorrect order would trigger a critical error. */
    }
  }

  /* Clear memory. */
  mpz_clear(d);
  mpz_clear(r);
  mpz_clear(alpha);
  mpz_clear(j);

  mpz_clear(tmp);

  mpz_clear(pow2ml);

  random_close(&random_state);
  parameters_clear(&parameters);
}

void test_sample_j_k_from_alpha_d() {
  printf("Testing sample_j_k_from_alpha_d()...\n");

  const uint32_t m = 2048;
  const uint32_t l = m / 8;
  const uint32_t t = 30;

  Random_State random_state;
  random_init(&random_state);

  mpz_t d;
  mpz_init(d);

  mpz_t r;
  mpz_init(r);

  mpz_t alpha;
  mpz_init(alpha);

  mpz_t j;
  mpz_init(j);

  mpz_t k;
  mpz_init(k);

  mpz_t tmp;
  mpz_init(tmp);

  mpz_t tmp2;
  mpz_init(tmp2);

  mpz_t pow2m;
  mpz_init_set_ui(pow2m, 0);
  mpz_setbit(pow2m, m);

  mpz_t pow2ml;
  mpz_init_set_ui(pow2ml, 0);
  mpz_setbit(pow2ml, m + l);

  mpz_t pow2mlm1;
  mpz_init_set_ui(pow2mlm1, 0);
  mpz_setbit(pow2mlm1, m + l - 1);

  parameters_selection_random_d_or_r(d, m);
  mpz_set(r, d); /* by convention */

  Parameters parameters;
  parameters_init(&parameters);
  parameters_explicit_m_l(&parameters, d, r, m, l, t);

  for (uint32_t i = m - 30; i < m + 30; i++) {
    for (uint32_t i2 = 0; i2 < 10; i2++) {
      double min_log_alpha = 10 * i + i2;     min_log_alpha /= 10;
      double max_log_alpha = 10 * i + i2 + 1; max_log_alpha /= 10;

      /* Standard order. Positive. */
      sample_alpha_from_region(
        alpha, min_log_alpha, max_log_alpha, kappa(d), &random_state);
      sample_j_k_from_alpha_d(j, k, alpha, &parameters, &random_state);

      mpz_mul(tmp, j, d);
      mpz_mul(tmp2, k, pow2m);
      mpz_add(tmp, tmp, tmp2);
      mpz_mod(tmp, tmp, pow2ml);
      if (mpz_cmp(tmp, pow2mlm1) >= 0) {
        mpz_sub(tmp, tmp, pow2ml);
      }

      if (0 != mpz_cmp(tmp, alpha)) {
        critical("test_sample_j_k_from_alpha_d(): "
          "Failed to correctly sample (j, k).");
      }

      /* Standard order. Negative. */
      sample_alpha_from_region(
        alpha, -min_log_alpha, -max_log_alpha, kappa(d), &random_state);
      sample_j_k_from_alpha_d(j, k, alpha, &parameters, &random_state);

      mpz_mul(tmp, j, d);
      mpz_mul(tmp2, k, pow2m);
      mpz_add(tmp, tmp, tmp2);
      mpz_mod(tmp, tmp, pow2ml);
      if (mpz_cmp(tmp, pow2mlm1) >= 0) {
        mpz_sub(tmp, tmp, pow2ml);
      }

      if (0 != mpz_cmp(tmp, alpha)) {
        critical("test_sample_j_k_from_alpha_d(): "
          "Failed to correctly sample (j, k).");
      }

      /* Note: Calling with incorrect order would trigger a critical error. */
    }
  }

  /* Clear memory. */
  mpz_clear(d);
  mpz_clear(r);
  mpz_clear(alpha);
  mpz_clear(j);
  mpz_clear(k);

  mpz_clear(tmp);
  mpz_clear(tmp2);

  mpz_clear(pow2m);
  mpz_clear(pow2ml);
  mpz_clear(pow2mlm1);

  random_close(&random_state);
  parameters_clear(&parameters);
}

void test_sample_j_k_from_alpha_d_r() {
  printf("Testing sample_j_k_from_alpha_d_r()...\n");

  const uint32_t m = 2048;
  const uint32_t l = m / 8;
  const uint32_t t = 30;

  Random_State random_state;
  random_init(&random_state);

  mpz_t d;
  mpz_init(d);

  mpz_t r;
  mpz_init(r);

  mpz_t alpha_d;
  mpz_init(alpha_d);

  mpz_t alpha_r;
  mpz_init(alpha_r);

  mpz_t original_alpha_d;
  mpz_init(original_alpha_d);

  mpz_t original_alpha_r;
  mpz_init(original_alpha_r);

  mpz_t j;
  mpz_init(j);

  mpz_t k;
  mpz_init(k);

  mpz_t tmp;
  mpz_init(tmp);

  mpz_t tmp2;
  mpz_init(tmp2);

  mpz_t pow2m;
  mpz_init_set_ui(pow2m, 0);
  mpz_setbit(pow2m, m);

  mpz_t pow2ml;
  mpz_init_set_ui(pow2ml, 0);
  mpz_setbit(pow2ml, m + l);

  mpz_t pow2mlm1;
  mpz_init_set_ui(pow2mlm1, 0);
  mpz_setbit(pow2mlm1, m + l - 1);

  mpz_t bound_diff_alpha;
  mpz_init_set_ui(bound_diff_alpha, 0);
  mpz_setbit(bound_diff_alpha, m / 2 + 24);

  mpfr_t alpha_f;
  mpfr_init2(alpha_f, PRECISION);

  parameters_selection_random_d_and_r(d, r, m);

  Parameters parameters;
  parameters_init(&parameters);
  parameters_explicit_m_l(&parameters, d, r, m, l, t);

  Distribution distribution;
  distribution_init(&distribution, &parameters, 1); /* for the lattice */

  for (uint32_t id = m - 30; id < m + 30; id++) {
  for (uint32_t id2 = 0; id2 < 10; id2++) {
  for (uint32_t id3 = 0; id3 <= 1; id3++) {
    double min_log_alpha_d = 10 * id + id2;     min_log_alpha_d /= 10;
    double max_log_alpha_d = 10 * id + id2 + 1; max_log_alpha_d /= 10;

    if (1 == id3) {
      min_log_alpha_d *= -1; max_log_alpha_d *= -1;
    }

    for (uint32_t ir = m - 30; ir < m + 30; ir++) {
    for (uint32_t ir2 = 0; ir2 < 10; ir2++) {
    for (uint32_t ir3 = 0; ir3 <= 1; ir3++) {
      double min_log_alpha_r = 10 * ir + ir2;     min_log_alpha_r /= 10;
      double max_log_alpha_r = 10 * ir + ir2 + 1; max_log_alpha_r /= 10;

      if (1 == ir3) {
        min_log_alpha_r *= -1; max_log_alpha_r *= -1;
      }

      /* Standard order. Positive. */
      sample_alpha_from_region(
        alpha_d, min_log_alpha_d, max_log_alpha_d, kappa(d), &random_state);
      sample_alpha_from_region(
        alpha_r, min_log_alpha_r, max_log_alpha_r, kappa(r), &random_state);

      mpz_set(original_alpha_d, alpha_d);
      mpz_set(original_alpha_r, alpha_r);

      /* Map the argument pair to the closest admissible argument pair. */
      lattice_alpha_map(
        alpha_d,
        alpha_r,
        &(distribution.lattice_alpha),
        &parameters);

      mpz_sub(tmp, original_alpha_d, alpha_d);
      mpz_abs(tmp, tmp);
      if (mpz_cmp(tmp, bound_diff_alpha) > 0) {
        critical("test_sample_j_k_from_alpha_d_r(): "
          "The mapped alpha_d is too far from the original alpha_d.");
      }

      mpz_sub(tmp, original_alpha_r, alpha_r);
      mpz_abs(tmp, tmp);
      if (mpz_cmp(tmp, bound_diff_alpha) > 0) {
        critical("test_sample_j_k_from_alpha_d_r(): "
          "The mapped alpha_r is too far from the original alpha_r.");
      }

      sample_j_k_from_alpha_d_r(
        j, k, alpha_d, alpha_r, &parameters, &random_state);

      /* Recompute alpha_d from (j, k). */
      mpz_mul(tmp, j, d);
      mpz_mul(tmp2, k, pow2m);
      mpz_add(tmp, tmp, tmp2);
      mpz_mod(tmp, tmp, pow2ml);
      if (mpz_cmp(tmp, pow2mlm1) >= 0) {
        mpz_sub(tmp, tmp, pow2ml);
      }

      if (0 != mpz_cmp(tmp, alpha_d)) {
        critical("test_sample_j_k_from_alpha_d_r(): Failed to sample (j, k).");
      }

      /* Recompute alpha_r from j. */
      mpz_mul(tmp, j, r);
      mpz_mod(tmp, tmp, pow2ml);
      if (mpz_cmp(tmp, pow2mlm1) >= 0) {
        mpz_sub(tmp, tmp, pow2ml);
      }

      if (0 != mpz_cmp(tmp, alpha_r)) {
        critical("test_sample_j_k_from_alpha_d_r(): Failed to sample j.");
      }


      /* Check the size and the compatibility requirement for alpha_d. */
      mpfr_set_z(alpha_f, alpha_d, MPFR_RNDN);
      mpfr_abs(alpha_f, alpha_f, MPFR_RNDN);
      mpfr_log2(alpha_f, alpha_f, MPFR_RNDN);

      double log_alpha = mpfr_get_d(alpha_f, MPFR_RNDN);

      if ((log_alpha < abs_d(min_log_alpha_d)) ||
          (log_alpha > abs_d(max_log_alpha_d)))
      {
        critical("test_sample_j_k_from_alpha_d_r(): "
          "Failed to sample alpha_d. Incorrect magnitude.");
      }

      mpz_mod_ui(tmp, alpha_d, 1UL << kappa(d));
      if (0 != mpz_cmp_ui(tmp, 0)) {
        critical("test_sample_j_k_from_alpha_d_r(): "
          "Failed to sample alpha_d. Not divisible by 2^kappa_d.");
      }


      /* Check the size and the compatibility requirement for alpha_r. */
      mpfr_set_z(alpha_f, alpha_r, MPFR_RNDN);
      mpfr_abs(alpha_f, alpha_f, MPFR_RNDN);
      mpfr_log2(alpha_f, alpha_f, MPFR_RNDN);

      log_alpha = mpfr_get_d(alpha_f, MPFR_RNDN);

      if ((log_alpha < abs_d(min_log_alpha_r)) ||
          (log_alpha > abs_d(max_log_alpha_r)))
      {
        critical("test_sample_j_k_from_alpha_d_r(): "
          "Failed to sample alpha_r. Incorrect magnitude.");
      }

      mpz_mod_ui(tmp, alpha_r, 1UL << kappa(r));
      if (0 != mpz_cmp_ui(tmp, 0)) {
        critical("test_sample_j_k_from_alpha_d_r(): "
          "Failed to sample alpha_r. Not divisible by 2^kappa_r.");
      }
    }
    }
    }
  }
  }
  }

  /* Clear memory. */
  mpz_clear(d);
  mpz_clear(r);
  mpz_clear(alpha_d);
  mpz_clear(alpha_r);
  mpz_clear(original_alpha_d);
  mpz_clear(original_alpha_r);
  mpz_clear(j);
  mpz_clear(k);

  mpz_clear(tmp);
  mpz_clear(tmp2);

  mpz_clear(pow2m);
  mpz_clear(pow2ml);
  mpz_clear(pow2mlm1);

  mpz_clear(bound_diff_alpha);

  mpfr_clear(alpha_f);

  random_close(&random_state);
  parameters_clear(&parameters);
  distribution_clear(&distribution);
}

void test_sample_j_from_diagonal_alpha_r() {
  printf("Testing sample_j_from_diagonal_alpha_r()...\n");

  const uint32_t t = 30;
  const uint32_t eta_bound = 0;
  const uint32_t m = 2048;
  const uint32_t sigma = 5;
  const uint32_t l = m + sigma;

  Random_State random_state;
  random_init(&random_state);

  mpz_t d;
  mpz_init(d);

  mpz_t r;
  mpz_init(r);

  mpz_t alpha_r;
  mpz_init(alpha_r);

  mpz_t j;
  mpz_init(j);

  mpz_t tmp;
  mpz_init(tmp);

  mpz_t pow2msigma;
  mpz_init_set_ui(pow2msigma, 0);
  mpz_setbit(pow2msigma, m + sigma);

  parameters_selection_random_d_and_r(d, r, m);

  Diagonal_Parameters parameters;
  diagonal_parameters_init(&parameters);
  diagonal_parameters_explicit_m_l(
    &parameters, d, r, m, sigma, l, eta_bound, t);

  for (uint32_t i = m - 30; i < m + sigma - 1; i++) {
    for (uint32_t i2 = 0; i2 < 10; i2++) {
    for (uint32_t i3 = 0; i3 <= 1; i3++) {
      double min_log_alpha_r = 10 * i + i2;     min_log_alpha_r /= 10;
      double max_log_alpha_r = 10 * i + i2 + 1; max_log_alpha_r /= 10;

      if (1 == i3) {
        min_log_alpha_r *= -1; max_log_alpha_r *= -1;
      }

      sample_alpha_from_region(
        alpha_r, min_log_alpha_r, max_log_alpha_r, kappa(r), &random_state);

      sample_j_from_diagonal_alpha_r(
        j, alpha_r, &parameters, &random_state);

      mpz_mul(tmp, j, r);
      mod_reduce(tmp, pow2msigma);

      if (0 != mpz_cmp(tmp, alpha_r)) {
        critical("test_sample_j_from_diagonal_alpha_r(): "
          "Failed to correctly sample j.");
      }

      /* Note: Calling with incorrect order would trigger a critical error. */
    }
    }
  }

  /* Clear memory. */
  mpz_clear(d);
  mpz_clear(r);
  mpz_clear(alpha_r);
  mpz_clear(j);

  mpz_clear(tmp);

  mpz_clear(pow2msigma);

  random_close(&random_state);
  diagonal_parameters_clear(&parameters);
}

void test_sample_k_from_diagonal_j_eta_pivot_kat() {
  printf("Testing test_sample_k_from_diagonal_j_eta_pivot_kat() via KAT...\n");

  /* As set in main() in main_generate_distribution.cpp. */
  mpfr_set_default_prec(PRECISION);

  const uint32_t m_s_entries[7][2] = {
    { 128, 10},
    { 256, 20},
    { 512, 30},
    {1024, 40},
    {2048, 50},
    {4096, 80},
    {8192, 80}
  };

  const uint32_t s_entries[16] = {
    1, 2, 3, 4, 5, 6, 7, 8, 10, 20, 30, 40, 50, 60, 70, 80
  };

  const uint32_t sigma_entries[6] = {
    0, 1, 2, 3, 4, 5
  };

  for (uint32_t m_index = 0; m_index < 7; m_index++) {
    const uint32_t m = m_s_entries[m_index][0];
    const uint32_t s_max = m_s_entries[m_index][1];

    for (uint32_t s_index = 0; s_index < 16; s_index++) {
      const uint32_t s = s_entries[s_index];
      if (s > s_max) {
        break;
      }

      for (uint32_t sigma_index = 0; sigma_index < 6; sigma_index++) {
        const uint32_t sigma = sigma_entries[sigma_index];
        const uint32_t l = ceil(((double)m) / ((double)s));

        const uint32_t precision = 2 * l;

        mpz_t d;
        mpz_init(d);

        mpz_t r;
        mpz_init(r);

        parameters_selection_deterministic_d_r(d, r, m);

        Diagonal_Parameters parameters;
        diagonal_parameters_init(&parameters);
        diagonal_parameters_explicit_m_l(
          &parameters,
          d,
          r,
          m,
          sigma,
          l,
          25, /* eta_bound = 25 */
          30); /* t = 30 */

        char path[MAX_BUFFER_SIZE];
        sprintf(path,
          "res/test-vectors/sample-k-from-diagonal-j-eta-pivot-"
            "m-%u-sigma-%u-s-%u.txt", m, sigma, s);
        printf(" Processing: %s\n", path);

        FILE * file = fopen(path, "rb");
        if (NULL == file) {
          critical("test_sample_k_from_diagonal_j_eta_pivot_kat(): "
            "Failed to open \"%s\".", path);
        }

        /* Declare variables. */
        mpz_t j;
        mpz_init(j);

        mpz_t k;
        mpz_init(k);

        mpz_t exp_k;
        mpz_init(exp_k);

        mpfr_t alpha_phi;
        mpfr_init2(alpha_phi, precision);

        mpfr_t exp_alpha_phi;
        mpfr_init2(exp_alpha_phi, precision);

        for (uint32_t iteration = 0; iteration < 25; iteration++) {
          /* Load j. */
          test_mpz_load(j, file);

          /* Load eta. */
          int32_t eta;
          if (1 != fscanf(file, "%d\n", &eta)) {
            critical("test_sample_k_from_diagonal_j_eta_pivot_kat(): "
              "Failed to load eta from file.");
          }

          /* Load pivot. */
          long double pivot;
          if (1 != fscanf(file, "%Lf\n", &pivot)) {
            critical("test_sample_k_from_diagonal_j_eta_pivot_kat(): "
              "Failed to load pivot from file.");
          }

          /* Load expected k. */
          test_mpz_load(exp_k, file);

          /* Load expected alpha_phi. */
          test_mpfr_load(exp_alpha_phi, file);

          /* Compute k and alpha_phi. */
          sample_k_from_diagonal_j_eta_pivot(
            &parameters,
            pivot,
            j,
            eta,
            UINT32_MAX, /* = delta_bound */
            k,
            alpha_phi);

          if (0 != mpz_cmp(exp_k, k)) {
            critical("test_sample_k_from_diagonal_j_eta_pivot_kat(): "
              "Failed to verify k.");
          }

          const long double alpha_phi_ld =
            mpfr_get_ld(alpha_phi, MPFR_RNDN);
          const long double exp_alpha_phi_ld =
            mpfr_get_ld(exp_alpha_phi, MPFR_RNDN);

          if (TRUE != test_cmp_ld(exp_alpha_phi_ld, alpha_phi_ld)) {
            critical("test_sample_k_from_diagonal_j_eta_pivot_kat(): "
              "Failed to verify alpha_phi.");
          }
        }

        /* Close the file. */
        fclose(file);
        file = NULL;

        /* Clear memory. */
        mpz_clear(d);
        mpz_clear(r);

        mpz_clear(j);
        mpz_clear(k);
        mpz_clear(exp_k);

        mpfr_clear(alpha_phi);
        mpfr_clear(exp_alpha_phi);

        diagonal_parameters_clear(&parameters);
      }
    }
  }
}

void test_sample() {
  test_sample_approximate_alpha_from_region();
  test_sample_alpha_from_region();

  test_sample_j_from_alpha_r();
  test_sample_j_k_from_alpha_d();
  test_sample_j_k_from_alpha_d_r();
  test_sample_j_from_diagonal_alpha_r();
}