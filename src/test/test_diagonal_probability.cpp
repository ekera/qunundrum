/*!
 * \file    test/test_diagonal_probability.cpp
 * \ingroup unit_tests_diagonal_distribution_probability
 *
 * \brief   The definition of unit tests for the probability functions in the
 *          \ref diagonal_distribution_probability module.
 */

#include "test_diagonal_probability.h"

#include "test_common.h"

#include "../diagonal_parameters.h"
#include "../diagonal_probability.h"
#include "../errors.h"
#include "../math.h"
#include "../parameters_selection.h"

#include <gmp.h>
#include <mpfr.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>

/*!
 * \brief   The maximum buffer size in bytes.
 */
#define MAX_BUFFER_SIZE 8192

void test_diagonal_probability_f_eta_approx_kat() {
  printf("Testing diagonal_probability_f_eta_approx() via KAT...\n");

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

  for (uint32_t i = 0; i < 7; i++) {
    const uint32_t m = m_s_entries[i][0];
    const uint32_t s_max = m_s_entries[i][1];

    for (uint32_t j = 0; j < 16; j++) {
      const uint32_t s = s_entries[j];
      if (s > s_max) {
        break;
      }

      for (uint32_t k = 0; k < 6; k++) {
        const uint32_t sigma = sigma_entries[k];
        const uint32_t l = ceil(((double)m) / ((double)s));

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

        mpz_t alpha_r;
        mpz_init(alpha_r);

        mpfr_t theta_r;
        mpfr_init2(theta_r, PRECISION);

        mpfr_t scale_factor;
        mpfr_init2(scale_factor, PRECISION);

        mpfr_t norm;
        mpfr_init2(norm, PRECISION);

        mpfr_t exp_norm;
        mpfr_init2(exp_norm, PRECISION);

        mpfr_t tmp;
        mpfr_init2(tmp, PRECISION);

        /* Compute the scale factor. */
        mpfr_const_pi(scale_factor, MPFR_RNDN);
        mpfr_mul_ui(scale_factor, scale_factor, 2, MPFR_RNDN);

        mpfr_set_ui_2exp(tmp, 1, (mpfr_exp_t)(m + sigma), MPFR_RNDN);
        mpfr_div(scale_factor, scale_factor, tmp, MPFR_RNDN);

        char path[MAX_BUFFER_SIZE];
        sprintf(path,
          "res/test-vectors/diagonal-probabilities-f-eta-"
            "det-m-%u-sigma-%u-s-%u.txt",
              m, sigma, s);
        printf(" Processing: %s\n", path);

        FILE * file = fopen(path, "rb");
        if (NULL == file) {
          critical("test_diagonal_probability_f_eta_approx_kat(): "
            "Failed to open \"%s\".", path);
        }

        /* Load alpha_r. */
        test_mpz_load(alpha_r, file);

        /* Compute theta_r. */
        mpfr_mul_z(theta_r, scale_factor, alpha_r, MPFR_RNDN);

        for (int32_t eta = -25; eta <= 25; eta++) {
          diagonal_probability_approx_f_eta(norm, theta_r, eta, &parameters);

          test_mpfr_load(exp_norm, file);

          const long double norm_ld = mpfr_get_ld(norm, MPFR_RNDN);
          const long double exp_norm_ld = mpfr_get_ld(exp_norm, MPFR_RNDN);

          if ((exp_norm_ld <= 0) || (norm_ld <= 0)) {
            critical("test_diagonal_probability_f_eta_approx_kat(): "
              "Failed to verify probability estimate.");
          }

          if (TRUE != test_cmp_ld(norm_ld, exp_norm_ld)) {
            critical("test_diagonal_probability_f_eta_approx_kat(): "
              "Failed to verify probability estimate.");
          }
        }

        /* Close the file. */
        fclose(file);
        file = NULL;

        /* Clear memory. */
        mpz_clear(d);
        mpz_clear(r);

        mpz_clear(alpha_r);

        mpfr_clear(theta_r);
        mpfr_clear(scale_factor);
        mpfr_clear(norm);
        mpfr_clear(exp_norm);
        mpfr_clear(tmp);

        diagonal_parameters_clear(&parameters);
      }
    }
  }
}

void test_diagonal_probability_h_approx_kat() {
  printf("Testing diagonal_probability_h_approx() via KAT...\n");

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

  for (uint32_t i = 0; i < 7; i++) {
    const uint32_t m = m_s_entries[i][0];
    const uint32_t s_max = m_s_entries[i][1];

    for (uint32_t j = 0; j < 16; j++) {
      const uint32_t s = s_entries[j];
      if (s > s_max) {
        break;
      }

      for (uint32_t k = 0; k < 6; k++) {
        const uint32_t sigma = sigma_entries[k];
        const uint32_t l = ceil(((double)m) / ((double)s));

        const uint32_t precision = 3 * l;

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
          "res/test-vectors/diagonal-probabilities-h-"
            "det-m-%u-sigma-%u-s-%u.txt",
              m, sigma, s);
        printf(" Processing: %s\n", path);

        FILE * file = fopen(path, "rb");
        if (NULL == file) {
          critical("test_diagonal_probability_h_approx_kat(): "
            "Failed to open \"%s\".", path);
        }

        /* Declare variables. */
        mpfr_t phi;
        mpfr_init2(phi, precision);

        mpfr_t norm;
        mpfr_init2(norm, precision);

        mpfr_t exp_norm;
        mpfr_init2(exp_norm, precision);

        for (uint32_t iteration = 0; iteration < 25; iteration++) {
          /* Load phi. */
          test_mpfr_load(phi, file);

          /* Load expected norm. */
          test_mpfr_load(exp_norm, file);

          /* Compute norm. */
          diagonal_probability_approx_h(norm, phi, &parameters);

          const long double norm_ld = mpfr_get_ld(norm, MPFR_RNDN);
          const long double exp_norm_ld = mpfr_get_ld(exp_norm, MPFR_RNDN);

          if ((exp_norm_ld <= 0) || (norm_ld <= 0)) {
            critical("test_diagonal_probability_h_approx_kat(): "
              "Failed to verify probability estimate.");
          }

          if (TRUE != test_cmp_ld(norm_ld, exp_norm_ld)) {
            critical("test_diagonal_probability_h_approx_kat(): "
              "Failed to verify probability estimate.");
          }
        }

        /* Close the file. */
        fclose(file);
        file = NULL;

        /* Clear memory. */
        mpz_clear(d);
        mpz_clear(r);

        mpfr_clear(phi);

        mpfr_clear(norm);
        mpfr_clear(exp_norm);

        diagonal_parameters_clear(&parameters);
      }
    }
  }
}

void test_diagonal_probability() {
  test_diagonal_probability_f_eta_approx_kat();
  test_diagonal_probability_h_approx_kat();
}
