/*!
 * \file    test/test_diagonal_probability.cpp
 * \ingroup unit_tests_diagonal_distribution_probability
 * 
 * \brief   The definition of unit tests for the probability functions in the 
 *          \ref diagonal_distribution_probability module.
 */

#include "test_diagonal_probability.h"

#include "test_common.h"

#include "../diagonal_probability.h"
#include "../parameters_selection.h"
#include "../parameters.h"
#include "../math.h"
#include "../errors.h"

#include <gmp.h>
#include <mpfr.h>

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

/*!
 * \brief   The maximum buffer size in bytes.
 */
#define MAX_BUFFER_SIZE 8192

void test_diagonal_probability_approx_kat() {
  printf("Testing diagonal_probability() via KAT...\n");

  /* As set in main() in main_generate_distribution.cpp. */
  mpfr_set_default_prec(PRECISION);

  /* As set in main() in main_generate_distribution.cpp. */
  const static uint32_t m_entries[7] = {
    128, 256, 512, 1024, 2048, 4096, 8192
  };
  
  const static uint32_t l_entries[6] = {
    0, 1, 2, 3, 4, 5
  };

  for (uint32_t i = 0; i < 7; i++) {
    const uint32_t m = m_entries[i];

    for (uint32_t j = 0; j < 6; j++) {
      const uint32_t l = l_entries[j];

      const uint32_t precision = 2 * (m + l); /* recommended setting */

      mpz_t d;
      mpz_init(d);

      mpz_t r;
      mpz_init(r);
      
      parameters_selection_deterministic_d_r(d, r, m);

      Parameters parameters;
      parameters_init(&parameters);
      parameters_explicit_m_l(&parameters, d, r, m, l, 30); /* t = 30 */

      mpz_t alpha_d;
      mpz_init(alpha_d);

      mpz_t alpha_r;
      mpz_init(alpha_r);

      mpfr_t theta_d;
      mpfr_init2(theta_d, precision);

      mpfr_t theta_r;
      mpfr_init2(theta_r, precision);

      mpfr_t scale_factor;
      mpfr_init2(scale_factor, precision);

      mpfr_const_pi(scale_factor, MPFR_RNDN);
      mpfr_mul_ui(scale_factor, scale_factor, 2, MPFR_RNDN);

      {
        mpfr_t tmp;
        mpfr_init2(tmp, precision);
        
        mpfr_set_ui_2exp(tmp, 1, (mpfr_exp_t)(m + l), MPFR_RNDN);
        mpfr_div(scale_factor, scale_factor, tmp, MPFR_RNDN);

        mpfr_clear(tmp);
      }

      mpfr_t norm;
      mpfr_init2(norm, PRECISION);

      mpfr_t exp_norm;
      mpfr_init2(exp_norm, PRECISION);

      char path[MAX_BUFFER_SIZE];
      sprintf(path,
        "res/test-vectors/diagonal-probabilities-det-m-%u-l-%u.txt", m, l);
      printf(" Processing: %s\n", path);

      FILE * file = fopen(path, "rb");
      if (NULL == file) {
        critical("test_diagonal_probability_kat(): "
          "Failed to open \"%s\".", path);
      }

      for (uint32_t t = m - 30; t < m + l - 1; t++) {
        /* Read alpha_r and compute_d. */
        test_mpz_load(alpha_r, file);
        mpfr_mul_z(theta_r, scale_factor, alpha_r, MPFR_RNDN);

        for (int32_t delta = -20; delta <= 20; delta++) {
          /* Compute alpha_d. */
          mpfr_set_z(theta_d, alpha_r, MPFR_RNDN);
          mpfr_mul_z(theta_d, theta_d, parameters.d, MPFR_RNDN);
          mpfr_div_z(theta_d, theta_d, parameters.r, MPFR_RNDN);
          mpfr_round(theta_d, theta_d);
          mpfr_get_z(alpha_d, theta_d, MPFR_RNDN);

          if (delta < 0) {
            mpz_sub_ui(alpha_d, alpha_d, (uint32_t)(-delta));
          } else {
            mpz_add_ui(alpha_d, alpha_d, (uint32_t)( delta));
          }

          /* Compute theta_d. */
          mpfr_mul_z(theta_d, scale_factor, alpha_d, MPFR_RNDN);

          /* Compute and compare. */
          diagonal_probability_approx(norm, theta_d, theta_r, &parameters);

          test_mpfr_load(exp_norm, file);

          const long double norm_ld = mpfr_get_ld(norm, MPFR_RNDN);
          const long double exp_norm_ld = mpfr_get_ld(exp_norm, MPFR_RNDN);

          if ((exp_norm_ld <= 0) || (norm_ld <= 0)) {
            critical("test_diagonal_probability_kat(): "
              "Failed to verify probability estimate.");
          }

          if (!test_cmp_ld(norm_ld, exp_norm_ld)) {
            critical("test_diagonal_probability_kat(): "
              "Failed to verify probability estimate.");
          }
        }
      }

      fclose(file);
      file = NULL;

      mpz_clear(d);
      mpz_clear(r);

      mpz_clear(alpha_d);
      mpz_clear(alpha_r);
      
      mpfr_clear(theta_d);
      mpfr_clear(theta_r);
      mpfr_clear(scale_factor);
      mpfr_clear(norm);
      mpfr_clear(exp_norm);
      
      parameters_clear(&parameters);
    }
  }
}

void test_diagonal_probability() {
  test_diagonal_probability_approx_kat();
}
