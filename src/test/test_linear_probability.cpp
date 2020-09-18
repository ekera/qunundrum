/*!
 * \file    test/test_linear_probability.cpp
 * \ingroup unit_tests_linear_distribution_probability
 * 
 * \brief   The definition of unit tests for the probability functions in the 
 *          \ref linear_distribution_probability module.
 */

#include "test_linear_probability.h"

#include "test_common.h"

#include "../linear_probability.h"
#include "../parameters_selection.h"
#include "../parameters.h"
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

void test_linear_probability_d_kat() {
  printf("Testing linear_probability_d() via KAT...\n");

  /* As set in main() in main_generate_distribution.cpp. */
  mpfr_set_default_prec(PRECISION);

  /* As set in main() in main_generate_distribution.cpp. */
  const static uint32_t m_s_entries[7][2] = {
    { 128, 10}, { 256, 20}, { 512, 30},
    {1024, 40}, {2048, 50}, {4096, 80}, {8192, 80}
  };
  
  const static uint32_t s_entries[16] = {
    1, 2, 3, 4, 5, 6, 7, 8, 10, 20, 30, 40, 50, 60, 70, 80
  };

  for (uint32_t i = 0; i < 7; i++) {
    const uint32_t m = m_s_entries[i][0];
    const uint32_t s_max = m_s_entries[i][1];

    for (uint32_t j = 0; j < 16; j++) {
      const uint32_t s = s_entries[j];
      if (s > s_max) {
        break;
      }

      mpz_t d;
      mpz_init(d);

      mpz_t r;
      mpz_init(r);
      
      parameters_selection_deterministic_d_r(d, r, m);
      mpz_set(r, d); /* by convention */
    
      const uint32_t t = 30;
  
      Parameters parameters;
      parameters_init(&parameters);
      parameters_explicit_m_s(&parameters, d, r, m, s, t);

      mpfr_t norm;
      mpfr_init2(norm, PRECISION);

      mpfr_t theta_d;
      mpfr_init2(theta_d, PRECISION);

      mpfr_t exp_norm;
      mpfr_init2(exp_norm, PRECISION);

      char path[MAX_BUFFER_SIZE];
      sprintf(path, 
        "res/test-vectors/linear-probabilities-det-d-m-%u-s-%u.txt", m, s);
      printf(" Processing: %s\n", path);

      FILE * file = fopen(path, "rb");
      if (NULL == file) {
        critical("test_linear_probability_d_kat(): "
          "Failed to open \"%s\".", path);
      }

      for (uint32_t k = 0; k < 61; k++) {
        test_mpfr_load(theta_d, file);
        
        linear_probability_d(norm, theta_d, &parameters);

        test_mpfr_load(exp_norm, file);

        const long double norm_ld = mpfr_get_ld(norm, MPFR_RNDN);
        const long double exp_norm_ld = mpfr_get_ld(exp_norm, MPFR_RNDN);

        if ((exp_norm_ld < 0) || (norm_ld < 0)) {
          critical("test_linear_probability_d_kat(): "
            "Failed to verify probability estimate.");
        }

        if (!test_cmp_ld(norm_ld, exp_norm_ld)) {
          critical("test_linear_probability_d_kat(): "
            "Failed to verify probability estimate.");
        }
      }

      /* Close the file. */
      fclose(file);
      file = NULL;

      /* Clear memory. */
      mpfr_clear(norm);
      mpfr_clear(theta_d);
      mpfr_clear(exp_norm);
      
      mpz_clear(d);
      mpz_clear(r);
      
      parameters_clear(&parameters);
    }
  }
}

void test_linear_probability_r_kat() {
  printf("Testing linear_probability_r() via KAT...\n");

  /* As set in main() in main_generate_distribution.cpp. */
  mpfr_set_default_prec(PRECISION);

  /* As set in main() in main_generate_distribution.cpp. */
  const static uint32_t m_s_entries[7][2] = {
    { 128, 10}, { 256, 20}, { 512, 30},
    {1024, 40}, {2048, 50}, {4096, 80}, {8192, 80}
  };
  
  const static uint32_t s_entries[16] = {
    1, 2, 3, 4, 5, 6, 7, 8, 10, 20, 30, 40, 50, 60, 70, 80
  };

  for (uint32_t i = 0; i < 7; i++) {
    const uint32_t m = m_s_entries[i][0];
    const uint32_t s_max = m_s_entries[i][1];

    for (uint32_t j = 0; j < 16; j++) {
      const uint32_t s = s_entries[j];
      if (s > s_max) {
        break;
      }

      mpz_t d;
      mpz_init(d);

      mpz_t r;
      mpz_init(r);
      
      parameters_selection_deterministic_d_r(d, r, m);
      mpz_set(d, r); /* by convention */
    
      const uint32_t t = 30;
  
      Parameters parameters;
      parameters_init(&parameters);
      parameters_explicit_m_s(&parameters, d, r, m, s, t);

      mpfr_t norm;
      mpfr_init2(norm, PRECISION);

      mpfr_t theta_r;
      mpfr_init2(theta_r, PRECISION);

      mpfr_t exp_norm;
      mpfr_init2(exp_norm, PRECISION);

      char path[MAX_BUFFER_SIZE];
      sprintf(path,
        "res/test-vectors/linear-probabilities-det-r-m-%u-s-%u.txt", m, s);
      printf(" Processing: %s\n", path);

      FILE * file = fopen(path, "rb");
      if (NULL == file) {
        critical("test_linear_probability_r_kat(): "
          "Failed to open \"%s\".", path);
      }

      for (uint32_t k = 0; k < 61; k++) {
        test_mpfr_load(theta_r, file);
        
        linear_probability_r(norm, theta_r, &parameters);

        test_mpfr_load(exp_norm, file);
        mpfr_set(norm, norm, MPFR_RNDN);

        const long double norm_ld = mpfr_get_ld(norm, MPFR_RNDN);
        const long double exp_norm_ld = mpfr_get_ld(exp_norm, MPFR_RNDN);

        if ((exp_norm_ld <= 0) || (norm_ld <= 0)) {
          critical("test_linear_probability_r_kat(): "
            "Failed to verify probability estimate.");
        }

        if (!test_cmp_ld(norm_ld, exp_norm_ld)) {
          critical("test_linear_probability_r_kat(): "
            "Failed to verify probability estimate.");
        }
      }

      /* Close the file. */
      fclose(file);
      file = NULL;

      /* Clear memory. */
      mpfr_clear(norm);
      mpfr_clear(theta_r);
      mpfr_clear(exp_norm);
      
      mpz_clear(d);
      mpz_clear(r);
      
      parameters_clear(&parameters);
    }
  }
}

void test_linear_probability() {
  test_linear_probability_d_kat();
  test_linear_probability_r_kat();
}