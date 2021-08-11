/*!
 * \file    test/test_probability.cpp
 * \ingroup unit_tests_two_dimensional_distribution
 * 
 * \brief   The definition of unit tests for the probability functions in the 
 *          \ref two_dimensional_distribution_probability module.
 */

#include "test_probability.h"

#include "test_common.h"

#include "../errors.h"
#include "../math.h"
#include "../parameters.h"
#include "../parameters_selection.h"
#include "../probability.h"

#include <gmp.h>
#include <mpfr.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>

/*! 
 * \brief   The maximum buffer size in bytes.
 */
#define MAX_BUFFER_SIZE 8192

void test_probability_approx_heuristic_sigma_kat() {
  printf("Testing probability_approx() for heuristic sigma via KAT...\n");

  /* As set in main() in main_generate_distribution.cpp. */
  mpfr_set_default_prec(PRECISION);

  /* As set in main() in main_generate_distribution.cpp. */
  const uint32_t m_s_entries[8][2] = {
    { 128,  2}, { 256,  4}, { 384,  6}, { 512,  8},
    {1024, 10}, {2048, 30}, {4096, 50}, {8192, 80}
  };
  
  const uint32_t s_entries[14] = {
    1, 2, 3, 4, 5, 6, 7, 8, 10, 20, 30, 40, 50, 80
  };

  for (uint32_t i = 0; i < 8; i++) {
    const uint32_t m = m_s_entries[i][0];
    const uint32_t s_max = m_s_entries[i][1];

    for (uint32_t j = 0; j < 14; j++) {
      const uint32_t s = s_entries[j];
      if (s > s_max) {
        break;
      }

      mpz_t d;
      mpz_init(d);

      mpz_t r;
      mpz_init(r);
      
      parameters_selection_deterministic_d_r(d, r, m);

      const uint32_t t = 30;
  
      Parameters parameters;
      parameters_init(&parameters);
      parameters_explicit_m_s(&parameters, d, r, m, s, t);
        
      /* Select sigma heuristically as in generate_distribution. */
      const uint32_t t_sigma = 11;

      /* Note that we could extract t from the parameters above. However, we 
       * have also imposed a hard limit on t not growing past 10 when generating
       * the distribution, and this limit is not accounted for when reading from
       * the parameters. Note also that we set t to the maximum possible value 
       * above. It may be possible to do even better by varying t. */

      const uint32_t sigma =
        round(((float)parameters.l + t_sigma + 4 - 1.6515f) / 2.0f);
      
      mpfr_t norm;
      mpfr_init2(norm, PRECISION);

      mpfr_t error;
      mpfr_init2(error, PRECISION);

      mpfr_t theta_d;
      mpfr_init2(theta_d, PRECISION);

      mpfr_t theta_r;
      mpfr_init2(theta_r, PRECISION);

      mpfr_t exp_norm;
      mpfr_init2(exp_norm, PRECISION);

      mpfr_t exp_error;
      mpfr_init2(exp_error, PRECISION);

      mpfr_t quick_norm;
      mpfr_init2(quick_norm, PRECISION);

      char path[MAX_BUFFER_SIZE];
      sprintf(path, "res/test-vectors/probabilities-det-m-%u-s-%u.txt", m, s);
      printf(" Processing: %s\n", path);

      FILE * file = fopen(path, "rb");
      if (NULL == file) {
        critical("test_probability_approx_heuristic_sigma_kat(): "
          "Failed to open \"%s\".", path);
      }

      for (uint32_t k = 0; k < 1764; k++) {
        test_mpfr_load(theta_d, file);
        test_mpfr_load(theta_r, file);
        
        probability_approx(
          norm,
          error,
          sigma,
          theta_d,
          theta_r,
          &parameters);

        test_mpfr_load(exp_norm, file);
        test_mpfr_load(exp_error, file);

        const long double norm_ld = mpfr_get_ld(norm, MPFR_RNDN);
        const long double exp_norm_ld = mpfr_get_ld(exp_norm, MPFR_RNDN);

        const long double error_ld = mpfr_get_ld(error, MPFR_RNDN);
        const long double exp_error_ld = mpfr_get_ld(exp_error, MPFR_RNDN);

        if ((exp_norm_ld <= 0) || (norm_ld <= 0)) {
          critical("test_probability_approx_heuristic_sigma_kat(): "
            "Failed to verify probability estimate.");
        }

        if (!test_cmp_ld(norm_ld, exp_norm_ld)) {
          critical("test_probability_approx_heuristic_sigma_kat(): "
            "Failed to verify probability estimate.");
        }

        if ((exp_error_ld <= 0) || (error_ld <= 0)) {
          critical("test_probability_approx_heuristic_sigma_kat(): "
            "Failed to verify error bound estimate.");
        }

        if (!test_cmp_ld(error_ld, exp_error_ld)) {
          critical("test_probability_approx_heuristic_sigma_kat(): "
            "Failed to verify error bound estimate.");
        }

        /* Also test the approximate norm. */
        probability_approx_quick(quick_norm, theta_d, theta_r, &parameters);

        const long double quick_norm_ld = mpfr_get_ld(quick_norm, MPFR_RNDN);

        if (quick_norm_ld <= 0) {
          critical("test_probability_approx_heuristic_sigma_kat(): "
            "Failed to verify probability estimate.");
        }
        
        if (!test_cmp_tol_ld(norm_ld, quick_norm_ld, 1e-4)) {
          critical("test_probability_approx_heuristic_sigma_kat(): "
            "Failed to verify ** quick ** probability estimate.");
        }
      }

      /* Close the file. */
      fclose(file);
      file = NULL;

      /* Clear memory. */
      mpfr_clear(norm);
      mpfr_clear(error);
      mpfr_clear(theta_d);
      mpfr_clear(theta_r);
      mpfr_clear(exp_norm);
      mpfr_clear(exp_error);
      mpfr_clear(quick_norm);
      
      mpz_clear(d);
      mpz_clear(r);
      
      parameters_clear(&parameters);
    }
  }
}

void test_probability() {
  test_probability_approx_heuristic_sigma_kat();
}