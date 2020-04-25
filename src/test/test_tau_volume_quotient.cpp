/*!
 * \file    test/test_tau_volume_quotient.cpp
 * \ingroup unit_tests_estimating_volume_quotients
 * 
 * \brief   The definition of unit tests for the 
 *          \ref estimating_volume_quotients module.
 */

#include "test_tau_volume_quotient.h"

#include "test_common.h"

#include "../tau_volume_quotient.h"
#include "../parameters_selection.h"
#include "../parameters.h"
#include "../errors.h"
#include "../math.h"

#include <gmp.h>
#include <mpfr.h>

#include <math.h>

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

void test_tau_volume_quotient_kat() {
  printf("Testing tau_volume_quotient() via KAT...\n");

  uint32_t m;
  uint32_t s;
  uint32_t n;

  long double tau;

  mpfr_t v;
  mpfr_init2(v, PRECISION);

  mpfr_t exp_v;
  mpfr_init2(exp_v, PRECISION);

  mpz_t d;
  mpz_init(d);

  mpz_t r;
  mpz_init(r);

  const static char * const path = "res/test-vectors/tau-volume-quotients.txt";

  FILE * file = fopen(path, "rb");
  if (NULL == file) {
    critical("test_tau_volume_quotient_kat(): Failed to open \"%s\".", path);
  }

  uint32_t count;

  if (1 != fscanf(file, "%u\n", &count)) {
    critical("test_tau_volume_quotient_kat(): Failed to parse count.");
  }
  
  for (uint32_t i = 0; i < count; i++) {
    /* Read parameters. */
    if (4 != fscanf(file, "%u %u %u %Lf\n", &m, &s, &n, &tau)) {
      critical("test_tau_volume_quotient_kat(): "
        "Failed to parse m, s, n and tau.");
    }

    /* Read the expected volume quotient. */
    test_mpfr_load(exp_v, file);

    /* Setup parameters. */   
    parameters_selection_deterministic_d_r(d, r, m);
    const uint32_t l = ceil(((double)m) / ((double)s));

    /* Compute the volume quotient v. */
    tau_volume_quotient(m, l, n, d, tau, v);

    /* Verify v by comparing it to v_exp. */
    if (!test_cmp_tol_ld(mpfr_get_ld(v, MPFR_RNDN), 
                         mpfr_get_ld(exp_v, MPFR_RNDN),
                         1e-3))
    {
      critical("test_tau_volume_quotient_kat(): "
        "Failed to verify volume quotient for (m, s, n) = (%u, %u, %u).",
          m, s, n);
    }
  }

  /* Clean up */
  fclose(file);
  file = NULL;

  mpz_clear(d);
  mpz_clear(r);
  
  mpfr_clear(v);
  mpfr_clear(exp_v);
}

void test_tau_volume_quotient() {
  test_tau_volume_quotient_kat();
}
