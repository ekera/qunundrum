/*!
 * \file    continued_fractions.cpp
 * \ingroup continued_fractions
 *
 * \brief   The definition of functions for computing continued fractions.
 *
 * These functions are used to implement Shor's solver for order finding.
 */

#include "continued_fractions.h"

#include "parameters.h"
#include "errors.h"

#include <mpfr.h>
#include <gmp.h>  

#include <stdint.h>
#include <stdlib.h>

bool continued_fractions_solve(
  const mpz_t j,
  const uint32_t search_bound_cofactor,
  const Parameters * const parameters)
{
  /* Setup the precision. */
  const uint32_t precision = 2 * (parameters->m + parameters->l);

  /* Setup constants. */
  mpfr_t one;
  mpfr_init2(one, precision);
  mpfr_set_ui(one, 1, MPFR_RNDN);
  
  mpz_t pow2m;
  mpz_init_set_ui(pow2m, 0);
  mpz_setbit(pow2m, parameters->m);

  mpz_t pow2ml;
  mpz_init_set_ui(pow2ml, 0);
  mpz_setbit(pow2ml, parameters->m + parameters->l);

  mpfr_t fraction;
  mpfr_init2(fraction, precision);
  mpfr_set_z(fraction, j, MPFR_RNDN);
  mpfr_div_z(fraction, fraction, pow2ml, MPFR_RNDN);

  mpfr_t integer_part;
  mpfr_init2(integer_part, precision);

  mpz_t denominator;
  mpz_init(denominator);

  mpz_t km1;
  mpz_init_set_ui(km1, 0);

  mpz_t km2;
  mpz_init_set_ui(km2, 1);

  mpz_t tmp;
  mpz_init(tmp);

  /* We recursively compute the continued fractions coefficients 
   * 
   *          fraction = [ a_0; a_1, a_2, .. ]
   * 
   * below by letting f_0 = fraction, a_i = floor(f_i), and recursively letting
   * f_{i+1} = 1 / (f_i - floor(f_i)) = 1 / (f_i - a_i) provided f_i ≠ a_i. If 
   * it is the case that f_i = a_i we have obtained the best approximation.
   * 
   * From the coefficients we can construct increasingly good quotients that
   * approximate the fraction. The first quotient is a_0, the second a_0 + 1/a_1
   * and so forth, recursively, as we have that
   * 
   *                            
   *                                         1
   *           fraction = a_0 + ---------------------------
   *                                            1
   *                            a_1 + ---------------------
   *                                               1
   *                                  a_2 + ---------------
   *                                                  1
   *                                        a_3 + --------- .
   *                                              a_4 + ...
   * 
   * We may however compute these quotients, called convergents, considerably 
   * more efficiently than by using the above formula. Specifically:
   * 
   * Let h_{-1} = 1, h_{-2} = 0 and recursively h_i = a_n h_{i-1} + h_{i-2}.
   * 
   * Let k_{-1} = 0, k_{-2} = 1 and recursively k_i = a_n k_{i-1} + k_{i-2}.
   * 
   * Then the i:th convergent is given by the quotient h_i / k_i for i >= 0.
   * 
   * Now, since we are only interested in the denominator, not the numerator, it
   * suffices to recursivelty compute a_0, k_0, a_1, k_1, .. keeping a_{i-1}, 
   * k_{i-1} and k_{i-2} in memory when computing a_i and k_i. */

  bool found = FALSE;

  for (uint32_t depth = 1; ; depth++) {
    mpfr_floor(integer_part, fraction); /* next cofficient */

    /* Construct the next denominator. */
    mpfr_get_z(denominator, integer_part, MPFR_RNDN);
    mpz_mul(denominator, denominator, km1);
    mpz_add(denominator, denominator, km2);
    
    if (mpz_cmp(denominator, pow2m) >= 0) {
      break; /* The denominator is larger than r. */
    }

    /* Test if the denominator is r, or a multiple of r, before proceeding.
     * 
     * As was originally pointed out by Odlyzko, we need to check multiples, 
     * since we could have z / r, where z and r share a common factor z, 
     * leading us to obtain (z / c) / (r / c). We solve this by checking if 
     * this holds for c within the search bound on the cofactor. */
    mpz_mod(tmp, parameters->r, denominator);
    if (mpz_cmp_ui(tmp, 0) == 0) {
      /* The denominator devides r. Set tmp = c and test the size of c. */
      mpz_div(tmp, parameters->r, denominator);

      /* Check whether the cofactor is within the search bound. */
      if (mpz_cmp_ui(tmp, search_bound_cofactor) <= 0) {
        found = TRUE;
        break;
      }
    }

    /* Update the recursion. */
    mpz_set(km2, km1);
    mpz_set(km1, denominator);

    /* Update the fraction: The next fraction is 1 / ( f - floor(f) ). */
    mpfr_sub(fraction, fraction, integer_part, MPFR_RNDN);
    if (mpfr_cmp_ui(fraction, 0) == 0) {
      break; /* The best approximation has been reached. */
    }
    mpfr_div(fraction, one, fraction, MPFR_RNDN);
  }

  /* Clean up. */
  mpfr_clear(one);
  mpz_clear(pow2m);
  mpz_clear(pow2ml);

  mpfr_clear(fraction);
  mpfr_clear(integer_part);

  mpz_clear(denominator);

  mpz_clear(km1);
  mpz_clear(km2);

  mpz_clear(tmp);

  /* Signal result. */
  return found;
}
