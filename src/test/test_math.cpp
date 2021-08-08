/*!
 * \file    test/test_math.cpp
 * \ingroup unit_tests_math
 * 
 * \brief   The definition of unit tests for the mathematical functions in the 
 *          \ref math module.
 */

#include "test_math.h"

#include "../math.h"
#include "../errors.h"

#include <stdint.h>
#include <string.h>
#include <stdio.h>

void test_math_mod_reduce() {
  printf("Testing mod_reduce...\n");

  mpz_t x;
  mpz_init(x);

  mpz_t n;
  mpz_init(n);

  /* Tests for n = 7. */
  static const int32_t expected_7[7] = {0, 1, 2, 3, -3, -2, -1};

  mpz_set_ui(n, 7);
  
  for (int32_t i = -3 * 7; i <= 3 * 7; i++) {
    mpz_set_si(x, i);
    mod_reduce(x, n);

    /* Let j = i mod 7 on [0, 7). */
    int32_t j = i;
    while (j <  0) { j += 7; };
    while (j >= 7) { j -= 7; };

    if (0 != mpz_cmp_si(x, expected_7[j])) {
      critical("test_math_mod_reduce(): Comparison failed. "
        "Found %d but expected %d.", mpz_get_si(x), expected_7[j]);
    }
  }

  /* Tests for n = 8. */
  static const int32_t expected_8[8] = {0, 1, 2, 3, -4, -3, -2, -1};
  
  mpz_set_ui(n, 8);

  for (int32_t i = -3 * 8; i <= 3 * 8; i++) {
    mpz_set_si(x, i);
    mod_reduce(x, n);

    /* Let j = i mod 8 on [0, 8). */
    int32_t j = i;
    while (j <  0) { j += 8; };
    while (j >= 8) { j -= 8; };

    if (0 != mpz_cmp_si(x, expected_8[j])) {
      critical("test_math_mod_reduce(): Comparison failed. "
        "Found %d but expected %d.", mpz_get_si(x), expected_8[j]);
    }
  }

  /* Clear memory. */
  mpz_clear(x);
  mpz_clear(n);
}

void test_math() {
  test_math_mod_reduce();
}
