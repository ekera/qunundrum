/*!
 * \file    test/test_math.cpp
 * \ingroup unit_tests_math
 *
 * \brief   The definition of unit tests for the basic mathematical functions
 *          in the \ref math module.
 */

#include "test_math.h"

#include "../errors.h"
#include "../math.h"

#include <gmp.h>

#include <stdint.h>
#include <stdio.h>

void test_abs() {
  printf("Testing abs...\n");

  /* Integers. */
  if ((0 != abs_i(0)) || (2 != abs_i(2)) || (2 != abs_i(-2))) {
    critical("test_abs(): Test failed for abs_i().");
  }

  /* Doubles. */
  if ((0 != abs_d(0)) || (2 != abs_d(2)) || (2 != abs_d(-2))) {
    critical("test_abs(): Test failed for abs_d().");
  }

  /* Long doubles. */
  if ((0 != abs_ld(0)) || (2 != abs_ld(2)) || (2 != abs_ld(-2))) {
    critical("test_abs(): Test failed for abs_ld().");
  }
}

void test_sgn() {
  printf("Testing sgn...\n");

  /* Integers. */
  if ((1 != sgn_i(0)) || (1 != sgn_i(2)) || (-1 != sgn_i(-2))) {
    critical("test_sgn(): Test failed for sgn_i().");
  }

  /* Doubles. */
  if ((1 != sgn_d(0)) || (1 != sgn_d(2)) || (-1 != sgn_d(-2))) {
    critical("test_sgn(): Test failed for sgn_d().");
  }

  /* Long doubles. */
  if ((1 != sgn_ld(0)) || (1 != sgn_ld(2)) || (-1 != sgn_ld(-2))) {
    critical("test_sgn(): Test failed for sgn_ld().");
  }
}

void test_min() {
  printf("Testing min...\n");

  /* Signed integers. */
  if ((0 != min_i(0, 1)) || (0 != min_i(1, 0))) {
    critical("test_min(): Test failed for min_i().");
  }

  if ((-1 != min_i(0, -1)) || (-1 != min_i(-1, 0))) {
    critical("test_min(): Test failed for min_i().");
  }

  if ((-1 != min_i(1, -1)) || (-1 != min_i(-1, 1))) {
    critical("test_min(): Test failed for min_i().");
  }

  /* Unsigned integers. */
  if ((0 != min_ui(0, 1)) || (0 != min_ui(1, 0))) {
    critical("test_min(): Test failed for min_ui().");
  }

  /* Doubles. */
  if ((0 != min_d(0, 1)) || (0 != min_d(1, 0))) {
    critical("test_min(): Test failed for min_d().");
  }

  if ((-1 != min_d(0, -1)) || (-1 != min_d(-1, 0))) {
    critical("test_min(): Test failed for min_d().");
  }

  if ((-1 != min_d(1, -1)) || (-1 != min_d(-1, 1))) {
    critical("test_min(): Test failed for min_d().");
  }

  /* Long doubles. */
  if ((0 != min_ld(0, 1)) || (0 != min_ld(1, 0))) {
    critical("test_min(): Test failed for min_ld().");
  }

  if ((-1 != min_ld(0, -1)) || (-1 != min_ld(-1, 0))) {
    critical("test_min(): Test failed for min_ld().");
  }

  if ((-1 != min_ld(1, -1)) || (-1 != min_ld(-1, 1))) {
    critical("test_min(): Test failed for min_ld().");
  }
}

void test_max() {
  printf("Testing max...\n");

  /* Signed integers. */
  if ((1 != max_i(0, 1)) || (1 != max_i(1, 0))) {
    critical("test_max(): Test failed for max_i().");
  }

  if ((0 != max_i(0, -1)) || (0 != max_i(-1, 0))) {
    critical("test_max(): Test failed for max_i().");
  }

  if ((1 != max_i(1, -1)) || (1 != max_i(-1, 1))) {
    critical("test_max(): Test failed for max_i().");
  }

  /* Unsigned integers. */
  if ((1 != max_ui(0, 1)) || (1 != max_ui(1, 0))) {
    critical("test_max(): Test failed for max_ui().");
  }

  /* Doubles. */
  if ((1 != max_d(0, 1)) || (1 != max_d(1, 0))) {
    critical("test_max(): Test failed for max_d().");
  }

  if ((0 != max_d(0, -1)) || (0 != max_d(-1, 0))) {
    critical("test_max(): Test failed for max_d().");
  }

  if ((1 != max_d(1, -1)) || (1 != max_d(-1, 1))) {
    critical("test_max(): Test failed for max_d().");
  }

  /* Long doubles. */
  if ((1 != max_ld(0, 1)) || (1 != max_ld(1, 0))) {
    critical("test_max(): Test failed for max_ld().");
  }

  if ((0 != max_ld(0, -1)) || (0 != max_ld(-1, 0))) {
    critical("test_max(): Test failed for max_ld().");
  }

  if ((1 != max_ld(1, -1)) || (1 != max_ld(-1, 1))) {
    critical("test_max(): Test failed for max_ld().");
  }
}

void test_math_mod_reduce() {
  printf("Testing mod_reduce...\n");

  mpz_t x;
  mpz_init(x);

  mpz_t n;
  mpz_init(n);

  /* Tests for n = 7. */
  const int32_t expected_7[7] = {0, 1, 2, 3, -3, -2, -1};

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
  const int32_t expected_8[8] = {0, 1, 2, 3, -4, -3, -2, -1};

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

void test_kappa() {
  printf("Testing kappa...\n");

  const uint32_t expected_kappa[9] {
    0, /* = kappa( 0) */
    0, /* = kappa(±1) */
    1, /* = kappa(±2) */
    0, /* = kappa(±3) */
    2, /* = kappa(±4) */
    0, /* = kappa(±5) */
    1, /* = kappa(±6) */
    0, /* = kappa(±7) */
    3  /* = kappa(±8) */
  };

  mpz_t x;
  mpz_init(x);

  for (uint32_t i = 0; i <= 8; i++) {
    mpz_set_ui(x, i);

    if (expected_kappa[i] != kappa(x)) {
      critical("test_kappa(): Failed to correctly compute kappa(%u).", i);
    }

    if (0 != i) {
      mpz_neg(x, x);

      if (expected_kappa[i] != kappa(x)) {
        critical("test_kappa(): Failed to correctly compute kappa(-%u).", i);
      }
    }
  }

  /* Clear memory. */
  mpz_clear(x);
}

void test_is_pow2() {
  printf("Testing is_pow2...\n");

  /* Powers of two. */
  for (uint32_t i = 0; i < 32; i++) {
    const uint32_t x = 1 << i;

    if (!is_pow2(x)) {
      critical("test_is_pow2(): Failed to correctly detect powers of two.");
    }
  }

  /* Non-powers of two. */
  const uint32_t non_powers_of_two[8] = {
    0, 3, 5, 0xffffffff, 0x579e598e, 0xbe6cb578, 0x20830ea6, 0xf71bb829
  };

  for (uint32_t i = 0; i < 8; i++) {
    if (is_pow2(non_powers_of_two[i])) {
      critical("test_is_pow2(): Failed to correctly detect powers of two.");
    }
  }
}

void test_math() {
  test_abs();
  test_sgn();
  test_min();
  test_max();

  test_math_mod_reduce();
  test_kappa();
  test_is_pow2();
}
