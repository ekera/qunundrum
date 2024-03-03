/*!
 * \file    test_common.h
 * \ingroup unit_tests
 *
 * \brief   The declaration of common functions used by the unit tests.
 */

/*!
 * \defgroup unit_tests Unit tests
 * \ingroup  unit_tests_debug_tracing
 *
 * \brief    A module for unit tests.
 */

#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <gmp.h>
#include <mpfr.h>

#include <stdio.h>

/*!
 * \brief   Loads a multiprecision floating point value from a file.
 *
 * \param[in, out] value    The variable in which to store the value parsed.
 * \param[in, out] file     The file from which to read the value.
 */
void test_mpfr_load(
  mpfr_t value,
  FILE * const file);

/*!
 * \brief   Loads a multiprecision integer value from a file.
 *
 * \param[in, out] value    The variable in which to store the value parsed.
 * \param[in, out] file     The file from which to read the value.
 */
void test_mpz_load(
  mpz_t value,
  FILE * const file);

/*!
 * \brief   Compares two long doubles for equality.
 *
 * Returns #TRUE if a is equal to b, or if abs(a - b) / min(a, b) < tolerance,
 * returns #FALSE otherwise.
 *
 * This function fails with a critical error if a and b are not both positive.
 *
 * \param[in] a   The long double a.
 * \param[in] b   The long double b.
 *
 * \param[in] tolerance   The tolerance. Defaults to 1e-6.
 *
 * \return  #TRUE if the value are equal or within tolerance, #FALSE otherwise.
 */
bool test_cmp_ld(
  const long double a,
  const long double b,
  const long double tolerance = 1e-6);

#endif /* TEST_COMMON_H */
