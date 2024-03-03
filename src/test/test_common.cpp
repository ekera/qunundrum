/*!
 * \file    test_common.cpp
 * \ingroup unit_tests
 *
 * \brief   The definition of common functions used by the unit tests.
 */

#include "test_common.h"

#include "../common.h"
#include "../errors.h"
#include "../math.h"

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>
#include <stdio.h>

/*!
 * \brief   The maximum buffer size in bytes.
 */
#define MAX_BUFFER_SIZE 32768

void test_mpfr_load(
  mpfr_t value,
  FILE * const file)
{
  char buffer[MAX_BUFFER_SIZE];

  for (uint32_t i = 0; i < MAX_BUFFER_SIZE; i++) {
    if (1 != fread(&buffer[i], 1, 1, file)) {
      critical("test_mpfr_load(): Failed to read from file.");
    }

    if ('\n' == buffer[i]) {
      buffer[i] = 0;
      break;
    }

    if ((i + 1) >= MAX_BUFFER_SIZE) {
      critical("test_mpfr_load(): Buffer length exceeded.\n");
    }
  }

  mpfr_set_str(value, buffer, 10, MPFR_RNDN);
}

void test_mpz_load(
  mpz_t value,
  FILE * const file)
{
  char buffer[MAX_BUFFER_SIZE];

  for (uint32_t i = 0; i < MAX_BUFFER_SIZE; i++) {
    if (1 != fread(&buffer[i], 1, 1, file)) {
      critical("test_mpz_load(): Failed to read from file.");
    }

    if ('\n' == buffer[i]) {
      buffer[i] = 0;
      break;
    }

    if ((i + 1) >= MAX_BUFFER_SIZE) {
      critical("test_mpz_load(): Buffer length exceeded.\n");
    }
  }

  if (0 != mpz_set_str(value, buffer, 10)) {
    critical("test_mpz_load(): Failed to deserialize integer.\n");
  }
}

bool test_cmp_ld(
  const long double a,
  const long double b,
  const long double tolerance)
{
  if ((a < 0) || (b < 0)) {
    critical("test_cmp_ld(): At least one value is negative.");
  }

  if (a == b) {
    return TRUE;
  }

  if (abs_ld(a - b) / min_ld(abs_ld(a), abs_ld(b)) < tolerance) {
    return TRUE;
  } else {
    return FALSE;
  }
}
