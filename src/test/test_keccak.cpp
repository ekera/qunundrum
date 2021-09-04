/*!
 * \file    test/test_keccak.cpp
 * \ingroup unit_tests_keccak
 *
 * \brief   The definition of unit tests for \ref keccak module.
 */

#include "test_keccak.h"

#include "../errors.h"
#include "../keccak.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

void test_keccak_f() {
  printf("Testing keccak...\n");

  uint64_t state[KECCAK_LANE_COUNT];
  memset(state, 0, sizeof(state));

  keccak_f(state);

  for (uint32_t i = 0; i < KECCAK_LANE_COUNT; i++) {
    if (state[i] != KECCAK_EXPECTED_FIRST[i]) {
      critical("test_keccak_f(): "
        "Expected %.16llx but found %.16llx at position %u.",
          KECCAK_EXPECTED_FIRST[i], state[i], i);
    }
  }

  keccak_f(state);

  for (uint32_t i = 0; i < KECCAK_LANE_COUNT; i++) {
    if (state[i] != KECCAK_EXPECTED_SECOND[i]) {
      critical("test_keccak_f(): "
        "Expected %.16llx but found %.16llx at position %u.",
          KECCAK_EXPECTED_SECOND[i], state[i], i);
    }
  }
}

void test_keccak() {
  test_keccak_f();
}
