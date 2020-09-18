/*!
 * \file    test/test_keccak_random.cpp
 * \ingroup unit_tests_keccak_random
 * 
 * \brief   The definition of unit tests for the \ref keccak_random module.
 */

#include "test_keccak_random.h"

#include "test_keccak.h"

#include "../keccak_random.h"
#include "../keccak.h"
#include "../errors.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void test_keccak_random() {
  printf("Testing keccak_random...\n");

  Keccak_Random_State state;

  uint8_t seed[KECCAK_RANDOM_SEED_LENGTH];
  memset(seed, 0, sizeof(seed));
  keccak_random_init_seed(&state, seed);

  if (0 != memcmp(seed, state.seed, KECCAK_RANDOM_SEED_LENGTH)) {
    critical("Incorrect seed stored in the state.");
  }

  /* When initialized with the zero seed, the state after initialization should
   * be equal to that resulting when keccak_f() is applied to the zero block. */
  for (uint32_t i = 0; i < KECCAK_LANE_COUNT; i++) {
    if (state.lanes[i] != KECCAK_EXPECTED_FIRST[i]) {
      critical("Expected %.16llx but found %.16llx at position %u.",
        KECCAK_EXPECTED_FIRST[i], state.lanes[i], i);
    }
  }

  /* If we proceed to generate output, we expect to obtain data equal to lanes 
   * 4 thru 24 in #KECCAK_EXPECTED_FIRST and #KECCAK_EXPECTED_SECOND. */
  uint8_t EXPECTED[2 * (8 * KECCAK_LANE_COUNT - KECCAK_RANDOM_SEED_LENGTH)];

  for (uint32_t i = KECCAK_RANDOM_SEED_LENGTH; i < 8 * KECCAK_LANE_COUNT; i++) {
    EXPECTED[i - KECCAK_RANDOM_SEED_LENGTH] = 
      (uint8_t)((KECCAK_EXPECTED_FIRST[i / 8] >> (56 - 8 * (i % 8))) & 0xff);
    EXPECTED[i + 8 * KECCAK_LANE_COUNT - 2 * KECCAK_RANDOM_SEED_LENGTH] = 
      (uint8_t)((KECCAK_EXPECTED_SECOND[i / 8] >> (56 - 8 * (i % 8))) & 0xff);
  }

  /* Read byte by byte from the state and verify. */
  const uint32_t length = 
    2 * (8 * KECCAK_LANE_COUNT - KECCAK_RANDOM_SEED_LENGTH);

  for (uint32_t i = 0; i < length; i++) {
    uint8_t byte;

    keccak_random_generate(&byte, 1, &state);

    if (EXPECTED[i] != byte) {
      critical("Expected %.2x but found %.2x at position %u.",
        EXPECTED[i], byte, i);
    }
  }

  /* Read the entire chunk in one go and compare. */
  memset(&state, 0, sizeof(state));
  keccak_random_init_seed(&state, seed);

  if (0 != memcmp(seed, state.seed, KECCAK_RANDOM_SEED_LENGTH)) {
    critical("Incorrect seed stored in the state.");
  }

  uint8_t * buffer = (uint8_t *)malloc(length * sizeof(uint8_t));
  if (NULL == buffer) {
    critical("Failed to allocate memory.");
  }

  memset(buffer, 0, length);
  keccak_random_generate(buffer, length, &state);
  if (0 != memcmp(buffer, EXPECTED, length)) {
    critical("The generated and expected sequences differ.");
  }

  /* Read in chunks of varying lengths on [0, 32) and compare. */
  memset(&state, 0, sizeof(state));
  keccak_random_init_seed(&state, seed);

  if (0 != memcmp(seed, state.seed, KECCAK_RANDOM_SEED_LENGTH)) {
    critical("Incorrect seed stored in the state.");
  }
  
  Keccak_Random_State read_state;
  keccak_random_init(&read_state);

  uint32_t offset = 0;

  memset(buffer, 0, length);

  while (offset < length) {
    uint8_t chunk_length;
    keccak_random_generate(&chunk_length, 1, &read_state);

    chunk_length %= 32;
    chunk_length %= (length + 1 - offset);

    keccak_random_generate(&buffer[offset], chunk_length, &state);
    offset += chunk_length;
  }

  if (0 != memcmp(buffer, EXPECTED, length)) {
    critical("The generated and expected sequences differ.");
  }

  /* Clear memory. */
  keccak_random_close(&state);
  keccak_random_close(&read_state);
  
  free(buffer);
  buffer = NULL;
}
