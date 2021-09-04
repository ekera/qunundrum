/*!
 * \file    keccak_random.c
 * \ingroup keccak_random
 *
 * \brief   The definition of functions for a simple deterministic random
 *          number generator based on Keccak-f.
 *
 * \remark  This random number generator has not been formally evaluated. It is
 *          not fit for, and must not be used for, any cryptographic purposes.
 */

#include "keccak_random.h"

#include "errors.h"
#include "keccak.h"

#include "debug_common.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

/*!
 * \brief   A random initialization canary value.
 */
#define KECCAK_RANDOM_CANARY    0xf19d61ae758d6f90ULL

void keccak_random_init(
  Keccak_Random_State * const state)
{
  FILE * file = fopen(KECCAK_RANDOM_DEVICE, "rb");
  if (NULL == file) {
    critical("keccak_random_init(): "
      "Failed to open \"%s\" for reading.", KECCAK_RANDOM_DEVICE);
  }

  uint8_t seed[KECCAK_RANDOM_SEED_LENGTH];
  if (1 != fread(seed, KECCAK_RANDOM_SEED_LENGTH, 1, file)) {
    critical("keccak_random_init(): "
      "Failed to read seed from \"%s\".", KECCAK_RANDOM_DEVICE);
  }

  fclose(file);
  file = NULL;

  #ifdef DEBUG_TRACE_RNG
  printf("keccak_random_init(): "
    "Debug: Read seed (from \"%s\" to initialize state %p):\n",
      KECCAK_RANDOM_DEVICE, state);
  debug_print_buffer(seed, KECCAK_RANDOM_SEED_LENGTH);
  #endif

  keccak_random_init_seed(state, seed);

  #ifdef DEBUG_TRACE_RNG
  printf("keccak_random_init(): "
    "Debug: Finished initializing state: %p\n", state);
  #endif
}

void keccak_random_init_seed(
  Keccak_Random_State * const state,
  const uint8_t * seed)
{
  for (uint32_t i = 0; i < KECCAK_RANDOM_SEED_LENGTH; i++) {
    state->seed[i] = seed[i];
  }

  for (uint32_t i = 0; i < KECCAK_LANE_COUNT; i++) {
    state->lanes[i] = 0;
  }

  for (uint32_t i = 0; i < KECCAK_RANDOM_SEED_LENGTH; i++) {
    state->lanes[i / 8] <<= 8;
    state->lanes[i / 8] ^= seed[i];
  }

  keccak_f(state->lanes);

  state->offset = KECCAK_RANDOM_SEED_LENGTH;
  state->canary = KECCAK_RANDOM_CANARY;

  #ifdef DEBUG_TRACE_RNG
  printf("keccak_random_init_seed(): "
    "Debug: Initialized state %p with seed:\n", state);
  debug_print_buffer(seed, KECCAK_RANDOM_SEED_LENGTH);
  #endif
}

void keccak_random_close(
  Keccak_Random_State * const state)
{
  memset(state, 0, sizeof(Keccak_Random_State));

  #ifdef DEBUG_TRACE_RNG
  printf("keccak_random_init(): Debug: Closed state: %p\n", state);
  #endif
}

void keccak_random_generate(
  uint8_t * dst,
  const uint32_t length,
  Keccak_Random_State * const state)
{
  if (KECCAK_RANDOM_CANARY != state->canary) {
    critical("keccak_random_generate(): The state is not initialized.");
  }

  uint32_t offset = state->offset;

  for (uint32_t i = 0; i < length; i++, offset++) {
    if (offset >= (8 * KECCAK_LANE_COUNT)) {
      keccak_f(state->lanes);
      offset = KECCAK_RANDOM_SEED_LENGTH;
    }

    const uint64_t lane = state->lanes[offset / 8];
    dst[i] = (uint8_t)((lane >> (56 - (8 * (offset % 8)))) & 0xff);
  }

  /* Update the offset in the state. */
  state->offset = offset;

  #ifdef DEBUG_TRACE_RNG
  printf("keccak_random_generate(): "
    "Debug: Generated %u byte(s) from state %p:\n", length, state);
  debug_print_buffer(dst, length);
  #endif
}
