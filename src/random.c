/*!
 * \file    random.c
 * \ingroup random
 *
 * \brief   The definition of functions for deterministic and non-deterministic
 *          random number generation.
 *
 * \remark  This random number generator has not been formally evaluated. It is
 *          not fit for, and must not used for, any cryptographic purposes.
 */

#include "random.h"

#include "keccak_random.h"
#include "errors.h"

#include "debug_common.h"

#include <gmp.h>

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

/*!
 * \brief   A random initialization canary value.
 */
#define RANDOM_CANARY           0x9635dddb82658d7eULL

/*!
 * \brief   The buffer size in bytes to use for generating random integers.
 */
#define RANDOM_MPZ_BUFFER_SIZE  16384

void random_init(
  Random_State * const state)
{
  memset(state, 0, sizeof(Random_State));

  keccak_random_init(&(state->keccak_state));

  state->random_device = NULL;

  state->canary = RANDOM_CANARY;
  
  #ifdef DEBUG_TRACE_RNG
  printf("random_init(): Debug: Finished initializing state: %p\n", state);
  #endif
}

void random_init_device(
  Random_State * const state,
  const char * const path)
{
  memset(state, 0, sizeof(Random_State));

  state->random_device = fopen(path, "rb");
  if (NULL == state->random_device) {
    critical("random_init_device(): "
      "Failed to open \"%s\" for reading.", path);
  }

  state->canary = RANDOM_CANARY;
    
  #ifdef DEBUG_TRACE_RNG
  printf("random_init_device(): "
    "Debug: Initialized state (with device \"%s\"): %p\n", path, state);
  #endif
}

void random_close(
  Random_State * const state)
{
  if (NULL != state->random_device) {
    fclose(state->random_device);
    state->random_device = NULL;
  } else {
    keccak_random_close(&(state->keccak_state));
  }

  memset(state, 0, sizeof(Random_State));

  #ifdef DEBUG_TRACE_RNG
  printf("random_close(): Debug: Closed state: %p\n", state);
  #endif
}

void random_generate(
  void * dst,
  const uint32_t length,
  Random_State * const state)
{
  if (RANDOM_CANARY != state->canary) {
    critical("random_generate(): The state is not initialized.");
  }

  if (0 == length) {
    return; /* Proceeding to read zero bytes below would produce an error. */
  }

  if (NULL == state->random_device) {
    keccak_random_generate(dst, length, &(state->keccak_state));
  } else {
    if (1 != fread(dst, length, 1, state->random_device)) {
      critical("random_generate(): Failed to read random data.");
    }
  }

  #ifdef DEBUG_TRACE_RNG
  printf("random_generate(): "
    "Debug: Generated %u byte(s) from state %p:\n", length, state);
  debug_print_buffer(dst, length);
  #endif
}

long double random_generate_pivot_inclusive(
  Random_State * const state)
{
  uint64_t value = 0;

  random_generate(&value, sizeof(uint64_t), state);
  
  long double result;
  
  result  = (long double)(value); /* on [0, 2^64) */
  result /= (long double)(0xffffffffffffffffLLU); /* divide by 2^64 - 1 */

  #ifdef DEBUG_TRACE_RNG
  printf("random_generate_pivot_inclusive(): "
    "Debug: Sampled pivot (from value = %llx): %Lf\n",
      (long long unsigned int)value, result);
  #endif

  return result;
}

long double random_generate_pivot_exclusive(
  Random_State * const state)
{
  uint64_t value = 0;

  random_generate(&value, sizeof(uint64_t), state);

  long double result;
  
  result  = (long double)(value & 0x7fffffffffffffffLLU); /* on [0, 2^63) */
  result /= (long double)(0x8000000000000000LLU); /* divide by 2^63 */

  #ifdef DEBUG_TRACE_RNG
  printf("random_generate_pivot_exclusive(): "
    "Debug: Sampled pivot (from value = %llx): %Lf\n",
      (long long unsigned int)value, result);
  #endif

  return result;
}

void random_generate_mpz(
  mpz_t value,
  const mpz_t modulus,
  Random_State * const state)
{
  size_t length = mpz_sizeinbase(modulus, 2) + 64;
  length = (length + 8) / 8;

  if (length > RANDOM_MPZ_BUFFER_SIZE) {
    critical("random_generate_mpz(): Generating a random integer for this "
      "modulus would exceed the buffer capacity.");
  }

  uint8_t buffer[RANDOM_MPZ_BUFFER_SIZE];
  random_generate(buffer, (uint32_t)length, state);

  mpz_import(value, length, 1, 1, 1, 0, buffer);
  mpz_mod(value, value, modulus);

  #ifdef DEBUG_TRACE_RNG
  gmp_printf("random_generate_mpz(): "
    "Debug: Generated integer (for modulus = %Zd): %Zd\n", modulus, value);
  #endif
}