/*!
 * \file    test/test_random.cpp
 * \ingroup unit_tests_random
 * 
 * \brief   The definition of unit tests for the \ref random module.
 */

#include "test_random.h"

#include "../keccak_random.h"
#include "../random.h"
#include "../errors.h"

#include <gmp.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/*!
 * \brief   The path to a temporary file used to hold random data.
 */
const static char * TEMPORARY_RANDOM_FILE = "/tmp/tmp-random.bin";

/*!
 * \brief   The path to a random device.
 */
const static char * RANDOM_DEVICE = "/dev/urandom";

/*!
 * \brief   The random buffer size in byte.
 */
#define RANDOM_BUFFER_SIZE  65536

void test_random_device() {
  printf("Testing random number generation via file or device...\n");

  /* Read from the random device. */
  FILE * random_device = fopen(RANDOM_DEVICE, "rb");
  if (NULL == random_device) {
    critical("Failed to open \"%s\" for reading.", RANDOM_DEVICE);
  }

  /* Read data into the buffer. */
  uint8_t * buffer = (uint8_t *)malloc(RANDOM_BUFFER_SIZE);
  if (NULL == buffer) {
    critical("Failed to allocate memory.");
  }

  memset(buffer, 0, RANDOM_BUFFER_SIZE);

  if (1 != fread(buffer, RANDOM_BUFFER_SIZE, 1, random_device)) {
    critical("Failed to read from \"%s\".", RANDOM_DEVICE);
  }

  fclose(random_device);
  random_device = NULL;

  /* Setup a random temporary file. */
  FILE * file = fopen(TEMPORARY_RANDOM_FILE, "wb");
  if (NULL == file) {
    critical("Failed to open \"%s\" for writing.", TEMPORARY_RANDOM_FILE);
  }

  if (1 != fwrite(buffer, RANDOM_BUFFER_SIZE, 1, file)) {
    critical("Failed to write to \"%s\".", TEMPORARY_RANDOM_FILE);
  }

  fclose(file);
  file = NULL;
  
  /* Setup a random state for reading back the file in one go. */
  Random_State random_state;
  random_init_device(&random_state, TEMPORARY_RANDOM_FILE);
  if (NULL == random_state.random_device) {
    critical("Failed to initialize the random state for reading a device.");
  }

  uint8_t * buffer2 = (uint8_t *)malloc(RANDOM_BUFFER_SIZE);
  if (NULL == buffer2) {
    critical("Failed to allocate memory.");
  }
  
  memset(buffer2, 0, RANDOM_BUFFER_SIZE);

  random_generate(buffer2, RANDOM_BUFFER_SIZE, &random_state);

  if (0 != memcmp(buffer, buffer2, RANDOM_BUFFER_SIZE)) {
    critical("Failed to correctly generate random output.");
  }
  
  random_close(&random_state);

  /* Setup a random state for reading back the file in increments. */
  random_init_device(&random_state, TEMPORARY_RANDOM_FILE);
  if (NULL == random_state.random_device) {
    critical("Failed to initialize the random state for reading a device.");
  }
  
  uint32_t offset = 0;

  memset(buffer2, 0, RANDOM_BUFFER_SIZE);
  
  Keccak_Random_State chunk_length_state;
  keccak_random_init(&chunk_length_state);

  while (offset < RANDOM_BUFFER_SIZE) {
    uint8_t chunk_length;
    keccak_random_generate(&chunk_length, 1, &chunk_length_state);

    chunk_length %= 32;
    chunk_length %= (RANDOM_BUFFER_SIZE + 1 - offset);

    random_generate(&buffer2[offset], chunk_length, &random_state);
    offset += chunk_length;
  }

  if (0 != memcmp(buffer, buffer2, RANDOM_BUFFER_SIZE)) {
    critical("Failed to correctly generate random output.");
  }
  
  keccak_random_close(&chunk_length_state);

  random_close(&random_state);

  /* Clear memory. */
  free(buffer);
  buffer = NULL;

  free(buffer2);
  buffer2 = NULL;
}

void test_random_keccak() {
  printf("Testing random number generation via keccak...\n");

  /* Setup a random state for reading back the file in one go. */
  Random_State random_state;
  random_init(&random_state);
  if (NULL != random_state.random_device) {
    critical("Failed to initialize the random state for using Keccak.");
  }

  /* Setup a Keccak state with the same seed. */
  Keccak_Random_State keccak_state;
  keccak_random_init_seed(&keccak_state, random_state.keccak_state.seed);
  
  if (0 != memcmp(&(random_state.keccak_state.lanes), &(keccak_state.lanes), 
    sizeof(KECCAK_LANE_COUNT * sizeof(uint64_t))))
  {
    critical("Failed to initialize the random state for using Keccak.");
  }

  if (0 != memcmp(&(random_state.keccak_state.seed), &(keccak_state.seed), 
    KECCAK_RANDOM_SEED_LENGTH))
  {
    critical("Failed to initialize the random state for using Keccak.");
  }

  if (random_state.keccak_state.offset != keccak_state.offset) {
    critical("Failed to initialize the random state for using Keccak.");
  }

  /* Check that reading through the wrapper produces the expected result. */
  uint8_t * buffer = (uint8_t *)malloc(RANDOM_BUFFER_SIZE);
  if (NULL == buffer) {
    critical("Failed to allocate memory.");
  }

  memset(buffer, 0, RANDOM_BUFFER_SIZE);

  keccak_random_generate(buffer, RANDOM_BUFFER_SIZE, &keccak_state);

  uint8_t * buffer2 = (uint8_t *)malloc(RANDOM_BUFFER_SIZE);
  if (NULL == buffer2) {
    critical("Failed to allocate memory.");
  }

  memset(buffer2, 0, RANDOM_BUFFER_SIZE);

  random_generate(buffer2, RANDOM_BUFFER_SIZE, &random_state);

  if (0 != memcmp(buffer, buffer2, RANDOM_BUFFER_SIZE)) {
    critical("Failed to correctly generate random output.");
  }

  memset(buffer2, 0, RANDOM_BUFFER_SIZE);

  Keccak_Random_State chunk_length_state;
  keccak_random_init(&chunk_length_state);

  uint32_t offset = 0;

  while (offset < RANDOM_BUFFER_SIZE) {
    uint8_t chunk_length;
    keccak_random_generate(&chunk_length, 1, &chunk_length_state);

    chunk_length %= 32;
    chunk_length %= (RANDOM_BUFFER_SIZE + 1 - offset);

    random_generate(&buffer2[offset], chunk_length, &random_state);
    offset += chunk_length;
  }

  keccak_random_close(&chunk_length_state);

  keccak_random_generate(buffer, RANDOM_BUFFER_SIZE, &keccak_state);

  if (0 != memcmp(buffer, buffer2, RANDOM_BUFFER_SIZE)) {
    critical("Failed to correctly generate random output.");
  }
  
  /* Clear memory. */
  keccak_random_close(&keccak_state);
  random_close(&random_state);

  free(buffer);
  buffer = NULL;

  free(buffer2);
  buffer2 = NULL;
}

void test_random_generate_pivot() {
  printf("Testing random_generate_pivot_*()...\n");

  /* Setup a random temporary file. */
  FILE * file = fopen(TEMPORARY_RANDOM_FILE, "wb");
  if (NULL == file) {
    critical("Failed to open \"%s\" for writing.", TEMPORARY_RANDOM_FILE);
  }

  uint8_t buffer[32] = {
    0xa8, 0xb0, 0x2d, 0xac, 0xc8, 0x79, 0x8f, 0x43,
    0x05, 0x0b, 0x0b, 0x5d, 0x1e, 0x1f, 0x27, 0x9d,
    0xbe, 0xc7, 0x27, 0x5a, 0xe1, 0xe9, 0x46, 0xc9,
    0x57, 0x8a, 0x9f, 0xf1, 0xbd, 0x7b, 0xe4, 0xc0 
  };

  if (1 != fwrite(buffer, 32, 1, file)) {
    critical("Failed to write to \"%s\".", TEMPORARY_RANDOM_FILE);
  }

  fclose(file);
  file = NULL;

  Random_State random_state;
  random_init_device(&random_state, TEMPORARY_RANDOM_FILE);

  long double tmp;

  long double expected;

  tmp = random_generate_pivot_exclusive(&random_state);
  expected = 
    (long double)0x438f79c8ac2db0a8ULL / (long double)0x8000000000000000ULL;
  if (tmp != expected) {
    critical("Failed to sample exclusive (1).");
  }

  tmp = random_generate_pivot_exclusive(&random_state);
  expected = 
    (long double)0x1d271f1e5d0b0b05ULL / (long double)0x8000000000000000ULL;
  if (tmp != expected) {
    critical("Failed to sample exclusive (2).");
  }

  tmp = random_generate_pivot_inclusive(&random_state);
  expected = 
    (long double)0xc946e9e15a27c7beULL / (long double)0xffffffffffffffffULL;
  if (tmp != expected) {
    critical("Failed to sample inclusive (1).");
  }

  tmp = random_generate_pivot_inclusive(&random_state);
  expected = 
    (long double)0xc0e47bbdf19f8a57ULL / (long double)0xffffffffffffffffULL;
  if (tmp != expected) {
    critical("Failed to sample inclusive (2).");
  }

  random_close(&random_state);
}

void test_random_generate_mpz() {
  printf("Testing random_generate_mpz()...\n");

  /* Setup a random temporary file. */
  FILE * file = fopen(TEMPORARY_RANDOM_FILE, "wb");
  if (NULL == file) {
    critical("Failed to open \"%s\" for writing.", TEMPORARY_RANDOM_FILE);
  }

  uint8_t buffer[33] = {
    0xa8, 0xb0, 0x2d, 0xac, 0xc8, 0x79, 0x8f, 0x43,
    0x05, 0x0b, 0x0b, 0x5d, 0x1e, 0x1f, 0x27, 0x9d,
    0xbe, 0xc7, 0x27, 0x5a, 0xe1, 0xe9, 0x46, 0xc9,
    0x57, 0x8a, 0x9f, 0xf1, 0xbd, 0x7b, 0xe4, 0xc0,
    0x87
  };

  if (1 != fwrite(buffer, 33, 1, file)) {
    critical("Failed to write to \"%s\".", TEMPORARY_RANDOM_FILE);
  }

  fclose(file);
  file = NULL;

  Random_State random_state;
  random_init_device(&random_state, TEMPORARY_RANDOM_FILE);

 
  mpz_t modulus;
  mpz_init_set_str(modulus,
    "a6cb8df7f7dbe3de7165ff93552c8edc73b6039daced4d4b", 16);

  mpz_t value;
  mpz_init(value);

  random_generate_mpz(value, modulus, &random_state);

  mpz_t expected;
  mpz_init_set_str(expected,
    "37c26eca1428c2c0234ab97064c9f631325e56856510c764", 16);
  
  if (0 != mpz_cmp(value, expected)) {
    critical("Failed to generate mpz.");
  }

  /* Close the random state. */
  random_close(&random_state);

  /* Clear memory. */
  mpz_clear(modulus);
  mpz_clear(value);
  mpz_clear(expected);
}

void test_random() {
  printf("Testing random...\n");

  test_random_device();
  test_random_keccak();
  test_random_generate_pivot();
  test_random_generate_mpz();
}