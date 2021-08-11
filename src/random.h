/*!
 * \file    random.h
 * \ingroup random
 *
 * \brief   The declaration of functions for deterministic and non-deterministic
 *          random number generation, and of data structures for random number
 *          generation.
 *
 * \remark  This random number generator has not been formally evaluated. It is
 *          not fit for, and must not used for, any cryptographic purposes.
 */

/*!
 * \defgroup  random Random number generation
 *
 * \brief A module for random number generation.
 */

#ifndef RANDOM_H
#define RANDOM_H

#include "keccak_random.h"

#include <gmp.h>

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief   A data structure representing the state of a random number
 *          generator that may either read randomness from file or expand
 *          a seed using Keccak.
 */
typedef struct {
  /*!
   * \brief   The state of the Keccak-based random number generator.
   */
  Keccak_Random_State keccak_state;

  /*!
   * \brief   A handle to a file from which to read random data.
   *
   * Set to NULL if a seed is instead to be read and expanded using the
   * Keccak-based deterministic random number generator.
   */
  FILE * random_device;

  /*!
   * \brief   An initialization canary.
   */
  uint64_t canary;
} Random_State;

/*!
 * \name Initialization and closing
 * \{
 */

/*!
 * \brief   Initializes a random generator state by reading a seed from
 *          #KECCAK_RANDOM_DEVICE and using it to setup a Keccak-based random
 *          number generator state.
 *
 * This function calls keccak_random_init().
 *
 * \param[in, out] state  The state to be initialized.
 */
void random_init(
  Random_State * const state);

/*!
 * \brief   Initializes a random number generator state by opening a file or
 *          random device from which raw random data will be read.
 *
 * \param[in, out] state  The state to be initialized.
 * \param[in] path        The path to the file or device to read.
 */
void random_init_device(
  Random_State * const state,
  const char * const path);

/*!
 * \brief   Closes the random number generator state.
 */
void random_close(
  Random_State * const state);

/*!
 * \}
 */

/*!
 * \name Generation
 * \{
 */

/*!
 * \brief   Generates output from an initialized random number generator state.
 *
 * This function will either read from the Keccak-based deterministic random
 * generator state, or from a file or device, depending on how the state was
 * initialized, see random_init() and random_init_device().
 *
 * \param[in, out] dst    The destination buffer.
 * \param[in] length      The number of bytes of random data to generate and
 *                        write to the destination buffer.
 * \param[in, out] state  The initialized random number generator state from
 *                        which to generate random data.
 */
void random_generate(
  void * dst,
  const uint32_t length,
  Random_State * const state);

/*!
 * \brief   Generates a coarsely uniformly sampled long double on [0, 1].
 *
 * This function generates an integer on [0, 2^64) uniformly at random and
 * divides the resulting integer by 2^64 - 1 to construct the long double.
 *
 * This function calls random_generate() to generate raw random data.
 *
 * \param[in, out] state  The initialized random number generator state from
 *                        which to generate random data.
 *
 * \return A coarsely uniformly sampled random long double on [0, 1].
 */
long double random_generate_pivot_inclusive(
  Random_State * const state);

/*!
 * \brief   Generates a coarsely uniformly sampled long double on [0, 1).
 *
 * This function generates an integer on [0, 2^63) uniformly at random and
 * divides the resulting integer by 2^63 to construct the long double.
 *
 * This function calls random_generate() to generate raw random data.
 *
 * \param[in, out] state  The initialized random number generator state from
 *                        which to generate random data.
 *
 * \return A coarsely uniformly sampled random long double on [0, 1).
 */
long double random_generate_pivot_exclusive(
  Random_State * const state);

/*!
 * \brief   Generates a uniformly random integer on [0, modulus).
 *
 * Let 2^(m-1) <= modulus < 2^m. This function then generates an integer on
 * [0, 2^(m + 64 + c)) uniformly at random and reduces it by the modulus,
 * for c a constant on (0, 8] such that m + 64 + t is a multiple of 8 bits.
 *
 * This function calls random_generate() to generate raw random data.
 *
 * \param[in, out] value  The destination value.
 * \param[in] modulus     The modulus.
 * \param[in, out] state  The initialized random number generator state from
 *                        which to generate random data.
 */
void random_generate_mpz(
  mpz_t value,
  const mpz_t modulus,
  Random_State * const state);

/*!
 * \}
 */

#ifdef __cplusplus
}
#endif

#endif /* RANDOM_H */
