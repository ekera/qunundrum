/*!
 * \file    keccak_random.h
 * \ingroup keccak_random
 *
 * \brief   The declaration of functions and data structures for a simple
 *          deterministic random number generator based on Keccak-f.
 *
 * The generator simply sets the first 256 bits of the state to a seed value,
 * zeroes the remainder of the state, and then squeezes the sponge by
 * iteratively applying Keccak-f to generate output.
 *
 * \remark  This random number generator has not been formally evaluated. It is
 *          not fit for, and must not be used for, cryptographic purposes.
 */

/*!
 * \defgroup keccak_random Keccak-based random number generation
 * \ingroup  random
 *
 * \brief    A module for deterministic random generation based on Keccak-f.
 */

#ifndef KECCAK_RANDOM_H
#define KECCAK_RANDOM_H

#include "keccak.h"

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief   The seed length in bytes.
 */
#define KECCAK_RANDOM_SEED_LENGTH   32

/*!
 * \brief   The random device to use to read seeds.
 */
static const char * const KECCAK_RANDOM_DEVICE = "/dev/urandom";


/*!
 * \brief   A data structure representing the state of the deterministic random
 *          number generator based on Keccak-f.
 */
typedef struct {
  /*!
   * \brief   The #KECCAK_LANE_COUNT lanes of the Keccak state.
   */
  uint64_t lanes[KECCAK_LANE_COUNT];

  /*!
   * \brief   The seed used to initialize the state.
   */
  uint8_t seed[KECCAK_RANDOM_SEED_LENGTH];

  /*!
   * \brief   The offset within the state.
   */
  uint32_t offset;

  /*!
   * \brief   An initialization canary.
   */
  uint64_t canary;
} Keccak_Random_State;


/*!
 * \name Initialization
 * \{
 */

/*!
 * \brief   Initializes the random state by reading a seed from
 *          #KECCAK_RANDOM_DEVICE.
 *
 * \param[in, out] state      The state to initialize.
 */
void keccak_random_init(
  Keccak_Random_State * const state);

/*!
 * \brief   Initializes the random state with a seed.
 *
 * \param[in, out] state      The state to initialize.
 * \param[in] seed            A seed of length #KECCAK_RANDOM_SEED_LENGTH bytes.
 */
void keccak_random_init_seed(
  Keccak_Random_State * const state,
  const uint8_t * seed);

/*!
 * \brief   Closes the random state.
 *
 * \param[in, out] state      The state to close.
 */
void keccak_random_close(
  Keccak_Random_State * const state);

/*!
 * \}
 */

/*!
 * \name Generation
 * \{
 */

/*!
 * \brief   Generates output from the random number generator state.
 *
 * \param[in, out] dst        The destination buffer in which to store the
 *                            output generated.
 * \param[in] length          The number of bytes of generate.
 * \param[in, out] state      The state from which to generate output.
 */
void keccak_random_generate(
  uint8_t * dst,
  const uint32_t length,
  Keccak_Random_State * const state);

/*!
 * \}
 */

#ifdef __cplusplus
}
#endif

#endif /* KECCAK_RANDOM_H */
