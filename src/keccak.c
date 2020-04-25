/*!
 * \file    keccak.c
 * \ingroup keccak
 *
 * \brief   The definition the Keccak-f permutation.
 *
 * \remark  This implementation of Keccak-f has not been formally evaluated. It
 *          is not fit for, and must not used for, any cryptographic purposes.
 *
 * This implementation is inspired by the Keccak team's reference implementation
 * and by the excellent tiny SHA-3 implementation by Saarinen.
 */

#include "rotate.h"

#include <stdint.h>

/*!
 * \brief   The number of rounds in Keccak-f.
 */
#define KECCAK_ROUNDS     24

/*!
 * \brief   The round constants for Keccak-f.
 */
static const uint64_t KECCAK_ROUND_CONSTANTS[KECCAK_ROUNDS] =
{
  0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL,
  0x8000000080008000ULL, 0x000000000000808bULL, 0x0000000080000001ULL,
  0x8000000080008081ULL, 0x8000000000008009ULL, 0x000000000000008aULL,
  0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
  0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL,
  0x8000000000008003ULL, 0x8000000000008002ULL, 0x8000000000000080ULL,
  0x000000000000800aULL, 0x800000008000000aULL, 0x8000000080008081ULL,
  0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL
};

/*!
 * \brief   The rotational constants for Keccak-f.
 */
static const uint32_t KECCAK_ROTATIONAL_CONSTANTS[KECCAK_ROUNDS] =
{
   1,  3,  6, 10, 15, 21, 28, 36, 45, 55,  2, 14,
  27, 41, 56,  8, 25, 43, 62, 18, 39, 61, 20, 44
};

/*!
 * \brief   The lane re-ordering constants for Keccak-f.
 */
static const uint32_t KECCAK_PI_LANE_INDICES[KECCAK_ROUNDS] =
{
  10,  7, 11, 17, 18, 3,  5, 16,  8, 21, 24, 4,
  15, 23, 19, 13, 12, 2, 20, 14, 22,  9,  6, 1
};

void keccak_f(
  uint64_t * const lanes)
{
  /* Temporary variabls. */
  uint64_t tmp_lanes[5];

  uint64_t tmp_lane;

  /* Iterate through the rounds. */
  for (uint32_t r = 0; r < KECCAK_ROUNDS; r++) {
    /* The theta transform. */
    for (uint32_t i = 0; i < 5; i++) {
      tmp_lanes[i] =
        lanes[i] ^ lanes[i + 5] ^ lanes[i + 10] ^ lanes[i + 15] ^ lanes[i + 20];
    }

    for (uint32_t i = 0; i < 5; i++) {
      tmp_lane = tmp_lanes[(i + 4) % 5] ^
        ROTATE_LEFT_UINT64(tmp_lanes[(i + 1) % 5], 1);

      for (uint32_t j = 0; j < 25; j += 5) {
        lanes[j + i] ^= tmp_lane;
      }
    }

    /* The rho and pi transforms combined. */
    tmp_lane = lanes[1];

    for (uint32_t i = 0; i < 24; i++) {
      const uint32_t j = KECCAK_PI_LANE_INDICES[i];

      const uint64_t tmp = lanes[j];
      lanes[j] = ROTATE_LEFT_UINT64(tmp_lane, KECCAK_ROTATIONAL_CONSTANTS[i]);
      tmp_lane = tmp;
    }

    /* The chi transform. */
    for (uint32_t j = 0; j < 25; j += 5) {
      for (uint32_t i = 0; i < 5; i++) {
        tmp_lanes[i] = lanes[j + i];
      }

      for (uint32_t i = 0; i < 5; i++) {
        lanes[j + i] ^= (~tmp_lanes[(i + 1) % 5]) & tmp_lanes[(i + 2) % 5];
      }
    }

    /* The iota transform. */
    lanes[0] ^= KECCAK_ROUND_CONSTANTS[r];
  }
}
