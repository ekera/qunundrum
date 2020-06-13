/*!
 * \file    keccak.h
 * \ingroup keccak
 *
 * \brief   The declaration of the Keccak-f permutation.
 *
 * \remark  This implementation of Keccak-f has not been formally evaluated. It
 *          is not fit for, and must not be used for, cryptographic purposes.
 *
 * This implementation is inspired by the Keccak team's reference implementation
 * and by the excellent tiny SHA-3 implementation by Saarinen.
 */

/*!
 * \defgroup keccak The Keccak-f permutation
 * \ingroup  keccak_random
 *
 * \brief    A module for the Keccak-f permutation.
 */

#ifndef KECCAK_H
#define KECCAK_H

#include <stdint.h>

/*!
 * \brief   The number of lanes in a state for Keccak-f.
 */
#define KECCAK_LANE_COUNT     25

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief   Applies the Keccak-f permutation to a state.
 *
 * The state is a 5 x 5 matrix with a 64 bit lane in each position.
 *
 * \param[in, out] lanes    The 25 lanes to which to apply Keccak-f in place.
 */
void keccak_f(
  uint64_t * const lanes);

#ifdef __cplusplus
}
#endif

#endif /* KECCAK_H */
