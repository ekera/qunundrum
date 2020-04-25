/*!
 * \file    rotate.h
 * \ingroup keccak
 *
 * \brief   The definition the macros for cyclically rotating unsigned integers.
 */

#ifndef ROTATE_H
#define ROTATE_H

/*!
 * \brief   Rotates a 64 bit integer s steps to the left cyclically.
 */
#define ROTATE_LEFT_UINT64(x, s) (((x) << (s)) | ((x) >> (64 - (s))));

/*!
 * \brief   Rotates a 64 bit integer s steps to the right cyclically.
 */
#define ROTATE_RIGHT_UINT64(x, s) (((x) >> (s)) | ((x) << (64 - (s))));

#endif /* ROTATE_H */
