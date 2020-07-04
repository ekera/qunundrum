/*!
 * \file    gmp_mpi.h
 * \ingroup mpi
 *
 * \brief   The declaration of functions for sending and receiving arbitrary
 *          precision integers over the message passing interface.
 */

/*!
 * \defgroup mpi Message passing
 * \ingroup  utility
 *
 * \brief    A module for conveniency functions for message passing.
 */

#ifndef GMP_MPI_H
#define GMP_MPI_H

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief   Sends an arbitrary precision integer to a node.
 *
 * \param[in] value       The integer to send.
 * \param[in] rank        The rank of the destination node to which to send.
 */
void mpz_send(
  const mpz_t value,
  const int rank);

/*!
 * \brief   Receives an arbitrary precision integer from a node.
 *
 * \param[in, out] value  The integer to be set to the value received.
 * \param[in] rank        The rank of the source node from which to receive.
 */
void mpz_recv(
  mpz_t value,
  const int rank);

/*!
 * \brief   Sends a broadcast of an arbitrary precision integer.
 *
 * \param[in] value       The integer to broadcast.
 * \param[in] root        The rank of the broadcasting node.
 */
void mpz_bcast_send(
  const mpz_t value,
  const int root);

/*!
 * \brief   Receives a broadcast of an arbitrary precision integer.
 *
 * \param[in, out] value  The integer to broadcast.
 * \param[in] root        The rank of the broadcasting node.
 */
void mpz_bcast_recv(
  mpz_t value,
  const int root);

#ifdef __cplusplus
}
#endif

#endif /* GMP_MPI_H */
