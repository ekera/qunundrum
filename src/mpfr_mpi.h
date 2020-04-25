/*!
 * \file    mpfr_mpi.h
 * \ingroup mpi
 *
 * \brief   The declaration of functions for sending and receiving arbitrary
 *          precision floating point values over the message passing interface.
 */

#ifndef MPFR_MPI_H
#define MPFR_MPI_H

#include <mpfr.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief   Sends an arbitrary precision integer to a node.
 *
 * \param[in] value       The floating point value to send.
 * \param[in] rank        The rank of the destination node to which to send.
 */
void mpfr_send(
  const mpfr_t value,
  const int rank);

/*!
 * \brief   Receives an arbitrary precision integer from a node.
 *
 * \param[in, out] value  The floating point value to be set to the value
 *                        received.
 * \param[in] rank        The rank of the source node from which to receive.
 */
void mpfr_recv(
  mpfr_t value,
  const int rank);

#ifdef __cplusplus
}
#endif

#endif /* MPFR_MPI_H */
