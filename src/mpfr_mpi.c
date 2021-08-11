/*!
 * \file    mpfr_mpi.c
 * \ingroup mpi
 *
 * \brief   The definition of functions for sending and receiving arbitrary
 *          precision floating point values over the message passing interface.
 */

#include "mpfr_mpi.h"

#include "common.h"
#include "errors.h"

#include <mpfr.h>

#include <mpi.h>

/*!
 * \brief   An MPI tag used to send the length of the serialized representation
 *          of an arbitrary precision floating point value.
 */
#define MPI_MPFR_LENGTH_TAG                       11021

/*!
 * \brief   An MPI tag used to send the the serialized representation of an
 *          arbitrary precision floating point value.
 */
#define MPI_MPFR_DATA_TAG                         11022

/*!
 * \brief   The length in bytes of the buffer used to serialize arbitrary 
 *          precision floating point integers.
 */
#define MPFR_LENGTH_BUFFER                        32768

void mpfr_send(
  const mpfr_t value,
  const int rank)
{
  char buffer[MPFR_LENGTH_BUFFER];
  mpfr_snprintf(buffer, MPFR_LENGTH_BUFFER, "%.200Rg", value);

  if (MPI_SUCCESS != MPI_Send(
      buffer,
      MPFR_LENGTH_BUFFER, /* count */
      MPI_CHAR,
      rank,
      MPI_MPFR_DATA_TAG,
      MPI_COMM_WORLD))
  {
    critical("mpfr_send(): Failed to send floating point value.");
  }
}

void mpfr_recv(
  mpfr_t value,
  const int rank)
{
  char buffer[MPFR_LENGTH_BUFFER];

  MPI_Status status;

  if (MPI_SUCCESS != MPI_Recv(
      buffer,
      MPFR_LENGTH_BUFFER, /* count */
      MPI_CHAR,
      rank,
      MPI_MPFR_DATA_TAG,
      MPI_COMM_WORLD,
      &status))
  {
    critical("mpfr_recv(): Failed to receive floating point value.");
  }

  mpfr_set_str(value, buffer, 10, MPFR_RNDN);
}
