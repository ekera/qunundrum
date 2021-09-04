/*!
 * \file    gmp_mpi.c
 * \ingroup mpi
 *
 * \brief   The definition of functions for sending and receiving arbitrary
 *          precision integers over the message passing interface.
 */

#include "gmp_mpi.h"

#include "common.h"
#include "errors.h"

#include <gmp.h>

#include <mpi.h>

#include <stdint.h>
#include <stdlib.h>

/*!
 * \brief   An MPI tag used to send the length of the serialized representation
 *          of an arbitrary precision integer.
 */
#define MPI_MPZ_LENGTH_TAG                        11011

/*!
 * \brief   An MPI tag used to send the the serialized representation of an
 *          arbitrary precision integer.
 */
#define MPI_MPZ_DATA_TAG                          11012

void mpz_send(
  const mpz_t value,
  const int rank)
{
  if (mpz_cmp_ui(value, 0) < 0) {
    critical("mpz_send(): Negative GMP integers are not supported.");
  }

  size_t length;
  uint8_t * buffer = (uint8_t *)mpz_export(NULL, &length, 1, 1, 1, 0, value);

  uint32_t length_u32 = (uint32_t)length;

  if (MPI_SUCCESS != MPI_Send(
      &length_u32,
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_MPZ_LENGTH_TAG,
      MPI_COMM_WORLD))
  {
    critical("mpz_send(): Failed to send GMP integer.");
  }

  if (length > 0) {
    if (MPI_SUCCESS != MPI_Send(
        buffer,
        length,
        MPI_BYTE,
        rank,
        MPI_MPZ_DATA_TAG,
        MPI_COMM_WORLD))
    {
      critical("mpz_send(): Failed to send GMP integer.");
    }
  }

  /* Clear memory. */
  free(buffer);
  buffer = NULL;
}

void mpz_recv(
  mpz_t value,
  const int rank)
{
  uint32_t length;

  MPI_Status status;

  if (MPI_SUCCESS != MPI_Recv(
      &length,
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_MPZ_LENGTH_TAG,
      MPI_COMM_WORLD,
      &status))
  {
    critical("mpz_recv(): Failed to receive GMP integer.");
  }

  if (0 == length) {
    /* If the value is zero then the value of the operand is zero. */
    mpz_set_ui(value, 0);
    return;
  }

  uint8_t * buffer = malloc(length * sizeof(uint8_t));
  if (NULL == buffer) {
    critical("mpz_recv(): Failed to allocate memory.");
  }

  if (MPI_SUCCESS != MPI_Recv(
      buffer,
      length,
      MPI_BYTE,
      rank,
      MPI_MPZ_DATA_TAG,
      MPI_COMM_WORLD,
      &status))
  {
    critical("mpz_recv(): Failed to receive GMP integer.");
  }

  mpz_import(value, length, 1, 1, 1, 0, buffer);

  /* Clear memory. */
  free(buffer);
  buffer = NULL;
}

void mpz_bcast_send(
  const mpz_t value,
  const int root)
{
  if (mpz_cmp_ui(value, 0) < 0) {
    critical("mpz_bcast_send(): Negative GMP integers are not supported.");
  }

  size_t length;
  uint8_t * buffer = (uint8_t *)mpz_export(NULL, &length, 1, 1, 1, 0, value);
  if (NULL == buffer) {
    critical("mpz_bcast_send(): Failed to allocate memory.");
  }

  if (MPI_SUCCESS != MPI_Bcast(
      &length,
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("mpz_bcast_send(): Failed to send GMP integer.");
  }

  if (length > 0) {
    if (MPI_SUCCESS != MPI_Bcast(
        buffer,
        length,
        MPI_BYTE,
        root,
        MPI_COMM_WORLD))
    {
      critical("mpz_bcast_send(): Failed to send GMP integer.");
    }
  }

  /* Clear memory. */
  free(buffer);
  buffer = NULL;
}

void mpz_bcast_recv(
  mpz_t value,
  const int root)
{
  uint32_t length;

  if (MPI_SUCCESS != MPI_Bcast(
      &length,
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("mpz_bcast_recv(): Failed to receive GMP integer.");
  }

  if (0 == length) {
    /* If the value is zero then the value of the operand is zero. */
    mpz_set_ui(value, 0);
    return;
  }

  uint8_t * buffer = malloc(length * sizeof(uint8_t));
  if (NULL == buffer) {
    critical("mpz_bcast_recv(): Failed to allocate memory.");
  }

  if (MPI_SUCCESS != MPI_Bcast(
      buffer,
      length,
      MPI_BYTE,
      root,
      MPI_COMM_WORLD))
  {
    critical("mpz_bcast_recv(): Failed to receive GMP integer.");
  }

  mpz_import(value, length, 1, 1, 1, 0, buffer);

  /* Clear memory. */
  free(buffer);
  buffer = NULL;
}
