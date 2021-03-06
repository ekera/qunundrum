/*!
 * \file    diagonal_distribution_mpi.cpp
 * \ingroup diagonal_distribution
 *
 * \brief   The definition of functions for sending, receiving and broadcasting
 *          diagonal probability distributions.
 */

#include "diagonal_distribution.h"

#include "diagonal_distribution_slice.h"
#include "parameters.h"
#include "errors.h"

#include <mpi.h>
#include <stdint.h>
#include <stdlib.h>

void diagonal_distribution_init_bcast_recv(
  Diagonal_Distribution * const distribution,
  const int root)
{
  /* Broadcast the parameters. */
  parameters_init(&(distribution->parameters));
  parameters_bcast_recv(&(distribution->parameters), root);

  /* Broadcast the precision. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(distribution->precision),
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_init_bcast_recv(): "
      "Failed to broadcast precision.");
  }

  /* Broadcast the flags. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(distribution->flags),
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_init_bcast_recv(): "
      "Failed to broadcast flags.");
  }

  /* Broadcast the capacity. */
  uint32_t capacity;

  if (MPI_SUCCESS != MPI_Bcast(
      &capacity,
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_init_bcast_recv(): "
      "Failed to broadcast capacity.");
  }

  /* Store the capacity and zeroize the count. */
  distribution->capacity = capacity;
  distribution->count = 0;

  /* Allocate space for slices according to the capacity. */
  distribution->slices =
    (Diagonal_Distribution_Slice **)malloc(
        sizeof(Diagonal_Distribution_Slice *) * capacity);
  if (NULL == distribution->slices) {
    critical("diagonal_distribution_init_bcast_recv(): "
      "Failed to allocate memory.");
  }

  for (uint32_t i = 0; i < capacity; i++) {
    distribution->slices[i] = NULL;
  }

  /* Broadcast the number of slices. */
  uint32_t count;

  if (MPI_SUCCESS != MPI_Bcast(
      &count,
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_init_bcast_recv(): "
      "Failed to broadcast count.");
  }

  /* Broadcast the slices. */
  distribution->total_probability = 0;
  distribution->total_error = 0;

  for (uint32_t i = 0; i < count; i++) {
    /* Broadcast slice. */
    Diagonal_Distribution_Slice * slice = diagonal_distribution_slice_alloc();
    diagonal_distribution_slice_init_bcast_recv(slice, root);

    /* Insert the slice into the distribution. */
    diagonal_distribution_insert_slice(distribution, slice);
  }

  /* Note: The total_probability and total_error will now have been set. */
}

void diagonal_distribution_bcast_send(
  const Diagonal_Distribution * const distribution,
  const int root)
{
  /* Broadcast the parameters. */
  parameters_bcast_send(&(distribution->parameters), root);

  /* Broadcast the precision. */
  uint32_t precision = distribution->precision;

  if (MPI_SUCCESS != MPI_Bcast(
      &precision,
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_bcast_send(): "
      "Failed to broadcast precision.");
  }

  /* Broadcast the flags. */
  uint32_t flags = distribution->flags;

  if (MPI_SUCCESS != MPI_Bcast(
      &flags,
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_bcast_send(): "
      "Failed to broadcast flags.");
  }

  /* Broadcast the capacity. */
  uint32_t capacity = distribution->capacity;

  if (MPI_SUCCESS != MPI_Bcast(
      &capacity,
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_bcast_send(): "
      "Failed to broadcast capacity.");
  }

  /* Broadcast the number of slices. */
  uint32_t count = distribution->count;

  if (MPI_SUCCESS != MPI_Bcast(
      &count,
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_bcast_send(): Failed to broadcast count.");
  }

  /* Broadcast the slices. */
  for (uint32_t i = 0; i < distribution->count; i++) {
    /* Broadcast slice. */
    diagonal_distribution_slice_bcast_send(distribution->slices[i], root);
  }
}
