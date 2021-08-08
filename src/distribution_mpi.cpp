/*!
 * \file    distribution_mpi.cpp
 * \ingroup two_dimensional_distribution
 *
 * \brief   The definition of functions for sending, receiving and broadcasting
 *          two-dimensional probability distributions.
 */

#include "distribution.h"

#include "distribution_slice.h"
#include "lattice_sample.h"
#include "parameters.h"
#include "errors.h"

#include <mpi.h>

#include <stdint.h>
#include <stdlib.h>

void distribution_init_bcast_recv(
  Distribution * const distribution,
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
    critical("distribution_init_bcast_recv(): "
      "Failed to broadcast the precision.");
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
    critical("distribution_init_bcast_recv(): "
      "Failed to broadcast the capacity.");
  }

  /* Store the capacity and zeroize the count. */
  distribution->capacity = capacity;
  distribution->count = 0;

  /* Allocate space for slices according to the capacity. */
  distribution->slices =
    (Distribution_Slice **)malloc(sizeof(Distribution_Slice *) * capacity);
  if (NULL == distribution->slices) {
    critical("distribution_init_bcast_recv(): Failed to allocate memory.");
  }

  for (uint32_t i = 0; i < capacity; i++) {
    distribution->slices[i] = NULL;
  }

  /* Broadcast the slice count. */
  uint32_t count;

  if (MPI_SUCCESS != MPI_Bcast(
      &count,
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_init_bcast_recv(): "
      "Failed to broadcast the slice count.");
  }

  /* Broadcast the slices. */
  distribution->total_probability = 0;
  distribution->total_error = 0;

  for (uint32_t i = 0; i < count; i++) {
    /* Broadcast slice. */
    Distribution_Slice * slice = distribution_slice_alloc();
    distribution_slice_init_bcast_recv(slice, root);

    /* Insert the slice into the distribution. */
    distribution_insert_slice(distribution, slice);
  }

  /* Compute the alpha lattice to prepare for future sampling requests. */
  lattice_alpha_init(
    &(distribution->lattice_alpha),
    &(distribution->parameters));
}

void distribution_bcast_send(
  const Distribution * const distribution,
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
    critical("distribution_bcast_send(): "
      "Failed to broadcast the precision.");
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
    critical("distribution_bcast_send(): "
      "Failed to broadcast the capacity.");
  }

  /* Broadcast the slice count. */
  uint32_t count = distribution->count;

  if (MPI_SUCCESS != MPI_Bcast(
      &count,
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_bcast_send(): "
      "Failed to broadcast the slice count.");
  }

  /* Broadcast the slices. */
  for (uint32_t i = 0; i < distribution->count; i++) {
    /* Broadcast slice. */
    distribution_slice_bcast_send(distribution->slices[i], root);
  }
}
