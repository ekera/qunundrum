/*!
 * \file    linear_distribution_slice_mpi.cpp
 * \ingroup linear_distribution_slice
 *
 * \brief   The definition of functions for sending, receiving and broadcasting
 *          slices in linear probability distributions.
 */

#include "linear_distribution_slice.h"

#include "parameters.h"
#include "errors.h"
#include "math.h"
#include "common.h"

#include <mpfr.h>
#include <gmp.h>
#include <mpi.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void linear_distribution_slice_recv_common(
  Linear_Distribution_Slice * const slice,
  const int rank)
{
  MPI_Status status;

  if (MPI_SUCCESS != MPI_Recv(
      &(slice->min_log_alpha),
      1, /* count */
      MPI_INT,
      rank,
      MPI_TAG_SLICE_ALPHAS,
      MPI_COMM_WORLD,
      &status))
  {
    critical("linear_distribution_slice_recv_common(): "
      "Failed to receive alpha.");
  }

  if (MPI_SUCCESS != MPI_Recv(
      &(slice->flags),
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_TAG_SLICE_FLAGS,
      MPI_COMM_WORLD,
      &status))
  {
    critical("linear_distribution_slice_recv_common(): "
      "Failed to receive flags.");
  }

  if (MPI_SUCCESS != MPI_Recv(
      slice->norm_vector,
      slice->dimension, /* count */
      MPI_LONG_DOUBLE,
      rank,
      MPI_TAG_SLICE_NORM_MATRIX,
      MPI_COMM_WORLD,
      &status))
  {
    critical("linear_distribution_slice_recv_common(): "
      "Failed to receive the norm vector.");
  }

  if (MPI_SUCCESS != MPI_Recv(
      &(slice->total_error),
      1, /* count */
      MPI_LONG_DOUBLE,
      rank,
      MPI_TAG_SLICE_TOTAL_ERROR,
      MPI_COMM_WORLD,
      &status))
  {
    critical("linear_distribution_slice_recv_common(): "
      "Failed to receive the total error.");
  }

  /* Compute the total probability. */
  slice->total_probability = 0;
  for (uint32_t i = 0; i < slice->dimension; i++) {
    slice->total_probability += slice->norm_vector[i];
  }
}

void linear_distribution_slice_recv(
  Linear_Distribution_Slice * const slice,
  const int rank)
{
  uint32_t dimension;

  MPI_Status status;

  if (MPI_SUCCESS != MPI_Recv(
      &dimension,
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_TAG_SLICE_DIMENSION,
      MPI_COMM_WORLD,
      &status))
  {
    critical("linear_distribution_slice_recv(): "
      "Failed to receive the dimension.");
  }

  if (dimension != slice->dimension) {
    linear_distribution_slice_clear(slice);
    linear_distribution_slice_init(slice, dimension);
  }

  linear_distribution_slice_recv_common(slice, rank);
}

void linear_distribution_slice_init_recv(
  Linear_Distribution_Slice * const slice,
  const int rank)
{
  uint32_t dimension;

  MPI_Status status;

  if (MPI_SUCCESS != MPI_Recv(
      &dimension,
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_TAG_SLICE_DIMENSION,
      MPI_COMM_WORLD,
      &status))
  {
    critical("linear_distribution_slice_init_recv(): "
      "Failed to receive the dimension.");
  }

  linear_distribution_slice_init(slice, dimension);

  linear_distribution_slice_recv_common(slice, rank);
}

void linear_distribution_slice_send(
  const Linear_Distribution_Slice * const slice,
  const int rank)
{
  if (MPI_SUCCESS != MPI_Send(
      &(slice->dimension),
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_TAG_SLICE_DIMENSION,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_send(): Failed to send the dimension.");
  }

  if (MPI_SUCCESS != MPI_Send(
      &(slice->min_log_alpha),
      1, /* count */
      MPI_INT,
      rank,
      MPI_TAG_SLICE_ALPHAS,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_send(): "
      "Failed to send alpha_d and alpha_r.");
  }

  if (MPI_SUCCESS != MPI_Send(
      &(slice->flags),
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_TAG_SLICE_FLAGS,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_send(): Failed to send flags.");
  }

  if (MPI_SUCCESS != MPI_Send(
      slice->norm_vector,
      slice->dimension, /* count */
      MPI_LONG_DOUBLE,
      rank,
      MPI_TAG_SLICE_NORM_MATRIX,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_send(): "
      "Failed to send the norm vector.");
  }

  if (MPI_SUCCESS != MPI_Send(
      &(slice->total_error),
      1, /* count */
      MPI_LONG_DOUBLE,
      rank,
      MPI_TAG_SLICE_TOTAL_ERROR,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_send(): "
      "Failed to send the total error.");
  }
}

void linear_distribution_slice_init_bcast_recv(
  Linear_Distribution_Slice * const slice,
  const int root)
{
  /* Broadcast dimension. */
  uint32_t dimension;

  if (MPI_SUCCESS != MPI_Bcast(
      &dimension,
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_init_bcast_recv(): "
      "Failed to broadcast the dimension.");
  }

  /* Store dimension. */
  slice->dimension = dimension;

  /* Allocate memory for the norm vector. */
  slice->norm_vector =
    (long double *)malloc(dimension * sizeof(long double));
  if (NULL == slice->norm_vector) {
    critical("linear_distribution_slice_init_bcast_recv(): "
      "Failed to allocate memory.");
  }

  /* Broadcast alpha_d and alpha_r. */
  int32_t min_log_alpha;

  if (MPI_SUCCESS != MPI_Bcast(
      &min_log_alpha,
      1, /* count */
      MPI_INT,
      root,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_init_bcast_recv(): "
      "Failed to broadcast alpha.");
  }

  /* Store alpha. */
  slice->min_log_alpha = min_log_alpha;

  /* Broadcast flags. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->flags),
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_init_bcast_recv(): "
      "Failed to broadcast flags.");
  }

  /* Broadcast norm vector. */
  if (MPI_SUCCESS != MPI_Bcast(
      slice->norm_vector,
      slice->dimension, /* count */
      MPI_LONG_DOUBLE,
      root,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_init_bcast_recv(): "
      "Failed to broadcast the norm vector.");
  }

  /* Broadcast total probability. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->total_probability),
      1, /* count */
      MPI_LONG_DOUBLE,
      root,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_init_bcast_recv(): "
      "Failed to broadcast the total probability.");
  }

  /* Broadcast total error. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->total_error),
      1, /* count */
      MPI_LONG_DOUBLE,
      root,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_init_bcast_recv(): "
      "Failed to broadcast the total error.");
  }
}

void linear_distribution_slice_bcast_send(
  Linear_Distribution_Slice * const slice,
  const int root)
{
  /* Broadcast dimension. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->dimension),
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_bcast_send(): "
      "Failed to broadcast the dimension.");
  }

  /* Broadcast alpha_d and alpha_r. */
  int32_t min_log_alpha = slice->min_log_alpha;

  if (MPI_SUCCESS != MPI_Bcast(
      &min_log_alpha,
      1, /* count */
      MPI_INT,
      root,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_bcast_send(): "
      "Failed to broadcast alpha.");
  }

  /* Broadcast flags. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->flags),
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_bcast_send(): "
      "Failed to broadcast flags.");
  }

  /* Broadcast norm vector. */
  if (MPI_SUCCESS != MPI_Bcast(
      slice->norm_vector,
      slice->dimension, /* count */
      MPI_LONG_DOUBLE,
      root,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_bcast_send(): "
      "Failed to broadcast the norm vector.");
  }

  /* Broadcast total probability. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->total_probability),
      1, /* count */
      MPI_LONG_DOUBLE,
      root,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_init_bcast_recv(): "
      "Failed to broadcast the total probability.");
  }

  /* Broadcast total error. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->total_error),
      1, /* count */
      MPI_LONG_DOUBLE,
      root,
      MPI_COMM_WORLD))
  {
    critical("linear_distribution_slice_bcast_send(): "
      "Failed to broadcast the total error.");
  }
}
