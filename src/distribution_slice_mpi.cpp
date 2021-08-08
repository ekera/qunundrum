/*!
 * \file    distribution_slice_mpi.cpp
 * \ingroup two_dimensional_distribution_slice
 *
 * \brief   The definition of functions for sending, receiving and broadcasting
 *          slices in two-dimensional probability distributions.
 */

#include "distribution_slice.h"

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

static void distribution_slice_recv_common(
  Distribution_Slice * const slice,
  const int rank)
{
  MPI_Status status;

  /* Receive min_log_alpha_d and min_log_alpha_r. */
  int32_t min_log_alpha[2];

  if (MPI_SUCCESS != MPI_Recv(
      min_log_alpha,
      2, /* count */
      MPI_INT,
      rank,
      MPI_TAG_SLICE_ALPHAS,
      MPI_COMM_WORLD,
      &status))
  {
    critical("distribution_slice_recv_common(): "
      "Failed to receive min_log_alpha_d and min_log_alpha_r.");
  }

  /* Store min_log_alpha_d and min_log_alpha_r. */
  slice->min_log_alpha_d = min_log_alpha[0];
  slice->min_log_alpha_r = min_log_alpha[1];

  /* Receive the flags. */
  if (MPI_SUCCESS != MPI_Recv(
      &(slice->flags),
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_TAG_SLICE_FLAGS,
      MPI_COMM_WORLD,
      &status))
  {
    critical("distribution_slice_recv_common(): "
      "Failed to receive the flags.");
  }

  /* Receive the norm matrix. */
  if (MPI_SUCCESS != MPI_Recv(
      slice->norm_matrix,
      slice->dimension * slice->dimension, /* count */
      MPI_LONG_DOUBLE,
      rank,
      MPI_TAG_SLICE_NORM_MATRIX,
      MPI_COMM_WORLD,
      &status))
  {
    critical("distribution_slice_recv_common(): "
      "Failed to receive the norm matrix.");
  }

  /* Receive the total error. */
  if (MPI_SUCCESS != MPI_Recv(
      &(slice->total_error),
      1, /* count */
      MPI_LONG_DOUBLE,
      rank,
      MPI_TAG_SLICE_TOTAL_ERROR,
      MPI_COMM_WORLD,
      &status))
  {
    critical("distribution_slice_recv(): "
      "Failed to receive the total error.");
  }

  /* Compute the total probability. */
  slice->total_probability = 0;
  for (uint32_t i = 0; i < (slice->dimension * slice->dimension); i++) {
    slice->total_probability += slice->norm_matrix[i];
  }
}

void distribution_slice_recv(
  Distribution_Slice * const slice,
  const int rank)
{
  MPI_Status status;

  /* Receive the dimension. */
  uint32_t dimension;

  if (MPI_SUCCESS != MPI_Recv(
      &dimension,
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_TAG_SLICE_DIMENSION,
      MPI_COMM_WORLD,
      &status))
  {
    critical("distribution_slice_recv(): "
      "Failed to receive the dimension.");
  }

  /* Re-initialize the slice if the dimension differs. */
  if (dimension != slice->dimension) {
    distribution_slice_clear(slice);
    distribution_slice_init(slice, dimension);
  }

  /* Call distribution_slice_recv_common(). */
  distribution_slice_recv_common(slice, rank);
}

void distribution_slice_init_recv(
  Distribution_Slice * const slice,
  const int rank)
{
  MPI_Status status;

  /* Receive the dimension. */
  uint32_t dimension;

  if (MPI_SUCCESS != MPI_Recv(
      &dimension,
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_TAG_SLICE_DIMENSION,
      MPI_COMM_WORLD,
      &status))
  {
    critical("distribution_slice_init_recv(): "
      "Failed to receive the dimension.");
  }

  /* Initialize the slice. */
  distribution_slice_init(slice, dimension);

  /* Call distribution_slice_recv_common(). */
  distribution_slice_recv_common(slice, rank);
}

void distribution_slice_send(
  const Distribution_Slice * const slice,
  const int rank)
{
  /* Send the dimension. */
  if (MPI_SUCCESS != MPI_Send(
      &(slice->dimension),
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_TAG_SLICE_DIMENSION,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_send(): Failed to send the dimension.");
  }

  /* Send min_log_alpha_d and min_log_alpha_r. */
  int32_t min_log_alpha[2] = {
    slice->min_log_alpha_d, slice->min_log_alpha_r
  };

  if (MPI_SUCCESS != MPI_Send(
      min_log_alpha,
      2, /* count */
      MPI_INT,
      rank,
      MPI_TAG_SLICE_ALPHAS,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_send(): "
      "Failed to send min_log_alpha_d and min_log_alpha_r.");
  }

  /* Send the flags. */
  if (MPI_SUCCESS != MPI_Send(
      &(slice->flags),
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_TAG_SLICE_FLAGS,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_send(): Failed to send the flags.");
  }

  /* Send the norm matrix. */
  if (MPI_SUCCESS != MPI_Send(
      slice->norm_matrix,
      slice->dimension * slice->dimension, /* count */
      MPI_LONG_DOUBLE,
      rank,
      MPI_TAG_SLICE_NORM_MATRIX,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_send(): Failed to send the norm matrix.");
  }

  /* Send the total error. */
  if (MPI_SUCCESS != MPI_Send(
      &(slice->total_error),
      1, /* count */
      MPI_LONG_DOUBLE,
      rank,
      MPI_TAG_SLICE_TOTAL_ERROR,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_send(): Failed to send the total error.");
  }
}

void distribution_slice_init_bcast_recv(
  Distribution_Slice * const slice,
  const int root)
{
  /* Broadcast the dimension. */
  uint32_t dimension;

  if (MPI_SUCCESS != MPI_Bcast(
      &dimension,
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_init_bcast_recv(): "
      "Failed to broadcast the dimension.");
  }

  /* Store the dimension. */
  slice->dimension = dimension;

  /* Allocate memory for the norm matrix. */
  slice->norm_matrix =
    (long double *)malloc(dimension * dimension * sizeof(long double));
  if (NULL == slice->norm_matrix) {
    critical("distribution_slice_init_bcast_recv(): "
      "Failed to allocate memory.");
  }

  /* Broadcast min_log_alpha_d and min_log_alpha_r. */
  int32_t min_log_alpha[2];

  if (MPI_SUCCESS != MPI_Bcast(
      min_log_alpha,
      2, /* count */
      MPI_INT,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_init_bcast_recv(): "
      "Failed to broadcast min_log_alpha_d and min_log_alpha_r.");
  }

  /* Store min_log_alpha_d and min_log_alpha_r. */
  slice->min_log_alpha_d = min_log_alpha[0];
  slice->min_log_alpha_r = min_log_alpha[1];

  /* Broadcast the flags. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->flags),
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_init_bcast_recv(): "
      "Failed to broadcast the flags.");
  }

  /* Broadcast the norm matrix. */
  if (MPI_SUCCESS != MPI_Bcast(
      slice->norm_matrix,
      slice->dimension * slice->dimension, /* count */
      MPI_LONG_DOUBLE,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_init_bcast_recv(): "
      "Failed to broadcast the norm matrix.");
  }

  /* Broadcast the total probability. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->total_probability),
      1, /* count */
      MPI_LONG_DOUBLE,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_init_bcast_recv(): "
      "Failed to broadcast the total probability.");
  }

  /* Broadcast the total error. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->total_error),
      1, /* count */
      MPI_LONG_DOUBLE,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_init_bcast_recv(): "
      "Failed to broadcast the total error.");
  }
}

void distribution_slice_bcast_send(
  Distribution_Slice * const slice,
  const int root)
{
  /* Broadcast the dimension. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->dimension),
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_bcast_send(): "
      "Failed to broadcast the dimension.");
  }

  /* Broadcast min_log_alpha_d and min_log_alpha_r. */
  int32_t min_log_alpha[2] = {
    slice->min_log_alpha_d, slice->min_log_alpha_r
  };

  if (MPI_SUCCESS != MPI_Bcast(
      min_log_alpha,
      2, /* count */
      MPI_INT,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_bcast_send(): "
      "Failed to broadcast min_log_alpha_d and min_log_alpha_r.");
  }

  /* Broadcast the flags. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->flags),
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_bcast_send(): Failed to broadcast the flags.");
  }

  /* Broadcast the norm matrix. */
  if (MPI_SUCCESS != MPI_Bcast(
      slice->norm_matrix,
      slice->dimension * slice->dimension, /* count */
      MPI_LONG_DOUBLE,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_bcast_send(): "
      "Failed to broadcast the norm matrix.");
  }

  /* Broadcast the total probability. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->total_probability),
      1, /* count */
      MPI_LONG_DOUBLE,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_bcast_send(): "
      "Failed to broadcast the total probability.");
  }

  /* Broadcast the total error. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->total_error),
      1, /* count */
      MPI_LONG_DOUBLE,
      root,
      MPI_COMM_WORLD))
  {
    critical("distribution_slice_bcast_send(): "
      "Failed to broadcast the total error.");
  }
}
