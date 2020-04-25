/*!
 * \file    diagonal_distribution_slice_mpi.cpp
 * \ingroup diagonal_distribution_slice
 *
 * \brief   The definition of functions for sending, receiving and broadcasting
 *          slices in diagonal probability distributions.
 */

#include "diagonal_distribution_slice.h"

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

static void diagonal_distribution_slice_recv_common(
  Diagonal_Distribution_Slice * const slice,
  const int rank)
{
  MPI_Status status;

  if (MPI_SUCCESS != MPI_Recv(
      &(slice->min_log_alpha_r),
      1, /* count */
      MPI_INT,
      rank,
      MPI_TAG_SLICE_MIN_LOG_ALPHA_R,
      MPI_COMM_WORLD,
      &status))
  {
    critical("diagonal_distribution_slice_recv_common(): "
      "Failed to receive min_log_alpha_r.");
  }

  if (MPI_SUCCESS != MPI_Recv(
      &(slice->offset_alpha_d),
      1, /* count */
      MPI_INT,
      rank,
      MPI_TAG_SLICE_OFFSET_ALPHA_D,
      MPI_COMM_WORLD,
      &status))
  {
    critical("diagonal_distribution_slice_recv_common(): "
      "Failed to receive offset_alpha_d.");
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
    critical("diagonal_distribution_slice_recv_common(): "
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
    critical("diagonal_distribution_slice_recv_common(): "
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
    critical("diagonal_distribution_slice_recv_common(): "
      "Failed to receive the total error.");
  }

  /* Compute the total probability. */
  slice->total_probability = 0;
  for (uint32_t i = 0; i < slice->dimension; i++) {
    slice->total_probability += slice->norm_vector[i];
  }
}

void diagonal_distribution_slice_recv(
  Diagonal_Distribution_Slice * const slice,
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
    critical("diagonal_distribution_slice_recv(): "
      "Failed to receive the dimension.");
  }

  if (dimension != slice->dimension) {
    diagonal_distribution_slice_clear(slice);
    diagonal_distribution_slice_init(slice, dimension);
  }

  diagonal_distribution_slice_recv_common(slice, rank);
}

void diagonal_distribution_slice_init_recv(
  Diagonal_Distribution_Slice * const slice,
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
    critical("diagonal_distribution_slice_init_recv(): "
      "Failed to receive the dimension.");
  }

  diagonal_distribution_slice_init(slice, dimension);

  diagonal_distribution_slice_recv_common(slice, rank);
}

void diagonal_distribution_slice_send(
  const Diagonal_Distribution_Slice * const slice,
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
    critical("diagonal_distribution_slice_send(): "
      "Failed to send the dimension.");
  }

  if (MPI_SUCCESS != MPI_Send(
      &(slice->min_log_alpha_r),
      1, /* count */
      MPI_INT,
      rank,
      MPI_TAG_SLICE_MIN_LOG_ALPHA_R,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_slice_send(): "
      "Failed to send min_log_alpha_r.");
  }

  if (MPI_SUCCESS != MPI_Send(
      &(slice->offset_alpha_d),
      1, /* count */
      MPI_INT,
      rank,
      MPI_TAG_SLICE_OFFSET_ALPHA_D,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_slice_send(): "
      "Failed to send offset_alpha_d.");
  }

  if (MPI_SUCCESS != MPI_Send(
      &(slice->flags),
      1, /* count */
      MPI_UNSIGNED,
      rank,
      MPI_TAG_SLICE_FLAGS,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_slice_send(): "
      "Failed to send flags.");
  }

  if (MPI_SUCCESS != MPI_Send(
      slice->norm_vector,
      slice->dimension, /* count */
      MPI_LONG_DOUBLE,
      rank,
      MPI_TAG_SLICE_NORM_MATRIX,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_slice_send(): "
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
    critical("diagonal_distribution_slice_send(): "
      "Failed to send the total error.");
  }
}

void diagonal_distribution_slice_init_bcast_recv(
  Diagonal_Distribution_Slice * const slice,
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
    critical("diagonal_distribution_slice_init_bcast_recv(): "
      "Failed to broadcast the dimension.");
  }

  /* Store dimension. */
  slice->dimension = dimension;

  /* Allocate memory for the norm vector. */
  slice->norm_vector =
    (long double *)malloc(dimension * sizeof(long double));
  if (NULL == slice->norm_vector) {
    critical("diagonal_distribution_slice_init_bcast_recv(): "
      "Failed to allocate memory.");
  }

  /* Broadcast min_log_alpha_r. */
  int32_t min_log_alpha_r;

  if (MPI_SUCCESS != MPI_Bcast(
      &min_log_alpha_r,
      1, /* count */
      MPI_INT,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_slice_init_bcast_recv(): "
      "Failed to broadcast min_log_alpha_r.");
  }

  /* Store alpha. */
  slice->min_log_alpha_r = min_log_alpha_r;

  /* Broadcast offset_alpha_d. */
  int32_t offset_alpha_d;

  if (MPI_SUCCESS != MPI_Bcast(
      &offset_alpha_d,
      1, /* count */
      MPI_INT,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_slice_init_bcast_recv(): "
      "Failed to broadcast offset_alpha_d.");
  }

  /* Store alpha. */
  slice->offset_alpha_d = offset_alpha_d;

  /* Broadcast flags. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->flags),
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_slice_init_bcast_recv(): "
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
    critical("diagonal_distribution_slice_init_bcast_recv(): "
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
    critical("diagonal_distribution_slice_init_bcast_recv(): "
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
    critical("diagonal_distribution_slice_init_bcast_recv(): "
      "Failed to broadcast the total error.");
  }
}

void diagonal_distribution_slice_bcast_send(
  Diagonal_Distribution_Slice * const slice,
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
    critical("diagonal_distribution_slice_bcast_send(): "
      "Failed to broadcast the dimension.");
  }

  /* Broadcast min_log_alpha_r. */
  int32_t min_log_alpha_r = slice->min_log_alpha_r;

  if (MPI_SUCCESS != MPI_Bcast(
      &min_log_alpha_r,
      1, /* count */
      MPI_INT,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_slice_bcast_send(): "
      "Failed to broadcast min_log_alpha_r.");
  }

  /* Broadcast offset_alpha_d. */
  int32_t offset_alpha_d = slice->offset_alpha_d;

  if (MPI_SUCCESS != MPI_Bcast(
      &offset_alpha_d,
      1, /* count */
      MPI_INT,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_slice_bcast_send(): "
      "Failed to broadcast offset_alpha_d.");
  }

  /* Broadcast flags. */
  if (MPI_SUCCESS != MPI_Bcast(
      &(slice->flags),
      1, /* count */
      MPI_UNSIGNED,
      root,
      MPI_COMM_WORLD))
  {
    critical("diagonal_distribution_slice_bcast_send(): "
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
    critical("diagonal_distribution_slice_bcast_send(): "
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
    critical("diagonal_distribution_slice_bcast_send(): "
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
    critical("diagonal_distribution_slice_bcast_send(): "
      "Failed to broadcast total error.");
  }
}
