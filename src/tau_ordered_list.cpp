/*!
 * \file    tau_ordered_list.cpp
 * \ingroup estimating_volume_quotients
 *
 * \brief   The definition of functions for manipulating truncated order lists
 *          of tau estimates.
 */

#include "tau_ordered_list.h"

#include "common.h"
#include "errors.h"

#include <mpi.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*!
 * \brief   An MPI tag used to send the item count for an ordered list.
 */
#define MPI_TAG_TAU_ORDERED_LIST_COUNT            17051

/*!
 * \brief   An MPI tag used to send the tau values in an ordered list.
 */
#define MPI_TAG_TAU_ORDERED_LIST_TAU              17052

void tau_ordered_list_init(
  Tau_Ordered_List * const list,
  const unsigned int capacity)
{
  list->capacity = capacity;
  list->size = 0;

  list->tau = (long double *)malloc(capacity * sizeof(long double));
  if (NULL == (list->tau)) {
    critical("tau_ordered_list_init(): Failed to allocate memory.");
  }

  memset(list->tau, 0, capacity * sizeof(long double));
}

void tau_ordered_list_clear(
    Tau_Ordered_List * const list)
{
  free(list->tau);
  list->tau = NULL;

  memset(list, 0, sizeof(Tau_Ordered_List));
}

void tau_ordered_list_print(
  const Tau_Ordered_List * const list)
{
  printf("List(size = %u, capacity = %u):\n", list->size, list->capacity);
  for (uint32_t i = 0; i < list->size; i++) {
    printf(" %u: %Lf\n", i, list->tau[i]);
  }
  printf("\n");
}

void tau_ordered_list_insert(
  Tau_Ordered_List * const list,
  const long double tau)
{
  /* Check if the list is empty and if so simply insert the value. */
  if (0 == list->size) {
    list->tau[0] = tau;
    list->size = 1;

    return;
  }

  /* Check if the list is at full capacity. */
  if (list->capacity == list->size) {
    /* If the value to be inserted is smaller than all elements in the list
     * then we may simply return at this point. */
    if (tau < list->tau[0]) {
      return;
    }

    /* Compare to the other elements in the list. */
    for (uint32_t i = 1; i < list->size; i++) {
      if (tau < list->tau[i]) {
        /* Make room by shifting the list to the left and inserting. */
        for (unsigned int j = 0; j < i - 1; j++) {
          list->tau[j] = list->tau[j + 1];
        }
        list->tau[i - 1] = tau;
        return;
      }
    }

    /* Shift the whole list one step to the left and insert. */
    for (uint32_t j = 1; j < list->size; j++) {
      list->tau[j - 1] = list->tau[j];
    }
    list->tau[list->size - 1] = tau;
  } else {
    /* Expand the list. */
    for (uint32_t i = 0; i < list->size; i++) {
      if (tau < list->tau[i]) {
        for (uint32_t j = list->size; j > i; j--) {
          list->tau[j] = list->tau[j - 1];
        }

        list->tau[i] = tau;
        list->size += 1;

        return;
      }
    }

    list->tau[list->size] = tau;
    list->size += 1;
  }
}

void tau_ordered_list_send_merge(
  const Tau_Ordered_List * const list,
  const int rank)
{
  if (MPI_SUCCESS != MPI_Send(
    &(list->size),
    1,
    MPI_UNSIGNED,
    rank,
    MPI_TAG_TAU_ORDERED_LIST_COUNT,
    MPI_COMM_WORLD))
  {
    critical("tau_ordered_list_send_merge(): Failed to send the list count.");
  }

  if (0 == list->size) {
    return;
  }

  if (MPI_SUCCESS != MPI_Send(
    list->tau,
    list->size,
    MPI_LONG_DOUBLE,
    rank,
    MPI_TAG_TAU_ORDERED_LIST_TAU,
    MPI_COMM_WORLD))
  {
    critical("tau_ordered_list_send_merge(): "
      "Failed to send the list of tau values.");
  }
}

void tau_ordered_list_recv_merge(
  Tau_Ordered_List * const list,
  const int rank)
{
  unsigned int count;

  MPI_Status status;

  if (MPI_SUCCESS != MPI_Recv(
    &count,
    1,
    MPI_UNSIGNED,
    rank,
    MPI_TAG_TAU_ORDERED_LIST_COUNT,
    MPI_COMM_WORLD,
    &status))
  {
    critical("tau_ordered_list_recv_merge(): "
      "Failed to receive the list count.");
  }

  if (0 == count) {
    return;
  }

  /* Allocate space for the elements to be inserted. */
  long double * tau = (long double *)malloc(count * sizeof(long double));
  if (NULL == tau) {
    critical("tau_ordered_list_recv_merge(): "
      "Failed to allocate memory.");
  }

  if (MPI_SUCCESS != MPI_Recv(
    tau,
    count,
    MPI_LONG_DOUBLE,
    rank,
    MPI_TAG_TAU_ORDERED_LIST_TAU,
    MPI_COMM_WORLD,
    &status))
  {
    critical("tau_ordered_list_recv_merge(): "
      "Failed to receive the list of tau values.");
  }

  for (uint32_t i = 0; i < count; i++) {
    tau_ordered_list_insert(list, tau[i]);
  }

  /* Clear memory. */
  free(tau);
  tau = NULL;
}
