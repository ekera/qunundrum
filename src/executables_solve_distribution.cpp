/*!
 * \file    executables_solve_distribution.cpp
 * \ingroup solve_executable
 *
 * \brief   The definition of functions shared by the executables for
 *          post-processing simulated output.
 */

#include "executables_solve_distribution.h"

#include "timer.h"
#include "common.h"

#include <stdint.h>
#include <stdio.h>

void solution_status_init(
  Solution_Status * const status,
  const uint32_t m,
  const uint32_t sigma,
  const uint32_t s,
  const uint32_t l,
  const uint32_t n,
  const bool has_sigma)
{
  /* Store m, sigma, s and l. */
  status->m = m;
  status->sigma = sigma;
  status->s = s;
  status->l = l;

  /* Store n. */
  status->n = n;

  /* Store information on sigma. */
  status->has_sigma = has_sigma;

  /* Initialize the timers.  */
  timer_statistics_init(&(status->statistics_prepare_system));
  timer_statistics_init(&(status->statistics_solve_system));

  /* Reset the solution status. */
  solution_status_reset(status);
}

void solution_status_clear(
  Solution_Status * const status)
{
  /* This function currently performs no operation. */
  (void)status;

  return;
}

void solution_status_reset(
  Solution_Status * const status)
{
  status->success_count = 0;
  status->fail_count = 0;
  status->fail_out_of_bounds_count = 0;
  status->issued_count = 0;

  timer_statistics_reset(&(status->statistics_solve_system));
  timer_statistics_reset(&(status->statistics_prepare_system));
}

void solution_status_print(
  FILE * const file,
  const Solution_Status * const status,
  const char * const suffix)
{
  /* Extract statistics. */
  double avg_time_prepare_system;

  timer_statistics_export(
    NULL, /* min */
    NULL, /* max */
    NULL, /* sum */
    NULL, /* count */
    &avg_time_prepare_system,
    &(status->statistics_prepare_system));

  uint64_t min_time_solve_system;
  uint64_t max_time_solve_system;
  double avg_time_solve_system;

  timer_statistics_export(
    &min_time_solve_system, /* min */
    &max_time_solve_system, /* max */
    NULL, /* sum */
    NULL, /* count */
    &avg_time_solve_system,
    &(status->statistics_solve_system));

  if (TRUE == (status->has_sigma)) {
    fprintf(file, "m: %u sigma: %u %c: %u n: %u -- success: %u -- fail: "
      "%u (%u) -- prepare: %9.3f ms solve: %9.3f ms [%9.3f, %9.3f] %s\n",
      status->m,
      status->sigma,
      (0 != status->s) ? 's' : 'l',
      (0 != status->s) ? status->s : status->l,
      status->n,
      status->success_count,
      status->fail_count,
      status->fail_out_of_bounds_count,
      avg_time_prepare_system / 1000, /* reduce us to ms */
      avg_time_solve_system / 1000, /* reduce us to ms */
      ((double)min_time_solve_system) / 1000, /* reduce us to ms */
      ((double)max_time_solve_system) / 1000, /* reduce us to ms */
      (NULL == suffix) ? "" : suffix);
  } else {
    fprintf(file, "m: %u %c: %u n: %u -- success: %u -- fail: %u (%u) "
      "-- prepare: %9.3f ms solve: %9.3f ms [%9.3f, %9.3f] %s\n",
      status->m,
      (0 != status->s) ? 's' : 'l',
      (0 != status->s) ? status->s : status->l,
      status->n,
      status->success_count,
      status->fail_count,
      status->fail_out_of_bounds_count,
      avg_time_prepare_system / 1000, /* reduce us to ms */
      avg_time_solve_system / 1000, /* reduce us to ms */
      ((double)min_time_solve_system) / 1000, /* reduce us to ms */
      ((double)max_time_solve_system) / 1000, /* reduce us to ms */
      (NULL == suffix) ? "" : suffix);
  }
}
