/*!
 * \file    executables_solve_distribution.cpp
 * \ingroup solve_executable
 *
 * \brief   The definition of functions shared by the executables for
 *          post-processing simulated output.
 */

#include "executables_solve_distribution.h"

#include "parameters.h"
#include "timer.h"

#include <stdint.h>
#include <stdio.h>

void solution_status_init(
  Solution_Status * const status,
  const Parameters * const parameters,
  const uint32_t n)
{
  status->m = parameters->m;
  status->s = parameters->s;
  status->l = parameters->l;

  status->n = n;

  timer_statistics_init(&(status->statistics_prepare_system));
  timer_statistics_init(&(status->statistics_solve_system));

  solution_status_reset(status);
}

void solution_status_clear(
  Solution_Status * const status)
{
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
