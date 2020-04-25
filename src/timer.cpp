/*!
 * \file    timer.cpp
 * \ingroup timer
 *
 * \brief   The definition of functions for manipulating timers and collecting
 *          timing statistics.
 */

#include "timer.h"

#include "errors.h"

#include <time.h>

#include <stdint.h>
#include <string.h>

void timer_start(
  Timer * const timer)
{
  clock_gettime(CLOCK_MONOTONIC_RAW, &(timer->start));
  memcpy(&(timer->stop), &(timer->start), sizeof(struct timespec));
}

uint64_t timer_stop(
  Timer * const timer)
{
  clock_gettime(CLOCK_MONOTONIC_RAW, &(timer->stop));

  int64_t delta;

  delta  = (int64_t)(timer->stop.tv_sec) -
              (int64_t)(timer->start.tv_sec);

  delta *= 1000 * 1000; /* compensate for us resolution */

  delta += ((int64_t)(timer->stop.tv_nsec) -
              (int64_t)(timer->start.tv_nsec)) / 1000; /* reduce ns to us */

  if (delta < 0) {
    delta = 0;
  }

  return (uint64_t)delta;
}

void timer_statistics_init(
  Timer_Statistics * const statistics)
{
  timer_statistics_reset(statistics);
}

void timer_statistics_reset(
  Timer_Statistics * const statistics)
{
  statistics->min = (uint64_t)(-1);
  statistics->max = 0;
  statistics->sum = 0;
  statistics->count = 0;
}

void timer_statistics_insert(
  Timer_Statistics * const statistics,
  const uint64_t t)
{
  if (statistics->min > t) {
    statistics->min = t;
  }

  if (statistics->max < t) {
    statistics->max = t;
  }

  statistics->sum += t;
  statistics->count += 1;
}

void timer_statistics_export(
  uint64_t * const min,
  uint64_t * const max,
  uint64_t * const sum,
  uint32_t * const count,
  double * const avg,
  const Timer_Statistics * const statistics)
{
  if (NULL != min) {
    if (statistics->count == 0) {
      (*min) = 0;
    } else {
      (*min) = statistics->min;
    }
  }

  if (NULL != max) {
    if (statistics->count == 0) {
      (*max) = 0;
    } else {
      (*max) = statistics->max;
    }
  }

  if (NULL != sum) {
    if (statistics->count == 0) {
      (*sum) = 0;
    } else {
      (*sum) = statistics->sum;
    }
  }

  if (NULL != count) {
    (*count) = statistics->count;
  }

  if (NULL != avg) {
    if (statistics->count == 0) {
      (*avg) = 0;
    } else {
      (*avg) = ((double)statistics->sum) / ((double)statistics->count);
    }
  }
}
