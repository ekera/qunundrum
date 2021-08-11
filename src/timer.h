/*!
 * \file    timer.h
 * \ingroup timer
 *
 * \brief   The declaration of data structures for implementing timers and
 *          timing statistics, and of functions for manipulating such timers and
 *          collecting timing statistics.
 */

/*!
 * \defgroup timer  Timers and timing statistics
 * \ingroup  utility
 *
 * \brief    A module for timers and for collecting timing statistics.
 */

#ifndef TIMER_H
#define TIMER_H

#include <stdint.h>
#include <time.h>

/*!
 * \brief   A data structure for a high resolution timer.
 */
typedef struct {
  /*!
   * \brief   The start time.
   */
  struct timespec start;

  /*!
   * \brief   The stop time.
   */
  struct timespec stop;
} Timer;

/*!
 * \brief   A data structure for holding high resolution timing statistics.
 */
typedef struct {
  /*!
   * \brief   The minimum time sample.
   */
  uint64_t min;

  /*!
   * \brief   The maximum time sample.
   */
  uint64_t max;

  /*!
   * \brief   The sum of all time samples added to this data structure.
   */
  uint64_t sum;

  /*!
   * \brief   The number of time samples added to this data structure.
   */
  uint32_t count;
} Timer_Statistics;

/*!
 * \brief Starts a timer.
 *
 * \param[in, out] timer  The timer to start.
 */
void timer_start(
  Timer * const timer);

/*!
 * \brief Stops a timer.
 *
 * \param[in, out] timer  The timer to stop.
 *
 * \return The number of us elapsed inbetween the point in time when the timer
 *         was started and the time it was stopped.
 */
uint64_t timer_stop(
  Timer * const timer);

/*!
 * \brief Initializes a timing statistics data structure.
 *
 * \param[in, out] statistics The timing statistics data structure to
 *                            initialize.
 */
void timer_statistics_init(
  Timer_Statistics * const statistics);

/*!
 * \brief Resets a timing statistics data structure.
 *
 * \param[in, out] statistics The timing statistics data structure to reset.
 */
void timer_statistics_reset(
  Timer_Statistics * const statistics);

/*!
 * \brief Inserts a time sample into a timing statistics data structure.
 *
 * \param[in, out] statistics The timing statistics data structure.
 * \param[in] t               The time sample.
 */
void timer_statistics_insert(
  Timer_Statistics * const statistics,
  const uint64_t t);

/*!
 * \brief Exports statistics from a timing statistics data structure.
 *
 * \param[out] min        A pointer to an integer to set to the minimum time
 *                        sample, or NULL if not of interest.
 * \param[out] max        A pointer to an integer to set to the maximum time
 *                        sample, or NULL if not of interest.
 * \param[out] sum        A pointer to an integer to set to the sum of all time
 *                        samples, or NULL if not of interest.
 * \param[out] count      A pointer to an integer to set to the number of time
 *                        samples, or NULL if not of interest.
 * \param[out] avg        A pointer to a double to set to the average time
 *                        sample, or NULL if not of interest.
 * \param[in] statistics  The timing statistics data structure.
 */
void timer_statistics_export(
  uint64_t * const min,
  uint64_t * const max,
  uint64_t * const sum,
  uint32_t * const count,
  double * const avg,
  const Timer_Statistics * const statistics);

#endif /* TIMER_H */
