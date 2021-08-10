/*!
 * \file    executables_solve_distribution.h
 * \ingroup solve_executable
 *
 * \brief   The declaration of data structures and functions shared by the
 *          executables for post-processing simulated output.
 */

/*!
 * \defgroup solve_executable Executables for post-processing simulated output
 * \ingroup  executable
 *
 * \brief    A module for executables for post-processing simulated output.
 */

#ifndef EXECUTABLES_SOLVE_DISTRIBUTION_H
#define EXECUTABLES_SOLVE_DISTRIBUTION_H

#include "timer.h"
#include "common.h"

#include <stdint.h>
#include <stdio.h>

/*!
 * \brief   An enumeration of states that an adaptive search may assume.
 *
 * For the initial n, the pivot state #SEARCH_STATE_PIVOT is assumed.
 *
 * If the number of failures is at or below the limit for the initial n, the
 * state will change to #SEARCH_STATE_INCREASE and n will be decreased in
 * decrements of one until the limit is exceeded.
 *
 * Analogously, if the number of failures is above the limit for the initial n,
 * the state will change to #SEARCH_STATE_DECREASE and n will be incremented in
 * increments of one until the limit is met or exceeded.
 */
typedef enum {

  /*!
   * \brief   Indicates that the search is in the initial pivot state.
   */
  SEARCH_STATE_PIVOT = 1,

  /*!
   * \brief   Indicates that the search is in the state where n is increased.
   */
  SEARCH_STATE_INCREASE = 2,

  /*!
   * \brief   Indicates that the search is in the state where n is decreased.
   */
  SEARCH_STATE_DECREASE = 0

} Adaptive_Search_State;

/*!
 * \brief   An enumeration of search strategies.
 */
typedef enum {

  /*!
   * \brief   Indicates that the number of pairs n included in the lattice
   *          should be adaptively increased on decreased to find the minimum
   *          number of pairs required to bound the failure rate below a limit.
   */
  SEARCH_STRATEGY_ADAPTIVE = 1,

  /*!
   * \brief   Indicates that the number of pairs n included in the lattice
   *          should be kept fixed to the initial value given on the command
   *          line, and that all sample sets should be solved.
   */
  SEARCH_STRATEGY_NON_ADAPTIVE = 2,

  /*!
   * \brief   Indicates that the number of pairs n included in the lattice
   *          should be kept fixed to the initial value given on the command
   *          line, and that the executable should abort early if the limit on
   *          the failure rate is exceeded.
   */
  SEARCH_STRATEGY_NON_ADAPTIVE_EARLY_ABORT = 3,

  /*!
   * \brief   Indicates that no search strategy has been specified.
   */
  SEARCH_STRATEGY_DEFAULT = 0

} Search_Strategy;

/*!
 * \brief   An enumeration of solution methods.
 */
typedef enum {
  /*!
   * \brief   Indicates that the vector v should be mapped to the closest
   *          lattice vector to find the vector u sought.
   */
  SOLUTION_METHOD_CLOSEST = 1,

  /*!
   * \brief   Indicates that all lattice vectors in a ball around a vector v
   *          should be enumerated to find the vector u sought.
   */
  SOLUTION_METHOD_ENUMERATE = 2,

  /*!
   * \brief   Indicates that no solution method has been specified.
   */
  SOLUTION_METHOD_DEFAULT = 0
} Solution_Method;

/*!
 * \brief   An enumeration of options for detecting if the order is smooth.
 */
typedef enum {
  /*!
   * \brief   Indicates that smooth orders should be detected and leveraged to 
   *          faciliate solving when possible.
   */
  DETECT_SMOOTH_ORDER_TRUE = 1,

  /*!
   * \brief   Indicates that smooth orders should not be detected, and hence 
   *          not leveraged to faciliate solving.
   */
  DETECT_SMOOTH_ORDER_FALSE = 2,

  /*!
   * \brief   Indicates that no option for detecting smooth orders has been 
   *          specified.
   */
  DETECT_SMOOTH_ORDER_DEFAULT = 0
} Detect_Smooth_Order_Option;

/*!
 * \brief   A data structure for holding the current solution status and for
 *          printing information and statistics to log files and to the console.
 */
typedef struct {
  /*!
   * \brief   The length m in bits of the logarithm d or order r.
   */
  uint32_t m;

  /*!
   * \brief   The parameter sigma.
   */
  uint32_t sigma;

  /*!
   * \brief   The tradeoff factor s.
   *
   * Set to zero if l is explicitly specified.
   */
  uint32_t s;

  /*!
   * \brief   The parameter l.
   */
  uint32_t l;

  /*!
   * \brief   The number of runs n for which this data structure is valid.
   *
   * Note that n as stored in this structure may differ from n passed to this
   * executable via the command line if an adaptive search strategy is used.
   */
  uint32_t n;

  /*!
   * \brief   A flag that is set to #TRUE if sigma is to be printed, and to 
   *          #FALSE otherwise.
   */
  bool has_sigma;

  /*!
   * \brief   The number of jobs issued for the current n.
   */
  uint32_t issued_count;

  /*!
   * \brief   The number of successfully completed jobs for the current n.
   */
  uint32_t success_count;

  /*!
   * \brief   The number of failed jobs.
   */
  uint32_t fail_count;

  /*!
   * \brief   The number of jobs that failed due to out of bounds sampling.
   *
   * This count is included in the count given in Solution_Status::fail_count.
   */
  uint32_t fail_out_of_bounds_count;

  /*!
   * \brief   Timing statistics for the time required to prepare the system.
   */
  Timer_Statistics statistics_prepare_system;

  /*!
   * \brief   Timing statistics for the time required to solve the system.
   *
   * This includes the runtime for the lattice basis reduction algorithms, GSO,
   * Babai's nearest plane algorithm and enumeration.
   */
  Timer_Statistics statistics_solve_system;
} Solution_Status;

/*!
 * \name Initialization
 * \{
 */

/*!
 * \brief   Initializes a solution status data structure.
 *
 * \param[in, out] status   The solution status data structure to initialize.
 * \param[in] m             The bit length m of d or r.
 * \param[in] sigma         The parameter sigma.
 * \param[in] s             The tradeoff factor s, or zero if not specified.
 * \param[in] l             The parameter l.
 * \param[in] n             The initial number of runs n for which to solve.
 * \param[in] has_sigma     A flag that should be set to #TRUE if sigma should
 *                          be printed, and to #FALSE otherwise.
 */
void solution_status_init(
  Solution_Status * const status,
  const uint32_t m,
  const uint32_t sigma,
  const uint32_t s,
  const uint32_t l,
  const uint32_t n,
  const bool has_sigma = FALSE);

/*!
 * \brief   Clears a solution status data structure.
 *
 * This function currently performs no operation. It is reserved for future
 * use, and should be called to ensure forward compatibility.
 * 
 * \param[in, out] status   The solution status data structure to clear.
 */
void solution_status_clear(
  Solution_Status * const status);

/*!
 * \}
 */

/*!
 * \name State management
 * \{
 */

/*!
 * \brief   Resets a solution status data structure.
 *
 * \param[in, out] status   The solution status data structure to reset.
 */
void solution_status_reset(
  Solution_Status * const status);

/*!
 * \}
 */

/*!
 * \name Printing
 * \{
 */

/*!
 * \brief   Prints a single line human-readable summary of a solution status
 *          data structure to a file.
 *
 * \param[in, out] file     The file to which to print.
 * \param[in] status        The solution status data structure to print.
 * \param[in] suffix        An optional suffix to print at the end of the line.
 */
void solution_status_print(
  FILE * const file,
  const Solution_Status * const status,
  const char * const suffix = NULL);

/*!
 * \}
 */

#endif /* EXECUTABLES_SOLVE_DISTRIBUTION_H */
