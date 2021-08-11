/*!
 * \file    linear_distribution_loader.h
 * \ingroup linear_distribution_loader
 *
 * \brief   The declaration of functions for loading linear probability
 *          distributions and the definition of data structures for handling
 *          such loading operations.
 */

/*!
 * \defgroup linear_distribution_loader Linear distribution loader
 * \ingroup  linear_distribution
 *
 * \brief    A module for a multi-threaded loader for linear probability
 *           distributions.
 */

#ifndef LINEAR_DISTRIBUTION_LOADER_H
#define LINEAR_DISTRIBUTION_LOADER_H

#include "linear_distribution.h"
#include "thread_pool.h"

#include <stdint.h>

#include <sys/types.h>

/*!
 * \brief   The maximum number of loaded linear distributions to simultaneously 
 *          reside in primary memory at any one time.
 * 
 * This cap primarily serves to prevent excessive use of primary memory.
 * 
 * Note that popped distributions also reside in primary memory until they are 
 * deallocated. They do not count in this context.
 */
#define LINEAR_DISTRIBUTION_LOADER_MAX_SIMULTANEOUSLY_LOADED_DISTRIBUTIONS 8

/*!
 * \brief   A data structure for handling multi-threaded loading of linear
 *          probability distributions.
 *
 * \ingroup linear_distribution_loader
 */
typedef struct {
  /*!
   * \brief   Pointers to the distributions.
   */
  Linear_Distribution ** distributions;

  /*!
   * \brief   The paths to the distributions to be loaded.
   */
  char ** paths;

  /*!
   * \brief   The total number of distributions to be loaded.
   */
  uint32_t count;

  /*!
   * \brief   The offset of the distribution currently being processed within
   *          Linear_Distribution_Loader::distributions.
   */
  uint32_t offset;

  /*!
   * \brief   A mutex used to lock this data structure.
   */
  pthread_mutex_t mutex;

  /*!
   * \brief   A pool of worker threads used to load distributions.
   */
  Thread_Pool pool;
} Linear_Distribution_Loader;

/*!
 * \name Initialization
 * \{
 */

/*!
 * \brief   Initializes the distribution loader with the paths to one or more
 *          distributions to be loaded.
 *
 * \param[in, out] loader   The distribution loader to be initialized.
 * \param[in] paths         The paths to the distributions to be loaded.
 * \param count             The number of distributions to be loaded.
 */
void linear_distribution_loader_init(
  Linear_Distribution_Loader * const loader,
  char ** paths,
  const uint32_t count);

/*!
 * \brief Clears the distribution loader.
 *
 * \param[in, out] loader   An initialized distribution loader to be cleared.
 */
void linear_distribution_loader_clear(
  Linear_Distribution_Loader * const loader);

/*!
 * \}
 */

/*!
 * \name Popping
 * \{
 */

/*!
 * \brief   Pops a distribution from the distribution loader.
 *
 * \param[in, out] loader   An initialized distribution loader.
 *
 * \return  The point to the distribution popped, or NULL, if there are no more
 *          distributions left to pop.
 */
Linear_Distribution * linear_distribution_loader_pop(
  Linear_Distribution_Loader * const loader);

/*!
 * \}
 */

#endif /* LINEAR_DISTRIBUTION_LOADER_H */
