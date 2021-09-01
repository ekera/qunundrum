/*!
 * \file    distribution_loader.h
 * \ingroup two_dimensional_distribution_loader
 *
 * \brief   The declaration of functions for loading two-dimensional probability
 *          distributions, and the definition of data structures for handling
 *          such loading operations.
 */

/*!
 * \defgroup  two_dimensional_distribution_loader \
 *            Two-dimensional distribution loader
 * \ingroup   two_dimensional_distribution
 *
 * \brief   A module for a multi-threaded loader for two-dimensional
 *          probability distributions.
 */

#ifndef DISTRIBUTION_LOADER_H
#define DISTRIBUTION_LOADER_H

#include "distribution.h"
#include "thread_pool.h"

#include <stdint.h>

#include <sys/types.h>

/*!
 * \brief   The maximum number of loaded two-dimensional distributions to
 *          simultaneously reside in primary memory at any one time.
 *
 * This cap primarily serves to prevent excessive use of primary memory.
 *
 * Note that popped distributions also reside in primary memory until they are
 * deallocated. They do not count in this context.
 */
#define DISTRIBUTION_LOADER_MAX_SIMULTANEOUSLY_LOADED_DISTRIBUTIONS 3

/*!
 * \brief   A data structure for handling multi-threaded loading of
 *          two-dimensional probability distributions.
 *
 * \ingroup two_dimensional_distribution_loader
 */
typedef struct {
  /*!
   * \brief   Pointers to the distributions.
   */
  Distribution ** distributions;

  /*!
   * \brief   The paths to the distributions to be loaded.
   */
  char ** paths;

  /*!
   * \brief   The total number of distributions to be loaded.
   */
  uint32_t count;

  /*!
   * \brief   The offset of the distribution currently being loaded within
   *          Distribution_Loader::distributions.
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
} Distribution_Loader;

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
 * \param[in] count         The number of distributions to be loaded.
 */
void distribution_loader_init(
  Distribution_Loader * const loader,
  char ** paths,
  const uint32_t count);

/*!
 * \brief Clears the distribution loader.
 *
 * \param[in, out] loader   An initialized distribution loader to be cleared.
 */
void distribution_loader_clear(
  Distribution_Loader * const loader);

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
Distribution * distribution_loader_pop(
  Distribution_Loader * const loader);

/*!
 * \}
 */

#endif /* DISTRIBUTION_LOADER_H */
