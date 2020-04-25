/*!
 * \file    thread_pool.h
 * \ingroup thread_pool
 *
 * \brief   The declaration of data structures for implementing thread pools,
 *          and of functions for manipulating such pools.
 */

/*!
 * \defgroup thread_pool  Thread pools
 * \ingroup  utility
 *
 * \brief    A module for thread pools.
 */

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <pthread.h>

#include <stdint.h>

/*!
 * \brief   The maximum capacity of the thread pool in threads.
 */
#define THREAD_POOL_CAPACITY                        16

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief   A data structure for representing a thread pool.
 */
typedef struct {
  /*!
   * \brief   Handles for the threads in the pool.
   */
  pthread_t threads[THREAD_POOL_CAPACITY];

  /*!
   * \brief   Indices for the threads in the pool.
   */
  uint32_t indices[THREAD_POOL_CAPACITY];

  /*!
   * \brief   A counter used to assign a running index to each thread spawned.
   *
   * Initially assigned index one and incremented for each thread spawned by
   * the thread_pool_spawn(). If the counter overflows, this will be result in
   * a critical error.
   */
  uint32_t index;
} Thread_Pool;

/*!
 * \name Initialization
 * \{
 */

/*!
 * \brief   Initializes the thread pool.
 *
 * This function assigns all thread handles an index of zero.
 *
 * \param[in, out] pool       The thread pool to initialize.
 */
void thread_pool_init(
  Thread_Pool * const pool);

/*!
 * \brief   Clears the thread pool.
 *
 * This function calls thread_pool_join() as a part of clearing the pool.
 *
 * \param[in, out] pool       The thread pool to clear.
 */
void thread_pool_clear(
  Thread_Pool * const pool);

/*!
 * \}
 */

/*!
 * \name Spawning and joining
 * \{
 */

/*!
 * \brief   Spawns a new thread from the pool and calls a function with a given
 *          argument from the thread spawned.
 *
 * The first available handle is used to spawn the thread. If no thread handle
 * has a zero index, the thread handle with the lowest index is first joined
 * with the thread used to call this function and then re-used to spawn a new
 * thread.
 *
 * This is implemented as follows: All thread handles are initially assigned a
 * zero index by thread_pool_init(). This function iterates through the list of
 * handles in the pool and uses the first handle with a zero index to spawn a
 * new thread. If no thread handle has a zero index, the thread handle with the
 * lowest index is first joined with the thread used to call this function and
 * then re-used to spawn a new thread.
 *
 * The thread handle used to spawned a new thread is assigned an index equal to
 * the counter Thread_Pool::index. The counter is the incremented. If the
 * counter overflows, this will be result in a critical error.
 *
 * \param[in, out] pool       The thread pool.
 * \param[in] start_routine   The function to call from the worker thread.
 * \param[in] arg             A pointer that the worker thread should pass to
 *                            the function.
 */
void thread_pool_spawn(
  Thread_Pool * const pool,
  void * (* start_routine)(void *),
  void * arg);

/*!
 * \brief   Joins all threads in the pool with the main thread.
 *
 * \param[in, out] pool       The thread pool.
 */
void thread_pool_join(
  Thread_Pool * const pool);

/*!
 * \}
 */

#ifdef __cplusplus
}
#endif

#endif /* THREAD_POOL_H */
