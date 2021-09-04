/*!
 * \file    thread_pool.c
 * \ingroup thread_pool
 *
 * \brief   The definition of functions for implementing thread pools.
 */

#include "thread_pool.h"

#include "errors.h"

#include <stdint.h>
#include <string.h>

#include <pthread.h>

void thread_pool_init(
  Thread_Pool * const pool)
{
  for (uint32_t i = 0; i < THREAD_POOL_CAPACITY; i++) {
    pool->indices[i] = 0;
  }

  pool->index = 1;
}

void thread_pool_clear(
  Thread_Pool * const pool)
{
  thread_pool_join(pool);

  memset(pool, 0, sizeof(Thread_Pool));
}

void thread_pool_spawn(
  Thread_Pool * const pool,
  void * (* start_routine)(void *),
  void * arg)
{
  for (uint32_t i = 0; i < THREAD_POOL_CAPACITY; i++) {
    if (0 == pool->indices[i]) {
      if (0 != pthread_create(
        &(pool->threads[i]),
        NULL, /* attr */
        start_routine,
        arg))
      {
        critical("thread_pool_spawn(): Failed to spawn thread.");
      }

      pool->indices[i] = pool->index;

      pool->index++;
      if (0 == pool->index) {
        critical("thread_pool_spawn(): The thread index counter is exhausted.");
      }

      return;
    }
  }

  /* The thread pool is at maximum capacity. Join the oldest train in the thread
   * pool with the main thread to free up a slot in the pool. */

  /* Find the oldest thread in the pool. */
  uint32_t min_index = 0xFFFFFFFFU;
  uint32_t min_i = 0;

  for (uint32_t i = 0; i < THREAD_POOL_CAPACITY; i++) {
    if ((0 != pool->indices[i]) && (min_index >= pool->indices[i])) {
      min_index = pool->indices[i];
      min_i = i;
    }
  }

  /* Join the oldest thread in the pool with the main thread. */
  if (0 != pthread_join(pool->threads[min_i], NULL)) {
    critical("thread_pool_spawn(): Failed to join threads.");
  }

  pool->indices[min_i] = 0;

  /* Create a new thread. */
  if (0 != pthread_create(
    &(pool->threads[min_i]),
    NULL, /* attr */
    start_routine,
    arg))
  {
    critical("thread_pool_spawn(): Failed to spawn thread.");
  }

  pool->indices[min_i] = pool->index;

  pool->index++;
  if (0 == pool->index) {
    critical("thread_pool_spawn(): The thread index counter is exhausted.");
  }
}

void thread_pool_join(
  Thread_Pool * const pool)
{
  for (uint32_t i = 0; i < THREAD_POOL_CAPACITY; i++) {
    if (0 != pool->indices[i]) {
      if (0 != pthread_join(pool->threads[i], NULL)) {
        critical("thread_pool_join(): Failed to join threads.");
      }

      pool->indices[i] = 0;
    }
  }
}
