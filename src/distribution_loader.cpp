/*!
 * \file    distribution_loader.cpp
 * \ingroup two_dimensional_distribution_loader
 *
 * \brief   The definition of functions for loading two-dimensional probability
 *          distributions.
 */

#include "distribution_loader.h"
#include "distribution.h"
#include "thread_pool.h"
#include "errors.h"
#include "string_utilities.h"

#include <pthread.h>

#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static void * main_load_distributions(void * ptr)
{
  /* Map the input. */
  Distribution_Loader * const loader = (Distribution_Loader *)ptr;

  /* Load distributions. */
  const uint32_t count = loader->count;

  for (uint32_t i = 0; i < count; i++) {
    /* Check if we should wait for distributions to be popped before loading. */
    while (TRUE) {
      pthread_mutex_lock(&(loader->mutex));
      const uint32_t offset = loader->offset;
      pthread_mutex_unlock(&(loader->mutex));
      
      if ((i - offset) < 
        DISTRIBUTION_LOADER_MAX_SIMULTANEOUSLY_LOADED_DISTRIBUTIONS)
      {
        break;
      }
      
      sleep(1);
    }

    printf("Loading distribution \"%s\"...\n", loader->paths[i]);

    Distribution * distribution = distribution_alloc();

    FILE * file = fopen(loader->paths[i], "rb");
    if (NULL == file) {
      critical("main_load_distributions(): "
        "Failed to read distribution from file \"%s\".",
          loader->paths[i]);
    }

    distribution_init_import(distribution, file);

    fclose(file);
    file = NULL;

    printf("Finished loading distribution from file ""\"%s\".\n",
      loader->paths[i]);

    /* Insert the distribution into the distribution loader. */
    pthread_mutex_lock(&(loader->mutex));
    loader->distributions[i] = distribution;
    pthread_mutex_unlock(&(loader->mutex));

    distribution = NULL;
  }

  /* Finished. */
  return NULL;
}

void distribution_loader_init(
  Distribution_Loader * const loader,
  char ** paths,
  const uint32_t count)
{
  /* Start off by zeroizing the entire data structure. */
  memset(loader, 0, sizeof(Distribution_Loader));

  /* Allocate memory for the paths array. */
  loader->paths = (char **)malloc(count * sizeof(char *));
  if (NULL == loader->paths) {
    critical("distribution_loader_init(): Failed to allocate memory.");
  }

  /* Copy the paths array. */
  for (uint32_t i = 0; i < count; i++) {
    size_t length = strlen(paths[i]);

    loader->paths[i] = (char *)malloc((length + 1) * sizeof(char));
    if (NULL == loader->paths[i]) {
      critical("distribution_loader_init(): Failed to allocate memory.");
    }

    safe_strlcpy(loader->paths[i], paths[i], length + 1);
  }

  /* Allocate memory for the distributions array. */
  loader->distributions =
    (Distribution **)malloc(count * sizeof(Distribution *));
  if (NULL == loader->distributions) {
    critical("distribution_loader_init(): Failed to allocate memory.");
  }

  /* Initialize the distributions array. */
  for (uint32_t i = 0; i < count; i++) {
    loader->distributions[i] = NULL;
  }

  /* Store the count and setup other fields. */
  loader->count = count;
  loader->offset = 0;

  /* Initialize the mutex. */
  pthread_mutex_init(&(loader->mutex), NULL);

  /* Initialize the thread pool. */
  thread_pool_init(&(loader->pool));

  /* Spawn a worker thread to load the distributions. */
  thread_pool_spawn(&(loader->pool), main_load_distributions, loader);
}

Distribution * distribution_loader_pop(
  Distribution_Loader * const loader)
{
  while (TRUE) {
    pthread_mutex_lock(&(loader->mutex));

    if ((loader->offset) >= (loader->count)) {
      /* The distribution loader has been exhausted. */
      pthread_mutex_unlock(&(loader->mutex));
      return NULL;
    }

    Distribution * distribution = loader->distributions[loader->offset];
    if (NULL != distribution) {
      printf("Popped distribution \"%s\".\n", loader->paths[loader->offset]);

      loader->distributions[loader->offset] = NULL;
      loader->offset += 1;
    }

    pthread_mutex_unlock(&(loader->mutex));

    if (NULL != distribution) {
      return distribution;
    }

    sleep(1); /* Sleep for one second before trying again... */
  }
}

void distribution_loader_clear(
  Distribution_Loader * const loader)
{
  /* Join the worker thread to the main thread and clear the thread pool. */
  thread_pool_join(&(loader->pool));
  thread_pool_clear(&(loader->pool));

  /* Destroy the mutex. */
  pthread_mutex_destroy(&(loader->mutex));

  for (uint32_t i = 0; i < (loader->count); i++) {
    if (NULL != loader->distributions[i]) {
      distribution_clear(loader->distributions[i]);
      distribution_dealloc(&(loader->distributions[i]));
    }

    if (NULL != loader->paths[i]) {
      free(loader->paths[i]);
      loader->paths[i] = NULL;
    }
  }

  free(loader->distributions);
  loader->distributions = NULL;

  free(loader->paths);
  loader->paths = NULL;

  memset(loader, 0, sizeof(Distribution_Loader));
}
