/*!
 * \file    diagonal_distribution_loader.cpp
 * \ingroup diagonal_distribution_loader
 *
 * \brief   The definition of functions for loading diagonal probability
 *          distributions.
 */

#include "diagonal_distribution_loader.h"

#include "common.h"
#include "diagonal_distribution.h"
#include "errors.h"
#include "string_utilities.h"
#include "thread_pool.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>
#include <unistd.h>

static void * main_load_diagonal_distributions(void * ptr)
{
  /* Map the input. */
  Diagonal_Distribution_Loader * const loader = 
    (Diagonal_Distribution_Loader *)ptr;

  /* Load distributions. */
  const uint32_t count = loader->count;

  for (uint32_t i = 0; i < count; i++) {
    /* Check if we should wait for distributions to be popped before loading. */
    while (TRUE) {
      pthread_mutex_lock(&(loader->mutex));
      const uint32_t offset = loader->offset;
      pthread_mutex_unlock(&(loader->mutex));
      
      if ((i - offset) < 
        DIAGONAL_DISTRIBUTION_LOADER_MAX_SIMULTANEOUSLY_LOADED_DISTRIBUTIONS)
      {
        break;
      }
      
      sleep(1);
    }
    
    printf("Loading distribution \"%s\"...\n", loader->paths[i]);

    Diagonal_Distribution * distribution = diagonal_distribution_alloc();

    FILE * file = fopen(loader->paths[i], "rb");
    if (NULL == file) {
      critical("main_load_diagonal_distributions(): "
        "Failed to read distribution from file \"%s\".",
          loader->paths[i]);
    }

    diagonal_distribution_init_import(distribution, file);

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

void diagonal_distribution_loader_init(
  Diagonal_Distribution_Loader * const loader,
  char ** paths,
  const uint32_t count)
{
  /* Start off by zeroizing the entire data structure. */
  memset(loader, 0, sizeof(Diagonal_Distribution_Loader));

  /* Allocate memory for the paths array. */
  loader->paths = (char **)malloc(count * sizeof(char *));
  if (NULL == loader->paths) {
    critical("diagonal_distribution_loader_init(): Failed to allocate memory.");
  }

  /* Copy the paths array. */
  for (uint32_t i = 0; i < count; i++) {
    size_t length = strlen(paths[i]);

    loader->paths[i] = (char *)malloc((length + 1) * sizeof(char));
    if (NULL == loader->paths[i]) {
      critical("diagonal_distribution_loader_init(): "
        "Failed to allocate memory.");
    }

    safe_strlcpy(loader->paths[i], paths[i], length + 1);
  }

  /* Allocate memory for the distributions array. */
  loader->distributions =
    (Diagonal_Distribution **)malloc(count * sizeof(Diagonal_Distribution *));
  if (NULL == loader->distributions) {
    critical("diagonal_distribution_loader_init(): Failed to allocate memory.");
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
  thread_pool_spawn(&(loader->pool), main_load_diagonal_distributions, loader);
}

Diagonal_Distribution * diagonal_distribution_loader_pop(
  Diagonal_Distribution_Loader * const loader)
{
  while (TRUE) {
    pthread_mutex_lock(&(loader->mutex));

    if ((loader->offset) >= (loader->count)) {
      /* The distribution loader has been exhausted. */
      pthread_mutex_unlock(&(loader->mutex));
      return NULL;
    }

    Diagonal_Distribution * distribution = 
      loader->distributions[loader->offset];
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

void diagonal_distribution_loader_clear(
  Diagonal_Distribution_Loader * const loader)
{
  /* Join the worker thread to the main thread and clear the thread pool. */
  thread_pool_join(&(loader->pool));
  thread_pool_clear(&(loader->pool));

  /* Destroy the mutex. */
  pthread_mutex_destroy(&(loader->mutex));

  for (uint32_t i = 0; i < (loader->count); i++) {
    if (NULL != loader->distributions[i]) {
      diagonal_distribution_clear(loader->distributions[i]);
      diagonal_distribution_dealloc(&(loader->distributions[i]));
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

  memset(loader, 0, sizeof(Diagonal_Distribution_Loader));
}
