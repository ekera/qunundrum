/*!
 * \file    main_generate_linear_distribution_rsa.cpp
 * \ingroup generate_linear_distribution_rsa_exe
 *
 * \brief   The definition of the main entry point to the
 *          generate_linear_distribution_rsa executable, and of associated
 *          functions.
 */

/*!
 * \defgroup generate_linear_distribution_rsa_exe \
 *           The generate_linear_distribution_rsa executable
 * \ingroup  generate_executable
 *
 * \brief    A module for the generate_linear_distribution_rsa executable.
 */

#include "executables.h"
#include "executables_generate_distribution.h"

#include "linear_distribution.h"
#include "linear_distribution_slice.h"
#include "linear_distribution_enumerator.h"
#include "linear_probability.h"
#include "parameters.h"
#include "parameters_selection.h"
#include "thread_pool.h"
#include "math.h"
#include "errors.h"
#include "string_utilities.h"
#include "rsa.h"

#include <gmp.h>
#include <mpfr.h>
#include <mpi.h>

#include <pthread.h>

#include <memory.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <unistd.h>

#include <sys/stat.h>
#include <sys/types.h>

/*!
 * \brief   A data structure representing argument entries in the form of
 *          parsed \<n\> values from the command line arguments.
 *
 * \ingroup generate_linear_distribution_rsa_exe
 */
typedef struct {
  /*!
   * \brief   The parameter m.
   */
  uint32_t m;

  /*!
   * \brief   The parameter l.
   */
  uint32_t l;

  /*!
   * \brief   The parameter l.
   */
  uint32_t n;

  /*!
   * \brief   The value of d for the selection method.
   *
   * The selection method is common for all argument tuples, see
   * Generate_Linear_Distribution_RSA_Arguments::selection_method.
   */
  mpz_t d;

  /*!
   * \brief   The value of p for the selection method, if applicable, or zero.
   *
   * The selection method is common for all argument tuples, see
   * Generate_Linear_Distribution_RSA_Arguments::selection_method.
   */
  mpz_t p;

  /*!
   * \brief   The value of q for the selection method, if applicable, or zero.
   *
   * The selection method is common for all argument tuples, see
   * Generate_Linear_Distribution_RSA_Arguments::selection_method.
   */
  mpz_t q;

  /*!
   * \brief   The path to which to export this distribution.
   */
  char path[MAX_SIZE_PATH_BUFFER];
} Generate_Linear_Distribution_RSA_Arguments_Entry;

/*!
 * \brief   A data structure representing parsed command line arguments.
 *
 * \ingroup generate_linear_distribution_rsa_exe
 */
typedef struct {
  /*!
   * \brief   The selection method for the logarithm d or order r.
   */
  Selection_Method selection_method;

  /*!
   * \brief   The value of p for the selection method, if applicable, or zero.
   */
  mpz_t p;

  /*!
   * \brief   The value of q for the selection method, if applicable, or zero.
   */
  mpz_t q;

  /*!
   * \brief   The value of the parameter delta.
   */
  uint32_t delta;

  /*!
   * \brief   The dimension of the slices in the distribution.
   */
  uint32_t dimension;

  /*!
   * \brief   The number of distributions to generate.
   */
  uint32_t count;

  /*!
   * \brief   A vector of entries for each distribution to generate.
   */
  Generate_Linear_Distribution_RSA_Arguments_Entry * entries;
} Generate_Linear_Distribution_RSA_Arguments;

/*!
 * \brief   A data structure representing an export job.
 *
 * \ingroup generate_linear_distribution_rsa_exe
 */
typedef struct {
  /*!
   * \brief   A pointer to the linear distribution to export.
   */
  Linear_Distribution * distribution;

  /*!
   * \brief   The path to which to export the distribution.
   */
  char path[MAX_SIZE_PATH_BUFFER];
} Linear_Distribution_RSA_Export_Job;

/*!
 * \brief   Parses the command line arguments.
 *
 * \param[in, out] arguments   The argument data structure in which to store
 *                             the parsed command line arguments.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \remark   So as to allow parallelized executables to be shut down gracefully,
 *           this function does not call critical() to signal a critical error.
 *           Rather it prints an informative error message to stderr and returns
 *           False. When returning False, this function expects the caller to
 *           terminate the executable.
 *
 * \return   True if the command line arguments were successfully parsed, False
 *           otherwise. If False is returned, the data structure may be only
 *           partially initialized and memory only partially allocated.
 */
static bool arguments_init_parse_command_line(
  Generate_Linear_Distribution_RSA_Arguments * const arguments,
  const int argc,
  char ** argv)
{
  /* Initialize the arguments data structure. */
  arguments->selection_method = SELECTION_METHOD_UNDEFINED;
  arguments->delta = 0;
  arguments->dimension = 0;
  arguments->entries = NULL;
  arguments->count = 0;

  mpz_init(arguments->p);
  mpz_set_ui(arguments->p, 0);
  mpz_init(arguments->q);
  mpz_set_ui(arguments->q, 0);

  /* Iterate over the command line arguments. */
  int i;

  for (i = 1; i < argc; i++) {
    /* Parse the selection method. */
    Selection_Method selection_method = SELECTION_METHOD_UNDEFINED;

    if (0 == strcmp(argv[i], "-max")) {
      selection_method = SELECTION_METHOD_MAXIMAL;
    } else if (0 == strcmp(argv[i], "-rnd")) {
      selection_method = SELECTION_METHOD_RANDOM;
    } else if (0 == strcmp(argv[i], "-exp")) {
      selection_method = SELECTION_METHOD_EXPLICIT;
    }

    if (SELECTION_METHOD_UNDEFINED != selection_method) {
      /* Check that a selection method has not already been specified. */
      if (SELECTION_METHOD_UNDEFINED != (arguments->selection_method)) {
        fprintf(stderr, "Error: The selection method cannot be twice "
          "specified.\n");
        return FALSE;
      }

      /* Store the selection method. */
      arguments->selection_method = selection_method;

      if (SELECTION_METHOD_EXPLICIT == selection_method) {
        if ((i + 2) >= argc) {
          fprintf(stderr, "Error: Expected <p> <q> to follow after -exp.\n");
          return FALSE;
        }

        /* Parse p. */
        if (0 != mpz_set_str(arguments->p, argv[i + 1], 10)) {
          fprintf(stderr, "Error: Failed to parse <p> after -exp as an "
            "integer.\n");
          return FALSE;
        }

        if (mpz_cmp_ui(arguments->p, 0) <= 0) {
          fprintf(stderr, "Error: The value of <p> after -exp must be "
            "positive.\n");
          return FALSE;
        }

        /* Parse q. */
        if (0 != mpz_set_str(arguments->q, argv[i + 2], 10)) {
          fprintf(stderr, "Error: Failed to parse <q> after -exp as an "
            "integer.\n");
          return FALSE;
        }

        if (mpz_cmp_ui(arguments->q, 0) <= 0) {
          fprintf(stderr, "Error: The value of <q> after -exp must be "
            "positive.\n");
          return FALSE;
        }

        /* Check that p and q are distinct. */
        if (mpz_cmp(arguments->p, arguments->q) == 0) {
          fprintf(stderr,
            "Error: The values of <q> and <q> after -exp must be dinstict.\n");
          return FALSE;
        }

        /* Check that p and q are of the same length in bits. */
        if (mpz_sizeinbase(arguments->p, 2) !=
            mpz_sizeinbase(arguments->q, 2))
        {
          fprintf(stderr,
            "Error: The values of <p> and <q> after -exp must be of the same "
            "length in bits when using this convenience utility. For other "
            "<p> and <q>, use the general linear distribution generator.\n");
          return FALSE;
        }

        /* Check that p and q are (probable) prime numbers. */
        int result;

        result = mpz_probab_prime_p(arguments->p, RSA_MILLER_RABIN_ITERATIONS);
        /* Note: Returns 2 if n is prime, return 1 if n is probably prime 
         * (without being certain), or return 0 if n is definitely non-prime. */
        if ((1 != result) && (2 != result)) {
          fprintf(stderr,
             "Error: The value of <p> after -exp must be prime.\n");
          return FALSE;
        }

        result = mpz_probab_prime_p(arguments->q, RSA_MILLER_RABIN_ITERATIONS);
        /* Note: Returns 2 if n is prime, return 1 if n is probably prime 
         * (without being certain), or return 0 if n is definitely non-prime. */
        if ((1 != result) && (2 != result)) {
          fprintf(stderr,
             "Error: The value of <q> after -exp must be prime.\n");
          return FALSE;
        }

        i += 2;
      }

      continue;
    }

    /* Parse delta. */
    if (0 == strcmp(argv[i], "-delta")) {
      /* Check that delta has not already been specified. */
      if (0 != arguments->delta) {
        fprintf(stderr, "Error: The delta cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr, "Error: Expected value to follow after -delta.\n");
        return FALSE;
      }

      int x = atoi(argv[i+1]);
      if (x <= 0) {
        fprintf(stderr, "Error: The value of <delta> must be positive.\n");
        return FALSE;
      }

      /* Store the delta. */
      arguments->delta = (uint32_t)x;

      i++;

      continue;
    }

    /* Parse the dimension. */
    if (0 == strcmp(argv[i], "-dim")) {
      /* Check that a dimension has not already been specified. */
      if (0 != arguments->dimension) {
        fprintf(stderr, "Error: The dimension cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr, "Error: Expected <dimension> to follow after -dim.\n");
        return FALSE;
      }

      int x = atoi(argv[i+1]);
      if ((x < 16) || (x > 8192) || (!is_pow2((uint32_t)x))) {
        fprintf(stderr, "Error: The <dimension> passed to -dim must be a power "
          "of two on the interval [16, 8192].\n");
        return FALSE;
      }

      /* Store the dimension. */
      arguments->dimension = (uint32_t)x;

      i++;

      continue;
    }

    /* Stop parsing. */
    break;
  }

  /* Set default parameters if arguments where not explicitly specified. */
  if (0 == arguments->delta) {
    arguments->delta = 20;
  }

  if (0 == arguments->dimension) {
    arguments->dimension = 2048;
  }

  if (SELECTION_METHOD_UNDEFINED == arguments->selection_method) {
    arguments->selection_method = SELECTION_METHOD_MAXIMAL;
  }

  /* Parse entries { <n> }. */
  if ((argc - i) <= 0) {
    fprintf(stderr, "Error: Incorrect command line arguments; expected tuples "
      "but found an odd number of arguments.\n");
    return FALSE;
  }

  arguments->count = (uint32_t)(argc - i);

  arguments->entries =
    (Generate_Linear_Distribution_RSA_Arguments_Entry *)malloc(
      (arguments->count) * sizeof(
        Generate_Linear_Distribution_RSA_Arguments_Entry));
  if (NULL == arguments->entries) {
    fprintf(stderr, "Error: Failed to allocate memory for argument entries.\n");
    return FALSE;
  }

  for (uint32_t j = 0; j < arguments->count; j++) {
    mpz_init(arguments->entries[j].d);
    mpz_set_ui(arguments->entries[j].d, 0);

    mpz_init(arguments->entries[j].p);
    mpz_set_ui(arguments->entries[j].p, 0);

    mpz_init(arguments->entries[j].q);
    mpz_set_ui(arguments->entries[j].q, 0);
  }

  /* Check that the output directory exists and is accessible. */
  if (0 != access(DISTRIBUTIONS_DIRECTORY, F_OK)) {
    if (0 != mkdir(DISTRIBUTIONS_DIRECTORY, DIRECTORY_PERMISSIONS)) {
      fprintf(stderr, "Error: The output directory \"%s\" does not exist.\n",
        DISTRIBUTIONS_DIRECTORY);
      return FALSE;
    }
  }

  if (0 != access(DISTRIBUTIONS_DIRECTORY, R_OK | W_OK | X_OK)) {
    fprintf(stderr, "Error: The output directory \"%s\" is not readable, "
      "writeable and executable.\n", DISTRIBUTIONS_DIRECTORY);
    return FALSE;
  }

  /* Iterate over the tuples. */
  for (uint32_t j = 0; j < arguments->count; j++, i++) {
    int x = atoi(argv[i]);
    if (x <= 0) {
      fprintf(stderr, "Error: Failed to parse <n>.\n");
      return FALSE;
    }

    const uint32_t n = (uint32_t)x;

    if (((n % 2) != 0) || (n < 128)) {
      fprintf(stderr, "Error: The value of <n> must be even and >= 128.\n");
      return FALSE;
    }

    if (n < 2 * (arguments->delta + 1)) {
      fprintf(stderr, "Error: The value of <n> must be > 2(delta + 1).\n");
      return FALSE;
    }

    arguments->entries[j].n = n;
    arguments->entries[j].m = n / 2 - 1;
    arguments->entries[j].l = arguments->entries[j].m - arguments->delta;

    char suffix[MAX_SIZE_PATH_BUFFER];
    safe_snprintf(suffix,
      MAX_SIZE_PATH_BUFFER,
      "dim-%u-rsa-n-%u-delta-%u-m-%u-l-%u.txt",
      arguments->dimension,
      arguments->entries[j].n,
      arguments->delta,
      arguments->entries[j].m,
      arguments->entries[j].l);

    /* Set the value. */
    if (SELECTION_METHOD_MAXIMAL == arguments->selection_method) {
      mpz_setbit(arguments->entries[j].d, arguments->entries[j].m);
      mpz_sub_ui(arguments->entries[j].d, arguments->entries[j].d, 1);

      safe_snprintf(arguments->entries[j].path, MAX_SIZE_PATH_BUFFER,
        "%s/linear-distribution-max-%s", DISTRIBUTIONS_DIRECTORY, suffix);
    } else if (SELECTION_METHOD_RANDOM == arguments->selection_method) {
      parameters_selection_random_rsa(
        arguments->entries[j].d,
        arguments->entries[j].p,
        arguments->entries[j].q,
        arguments->entries[j].n);

      safe_snprintf(arguments->entries[j].path, MAX_SIZE_PATH_BUFFER,
        "%s/linear-distribution-rnd-%s", DISTRIBUTIONS_DIRECTORY, suffix);
    } else if (SELECTION_METHOD_EXPLICIT == arguments->selection_method) {
      if ((arguments->entries[j].n != 2 * mpz_sizeinbase(arguments->p, 2)) ||
          (arguments->entries[j].n != 2 * mpz_sizeinbase(arguments->q, 2)))
      {
        fprintf(stderr, "Error: The value of <n> is incompatible with the "
          "explicitly specified values of <p> and <q> using the -exp flag.\n");
        return FALSE;
      }

      mpz_set(arguments->entries[j].p, arguments->p);
      mpz_set(arguments->entries[j].q, arguments->q);

      parameters_selection_explicit_rsa(
        arguments->entries[j].d,
        arguments->entries[j].p,
        arguments->entries[j].q,
        arguments->entries[j].n);

      safe_snprintf(arguments->entries[j].path, MAX_SIZE_PATH_BUFFER,
        "%s/linear-distribution-exp-%s", DISTRIBUTIONS_DIRECTORY, suffix);
    } else {
      fprintf(stderr, "Error: Internal error: Unknown selection method.\n");
      return FALSE;
    }

    /* Check that the output file does not already exist. */
    if (0 == access(arguments->entries[j].path, F_OK)) {
      fprintf(stderr, "Error: The distribution \"%s\" already exists.\n",
        arguments->entries[j].path);
      return FALSE;
    }
  }

  /* Signal success. */
  return TRUE;
}

/*!
 * \brief   Clears an initialized command line arguments data structure.
 * 
 * \param[in, out] arguments   The argument data structure to clear.
 */
static void arguments_clear(
  Generate_Linear_Distribution_RSA_Arguments * const arguments)
{
  mpz_clear(arguments->p);
  mpz_clear(arguments->q);

  if (NULL != arguments->entries) {
    for (uint32_t i = 0; i < arguments->count; i++) {
      mpz_clear(arguments->entries[i].d);
      mpz_clear(arguments->entries[i].p);
      mpz_clear(arguments->entries[i].q);
    }

    free(arguments->entries);
    arguments->entries = NULL;
  }

  memset(arguments, 0, sizeof(Generate_Linear_Distribution_RSA_Arguments));
}

/*!
 * \brief   Exports the computed probability distribution to file.
 *
 * This function is intended to run in a separate thread.
 *
 * \param[in, out] ptr    A pointer to a Linear_Distribution_RSA_Export_Job data
 *                        structure that describes the export job.
 *
 * \remark   The export job and the distribution it contains will be cleared
 *           and deallocated by this function.
 */
void * main_server_export_distribution(void * ptr)
{
  /* Map the input to data structures. */
  Linear_Distribution_RSA_Export_Job * job =
    (Linear_Distribution_RSA_Export_Job *)ptr;

  Linear_Distribution * distribution = job->distribution;
  const char * path = job->path;

  /* Sort the slices in the distribution. */
  printf("Sorting the slices in the distribution...\n");
  linear_distribution_sort_slices(distribution);

  /* Export, clear and deallocate the distribution. */
  printf("Exporting the distribution to \"%s\"...\n", path);

  FILE * file = fopen(path, "wb+");
  if (NULL == file) {
    critical("main_server_export_distribution(): Failed to open \"%s\".", path);
  }

  linear_distribution_export(distribution, file);
  fclose(file);

  linear_distribution_clear(distribution);
  linear_distribution_dealloc(&distribution);

  printf("Finished exporting the distribution to \"%s\".\n", path);

  /* Finished. */
  return NULL;
}

/*!
 * \brief   The main function on the server node.
 *
 * This function is called once by main() for each distribution to generate.
 *
 * \param[in] arguments   The parsed command line arguments.
 * \param[in] entry       The entry in the command lines arguments data
 *                        structure corresponding to the distribution to
 *                        generate.
 * \param[in] mpi_size    The number of nodes.
 * \param[in] pool        The thread pool from which to spawn threads for
 *                        export operations that run in the background.
 */
void main_server(
  const Generate_Linear_Distribution_RSA_Arguments * const arguments,
  const Generate_Linear_Distribution_RSA_Arguments_Entry * const entry,
  int mpi_size,
  Thread_Pool * const pool)
{
  /* Setup the parameters. */
  const uint32_t t = 30;

  Parameters parameters;
  parameters_init(&parameters);

  parameters_explicit_m_l(
    &parameters,
    entry->d, /* d */
    entry->d, /* r */
    entry->m,
    entry->l,
    t);

  /* Send broadcast of the parameters. */
  parameters_bcast_send(&parameters, MPI_RANK_ROOT);

  /* Send broadcast of the dimension. */
  uint32_t dimension = arguments->dimension;
  if (MPI_SUCCESS !=
    MPI_Bcast(&dimension,
      1, /* count */
      MPI_UNSIGNED,
      MPI_RANK_ROOT,
      MPI_COMM_WORLD))
  {
    critical("main_server(): Failed to send broadcast of the dimension.");
  }

  /* Setup the distribution enumerator. */
  const bool mirrored = TRUE;

  Linear_Distribution_Enumerator enumerator;
  linear_distribution_enumerator_init(&enumerator, &parameters, mirrored);

  /* Get the required distribution capacity in slices. */
  uint32_t capacity = linear_distribution_enumerator_count(&enumerator);

  if (mirrored) {
    capacity *= 2;
  }

  /* Setup the distribution. */
  uint32_t flags = LINEAR_DISTRIBUTION_FLAG_COMPUTED;
  flags |= LINEAR_DISTRIBUTION_FLAG_D;
  flags |= LINEAR_DISTRIBUTION_FLAG_RSA;

  if (SELECTION_METHOD_EXPLICIT == arguments->selection_method) {
    flags |= LINEAR_DISTRIBUTION_FLAG_EXPLICIT;
  } else if (SELECTION_METHOD_RANDOM == arguments->selection_method) {
    flags |= LINEAR_DISTRIBUTION_FLAG_RANDOM;
  } else if (SELECTION_METHOD_MAXIMAL == arguments->selection_method) {
    flags |= LINEAR_DISTRIBUTION_FLAG_MAXIMAL;
  } else {
    critical("main_server(): Unknown selection method.");
  }

  Linear_Distribution * distribution = linear_distribution_alloc();
  linear_distribution_init(distribution, &parameters, flags, capacity);

  /* Keep track of the number of nodes remaining. */
  int mpi_nodes_remaining = mpi_size;

  while (mpi_nodes_remaining > 1) {
    /* Listen for a notification. */
    MPI_Status status;

    uint32_t notification;

    if (MPI_SUCCESS != MPI_Recv(
        &notification,
        1, /* count */
        MPI_UNSIGNED,
        MPI_ANY_SOURCE, /* source */
        MPI_TAG_NOTIFY,
        MPI_COMM_WORLD,
        &status))
    {
      critical("main_server(): Failed to receive notification.");
    }

    if (MPI_NOTIFY_READY == notification) {
      int32_t min_log_alpha;

      bool found = linear_distribution_enumerator_next(
                      &min_log_alpha,
                      &enumerator);
      if (found) {
        printf("Processing slice: %u / %u (%d)\n",
          enumerator.offset,
          enumerator.count,
          min_log_alpha);

        uint32_t job = MPI_JOB_PROCESS_SLICE;

        if (MPI_SUCCESS != MPI_Send(
            &job,
            1, /* count */
            MPI_UNSIGNED,
            status.MPI_SOURCE, /* destination */
            MPI_TAG_JOB,
            MPI_COMM_WORLD))
        {
          critical("main_server(): Failed to send job.");
        }

        if (MPI_SUCCESS != MPI_Send(
            &min_log_alpha,
            1, /* count */
            MPI_INT,
            status.MPI_SOURCE, /* destination */
            MPI_TAG_SLICE_ALPHAS,
            MPI_COMM_WORLD))
        {
          critical("main_server(): Failed to send MPI_TAG_SLICE_ALPHAS.");
        }
      } else {
        uint32_t job = MPI_JOB_STOP;

        if (MPI_SUCCESS != MPI_Send(
            &job,
            1, /* count */
            MPI_UNSIGNED,
            status.MPI_SOURCE, /* destination */
            MPI_TAG_JOB,
            MPI_COMM_WORLD))
        {
          critical("main_server(): Failed to send job.");
        }

        mpi_nodes_remaining--;
      }
    } else if (MPI_NOTIFY_SLICE_DONE == notification) {
      /* The node has finished processing a slice. Receive the slice. */
      Linear_Distribution_Slice * slice = linear_distribution_slice_alloc();

      linear_distribution_slice_init_recv(slice, status.MPI_SOURCE);

      linear_distribution_insert_slice(distribution, slice);

      if (mirrored) {
        Linear_Distribution_Slice * mirrored_slice =
          linear_distribution_slice_alloc();
        linear_distribution_slice_init_copy(mirrored_slice, slice);

        mirrored_slice->min_log_alpha *= -1; /* Mirror alpha. */

        mirrored_slice->flags |= SLICE_FLAGS_MIRRORED;

        linear_distribution_insert_slice(distribution, mirrored_slice);
      }

      /* Print some feedback to the user. */
      printf("Slice total probability is: %.24LG\n", slice->total_probability);
      printf("Slice dimension is: %u\n", slice->dimension);

      printf("Total probability is now: %.24LG\n",
        distribution->total_probability);

      /* Note: The slice is inserted by reference and will be deallocated by
       *       the call to distribution_dealloc(). It must not be deallocated
       *       here as doing so would lead to memory corruption. */
    } else if (MPI_NOTIFY_SLICE_SKIP == notification) {
      /* Ignore. */
    } else {
      critical("main_server(): Unknown notification.");
    }
  }

  /* Clear memory. */
  linear_distribution_enumerator_clear(&enumerator);

  parameters_clear(&parameters);

  /* Setup an export job. */
  Linear_Distribution_RSA_Export_Job * export_job =
    (Linear_Distribution_RSA_Export_Job *)
      malloc(sizeof(Linear_Distribution_RSA_Export_Job));
  if (NULL == export_job) {
    critical("main_server(): Failed to allocate memory.");
  }

  export_job->distribution = distribution;
  safe_strlcpy(export_job->path, entry->path, MAX_SIZE_PATH_BUFFER);

  /* Spawn a thread to sort, export and deallocate the distribution. */
  thread_pool_spawn(pool, main_server_export_distribution, export_job);
}

/*!
 * \brief   The main function on the client node.
 *
 * This function is called once by main() for each distribution to generate.
 */
void main_client()
{
  Parameters parameters;
  parameters_init(&parameters);

  /* Receive broadcast of the distribution parameters. */
  parameters_bcast_recv(&parameters, MPI_RANK_ROOT);

  /* Receive broadcast of the dimension. */
  uint32_t dimension;
  if (MPI_SUCCESS !=
    MPI_Bcast(
      &dimension,
      1, /* count */
      MPI_UNSIGNED,
      MPI_RANK_ROOT,
      MPI_COMM_WORLD))
  {
    critical("main_client(): Failed to receive broadcast of the dimension.");
  }

  while (TRUE) {
    /* Notify the server that we are ready to receive a job. */
    uint32_t notification = MPI_NOTIFY_READY;

    if (MPI_SUCCESS != MPI_Send(
        &notification,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_TAG_NOTIFY,
        MPI_COMM_WORLD))
    {
      critical("main_client(): Failed to send notification.");
    }

    /* Receive the response. */
    uint32_t job;

    MPI_Status status;

    if (MPI_SUCCESS != MPI_Recv(
        &job,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_TAG_JOB,
        MPI_COMM_WORLD,
        &status))
    {
      critical("main_client(): Failed to receive job.");
    }

    if (MPI_JOB_STOP == job) {
      printf("Stopping.\n");
      break;
    } else if (MPI_JOB_PROCESS_SLICE != job) {
      critical("main_client(): Unknown job (%u).", job);
    }

    /* Receive the slice region. */
    uint32_t min_log_alpha;

    if (MPI_SUCCESS != MPI_Recv(
        &min_log_alpha,
        1, /* count */
        MPI_INT,
        MPI_RANK_ROOT,
        MPI_TAG_SLICE_ALPHAS,
        MPI_COMM_WORLD,
        &status))
    {
      critical("main_client(): Failed to receive MPI_TAG_SLICE_ALPHAS.");
    }

    const int32_t max_alpha = abs_i(min_log_alpha);

    const int32_t m = (int32_t)(parameters.m);

    /* Note that it is technically only necessary to skip in the two-dimensional 
     * case, as the error otherwise grows to great. We could remove the below 
     * restriction with a small gain in probability mass, and obtain a more 
     * symmetric distribution. */

    if (max_alpha > m + 10) {
      /* Notify the server that we are skipping this slice. */
      notification = MPI_NOTIFY_SLICE_SKIP;

      if (MPI_SUCCESS != MPI_Send(
          &notification,
          1, /* count */
          MPI_UNSIGNED,
          MPI_RANK_ROOT,
          MPI_TAG_NOTIFY,
          MPI_COMM_WORLD))
      {
        critical("main_client(): Failed to send notification.");
      }

      continue;
    }

    /* Compute slice. */
    Linear_Distribution_Slice slice;
    linear_distribution_slice_init(&slice, dimension);

    linear_distribution_slice_compute_richardson(
        &slice,
        &parameters,
        LINEAR_DISTRIBUTION_SLICE_COMPUTE_TARGET_D,
        min_log_alpha);

    /* Notify the server we are done computing the slice. */
    notification = MPI_NOTIFY_SLICE_DONE;

    if (MPI_SUCCESS != MPI_Send(
        &notification,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_TAG_NOTIFY,
        MPI_COMM_WORLD))
    {
      critical("main_client(): Failed to send notification.");
    }

    /* Send back the slice. */
    linear_distribution_slice_send(&slice, MPI_RANK_ROOT);

    /* Clear memory. */
    linear_distribution_slice_clear(&slice);
  }

  /* Clear memory. */
  parameters_clear(&parameters);
}

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  fprintf(file, "Synopsis: mpirun generate_linear_distribution_rsa \\\n"
          "   [ -delta <delta> ] [ -dim <dimension> ] \\\n"
          "      [ -max | -rnd | -exp <p> <q> ] <n> { <n> }\n");

  fprintf(file, "\n");
  fprintf(file, "Selection method for d : -- defaults to -max\n");

  fprintf(file, " -max  set d to 2^m - 1\n");
  fprintf(file,
    " -rnd  select d from random p and q where N = pq\n");
  fprintf(file,
      " -exp  explicitly set d given <p> and <q>\n");

  fprintf(file, "\n");
  fprintf(file,
    "Entries <n>: -- one distribution is generated for each entry\n");
  fprintf(file,
    " <n>   the even length n in bits of the RSA integer N = pq\n");

  fprintf(file, "\n");
  fprintf(file, "Delta: -- defaults to 20\n");
  fprintf(file,
    " -delta explicitly set delta to <delta>\n");

  fprintf(file, "\n");
  fprintf(file, "Dimension: -- defaults to 2048\n");
  fprintf(file,
    " -dim  explicitly set the slice dimension to <dimension>\n");

  fprintf(file,
    "\nNote: This implementation is optimized for small to medium kappa. If ");
  fprintf(file,
    "\nkappa is very large, it may fail to capture the probability mass.\n");
}

/*!
 * \brief The main entry point to the generate_linear_distribution_rsa
 *        executable.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \return Zero upon successful execution, non-zero otherwise.
 */
int main(int argc, char ** argv) {
  mpfr_set_default_prec(PRECISION);

  int provided;

  if (MPI_SUCCESS != MPI_Init_thread(
                        &argc,
                        &argv,
                        MPI_THREAD_FUNNELED, /* requested level of support */
                        &provided))
  {
    critical("main(): Failed to initialize MPI.");
  }

  int mpi_rank;
  int mpi_size;

  if (MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank)) {
    critical("main(): Failed to get the MPI world communicator rank.");
  }

  if (MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD, &mpi_size)) {
    critical("main(): Failed to get the MPI world communicator size.");
  }

  if (1 == argc) {
    if (MPI_RANK_ROOT == mpi_rank) {
      print_synopsis(stdout);
    }

    /* Exit gracefully. */
    MPI_Finalize();
    return 0;
  }

  if (mpi_size < 2) {
    if (MPI_RANK_ROOT == mpi_rank) {
      fprintf(stderr, "Incorrect environment: "
        "At least two MPI processors are required.\n");
    }

    /* Exit gracefully. */
    MPI_Finalize();
    return 0;
  }

  /* Parse the command line arguments and signal errors in a graceful way. */
  uint32_t result;

  Generate_Linear_Distribution_RSA_Arguments arguments;

  if (MPI_RANK_ROOT == mpi_rank) {
    result =
      (uint32_t)arguments_init_parse_command_line(&arguments, argc, argv);
  }

  if (MPI_SUCCESS != 
    MPI_Bcast(&result, 1, MPI_UNSIGNED, MPI_RANK_ROOT, MPI_COMM_WORLD))
  {
    critical("main(): "
      "Failed to broadcast the result of parsing the arguments.");
  };

  if (!result) {
    /* Exit gracefully if an error occurred. */
    MPI_Finalize();
    return 0;
  }

  /* Begin the generation process. */
  if (MPI_RANK_ROOT == mpi_rank) {
    printf("Number of MPI processors: %u\n", mpi_size);

    /* Report back information on the level of MPI thread support provided. */
    switch (provided) {
      case MPI_THREAD_SINGLE:
        printf("Level of thread support provided: "
          "MPI_THREAD_SINGLE\n");
        break;

      case MPI_THREAD_FUNNELED:
        printf("Level of thread support provided: "
          "MPI_THREAD_FUNNELED\n");
        break;

      case MPI_THREAD_SERIALIZED:
        printf("Level of thread support provided: "
          "MPI_THREAD_SERIALIZED\n");
        break;

      case MPI_THREAD_MULTIPLE:
        printf("Level of thread support provided: "
          "MPI_THREAD_MULTIPLE\n");
        break;

      default:
        printf("Level of thread support provided: "
          "Unknown (%d)\n", provided);
        break;
    }

    printf("\n");

    /* Broadcast the number of distributions. */
    if (MPI_SUCCESS !=
      MPI_Bcast(
        &(arguments.count),
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
    {
      critical("main(): Failed to broadcast the number of distributions.");
    }

    /* Process each distribution. */
    Thread_Pool pool;
    thread_pool_init(&pool);

    for (uint32_t i = 0; i < arguments.count; i++) {
      main_server(&arguments, &(arguments.entries[i]), mpi_size, &pool);
    }

    printf("Waiting for all export operations to finish...\n");

    thread_pool_join(&pool);
    thread_pool_clear(&pool);

    arguments_clear(&arguments);
  } else {
    uint32_t count;

    if (MPI_SUCCESS !=
      MPI_Bcast(
        &count,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
    {
      critical("main(): Failed to broadcast the number of distributions.");
    }

    for (uint32_t i = 0; i < count; i++) {
      main_client();
    }
  }

  if (MPI_SUCCESS != MPI_Finalize()) {
    critical("main(): Failed to finalize MPI.");
  }

  /* Signal success. */
  return 0;
}
