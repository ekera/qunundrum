/*!
 * \file    main_generate_linear_distribution.cpp
 * \ingroup generate_linear_distribution_exe
 *
 * \brief   The definition of the main entry point to the
 *          generate_linear_distribution executable, and of associated
 *          functions.
 */

/*!
 * \defgroup generate_linear_distribution_exe \
 *           The generate_linear_distribution executable
 * \ingroup  generate_executable
 *
 * \brief    A module for the generate_linear_distribution executable.
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
#include "random.h"
#include "string_utilities.h"

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
 *          parsed \<m\> \<s\> or \<m\> \<l\> tuples from the command line
 *          arguments.
 *
 * \ingroup generate_linear_distribution_exe
 */
typedef struct {
  /*!
   * \brief   The parameter m.
   */
  uint32_t m;

  /*!
   * \brief   The parameter s.
   */
  uint32_t s;

  /*!
   * \brief   The parameter l.
   */
  uint32_t l;

  /*!
   * \brief   The value of d or r for the combination of distribution type and
   *          selection method.
   *
   * The selection method and distribution type are common for all argument
   * tuples, see Generate_Linear_Distribution_::distribution_type and
   * Generate_Linear_Distribution_Arguments::selection_method.
   */
  mpz_t value;

  /*!
   * \brief   The path to which to export this distribution.
   */
  char path[MAX_SIZE_PATH_BUFFER];
} Generate_Linear_Distribution_Arguments_Entry;

/*!
 * \brief   A data structure representing parsed command line arguments.
 *
 * \ingroup generate_linear_distribution_exe
 */
typedef struct {
  /*!
   * \brief   The distribution type.
   */
  Distribution_Type distribution_type;

  /*!
   * \brief   The selection method for the logarithm d or order r.
   */
  Selection_Method selection_method;

  /*!
   * \brief   The tradeoff method.
   */
  Tradeoff_Method tradeoff_method;

  /*!
   * \brief   The explicitly specified value of the logarithm d or order r.
   *
   * If the value is not explicitly specified, this entry is set to zero.
   */
  mpz_t explicit_value;

  /*!
   * \brief   The dimension of the slices in the distribution.
   */
  uint32_t dimension;

  /*!
   * \brief   The number of \<m\> \<s\> or \<m\> \<l\> tuples passed to the
   *          executable.
   *
   * This number also corresponds to the number of distributions to generate.
   */
  uint32_t count;

  /*!
   * \brief   A vector of entries in the form of parsed tuples.
   *
   * Each entry corresponds to a distribution to generate.
   */
  Generate_Linear_Distribution_Arguments_Entry * entries;
} Generate_Linear_Distribution_Arguments;

/*!
 * \brief   A data structure representing an export job.
 *
 * \ingroup generate_linear_distribution_exe
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
} Linear_Distribution_Export_Job;

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
  Generate_Linear_Distribution_Arguments * const arguments,
  const int argc,
  char ** argv)
{
  /* Initialize the arguments data structure. */
  arguments->distribution_type = DISTRIBUTION_TYPE_UNDEFINED;
  arguments->selection_method = SELECTION_METHOD_UNDEFINED;
  arguments->tradeoff_method = TRADEOFF_METHOD_UNDEFINED;
  arguments->dimension = 0;
  arguments->entries = NULL;
  arguments->count = 0;

  mpz_init(arguments->explicit_value);
  mpz_set_ui(arguments->explicit_value, 0);

  /* Iterate over the command line arguments. */
  int i;

  for (i = 1; i < argc; i++) {
    /* Parse the selection method. */
    Selection_Method selection_method = SELECTION_METHOD_UNDEFINED;

    if (0 == strcmp(argv[i], "-det")) {
      selection_method = SELECTION_METHOD_DETERMINISTIC;
    } else if (0 == strcmp(argv[i], "-min")) {
      selection_method = SELECTION_METHOD_MINIMAL;
    } else if (0 == strcmp(argv[i], "-max")) {
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
        if ((i + 1) >= argc) {
          fprintf(stderr, "Error: Expected <value> to follow after -exp.\n");
          return FALSE;
        }

        if (0 != mpz_set_str(arguments->explicit_value, argv[i + 1], 10)) {
          fprintf(stderr, "Error: Failed to parse value after -exp as an "
            "integer.\n");
          return FALSE;
        }

        if (mpz_cmp_ui(arguments->explicit_value, 0) <= 0) {
          fprintf(stderr, "Error: The value after -exp must be positive.\n");
          return FALSE;
        }

        i++;
      }

      continue;
    }

    /* Parse the tradeoff method. */
    Tradeoff_Method tradeoff_method = TRADEOFF_METHOD_UNDEFINED;

    if (0 == strcmp(argv[i], "-s")) {
      tradeoff_method = TRADEOFF_METHOD_FACTOR;
    } else if (0 == strcmp(argv[i], "-l")) {
      tradeoff_method = TRADEOFF_METHOD_EXPLICIT;
    }

    if (TRADEOFF_METHOD_UNDEFINED != tradeoff_method) {
      /* Check that a tradeoff method has not already been specified. */
      if (TRADEOFF_METHOD_UNDEFINED != (arguments->tradeoff_method)) {
        fprintf(stderr, "Error: The tradeoff method cannot be twice "
          "specified.\n");
        return FALSE;
      }

      /* Store the distribution type. */
      arguments->tradeoff_method = tradeoff_method;

      continue;
    }


    /* Parse the distribution type. */
    Distribution_Type distribution_type = DISTRIBUTION_TYPE_UNDEFINED;

    if (0 == strcmp(argv[i], "-d")) {
      distribution_type = DISTRIBUTION_TYPE_LOGARITHM;
    } else if (0 == strcmp(argv[i], "-r")) {
      distribution_type = DISTRIBUTION_TYPE_ORDER;
    }

    if (DISTRIBUTION_TYPE_UNDEFINED != distribution_type) {
      /* Check that a distribution type has not already been specified. */
      if (DISTRIBUTION_TYPE_UNDEFINED != (arguments->distribution_type)) {
        fprintf(stderr, "Error: The distribution type cannot be twice "
          "specified.\n");
        return FALSE;
      }

      /* Store the distribution type. */
      arguments->distribution_type = distribution_type;

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
  if (0 == arguments->dimension) {
    arguments->dimension = 2048;
  }

  if (SELECTION_METHOD_UNDEFINED == arguments->selection_method) {
    arguments->selection_method = SELECTION_METHOD_MAXIMAL;
  }

  if (DISTRIBUTION_TYPE_UNDEFINED == arguments->distribution_type) {
    arguments->distribution_type = DISTRIBUTION_TYPE_LOGARITHM;
  }

  if (TRADEOFF_METHOD_UNDEFINED == arguments->tradeoff_method) {
    arguments->tradeoff_method = TRADEOFF_METHOD_FACTOR;
  }

  /* Parse tuples {<m> <s>} or {<m> <l>}. */
  if (((argc - i) <= 0) || (0 != ((argc - i) % 2))) {
    fprintf(stderr, "Error: Incorrect command line arguments; expected tuples "
      "but found an odd number of arguments.\n");
    return FALSE;
  }

  arguments->count = (uint32_t)((argc - i) / 2);

  arguments->entries =
    (Generate_Linear_Distribution_Arguments_Entry *)malloc(
      (arguments->count) * sizeof(
        Generate_Linear_Distribution_Arguments_Entry));
  if (NULL == arguments->entries) {
    fprintf(stderr, "Error: Failed to allocate memory for argument entries.\n");
    return FALSE;
  }

  for (uint32_t j = 0; j < arguments->count; j++) {
    mpz_init(arguments->entries[j].value);
    mpz_set_ui(arguments->entries[j].value, 0);
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
  for (uint32_t j = 0; j < arguments->count; j++, i += 2) {
    int m = atoi(argv[i]);
    if (m <= 1) {
      fprintf(stderr, "Error: Failed to parse <m>.\n");
      return FALSE;
    }

    arguments->entries[j].m = (uint32_t)m;

    if (TRADEOFF_METHOD_FACTOR == arguments->tradeoff_method) {
      int s = atoi(argv[i + 1]);
      if (s <= 0) {
        fprintf(stderr, "Error: Failed to parse <s>.\n");
        return FALSE;
      }

      arguments->entries[j].s = (uint32_t)s;
      arguments->entries[j].l = ceil((uint32_t)m / (uint32_t)s);
    } else {
      int l = atoi(argv[i + 1]);
      if (l <= 0) {
        fprintf(stderr, "Error: Failed to parse <l>.\n");
        return FALSE;
      }

      arguments->entries[j].s = 0;
      arguments->entries[j].l = (uint32_t)l;
    }

    char suffix[MAX_SIZE_PATH_BUFFER];
    if (TRADEOFF_METHOD_FACTOR == arguments->tradeoff_method) {
      safe_snprintf(suffix,
        MAX_SIZE_PATH_BUFFER,
        "dim-%u-%c-m-%u-s-%u.txt",
        arguments->dimension,
        (DISTRIBUTION_TYPE_ORDER == arguments->distribution_type) ? 'r' : 'd',
        arguments->entries[j].m,
        arguments->entries[j].s);
    } else {
      safe_snprintf(suffix,
        MAX_SIZE_PATH_BUFFER,
        "dim-%u-%c-m-%u-l-%u.txt",
        arguments->dimension,
        (DISTRIBUTION_TYPE_ORDER == arguments->distribution_type) ? 'r' : 'd',
        arguments->entries[j].m,
        arguments->entries[j].l);
    }

    /* Set the value. */
    if (SELECTION_METHOD_MAXIMAL == arguments->selection_method) {
      mpz_setbit(arguments->entries[j].value, arguments->entries[j].m);
      mpz_sub_ui(arguments->entries[j].value, arguments->entries[j].value, 1);

      safe_snprintf(arguments->entries[j].path, MAX_SIZE_PATH_BUFFER,
        "%s/linear-distribution-max-%s", DISTRIBUTIONS_DIRECTORY, suffix);
    } else if (SELECTION_METHOD_MINIMAL == arguments->selection_method) {
      mpz_setbit(arguments->entries[j].value, arguments->entries[j].m - 1);
      mpz_add_ui(arguments->entries[j].value, arguments->entries[j].value, 1);

      safe_snprintf(arguments->entries[j].path, MAX_SIZE_PATH_BUFFER,
        "%s/linear-distribution-min-%s", DISTRIBUTIONS_DIRECTORY, suffix);
    } else if (SELECTION_METHOD_DETERMINISTIC == arguments->selection_method) {
      if (m > 8192) {
        fprintf(stderr,
          "Error: When using the -det flag it is required that m <= 8192.\n");
        return FALSE;
      }

      {
        mpz_t r;
        mpz_init(r);

        mpz_t d;
        mpz_init(d);

        parameters_selection_deterministic_d_r(d, r, (uint32_t)m);

        if (DISTRIBUTION_TYPE_ORDER == arguments->distribution_type) {
          mpz_set(arguments->entries[j].value, r);
        } else {
          mpz_set(arguments->entries[j].value, d);
        }

        /* Clear memory. */
        mpz_clear(d);
        mpz_clear(r);
      }

      safe_snprintf(arguments->entries[j].path, MAX_SIZE_PATH_BUFFER,
        "%s/linear-distribution-det-%s", DISTRIBUTIONS_DIRECTORY, suffix);
    } else if (SELECTION_METHOD_RANDOM == arguments->selection_method) {
      if (j > 0) { /* Check if we should copy the previous value */
        if ((arguments->entries[j-1].m) == (arguments->entries[j].m)) {
          mpz_set(arguments->entries[j].value, arguments->entries[j-1].value);
        }
      }

      /* If we did not copy the previous value, select the value at random. */
      if (0 == mpz_cmp_ui(arguments->entries[j].value, 0)) {
        parameters_selection_random_d_or_r(
          arguments->entries[j].value, arguments->entries[j].m);
      }

      safe_snprintf(arguments->entries[j].path, MAX_SIZE_PATH_BUFFER,
        "%s/linear-distribution-rnd-%s", DISTRIBUTIONS_DIRECTORY, suffix);
    } else if (SELECTION_METHOD_EXPLICIT == arguments->selection_method) {
      /* Check that the explicitly specified value is less than 2^m. */
      mpz_setbit(arguments->entries[j].value, arguments->entries[j].m);

      if (mpz_cmp(arguments->explicit_value, arguments->entries[j].value) >= 0)
      {
        fprintf(stderr, "Error: The explicitly specified value with -exp is "
          "greater than or equal to 2^m for at least one <m>.\n");
        return FALSE;
      }

      if (DISTRIBUTION_TYPE_ORDER == arguments->distribution_type) {
        /* Check that the specified value is greater than 2^(m-1). */
        mpz_set_ui(arguments->entries[j].value, 0);
        mpz_setbit(arguments->entries[j].value, arguments->entries[j].m - 1);

        if (mpz_cmp(arguments->explicit_value, arguments->entries[j].value) 
          <= 0)
        {
          fprintf(stderr, "Error: The explicitly specified value with -exp is "
            "less than or equal to 2^(m-1) for at least one <m>.\n");
          return FALSE;
        }
      }

      /* Note: It has already been checked that 0 < value when parsing the
       *       command line arguments, so we need not do it again. */

      mpz_set(arguments->entries[j].value, arguments->explicit_value);

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
  Generate_Linear_Distribution_Arguments * const arguments)
{
  mpz_clear(arguments->explicit_value);

  if (NULL != arguments->entries) {
    for (uint32_t i = 0; i < arguments->count; i++) {
      mpz_clear(arguments->entries[i].value);
    }

    free(arguments->entries);
    arguments->entries = NULL;
  }

  memset(arguments, 0, sizeof(Generate_Linear_Distribution_Arguments));
}

/*!
 * \brief   Exports the computed probability distribution to file.
 *
 * This function is intended to run in a separate thread.
 *
 * \param[in, out] ptr    A pointer to a Linear_Distribution_Export_Job data
 *                        structure that describes the export job.
 *
 * \remark   The export job and the distribution it contains will be cleared
 *           and deallocated by this function.
 */
static void * main_server_export_distribution(void * ptr)
{
  /* Map the input to data structures. */
  Linear_Distribution_Export_Job * job = (Linear_Distribution_Export_Job *)ptr;

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

  free(job);
  job = NULL;
  ptr = NULL;

  /* Finished. */
  return NULL;
}

/*!
 * \brief   The main function on the server node.
 *
 * This function is called once by main() for each distribution to generate.
 *
 * \param[in] arguments   The parsed command line arguments.
 * \param[in] entry       The \<m\> \<s\> or \<m\> \<l\> tuple to process.
 * \param[in] mpi_size    The number of nodes.
 * \param[in] pool        The thread pool from which to spawn threads for
 *                        export operations that run in the background.
 */
static void main_server(
  const Generate_Linear_Distribution_Arguments * const arguments,
  const Generate_Linear_Distribution_Arguments_Entry * const entry,
  int mpi_size,
  Thread_Pool * const pool)
{
  /* Setup the parameters. */
  const uint32_t t = 30;

  Parameters parameters;
  parameters_init(&parameters);

  if (TRADEOFF_METHOD_FACTOR == arguments->tradeoff_method) {
    parameters_explicit_m_s(
      &parameters,
      entry->value, /* d */
      entry->value, /* r */
      entry->m,
      entry->s,
      t);
  } else {
    parameters_explicit_m_l(
      &parameters,
      entry->value, /* d */
      entry->value, /* r */
      entry->m,
      entry->l,
      t);
  }

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

  /* Send broadcast of the distribution type. */
  uint32_t tmp = (uint32_t)(arguments->distribution_type);

  if (MPI_SUCCESS !=
    MPI_Bcast(
      &tmp,
      1, /* count */
      MPI_UNSIGNED,
      MPI_RANK_ROOT,
      MPI_COMM_WORLD))
  {
    critical("main_server(): "
      "Failed to send broadcast of the distribution type.");
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

  if (DISTRIBUTION_TYPE_LOGARITHM == arguments->distribution_type) {
    flags |= LINEAR_DISTRIBUTION_FLAG_D;
  } else if (DISTRIBUTION_TYPE_ORDER == arguments->distribution_type) {
    flags |= LINEAR_DISTRIBUTION_FLAG_R;
  } else {
    critical("main_server(): Unknown distribution type.");
  }

  if (SELECTION_METHOD_DETERMINISTIC == arguments->selection_method) {
    flags |= LINEAR_DISTRIBUTION_FLAG_DETERMINISTIC;
  } else if (SELECTION_METHOD_EXPLICIT == arguments->selection_method) {
    flags |= LINEAR_DISTRIBUTION_FLAG_EXPLICIT;
  } else if (SELECTION_METHOD_RANDOM == arguments->selection_method) {
    flags |= LINEAR_DISTRIBUTION_FLAG_RANDOM;
  } else if (SELECTION_METHOD_MINIMAL == arguments->selection_method) {
    flags |= LINEAR_DISTRIBUTION_FLAG_MINIMAL;
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

        printf("Stopping node %d with %u node(s) remaining...\n",
          status.MPI_SOURCE,
          mpi_nodes_remaining);
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
       *       the call to linear_distribution_dealloc(). It must not be
       *       deallocated here as doing so would lead to memory corruption. */
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
  Linear_Distribution_Export_Job * export_job =
    (Linear_Distribution_Export_Job *)
      malloc(sizeof(Linear_Distribution_Export_Job));
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
static void main_client()
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

  /* Receive broadcast of the distribution type. */
  uint32_t tmp;

  if (MPI_SUCCESS !=
    MPI_Bcast(
      &tmp,
      1, /* count */
      MPI_UNSIGNED,
      MPI_RANK_ROOT,
      MPI_COMM_WORLD))
  {
    critical("main_client(): "
      "Failed to receive broadcast of the distribution type.");
  }

  const Distribution_Type distribution_type = (Distribution_Type)tmp;

  /* Process jobs. */
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
      break;
    } else if (MPI_JOB_PROCESS_SLICE != job) {
      critical("main_client(): Unknown job (%u).", job);
    }

    /* Receive the slice region. */
    int32_t min_log_alpha;

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

    /* Skip this slice if the region is out of bounds. */
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

    Linear_Distribution_Slice_Compute_Target target = 
      LINEAR_DISTRIBUTION_SLICE_COMPUTE_TARGET_D;

    if (DISTRIBUTION_TYPE_ORDER == distribution_type) {
      target = LINEAR_DISTRIBUTION_SLICE_COMPUTE_TARGET_R;
    } else if (DISTRIBUTION_TYPE_LOGARITHM != distribution_type) {
      critical("main_client(): Incorrect distribution type.");
    }

    linear_distribution_slice_compute_richardson(
        &slice,
        &parameters,
        target,
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
  fprintf(file, "Synopsis: mpirun generate_linear_distribution \\\n"
          "   [ -d | -r ] [ -dim <dimension> ] \\\n"
          "      [ -min | -max | -det | -rnd | -exp <value> ] \\\n"
          "         ( <m> <s> { <m> <s> } | -l <m> <l> { <m> <l> } )\n");

  fprintf(file, "\n");
  fprintf(file, "Distribution type: -- defaults to -d\n");
  fprintf(file,
    " -d    generate a distribution for a short discrete logarithm d\n");
  fprintf(file,
    " -r    generate a distribution for an order r\n");

  fprintf(file, "\n");
  fprintf(file, "Selection method for d or r: -- defaults to -max\n");

  fprintf(file, " -min  set d or r to 2^(m-1) + 1\n");
  fprintf(file, " -max  set d or r to 2^m - 1\n");
  fprintf(file,
    " -det  select d or r deterministically from Catalan's constant\n");
  fprintf(file,
    " -rnd  select d or r uniformly at random from (2^(m-1), 2^m)\n");
  fprintf(file,
      " -exp  explicitly set d or r to <value>\n");

  fprintf(file, "\n");
  fprintf(file,
    "Tuples <m> <s> or <m> <l>: -- one distribution is generated for each "
      "tuple\n");
  fprintf(file,
    " <m>   the length in bits of d or r\n");
  fprintf(file,
    " <s>   the tradeoff factor s; used to set l = ceil(m / s)\n");
  fprintf(file,
    " <l>   the length in bits of the prefix\n");

  fprintf(file, "\n");
  fprintf(file, "Dimension: -- defaults to 2048\n");
  fprintf(file,
    " -dim  explicitly set the slice dimension to <dimension>\n");

  fprintf(file,
    "\nNote: If -rnd is specified and m is kept constant for consecutive "
      "tuples\n");
  fprintf(file,
    "<m> <s> or <m> <l>, the same value of d or r will be re-used.\n");

  fprintf(file,
    "\nNote: This implementation is optimized for small to medium kappa. If ");
  fprintf(file,
    "\nkappa is very large, it may fail to capture the probability mass.\n");
}

/*!
 * \brief   The main entry point to the generate_linear_distribution executable.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \return  Zero upon successful execution, non-zero otherwise.
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

  Generate_Linear_Distribution_Arguments arguments;

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
