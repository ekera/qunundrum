/*!
 * \file    main_solve_diagonal_distribution.cpp
 * \ingroup solve_diagonal_distribution_exe
 *
 * \brief   The definition of the main entry point to the
 *          solve_diagonal_distribution executable, and of associated functions.
 */

/*!
 * \defgroup solve_diagonal_distribution_exe \
 *           The solve_diagonal_distribution executable
 * \ingroup  solve_executable
 *
 * \brief    A module for the solve_diagonal_distribution executable.
 */

#include "executables.h"
#include "executables_solve_distribution.h"

#include "common.h"
#include "diagonal_distribution.h"
#include "diagonal_distribution_loader.h"
#include "errors.h"
#include "lattice.h"
#include "lattice_enumerate.h"
#include "lattice_solve.h"
#include "log.h"
#include "math.h"
#include "random.h"
#include "sample.h"
#include "string_utilities.h"
#include "timer.h"

#include <gmp.h>
#include <mpfr.h>

#include <mpi.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include <sys/stat.h>

/*!
 * \brief   A data structure representing argument entries in the form of
 *          parsed \<distribution\> \<n\> tuples from the command line
 *          arguments.
 *
 * \ingroup solve_diagonal_distribution_exe
 */
typedef struct {
  /*!
   * \brief   The parameter n.
   */
  uint32_t n;

  /*!
   * \brief   The path to the distribution.
   */
  char * path;
} Solve_Diagonal_Distribution_Arguments_Entry;

/*!
 * \brief   A data structure representing parsed command line arguments.
 *
 * \ingroup solve_diagonal_distribution_exe
 */
typedef struct {
  /*!
   * \brief   The bound on the offset Delta from k_0 when sampling k given j.
   */
  uint32_t delta_bound;

  /*!
   * \brief   The bound on the peak index eta.
   */
  uint32_t eta_bound;

  /*!
   * \brief   A boolean flag that indicates whether the bound on the peak index
   *          eta has been specified. Otherwise, the eta bound used to setup
   *          each distribution is used.
   */
  bool eta_bound_specified;

  /*!
   * \brief   The search strategy.
   */
  Search_Strategy search_strategy;

  /*!
   * \brief   The solution method.
   */
  Solution_Method solution_method;

  /*!
   * \brief   The reduction algorithm.
   */
  Lattice_Reduction_Algorithm reduction_algorithm;

  /*!
   * \brief   The enumeration timeout.
   */
  uint32_t timeout;

  /*!
   * \brief   The number of distributions to process.
   */
  uint32_t count;

  /*!
   * \brief   Entries for each distribution to process.
   */
  Solve_Diagonal_Distribution_Arguments_Entry * entries;
} Solve_Diagonal_Distribution_Arguments;


/*!
 * \name Command line arguments
 * \{
 */

/*!
 * \brief   Parses the command line arguments.
 *
 * \param[in, out] arguments   The arguments data structure in which to store
 *                             the parsed command line arguments.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \remark   So as to allow parallelized executables to be shut down gracefully,
 *           this function does not call critical() to signal a critical error.
 *           Rather it prints an informative error message to stderr and returns
 *           #FALSE. When returning #FALSE, this function expects the caller to
 *           terminate the executable.
 *
 * \return   #TRUE if the command line arguments were successfully parsed,
 *           #FALSE otherwise. If #FALSE is returned, the data structure may be
 *           only partially initialized and memory only partially allocated.
 */
static bool arguments_init_parse_command_line(
  Solve_Diagonal_Distribution_Arguments * const arguments,
  const int argc,
  char ** argv)
{
  /* Initialize the arguments data structure. */
  arguments->delta_bound = BOUND_DELTA;
  arguments->eta_bound = 0;
  arguments->eta_bound_specified = FALSE;
  arguments->search_strategy = SEARCH_STRATEGY_DEFAULT;
  arguments->solution_method = SOLUTION_METHOD_DEFAULT;
  arguments->reduction_algorithm = REDUCTION_ALGORITHM_DEFAULT;
  arguments->timeout = 0;
  arguments->entries = NULL;
  arguments->count = 0;

  bool delta_bound_specified = FALSE;
  bool timeout_specified = FALSE;

  /* Iterate over the command line arguments. */
  int i;

  for (i = 1; i < argc; i++) {
    /* Parse the delta bound. */
    if (0 == strcmp(argv[i], "-delta-bound")) {
      /* Check that a delta bound has not already been specified. */
      if (FALSE != delta_bound_specified) {
        fprintf(stderr, "Error: The delta bound cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr, "Error: Expected <delta-bound> to follow after "
          "-delta-bound.\n");
        return FALSE;
      }

      const int x = atoi(argv[i + 1]);
      if (x < 0) {
        fprintf(stderr, "Error: The <delta-bound> passed to -delta-bound must "
          "be non-negative.\n");
        return FALSE;
      }

      /* Store the delta bound. */
      arguments->delta_bound = (uint32_t)x;
      delta_bound_specified = TRUE;

      i++;

      continue;
    }

    /* Parse the eta bound. */
    if (0 == strcmp(argv[i], "-eta-bound")) {
      /* Check that an eta bound has not already been specified. */
      if (FALSE != arguments->eta_bound_specified) {
        fprintf(stderr, "Error: The eta bound cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr, "Error: Expected <eta-bound> to follow after "
          "-eta-bound.\n");
        return FALSE;
      }

      const int x = atoi(argv[i + 1]);
      if ((x < 0) || (x > 256)) {
        fprintf(stderr, "Error: The <eta-bound> passed to -eta-bound must be "
          "on the interval [0, 256].\n");
        return FALSE;
      }

      /* Store the eta bound. */
      arguments->eta_bound = (uint32_t)x;
      arguments->eta_bound_specified = TRUE;

      i++;

      continue;
    }

    /* Parse the search strategy. */
    Search_Strategy search_strategy = SEARCH_STRATEGY_DEFAULT;

    if (0 == strcmp(argv[i], "-adaptive")) {
      search_strategy = SEARCH_STRATEGY_ADAPTIVE;
    } else if (0 == strcmp(argv[i], "-non-adaptive")) {
      search_strategy = SEARCH_STRATEGY_NON_ADAPTIVE;
    } else if (0 == strcmp(argv[i], "-non-adaptive-early-abort")) {
      search_strategy = SEARCH_STRATEGY_NON_ADAPTIVE_EARLY_ABORT;
    }

    if (SEARCH_STRATEGY_DEFAULT != search_strategy) {
      /* Check that a search strategy has not already been specified. */
      if (SEARCH_STRATEGY_DEFAULT != (arguments->search_strategy)) {
        fprintf(stderr, "Error: The search strategy cannot be twice "
          "specified.\n");
        return FALSE;
      }

      /* Store the search strategy. */
      arguments->search_strategy = search_strategy;

      continue;
    }

    /* Parse the solution method. */
    Solution_Method solution_method = SOLUTION_METHOD_DEFAULT;

    if (0 == strcmp(argv[i], "-closest")) {
      solution_method = SOLUTION_METHOD_CLOSEST;
    } else if (0 == strcmp(argv[i], "-enumerate")) {
      solution_method = SOLUTION_METHOD_ENUMERATE;
    }

    if (SOLUTION_METHOD_DEFAULT != solution_method) {
      /* Check that a selection method has not already been specified. */
      if (SOLUTION_METHOD_DEFAULT != (arguments->solution_method)) {
        fprintf(stderr, "Error: The solution method cannot be twice "
          "specified.\n");
        return FALSE;
      }

      /* Store the solution method. */
      arguments->solution_method = solution_method;

      continue;
    }

    /* Parse the reduction algorithm. */
    Lattice_Reduction_Algorithm reduction_algorithm =
      REDUCTION_ALGORITHM_DEFAULT;

    if (0 == strcmp(argv[i], "-lll")) {
      reduction_algorithm = REDUCTION_ALGORITHM_LLL;
    } else if (0 == strcmp(argv[i], "-bkz")) {
      reduction_algorithm = REDUCTION_ALGORITHM_BKZ;
    } else if (0 == strcmp(argv[i], "-lll-then-bkz")) {
      reduction_algorithm = REDUCTION_ALGORITHM_LLL_BKZ;
    } else if (0 == strcmp(argv[i], "-hkz")) {
      reduction_algorithm = REDUCTION_ALGORITHM_HKZ;
    }

    if (REDUCTION_ALGORITHM_DEFAULT != reduction_algorithm) {
      /* Check that a reduction algorithm has not already been specified. */
      if (REDUCTION_ALGORITHM_DEFAULT != (arguments->reduction_algorithm)) {
        fprintf(stderr, "Error: The reduction algorithm cannot be twice "
          "specified.\n");
        return FALSE;
      }

      /* Store the reduction algorithm. */
      arguments->reduction_algorithm = reduction_algorithm;

      continue;
    }

    /* Parse the timeout. */
    if (0 == strcmp(argv[i], "-timeout")) {
      /* Check that a timeout has not already been specified. */
      if (timeout_specified) {
        fprintf(stderr, "Error: The timeout cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr, "Error: Expected value to follow after -timeout.\n");
        return FALSE;
      }

      const int x = atoi(argv[i+1]);
      if (x <= 0) {
        fprintf(stderr, "Error: The timeout must be positive.\n");
        return FALSE;
      }

      /* Store the timeout. */
      arguments->timeout = (uint32_t)x;
      timeout_specified = TRUE;

      i++;

      continue;
    }

    /* Stop parsing. */
    break;
  }

  /* Set default parameters if arguments were not explicitly specified. */
  if (SOLUTION_METHOD_DEFAULT == arguments->solution_method) {
    arguments->solution_method = SOLUTION_METHOD_CLOSEST;
  }

  if (REDUCTION_ALGORITHM_DEFAULT == arguments->reduction_algorithm) {
    arguments->reduction_algorithm = REDUCTION_ALGORITHM_LLL_BKZ;
  }

  if (SEARCH_STRATEGY_DEFAULT == arguments->search_strategy) {
    arguments->search_strategy = SEARCH_STRATEGY_ADAPTIVE;
  }

  if (TRUE != timeout_specified) {
    arguments->timeout = 300;
  }

  /* Parse tuples { <distribution> <n> }. */
  if (((argc - i) <= 0) || (0 != ((argc - i) % 2))) {
    fprintf(stderr, "Error: Incorrect command line arguments; expected tuples "
      "but found an odd number of arguments.\n");
    return FALSE;
  }

  arguments->count = (uint32_t)((argc - i) / 2);

  arguments->entries =
    (Solve_Diagonal_Distribution_Arguments_Entry *)malloc(
      (arguments->count) * sizeof(Solve_Diagonal_Distribution_Arguments_Entry));
  if (NULL == arguments->entries) {
    fprintf(stderr, "Error: Failed to allocate memory for argument entries.\n");
    return FALSE;
  }

  for (uint32_t j = 0; j < arguments->count; j++) {
    arguments->entries[j].path = NULL;
  }

  /* Iterate over the tuples. */
  for (uint32_t j = 0; j < arguments->count; j++, i += 2) {
    /* Parse the path. */
    const uint32_t length = (uint32_t)strlen(argv[i]);

    arguments->entries[j].path = (char *)malloc((length + 1) * sizeof(char));
    if (NULL == arguments->entries[j].path) {
      critical("arguments_init_parse_command_line(): "
        "Failed to allocate memory.");
    }

    safe_strlcpy(arguments->entries[j].path, argv[i], length + 1);

    /* Check that the distribution exists. */
    if (0 != access(arguments->entries[j].path, F_OK)) {
      fprintf(stderr, "Error: The distribution \"%s\" does not exist.\n",
        arguments->entries[j].path);
      return FALSE;
    }

    if (0 != access(arguments->entries[j].path, R_OK)) {
      fprintf(stderr, "Error: The distribution \"%s\" is not readable.\n",
        arguments->entries[j].path);
      return FALSE;
    }

    /* Parse n. */
    const int n = atoi(argv[i + 1]);
    if (n < 1) {
      fprintf(stderr, "Error: Failed to parse <n>.\n");
      return FALSE;
    }

    arguments->entries[j].n = (uint32_t)n;
  }

  /* Signal success. */
  return TRUE;
}

/*!
 * \brief   Broadcasts the command line arguments to all other nodes.
 *
 * \param[in] arguments   The parsed command line arguments to broadcast.
 * \param[in] root        The rank of the root node.
 */
static void arguments_bcast_send(
  const Solve_Diagonal_Distribution_Arguments * const arguments,
  const int root)
{
  int result;

  uint32_t data[8];
  data[0] = arguments->delta_bound;
  data[1] = arguments->eta_bound;
  data[2] = (TRUE == arguments->eta_bound_specified) ? 1 : 0;

  data[3] = (uint32_t)(arguments->search_strategy);
  data[4] = (uint32_t)(arguments->solution_method);
  data[5] = (uint32_t)(arguments->reduction_algorithm);
  data[6] = (uint32_t)(arguments->timeout);
  data[7] = (uint32_t)(arguments->count);

  result = MPI_Bcast(data, 8, MPI_UNSIGNED, root, MPI_COMM_WORLD);
  if (MPI_SUCCESS != result) {
    critical("arguments_bcast_send(): Failed to send broadcast.");
  };

  for (uint32_t i = 0; i < (arguments->count); i++) {
    result = MPI_Bcast(
                &(arguments->entries[i].n),
                1,
                MPI_UNSIGNED,
                root,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("arguments_bcast_send(): Failed to send broadcast of n.");
    };

    uint32_t length = (uint32_t)strlen(arguments->entries[i].path);

    result = MPI_Bcast(
                &length,
                1,
                MPI_UNSIGNED,
                root,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("arguments_bcast_send(): Failed to send broadcast of the path "
        "length.");
    };

    result = MPI_Bcast(
                arguments->entries[i].path,
                length + 1,
                MPI_CHAR,
                root,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("arguments_bcast_send(): Failed to send broadcast of the path.");
    };
  }
}

/*!
 * \brief   Initializes the command line arguments by receiving a broadcast from
 *          a node.
 *
 * \param[in, out] arguments  The command line arguments to initialize.
 * \param[in] root            The rank of the node from which to receive.
 */
static void arguments_init_bcast_recv(
  Solve_Diagonal_Distribution_Arguments * const arguments,
  const int root)
{
  arguments->entries = NULL;
  arguments->count = 0;

  int result;

  uint32_t data[8];

  result = MPI_Bcast(data, 8, MPI_UNSIGNED, root, MPI_COMM_WORLD);
  if (MPI_SUCCESS != result) {
    critical("arguments_init_bcast_recv(): "
      "Failed to receive broadcast of collected metadata.");
  };

  arguments->delta_bound = data[0];
  arguments->eta_bound = data[1];
  arguments->eta_bound_specified = (1 == data[2]) ? TRUE : FALSE;
  arguments->search_strategy = (Search_Strategy)(data[3]);
  arguments->solution_method = (Solution_Method)(data[4]);
  arguments->reduction_algorithm = (Lattice_Reduction_Algorithm)(data[5]);
  arguments->timeout = data[6];
  arguments->count = data[7];

  arguments->entries =
    (Solve_Diagonal_Distribution_Arguments_Entry *)malloc(
      (arguments->count) * sizeof(Solve_Diagonal_Distribution_Arguments_Entry));
  if (NULL == arguments->entries) {
    critical("arguments_init_bcast_recv(): Failed to allocate memory.");
  }

  for (uint32_t i = 0; i < arguments->count; i++) {
    arguments->entries[i].path = NULL;
  }

  for (uint32_t i = 0; i < arguments->count; i++) {
    result = MPI_Bcast(
                &(arguments->entries[i].n),
                1,
                MPI_UNSIGNED,
                root,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("arguments_init_bcast_recv(): "
        "Failed to receive broadcast of n.");
    };

    uint32_t length;

    result = MPI_Bcast(
                &length,
                1,
                MPI_UNSIGNED,
                root,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("arguments_init_bcast_recv(): "
        "Failed to receive broadcast of the path length.");
    };

    arguments->entries[i].path = (char *)malloc((length + 1) * sizeof(char));
    if (NULL == arguments->entries[i].path) {
      critical("arguments_init_bcast_recv(): Failed to allocate memory.");
    }

    result = MPI_Bcast(
                arguments->entries[i].path,
                length + 1,
                MPI_CHAR,
                root,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("arguments_init_bcast_recv(): "
        "Failed to receive broadcast of the path.");
    };
  }
}

/*!
 * \brief   Prints the command line arguments.
 *
 * \param[in, out] file     Then file to which to print the arguments.
 * \param[in] arguments     The parsed command line arguments to print.
 * \param[in] distribution  The distribution being processed.
 */
static void arguments_fprintf(
  FILE * const file,
  const Solve_Diagonal_Distribution_Arguments * const arguments,
  const Diagonal_Distribution * const distribution)
{
  if (TRUE == arguments->eta_bound_specified) {
    mpfr_fprintf(file, "# Bounds: (eta = %u (%u), delta = %u)\n",
      arguments->eta_bound,
      distribution->parameters.eta_bound,
      arguments->delta_bound);
  } else {
    mpfr_fprintf(file, "# Bounds: (eta = <all> (%u), delta = %u)\n",
      distribution->parameters.eta_bound,
      arguments->delta_bound);
  }

  switch (arguments->search_strategy) {
    case SEARCH_STRATEGY_DEFAULT:
    case SEARCH_STRATEGY_ADAPTIVE:
      fprintf(file, "# Search strategy: Adaptive\n");
      break;

    case SEARCH_STRATEGY_NON_ADAPTIVE:
      fprintf(file, "# Search strategy: Non-adaptive\n");
      break;

    case SEARCH_STRATEGY_NON_ADAPTIVE_EARLY_ABORT:
      fprintf(file, "# Search strategy: Non-adaptive with early abort\n");
      break;
  }

  switch (arguments->solution_method) {
    case SOLUTION_METHOD_DEFAULT:
    case SOLUTION_METHOD_CLOSEST:
      fprintf(file, "# Solution method: Closest\n");
      break;

    case SOLUTION_METHOD_ENUMERATE:
      fprintf(file, "# Solution method: Enumerate (timeout: %u s)\n",
        arguments->timeout);
      break;
  }

  switch (arguments->reduction_algorithm) {
    case REDUCTION_ALGORITHM_DEFAULT:
    case REDUCTION_ALGORITHM_LLL_BKZ:
      fprintf(file, "# Reduction algorithm: LLL then BKZ\n");
      break;

    case REDUCTION_ALGORITHM_LLL:
      fprintf(file, "# Reduction algorithm: LLL\n");
      break;

    case REDUCTION_ALGORITHM_BKZ:
      fprintf(file, "# Reduction algorithm: BKZ\n");
      break;

    case REDUCTION_ALGORITHM_HKZ:
      fprintf(file, "# Reduction algorithm: HKZ\n");
      break;
  }

  log_timestamp_fprintf(file);
}

/*!
 * \brief   Clears an initialized command line arguments data structure.
 *
 * \param[in, out] arguments   The arguments data structure to clear.
 */
static void arguments_clear(
  Solve_Diagonal_Distribution_Arguments * const arguments)
{
  if (NULL != arguments->entries) {
    for (uint32_t i = 0; i < arguments->count; i++) {
      if (NULL != arguments->entries[i].path) {
        free(arguments->entries[i].path);
        arguments->entries[i].path = NULL;
      }
    }

    free(arguments->entries);
    arguments->entries = NULL;
  }

  memset(arguments, 0, sizeof(Solve_Diagonal_Distribution_Arguments));
}

/*!
 * \}
 */

/*!
 * \brief   The main function on the server node.
 *
 * This function is called once by main() for each distribution to process.
 *
 * \param[in] arguments    The parsed command line arguments.
 * \param[in] entry        The arguments entry corresponding to the distribution
 *                         to use for sampling problem instances to solve.
 * \param[in] distribution The distribution to use for sampling problem
 *                         instances to solve.
 * \param[in] mpi_size     The number of nodes.
 */
static void main_server(
  const Solve_Diagonal_Distribution_Arguments * const arguments,
  const Solve_Diagonal_Distribution_Arguments_Entry * const entry,
  const Diagonal_Distribution * const distribution,
  const int mpi_size)
{
  /* Create the log directory if it does not exist. */
  if (0 != access(LOGS_DIRECTORY, F_OK)) {
    if (0 != mkdir(LOGS_DIRECTORY, DIRECTORY_PERMISSIONS)) {
      critical("main_server(): Failed to create the directory \"%s\".",
        LOGS_DIRECTORY);
    }
  }

  /* Open the log file. */
  char log_path[MAX_SIZE_PATH_BUFFER];
  safe_snprintf(
    log_path, MAX_SIZE_PATH_BUFFER,
    "%s/solve-diagonal.txt", LOGS_DIRECTORY);

  FILE * log_file = fopen(log_path, "a+");
  if (NULL == log_file) {
    critical("main_server(): Failed to open \"%s\" for appending.", log_path);
  }

  fprintf(log_file, "\n# Processing: %s\n", truncate_path(entry->path));
  arguments_fprintf(log_file, arguments, distribution);

  /* Broadcast the distribution. */
  diagonal_distribution_bcast_send(distribution, MPI_RANK_ROOT);

  /* Setup a solution status data structure. */
  Solution_Status solution_status;
  solution_status_init(
    &solution_status,
    distribution->parameters.m,
    distribution->parameters.sigma,
    distribution->parameters.s,
    distribution->parameters.l,
    entry->n,
    TRUE); /* has_sigma */

  /* The number of currently running client nodes. */
  int nodes = mpi_size - 1;

  /* Declare status. */
  MPI_Status status;

  /* Declare state machine state. */
  Adaptive_Search_State state = SEARCH_STATE_PIVOT;

  while (TRUE) {
     /* Listen for a notification. */
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

    /* If a solution is done, receive the solution, otherwise proceed. */
    if (MPI_NOTIFY_SAMPLING_FAILED_OUT_OF_BOUNDS == notification) {
      uint32_t tmp_n;

      if (MPI_SUCCESS != MPI_Recv(
          &tmp_n,
          1, /* count */
          MPI_UNSIGNED,
          status.MPI_SOURCE,
          MPI_TAG_SAMPLE_COUNT,
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE))
      {
        critical("main_server(): Failed to receive n.");
      }

      if (tmp_n == solution_status.n) {
        /* Update the solution status data structure. */
        solution_status.fail_count++;
        solution_status.fail_out_of_bounds_count++;
      } else {
        /* Discard the notification as it is not for this n. */
      }
    } else if (MPI_NOTIFY_SOLUTION_DONE == notification) {
      uint32_t solution[6];

      if (MPI_SUCCESS != MPI_Recv(
          solution,
          6, /* count */
          MPI_UNSIGNED,
          status.MPI_SOURCE,
          MPI_TAG_SOLUTION,
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE))
      {
        critical("main_server(): Failed to receive solution.");
      }

      /* Check if the solution corresponds to the value of n currently being
       * processed. If n was recently increased or decreased, it may otherwise
       * be that we receive a solution to a problem for an old value of n. If
       * so, we simply ignore the solution and send out a new job. */
      if (solution[0] == solution_status.n) {
        /* Update the solution status data structure. */

        /* Update the timing statistics. */
        uint64_t time_prepare_system_us;

        time_prepare_system_us   = (uint64_t)solution[2]; /* high 32 bits */
        time_prepare_system_us <<= 32;
        time_prepare_system_us  ^= (uint64_t)solution[3]; /* low 32 bits */

        timer_statistics_insert(
          &(solution_status.statistics_prepare_system),
          time_prepare_system_us);

        uint64_t time_solve_system_us;

        time_solve_system_us   = (uint64_t)solution[4]; /* high 32 bits */
        time_solve_system_us <<= 32;
        time_solve_system_us  ^= (uint64_t)solution[5]; /* low 32 bits */

        timer_statistics_insert(
          &(solution_status.statistics_solve_system),
          time_solve_system_us);

        /* Update the success and failure counts. */
        const Lattice_Status_Recovery recovery_status =
          (Lattice_Status_Recovery)solution[1];

        switch (recovery_status) {
          case LATTICE_STATUS_RECOVERED_IMMEDIATE:
          case LATTICE_STATUS_RECOVERED_IMMEDIATE_SMOOTH:
          case LATTICE_STATUS_RECOVERED_SEARCH:
          case LATTICE_STATUS_RECOVERED_ENUMERATION:
          case LATTICE_STATUS_RECOVERED_ENUMERATION_SMOOTH:
           solution_status.success_count++;
            break;

          default:
            solution_status.fail_count++;
            break;
        }
      } else {
        /* Discard the notification as it is not for this n. */
      }
    } else if (MPI_NOTIFY_READY != notification) {
      critical("main_server(): Received unknown notification.");
    }

    /* Report status. */
    if (MPI_NOTIFY_READY != notification) {
      solution_status_print(stdout, &solution_status, NULL);
    }

    /* Check if a state change is in order. */
    if (solution_status.fail_count > 10) {
      if (SEARCH_STRATEGY_NON_ADAPTIVE_EARLY_ABORT ==
        arguments->search_strategy)
      {
        printf("** stopping\n");
        solution_status_print(log_file, &solution_status, "** stopping");
        break;
      }

      if (SEARCH_STATE_PIVOT == state) {
        if (SEARCH_STRATEGY_ADAPTIVE == arguments->search_strategy) {
          state = SEARCH_STATE_INCREASE;
        }
      }

      if (SEARCH_STATE_INCREASE == state) {
        printf("** incrementing n\n");
        solution_status_print(log_file, &solution_status, "** incrementing");
        solution_status_reset(&solution_status);

        solution_status.n++;
      } else if (SEARCH_STATE_DECREASE == state) {
        /* We have searched for a failure and found one so we abort. */
        printf("** stopping\n");
        solution_status_print(log_file, &solution_status, "** stopping");

        break;
      }
    }

    if ((solution_status.success_count + solution_status.fail_count) >= 1000) {
      if (SEARCH_STATE_PIVOT == state) {
        if (SEARCH_STRATEGY_ADAPTIVE == arguments->search_strategy) {
          state = SEARCH_STATE_DECREASE;
        } else {
          printf("** stopping\n");
          solution_status_print(log_file, &solution_status, "** stopping");

          break;
        }
      }

      if (SEARCH_STATE_DECREASE == state) {
        if (solution_status.n <= 1) {
          printf("** stopping\n");
          solution_status_print(log_file, &solution_status, "** stopping");

          break;
        }

        printf("** decrementing n\n");
        solution_status_print(log_file, &solution_status, "** decrementing");
        solution_status_reset(&solution_status);

        solution_status.n--;
      } else if (SEARCH_STATE_INCREASE == state) {
        printf("** stopping\n");
        solution_status_print(log_file, &solution_status, "** stopping");

        break;
      }
    }

    /* If the number of issued jobs is already at 1000, sleep the node so that
     * we can collect jobs from the nodes that are still processing problems. */
    if (solution_status.issued_count >= 1000) {
      /* Sleep the node whilst we wait for the other nodes to complete. */
      uint32_t job = MPI_JOB_SLEEP;

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
    } else {
      /* Issue a new job for the current n. */
      solution_status.issued_count++;

      /* Send a new job to the node. */
      uint32_t job = MPI_JOB_SOLVE_SYSTEM;

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

      /* Send n to the client */
      if (MPI_SUCCESS != MPI_Send(
          &(solution_status.n),
          1, /* count */
          MPI_UNSIGNED,
          status.MPI_SOURCE, /* destination */
          MPI_TAG_SAMPLE_COUNT,
          MPI_COMM_WORLD))
      {
        critical("main_server(): Failed to send n.");
      }
    }
  }

  /* Stop the node. */
  uint32_t job = MPI_JOB_STOP;

  if (MPI_SUCCESS != MPI_Send(
      &job,
      1, /* count */
      MPI_UNSIGNED,
      status.MPI_SOURCE,
      MPI_TAG_JOB,
      MPI_COMM_WORLD))
  {
    critical("main_server(): Failed to send job.");
  }

  nodes--; /* Decrement the node count. */

  printf("Stopping node %d with %u node(s) remaining...\n",
    status.MPI_SOURCE, nodes);

  /* Stop all remaining nodes as they report back. */
  while (nodes > 0) {
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

    if (MPI_NOTIFY_SOLUTION_DONE == notification) {
      uint32_t solution[6];

      if (MPI_SUCCESS != MPI_Recv(
          solution,
          6, /* count */
          MPI_UNSIGNED,
          status.MPI_SOURCE,
          MPI_TAG_SOLUTION,
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE))
      {
        critical("main_server(): Failed to receive solution.");
      }

      /* Ignore the solution as we are to shut down the node. */
    } else if (MPI_NOTIFY_SAMPLING_FAILED_OUT_OF_BOUNDS == notification) {
      uint32_t tmp_n;

      if (MPI_SUCCESS != MPI_Recv(
          &tmp_n,
          1, /* count */
          MPI_UNSIGNED,
          status.MPI_SOURCE,
          MPI_TAG_SAMPLE_COUNT,
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE))
      {
        critical("main_server(): Failed to receive n.");
      }

      (void)tmp_n;

      /* Ignore the sampling error as we are to shut down the node. */
    } else if (MPI_NOTIFY_READY != notification) {
      critical("main_server(): Received unknown notification.");
    }

    /* Stop the node. */
    job = MPI_JOB_STOP;

    if (MPI_SUCCESS != MPI_Send(
        &job,
        1, /* count */
        MPI_UNSIGNED,
        status.MPI_SOURCE,
        MPI_TAG_JOB,
        MPI_COMM_WORLD))
    {
      critical("main_server(): Failed to send job.");
    }

    nodes--; /* Decrement the node count. */

    printf("Stopping node %d with %u node(s) remaining...\n",
      status.MPI_SOURCE, nodes);
  }

  /* Close the file. */
  fclose(log_file);
  log_file = NULL;
}

/*!
 * \brief   The main function on the client node.
 *
 * This function is called once by main() for each distribution to process.
 *
 * \param[in] arguments     The parsed command line arguments.
 */
static void main_client(
  const Solve_Diagonal_Distribution_Arguments * const arguments)
{
  /* Declare variables. */
  mpz_t alpha;
  mpz_init(alpha);

  /* Receive broadcast of the distribution. */
  Diagonal_Distribution distribution;
  diagonal_distribution_init_bcast_recv(&distribution, MPI_RANK_ROOT);

  /* Setup a random state. */
  Random_State random_state;
  random_init(&random_state);

  /* Notify the server that we are ready to receive a new job. */
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

  while (TRUE) {
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
    } else if (MPI_JOB_SLEEP == job) {
      /* Sleep for one second. */
      sleep(1);

      /* Notify the server that we woke up and request a new job. */
      notification = MPI_NOTIFY_READY;

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
    } else if (MPI_JOB_SOLVE_SYSTEM != job) {
      critical("main_client(): Unknown job (%u).", job);
    }

    /* Receive the number n of pairs (j, k) to sample. */
    uint32_t n;

    if (MPI_SUCCESS != MPI_Recv(
        &n,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_TAG_SAMPLE_COUNT,
        MPI_COMM_WORLD,
        &status))
    {
      critical("main_client(): Failed to receive n.");
    }

    /* Sample. */
    Timer timer_prepare_system;
    timer_start(&timer_prepare_system);

    mpz_t * js = (mpz_t *)malloc(n * sizeof(mpz_t));
    mpz_t * ks = (mpz_t *)malloc(n * sizeof(mpz_t));
    int32_t * etas = (int32_t *)malloc(n * sizeof(int32_t));
    if ((NULL == js) || (NULL == ks) || (NULL == etas)) {
      critical("main_client(): Failed to allocate memory.");
    }

    for (uint32_t i = 0; i < n; i++) {
      mpz_init(js[i]);
      mpz_init(ks[i]);
    }

    bool result = TRUE;

    for (uint32_t i = 0; i < n; i++) {
      result = diagonal_distribution_sample_pair_j_k(
                  &distribution,
                  &random_state,
                  arguments->delta_bound,
                  js[i],
                  ks[i],
                  &etas[i],
                  NULL);
      if (TRUE != result) {
        break;
      }

      /* If a bound on eta was not explicitly specified, we accept all eta
       * contained in the distribution, and hence skip the below check. */
      if (TRUE == arguments->eta_bound_specified) {
        if (abs_i(etas[i]) > arguments->eta_bound) {
          /* Consider this a sampling error as if the distribution did not
           * contain the slice sampled. Break the loop. */

          result = FALSE;
          break;
        }
      }
    }

    if (TRUE != result) {
      /* Clear memory. */
      for (uint32_t i = 0; i < n; i++) {
        mpz_clear(js[i]);
        mpz_clear(ks[i]);
      }

      free(js);
      js = NULL;

      free(ks);
      ks = NULL;

      free(etas);
      etas = NULL;

      /* Notify the server that the attempt to sample was out of bounds. */
      notification = MPI_NOTIFY_SAMPLING_FAILED_OUT_OF_BOUNDS;

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

      if (MPI_SUCCESS != MPI_Send(
          &n,
          1, /* count */
          MPI_UNSIGNED,
          MPI_RANK_ROOT,
          MPI_TAG_SAMPLE_COUNT,
          MPI_COMM_WORLD))
      {
        critical("main_client(): Failed to send n.");
      }

      continue;
    }

    const uint64_t time_prepare_system_us = timer_stop(&timer_prepare_system);

    /* Solve system. */
    Timer timer_solve_system;
    timer_start(&timer_solve_system);

    Lattice_Status_Recovery recovery_status;

    if (SOLUTION_METHOD_CLOSEST == arguments->solution_method) {
      lattice_solve_for_d_given_r(
        &recovery_status,
        js,
        ks,
        etas,
        n,
        &(distribution.parameters),
        arguments->reduction_algorithm,
        2 * (distribution.parameters.m + distribution.parameters.sigma));
    } else if (SOLUTION_METHOD_ENUMERATE == arguments->solution_method) {
      lattice_enumerate_for_d_given_r(
        &recovery_status,
        js,
        ks,
        etas,
        n,
        &(distribution.parameters),
        arguments->reduction_algorithm,
        2 * (distribution.parameters.m + distribution.parameters.sigma),
        arguments->timeout);
    } else {
      critical("main_client(): Unknown solution method.");
    }

    const uint64_t time_solve_system_us = timer_stop(&timer_solve_system);

    /* Notify the server of the result and request a new job. */
    notification = MPI_NOTIFY_SOLUTION_DONE;

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

    uint32_t solution[6] = {
      n,
      (uint32_t)recovery_status,
      (uint32_t)(time_prepare_system_us >> 32), /* high 32 bits */
      (uint32_t)(time_prepare_system_us & 0xffffffffUL), /* low 32 bits */
      (uint32_t)(time_solve_system_us >> 32), /* high 32 bits */
      (uint32_t)(time_solve_system_us & 0xffffffffUL) /* low 32 bits */
    };

    if (MPI_SUCCESS != MPI_Send(
        solution,
        6, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_TAG_SOLUTION,
        MPI_COMM_WORLD))
    {
      critical("main_client(): Failed to send solution.");
    }

    /* Clear memory. */
    for (uint32_t i = 0; i < n; i++) {
      mpz_clear(js[i]);
      mpz_clear(ks[i]);
    }

    free(js);
    js = NULL;

    free(ks);
    ks = NULL;

    free(etas);
    etas = NULL;
  }

  /* Clear memory. */
  diagonal_distribution_clear(&distribution);

  random_close(&random_state);

  mpz_clear(alpha);
}

/*!
 * \name Synopsis
 * \{
 */

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  fprintf(file, "Synopsis: mpirun solve_diagonal_distribution \\\n"
          "   [ -delta-bound <delta-bound> ] [ -eta-bound <eta-bound> ] \\\n"
          "      [ -adaptive | -non-adaptive | -non-adaptive-early-abort ] \\\n"
          "         [ -closest | -enumerate ] [ -timeout <timeout> ] \\\n"
          "            [ -lll | -lll-then-bkz | -bkz | -hkz ] \\\n"
          "               <distribution> <n> { <distribution> <n> }\n");

  fprintf(file, "\n");
  fprintf(file, "Delta bound: -- defaults to %u\n", BOUND_DELTA);
  fprintf(file,
    " -delta-bound  explicitly set the delta bound to <delta-bound>\n");

  fprintf(file, "\n");
  fprintf(file, "Eta bound: -- "
    "defaults to the eta bound for the distribution\n");
  fprintf(file,
    " -eta-bound    explicitly set the eta bound to <eta-bound>\n");

  fprintf(file, "\n");
  fprintf(file, "Search strategy: -- defaults to adaptive\n");
  fprintf(file,
    " -adaptive      increment or decrement n to find the minimum\n"
    " -non-adaptive  attempt to solve only for the n specified\n"
    " -non-adaptive-early-abort  abort if too many failures\n");

  fprintf(file, "\n");
  fprintf(file, "Solution method: -- defaults to closest\n");
  fprintf(file,
    " -closest       use the closest vector to a given vector\n"
    " -enumerate     enumerate vectors around a given vector\n");

  fprintf(file, "\n");
  fprintf(file, "Reduction algorithm: -- defaults to LLL then BKZ\n");
  fprintf(file,
    " -lll           use Lenstra-Lenstra-Lovász (LLL)\n"
    " -bkz           use block Korkin-Zolotarev (BKZ)\n"
    " -hkz           use Hermite Korkin-Zolotarev (HKZ)\n"
    " -lll-then-bkz  first try with LLL; if it fails with BKZ\n");

  fprintf(file, "\n");
  fprintf(file, "Timeout: -- defaults to 300 s = 5 min\n");
  fprintf(file,
    " -timeout  explicitly set the timeout to <timeout> seconds\n");
}

/*!
 * \}
 */

/*!
 * \brief The main entry point to the solve_diagonal_distribution executable.
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

  Solve_Diagonal_Distribution_Arguments arguments;
  memset(&arguments, 0, sizeof(arguments));

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

    /* Broadcast the command line arguments. */
    arguments_bcast_send(&arguments, MPI_RANK_ROOT);

    /* Setup a path array for the distribution loader. */
    char ** paths = (char **)malloc(arguments.count * sizeof(char *));
    if (NULL == paths) {
      critical("main(): Failed to allocate memory.");
    }

    for (uint32_t i = 0; i < arguments.count; i++) {
      paths[i] = arguments.entries[i].path;
    }

    /* Setup the distribution loader. */
    Diagonal_Distribution_Loader loader;
    diagonal_distribution_loader_init(&loader, paths, arguments.count);

    /* Process the distributions. */
    for (uint32_t i = 0; i < (arguments.count); i++) {
      Diagonal_Distribution * distribution =
        diagonal_distribution_loader_pop(&loader);

      printf("Processing: %s\n", paths[i]);

      if (TRUE == arguments.eta_bound_specified) {
        if (arguments.eta_bound > distribution->parameters.eta_bound) {
          critical("main(): The eta bound specified via -eta-bound is greater "
            "than the eta bound used to generate the distribution.");
        }
      }

      main_server(
        &arguments,
        &(arguments.entries[i]),
        distribution,
        mpi_size);

      diagonal_distribution_clear(distribution);
    }

    diagonal_distribution_loader_clear(&loader);

    free(paths);
    paths = NULL;

    arguments_clear(&arguments);
  } else {
    /* Broadcast the command line arguments. */
    arguments_init_bcast_recv(&arguments, MPI_RANK_ROOT);

    /* Process the distributions. */
    for (uint32_t i = 0; i < (arguments.count); i++) {
      main_client(&arguments);
    }

    arguments_clear(&arguments);
  }

  if (MPI_SUCCESS != MPI_Finalize()) {
    critical("main(): Failed to finalize MPI.");
  }

  /* Signal success. */
  return 0;
}
