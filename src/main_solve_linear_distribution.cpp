/*!
 * \file    main_solve_linear_distribution.cpp
 * \ingroup solve_linear_distribution_exe
 *
 * \brief   The definition of the main entry point to the
 *          solve_linear_distribution executable, and of associated functions.
 */

/*!
 * \defgroup solve_linear_distribution_exe \
 *           The solve_linear_distribution executable
 * \ingroup  solve_executable
 *
 * \brief    A module for the solve_linear_distribution executable.
 */

#include "executables.h"
#include "executables_solve_distribution.h"

#include "linear_distribution.h"
#include "linear_distribution_loader.h"
#include "lattice.h"
#include "lattice_solve.h"
#include "lattice_enumerate.h"
#include "random.h"
#include "common.h"
#include "sample.h"
#include "timer.h"
#include "math.h"
#include "errors.h"
#include "string_utilities.h"

#include <gmp.h>
#include <mpfr.h>
#include <mpi.h>

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include <unistd.h>

#include <sys/stat.h>
#include <sys/types.h>

/*!
 * \brief   A data structure representing argument entries in the form of
 *          parsed \<distribution\> \<n\> tuples from the command line
 *          arguments.
 *
 * \ingroup solve_linear_distribution_exe
 */
typedef struct {
  /*!
   * \brief   The parameter m.
   */
  uint32_t n;

  /*!
   * \brief   The path to the distribution.
   */
  char * path;
} Solve_Linear_Distribution_Arguments_Entry;

/*!
 * \brief   A data structure representing parsed command line arguments.
 *
 * \ingroup solve_linear_distribution_exe
 */
typedef struct {
  /*!
   * \brief   The search strategy.
   */
  Search_Strategy search_strategy;

  /*!
   * \brief   The solution method.
   */
  Solution_Method solution_method;

  /*!
   * \brief   The option for detecting smooth orders.
   */
  Detect_Smooth_Order_Option detect_smooth_order;

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
  Solve_Linear_Distribution_Arguments_Entry * entries;
} Solve_Linear_Distribution_Arguments;


/*!
 * \name Command line arguments
 * \{
 */

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
  Solve_Linear_Distribution_Arguments * const arguments,
  const int argc,
  char ** argv)
{
  /* Initialize the arguments data structure. */
  arguments->search_strategy = SEARCH_STRATEGY_DEFAULT;
  arguments->solution_method = SOLUTION_METHOD_DEFAULT;
  arguments->detect_smooth_order = DETECT_SMOOTH_ORDER_DEFAULT;
  arguments->reduction_algorithm = REDUCTION_ALGORITHM_DEFAULT;
  arguments->timeout = 0;
  arguments->entries = NULL;
  arguments->count = 0;

  bool timeout_specified = FALSE;

  /* Iterate over the command line arguments. */
  int i;

  for (i = 1; i < argc; i++) {
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

    /* Parse the option for detecting smooth orders. */
    if (0 == strcmp(argv[i], "-detect-smooth-order")) {
      /* Check that a detection option has not already been specified. */
      if (DETECT_SMOOTH_ORDER_DEFAULT != (arguments->detect_smooth_order)) {
        fprintf(stderr, "Error: The flag for detecting smooth orders cannot be "
          "twice specified.\n");
        return FALSE;
      }

      /* Store the option for detecting smooth orders. */
      arguments->detect_smooth_order = DETECT_SMOOTH_ORDER_TRUE;

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

  /* Set default parameters if arguments where not explicitly specified. */
  if (SOLUTION_METHOD_DEFAULT == arguments->solution_method) {
    arguments->solution_method = SOLUTION_METHOD_CLOSEST;
  }

  if (DETECT_SMOOTH_ORDER_DEFAULT == arguments->detect_smooth_order) {
    arguments->detect_smooth_order = DETECT_SMOOTH_ORDER_FALSE;
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

  /* Parse tuples {<distribution> <n>}. */
  if (((argc - i) <= 0) || (0 != ((argc - i) % 2))) {
    fprintf(stderr, "Error: Incorrect command line arguments; expected tuples "
      "but found an odd number of arguments.\n");
    return FALSE;
  }

  arguments->count = (uint32_t)((argc - i) / 2);

  arguments->entries =
    (Solve_Linear_Distribution_Arguments_Entry *)malloc(
      (arguments->count) * sizeof(Solve_Linear_Distribution_Arguments_Entry));
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

static void arguments_bcast_send(
  const Solve_Linear_Distribution_Arguments * const arguments,
  const int root)
{
  int result;

  uint32_t data[6];
  data[0] = (uint32_t)(arguments->search_strategy);
  data[1] = (uint32_t)(arguments->solution_method);
  data[2] = (uint32_t)(arguments->detect_smooth_order);
  data[3] = (uint32_t)(arguments->reduction_algorithm);
  data[4] = (uint32_t)(arguments->timeout);
  data[5] = (uint32_t)(arguments->count);

  result = MPI_Bcast(data, 6, MPI_UNSIGNED, root, MPI_COMM_WORLD);
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

static void arguments_init_bcast_recv(
  Solve_Linear_Distribution_Arguments * const arguments,
  const int root)
{
  arguments->entries = NULL;
  arguments->count = 0;

  int result;

  uint32_t data[6];

  result = MPI_Bcast(data, 6, MPI_UNSIGNED, root, MPI_COMM_WORLD);
  if (MPI_SUCCESS != result) {
    critical("arguments_init_bcast_recv(): "
      "Failed to receive broadcast of collected metadata.");
  };

  arguments->search_strategy = (Search_Strategy)(data[0]);
  arguments->solution_method = (Solution_Method)(data[1]);
  arguments->detect_smooth_order = (Detect_Smooth_Order_Option)(data[2]);
  arguments->reduction_algorithm = (Lattice_Reduction_Algorithm)(data[3]);
  arguments->timeout = data[4];
  arguments->count = data[5];

  arguments->entries =
    (Solve_Linear_Distribution_Arguments_Entry *)malloc(
      (arguments->count) * sizeof(Solve_Linear_Distribution_Arguments_Entry));
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
 * \brief   Clears an initialized command line arguments data structure.
 * 
 * \param[in, out] arguments   The argument data structure to clear.
 */
static void arguments_clear(
  Solve_Linear_Distribution_Arguments * const arguments)
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

  memset(arguments, 0, sizeof(Solve_Linear_Distribution_Arguments));
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
  const Solve_Linear_Distribution_Arguments * const arguments,
  const Solve_Linear_Distribution_Arguments_Entry * const entry,
  const Linear_Distribution * const distribution,
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
    "%s/solve-linear.txt", LOGS_DIRECTORY);

  FILE * log_file = fopen(log_path, "a+");
  if (NULL == log_file) {
    critical("main_server(): Failed to open \"%s\" for appending.", log_path);
  }

  fprintf(log_file, "\n# Processing: %s\n", truncate_path(entry->path));

  /* Print a warning if the -detect-smooth-order flag is used with distributions
   * for short discrete logarithms as it then has no effect. */
  if (DETECT_SMOOTH_ORDER_TRUE == arguments->detect_smooth_order) {
    if ((distribution->flags & LINEAR_DISTRIBUTION_FLAG_R) != 0) {
      printf("Warning: The -detect-smooth-order flag has no effect for this "
        "distribution.\n");
    }
  }

  /* Broadcast the distribution. */
  linear_distribution_bcast_send(distribution, MPI_RANK_ROOT);

  /* Setup a solution status data structure. */
  Solution_Status solution_status;
  solution_status_init(
    &solution_status,
    distribution->parameters.m,
    0, /* = sigma */
    distribution->parameters.s,
    distribution->parameters.l,
    entry->n,
    FALSE); /* = has_sigma */

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
      uint32_t sample_count;

      if (MPI_SUCCESS != MPI_Recv(
          &sample_count,
          1, /* count */
          MPI_UNSIGNED,
          status.MPI_SOURCE,
          MPI_TAG_SAMPLE_COUNT,
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE))
      {
        critical("main_server(): Failed to receive the sample count.");
      }

      if (sample_count == solution_status.n) {
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
        critical("main_server(): Failed to send the sample count.");
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
      uint32_t samples_count;

      if (MPI_SUCCESS != MPI_Recv(
          &samples_count,
          1, /* count */
          MPI_UNSIGNED,
          status.MPI_SOURCE,
          MPI_TAG_SAMPLE_COUNT,
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE))
      {
        critical("main_server(): Failed to receive the sample count.");
      }

      (void)samples_count;

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
  const Solve_Linear_Distribution_Arguments * const arguments)
{
  /* Declare variables. */
  mpz_t alpha;
  mpz_init(alpha);

  /* Receive broadcast of the distribution. */
  Linear_Distribution distribution;
  linear_distribution_init_bcast_recv(&distribution, MPI_RANK_ROOT);

  const uint32_t flags = distribution.flags;

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

    /* Receive the number of regions n to sample. */
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
      critical("main_client(): Failed to receive the sample count.");
    }

    /* Sample. */
    Timer timer_prepare_system;
    timer_start(&timer_prepare_system);

    mpz_t * js = (mpz_t *)malloc(n * sizeof(mpz_t));
    mpz_t * ks = (mpz_t *)malloc(n * sizeof(mpz_t));
    if ((NULL == js) || (NULL == ks)) {
      critical("main_client(): Failed to allocate memory.");
    }

    for (uint32_t i = 0; i < n; i++) {
      mpz_init(js[i]);
      mpz_init(ks[i]);
    }

    bool result = TRUE;

    for (uint32_t i = 0; i < n; i++) {
      result = linear_distribution_sample_alpha(
                  &distribution,
                  &random_state,
                  alpha);
      if (TRUE != result) {
        break;
      }

      if ((flags & LINEAR_DISTRIBUTION_FLAG_R) != 0) {
        sample_j_from_alpha_r(
          js[i],
          alpha,
          &(distribution.parameters),
          &random_state
        );
      } else {
        sample_j_k_from_alpha_d(
          js[i],
          ks[i],
          alpha,
          &(distribution.parameters),
          &random_state);
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
        critical("main_client(): Failed to send the sample count.");
      }

      continue;
    }

    const uint64_t time_prepare_system_us = timer_stop(&timer_prepare_system);

    /* Solve system. */
    Timer timer_solve_system;
    timer_start(&timer_solve_system);

    Lattice_Status_Recovery recovery_status;

    if ((flags & LINEAR_DISTRIBUTION_FLAG_R) != 0) {
      /* If the flag for detecting smooth orders has been set, pass TRUE for 
       * the detect_smooth_r argument, otherwise pass FALSE. */
      const bool detect_smooth_r = 
        (DETECT_SMOOTH_ORDER_TRUE == arguments->detect_smooth_order) ? 
          TRUE : FALSE;

      if (SOLUTION_METHOD_CLOSEST == arguments->solution_method) {
        lattice_solve_for_r(
          &recovery_status,
          js,
          n,
          &(distribution.parameters),
          arguments->reduction_algorithm,
          detect_smooth_r);
      } else if (SOLUTION_METHOD_ENUMERATE == arguments->solution_method) {
        lattice_enumerate_for_r(
          &recovery_status,
          js,
          n,
          &(distribution.parameters),
          arguments->reduction_algorithm,
          distribution.parameters.m, /* = precision */
          detect_smooth_r,
          arguments->timeout);
      } else {
        critical("main_client(): Unknown solution method.");
      }
    } else {
      /* Note: Detection of smooth orders is not applicable to short discrete
       *       logarithms since the order does not enter into the equation. */

      if (SOLUTION_METHOD_CLOSEST == arguments->solution_method) {
        lattice_solve_for_d(
          &recovery_status,
          js,
          ks,
          n,
          &(distribution.parameters),
          arguments->reduction_algorithm,
          distribution.parameters.m, /* = precision */
          FALSE); /* = detect_smooth_r */
      } else if (SOLUTION_METHOD_ENUMERATE == arguments->solution_method) {
        lattice_enumerate_for_d(
          &recovery_status,
          js,
          ks,
          n,
          &(distribution.parameters),
          arguments->reduction_algorithm,
          distribution.parameters.m, /* = precision */
          FALSE, /* = detect_smooth_r */
          arguments->timeout);
      } else {
        critical("main_client(): Unknown solution method.");
      }
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
  }

  /* Clear memory. */
  linear_distribution_clear(&distribution);

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
  fprintf(file, "Synopsis: mpirun solve_linear_distribution \\\n"
          "   [ -adaptive | -non-adaptive | -non-adaptive-early-abort ] \\\n"
          "      [ -closest | -enumerate ] [ -timeout <timeout> ] \\\n"
          "         [ -detect-smooth-order ] \\\n"
          "            [ -lll | -lll-then-bkz | -bkz | -hkz ] \\\n"
          "               <distribution> <n> { <distribution> <n> }\n");

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
  fprintf(file, "Options for orders: -- defaults to no detection\n");
  fprintf(file,
    " -detect-smooth-order  detect if the order is smooth, and if so\n"
    "                       leverage its smoothness when solving\n");

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
 * \brief The main entry point to the solve_linear_distribution executable.
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

  Solve_Linear_Distribution_Arguments arguments;
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
    Linear_Distribution_Loader loader;
    linear_distribution_loader_init(&loader, paths, arguments.count);

    /* Process the distributions. */
    for (uint32_t i = 0; i < (arguments.count); i++) {
      Linear_Distribution * distribution =
        linear_distribution_loader_pop(&loader);

      printf("Processing: %s\n", paths[i]);
      main_server(
        &arguments,
        &(arguments.entries[i]),
        distribution,
        mpi_size);

      linear_distribution_clear(distribution);
    }

    linear_distribution_loader_clear(&loader);

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
