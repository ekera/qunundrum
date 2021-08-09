/*!
 * \file    main_solve_linear_distribution_shor.cpp
 * \ingroup solve_linear_distribution_shor_exe
 *
 * \brief   The definition of the main entry point to the
 *          solve_linear_distribution_shor executable, and of associated
 *          functions.
 */

/*!
 * \defgroup solve_linear_distribution_shor_exe \
 *           The solve_linear_distribution_shor executable
 * \ingroup  solve_executable
 *
 * \brief    A module for the solve_linear_distribution_shor executable.
 */

#include "executables.h"
#include "executables_solve_distribution.h"

#include "linear_distribution.h"
#include "linear_distribution_loader.h"
#include "continued_fractions.h"
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
 * \brief   A data structure representing parsed command line arguments.
 *
 * \ingroup solve_linear_distribution_shor_exe
 */
typedef struct {
  /*!
   * \brief   Search bound on integers added to or subtracted from j.
   */
  uint32_t search_bound_j;

  /*!
   * \brief   Search bound on cofactors when solving using continued fractions.
   */
  uint32_t search_bound_cofactors;

  /*!
   * \brief   The number of distributions to process.
   */
  uint32_t count;

  /*!
   * \brief   The path to each distribution to process.
   */
  char ** paths;
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
  arguments->search_bound_j = 256; /* = 2^8 */
  arguments->search_bound_cofactors = 65536; /* = 2^16 */
  arguments->paths = NULL;
  arguments->count = 0;

  bool search_bound_j_specified = FALSE;
  bool search_bound_cofactors_specified = FALSE;

  /* Iterate over the command line arguments. */
  int i;

  for (i = 1; i < argc; i++) {
    /* Parse the search bound on j. */
    if (0 == strcmp(argv[i], "-search-bound-j")) {
      /* Check that a search bound has not already been specified. */
      if (search_bound_j_specified) {
        fprintf(stderr,
          "Error: The search bound on j cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr, 
          "Error: Expected value to follow after -search-bound-j.\n");
        return FALSE;
      }

      const int x = atoi(argv[i+1]);
      if (x < 0) {
        fprintf(stderr, "Error: The search bound on j must be non-negative.\n");
        return FALSE;
      }

      /* Store the search bound. */
      arguments->search_bound_j = (uint32_t)x;
      search_bound_j_specified = TRUE;

      i++;

      continue;
    }

    /* Parse the search bound on cofactors. */
    if (0 == strcmp(argv[i], "-search-bound-cofactors")) {
      /* Check that a search bound has not already been specified. */
      if (search_bound_cofactors_specified) {
        fprintf(stderr,
          "Error: The search bound on cofactors cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr, 
          "Error: Expected value to follow after -search-bound-cofactors.\n");
        return FALSE;
      }

      const int x = atoi(argv[i+1]);
      if (x <= 0) {
        fprintf(stderr, 
          "Error: The search bound on cofactors must be positive.\n");
        return FALSE;
      }

      /* Store the search bound. */
      arguments->search_bound_cofactors = (uint32_t)x;
      search_bound_cofactors_specified = TRUE;

      i++;

      continue;
    }
  
    /* Stop parsing. */
    break;
  }

  /* Parse entries on the form {<distribution>}. */
  if ((argc - i) <= 0) {
    fprintf(stderr, "Error: Incorrect command line arguments; expected "
      "a least one <distribution>.\n");
    return FALSE;
  }

  arguments->count = (uint32_t)(argc - i);

  arguments->paths = (char **)malloc((arguments->count) * sizeof(char *));
  if (NULL == arguments->paths) {
    critical("arguments_init_parse_command_line(): Failed to allocate memory.");
  }

  for (uint32_t j = 0; j < arguments->count; j++) {
    arguments->paths[j] = NULL;
  }
  
  /* Iterate over the distributions. */
  for (uint32_t j = 0; j < arguments->count; j++, i++) {
    /* Parse the path. */
    const uint32_t length = (uint32_t)strlen(argv[i]);

    arguments->paths[j] = (char *)malloc((length + 1) * sizeof(char));
    if (NULL == arguments->paths[j]) {
      critical("arguments_init_parse_command_line(): "
        "Failed to allocate memory.");
    }

    safe_strlcpy(arguments->paths[j], argv[i], length + 1);

    /* Check that the distribution exists. */
    if (0 != access(arguments->paths[j], F_OK)) {
      fprintf(stderr, "Error: The distribution \"%s\" does not exist.\n",
        arguments->paths[j]);
      return FALSE;
    }

    if (0 != access(arguments->paths[j], R_OK)) {
      fprintf(stderr, "Error: The distribution \"%s\" is not readable.\n",
        arguments->paths[j]);
      return FALSE;
    }
  }

  /* Signal success. */
  return TRUE;
}

static void arguments_bcast_send(
  const Solve_Linear_Distribution_Arguments * const arguments,
  const int root)
{
  int result;

  uint32_t data[3];
  data[0] = (uint32_t)(arguments->search_bound_j);
  data[1] = (uint32_t)(arguments->search_bound_cofactors);
  data[2] = (uint32_t)(arguments->count);

  result = MPI_Bcast(data, 3, MPI_UNSIGNED, root, MPI_COMM_WORLD);
  if (MPI_SUCCESS != result) {
    critical("arguments_bcast_send(): "
      "Failed to send broadcast of collected metadata.");
  };

  for (uint32_t i = 0; i < (arguments->count); i++) {
    uint32_t length = (uint32_t)strlen(arguments->paths[i]);

    result = MPI_Bcast(
                &length,
                1,
                MPI_UNSIGNED,
                root,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("arguments_bcast_send(): "
        "Failed to send broadcast of the path length.");
    };

    result = MPI_Bcast(
                arguments->paths[i],
                length + 1,
                MPI_CHAR,
                root,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("arguments_bcast_send(): "
        "Failed to send broadcast of the path.");
    };
  }
}

static void arguments_init_bcast_recv(
  Solve_Linear_Distribution_Arguments * const arguments,
  const int root)
{
  arguments->paths = NULL;
  arguments->count = 0;

  int result;

  uint32_t data[3];

  result = MPI_Bcast(data, 3, MPI_UNSIGNED, root, MPI_COMM_WORLD);
  if (MPI_SUCCESS != result) {
    critical("arguments_init_bcast_recv(): "
      "Failed to receive broadcast of collected metadata.");
  };

  arguments->search_bound_j = data[0];
  arguments->search_bound_cofactors = data[1];
  arguments->count = data[2];

  arguments->paths = (char **)malloc((arguments->count) * sizeof(char *));
  if (NULL == arguments->paths) {
    critical("arguments_init_bcast_recv(): Failed to allocate memory.");
  }

  for (uint32_t i = 0; i < arguments->count; i++) {
    arguments->paths[i] = NULL;
  }

  for (uint32_t i = 0; i < arguments->count; i++) {
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

    arguments->paths[i] = (char *)malloc((length + 1) * sizeof(char));
    if (NULL == arguments->paths[i]) {
      critical("arguments_init_bcast_recv(): Failed to allocate memory.");
    }

    result = MPI_Bcast(
                arguments->paths[i],
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
  if (NULL != arguments->paths) {
    for (uint32_t i = 0; i < arguments->count; i++) {
      if (NULL != arguments->paths[i]) {
        free(arguments->paths[i]);
        arguments->paths[i] = NULL;
      }
    }

    free(arguments->paths);
    arguments->paths = NULL;
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
 * \param[in] distribution The distribution to use for sampling problem
 *                         instances to solve.
 * \param[in] path         The path to the distribution.
 * \param[in] mpi_size     The number of nodes.
 */
static void main_server(
  const Linear_Distribution * const distribution,
  const char * const path,
  const int mpi_size)
{
  const uint32_t flags = distribution->flags;
  if ((flags & LINEAR_DISTRIBUTION_FLAG_R) == 0) {
    critical("main_server(): Expected a distribution for order finding.");
  }

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
    "%s/solve-linear-shor.txt", LOGS_DIRECTORY);

  FILE * log_file = fopen(log_path, "a+");
  if (NULL == log_file) {
    critical("main_server(): Failed to open \"%s\" for appending.", log_path);
  }

  fprintf(log_file, "\n# Processing: %s\n", truncate_path(path));

  /* Broadcast the distribution. */
  linear_distribution_bcast_send(distribution, MPI_RANK_ROOT);

  /* Setup a solution status data structure. */
  Solution_Status solution_status;
  solution_status_init(&solution_status, &(distribution->parameters), 1);

  /* The number of currently running client nodes. */
  int nodes = mpi_size - 1;

  /* Declare status. */
  MPI_Status status;

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
      /* Update the solution status data structure. */
      solution_status.fail_count++;
      solution_status.fail_out_of_bounds_count++;
    } else if (MPI_NOTIFY_SOLUTION_DONE == notification) {
      uint32_t solution[5];

      if (MPI_SUCCESS != MPI_Recv(
          solution,
          5, /* count */
          MPI_UNSIGNED,
          status.MPI_SOURCE,
          MPI_TAG_SOLUTION,
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE))
      {
        critical("main_server(): Failed to receive solution.");
      }

      /* Update the solution status data structure. */

      /* Update the timing statistics. */
      uint64_t time_prepare_system_us;

      time_prepare_system_us   = (uint64_t)solution[1]; /* high 32 bits */
      time_prepare_system_us <<= 32;
      time_prepare_system_us  ^= (uint64_t)solution[2]; /* low 32 bits */

      timer_statistics_insert(
        &(solution_status.statistics_prepare_system),
        time_prepare_system_us);

      uint64_t time_solve_system_us;

      time_solve_system_us   = (uint64_t)solution[3]; /* high 32 bits */
      time_solve_system_us <<= 32;
      time_solve_system_us  ^= (uint64_t)solution[4]; /* low 32 bits */

      timer_statistics_insert(
        &(solution_status.statistics_solve_system),
        time_solve_system_us);

      const bool found = solution[0];

      if (found) {
        solution_status.success_count++;
      } else {
        solution_status.fail_count++;
      }
    } else if (MPI_NOTIFY_READY != notification) {
      critical("main_server(): Received unknown notification.");
    }

    /* Report status. */
    if (MPI_NOTIFY_READY != notification) {
      solution_status_print(stdout, &solution_status, NULL);
    }

    if ((solution_status.success_count + solution_status.fail_count) >= 1000) {
      solution_status_print(log_file, &solution_status, NULL);
      break;
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
      uint32_t solution[5];

      if (MPI_SUCCESS != MPI_Recv(
          solution,
          5, /* count */
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
  if ((flags & LINEAR_DISTRIBUTION_FLAG_R) == 0) {
    critical("main_client(): Expected a distribution for order finding.");
  }

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

  mpz_t j;
  mpz_init(j);

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

    /* Sample. */
    Timer timer_prepare_system;
    timer_start(&timer_prepare_system);

    bool result;

    result = linear_distribution_sample_alpha(
                &distribution,
                &random_state,
                alpha);
    if (TRUE != result) {
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

      continue;
    }

    sample_j_from_alpha_r(
      j,
      alpha,
      &(distribution.parameters),
      &random_state);

    const uint64_t time_prepare_system_us = timer_stop(&timer_prepare_system);

    /* Solve system. */
    Timer timer_solve_system;
    timer_start(&timer_solve_system);

    bool found = FALSE;

    mpz_t candidate_j;
    mpz_init(candidate_j);
    
    for (uint32_t t = 0; t <= arguments->search_bound_j; t++) {
      mpz_add_ui(candidate_j, j, t);

      found = continued_fractions_solve(
                candidate_j,
                arguments->search_bound_cofactors,
                &(distribution.parameters));
      if (TRUE == found) {
        break;
      }

      if (0 == t) {
        continue;
      }

      mpz_sub_ui(candidate_j, j, t);

      found = continued_fractions_solve(
                candidate_j,
                arguments->search_bound_cofactors,
                &(distribution.parameters));
      if (TRUE == found) {
        break;
      }
    }

    mpz_clear(candidate_j);

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
      critical("main_client(): Failed to send MPI_NOTIFY_SOLUTION_DONE.");
    }

    uint32_t solution[5] = {
      (uint32_t)found,
      (uint32_t)(time_prepare_system_us >> 32), /* high 32 bits */
      (uint32_t)(time_prepare_system_us & 0xffffffffUL), /* low 32 bits */
      (uint32_t)(time_solve_system_us >> 32), /* high 32 bits */
      (uint32_t)(time_solve_system_us & 0xffffffffUL) /* low 32 bits */
    };

    if (MPI_SUCCESS != MPI_Send(
        solution,
        5, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_TAG_SOLUTION,
        MPI_COMM_WORLD))
    {
      critical("main_client(): Failed to send solution.");
    }
  }

  /* Clear memory. */
  linear_distribution_clear(&distribution);

  random_close(&random_state);

  mpz_clear(alpha);
  mpz_clear(j);
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
  fprintf(file, "Synopsis: mpirun solve_linear_distribution_shor\n"
            "   [ -search-bound-j <bound> ] "
              "[ -search-bound-cofactors <bound> ]\n"
            "      <distribution> { <distribution> }\n");
  fprintf(file, "\n");
  fprintf(file, "Search bound on j: -- defaults to 2^8\n");
  fprintf(file,
    " -search-bound-j          sets the search <bound> on |t| where t is\n"
    "                          we try to solve for all integers j + t\n");
  fprintf(file, "\n");
  fprintf(file, "Search bound on cofactors: -- defaults to 2^16\n");
  fprintf(file,
    " -search-bound-cofactors  sets the search <bound> on the cofactor\n"
    "                          between the order r and denominator when\n"
    "                          taking continued fractions\n");
}

/*!
 * \}
 */

/*!
 * \brief The main entry point to the solve_linear_distribution_shor executable.
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

    /* Setup the distribution loader. */
    Linear_Distribution_Loader loader;
    linear_distribution_loader_init(
        &loader, arguments.paths, arguments.count);

    /* Process the distributions. */
    for (uint32_t i = 0; i < (arguments.count); i++) {
      Linear_Distribution * distribution =
        linear_distribution_loader_pop(&loader);

      printf("Processing: %s\n", arguments.paths[i]);
      main_server(
        distribution,
        arguments.paths[i],
        mpi_size);

      linear_distribution_clear(distribution);
    }

    linear_distribution_loader_clear(&loader);
  
    arguments_clear(&arguments);
  } else {
    /* Broadcast the command line arguments. */
    arguments_init_bcast_recv(&arguments, MPI_RANK_ROOT);

    /* Process the distributions. */
    for (uint32_t i = 0; i < arguments.count; i++) {
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