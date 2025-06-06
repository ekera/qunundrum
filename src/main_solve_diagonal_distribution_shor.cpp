/*!
 * \file    main_solve_diagonal_distribution_shor.cpp
 * \ingroup solve_diagonal_distribution_shor_exe
 *
 * \brief   The definition of the main entry point to the
 *          solve_diagonal_distribution_shor executable, and of associated
 *          functions.
 */

/*!
 * \defgroup solve_diagonal_distribution_shor_exe \
 *           The solve_diagonal_distribution_shor executable
 * \ingroup  solve_executable
 *
 * \brief    A module for the solve_diagonal_distribution_shor executable.
 */

#include "executables.h"
#include "executables_solve_distribution.h"

#include "common.h"
#include "continued_fractions.h"
#include "diagonal_distribution.h"
#include "diagonal_distribution_loader.h"
#include "errors.h"
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
 * \brief   A data structure representing parsed command line arguments.
 *
 * \ingroup solve_diagonal_distribution_shor_exe
 */
typedef struct {
  /*!
   * \brief   An upper bound on t.
   */
  uint32_t t_bound;

  /*!
   * \brief   An upper bound on eta.
   */
  uint32_t eta_bound;

  /*!
   * \brief   The number of distributions to process.
   */
  uint32_t count;

  /*!
   * \brief   The path to each distribution to process.
   */
  char ** paths;
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
  arguments->t_bound = 0;
  arguments->eta_bound = 0;
  arguments->paths = NULL;
  arguments->count = 0;

  bool t_bound_specified = FALSE;
  bool eta_bound_specified = FALSE;

  /* Iterate over the command line arguments. */
  int i;

  for (i = 1; i < argc; i++) {
    /* Parse the t bound. */
    if (0 == strcmp(argv[i], "-t-bound")) {
      /* Check that a t bound has not already been specified. */
      if (FALSE != t_bound_specified) {
        fprintf(stderr, "Error: The t bound cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr,
          "Error: Expected value to follow after -t-bound.\n");
        return FALSE;
      }

      const int x = atoi(argv[i + 1]);
      if (x < 0) {
        fprintf(stderr, "Error: The t bound must be non-negative.\n");
        return FALSE;
      }

      /* Store the t bound. */
      arguments->t_bound = (uint32_t)x;
      t_bound_specified = TRUE;

      i++;

      continue;
    }

    if (0 == strcmp(argv[i], "-eta-bound")) {
      /* Check that a eta bound has not already been specified. */
      if (FALSE != eta_bound_specified) {
        fprintf(stderr, "Error: The eta bound cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr,
          "Error: Expected value to follow after -eta-bound.\n");
        return FALSE;
      }

      const int x = atoi(argv[i + 1]);
      if (x < 0) {
        fprintf(stderr, "Error: The eta bound must be non-negative.\n");
        return FALSE;
      }

      /* Store the eta bound. */
      arguments->eta_bound = (uint32_t)x;
      eta_bound_specified = TRUE;

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

  uint32_t data[3];
  data[0] = (uint32_t)(arguments->t_bound);
  data[1] = (uint32_t)(arguments->eta_bound);
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
  arguments->paths = NULL;
  arguments->count = 0;

  int result;

  uint32_t data[3];

  result = MPI_Bcast(data, 3, MPI_UNSIGNED, root, MPI_COMM_WORLD);
  if (MPI_SUCCESS != result) {
    critical("arguments_init_bcast_recv(): "
      "Failed to receive broadcast of collected metadata.");
  };

  arguments->t_bound = data[0];
  arguments->eta_bound = data[1];
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
  mpfr_fprintf(file, "# Bounds: (eta = %u (%u), t = %u)\n",
    arguments->eta_bound,
    distribution->parameters.eta_bound,
    arguments->t_bound);

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
 * \param[in] arguments     The parsed command line arguments.
 * \param[in] path          The path to the distribution.
 * \param[in] distribution  The distribution to use for sampling problem
 *                          instances to solve.
 * \param[in] mpi_size      The number of nodes.
 */
static void main_server(
  const Solve_Diagonal_Distribution_Arguments * const arguments,
  const char * const path,
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
    "%s/solve-diagonal-shor.txt", LOGS_DIRECTORY);

  FILE * log_file = fopen(log_path, "a+");
  if (NULL == log_file) {
    critical("main_server(): Failed to open \"%s\" for appending.", log_path);
  }

  fprintf(log_file, "\n# Processing: %s\n", truncate_path(path));
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
    1, /* = n */
    TRUE); /* = has_sigma */

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
  const Solve_Diagonal_Distribution_Arguments * const arguments)
{
  /* Receive broadcast of the distribution. */
  Diagonal_Distribution distribution;
  diagonal_distribution_init_bcast_recv(&distribution, MPI_RANK_ROOT);

  /* Retrieve the t and eta bounds. */
  const uint32_t t_bound = arguments->t_bound;
  const uint32_t eta_bound = arguments->eta_bound;

  /* Declare constants. */
  mpz_t pow2_msigma;
  mpz_init(pow2_msigma);
  mpz_set_ui(pow2_msigma, 0);
  mpz_setbit(pow2_msigma,
    distribution.parameters.m + distribution.parameters.sigma);

  mpz_t pow2_l;
  mpz_init(pow2_l);
  mpz_set_ui(pow2_l, 0);
  mpz_setbit(pow2_l, distribution.parameters.l);

  /* Declare variables. */
  const uint32_t precision =
    2 * (max_ui(distribution.parameters.m +
      distribution.parameters.sigma, PRECISION));

  mpz_t tmp;
  mpz_init(tmp);

  mpfr_t tmp_f;
  mpfr_init2(tmp_f, precision);

  mpz_t j;
  mpz_init(j);

  mpz_t k;
  mpz_init(k);

  mpz_t z;
  mpz_init(z);

  mpz_t z_eta;
  mpz_init(z_eta);

  mpz_t rp;
  mpz_init(rp);

  mpz_t tau;
  mpz_init(tau);

  mpz_t target_d;
  mpz_init(target_d);

  mpz_t candidate_d;
  mpz_init(candidate_d);

  mpz_t rk_div_pow2_l;
  mpz_init(rk_div_pow2_l);

  /* Compute delta_bound given t_bound. We require that
   *
   *   round ( r (delta_bound + 1/2) / 2^l ) >= t_bound
   *
   * to ensure that the sampling functions explore all valid values of Delta for
   * the range of values of t that we are to search.
   *
   * This requirement may be re-written as
   *
   *   round ( r (delta_bound + delta_b) / 2^l ) >= t_bound
   *   r (delta_bound + delta_b) / 2^l - 1 >= t_bound
   *   r (delta_bound + delta_b) / 2^l >= t_bound + 1
   *   delta_bound + delta_b >= 2^l (t_bound + 1) / r
   *   delta_bound >= 2^l (t_bound + 1) / r - delta_b
   *
   * which in turn can be simplified to
   *
   *   delta_bound >= 2^l (t_bound + 1) / r + 1/2,
   *
   * and so we pick delta_bound = ceil(2^l (t_bound + 1) / r + 1/2). */

  mpz_set_ui(tmp, t_bound); /* tmp = t_bound */
  mpz_add_ui(tmp, tmp, 1); /* tmp = t_bound + 1 */
  mpfr_set_z(tmp_f, tmp, MPFR_RNDN);  /* tmp_f = tmp = t_bound + 1 */

  mpz_set_ui(tmp, 0); /* tmp = 0 */
  mpz_setbit(tmp, distribution.parameters.l); /* tmp = 2^l */
  mpfr_mul_z(tmp_f, tmp_f, tmp, MPFR_RNDN); /* tmp_f = 2^l (t_bound + 1) */
  mpfr_div_z(tmp_f, tmp_f, distribution.parameters.r, MPFR_RNDN);
    /* tmp_f = 2^l (t_bound + 1) / r */
  mpfr_add_d(tmp_f, tmp_f, 0.5, MPFR_RNDN);
    /* tmp_f = 2^l (t_bound + 1) / r + 1/2 */

  mpfr_ceil(tmp_f, tmp_f); /* tmp_f = ceil(2^l (t_bound + 1) / r + 1/2) */
  mpfr_get_z(tmp, tmp_f, MPFR_RNDN);
    /* tmp = tmp_f = ceil(2^l (t_bound + 1) / r + 1/2) */

  /* Sanity check. */
  if (mpz_cmp_ui(tmp, UINT32_MAX) > 0)
  {
    critical("main_client(): "
      "Overflow when computing delta_bound from t_bound: This may be because "
        "2^l is large in relation to r, or because t_bound is large.");
  }

  /* Set the delta bound. */
  const uint32_t delta_bound = mpz_get_ui(tmp);

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

    /* Sample. */
    Timer timer_prepare_system;
    timer_start(&timer_prepare_system);

    int32_t eta;

    bool result;

    result = diagonal_distribution_sample_pair_j_k(
                &distribution,
                &random_state,
                delta_bound,
                j,
                k,
                &eta, /* Ignore. */
                NULL);
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

    const uint64_t time_prepare_system_us = timer_stop(&timer_prepare_system);

    /* Solve system. */
    Timer timer_solve_system;
    timer_start(&timer_solve_system);

    /* Set z = (rj - {rj}_{2^(m + l)}) / 2^(m + l). */
    mpz_mul(z, distribution.parameters.r, j); /* z = rj */
    mpz_set(tmp, z); /* tmp = z = rj */
    mod_reduce(tmp, pow2_msigma); /* tmp = {rj}_{2^(m + sigma)} */
    mpz_sub(z, z, tmp); /* z = rj - {rj}_{2^(m + sigma)} */
    mpz_div(z, z, pow2_msigma);
      /* z = (rj - {rj}_{2^(m + sigma)}) / 2^(m + sigma) */

    mpfr_set_z(tmp_f, distribution.parameters.r, MPFR_RNDN);
      /* tmp_f = r */
    mpfr_mul_z(tmp_f, tmp_f, k, MPFR_RNDN);
      /* tmp_f = r k */
    mpfr_div_z(tmp_f, tmp_f, pow2_l, MPFR_RNDN);
      /* tmp_f = r k / 2^l */
    mpfr_round(tmp_f, tmp_f);
      /* tmp_f = round(r k / 2^l) */
    mpfr_get_z(rk_div_pow2_l, tmp_f, MPFR_RNDN);
      /* rk_div_pow2l = tmp_f = round(r k / 2^l) */

    /* Search for the candidate d. */
    bool found = FALSE;

    for (uint32_t eta_abs = 0;
          (eta_abs <= eta_bound) && (FALSE == found); eta_abs++) {
      for (int32_t eta_sgn = 1;
            (eta_sgn >= -1) && (FALSE == found); eta_sgn -= 2) {

        if ((0 == eta_abs) && (-1 == eta_sgn)) {
          continue;
        }

        if (1 == eta_sgn) {
          mpz_add_ui(z_eta, z, eta_abs);
            /* z_eta = z + eta_sgn * eta_abs = z + eta_abs */
        } else {
          mpz_sub_ui(z_eta, z, eta_abs);
            /* z_eta = z + eta_sgn * eta_abs = z - eta_abs */
        }

        /* Let r' = r / tau for the least tau such that gcd(r', z_eta) = 1. */
        mpz_set(rp, distribution.parameters.r);
          /* r' = r */
        mpz_set_ui(tau, 1);
          /* tau = 1 */

        while (TRUE) {
            mpz_gcd(tmp, rp, z_eta); /* tmp = gcd(r', z_eta) */
            if (mpz_cmp_ui(tmp, 1) == 0) {
              break;
            }

            mpz_mul(tau, tau, tmp); /* tau = tau * tmp */
            mpz_div(rp, rp, tmp); /* r' = r' / tmp */
        }

        /* Set z_eta = z_eta^-1 (mod r'). */
        mpz_mod(z_eta, z_eta, rp);
        if (0 == mpz_invert(z_eta, z_eta, rp)) {
          /* Note: This should never occur given how we select r'. */
          critical("main_client(): The inverse of z_eta does not exist.");
        }

        /* Set target_d = d mod r'. */
        mpz_mod(target_d, distribution.parameters.d, rp);

        for (uint32_t t = 0; t <= t_bound; t++) {
          /* Non-negative offset in t. */
          mpz_sub_ui(candidate_d, rk_div_pow2_l, t);
            /* candidate_d = round(r k / 2^l) - t */
          mpz_neg(candidate_d, candidate_d);
            /* candidate_d = t - round(r k / 2^l) */
          mpz_mul(candidate_d, candidate_d, z_eta);
            /* candidate_d = (t - round(r k / 2^l)) ((z + eta)^-1 mod r) */
          mpz_mod(candidate_d, candidate_d, rp);
            /* candidate_d = (t - round(r k / 2^l)) (z + eta)^-1 mod r */

          if (mpz_cmp(candidate_d, target_d) == 0) {
            found = TRUE;
            break;
          }

          if (0 == t) {
            continue; /* The offset in t is zero. */
          }

          /* Negative offset in t. */
          mpz_add_ui(candidate_d, rk_div_pow2_l, t);
            /* candidate_d = round(r k / 2^l) + t */
          mpz_neg(candidate_d, candidate_d);
            /* candidate_d = -t - round(r k / 2^l) */
          mpz_mul(candidate_d, candidate_d, z_eta);
            /* candidate_d = (-t - round(r k / 2^l)) ((z + eta)^-1 mod r) */
          mpz_mod(candidate_d, candidate_d, rp);
            /* candidate_d = (-t - round(r k / 2^l)) (z + eta)^-1 mod r */

          if (mpz_cmp(candidate_d, target_d) == 0) {
            found = TRUE;
            break;
          }
        }

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
  diagonal_distribution_clear(&distribution);

  random_close(&random_state);

  mpz_clear(j);
  mpz_clear(k);

  mpz_clear(tmp);
  mpfr_clear(tmp_f);

  mpz_clear(pow2_msigma);
  mpz_clear(pow2_l);
  mpz_clear(rk_div_pow2_l);

  mpz_clear(z);
  mpz_clear(z_eta);
  mpz_clear(rp);
  mpz_clear(tau);
  mpz_clear(target_d);
  mpz_clear(candidate_d);
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
  fprintf(file, "Synopsis: mpirun solve_diagonal_distribution_shor \\\n"
            "   [ -t-bound <t-bound> ] [ -eta-bound <eta-bound> ] \\\n"
            "      <distribution> { <distribution> }\n");

  fprintf(file, "\n");
  fprintf(file, "t bound: -- defaults to 0\n");
  fprintf(file,
    " -t-bound    explicitly set the t bound to <t-bound>\n");

  fprintf(file, "\n");
  fprintf(file, "Eta bound: -- defaults to 0\n");
  fprintf(file,
    " -eta-bound  explicitly set the eta bound to <eta-bound>\n");
}

/*!
 * \}
 */

/*!
 * \brief The main entry point to the solve_diagonal_distribution_shor
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

    /* Setup the distribution loader. */
    Diagonal_Distribution_Loader loader;
    diagonal_distribution_loader_init(
        &loader, arguments.paths, arguments.count);

    /* Process the distributions. */
    for (uint32_t i = 0; i < arguments.count; i++) {
      Diagonal_Distribution * distribution =
        diagonal_distribution_loader_pop(&loader);

      printf("Processing: %s\n", arguments.paths[i]);

      if (arguments.eta_bound > distribution->parameters.eta_bound) {
        critical("main(): The eta bound specified via -eta-bound is greater "
          "than the eta bound used to generate the distribution.");
      }

      main_server(
        &arguments,
        arguments.paths[i],
        distribution,
        mpi_size);

      diagonal_distribution_clear(distribution);
    }

    diagonal_distribution_loader_clear(&loader);

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
