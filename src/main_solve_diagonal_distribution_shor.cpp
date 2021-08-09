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

#include "diagonal_distribution.h"
#include "diagonal_distribution_loader.h"
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
 * \ingroup solve_diagonal_distribution_shor_exe
 */
typedef struct {
  /*!
   * \brief   A bound B on t, to use the notation from the paper, that relates
   *          to the maximum offset in alpha_d searched.
   */
  uint32_t search_bound;

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
  Solve_Diagonal_Distribution_Arguments * const arguments,
  const int argc,
  char ** argv)
{
  /* Initialize the arguments data structure. */
  arguments->search_bound = 0;
  arguments->paths = NULL;
  arguments->count = 0;

  bool search_bound_specified = FALSE;

  /* Iterate over the command line arguments. */
  int i;

  for (i = 1; i < argc; i++) {
    /* Parse the search bound. */
    if (0 == strcmp(argv[i], "-search-bound")) {
      /* Check that a search bound has not already been specified. */
      if (search_bound_specified) {
        fprintf(stderr, "Error: The search bound cannot be twice specified.\n");
        return FALSE;
      }

      if ((i + 1) >= argc) {
        fprintf(stderr, 
          "Error: Expected value to follow after -search-bound.\n");
        return FALSE;
      }

      const int x = atoi(argv[i+1]);
      if (x < 0) {
        fprintf(stderr, "Error: The search bound must be non-negative.\n");
        return FALSE;
      }

      /* Store the search bound. */
      arguments->search_bound = (uint32_t)x;
      search_bound_specified = TRUE;

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
  const Solve_Diagonal_Distribution_Arguments * const arguments,
  const int root)
{
  int result;

  uint32_t data[2];
  data[0] = (uint32_t)(arguments->search_bound);
  data[1] = (uint32_t)(arguments->count);

  result = MPI_Bcast(data, 2, MPI_UNSIGNED, root, MPI_COMM_WORLD);
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

static void arguments_init_bcast_recv(
  Solve_Diagonal_Distribution_Arguments * const arguments,
  const int root)
{
  arguments->paths = NULL;
  arguments->count = 0;

  int result;

  uint32_t data[2];

  result = MPI_Bcast(data, 2, MPI_UNSIGNED, root, MPI_COMM_WORLD);
  if (MPI_SUCCESS != result) {
    critical("arguments_init_bcast_recv(): "
      "Failed to receive broadcast of collected metadata.");
  };

  arguments->search_bound = data[0];
  arguments->count = data[1];

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
 * \param[in] distribution  The distribution to use for sampling problem
 *                          instances to solve.
 * \param[in] path          The path to the distribution.
 * \param[in] mpi_size      The number of nodes.
 */
static void main_server(
  const Diagonal_Distribution * const distribution,
  const char * const path,
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

  /* Broadcast the distribution. */
  diagonal_distribution_bcast_send(distribution, MPI_RANK_ROOT);

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
  const Solve_Diagonal_Distribution_Arguments * const arguments)
{
  /* Receive broadcast of the distribution. */
  Diagonal_Distribution distribution;
  diagonal_distribution_init_bcast_recv(&distribution, MPI_RANK_ROOT);

  /* Retrieve the search bound. */
  const uint32_t search_bound = arguments->search_bound;

  /* Declare constants. */
  mpz_t pow2_lm;
  mpz_init(pow2_lm);
  mpz_set_ui(pow2_lm, 0);
  mpz_setbit(pow2_lm, distribution.parameters.l + distribution.parameters.m);

  /* Declare variables. */
  mpz_t tmp;
  mpz_init(tmp);

  mpfr_t tmp_f;
  mpfr_init2(tmp_f,
      2 * (distribution.parameters.l + distribution.parameters.m));

  mpz_t j;
  mpz_init(j);

  mpz_t k;
  mpz_init(k);

  mpz_t z;
  mpz_init(z);

  mpz_t rp;
  mpz_init(rp);

  mpz_t tau;
  mpz_init(tau);

  mpz_t target_d;
  mpz_init(target_d);

  mpz_t candidate_d;
  mpz_init(candidate_d);

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

    bool result;

    result = diagonal_distribution_sample_pair_j_k(
                &distribution,
                &random_state,
                j,
                k);
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
    mod_reduce(tmp, pow2_lm); /* tmp = {rj}_{2^(m + l)} */
    mpz_sub(z, z, tmp); /* z = rj - {rj}_{2^(m + l)} */
    mpz_div(z, z, pow2_lm); /* z = (rj - {rj}_{2^(m + l)}) / 2^(m + l) */
  
    /* Compute r' = r / tau for the smallest tau such that gcd(r', z) = 1. */
    mpz_set(rp, distribution.parameters.r);
    mpz_set_ui(tau, 1);

    while (TRUE) {
        mpz_gcd(tmp, rp, z);
        if (mpz_cmp_ui(tmp, 1) == 0) {
          break;
        }

        mpz_mul(tau, tau, tmp);
        mpz_div(rp, rp, tmp);
    }
  
    /* Set z = z^-1 (mod r'). */
    if (0 == mpz_invert(z, z, rp)) {
      /* Note: This should never occur given how we select r'. */
      critical("main_client(): The inverse of z does not exist.");
    }

    /* Set target_d = d mod r'. */
    mpz_mod(target_d, distribution.parameters.d, rp);

    mpfr_set_z(tmp_f, distribution.parameters.r, MPFR_RNDN);
    mpfr_mul_z(tmp_f, tmp_f, k, MPFR_RNDN);
    mpfr_div_z(tmp_f, tmp_f, pow2_lm, MPFR_RNDN);
    mpfr_round(tmp_f, tmp_f);
    mpfr_get_z(tmp, tmp_f, MPFR_RNDN);
    mpz_neg(tmp, tmp); /* tmp = -round(r * k / 2^(m + l)) */
  
    /* Search for the candidate d. */
    bool found = FALSE;

    for (uint32_t t = 0; t <= search_bound; t++) {
      /* Non-negative offset in t. */
      mpz_add_ui(candidate_d, tmp, t);
      mpz_mul(candidate_d, candidate_d, z);
      mpz_mod(candidate_d, candidate_d, rp);

      if (mpz_cmp(candidate_d, target_d) == 0) {
        found = TRUE;
        break;
      }

      if (0 == t) {
        continue; /* The offset in t is zero. */
      }

      /* Negative offset in t. */
      mpz_sub_ui(candidate_d, tmp, t);
      mpz_mul(candidate_d, candidate_d, z);
      mpz_mod(candidate_d, candidate_d, rp);

      if (mpz_cmp(candidate_d, target_d) == 0) {
        found = TRUE;
        break;
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
  mpz_clear(pow2_lm);
  mpz_clear(z);
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
            "   [ -search-bound <bound> ] \\\n"
            "      <distribution> { <distribution> }\n");
  fprintf(file, "\n");
  fprintf(file, "Search bound: -- defaults to zero\n");
  fprintf(file,
    " -search-bound     sets the search <bound> on |t| related to the\n"
    "                   maximum offset in alpha_d searched\n");
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
      main_server(
        distribution,
        arguments.paths[i],
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
