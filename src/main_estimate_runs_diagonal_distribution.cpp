/*!
 * \file    main_estimate_runs_diagonal_distribution.cpp
 * \ingroup estimate_runs_diagonal_distribution_exe
 *
 * \brief   The definition of the main entry point to the
 *          estimate_runs_diagonal_distribution executable, and of associated
 *          functions.
 */

/*!
 * \defgroup estimate_runs_diagonal_distribution_exe \
 *           The estimate_runs_diagonal_distribution executable
 * \ingroup  volume_executable
 *
 * \brief    A module for the estimate_runs_diagonal_distribution executable.
 */

#include "executables.h"
#include "executables_estimate_runs_distribution.h"

#include "common.h"
#include "errors.h"
#include "diagonal_distribution.h"
#include "diagonal_distribution_loader.h"
#include "log.h"
#include "random.h"
#include "string_utilities.h"
#include "tau_estimate.h"
#include "tau_ordered_list.h"
#include "tau_volume_quotient.h"

#include <mpfr.h>

#include <mpi.h>

#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <unistd.h>
#include <sys/stat.h>

/*!
 * \brief   A data structure representing parsed command line arguments.
 *
 * \ingroup estimate_runs_diagonal_distribution_exe
 */
typedef struct {
  /*!
   * \brief   The bound on the volume quotient.
   */
  mpfr_t v_bound;

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
   * \brief   The number of distributions to process.
   */
  uint32_t count;

  /*!
   * \brief   Paths to the distributions to process.
   */
  char ** paths;
} Estimate_Runs_Diagonal_Distribution_Arguments;


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
  Estimate_Runs_Diagonal_Distribution_Arguments * const arguments,
  const int argc,
  char ** argv)
{
  /* Initialize the arguments data structure. */
  arguments->delta_bound = BOUND_DELTA;
  arguments->eta_bound = 0;
  arguments->eta_bound_specified = FALSE;
  arguments->paths = NULL;
  arguments->count = 0;

  bool delta_bound_specified = FALSE;

  mpfr_init2(arguments->v_bound, PRECISION);
  mpfr_set_ui(arguments->v_bound, 0, MPFR_RNDN);

  /* Iterate over the command line arguments. */
  int i;

  for (i = 1; i < argc; i++) {
    if (0 == strcmp(argv[i], "-v-bound")) {
      if ((i + 1) >= argc) {
        fprintf(stderr, "Error: Expected <v-bound> to follow after "
          "-v-bound.\n");
        return FALSE;
      }

      if (mpfr_cmp_ui(arguments->v_bound, 0) != 0) {
        fprintf(stderr,
          "Error: The volume quotient bound cannot be twice specified.\n");
        return FALSE;
      }

      if (mpfr_set_str(arguments->v_bound, argv[i + 1], 10, MPFR_RNDN) != 0) {
        fprintf(stderr, "Error: Failed to parse <v-bound> after -v-bound.\n");
        return FALSE;
      }

      if (mpfr_cmp_ui(arguments->v_bound, 0) <= 0) {
        fprintf(stderr, "Error: The <v-bound> after -v-bound must be "
          "positive.\n");
        return FALSE;
      }

      i++;

      continue;
    }

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

    /* Stop parsing. */
    break;
  }

  /* Set default parameters if arguments were not explicitly specified. */
  if (mpfr_cmp_ui(arguments->v_bound, 0) == 0) {
    mpfr_set_ui(arguments->v_bound, 2, MPFR_RNDN);
  }

  /* Parse paths to distributions. */
  if ((argc - i) <= 0) {
    fprintf(stderr, "Error: Incorrect command line arguments; expected at "
      "least one <distribution>.\n");
    return FALSE;
  }

  arguments->count = (uint32_t)(argc - i);
  arguments->paths = &(argv[i]);

  /* Iterate over the paths. */
  for (uint32_t j = 0; j < (arguments->count); j++) {
    /* Check that the distribution exists. */
    if (0 != access(arguments->paths[j], F_OK)) {
      fprintf(stderr, "Error: The distribution \"%s\" does not exists.\n",
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
 * \brief   Prints the command line arguments.
 *
 * \param[in, out] file     Then file to which to print the arguments.
 * \param[in] arguments     The arguments data structure to print to the file.
 * \param[in] distribution  The distribution being processed.
 */
static void arguments_fprintf(
  FILE * const file,
  const Estimate_Runs_Diagonal_Distribution_Arguments * const arguments,
  const Diagonal_Distribution * const distribution)
{
  if (TRUE == arguments->eta_bound_specified) {
    mpfr_fprintf(file, "# Bounds: (eta = %u (%u), delta = %u, v = %Rg)\n",
      arguments->eta_bound,
      distribution->parameters.eta_bound,
      arguments->delta_bound,
      arguments->v_bound);
  } else {
    mpfr_fprintf(file, "# Bounds: (eta = <all> (%u), delta = %u, v = %Rg)\n",
      distribution->parameters.eta_bound,
      arguments->delta_bound,
      arguments->v_bound);
  }

  log_timestamp_fprintf(file);
}

/*!
 * \brief   Clears an initialized command line arguments data structure.
 *
 * \param[in, out] arguments   The arguments data structure to clear.
 */
static void arguments_clear(
  Estimate_Runs_Diagonal_Distribution_Arguments * const arguments)
{
  mpfr_clear(arguments->v_bound);

  memset(arguments, 0, sizeof(Estimate_Runs_Diagonal_Distribution_Arguments));
}

/*!
 * \}
 */

/*!
 * \brief   The main function on the client node.
 *
 * This function is called once by main() for each distribution to process.
 */
static void main_client()
{
  /* Receive broadcast of the distribution. */
  Diagonal_Distribution distribution;
  diagonal_distribution_init_bcast_recv(&distribution, MPI_RANK_ROOT);

  /* Receive broadcast the delta and eta bounds. */
  uint32_t bounds[2];

  if (MPI_SUCCESS !=
    MPI_Bcast(bounds, 2, MPI_UNSIGNED, MPI_RANK_ROOT, MPI_COMM_WORLD))
  {
    critical("main_server(): "
      "Failed to receive broadcast of the delta and eta bounds.");
  };

  const uint32_t delta_bound = bounds[0];
  const uint32_t eta_bound = bounds[1];

  /* Setup a random state. */
  Random_State random_state;
  random_init(&random_state);

  while (TRUE) {
    /* Setup the ordered lists. */
    Tau_Ordered_List tau_ordered_list;
    tau_ordered_list_init(
        &tau_ordered_list,
        TAU_CHUNK_COUNT * TAU_CHUNK_SIZE / 100); /* Keep top one percent. */

    /* Receive broadcast of n. */
    uint32_t n;

    if (MPI_SUCCESS != MPI_Bcast(
          &n,
          1, /* count */
          MPI_UNSIGNED,
          MPI_RANK_ROOT,
          MPI_COMM_WORLD))
    {
      critical("main_client(): Failed to receive broadcast of n.");
    }

    /* If n is set to zero, stop processing. */
    if (0 == n) {
      break;
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

      /* Receive a job from the server. */
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
        tau_ordered_list_send_merge(&tau_ordered_list, MPI_RANK_ROOT);

        break;
      } else if (MPI_JOB_ESTIMATE_TAU != job) {
        critical("main_client(): Unknown job (%u).", job);
      }

      /* Receive the chunk identifier. Use it setup deterministic randomness. */
      uint32_t chunk;

      if (MPI_SUCCESS != MPI_Recv(
          &chunk,
          1, /* count */
          MPI_UNSIGNED,
          MPI_RANK_ROOT,
          MPI_TAG_CHUNK,
          MPI_COMM_WORLD,
          &status))
      {
        critical("main_client(): Failed to receive MPI_TAG_CHUNK.");
      }

      for (uint32_t i = 0; i < TAU_CHUNK_SIZE; i++) {
        long double tau;

        tau_estimate_diagonal(
          &distribution,
          &random_state,
          n,
          delta_bound,
          eta_bound,
          tau);

        tau_ordered_list_insert(&tau_ordered_list, tau);
      }
    }

    /* Clear memory. */
    tau_ordered_list_clear(&tau_ordered_list);
  }

  /* Close the random state. */
  random_close(&random_state);
}

/*!
 * \brief   A function on the server node for stopping the client nodes.
 *
 * This function is called by main_server() one times for each distribution to
 * stop the client nodes from processing the distribution. In response to this
 * signal, the client nodes will either proceed to process the next
 * distribution, or shut down, if there is no more distribution to process.
 */
static void main_server_stop()
{
  /* Broadcast n. */
  uint32_t n = 0;

  if (MPI_SUCCESS != MPI_Bcast(
      &n,
      1, /* count */
      MPI_UNSIGNED,
      MPI_RANK_ROOT,
      MPI_COMM_WORLD))
  {
    critical("main_server_stop(): Failed to broadcast n.");
  }
}

/*!
 * \brief   A function on the server node for estimating a volume quotient.
 *
 * This function is called by main_server() one or more times for each
 * distribution to generate estimates for various n.
 *
 * Note that this function only issues jobs to client nodes where the actual
 * estimate is computed, for detail see main_client().
 *
 * \param[in, out] v       An arbitrary precision floating point value in which
 *                         to store the volume quotient computed.
 * \param[in] n            The value of n for which to compute the estimate.
 * \param[in] distribution The distribution to sample to compute the estimate.
 * \param[in] mpi_size     The number of nodes.
 *
 * \return  Returns #TRUE if the estimate was successfully computed. Otherwise,
 *          if the estimate computed is out of bounds, #FALSE is returned.
 */
static bool main_server_estimate_volume_quotient(
  mpfr_t v,
  uint32_t n,
  const Diagonal_Distribution * const distribution,
  const int mpi_size)
{
  /* Setup the ordered lists. */
  Tau_Ordered_List tau_ordered_list;
  tau_ordered_list_init(
      &tau_ordered_list,
      TAU_CHUNK_COUNT * TAU_CHUNK_SIZE / 100); /* Keep top one percent. */

  /* Broadcast n. */
  if (MPI_SUCCESS != MPI_Bcast(
        &n,
        1, /* count */
        MPI_UNSIGNED,
        MPI_RANK_ROOT,
        MPI_COMM_WORLD))
  {
    critical("main_server_estimate_volume_quotient(): Failed to broadcast n.");
  }

  printf("Processing n = %u...\n", n);

  uint32_t chunks_submitted = 0;

  uint32_t active_nodes = mpi_size - 1;

  while (active_nodes > 0) {
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
      critical("main_server_estimate_volume_quotient(): "
        "Failed to receive notification.");
    }

    if (MPI_NOTIFY_READY != notification) {
      critical("main_server_estimate_volume_quotient(): "
        "Unexpected notification (%u).", notification);
    }

    /* Send out a new job request or stop the node. */
    uint32_t job;

    if (chunks_submitted < TAU_CHUNK_COUNT) {
      /* Send out a new job request. */
      job = MPI_JOB_ESTIMATE_TAU;

      if (MPI_SUCCESS != MPI_Send(
          &job,
          1, /* count */
          MPI_UNSIGNED,
          status.MPI_SOURCE,
          MPI_TAG_JOB,
          MPI_COMM_WORLD))
      {
        critical("main_server_estimate_volume_quotient(): "
          "Failed to send job.");
      }

      if (MPI_SUCCESS != MPI_Send(
          &chunks_submitted,
          1, /* count */
          MPI_UNSIGNED,
          status.MPI_SOURCE,
          MPI_TAG_CHUNK,
          MPI_COMM_WORLD))
      {
        critical("main_server_estimate_volume_quotient(): "
          "Failed to send MPI_TAG_CHUNK.");
      }

      /* Increment the counter. */
      chunks_submitted++;

      printf("Processing %u / %u ...\n", chunks_submitted, TAU_CHUNK_COUNT);
    } else {
      /* Stop the node and collect its partial ordered lists of tau values. */
      job = MPI_JOB_STOP;

      if (MPI_SUCCESS != MPI_Send(
          &job,
          1, /* count */
          MPI_UNSIGNED,
          status.MPI_SOURCE,
          MPI_TAG_JOB,
          MPI_COMM_WORLD))
      {
        critical("main_server_estimate_volume_quotient(): "
          "Failed to send job.");
      }

      /* Fetch the lists from the node and merge with the server list. */
      tau_ordered_list_recv_merge(&tau_ordered_list, status.MPI_SOURCE);

      active_nodes--;

      printf("Stopping node %u / %u with %u nodes in progress...\n",
        status.MPI_SOURCE,
        mpi_size - 1,
        active_nodes);
    }
  }

  /* Open the logfile. */
  char log_path[MAX_SIZE_PATH_BUFFER];
  safe_snprintf(log_path, MAX_SIZE_PATH_BUFFER,
    "%s/estimate-runs-diagonal.txt", LOGS_DIRECTORY);

  FILE * logfile = fopen(log_path, "a+");
  if (NULL == logfile) {
    critical("main_server_estimate_volume_quotient(): "
      "Failed to open \"%s\" for appending.", log_path);
  }

  /* Extract constants from the distribution parameters. */
  const uint32_t m = distribution->parameters.m;
  const uint32_t sigma = distribution->parameters.sigma;
  const uint32_t l = distribution->parameters.l;
  const uint32_t s = distribution->parameters.s;

  /* Extract the tau value. */
  long double tau = tau_ordered_list.tau[0];

  /* Check if the values of tau returned is out of bounds. */
  if (DBL_MAX == tau) {
    /* The tau value returned was sampled out of bounds of the distribution so
     * we have no data from which to compute the quotients. */
     printf("Warning: The tau returned is out of bounds.\n");

    /* Log a note to this effect... */
    mpfr_fprintf(logfile, "m: %u sigma: %u %c: %u n: %u -- *** "
      "<-- out of bounds\n",
      m,
      sigma,
      (0 != s) ? 's' : 'l',
      (0 != s) ?  s  :  l,
      n);
    mpfr_printf("m: %u sigma: %u %c: %u n: %u -- *** <-- out of bounds\n",
      m,
      sigma,
      (0 != s) ? 's' : 'l',
      (0 != s) ?  s  :  l,
      n);

    fclose(logfile);
    logfile = NULL;

    return FALSE;
  }

  /* Check and print the number of values of tau that are out of bounds. */
  uint32_t i;

  for (i = 0; i < tau_ordered_list.size; i++) {
    if (DBL_MAX != tau_ordered_list.tau[tau_ordered_list.size - 1 - i]) {
      break;
    }
  }

  /* Compute the volume quotient. */
  tau_volume_quotient_diagonal(
    m,
    sigma,
    l,
    n,
    distribution->parameters.d,
    tau,
    v);

  /* Log the tau value and associated volume quotient. */
  mpfr_fprintf(logfile, "m: %u sigma: %u %c: %u n: %u -- "
    "tau: %Lf v: %RG <%u>\n",
    m,
    sigma,
    (0 != s) ? 's' : 'l',
    (0 != s) ?  s  :  l,
    n,
    tau, v, i);
  mpfr_printf("m: %u sigma: %u %c: %u n: %u -- "
    "tau: %Lf v: %RG <%u>\n",
    m,
    sigma,
    (0 != s) ? 's' : 'l',
    (0 != s) ?  s  :  l,
    n,
    tau, v, i);

  fclose(logfile);
  logfile = NULL;

  /* Clear memory. */
  tau_ordered_list_clear(&tau_ordered_list);

  /* Signal success. */
  return TRUE;
}

/*!
 * \brief   The main function on the server node.
 *
 * This function is called once by main() for each distribution to process.
 *
 * \param[in] arguments    The parsed command line arguments.
 * \param[in] distribution The distribution to sample to compute the estimate.
 * \param[in] path         The path to the distribution.
 * \param[in] mpi_size     The number of nodes.
 */
static void main_server(
  const Estimate_Runs_Diagonal_Distribution_Arguments * const arguments,
  const Diagonal_Distribution * const distribution,
  char * path,
  int mpi_size)
{
  /* Declare variables. */
  mpfr_t v_0;
  mpfr_init2(v_0, PRECISION);

  mpfr_t v_1;
  mpfr_init2(v_1, PRECISION);

  mpfr_t v_n;
  mpfr_init2(v_n, PRECISION);

  mpfr_t factor;
  mpfr_init2(factor, PRECISION);

  /* Broadcast the distribution. */
  diagonal_distribution_bcast_send(distribution, MPI_RANK_ROOT);

  /* Broadcast the delta and eta bounds. */
  uint32_t bounds[2];

  bounds[0] = arguments->delta_bound;

  if (TRUE == arguments->eta_bound_specified) {
    bounds[1] = arguments->eta_bound;
  } else {
    bounds[1] = distribution->parameters.eta_bound;
  }

  if (MPI_SUCCESS !=
    MPI_Bcast(bounds, 2, MPI_UNSIGNED, MPI_RANK_ROOT, MPI_COMM_WORLD))
  {
    critical("main_server(): Failed to broadcast the delta and eta bounds.");
  };

  /* Create the log directory if it does not exist. */
  if (0 != access(LOGS_DIRECTORY, F_OK)) {
    if (0 != mkdir(LOGS_DIRECTORY, DIRECTORY_PERMISSIONS)) {
      critical("main_server(): Failed to create the directory \"%s\".",
        LOGS_DIRECTORY);
    }
  }

  /* Write an entry in the logfile. */
  char log_path[MAX_SIZE_PATH_BUFFER];
  safe_snprintf(log_path, MAX_SIZE_PATH_BUFFER,
    "%s/estimate-runs-diagonal.txt", LOGS_DIRECTORY);

  FILE * logfile = fopen(log_path, "a+");
  if (NULL == logfile) {
    critical("main_server(): Failed to open \"%s\" for appending.", log_path);
  }

  fprintf(logfile, "\n# Processing: %s\n", truncate_path(path));
  arguments_fprintf(logfile, arguments, distribution);

  fclose(logfile);
  logfile = NULL;

  /* Extract s from the distributions parameters. */
  const uint32_t s = (distribution->parameters.s < 1) ?
    1 : distribution->parameters.s;

  uint32_t n;

  /* Compute the maximum volume quotient for n = s. */
  bool result;

  result = main_server_estimate_volume_quotient(v_0, s, distribution, mpi_size);
  if (FALSE == result) {
    main_server_stop();
    goto main_server_clear;
  }
  if (mpfr_cmp(v_0, arguments->v_bound) < 0) {
    main_server_stop();
    goto main_server_clear;
  }

  /* Compute the maximum volume quotient for n = s + 1. */
  result =
    main_server_estimate_volume_quotient(v_1, s + 1, distribution, mpi_size);
  if (FALSE == result) {
    main_server_stop();
    goto main_server_clear;
  }
  if (mpfr_cmp(v_1, arguments->v_bound) < 0) {
    main_server_stop();
    goto main_server_clear;
  }

  /* Interpolate and jump ahead. */
  mpfr_div(factor, v_0, v_1, MPFR_RNDN);
  mpfr_div(v_n, v_1, factor, MPFR_RNDN);

  for (n = s + 2 ; n <= MAX_N; n++) {
    if (mpfr_cmp(v_n, arguments->v_bound) < 0) {
      break;
    }

    mpfr_div(v_n, v_n, factor, MPFR_RNDN);
  }

  if (MAX_N == n) {
    critical("main_server(): Failed to interpolate in n.");
  }

  /* Compute for n. */
  result = main_server_estimate_volume_quotient(v_n, n, distribution, mpi_size);
  if (FALSE == result) {
    main_server_stop();
    goto main_server_clear;
  }

  if (mpfr_cmp(v_n, arguments->v_bound) < 0) {
    /* Decrease n by one until an unacceptable n is found. We require n >= s and
     * we have already computed v for n = s and n = s + 1. */

    for (n--; n > (s + 1); n--) {
      result =
        main_server_estimate_volume_quotient(v_n, n, distribution, mpi_size);
      if (FALSE == result) {
        break;
      }
      if (mpfr_cmp(v_n, arguments->v_bound) >= 0) {
        break;
      }
    }
  } else {
    /* Increase n by one until an acceptable n is found. */
    for (n++; n <= MAX_N; n++) {
      result =
        main_server_estimate_volume_quotient(v_n, n, distribution, mpi_size);
      if (FALSE == result) {
        break;
      }
      if (mpfr_cmp(v_n, arguments->v_bound) < 0) {
        break;
      }
    }

    if (MAX_N == n) {
      critical("main_server(): Maximum n reached.");
    }
  }

  main_server_stop();

main_server_clear:

  /* Clear memory. */
  mpfr_clear(v_0);
  mpfr_clear(v_1);
  mpfr_clear(v_n);

  mpfr_clear(factor);
}

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  fprintf(file, "Synopsis: mpirun estimate_runs_diagonal_distribution \\\n"
          "   [ -v-bound <v-bound> ] [ -delta-bound <delta-bound> ] \\\n"
          "      [ -eta-bound <eta-bound> ] \\\n"
          "         <distribution> { <distribution> }\n");

  fprintf(file, "\n");
  fprintf(file, "Volume quotient bound: -- defaults to 2\n");
  fprintf(file,
    " -v-bound      explicitly set the volume quotient bound to <v-bound>\n");

  fprintf(file, "\n");
  fprintf(file, "Delta bound: -- defaults to %u\n", BOUND_DELTA);
  fprintf(file,
    " -delta-bound  explicitly set the delta bound to <delta-bound>\n");

  fprintf(file, "\n");
  fprintf(file, "Eta bound: -- "
    "defaults to the eta bound for the distribution\n");
  fprintf(file,
    " -eta-bound    explicitly set the eta bound to <eta-bound>\n");
}

/*!
 * \brief   The main entry point to the estimate_runs_diagonal_distribution
 *          executable.
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

  Estimate_Runs_Diagonal_Distribution_Arguments arguments;

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

    Diagonal_Distribution_Loader loader;
    diagonal_distribution_loader_init(
      &loader, arguments.paths, arguments.count);

    for (uint32_t i = 0; i < arguments.count; i++) {
      Diagonal_Distribution * distribution =
        diagonal_distribution_loader_pop(&loader);

      printf("Processing: %s\n", arguments.paths[i]);

      if (TRUE == arguments.eta_bound_specified) {
        if (arguments.eta_bound > distribution->parameters.eta_bound) {
          critical("main(): The eta bound specified via -eta-bound is greater "
            "than the eta bound used to generate the distribution.");
        }
      }

      main_server(&arguments, distribution, arguments.paths[i], mpi_size);

      diagonal_distribution_clear(distribution);
    }

    diagonal_distribution_loader_clear(&loader);

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
