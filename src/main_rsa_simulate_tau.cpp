/*!
 * \file    main_rsa_simulate_tau.cpp
 * \ingroup rsa_simulate_tau_exe
 *
 * \brief   The definition of the main entry point to the rsa_simulate_tau
 *          executable, and of associated functions.
 */

/*!
 * \defgroup rsa_simulate_tau_exe The rsa_simulate_tau executable
 * \ingroup  executable
 *
 * \brief    A module for the rsa_simulate_tau executable.
 */

#include "executables.h"

#include "common.h"
#include "errors.h"
#include "gmp_mpi.h"
#include "random.h"
#include "rsa.h"
#include "string_utilities.h"

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
 * \brief   The maximum tau for which to collect statistics.
 *
 * Statistics is collected for tau = 0, 1, .., #MAX_TAU.
 */
#define MAX_TAU 30

/*!
 * \brief   The maximum prime in the factor basis.
 */
#define MAX_PRIME (1ULL << 24)

/*!
 * \brief   The number of samples to collect per record.
 */
#define SAMPLES_PER_RECORD 10

/*!
 * \brief   The main function on the client node.
 *
 * This function is called once by main() on each client node.
 */
static void main_client()
{
  /* Allocate memory. */
  mpz_t p;
  mpz_init(p);

  mpz_t q;
  mpz_init(q);

  mpz_t N;
  mpz_init(N);

  mpz_t g;
  mpz_init(g);

  mpz_t r;
  mpz_init(r);

  mpz_t d;
  mpz_init(d);

  mpz_t tmp;
  mpz_init(tmp);

  uint32_t results[MAX_TAU + 1];

  Random_State random_state;
  random_init(&random_state);

  uint32_t notification;

  int result;

  /* Receive broadcast of the modulus length. */
  uint32_t modulus_length;

  if (MPI_SUCCESS != MPI_Bcast(
      &modulus_length, 1, MPI_UNSIGNED, MPI_RANK_ROOT, MPI_COMM_WORLD))
  {
    critical("main_client(): "
      "Failed to receive broadcast of the modulus length.");
  };

  uint32_t tmp2;

  if (MPI_SUCCESS != MPI_Bcast(
      &tmp2, 1, MPI_UNSIGNED, MPI_RANK_ROOT, MPI_COMM_WORLD))
  {
    critical("main_client(): "
      "Failed to broadcast of the check modulus size flag.");
  };

  const bool check_modulus_size = (1 == tmp2) ? TRUE : FALSE;

  /* Notify the server that we are ready. */
  notification = MPI_NOTIFY_READY;

  result = MPI_Send(
              &notification,
              1,
              MPI_UNSIGNED,
              MPI_RANK_ROOT,
              MPI_TAG_NOTIFY,
              MPI_COMM_WORLD);
  if (MPI_SUCCESS != result) {
    critical("main_client(): Failed to send notification.");
  }


  /* Generate all primes up to B using the sieve of Eratosthenes. */
  uint32_t * primes = (uint32_t *)malloc(MAX_PRIME * sizeof(uint32_t));
  if (NULL == primes) {
    critical("main_client(): Failed to allocate memory.");
  }

  uint32_t primes_count = 0;

  uint32_t * sieve = (uint32_t *)malloc(MAX_PRIME * sizeof(uint32_t));
  if (NULL == sieve) {
      critical("main_client(): Failed to allocate memory.");
  }

  for (uint32_t i = 2; i < MAX_PRIME; i++) {
    sieve[i] = i;
  }

  for (uint32_t i = 2; i < MAX_PRIME; i++) {
    if (0 == sieve[i]) {
        continue;
    }

    primes[primes_count++] = sieve[i];

    for (uint32_t j = i + sieve[i]; j < MAX_PRIME; j += sieve[i]) {
      sieve[j] = 0;
    }
  }

  free(sieve);
  sieve = NULL;

  while (TRUE) {
    /* Receive job. */
    uint32_t job;

    MPI_Status status;

    result = MPI_Recv(
                &job,
                1,
                MPI_UNSIGNED,
                MPI_RANK_ROOT,
                MPI_TAG_JOB,
                MPI_COMM_WORLD,
                &status);
    if (MPI_SUCCESS != result) {
      critical("main_client(): Failed to receive job.");
    }

    if (MPI_JOB_STOP == job) {
      break;
    } else if (MPI_JOB_SOLVE_SYSTEM != job) {
      critical("main_client(): Received an unknown job.");
    }

    for (uint32_t tau = 0; tau <= MAX_TAU; tau++) {
      results[tau] = 0;
    }

    for (uint32_t i = 0; i < SAMPLES_PER_RECORD; i++) {
      /* Generate an RSA key. */
      rsa_generate_modulus(
        p, q, modulus_length, check_modulus_size, &random_state);
      mpz_mul(N, p, q);

      /* Pick a generator. */
      while (TRUE) {
        random_generate_mpz(g, N, &random_state);
        mpz_gcd(d, g, N);
        if (0 == mpz_cmp_ui(d, 1)) {
          break;
        }
      }

      /* Compute the maximal order. */
      mpz_sub_ui(r, p, 1);
      mpz_sub_ui(tmp, q, 1);
      mpz_gcd(d, r, tmp);
      mpz_mul(r, r, tmp);
      mpz_div(r, r, d);

      /* Factor out small divisor in the order. */
      for (uint32_t j = 0; j < primes_count; ) {
        mpz_mod_ui(tmp, r, primes[j]);
        if (0 != mpz_cmp_ui(tmp, 0)) {
          j++; /* Process next prime. */

          continue;
        }

        mpz_div_ui(tmp, r, primes[j]);
        mpz_powm(tmp, g, tmp, N);

        if (0 == mpz_cmp_ui(tmp, 1)) {
          mpz_div_ui(r, r, primes[j]); /* Remove prime from r. */
        } else {
          j++; /* Process next prime. */
        }
      }

      /* Test the size of the order. */
      mpz_sub_ui(d, p, 1);
      mpz_sub_ui(tmp, q, 1);
      mpz_mul(tmp, tmp, d); /* tmp = phi(N) = (p - 1) * (q - 1) */
      mpz_div(tmp, tmp, r); /* tmp = phi(N) / r */

      for (uint32_t tau = 0; tau <= MAX_TAU; tau++) {
        /* r >= phi(N)/2^tau <=> 2^tau >= phi(N)/r <=> phi(N)/r <= 2^tau */
        if (mpz_cmp_ui(tmp, 1 << tau) <= 0) {
          results[tau] += 1;
        }
      }
    }

    /* Send the result to the server. */
    notification = MPI_NOTIFY_SOLUTION_DONE;

    result = MPI_Send(
                &notification,
                1,
                MPI_UNSIGNED,
                MPI_RANK_ROOT,
                MPI_TAG_NOTIFY,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("main_client(): Failed to send notification.");
    }

    result = MPI_Send(
                results,
                MAX_TAU + 1,
                MPI_UNSIGNED,
                MPI_RANK_ROOT,
                MPI_TAG_SOLUTION,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("main_client(): Failed to send results.");
    }
  }

  /* Clear memory. */
  random_close(&random_state);

  mpz_clear(p);
  mpz_clear(q);
  mpz_clear(N);
  mpz_clear(g);
  mpz_clear(r);
  mpz_clear(d);
  mpz_clear(tmp);

  free(primes);
  primes = NULL;
}

/*!
 * \brief   The main function on the server node.
 *
 * This function is called once by main().
 *
 * \param[in] modulus_length      The modulus length in bits.
 * \param[in] check_modulus_size  A flag that should be set to #TRUE if the
 *                                modulus must be of length exactly n bit, and
 *                                to #FALSE otherwise.
 * \param[in] records             The number of records to generate.
 * \param[in] mpi_size            The number of nodes.
 */
static void main_server(
  uint32_t modulus_length,
  const bool check_modulus_size,
  const uint32_t records,
  const int mpi_size)
{
  /* Allocate memory. */
  uint32_t results[MAX_TAU + 1];

  for (uint32_t tau = 0; tau <= MAX_TAU; tau++) {
    results[tau] = 0;
  }

  MPI_Status status;

  MPI_Status secondary_status;

  uint32_t count = 0;
  uint32_t issued = 0;

  /* Send broadcast of the modulus length. */
  if (MPI_SUCCESS != MPI_Bcast(
      &modulus_length, 1, MPI_UNSIGNED, MPI_RANK_ROOT, MPI_COMM_WORLD))
  {
    critical("main_server(): Failed to send broadcast of the modulus length.");
  };

  /* Send broadcast of the check modulus size flag. */
  uint32_t tmp = (TRUE == check_modulus_size) ? 1 : 0;

  if (MPI_SUCCESS != MPI_Bcast(
      &tmp, 1, MPI_UNSIGNED, MPI_RANK_ROOT, MPI_COMM_WORLD))
  {
    critical("main_server(): "
      "Failed to send broadcast of the check modulus size flag.");
  };

  /* Compute the number of client nodes. */
  if (mpi_size < 2) {
    /* Note: This check is unreachable. */
    critical("main_server(): At least two MPI nodes are required.");
  }

  uint32_t clients = (uint32_t)(mpi_size - 1);

  /* Create the log directory if it does not exist. */
  if (0 != access(LOGS_DIRECTORY, F_OK)) {
    if (0 != mkdir(LOGS_DIRECTORY, DIRECTORY_PERMISSIONS)) {
      critical("main_server(): Failed to create the directory \"%s\".",
        LOGS_DIRECTORY);
    }
  }

  /* Open the log file for appending, creating it if it does not exist. */
  char log_path[MAX_SIZE_PATH_BUFFER];
  safe_snprintf(log_path, MAX_SIZE_PATH_BUFFER,
    "%s/rsa-simulate-tau-%u.txt", LOGS_DIRECTORY, modulus_length);

  FILE * log_file = fopen(log_path, "a+");
  if (NULL == log_file) {
    critical("main_server(): Failed to open \"%s\" for appending.", log_path);
  }

  /* Write a header to the log file. */
  fprintf(log_file, "Modulus length: %u bits\n", modulus_length);
  fprintf(log_file, "Records: %u\n", records);
  fprintf(log_file, "Samples: %u\n\n", records * SAMPLES_PER_RECORD);

  /* Process the records. */
  while (TRUE) {
    /* Receive notification. */
    uint32_t notification;

    int result;

    result = MPI_Recv(
                &notification,
                1,
                MPI_UNSIGNED,
                MPI_ANY_SOURCE,
                MPI_TAG_NOTIFY,
                MPI_COMM_WORLD,
                &status);
    if (MPI_SUCCESS != result) {
      critical("main_server(): Failed to receive a notification.");
    }

    if (MPI_NOTIFY_SOLUTION_DONE == notification) {
      uint32_t client_results[MAX_TAU + 1];

      result = MPI_Recv(
                  client_results,
                  MAX_TAU + 1,
                  MPI_UNSIGNED,
                  status.MPI_SOURCE,
                  MPI_TAG_SOLUTION,
                  MPI_COMM_WORLD,
                  &secondary_status);
      if (MPI_SUCCESS != result) {
        critical("main_server(): Failed to receive results.");
      }

      for (uint32_t tau = 0; tau <= MAX_TAU; tau++) {
        results[tau] += client_results[tau];
      }

      /* Update the count. */
      count++;
    } else if (MPI_NOTIFY_READY != notification) {
      critical("main_server(): Received an unknown notification.");
    }

    if (count > 0) {
      printf("Total Count: %u / %u\n",
        count * SAMPLES_PER_RECORD, records * SAMPLES_PER_RECORD);

      for (uint32_t tau = 0; tau <= MAX_TAU; tau++) {
        printf(" Tau: %u Count: %u\n", tau, results[tau]);
      }

      printf("\n");
      fflush(stdout);
    }

    /* Send out a new job. */
    uint32_t job = MPI_JOB_SOLVE_SYSTEM;

    if (issued >= records) {
      job = MPI_JOB_STOP; /* Stop this node. */
      clients--;

      printf("Stopping node %u with %u node(s) remaining...\n",
        status.MPI_SOURCE, clients);
    }

    result = MPI_Send(
                &job,
                1,
                MPI_UNSIGNED,
                status.MPI_SOURCE,
                MPI_TAG_JOB,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("main_server(): Failed to send job.");
    }

    issued++;

    if (count == records) {
      break;
    }
  }

  /* Write final statistics to the log file. */
  for (uint32_t tau = 0; tau <= MAX_TAU; tau++) {
    fprintf(log_file, " Tau: %u Count: %u\n", tau, results[tau]);
  }
  fprintf(log_file, "\n");

  /* Stop all remaining nodes. */
  while (clients > 0) {
    /* Receive notification. */
    uint32_t notification;

    int result;

    result = MPI_Recv(
                &notification,
                1,
                MPI_UNSIGNED,
                MPI_ANY_SOURCE,
                MPI_TAG_NOTIFY,
                MPI_COMM_WORLD,
                &status);
    if (MPI_SUCCESS != result) {
      critical("main_server(): Failed to receive notification.");
    }

    if (MPI_NOTIFY_READY != notification) {
      critical("main_server(): Expected notification MPI_NOTIFY_READY (%u) but "
        "received  notification (%u).", MPI_NOTIFY_READY, notification);
    }

    /* Send out a stop request. */
    uint32_t job = MPI_JOB_STOP;

    result = MPI_Send(
                &job,
                1,
                MPI_UNSIGNED,
                status.MPI_SOURCE,
                MPI_TAG_JOB,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("main_server(): Failed to send job.");
    }

    clients--;

    printf("Stopping node %u with %u node(s) remaining...\n",
      status.MPI_SOURCE, clients);
  }

  /* Close the file. */
  fclose(log_file);
  log_file = NULL;
}

/*!
 * \brief   Parses the command line arguments.
 *
 * \param[in, out] modulus_length     The modulus length in bits.
 * \param[in, out] check_modulus_size A flag that is set to #TRUE if the modulus
 *                                    must be of length exactly n bit, and to
 *                                    #FALSE otherwise.
 * \param[in, out] records            The number of records.
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
 *           #FALSE otherwise. If #FALSE is returned, the output variables may
 *           be only partially initialized.
 */
static bool parse_command_line(
  uint32_t * const modulus_length,
  bool * const check_modulus_size,
  uint32_t * const records,
  const int argc,
  char ** argv)
{
  if ((3 != argc) && (4 != argc)) {
    fprintf(stderr, "Error: Incorrect command line arguments.\n");
    return FALSE;
  }

  uint32_t argv_offset = 1;

  if (4 == argc) {
    if (0 != strcmp("-check-modulus-size", argv[argv_offset++])) {
      fprintf(stderr, "Error: Expected the -check-modulus-size flag but found "
        "a different first command line argument..\n");
      return FALSE;
    }

    (*check_modulus_size) = TRUE;
  } else {
    (*check_modulus_size) = FALSE;
  }

  if (1 != sscanf(argv[argv_offset++], "%u", modulus_length)) {
    fprintf(stderr, "Error: Failed to parse <modulus_length>.\n");
    return FALSE;
  }

  if (((*modulus_length) <= 0) || (0 != ((*modulus_length) % 2))) {
    fprintf(stderr,
      "Error: The <modulus_length> must be a positive even integer.\n");
    return FALSE;
  }

  if (1 != sscanf(argv[argv_offset++], "%u", records)) {
    fprintf(stderr, "Error: Failed to parse <records>.\n");
    return FALSE;
  }

  if (0 == (*records)) {
    fprintf(stderr, "Error: <records> must be a non-zero positive integer.\n");
    return FALSE;
  }

  /* Signal success. */
  return TRUE;
}

/*!
 * \brief   Prints the command line synopsis.
 *
 * \param[in, out] file   The file to which to print the synopsis.
 */
static void print_synopsis(
  FILE * const file)
{
  fprintf(file, "Synopsis: mpirun rsa_simulate_tau\n");
  fprintf(file, "   [ -check-modulus-size ] <modulus_length> <records>\n");
  fprintf(file, "\n");
  fprintf(file, "Check modulus size: -- defaults to false\n");
  fprintf(file, " --check-modulus-size  "
    "Only accept moduli N < 2^n.\n");
  fprintf(file, "\n");
  fprintf(file, "The length n in bits of p, q in N = pq is "
    "<modulus_length>/2. You must \n");
  fprintf(file, "specify an even value for <modulus_length>. A total of "
    "%u * <records>\n", SAMPLES_PER_RECORD);
  fprintf(file, "random moduli N are generated when statistics is "
    "collected.\n");
}

/*!
 * \brief The main entry point to the rsa_simulate_tau executable.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \return Zero upon successful execution, non-zero otherwise.
 */
int main(int argc, char ** argv) {
  mpfr_set_default_prec(PRECISION);

  if (MPI_SUCCESS != MPI_Init(&argc, &argv))
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
  uint32_t result = TRUE;

  uint32_t modulus_length;
  bool check_modulus_size = FALSE; /* needed to avoid compiler warning */
  uint32_t records;

  if (MPI_RANK_ROOT == mpi_rank) {
    result =
      (uint32_t)parse_command_line(
        &modulus_length, &check_modulus_size, &records, argc, argv);
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

    printf("\n");

    main_server(modulus_length, check_modulus_size, records, mpi_size);
  } else {
    main_client();
  }

  if (MPI_SUCCESS != MPI_Finalize()) {
    critical("main(): Failed to finalize MPI.");
  }

  return 0;
}
