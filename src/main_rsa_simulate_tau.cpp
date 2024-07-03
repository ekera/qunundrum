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
#define MAX_TAU 41

/*!
 * \brief   The maximum delta for which to collect statistics.
 *
 * Statistics is collected for delta = 0, 1, .., #MAX_DELTA.
 */
#define MAX_DELTA 41

/*!
 * \brief   The number of samples to collect per record.
 */
#define SAMPLES_PER_RECORD 100

/*!
 * \brief   The maximum tau for which a lower bound on the success probablity
 *          is available, giving raise to an expected lower bound.
 */
#define MAX_TAU_EXPECTED 22

/*!
 * \brief   Lower bounds on the success probability from tau equal to zero up to
 *          and including #MAX_TAU_EXPECTED.
 */
static long double lower_bounded_probability[MAX_TAU_EXPECTED + 1] =
{
  0, /* tau = 0 */
  0, /* tau = 1 */
  0, /* tau = 2 */
  0.11324380721014293063, /* tau =  3 */
  0.46634418467301933811, /* tau =  4 */
  0.67263910856031188914, /* tau =  5 */
  0.80112799120765551625, /* tau =  6 */
  0.88037408124863013370, /* tau =  7 */
  0.92882559994986759249, /* tau =  8 */
  0.95817480365440429994, /* tau =  9 */
  0.97567972465820551214, /* tau = 10 */
  0.98598795065912947759, /* tau = 11 */
  0.99200101938974547411, /* tau = 12 */
  0.99547227210972206611, /* tau = 13 */
  0.99745366269829229602, /* tau = 14 */
  0.99857655881236983978, /* tau = 15 */
  0.99920861704771607966, /* tau = 16 */
  0.99956226369547748292, /* tau = 17 */
  0.99975901422875727562, /* tau = 18 */
  0.99986789588772869315, /* tau = 19 */
  0.99992786545942676673, /* tau = 20 */
  0.99996075305442676241, /* tau = 21 */
  0.99997871733359514075  /* tau = 22 */
};

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

  mpz_t r;
  mpz_init(r);

  mpz_t r_mod_p;
  mpz_init(r_mod_p);

  mpz_t r_mod_q;
  mpz_init(r_mod_q);

  mpz_t d;
  mpz_init(d);

  mpz_t phi;
  mpz_init(phi);

  mpz_t tmp;
  mpz_init(tmp);

  mpz_t tmp2;
  mpz_init(tmp2);

  mpq_t tmp_q;
  mpq_init(tmp_q);

  mpq_t tmp2_q;
  mpq_init(tmp2_q);

  uint32_t results_tau[MAX_TAU + 1];
  uint32_t results_delta[MAX_DELTA + 1];

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

  uint32_t tmp_ui32;

  if (MPI_SUCCESS != MPI_Bcast(
      &tmp_ui32, 1, MPI_UNSIGNED, MPI_RANK_ROOT, MPI_COMM_WORLD))
  {
    critical("main_client(): "
      "Failed to broadcast the check modulus size flag.");
  };

  const bool check_modulus_size = (1 == tmp_ui32) ? TRUE : FALSE;

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
    } else if (MPI_JOB_TABULATE_BOUND != job) {
      critical("main_client(): Received an unknown job.");
    }

    /* Inititalize results tables. */
    for (uint32_t tau = 0; tau <= MAX_TAU; tau++) {
      results_tau[tau] = 0;
    }

    for (uint32_t delta = 0; delta <= MAX_DELTA; delta++) {
      results_delta[delta] = 0;
    }

    for (uint32_t i = 0; i < SAMPLES_PER_RECORD; i++) {
      /* Generate an RSA key. */
      rsa_generate_modulus(
        p,
        q,
        modulus_length,
        check_modulus_size,
        &random_state);

      /* Compute the order of g modulo p. */
      mpz_sub_ui(r_mod_p, p, 1);
      mpz_set(phi, r_mod_p);
      random_generate_mpz(tmp, r_mod_p, &random_state);
      mpz_gcd(d, tmp, r_mod_p);
      mpz_div(r_mod_p, r_mod_p, d);

      /* Compute the order of g modulo q. */
      mpz_sub_ui(r_mod_q, q, 1);
      mpz_mul(phi, phi, r_mod_q);
      random_generate_mpz(tmp, r_mod_q, &random_state);
      mpz_gcd(d, tmp, r_mod_q);
      mpz_div(r_mod_q, r_mod_q, d);

      /* Compute the order of g. */
      mpz_mul(r, r_mod_p, r_mod_q);
      mpz_gcd(d, r_mod_p, r_mod_q);
      mpz_div(r, r, d); /* r = lcm(r_mod_p, r_mod_q) */

      /* Test if r >= phi(N) / 2^tau for tau = 0, 1, .., MAX_TAU. */
      for (uint32_t tau = 0; tau <= MAX_TAU; tau++) {
        mpq_set_z(tmp_q, phi); /* tmp_q = phi(N) */

        mpz_set_ui(tmp, 0); /* tmp = 0 */
        mpz_setbit(tmp, tau); /* tmp = 2^tau */
        mpq_set_z(tmp2_q, tmp); /* tmp2_q = 2^tau */

        mpq_div(tmp_q, tmp_q, tmp2_q); /* tmp_q = phi(N) / 2^tau */

        /* Check if r >= phi(N) / 2^tau. */
        mpq_set_z(tmp2_q, r); /* tmp2_q = r */
        if (mpq_cmp(tmp2_q, tmp_q) >= 0) {
          results_tau[tau] += 1;
        }
      }

      /* Test if r >= 2^(m + ell) + (2^ell - 1) d for ell = m - delta where
       * delta = 0, 1, .., MAX_DELTA, and for m = l - 1 where l = n / 2. */
      const uint32_t l = modulus_length / 2;
      const uint32_t m = l - 1;

      mpz_add(d, p, q); /* d = p + q */
      mpz_sub_ui(d, d, 2); /* d = p + q - 2 */
      mpz_div_ui(d, d, 2); /* d = (p + q - 2) / 2 */
      mpz_set_ui(tmp, 0); /* tmp = 0 */
      mpz_setbit(tmp, l - 1); /* tmp = 2^(l - 1) */
      mpz_sub(d, d, tmp); /* d = (p + q - 2) / 2 - 2^(l - 1) */

      for (uint32_t delta = 0; delta <= MAX_DELTA; delta++) {
        const uint32_t ell = m - delta;

        mpz_set_ui(tmp, 0); /* tmp = 0 */
        mpz_setbit(tmp, ell); /* tmp = 2^ell */
        mpz_sub_ui(tmp, tmp, 1); /* tmp = 2^ell - 1 */
        mpz_mul(tmp, tmp, d); /* tmp = (2^ell - 1) d */

        mpz_set_ui(tmp2, 0); /* tmp2 = 0 */
        mpz_setbit(tmp2, m + ell); /* tmp2 = 2^(m + ell) */
        mpz_add(tmp, tmp2, tmp); /* tmp = 2^(m + ell) + (2^ell - 1) d */

        /* Check if r >= 2^(m + ell) + (2^ell - 1) d. */
        if (mpz_cmp(r, tmp) >= 0) {
          results_delta[delta] += 1;
        }
      }

    }

    /* Send the result to the server. */
    notification = MPI_NOTIFY_JOB_DONE;

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
                results_tau,
                MAX_TAU + 1,
                MPI_UNSIGNED,
                MPI_RANK_ROOT,
                MPI_TAG_SOLUTION_TAU,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("main_client(): Failed to send results for tau.");
    }

    result = MPI_Send(
                results_delta,
                MAX_DELTA + 1,
                MPI_UNSIGNED,
                MPI_RANK_ROOT,
                MPI_TAG_SOLUTION_DELTA,
                MPI_COMM_WORLD);
    if (MPI_SUCCESS != result) {
      critical("main_client(): Failed to send results for delta.");
    }
  }

  /* Clear memory. */
  random_close(&random_state);

  mpz_clear(p);
  mpz_clear(q);
  mpz_clear(r);
  mpz_clear(r_mod_p);
  mpz_clear(r_mod_q);
  mpz_clear(d);
  mpz_clear(phi);

  mpz_clear(tmp);
  mpz_clear(tmp2);

  mpq_clear(tmp_q);
  mpq_clear(tmp2_q);
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
  uint32_t results_tau[MAX_TAU + 1];

  for (uint32_t tau = 0; tau <= MAX_TAU; tau++) {
    results_tau[tau] = 0;
  }

  uint32_t results_delta[MAX_DELTA + 1];

  for (uint32_t delta = 0; delta <= MAX_DELTA; delta++) {
    results_delta[delta] = 0;
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
    "%s/rsa-simulate-tau.txt", LOGS_DIRECTORY);

  FILE * log_file = fopen(log_path, "a+");
  if (NULL == log_file) {
    critical("main_server(): Failed to open \"%s\" for appending.", log_path);
  }

  /* Write a header to the log file. */
  fprintf(log_file, "Modulus length: %u bits\n", modulus_length);
  fprintf(log_file, "Check modulus size: %s\n",
    (check_modulus_size) ? "True" : "False");
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

    if (MPI_NOTIFY_JOB_DONE == notification) {
      uint32_t client_results_tau[MAX_TAU + 1];

      result = MPI_Recv(
                  client_results_tau,
                  MAX_TAU + 1,
                  MPI_UNSIGNED,
                  status.MPI_SOURCE,
                  MPI_TAG_SOLUTION_TAU,
                  MPI_COMM_WORLD,
                  &secondary_status);
      if (MPI_SUCCESS != result) {
        critical("main_server(): Failed to receive results.");
      }

      for (uint32_t tau = 0; tau <= MAX_TAU; tau++) {
        results_tau[tau] += client_results_tau[tau];
      }

      uint32_t client_results_delta[MAX_DELTA + 1];

      result = MPI_Recv(
                  client_results_delta,
                  MAX_DELTA + 1,
                  MPI_UNSIGNED,
                  status.MPI_SOURCE,
                  MPI_TAG_SOLUTION_DELTA,
                  MPI_COMM_WORLD,
                  &secondary_status);
      if (MPI_SUCCESS != result) {
        critical("main_server(): Failed to receive results.");
      }

      for (uint32_t delta = 0; delta <= MAX_DELTA; delta++) {
        results_delta[delta] += client_results_delta[delta];
      }

      /* Update the count. */
      count++;
    } else if (MPI_NOTIFY_READY != notification) {
      critical("main_server(): Received an unknown notification.");
    }

    if (count > 0) {
      printf("Total Count: %u / %u\n",
        count * SAMPLES_PER_RECORD, records * SAMPLES_PER_RECORD);

      /* Tabulate in tau. */
      for (uint32_t tau = 0; tau <= MAX_TAU; tau++) {
        const long double probability =
          (long double)results_tau[tau] /
              (long double)(SAMPLES_PER_RECORD * count);

        char marker = ' ';

        if (tau <= MAX_TAU_EXPECTED) {
          if (probability < lower_bounded_probability[tau]) {
            marker = '!';
          } else {
            marker = 'o';
          }
        }

        printf(" Tau: %u Count: %u (%.16Lf) %c\n",
          tau, results_tau[tau], probability, marker);
      }

      printf("\n");
      fflush(stdout);

      /* Tabulate in delta. */
      for (uint32_t delta = 0; delta <= MAX_DELTA; delta++) {
        const long double probability =
          (long double)results_delta[delta] /
              (long double)(SAMPLES_PER_RECORD * count);

        printf(" Delta: %u Count: %u (%.16Lf)\n",
          delta, results_delta[delta], probability);
      }

      printf("\n");
      fflush(stdout);
    }

    /* Send out a new job. */
    uint32_t job = MPI_JOB_TABULATE_BOUND;

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

  /* Tabulate in tau. */
  for (uint32_t tau = 0; tau <= MAX_TAU; tau++) {
    const long double probability =
      (long double)results_tau[tau] /
        (long double)(SAMPLES_PER_RECORD * count);

    char marker = ' ';

    if (tau <= MAX_TAU_EXPECTED) {
      if (probability < lower_bounded_probability[tau]) {
        marker = '!';
      } else {
        marker = 'o';
      }
    }

    fprintf(log_file, " Tau: %u Count: %u (%.16Lf) %c\n",
          tau, results_tau[tau], probability, marker);
  }

  fprintf(log_file, "\n");

  /* Tabulate in delta. */
  for (uint32_t delta = 0; delta <= MAX_DELTA; delta++) {
    const long double probability =
      (long double)results_delta[delta] /
          (long double)(SAMPLES_PER_RECORD * count);

    fprintf(log_file, " Delta: %u Count: %u (%.16Lf)\n",
      delta, results_delta[delta], probability);
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
        "a different first command line argument.\n");
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

  if (((*modulus_length) < 128) || (0 != ((*modulus_length) % 2))) {
    fprintf(stderr,
      "Error: The <modulus_length> must be a positive even integer >= 128.\n");
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
  fprintf(file, "Synopsis: mpirun rsa_simulate_tau \\\n");
  fprintf(file, "   [ -check-modulus-size ] <modulus_length> <records>\n");
  fprintf(file, "\n");
  fprintf(file, "Check modulus size: -- defaults to false\n");
  fprintf(file, " -check-modulus-size   "
    "Only accept moduli N < 2^n.\n");
  fprintf(file, "\n");
  fprintf(file, "The length n in bits of p, q in N = pq is "
    "<modulus_length>/2. Specify \n");
  fprintf(file, "an even <modulus_length> that is >= 128. A total of "
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
