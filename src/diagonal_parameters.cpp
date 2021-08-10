/*!
 * \file    diagonal_parameters.cpp
 * \ingroup diagonal_parameters
 *
 * \brief   The definition of functions for manipulating data structures
 *          representing parameters for diagonal probability distributions.
 */

#include "diagonal_parameters.h"

#include "gmp_mpi.h"
#include "errors.h"

#include <gmp.h>
#include <mpfr.h>

#include <mpi.h>

#include <math.h>

#include <stdio.h>
#include <stdint.h>
#include <string.h>

/*!
 * \brief   Computes min_alpha_r and max_alpha_r as functions of m, sigma, and t
 *          as described in Diagonal_Parameters::t.
 *
 * \param[in, out] parameters   The parameters to setup.
 */
static void diagonal_parameters_setup_regions(
  Diagonal_Parameters * const parameters)
{
  /* Setup the region in alpha_r. */
  if (parameters->t > parameters->m) {
    parameters->min_alpha_r = 0;
  } else {
    parameters->min_alpha_r = parameters->m - parameters->t;
  }

  if (parameters->t >= parameters->sigma) {
    /* The maximum alpha is 2^(m+sigma-1) - 1. Hence this maximum is 
     * 2^(m+sigma-2), as we add up to one to the exponent during the search. */
    parameters->max_alpha_r = parameters->m + parameters->sigma - 2;
  } else {
    parameters->max_alpha_r = parameters->m + parameters->t - 1;
  }
}

void diagonal_parameters_init(
  Diagonal_Parameters * const parameters)
{
  memset(parameters, 0, sizeof(Diagonal_Parameters));

  mpz_init(parameters->r);
  mpz_init(parameters->d);
}

void diagonal_parameters_clear(
  Diagonal_Parameters * const parameters)
{
  mpz_clear(parameters->r);
  mpz_clear(parameters->d);
}

void diagonal_parameters_explicit_m_s(
  Diagonal_Parameters * const parameters,
  const mpz_t d,
  const mpz_t r,
  const uint32_t m,
  const uint32_t sigma,
  const uint32_t s,
  const uint32_t t)
{
  /* Store m, sigma and s, and compute l. */
  parameters->m = m;
  parameters->sigma = sigma;
  parameters->s = s;
  parameters->l = (uint32_t)ceil((double)(m + sigma) / (double)s);

  /* Store t. */
  parameters->t = t;

  /* Store d and r. */
  mpz_set(parameters->d, d);
  mpz_set(parameters->r, r);

  /* Setup regions. */
  diagonal_parameters_setup_regions(parameters);
}

void diagonal_parameters_explicit_m_l(
  Diagonal_Parameters * const parameters,
  const mpz_t d,
  const mpz_t r,
  const uint32_t m,
  const uint32_t sigma,
  const uint32_t l,
  const uint32_t t)
{
  /* Store m, sigma, s and l. */
  parameters->m = m;
  parameters->sigma = sigma;
  parameters->s = 0;
  parameters->l = l;

  /* Store t. */
  parameters->t = t;

  /* Store d and r. */
  mpz_set(parameters->d, d);
  mpz_set(parameters->r, r);

  /* Setup regions. */
  diagonal_parameters_setup_regions(parameters);
}

void diagonal_parameters_copy(
  Diagonal_Parameters * const dst,
  const Diagonal_Parameters * const src)
{
  /* Copy m, sigma, s and l. */
  dst->m = src->m;
  dst->sigma = src->sigma;
  dst->s = src->s;
  dst->l = src->l;

  /* Copy t. */
  dst->t = src->t;

  /* Copy the region. */
  dst->min_alpha_r = src->min_alpha_r;
  dst->max_alpha_r = src->max_alpha_r;

  /* Copy r and d. */
  mpz_set(dst->r, src->r);
  mpz_set(dst->d, src->d);
}

void diagonal_parameters_bcast_send(
  const Diagonal_Parameters * const parameters,
  const int root)
{
  /* Send broadcast of all integer parameters. */
  uint32_t data[7];

  /* Get m, sigma, s and l. */
  data[0] = parameters->m;
  data[1] = parameters->sigma;
  data[2] = parameters->s;
  data[3] = parameters->l;

  /* Get t. */
  data[4] = parameters->t;

  /* Get the region. */
  data[5] = parameters->min_alpha_r;
  data[6] = parameters->max_alpha_r;

  if (MPI_SUCCESS != MPI_Bcast(
    data,
    7, /* count */
    MPI_UNSIGNED,
    root,
    MPI_COMM_WORLD))
  {
    critical("diagonal_parameters_bcast_send(): "
      "Failed to send broadcast of parameters.");
  }

  /* Send broadcast of r and d. */
  mpz_bcast_send(parameters->r, root);
  mpz_bcast_send(parameters->d, root);
}

void diagonal_parameters_bcast_recv(
  Diagonal_Parameters * const parameters,
  const int root)
{
  /* Receive broadcast of all integer parameters. */
  uint32_t data[7];

  if (MPI_SUCCESS != MPI_Bcast(
    data,
    7, /* count */
    MPI_UNSIGNED,
    root,
    MPI_COMM_WORLD))
  {
    critical("diagonal_parameters_bcast_recv(): "
      "Failed to receive broadcast of parameters.");
  }

  /* Store m, sigma, s and l. */
  parameters->m = data[0];
  parameters->sigma = data[1];
  parameters->s = data[2];
  parameters->l = data[3];

  /* Store t. */
  parameters->t = data[4];

  /* Store the region. */
  parameters->min_alpha_r = data[5];
  parameters->max_alpha_r = data[6];

  /* Receive broadcast of r and d. */
  mpz_bcast_recv(parameters->r, root);
  mpz_bcast_recv(parameters->d, root);
}

void diagonal_parameters_export(
  const Diagonal_Parameters * const parameters,
  FILE * const file)
{
  /* Export m, sigma, s and l. */
  fprintf(file, "%u\n", parameters->m);
  fprintf(file, "%u\n", parameters->sigma);
  fprintf(file, "%u\n", parameters->s);
  fprintf(file, "%u\n", parameters->l);

  /* Export r and d. */
  gmp_fprintf(file, "%Zd\n", parameters->r);
  gmp_fprintf(file, "%Zd\n", parameters->d);

  /* Export t. */
  fprintf(file, "%u\n", parameters->t);

  /* Export the region. */
  fprintf(file, "%u\n", parameters->min_alpha_r);
  fprintf(file, "%u\n", parameters->max_alpha_r);
}

void diagonal_parameters_import(
  Diagonal_Parameters * const parameters,
  FILE * const file)
{
  /* Import m, sigma, s and l. */
  if (1 != fscanf(file, "%u\n", &(parameters->m))) {
    critical("diagonal_parameters_import(): Failed to import m.");
  }

  if (1 != fscanf(file, "%u\n", &(parameters->sigma))) {
    critical("diagonal_parameters_import(): Failed to import sigma.");
  }

  if (1 != fscanf(file, "%u\n", &(parameters->s))) {
    critical("diagonal_parameters_import(): Failed to import s.");
  }

  if (1 != fscanf(file, "%u\n", &(parameters->l))) {
    critical("diagonal_parameters_import(): Failed to import l.");
  }


  /* Import r and d. */
  if (1 != gmp_fscanf(file, "%Zd\n", parameters->r)) {
    critical("diagonal_parameters_import(): Failed to import r.");
  }

  if (1 != gmp_fscanf(file, "%Zd\n", parameters->d)) {
    critical("diagonal_parameters_import(): Failed to import d.");
  }


  /* Import t. */
  if (1 != fscanf(file, "%u\n", &(parameters->t))) {
    critical("diagonal_parameters_import(): Failed to import t.");
  }


  /* Import the region. */
  if (1 != fscanf(file, "%u\n", &(parameters->min_alpha_r))) {
    critical("diagonal_parameters_import(): Failed to import min_alpha_r.");
  }

  if (1 != fscanf(file, "%u\n", &(parameters->max_alpha_r))) {
    critical("diagonal_parameters_import(): Failed to import max_alpha_r.");
  }
}
