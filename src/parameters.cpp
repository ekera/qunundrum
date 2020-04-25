/*!
 * \file    parameters.cpp
 * \ingroup parameters
 *
 * \brief   The definition of functions for manipulating parameters for
 *          probability distributions.
 */

#include "parameters.h"

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
 * \brief   Computes min_alpha_d, max_alpha_d, min_alpha_r and max_alpha_r as
 *          functions of m, l, and t as described in Parameters::t.
 *
 * \param[in, out] parameters   The parameters to setup.
 */
static void parameters_setup_regions(
  Parameters * const parameters)
{
  /* Setup the region in alpha_d. */
  if (parameters->t > parameters->m) {
    parameters->min_alpha_d = 0;
  } else {
    parameters->min_alpha_d = parameters->m - parameters->t;
  }

  if (parameters->t >= parameters->l) {
    /* The maximum alpha is 2^(m+l-1) - 1. Hence this maximum is 2^(m+l-2), as
     * we add up to one to the exponent during the search. */
    parameters->max_alpha_d = parameters->m + parameters->l - 2;
  } else {
    parameters->max_alpha_d = parameters->m + parameters->t - 1;
  }

  /* Mirror the region in alpha_r. */
  parameters->min_alpha_r = parameters->min_alpha_d;
  parameters->max_alpha_r = parameters->max_alpha_d;
}

void parameters_init(
  Parameters * const parameters)
{
  memset(parameters, 0, sizeof(Parameters));

  mpz_init(parameters->r);
  mpz_init(parameters->d);
}

void parameters_clear(
  Parameters * const parameters)
{
  mpz_clear(parameters->r);
  mpz_clear(parameters->d);
}

void parameters_explicit_m_s(
  Parameters * const parameters,
  const mpz_t d,
  const mpz_t r,
  const uint32_t m,
  const uint32_t s,
  const uint32_t t)
{
  /* Store m and l in the parameters. */
  parameters->m = m;
  parameters->s = s;
  parameters->l = (uint32_t)ceil((double)m / (double)s);

  parameters->t = t;

  mpz_set(parameters->d, d);
  mpz_set(parameters->r, r);

  parameters_setup_regions(parameters);
}

void parameters_explicit_m_l(
  Parameters * const parameters,
  const mpz_t d,
  const mpz_t r,
  const uint32_t m,
  const uint32_t l,
  const uint32_t t)
{
  /* Store m and l in the parameters. */
  parameters->m = m;
  parameters->s = 0;
  parameters->l = l;

  parameters->t = t;

  mpz_set(parameters->d, d);
  mpz_set(parameters->r, r);

  parameters_setup_regions(parameters);
}

void parameters_copy(
  Parameters * const dst,
  const Parameters * const src)
{
  dst->m = src->m;
  dst->s = src->s;
  dst->l = src->l;

  dst->t = src->t;

  dst->min_alpha_d = src->min_alpha_d;
  dst->max_alpha_d = src->max_alpha_d;
  dst->min_alpha_r = src->min_alpha_r;
  dst->max_alpha_r = src->max_alpha_r;

  mpz_set(dst->r, src->r);
  mpz_set(dst->d, src->d);
}

void parameters_bcast_send(
  const Parameters * const parameters,
  const int root)
{
  uint32_t data[8];

  data[0] = parameters->m;
  data[1] = parameters->s;
  data[2] = parameters->l;

  data[3] = parameters->t;

  data[4] = parameters->min_alpha_d;
  data[5] = parameters->max_alpha_d;
  data[6] = parameters->min_alpha_r;
  data[7] = parameters->max_alpha_r;

  if (MPI_SUCCESS != MPI_Bcast(
    data,
    8, /* count */
    MPI_UNSIGNED,
    root,
    MPI_COMM_WORLD))
  {
    critical("parameters_bcast_send(): "
      "Failed to send broadcast of parameters.");
  }

  /* Send broadcast of r and d. */
  mpz_bcast_send(parameters->r, root);
  mpz_bcast_send(parameters->d, root);
}

void parameters_bcast_recv(
  Parameters * const parameters,
  const int root)
{
  uint32_t data[8];

  if (MPI_SUCCESS != MPI_Bcast(
    data,
    8, /* count */
    MPI_UNSIGNED,
    root,
    MPI_COMM_WORLD))
  {
    critical("parameters_bcast_recv(): "
      "Failed to receive broadcast of parameters.");
  }

  parameters->m = data[0];
  parameters->s = data[1];
  parameters->l = data[2];

  parameters->t = data[3];

  parameters->min_alpha_d = data[4];
  parameters->max_alpha_d = data[5];
  parameters->min_alpha_r = data[6];
  parameters->max_alpha_r = data[7];

  /* Receive broadcast of r and d. */
  mpz_bcast_recv(parameters->r, root);
  mpz_bcast_recv(parameters->d, root);
}

void parameters_export(
  const Parameters * const parameters,
  FILE * const file)
{
  fprintf(file, "%u\n", parameters->m);
  fprintf(file, "%u\n", parameters->s);
  fprintf(file, "%u\n", parameters->l);

  gmp_fprintf(file, "%Zd\n", parameters->r);
  gmp_fprintf(file, "%Zd\n", parameters->d);

  fprintf(file, "%u\n", parameters->t);

  fprintf(file, "%u\n", parameters->min_alpha_d);
  fprintf(file, "%u\n", parameters->max_alpha_d);
  fprintf(file, "%u\n", parameters->min_alpha_r);
  fprintf(file, "%u\n", parameters->max_alpha_r);
}

void parameters_import(
  Parameters * const parameters,
  FILE * const file)
{
  /* Import m, l and t. */
  if (1 != fscanf(file, "%u\n", &(parameters->m))) {
    critical("parameters_import(): Failed to import m.");
  }

  if (1 != fscanf(file, "%u\n", &(parameters->s))) {
    critical("parameters_import(): Failed to import s.");
  }

  if (1 != fscanf(file, "%u\n", &(parameters->l))) {
    critical("parameters_import(): Failed to import l.");
  }


  /* Import r, d and e. */
  if (1 != gmp_fscanf(file, "%Zd\n", parameters->r)) {
    critical("parameters_import(): Failed to import r.");
  }

  if (1 != gmp_fscanf(file, "%Zd\n", parameters->d)) {
    critical("parameters_import(): Failed to import d.");
  }


  /* Import minimum and maximum alpha_d and alpha_r values. */
  if (1 != fscanf(file, "%u\n", &(parameters->t))) {
    critical("parameters_import(): Failed to import t.");
  }

  if (1 != fscanf(file, "%u\n", &(parameters->min_alpha_d))) {
    critical("parameters_import(): Failed to import min_alpha_d.");
  }

  if (1 != fscanf(file, "%u\n", &(parameters->max_alpha_d))) {
    critical("parameters_import(): Failed to import max_alpha_d.");
  }

  if (1 != fscanf(file, "%u\n", &(parameters->min_alpha_r))) {
    critical("parameters_import(): Failed to import min_alpha_r.");
  }

  if (1 != fscanf(file, "%u\n", &(parameters->max_alpha_r))) {
    critical("parameters_import(): Failed to import max_alpha_r.");
  }
}
