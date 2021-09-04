/*!
 * \file    diagonal_parameters.h
 * \ingroup diagonal_parameters
 *
 * \brief   The definition of data structures representing parameters for
 *          diagonal probability distributions, and the declarations of
 *          functions for manipulating such parameters.
 */

/*!
 * \defgroup diagonal_parameters Diagonal Parameters
 * \ingroup  parameters
 *
 * \brief    A module for parameters for diagonal probability distributions.
 */

#ifndef DIAGONAL_PARAMETERS_H
#define DIAGONAL_PARAMETERS_H

#include "common.h"

#include <gmp.h>

#include <stdint.h>
#include <stdio.h>

/*!
 * \brief   A data structure that holds the parameters for diagonal probability
 *          distributions.
 * \ingroup diagonal_parameters
 */
typedef struct {
  /*!
   * \brief The length m in bits of the order r.
   */
  uint32_t m;

  /*!
   * \brief The padding length sigma.
   */
  uint32_t sigma;

  /*!
   * \brief The parameter l.
   *
   * If s is explicitly specified, the parameter l is set to
   * ceil((m + sigma) / s).
   */
  uint32_t l;

  /*!
   * \brief The tradeoff factor s.
   *
   * If l is explicitly specified, the tradeoff factor s is set to zero.
   */
  uint32_t s;

  /*!
   * \brief The parameter t.
   *
   * This parameter controls how large a part of the distribution is computed.
   *
   * We set
   *
   *    min_alpha_r = max(m - t, 0)
   *
   * and
   *
   *    max_alpha_r = min(m + t - 1, m + sigma - 2).
   *
   * \remark Further restrictions may be imposed on which regions are actually
   *         computed and included in the distribution. For details, see the
   *         skip notifications in the MPI programs.
   */
  uint32_t t;

  /*!
   * \brief The order r.
   */
  mpz_t r;

  /*!
   * \brief The logarithm d.
   */
  mpz_t d;

  /*!
   * \brief The minimum logarithmic unsigned alpha_r value.
   *
   * The region min_alpha_r <= log_alpha_r <= max_alpha_r + 1 is considered
   * on either side of the alpha_r-axis.
   */
  uint32_t min_alpha_r;

  /*!
   * \brief The maximum logarithmic unsigned alpha_r value.
   *
   * The region min_alpha_r <= log_alpha_r <= max_alpha_r + 1 is considered
   * on either side of the alpha_r-axis.
   */
  uint32_t max_alpha_r;
} Diagonal_Parameters;

/*!
 * \name Initialization
 * \{
 */

/*!
 * \brief   Initializes the parameters.
 *
 * \param[in, out] parameters   The parameters to be initialized.
 */
void diagonal_parameters_init(
  Diagonal_Parameters * const parameters);

/*!
 * \brief   Clears the parameters.
 *
 * \param[in, out] parameters   The parameters to be cleared.
 */
void diagonal_parameters_clear(
  Diagonal_Parameters * const parameters);

/*!
 * \}
 */

/*!
 * \name Setup
 * \{
 */

/*!
 * \brief   Sets up the parameters for explicitly specified d and r, and for
 *          explicitly specified values of m and s.
 *
 * This function sets l = ceil((m + sigma) / s).
 *
 * \param[in, out] parameters   The parameters to setup.
 * \param[in] d                 The logarithm d.
 * \param[in] r                 The order r.
 * \param[in] m                 The length in bits of r.
 * \param[in] sigma             The padding length sigma.
 * \param[in] s                 The tradeoff factor s.
 * \param[in] t                 The parameter t, see Parameters::t for details.
 */
void diagonal_parameters_explicit_m_s(
  Diagonal_Parameters * const parameters,
  const mpz_t d,
  const mpz_t r,
  const uint32_t m,
  const uint32_t sigma,
  const uint32_t s,
  const uint32_t t);

/*!
 * \brief   Sets up the parameters for explicitly specified d and r, and for
 *          explicitly specified values of m and l.
 *
 * This function sets s = 0.
 *
 * \param[in, out] parameters   The parameters to setup.
 * \param[in] d                 The logarithm d.
 * \param[in] r                 The order r.
 * \param[in] m                 The length in bits of r.
 * \param[in] sigma             The padding length sigma.
 * \param[in] l                 The parameter l.
 * \param[in] t                 The parameter t, see Parameters::t for details.
 */
void diagonal_parameters_explicit_m_l(
  Diagonal_Parameters * const parameters,
  const mpz_t d,
  const mpz_t r,
  const uint32_t m,
  const uint32_t sigma,
  const uint32_t l,
  const uint32_t t);

/*!
 * \}
 */

/*!
 * \name Copying
 * \{
 */

/*!
 * \brief   Copies parameters.
 *
 * \param[in, out] dst_parameters   The destination parameters.
 * \param[in] src_parameters        The source parameters.
 */
void diagonal_parameters_copy(
  Diagonal_Parameters * const dst_parameters,
  const Diagonal_Parameters * const src_parameters);

/*!
 * \}
 */

/*!
 * \name Importing and exporting
 * \{
 */

/*!
 * \brief   Imports parameters from file.
 *
 * \param[in, out] parameters   The parameters to which to import.
 * \param[in, out] file         The file from which to read the parameters.
 */
void diagonal_parameters_import(
  Diagonal_Parameters * const parameters,
  FILE * const file);

/*!
 * \brief   Exports parameters to file.
 *
 * \param[in, out] parameters   The parameters to export.
 * \param[in, out] file         The file to which to write the parameters.
 */
void diagonal_parameters_export(
  const Diagonal_Parameters * const parameters,
  FILE * const file);

/*!
 * \}
 */

/*!
 * \name Message passing
 * \{
 */

/*!
 * \brief   Sends a broadcast of parameters.
 *
 * \param[in, out] parameters   The parameters to send.
 * \param[in] root              The rank of the node broadcasting the
 *                              parameters.
 */
void diagonal_parameters_bcast_send(
  const Diagonal_Parameters * const parameters,
  const int root);

/*!
 * \brief   Receives a broadcast of parameters.
 *
 * \param[in, out] parameters   The parameters.
 * \param[in] root              The rank of the node broadcasting the
 *                              parameters.
 */
void diagonal_parameters_bcast_recv(
  Diagonal_Parameters * const parameters,
  const int root);

/*!
 * \}
 */

#endif /* DIAGONAL_PARAMETERS_H */
