/*!
 * \file    parameters.h
 * \ingroup parameters
 *
 * \brief   The definition of data structures representing parameters for
 *          probability distributions, and the declarations of functions for
 *          manipulating such parameters.
 */

/*!
 * \defgroup parameters Parameters
 * \ingroup  distribution
 *
 * \brief    A module for parameters for probability distributions.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "random.h"
#include "common.h"

#include <gmp.h>
#include <mpfr.h>

#include <stdint.h>

/*!
 * \brief   A data structure that holds the parameters for probability
 *          distributions.
 * \ingroup parameters
 */
typedef struct {
  /*!
   * \brief The length m in bits of the logarithm d or order r.
   */
  uint32_t m;

  /*!
   * \brief The parameter l.
   *
   * If s is explicitly specified, the parameter l is set to ceil(m / s).
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
   * In the two-dimensional case, we set
   *
   *    min_alpha_d = min_alpha_r = max(m - t, 0)
   *
   * and
   *
   *    max_alpha_d = max_alpha_r = min(m + t - 1, m + l - 2).
   *
   * \remark Further restrictions may be imposed on which regions are actually
   *         computed and included in the distribution, especially in the
   *         two-dimensional case, to keep the error bound under control. For
   *         details, see the skip notifications in the MPI programs.
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
   * \brief The minimum logarithmic unsigned alpha_d value.
   *
   * The region min_alpha_d <= log_alpha_d <= max_alpha_d + 1 is considered
   * on either side of the alpha_d-axis.
   */
  uint32_t min_alpha_d;

  /*!
   * \brief The maximum logarithmic unsigned alpha_d value.
   *
   * The region min_alpha_d <= log_alpha_d <= max_alpha_d + 1 is considered
   * on either side of the alpha_d-axis.
   */
  uint32_t max_alpha_d;

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
} Parameters;

/*!
 * \name Initialization
 * \{
 */

/*!
 * \brief   Initializes the parameters.
 *
 * \param[in, out] parameters   The parameters to be initialized.
 */
void parameters_init(
  Parameters * const parameters);

/*!
 * \brief   Clears the parameters.
 *
 * \param[in, out] parameters   The parameters to be cleared.
 */
void parameters_clear(
  Parameters * const parameters);

/*!
 * \}
 */

/*!
 * \name Setup
 * \{
 */

/*!
 * \brief   Sets up the parameters for an explicitly specified d or r, and 
 *          explicitly specified values of m and s.
 *
 * This function sets l = ceil(m / s).
 *
 * \param[in, out] parameters   The parameters to setup.
 * \param[in] d                 The logarithm d.
 * \param[in] r                 The order r.
 * \param[in] m                 An upper bound on the length in bits of d or r.
 * \param[in] s                 The tradeoff factor s.
 * \param[in] t                 The parameter t, see Parameters::t for details.
 */
void parameters_explicit_m_s(
  Parameters * const parameters,
  const mpz_t d,
  const mpz_t r,
  const uint32_t m,
  const uint32_t s,
  const uint32_t t);

/*!
 * \brief   Sets up the parameters for an explicitly specified d or r, and 
 *          explicit specified values of m and l.
 *
 * This function sets s = 0.
 *
 * \param[in, out] parameters   The parameters to setup.
 * \param[in] d                 The logarithm d.
 * \param[in] r                 The order r.
 * \param[in] m                 An upper bound on the length in bits of d or r.
 * \param[in] l                 The parameter l.
 * \param[in] t                 The parameter t, see Parameters::t for details.
 */
void parameters_explicit_m_l(
  Parameters * const parameters,
  const mpz_t d,
  const mpz_t r,
  const uint32_t m,
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
void parameters_copy(
  Parameters * const dst_parameters,
  const Parameters * const src_parameters);

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
void parameters_import(
  Parameters * const parameters,
  FILE * const file);

/*!
 * \brief   Exports parameters to file.
 *
 * \param[in, out] parameters   The parameters to export.
 * \param[in, out] file         The file to which to write the parameters.
 */
void parameters_export(
  const Parameters * const parameters,
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
void parameters_bcast_send(
  const Parameters * const parameters,
  const int root);

/*!
 * \brief   Receives a broadcast of parameters.
 *
 * \param[in, out] parameters   The parameters.
 * \param[in] root              The rank of the node broadcasting the
 *                              parameters.
 */
void parameters_bcast_recv(
  Parameters * const parameters,
  const int root);

/*!
 * \}
 */

#endif /* PARAMETERS_H */
