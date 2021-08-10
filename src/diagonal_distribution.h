/*!
 * \file    diagonal_distribution.h
 * \ingroup diagonal_distribution
 *
 * \brief   The declaration of data structures representing diagonal probability
 *          distributions, and of functions for manipulating such distributions.
 *
 * Diagonal probability distributions are used to represent the probability
 * distribution induced by Shor's algorithm for computing general discrete
 * logarithms [1] when modified as in [2] to start with a uniform superposition.
 * 
 * [1] Shor, P.W.: Polynomial-time algorithms for prime factorization and 
 * discrete logarithms on a quantum computer. In: SIAM Journal on Scientific 
 * Computing (SISC), volume 26(5), pp. 1484 (1997). 
 * 
 * [2] Ekerå, M.: Revisiting Shor's quantum algorithm for computing general 
 * discrete logarithms. In: ArXiv Pre-Print 1905.09084v2.
 */

/*!
 * \defgroup  diagonal_distribution Diagonal probability distributions
 * \ingroup   distribution
 *
 * \brief     A module for diagonal probability distributions.
 *
 * Diagonal probability distributions are used to represent the probability
 * distribution induced by Shor's algorithm for computing general discrete
 * logarithms [1] when modified as in [2] to start with a uniform superposition.
 * 
 * [1] Shor, P.W.: Polynomial-time algorithms for prime factorization and 
 * discrete logarithms on a quantum computer. In: SIAM Journal on Scientific 
 * Computing (SISC), volume 26(5), pp. 1484 (1997). 
 * 
 * [2] Ekerå, M.: Revisiting Shor's quantum algorithm for computing general 
 * discrete logarithms. In: ArXiv Pre-Print 1905.09084v2.
 */

#ifndef DIAGONAL_DISTRIBUTION_H
#define DIAGONAL_DISTRIBUTION_H

#include "diagonal_distribution_slice.h"
#include "random.h"
#include "diagonal_parameters.h"
#include "common.h"

#include <mpfr.h>

#include <stdio.h>
#include <stdint.h>

/*!
 * \brief   An enumeration of flags indicating how the logarithm d or order r
 *          were selected.
 */
enum {
  /*!
   * \brief   A flag indicating that the distribution was constructed by
   *          deterministically setting d or r from Catalan's constant.
   *
   * For more information, see parameters_deterministic().
   */
  DIAGONAL_DISTRIBUTION_FLAG_DETERMINISTIC = 0,

  /*!
   * \brief   A flag indicating that the distribution was constructed by
   *          setting a minimal r = 2^(m - 1) + 1 and random d.
   * 
   * \note    This flag is currently not used.
   */
  DIAGONAL_DISTRIBUTION_FLAG_MINIMAL = 4,

  /*!
   * \brief   A flag indicating that the distribution was constructed by
   *          setting a maximal r = 2^m - 1 and random d.
   * 
   * \note    This flag is currently not used.
   */
  DIAGONAL_DISTRIBUTION_FLAG_MAXIMAL = 8,

  /*!
   * \brief   A flag indicating that the distribution was constructed by
   *          randomly selecting d or r on (2^(m-1), 2^m).
   *
   * For more information, see parameters_selection_random_d_or_r().
   */
  DIAGONAL_DISTRIBUTION_FLAG_RANDOM = 16,

  /*!
   * \brief   A flag indicating that the distribution was constructed by
   *          explicitly setting d or r to a value provided by the caller.
   */
  DIAGONAL_DISTRIBUTION_FLAG_EXPLICIT = 32
};

/*!
 * \brief   A data structure representing a diagonal probability distribution.
 *
 * \ingroup diagonal_distribution
 */
typedef struct {
  /*! \brief  A vector of pointers to slices. */
  Diagonal_Distribution_Slice ** slices;

  /*! \brief  The total number of pointers to slices allocated. */
  uint32_t capacity;

  /*! \brief  The number of slices inserted into the distribution. */
  uint32_t count;

  /*! \brief  The distribution parameters. */
  Diagonal_Parameters parameters;

  /*! \brief  The total probability summed over all slices. */
  long double total_probability;

  /*! \brief  The total error summed over all slices. */
  long double total_error;

  /*! \brief   The distribution flags. */
  uint32_t flags;

  /*!
   * \brief   The precision using which the distribution was computed.
   *
   * \remark  Incorrect precisions may be reported for distributions imported
   *          and exported using binaries compiled to use different precision
   *          levels, especially if the distribution is modified by computing
   *          or recomputing slices inbetween import and export operations. */
  uint32_t precision;
} Diagonal_Distribution;

/*!
 * \name Memory allocation
 * \{
 */

/*!
 * \brief   Allocates memory for a distribution.
 *
 * \return  The pointer to the memory allocated for the distribution.
 */
Diagonal_Distribution * diagonal_distribution_alloc();

/*!
 * \brief   Deallocates the memory previously allocated for a distribution.
 *
 * \param[in, out] distribution   A pointer to a pointer to the memory allocated
 *                                for the distribution. This function assigns
 *                                NULL to the pointer as it deallocates the
 *                                memory to which it refers.
 */
void diagonal_distribution_dealloc(
  Diagonal_Distribution ** const distribution);

/*!
 * \}
 */

/*!
 * \name Initialization
 * \{
 */

/*!
 * \brief   Initializes a distribution with a given set of parameters.
 *
 * This function allocates memory for the slices in the distribution.
 *
 * \param[in, out] distribution   The distribution to initialize.
 * \param[in] parameters          The parameters to use to initialize the
 *                                distribution.
 * \param[in] flags               The distribution flags.
 * \param[in] capacity            The distribution capacity in slices.
 *
 * \remark  The parameter data structure is copied into the distribution data
 *          structure, rather than being referenced.
 */
void diagonal_distribution_init(
  Diagonal_Distribution * const distribution,
  const Diagonal_Parameters * const parameters,
  const uint32_t flags,
  const uint32_t capacity);

/*!
 * \brief   Initializes a distribution by reading from an already initialized
 *          source distribution and scaling down its dimension.
 *
 * This function allocates memory for the slices in the distribution.
 *
 * \param[in, out] dst_distribution   The distribution to initialize.
 * \param[in] src_distribution        The source distribution to be scaled.
 * \param[in] dimension               The dimension to be used for the
 *                                    destination distribution. Must divide the
 *                                    dimension of the source distribution.
 *
 * \remark  Used by the plotting functions to scale down the resolution in the
 *          probability distributions.
 */
void diagonal_distribution_init_copy_scale(
  Diagonal_Distribution * const dst_distribution,
  const Diagonal_Distribution * const src_distribution,
  const uint32_t dimension);

/*!
 * \brief   Clears a distribution.
 *
 * This function deallocates memory for slices in the distribution.
 *
 * \param[in, out] distribution   The distribution to clear.
 */
void diagonal_distribution_clear(
  Diagonal_Distribution * const distribution);

/*!
 * \}
 */

/*!
 * \name Importing and exporting
 * \{
 */

/*!
 * \brief   Initializes a distribution by importing it from a file.
 *
 * \param[in, out] distribution   The distribution to initialize.
 * \param[in, out] file           The file from which to read the distribution.
 */
void diagonal_distribution_init_import(
  Diagonal_Distribution * const distribution,
  FILE * const file);

/*!
 * \brief   Exports a probability distribution to file.
 *
 * \param[in] distribution        The distribution to export to file.
 * \param[in, out] file           The file to which to export the distribution.
 */
void diagonal_distribution_export(
  const Diagonal_Distribution * const distribution,
  FILE * const file);

/*!
 * \brief   Exports information about a probability distribution to file.
 *
 * \param[in] distribution        The distribution to export to file.
 * \param[in, out] file           The file to which to export information.
 * \param[in] export_slices       Export verbose information on slices.
 * \param[in] export_regions      Export verbose information on regions.
 */
void diagonal_distribution_export_info(
  const Diagonal_Distribution * const distribution,
  FILE * const file,
  const bool export_slices = TRUE,
  const bool export_regions = TRUE);

/*!
 * \}
 */

/*!
 * \name Message passing
 * \{
 */

/*!
 * \brief   Initializes a distribution by receiving a broadcast from a node.
 *
 * \param[in, out] distribution   The distribution to initialize.
 * \param[in] root                The rank of the node from which to receive.
 */
void diagonal_distribution_init_bcast_recv(
  Diagonal_Distribution * const distribution,
  const int root);

/*!
 * \brief   Broadcasts a distribution to all other nodes.
 *
 * \param[in] distribution        The distribution to broadcast.
 * \param[in] root                The rank of the root node.
 */
void diagonal_distribution_bcast_send(
  const Diagonal_Distribution * const distribution,
  const int root);

/*!
 * \}
 */

/*!
 * \name Managing slices
 * \{
 */

/*!
 * \brief   Inserts a slice into the distribution.
 *
 * \param[in, out] distribution   The distribution in which to insert the slice.
 * \param[in, out] slice          The slice to insert into the distribution.
 *
 * \remark  The slice data structure is referenced from the distribution data
 *          structure, rather than being copied into the structure. You may
 *          not modify or deallocate the slice after calling this function.
 */
void diagonal_distribution_insert_slice(
  Diagonal_Distribution * const distribution,
  Diagonal_Distribution_Slice * const slice);

/*!
 * \brief   Sorts the slices in the distribution in decreasing order with
 *          respect to the total probability of observing the slice.
 *
 * \remark  This function should be called prior to sampling the distribution.
 *          As the distribution slices are traversed in order when sampling,
 *          sampling a distribution in which the slices are sorted in decreasing
 *          order with respect to the total probability of observing the slice
 *          is much faster.
 *
 * \param[in, out] distribution   The distribution in which to sort the slices.
 */
void diagonal_distribution_sort_slices(
  Diagonal_Distribution * const distribution);

/*!
 * \}
 */

/*!
 * \name Sampling
 * \{
 */

/*!
 * \brief   Samples a slice from the distribution.
 *
 * \param[in] distribution        The distribution to sample.
 * \param[in, out] random_state   The random state to use when sampling.
 *
 * \return    A pointer to the slice sampled, or NULL if the slice sampled
 *            is outside the range of the distribution.
 */
const Diagonal_Distribution_Slice * diagonal_distribution_sample_slice(
  const Diagonal_Distribution * const distribution,
  Random_State * const random_state);

/*!
 * \brief   Samples a region on the alpha axis from the probability distribution
 *          and returns its bounds.
 *
 * \remark  The region bounds are represented as signed logarithmic alpha values
 *          log_alpha such that alpha = sgn(log_alpha) 2^(abs(log_alpha)).
 *
 * \param[in] distribution        The distribution to sample.
 * \param[in, out] random_state   The random state to use when sampling.
 *
 * \param[out] min_log_alpha_r    The minimum signed logarithmic alpha_r.
 * \param[out] max_log_alpha_r    The maximum signed logarithmic alpha_r.
 *
 * \return  Returns #TRUE if a region was successfully sampled, #FALSE if the
 *          region sampled is outside the range of the distribution.
 */
bool diagonal_distribution_sample_region(
  const Diagonal_Distribution * const distribution,
  Random_State * const random_state,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r);

/*!
 * \brief   Samples an argument alpha_r from the probability distribution.
 * 
 * The argument alpha_r is sampled with high resolution, and is guaranteed to be
 * admissible.
 *
 * \param[in] distribution        The distribution to sample.
 * \param[in, out] random_state   The random state to use when sampling.
 * 
 * \param[in, out] alpha_r        The argument alpha_r.
 *
 * \return  Returns #TRUE if the argument alpha was successfully sampled, #FALSE
 *          if the argument sampled is outside the range of the distribution.
 */
bool diagonal_distribution_sample_alpha_r(
  const Diagonal_Distribution * const distribution,
  Random_State * const random_state,
  mpz_t alpha_r);

/*!
 * \brief   Samples an integer pair (j, k) from the probability distribution and
 *          returns the components of the pair.
 *
 * \param[in] distribution        The distribution to sample.
 * \param[in, out] random_state   The random state to use when sampling.
 * 
 * \param[in, out] j              The integer j.
 * \param[in, out] k              The integer k.
 */
bool diagonal_distribution_sample_pair_j_k(
  const Diagonal_Distribution * const distribution,
  Random_State * const random_state,
  mpz_t j,
  mpz_t k);

/*!
 * \}
 */

#endif /* DIAGONAL_DISTRIBUTION_H */
