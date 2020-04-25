/*!
 * \file    distribution.h
 * \ingroup two_dimensional_distribution
 *
 * \brief   The declaration of data structures representing two-dimensional
 *          probability distributions, and of functions for manipulating such
 *          distributions.
 *
 * Two-dimensional probability distributions in the signed 
 * (alpha_d, alpha_r)-plane are used to represent the probability distribution 
 * induced by Ekerå's algorithm [1] for computing general discrete logarithms 
 * and orders with tradeoffs.
 *
 * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms 
 * and orders with tradeoffs. In: IACR ePrint Archive, 2018/797.
 */

/*!
 * \defgroup  two_dimensional_distribution \
 *            Two-dimensional probability distributions
 * \ingroup   distribution
 *
 * \brief   A module for two-dimensional probability distributions.
 *
 * Two-dimensional probability distributions in the signed 
 * (alpha_d, alpha_r)-plane are used to represent the probability distribution 
 * induced by Ekerå's algorithm [1] for computing general discrete logarithms 
 * and orders with tradeoffs.
 *
 * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms 
 * and orders with tradeoffs. In: IACR ePrint Archive, 2018/797.
 */

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include "distribution_slice.h"
#include "lattice_sample.h"
#include "random.h"
#include "parameters.h"
#include "common.h"

#include <mpfr.h>
#include <stdint.h>

/*!
 * \brief   A data structure representing a two-dimensional probability
 *          distribution.
 *
 * \ingroup two_dimensional_distribution
 */
typedef struct {
  /*! \brief  A vector of pointers to slices. */
  Distribution_Slice ** slices;

  /*! \brief  The total number of pointers to slices allocated. */
  uint32_t capacity;

  /*! \brief  The number of slices inserted into the distribution. */
  uint32_t count;

  /*! \brief  The distribution parameters. */
  Parameters parameters;

  /*! \brief  The total probability summed over all slices. */
  long double total_probability;

  /*! \brief  The total error summed over all slices. */
  long double total_error;

  /*!
   * \brief   The precision using which the distribution was computed.
   *
   * \remark  Incorrect precisions may be reported for distributions imported
   *          and exported using binaries compiled to use different precision
   *          levels, especially if the distribution is modified by computing
   *          or recomputing slices inbetween import and export operations. */
  uint32_t precision;

  /*! \brief  The alpha lattice for this distribution. */
  Lattice_Alpha lattice_alpha;
} Distribution;

/*!
 * \name Memory allocation
 * \{
 */

/*!
 * \brief   Allocates memory for a distribution.
 *
 * \return  The pointer to the memory allocated for the distribution.
 */
Distribution * distribution_alloc();

/*!
 * \brief   Deallocates the memory previously allocated for a distribution.
 *
 * \param[in, out] distribution   A pointer to a pointer to the memory allocated
 *                                for the distribution. This function assigns
 *                                NULL to the pointer as it deallocates the
 *                                memory to which it refers.
 */
void distribution_dealloc(
  Distribution ** const distribution);

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
 * \param[in] capacity            The distribution capacity in slices.
 *
 * \remark  The parameter data structure is copied into the distribution data
 *          structure, rather than being referenced.
 */
void distribution_init(
  Distribution * const distribution,
  const Parameters * const parameters,
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
 * \remark  Used by the plotting functions to scale down the resolution in
 *          the probability distributions.
 */
void distribution_init_copy_scale(
  Distribution * const dst_distribution,
  const Distribution * const src_distribution,
  const uint32_t dimension);

/*!
 * \brief   Clears a distribution.
 *
 * This function deallocates memory for slices in the distribution.
 *
 * \param[in, out] distribution   The distribution to clear.
 */
void distribution_clear(
  Distribution * const distribution);

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
void distribution_init_import(
  Distribution * const distribution,
  FILE * const file);

/*!
 * \brief   Exports a probability distribution to file.
 *
 * \param[in, out] distribution   The distribution to export to file.
 * \param[in, out] file           The file to which to export the distribution.
 */
void distribution_export(
  const Distribution * const distribution,
  FILE * const file);

/*!
 * \brief   Exports a probability distribution to file, clearing the
 *          distribution in an interleaved manner.
 *
 * Calling this function is equivalent to first calling distribution_export()
 * and then calling distribution_clear() and distribution_dealloc() but this
 * function minimizes the usage of primary memory over time.
 *
 * \param[in, out] distribution   The distribution to export to file.
 * \param[in, out] file           The file to which to export the distribution.
 */
void distribution_export_clear_dealloc(
  Distribution * const distribution,
  FILE * const file);

/*!
 * \brief   Exports information about a probability distribution to file.
 *
 * \param[in, out] distribution   The distribution to export to file.
 * \param[in, out] file           The file to which to export information.
 * \param[in] export_slices       Export verbose information on slices.
 */
void distribution_export_info(
  const Distribution * const distribution,
  FILE * const file,
  const bool export_slices = TRUE);

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
void distribution_init_bcast_recv(
  Distribution * const distribution,
  const int root);

/*!
 * \brief   Broadcasts a distribution to all other nodes.
 *
 * \param[in] distribution        The distribution to broadcast.
 * \param[in] root                The rank of the root node.
 */
void distribution_bcast_send(
  const Distribution * const distribution,
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
void distribution_insert_slice(
  Distribution * const distribution,
  Distribution_Slice * const slice);

/*!
 * \brief   Removes a distribution slice by index.
 *
 * \param[in, out] distribution The distribution from which to remove the slice.
 * \param[in] index             The index of the slice to remove.
 */
void distribution_remove_slice(
  Distribution * const distribution,
  const uint32_t index);

/*!
 * \brief   Replace a distribution slice by (alpha_d, alpha_r) coordinate.
 *
 * \param[in, out] distribution   The distribution in which to replace the
 *                                slice.
 * \param[in, out] slice          The slice to insert, replacing the existing
 *                                slice with the same (alpha_d, alpha_r)
 *                                coordinate pair as this slice.
 *
 * \remark  If this functions replaces a slices in the distribution, it will
 *          deallocate memory for the slice removed by replacement operation.
 *
 * \remark  The slice data structure is referenced from the distribution data
 *          structure, rather than being copied into the structure. You may
 *          not modify or deallocate the slice after calling this function.
 *
 * \return  Return #TRUE if a slice with coordinates (alpha_d, alpha_r) was
 *          found and replaced, #FALSE otherwise.
 */
bool distribution_replace_slice(
  Distribution * const distribution,
  Distribution_Slice * const slice);

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
void distribution_sort_slices(
  Distribution * const distribution);

/*!
 * \brief   Filters the slices in the distribution by removing all slices such
 *          that p = 0 or e / p > 0.01, for p the total probability and e the 
 *          total error in the slice.
 * 
 * \param[in, out] distribution   The distribution in which to filter the 
 *                                slices.
 */
void distribution_filter_slices(
  Distribution * const distribution);

/*!
 * \brief   Tests if filtering the distribution would leave it intact.
 * 
 * \param[in] distribution        The distribution for which to test if 
 *                                filtering would leave it intact.
 * 
 * \return  Returns #TRUE if filtering the distribution would not remove any 
 *          slices, #FALSE otherwise.
 */
bool distribution_is_filtered(
  const Distribution * const distribution);

/*!
 * \}
 */

/*!
 * \name Sampling
 * \{
 */

/*!
 * \brief   Samples a slice from the probability distribution.
 *
 * \param[in] distribution        The distribution to sample.
 * \param[in, out] random_state   The random state to use to sample.
 *
 * \return    A pointer to the slice sampled, or NULL if the slice sampled
 *            is outside the range of the distribution.
 */
const Distribution_Slice * distribution_sample_slice(
  const Distribution * const distribution,
  Random_State * const random_state);

/*!
 * \brief   Samples a rectangular region in the (alpha_d, alpha_r)-argument
 *          plane from the probability distribution and returns its bounds.
 *
 * \remark  The region bounds are represented as signed logarithmic alpha values
 *          log_alpha such that alpha = sgn(log_alpha) 2^(abs(log_alpha)).
 *
 * \param[in] distribution          The distribution to sample.
 * \param[in, out] random_state     The random state to use to sample pivots.
 *
 * \param[out] min_log_alpha_d      The minimum signed logarithmic alpha_d.
 * \param[out] max_log_alpha_d      The maximum signed logarithmic alpha_d.
 * \param[out] min_log_alpha_r      The minimum signed logarithmic alpha_r.
 * \param[out] max_log_alpha_r      The maximum signed logarithmic alpha_r.
 *
 * \return  Returns #TRUE if a region was successfully sampled, #FALSE if the
 *          region sampled is outside the range of the distribution.
 */
bool distribution_sample_region(
  const Distribution * const distribution,
  Random_State * const random_state,
  double * const min_log_alpha_d,
  double * const max_log_alpha_d,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r);

/*!
 * \brief   Samples an approximate argument pair (alpha_d, alpha_r) from the
 *          probability distribution and returns the components of the pair.
 *
 * The pair is coarsely sampled, and is not guaranteed to be admissible.
 *
 * \param[in] distribution          The distribution to sample.
 * \param[in, out] random_state     The random state to use to sample pivots.
 * \param[in, out] alpha_d          The approximate argument alpha_d.
 * \param[in, out] alpha_r          The approximate argument alpha_r.
 *
 * \return  Returns #TRUE if an argument pair was successfully sampled, #FALSE
 *          if the pair sampled is outside the range of the distribution.
 */
bool distribution_sample_approximate_alpha_d_r(
  const Distribution * const distribution,
  Random_State * const random_state,
  mpfr_t alpha_d,
  mpfr_t alpha_r);

/*!
 * \brief   Samples an argument pair (alpha_d, alpha_r) from the probability
 *          distribution and returns the components of the pair.
 *
 * The pair is sampled with high resolution, and is guaranteed to be admissible.
 *
 * To map outputs from this function to pairs (j, k), see the function 
 * sample_j_k_from_alpha_d_r(). See also the conveniency function
 * distribution_sample_pair_j_k() in this module.
 *
 * \param[in] distribution        The distribution to sample.
 * \param[in, out] random_state   The random state the use to sample pivots.
 * \param[in, out] alpha_d        The argument alpha_d.
 * \param[in, out] alpha_r        The argument alpha_r.
 *
 * \return  Returns #TRUE if an argument pair was successfully sampled, #FALSE
 *          if the pair sampled is outside the range of the distribution.
 */
bool distribution_sample_alpha_d_r(
  const Distribution * const distribution,
  Random_State * const random_state,
  mpz_t alpha_d,
  mpz_t alpha_r);

/*!
 * \brief   Samples an integer pair (j, k) from the probability distribution and
 *          returns the components of the pair.
 *
 * This is a conveniency function. Calling this function is equivalent to first
 * calling distribution_sample_alpha_d_r() and then calling the function
 * sample_j_k_from_alpha_d_r().
 *
 * \param[in] distribution        The distribution to sample.
 * \param[in, out] random_state   The random state the use to sample pivots.
 * \param[in, out] j              The integer j.
 * \param[in, out] k              The integer k.
 */
bool distribution_sample_pair_j_k(
  const Distribution * const distribution,
  Random_State * const random_state,
  mpz_t j,
  mpz_t k);

/*!
 * \}
 */

#endif /* DISTRIBUTION_H */
