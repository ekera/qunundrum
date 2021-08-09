/*!
 * \file    linear_distribution_slice.h
 * \ingroup linear_distribution_slice
 *
 * \brief   The declaration of data structures representing slices in linear
 *          probability distributions, and of functions for manipulating such
 *          slices.
 */

/*!
 * \defgroup linear_distribution_slice Linear distribution slices
 * \ingroup  linear_distribution
 *
 * \brief    A module for linear probability distribution slices.
 */

#ifndef LINEAR_DISTRIBUTION_SLICE_H
#define LINEAR_DISTRIBUTION_SLICE_H

#include "parameters.h"

#include <stdint.h>

/*!
 * \brief   An enumeration of targets when computing slices in linear 
 *          probability distributions.
 */
typedef enum {
  /*!
   * \brief Computes the slice for a short discrete logarithm d.
   */
  LINEAR_DISTRIBUTION_SLICE_COMPUTE_TARGET_D = 0,

  /*!
   * \brief Computes the slice for an order r.
   */
  LINEAR_DISTRIBUTION_SLICE_COMPUTE_TARGET_R = 1
} Linear_Distribution_Slice_Compute_Target;

/*!
 * \brief   A data structure representing slices in linear probability
 *          distributions.
 *
 * \ingroup linear_distribution
 */
typedef struct {
  /*!
   * \brief The dimension of the slice.
   *
   * A slice of dimension d contains d regions.
   */
  uint32_t dimension;

  /*!
   * \brief The minimum signed logarithmic alpha-coordinate of the slice.
   *
   *
   * The slice covers the interval
   * 
   *    [sgn abs(min_log_alpha), sgn (abs(min_log_alpha) + 1))
   * 
   * where sgn = sgn(min_log_alpha).
   */
  int32_t min_log_alpha;

  /*!
   * \brief The total probability mass summed over all regions in the slice.
   */
  long double total_probability;

  /*!
   * \brief The total error bound summed over all regions in the slice.
   */
  long double total_error;

  /*!
   * \brief Flags describing how this slice was created.
   */
  uint32_t flags;

  /*!
   * \brief The norm vector.
   *
   * For d the dimension of the slice, the norm vector contains d entries
   * corresponding to the probability mass integrated over each of the d
   * regions in the slice.
   */
  long double * norm_vector;
} Linear_Distribution_Slice;

/*!
 * \name Memory allocation
 * \{
 */

/*!
 * \brief   Allocates memory for a slice.
 *
 * \return  The pointer to the memory allocated for the slice.
 */
Linear_Distribution_Slice * linear_distribution_slice_alloc();

/*!
 * \brief   Deallocates the memory previously allocated for a slice.
 *
 * \param[in, out] slice   A pointer to a pointer to the memory allocated for
 *                         the slice. This function assigns NULL to the pointer
 *                         as it deallocates the memory to which it refers.
 */
void linear_distribution_slice_dealloc(
  Linear_Distribution_Slice ** const slice);

/*!
 * \}
 */

/*!
 * \name Initialization
 * \{
 */

/*!
 * \brief   Initializes a slice.
 *
 * This function allocates memory for the norm vector in the slice.
 *
 * \param[in, out] slice  The slice to initialize.
 * \param[in] dimension   The dimension of the slice.
 */
void linear_distribution_slice_init(
  Linear_Distribution_Slice * const slice,
  const uint32_t dimension);

/*!
 * \brief   Clears a slice.
 *
 * This function deallocates memory allocated for the norm vector in the slice.
 *
 * \param[in, out] slice  The slice to clear.
 */
void linear_distribution_slice_clear(
  Linear_Distribution_Slice * const slice);

/*!
 * \}
 */

/*!
 * \name Copying and scaling slices
 * \{
 */

/*!
 * \brief   Initializes a slice by reading from an already initialized source
 *          slice.
 *
 * \param[in, out] dst_slice  The destination slice to initialize.
 * \param[in] src_slice       The source slice from which to read.
 */
void linear_distribution_slice_init_copy(
  Linear_Distribution_Slice * const dst_slice,
  const Linear_Distribution_Slice * const src_slice);

/*!
 * \brief   Copies a slice from a source slice to a destination slice.
 *
 * \param[in, out] dst_slice  The destination slice.
 * \param[in] src_slice       The source slice.
 */
void linear_distribution_slice_copy(
  Linear_Distribution_Slice * const dst_slice,
  const Linear_Distribution_Slice * const src_slice);

/*!
 * \brief   Copies a slice, scaling the source slice to the dimension of the
 *          destination slice whilst copying.
 *
 * \param[in, out] dst_slice  The destination slice.
 * \param[in] src_slice       The source slice.
 */
void linear_distribution_slice_copy_scale(
  Linear_Distribution_Slice * const dst_slice,
  const Linear_Distribution_Slice * const src_slice);

/*!
 * \}
 */

/*!
 * \name Importing and exporting
 * \{
 */

/*!
 * \brief   Initializes a slice by importing it from a file.
 *
 * \param[in, out] slice  The slice to initialize.
 * \param[in, out] file   The file from which to read the slice.
 */
void linear_distribution_slice_init_import(
 Linear_Distribution_Slice * const slice,
 FILE * const file);

/*!
 * \brief   Imports a slice from file into an already initialized slice.
 *
 * If the dimension of the slice imported does not match the dimension for
 * which the destination slice was initialized, this function will
 * re-initialized the destination slice. Hence, calling this function may
 * result in memory being deallocated and allocated.
 *
 * \param[in, out] slice  The destination slice.
 * \param[in, out] file   The file from which to read the source slice.
 */
void linear_distribution_slice_import(
  Linear_Distribution_Slice * const slice,
  FILE * const file);

/*!
 * \brief   Export a slice to file.
 *
 * \param[in] slice       The slice to export to file.
 * \param[in, out] file   The file from which to export the slice.
 */
void linear_distribution_slice_export(
  const Linear_Distribution_Slice * const slice,
  FILE * const file);

/*!
 * \}
 */

/*!
 * \name Message passing
 * \{
 */

/*!
 * \brief   Initializes a slice by receiving it from another node.
 *
 * \param[in, out] slice  The slice to initialize.
 * \param[in] rank        The rank of the node from which to receive the slice.
 */
void linear_distribution_slice_init_recv(
  Linear_Distribution_Slice * const slice,
  const int rank);

/*!
 * \brief   Initializes a slice by receiving a broadcast.
 *
 * \param[in, out] slice  The slice to initialize.
 * \param[in] root        The rank of the node broadcasting the slice.
 */
void linear_distribution_slice_init_bcast_recv(
  Linear_Distribution_Slice * const slice,
  const int root);

/*!
 * \brief   Sends a slice to another node.
 *
 * \param[in] slice       The slice to send.
 * \param[in] rank        The rank of the node to which to send the slice.
 */
void linear_distribution_slice_send(
  const Linear_Distribution_Slice * const slice,
  const int rank);

/*!
 * \brief   Receives a slice from another node.
 *
 * \param[in, out] slice  The slice in which to store the slice received.
 * \param[in] rank        The rank of the node from which to receive the slice.
 */
void linear_distribution_slice_recv(
  Linear_Distribution_Slice * const slice,
  const int rank);

/*!
 * \brief   Sends a broadcast of a slice.
 *
 * \param[in, out] slice  The slice to send.
 * \param[in] root        The rank of the node broadcasting the slice.
 */
void linear_distribution_slice_bcast_send(
  Linear_Distribution_Slice * const slice,
  const int root);

/*!
 * \}
 */

/*!
 * \name Computing
 * \{
 */

/*!
 * \brief Computes a slice using Simpson's method of numerical integration.
 *
 * \param[in, out] slice      The slice to compute.
 * \param[in] parameters      The parameters for which to compute the slice.
 * \param[in] target          The target. Either a short discrete logarithm d or
 *                            an order r.
 * \param[in] min_log_alpha   The signed logarithmic alpha-coordinate of the 
 *                            slice to compute.
 */
void linear_distribution_slice_compute(
  Linear_Distribution_Slice * const slice,
  const Parameters * const parameters,
  const Linear_Distribution_Slice_Compute_Target target,
  const int32_t min_log_alpha);

/*!
 * \brief Computes a slice using Simpson's method of numerical integration, 
 *        followed by Richardson-extrapolation to cancel linear error terms.
 *
 * \param[in, out] slice      The slice to compute.
 * \param[in] parameters      The parameters for which to compute the slice.
 * \param[in] target          The target. Either a short discrete logarithm d or
 *                            an order r.
 * \param[in] min_log_alpha   The signed logarithmic alpha-coordinate of the 
 *                            slice to compute.
 */
void linear_distribution_slice_compute_richardson(
  Linear_Distribution_Slice * const slice,
  const Parameters * const parameters,
  const Linear_Distribution_Slice_Compute_Target target,
  const int32_t min_log_alpha);

/*!
 * \}
 */

/*!
 * \name Coordinates
 * \{
 */

/*!
 * \brief   Returns the coordinates on the signed logarithmic argument axis.
 *
 * \param[in] slice           The slice for which to return the coordinates.
 *
 * \param[out] min_log_alpha  The minimum signed logarithmic alpha_d or alpha_r.
 *                            May be set to NULL if this coordinate is not 
 *                            needed.
 * \param[out] max_log_alpha  The maximum signed logarithmic alpha_d or alpha_r.
 *                            May be set to NULL if this coordinate is not 
 *                            needed.
 */
void linear_distribution_slice_coordinates(
  const Linear_Distribution_Slice * const slice,
  double * const min_log_alpha,
  double * const max_log_alpha);

/*!
 * \brief   Returns the coordinates of j:th region in the slice in the signed
 *          logarithmic argument axis.
 *
 * \param[in] slice           The slice for which to return the coordinates.
 * \param[in] j               The index of the region corresponding to the
 *                            index in the norm vector stored in the data
 *                            structure for the slice.
 *
 * \param[out] min_log_alpha  The minimum signed logarithmic alpha_d or alpha_r.
 *                            May be set to NULL if this coordinate is not 
 *                            needed.
 * \param[out] max_log_alpha  The maximum signed logarithmic alpha_d or alpha_r.
 *                            May be set to NULL if this coordinate is not 
 *                            needed.
 */
void linear_distribution_slice_region_coordinates(
  const Linear_Distribution_Slice * const slice,
  const uint32_t j,
  double * const min_log_alpha,
  double * const max_log_alpha);

/*!
 * \}
 */

/*!
 * \name Sampling
 * \{
 */

/*!
 * \brief   Samples a region from a slice.
 *
 * \param[in] slice             The slice from which to sample.
 * \param[in, out] random_state The random state from which to read random data.
 *
 * \param[out] min_log_alpha    The minimum signed logarithmic of alpha_d or 
 *                              alpha_r.
 * \param[out] max_log_alpha    The maximum signed logarithmic of alpha_d or 
 *                              alpha_r.
 */
void linear_distribution_slice_sample_region(
  const Linear_Distribution_Slice * const slice,
  Random_State * const random_state,
  double * const min_log_alpha,
  double * const max_log_alpha);

/*!
 * \}
 */

#endif /* LINEAR_DISTRIBUTION_SLICE_H */
