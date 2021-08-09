/*!
 * \file    diagonal_distribution_slice.h
 * \ingroup diagonal_distribution_slice
 *
 * \brief   The declaration of data structures representing slices in diagonal
 *          probability distributions, and of functions for manipulating such
 *          slices.
 */

/*!
 * \defgroup  diagonal_distribution_slice Diagonal distribution slices
 * \ingroup   diagonal_distribution
 *
 * \brief   A module for diagonal probability distribution slices.
 */

#ifndef DIAGONAL_DISTRIBUTION_SLICE_H
#define DIAGONAL_DISTRIBUTION_SLICE_H

#include "parameters.h"

#include <stdint.h>

/*!
 * \brief   A data structure representing slices in diagonal probability
 *          distributions.
 *
 * \ingroup diagonal_distribution
 */
typedef struct {
  /*!
   * \brief   The dimension of the slice.
   *
   * A slice of dimension d contains d regions.
   */
  uint32_t dimension;

  /*!
   * \brief   The minimum signed logarithmic alpha_r-coordinate of the slice.
   *
   * The slice covers the interval
   * 
   *    [sgn abs(min_log_alpha_r), sgn (abs(min_log_alpha_r) + 1))
   * 
   * where sgn = sgn(min_log_alpha_r).
   */
  int32_t min_log_alpha_r;

  /*!
   * \brief   The offset from the optimal value of alpha_d.
   *
   * The argument alpha_d(alpha_r) = round( (d/r) alpha_r ) + offset_alpha_d.
   */
  int32_t offset_alpha_d;

  /*!
   * \brief   The total probability mass summed over all regions in the slice.
   */
  long double total_probability;

  /*!
   * \brief   The total error bound summed over all regions in the slice.
   *
   * \remark  Currently not used. Set to zero.
   */
  long double total_error;

  /*!
   * \brief   Flags describing how this slice was created.
   *
   * \remark  Currently not used. Set to zero.
   */
  uint32_t flags;

  /*!
   * \brief   The norm vector.
   *
   * For d the dimension of the slice, the norm vector contains d entries
   * corresponding to the probability mass integrated over each of the d
   * regions in the slice.
   */
  long double * norm_vector;
} Diagonal_Distribution_Slice;

/*!
 * \name Memory allocation
 * \{
 */

/*!
 * \brief   Allocates memory for a slice.
 *
 * \return  The pointer to the memory allocated for the slice.
 */
Diagonal_Distribution_Slice * diagonal_distribution_slice_alloc();

/*!
 * \brief   Deallocates the memory previously allocated for a slice.
 *
 * \param[in, out] slice   A pointer to a pointer to the memory allocated for
 *                         the slice. This function assigns NULL to the pointer
 *                         as it deallocates the memory to which it refers.
 */
void diagonal_distribution_slice_dealloc(
  Diagonal_Distribution_Slice ** const slice);

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
void diagonal_distribution_slice_init(
  Diagonal_Distribution_Slice * const slice,
  const uint32_t dimension);

/*!
 * \brief   Clears a slice.
 *
 * This function deallocates memory allocated for the norm vector in the slice.
 *
 * \param[in, out] slice  The slice to clear.
 */
void diagonal_distribution_slice_clear(
  Diagonal_Distribution_Slice * const slice);

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
void diagonal_distribution_slice_init_copy(
  Diagonal_Distribution_Slice * const dst_slice,
  const Diagonal_Distribution_Slice * const src_slice);

/*!
 * \brief   Copies a slice from a source slice to a destination slice.
 *
 * \param[in, out] dst_slice  The destination slice.
 * \param[in] src_slice       The source slice.
 */
void diagonal_distribution_slice_copy(
  Diagonal_Distribution_Slice * const dst_slice,
  const Diagonal_Distribution_Slice * const src_slice);

/*!
 * \brief   Copies a slice, scaling the source slice to the dimension of the
 *          destination slice whilst copying.
 *
 * \param[in, out] dst_slice  The destination slice.
 * \param[in] src_slice       The source slice.
 */
void diagonal_distribution_slice_copy_scale(
  Diagonal_Distribution_Slice * const dst_slice,
  const Diagonal_Distribution_Slice * const src_slice);

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
void diagonal_distribution_slice_init_import(
 Diagonal_Distribution_Slice * const slice,
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
void diagonal_distribution_slice_import(
  Diagonal_Distribution_Slice * const slice,
  FILE * const file);

/*!
 * \brief   Export a slice to file.
 *
 * \param[in, out] slice  The slice to export to file.
 * \param[in, out] file   The file from which to export the slice.
 */
void diagonal_distribution_slice_export(
  const Diagonal_Distribution_Slice * const slice,
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
void diagonal_distribution_slice_init_recv(
  Diagonal_Distribution_Slice * const slice,
  const int rank);

/*!
 * \brief   Initializes a slice by receiving a broadcast.
 *
 * \param[in, out] slice  The slice to initialize.
 * \param[in] root        The rank of the node broadcasting the slice.
 */
void diagonal_distribution_slice_init_bcast_recv(
  Diagonal_Distribution_Slice * const slice,
  const int root);

/*!
 * \brief   Sends a slice to another node.
 *
 * \param[in, out] slice  The slice to send.
 * \param[in] rank        The rank of the node to which to send the slice.
 */
void diagonal_distribution_slice_send(
  const Diagonal_Distribution_Slice * const slice,
  const int rank);

/*!
 * \brief   Receives a slice from another node.
 *
 * \param[in, out] slice  The slice in which to store the slice received.
 * \param[in] rank        The rank of the node from which to receive the slice.
 */
void diagonal_distribution_slice_recv(
  Diagonal_Distribution_Slice * const slice,
  const int rank);

/*!
 * \brief   Sends a broadcast of a slice.
 *
 * \param[in, out] slice  The slice to send.
 * \param[in] root        The rank of the node broadcasting the slice.
 */
void diagonal_distribution_slice_bcast_send(
  Diagonal_Distribution_Slice * const slice,
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
 * \remark  This function is still a work in progress. It is quite sensitive
 *          precision-wise, in that small changes in alpha_d have large impact
 *          on the probability obtained for a given alpha_r.
 * 
 * \param[in, out] slice      The slice to compute.
 * \param[in] parameters      The parameters for which to compute the slice.
 * \param[in] min_log_alpha_r The signed logarithmic alpha_r-coordinate of the
 *                            slice to compute.
 * \param[in] offset_alpha_d  The offset from the optimal value of alpha_d.
 */
void diagonal_distribution_slice_compute(
  Diagonal_Distribution_Slice * const slice,
  const Parameters * const parameters,
  const int32_t min_log_alpha_r,
  const int32_t offset_alpha_d);

/*!
 * \brief Computes a slice using Simpson's method of numerical integration, 
 *        followed by Richardson-extrapolation to cancel linear error terms.
 *
 * \param[in, out] slice      The slice to compute.
 * \param[in] parameters      The parameters for which to compute the slice.
 * \param[in] min_log_alpha_r The signed logarithmic alpha_r-coordinate of the 
 *                            slice to compute.
 * \param[in] offset_alpha_d  The offset from the optimal value of alpha_d.
 */
void diagonal_distribution_slice_compute_richardson(
  Diagonal_Distribution_Slice * const slice,
  const Parameters * const parameters,
  const int32_t min_log_alpha_r,
  const int32_t offset_alpha_d);

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
 * \param[in] slice             The slice for which to return the coordinates.
 *
 * \param[out] min_log_alpha_r  The minimum signed logarithmic alpha_r. May be 
 *                              set to NULL if this coordinate is not needed.
 * \param[out] max_log_alpha_r  The maximum signed logarithmic alpha_r. May be 
 *                              set to NULL if this coordinate is not needed.
 * \param[out] offset_alpha_d   The offset from the optimal value of alpha_d.
 *                              May be set to NULL if the offset is not needed.
 */
void diagonal_distribution_slice_coordinates(
  const Diagonal_Distribution_Slice * const slice,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r,
  int32_t * const offset_alpha_d);

/*!
 * \brief   Returns the coordinates of j:th region in the slice in the signed
 *          logarithmic argument axis.
 *
 * \param[in] slice             The slice for which to return the coordinates.
 * \param[in] j                 The index of the region corresponding to the
 *                              index in the norm vector stored in the data
 *                              structure for the slice.
 *
 * \param[out] min_log_alpha_r  The minimum signed logarithmic alpha_r. May be 
 *                              set to NULL if this coordinate is not needed.
 * \param[out] max_log_alpha_r  The maximum signed logarithmic alpha_r. May be 
 *                              set to NULL if this coordinate is not needed.
 * \param[out] offset_alpha_d   The offset from the optimal value of alpha_d.
 *                              May be set to NULL if the offset is not needed.
 */
void diagonal_distribution_slice_region_coordinates(
  const Diagonal_Distribution_Slice * const slice,
  const uint32_t j,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r,
  int32_t * const offset_alpha_d);

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
 * \param[out] min_log_alpha_r  The minimum signed logarithmic alpha_r.
 * \param[out] max_log_alpha_r  The maximum signed logarithmic alpha_r.
 * \param[out] offset_alpha_d  The offset from the optimal value of alpha_d.
 */
void diagonal_distribution_slice_sample_region(
  const Diagonal_Distribution_Slice * const slice,
  Random_State * const random_state,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r,
  int32_t * const offset_alpha_d);

/*!
 * \}
 */

#endif /* DIAGONAL_DISTRIBUTION_SLICE_H */
