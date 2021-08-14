/*!
 * \file    distribution_slice.h
 * \ingroup two_dimensional_distribution_slice
 *
 * \brief   The declaration of data structures representing slices in
 *          two-dimensional probability distributions, and of functions for
 *          manipulating such slices.
 */

/*!
 * \defgroup  two_dimensional_distribution_slice \
 *            Two-dimensional distribution slices
 * \ingroup   two_dimensional_distribution
 *
 * \brief   A module for two-dimensional probability distribution slices.
 */

#ifndef DISTRIBUTION_SLICE_H
#define DISTRIBUTION_SLICE_H

#include "parameters.h"
#include "random.h"

#include <stdint.h>
#include <stdio.h>

/*!
 * \brief   An enumeration of methods used to compute slices in two-dimensional
 *          probability distributions.
 */
typedef enum {
  /*!
   * \brief Computes the slice, with Richardson-extrapolation, with the
   *        error-bounded probability approximation by Ekerå, with heuristically
   *        globally selected sigma.
   * 
   * The error-bounded approximation in the paper by Ekerå [1] on computing 
   * general discrete logarithms and orders with tradeoffs is used to compute 
   * the distribution slices. The approximation is parameterized in sigma.
   *
   * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms
   * and orders with tradeoffs. In: IACR ePrint Archive, 2018/797.
   * 
   * Sigma is selected using the heuristic in [1] for fixed tau = 11.
   */
  DISTRIBUTION_SLICE_COMPUTE_METHOD_HEURISTIC_SIGMA = 0,

  /*!
   * \brief Computes the slice, with Richardson-extrapolation, with the
   *        error-bounded probability approximation by Ekerå, with adaptively
   *        selected sigma to locally minimize the error bound.
   * 
   * The error-bounded approximation in the paper by Ekerå [1] on computing 
   * general discrete logarithms and orders with tradeoffs is used to compute 
   * the distribution slices. The approximation is parameterized in sigma.
   *
   * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms
   * and orders with tradeoffs. In: IACR ePrint Archive, 2018/797.
   * 
   * Sigma is selected optimally at the outset, and adapted locally in each step
   * by seeing if increasing or decreasing sigma yields a lower error bound.
   */
  DISTRIBUTION_SLICE_COMPUTE_METHOD_OPTIMAL_LOCAL_SIGMA = 1,

  /*!
   * \brief Computes the slice with the quick and dirty heuristic approximation 
   *        by Ekerå that lacks an error bound.
   * 
   * The heuristic approximation in the paper by Ekerå [1] on computing general
   * discrete logarithms and orders with tradeoffs is used to compute the 
   * distribution slices. Note that there is currently no error bound available 
   * for this approximation.
   * 
   * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms
   * and orders with tradeoffs. In: IACR ePrint Archive, 2018/797.
   */
  DISTRIBUTION_SLICE_COMPUTE_METHOD_QUICK = 2
} Distribution_Slice_Compute_Method;

/*!
 * \brief   A data structure representing slices in two-dimensional probability
 *          distributions.
 *
 * \ingroup two_dimensional_distribution
 */
typedef struct {
  /*!
   * \brief The dimension of the slice.
   *
   * A slice of dimension d contains d x d regions.
   */
  uint32_t dimension;

  /*!
   * \brief The minimum signed logarithmic alpha_d-coordinate of the slice.
   *
   * The slice covers the interval
   * 
   *    [sgn abs(min_log_alpha_d), sgn (abs(min_log_alpha_d) + 1))
   * 
   * where sgn = sgn(min_log_alpha_d).
   */
  int32_t min_log_alpha_d;

  /*!
   * \brief The minimum signed logarithmic alpha_r-coordinate of the slice.
   *
   * The slice covers the interval
   * 
   *    [sgn abs(min_log_alpha_r), sgn (abs(min_log_alpha_r) + 1))
   * 
   * where sgn = sgn(min_log_alpha_r).
   */
  int32_t min_log_alpha_r;

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
   * \brief The norm matrix.
   *
   * For d the dimension of the slice, the norm matrix contains d^2 entries
   * corresponding to the probability mass integrated over each of the d^2
   * regions in the slice.
   */
  long double * norm_matrix;
} Distribution_Slice;


/*!
 * \name Memory allocation
 * \{
 */

/*!
 * \brief   Allocates memory for a slice.
 *
 * \return  The pointer to the memory allocated for the slice.
 */
Distribution_Slice * distribution_slice_alloc();

/*!
 * \brief   Deallocates the memory previously allocated for a slice.
 *
 * \param[in, out] slice   A pointer to a pointer to the memory allocated for
 *                         the slice. This function assigns NULL to the pointer
 *                         as it deallocates the memory to which it refers.
 */
void distribution_slice_dealloc(
  Distribution_Slice ** const slice);

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
 * This function allocates memory for the norm matrix in the slice.
 *
 * \param[in, out] slice  The slice to initialize.
 * \param[in] dimension   The dimension of the slice.
 */
void distribution_slice_init(
  Distribution_Slice * const slice,
  const uint32_t dimension);

/*!
 * \brief   Clears a slice.
 *
 * This function deallocates memory allocated for the norm matrix in the slice.
 *
 * \param[in, out] slice  The slice to clear.
 */
void distribution_slice_clear(
  Distribution_Slice * const slice);

/*!
 * \}
 */

/*!
 * \name Copying and scaling
 * \{
 */

/*!
 * \brief    Initializes a slice by reading from an already initialized source
 *           slice.
 *
 * \param[in, out] dst_slice The destination slice to initialize.
 * \param[in] src_slice      The source slice from which to read.
 */
void distribution_slice_init_copy(
 Distribution_Slice * const dst_slice,
 const Distribution_Slice * const src_slice);

/*!
 * \brief   Copies a slice from a source slice to a destination slice.
 *
 * \param[in, out] dst_slice  The destination slice.
 * \param[in] src_slice       The source slice.
 */
void distribution_slice_copy(
  Distribution_Slice * const dst_slice,
  const Distribution_Slice * const src_slice);

/*!
 * \brief   Copies a slice, scaling the source slice to the dimension of the
 *          destination slice whilst copying.
 *
 * \param[in, out] dst_slice  The destination slice.
 * \param[in] src_slice       The source slice.
 */
void distribution_slice_copy_scale(
  Distribution_Slice * const dst_slice,
  const Distribution_Slice * const src_slice);

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
void distribution_slice_init_import(
  Distribution_Slice * const slice,
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
void distribution_slice_import(
  Distribution_Slice * const slice,
  FILE * const file);

/*!
 * \brief   Export a slice to file.
 *
 * \param[in, out] slice  The slice to export to file.
 * \param[in, out] file   The file from which to export the slice.
 */
void distribution_slice_export(
  const Distribution_Slice * const slice,
  FILE * const file);

/*!
 * \}
 */


/*!
 * \name Message passing
 * \{
 */

/*!
 * \brief   Initializes a slice by receiving a message.
 *
 * \param[in, out] slice  The slice to initialize.
 * \param[in] rank        The rank of the node from which to receive the slice.
 */
void distribution_slice_init_recv(
  Distribution_Slice * const slice,
  const int rank);

/*!
 * \brief   Initializes a slice by receiving a broadcast.
 *
 * \param[in, out] slice  The slice to initialize.
 * \param[in] root        The rank of the node broadcasting the slice.
 */
void distribution_slice_init_bcast_recv(
  Distribution_Slice * const slice,
  const int root);

/*!
 * \brief   Sends a slice to another node.
 *
 * \param[in] slice       The slice to send.
 * \param[in] rank        The rank of the node to which to send the slice.
 */
void distribution_slice_send(
  const Distribution_Slice * const slice,
  const int rank);

/*!
 * \brief   Receives a slice from another node.
 *
 * \param[in, out] slice  The slice in which to store the slice received.
 * \param[in] rank        The rank of the node from which to receive the slice.
 */
void distribution_slice_recv(
  Distribution_Slice * const slice,
  const int rank);

/*!
 * \brief   Sends a broadcast of a slice.
 *
 * \param[in, out] slice  The slice to send.
 * \param[in] root        The rank of the node broadcasting the slice.
 */
void distribution_slice_bcast_send(
  Distribution_Slice * const slice,
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
 * \param[in, out] slice        The slice to compute.
 * \param[in] parameters        The parameters for which to compute the slice.
 * \param[in] method            The method to use to compute the size, with 
 *                              respect to which approximation method to use 
 *                              and how to select parameters. See the 
 *                              Distribution_Slice_Compute_Method enumeration 
 *                              for more information on the available options.
 * \param[in] min_log_alpha_d   The signed logarithmic alpha_d-coordinate of the
 *                              slice to compute.
 * \param[in] min_log_alpha_r   The signed logarithmic alpha_r-coordinate of the
 *                              slice to compute.
 */
void distribution_slice_compute(
  Distribution_Slice * const slice,
  const Parameters * const parameters,
  const Distribution_Slice_Compute_Method method,
  const int32_t min_log_alpha_d,
  const int32_t min_log_alpha_r);

/*!
 * \brief Computes a slice using Simpson's method of numerical integration, 
 *        followed by Richardson-extrapolation to cancel linear error terms.
 * 
 * \param[in, out] slice        The slice to compute.
 * \param[in] parameters        The parameters for which to compute the slice.
 * \param[in] method            The method to use to compute the size, with 
 *                              respect to which approximation method to use 
 *                              and how to select parameters. See the 
 *                              Distribution_Slice_Compute_Method enumeration 
 *                              for more information on the available options.
 * \param[in] min_log_alpha_d   The signed logarithmic alpha_d-coordinate of the
 *                              slice to compute.
 * \param[in] min_log_alpha_r   The signed logarithmic alpha_r-coordinate of the
 *                              slice to compute.
 */
void distribution_slice_compute_richardson(
  Distribution_Slice * const slice,
  const Parameters * const parameters,
  const Distribution_Slice_Compute_Method method,
  const int32_t min_log_alpha_d,
  const int32_t min_log_alpha_r);

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
 * \param[out] min_log_alpha_d  The minimum signed logarithmic alpha_d.
 * \param[out] max_log_alpha_d  The maximum signed logarithmic alpha_d.
 * \param[out] min_log_alpha_r  The minimum signed logarithmic alpha_r.
 * \param[out] max_log_alpha_r  The maximum signed logarithmic alpha_r.
 */
void distribution_slice_sample_region(
  const Distribution_Slice * const slice,
  Random_State * const random_state,
  double * const min_log_alpha_d,
  double * const max_log_alpha_d,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r);

/*!
 * \}
 */

/*!
 * \name Coordinates
 * \{
 */

/*!
 * \brief   Returns the coordinates in the signed logarithmic argument plane.
 *
 * \param[in] slice             The slice for which to return the coordinates.
 *
 * \param[out] min_log_alpha_d  The minimum signed logarithmic alpha_d. May be 
 *                              set to NULL if this coordinate is not needed.
 * \param[out] max_log_alpha_d  The maximum signed logarithmic alpha_d. May be 
 *                              set to NULL if this coordinate is not needed.
 * \param[out] min_log_alpha_r  The minimum signed logarithmic alpha_r. May be 
 *                              set to NULL if this coordinate is not needed.
 * \param[out] max_log_alpha_r  The maximum signed logarithmic alpha_r. May be 
 *                              set to NULL if this coordinate is not needed.
 */
void distribution_slice_coordinates(
  const Distribution_Slice * const slice,
  double * const min_log_alpha_d,
  double * const max_log_alpha_d,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r);

/*!
 * \brief   Returns the coordinates of j:th region in the slice in the signed
 *          logarithmic argument plane.
 *
 * \param[in] slice             The slice for which to return the coordinates.
 * \param[in] j                 The index of the region corresponding to the
 *                              index in the norm matrix stored in the data
 *                              structure for the slice.
 *
 * \param[out] min_log_alpha_d  The minimum signed logarithm of alpha_d. May be 
 *                              set to NULL if this coordinate is not needed.
 * \param[out] max_log_alpha_d  The maximum signed logarithm of alpha_d. May be 
 *                              set to NULL if this coordinate is not needed.
 * \param[out] min_log_alpha_r  The minimum signed logarithm of alpha_r. May be 
 *                              set to NULL if this coordinate is not needed.
 * \param[out] max_log_alpha_r  The maximum signed logarithm of alpha_r. May be 
 *                              set to NULL if this coordinate is not needed.
 */
void distribution_slice_region_coordinates(
  const Distribution_Slice * const slice,
  const uint32_t j,
  double * const min_log_alpha_d,
  double * const max_log_alpha_d,
  double * const min_log_alpha_r,
  double * const max_log_alpha_r);

/*!
 * \}
 */

#endif /* DISTRIBUTION_SLICE_H */
