/*!
 * \file    common.h
 *
 * \brief   The definition of various commonly used constants.
 */

#ifndef COMMON_H
#define COMMON_H

/*!
 * \brief   The default precision used in floating point arithmetic.
 */
#define PRECISION                                 192

/*!
 * \name Booleans
 * \{
 */

/*!
 * \brief   The boolean value true.
 */
#ifndef TRUE
#define TRUE                                      1
#endif

/*!
 * \brief  The boolean value false.
 */
#ifndef FALSE
#define FALSE                                     0
#endif

/*!
 * \}
 */

/*!
 * \name MPI
 * \{
 */

/*!
 * \brief   The MPI rank of the standard root node.
 */
#define MPI_RANK_ROOT                             0

/*!
 * \}
 */

/*!
 * \name MPI notifications
 * \{
 */

/*!
 * \brief   The definition of an MPI tag used to send notifications.
 */
#define MPI_TAG_NOTIFY                            17010

/*!
 * \brief   A notification code used with the MPI_TAG_NOTIFY tag to indicate
 *          that the node is ready.
 */
#define MPI_NOTIFY_READY                          17011

/*!
 * \brief   A notification code used with the MPI_TAG_NOTIFY tag to indicate
 *          that the node is has computed a slice.
 */
#define MPI_NOTIFY_SLICE_DONE                     17112

/*!
 * \brief   A notification code used with the MPI_TAG_NOTIFY tag to indicate
 *          that the node has computed a solution.
 */
#define MPI_NOTIFY_SOLUTION_DONE                  17113

/*!
 * \brief   A notification code used with the MPI_TAG_NOTIFY tag to indicate
 *          that the node skips over computing a slice.
 */
#define MPI_NOTIFY_SLICE_SKIP                     17114

/*!
 * \brief   A notification code used with the MPI_TAG_NOTIFY tag to indicate
 *          that a sampling operation was successful.
 */
#define MPI_NOTIFY_SAMPLING_SUCCESSFUL            17115

/*!
 * \brief   A notification code used with the MPI_TAG_NOTIFY tag to indicate
 *          that a sampling operation failed due to the sampled value being
 *          outside of the range covered by the stored probability distribution.
 */
#define MPI_NOTIFY_SAMPLING_FAILED_OUT_OF_BOUNDS  17116

/*!
 * \}
 */

/*!
 * \name MPI tags
 * \{
 */

/*!
 * \brief   An MPI tag used to send and receive min_log_alpha_d, min_log_alpha_r
 *          or the pair (min_log_alpha_d, min_log_alpha_r) for a slice.
 */
#define MPI_TAG_SLICE_MIN_LOG_ALPHA               17021

/*!
 * \brief   An MPI tag used to send and receive slice norm vectors or matrices.
 */
#define MPI_TAG_SLICE_NORM_MATRIX                 17022

/*!
 * \brief   An MPI tag used to send and receive the total error for a slice.
 */
#define MPI_TAG_SLICE_TOTAL_ERROR                 17026

/*!
 * \brief   An MPI tag used to send and receive the flags for a slice.
 */
#define MPI_TAG_SLICE_FLAGS                       17023

/*!
 * \brief   An MPI tag used to send and receive the dimension of a slice.
 */
#define MPI_TAG_SLICE_DIMENSION                   17024

/*!
 * \brief   An MPI tag used to send and receive the number of samples n that
 *          are to be used to setup the problem instance.
 */
#define MPI_TAG_SAMPLE_COUNT                      17025

/*!
 * \brief   An MPI tag used to send and receive information about the solution
 *          of a problem instance.
 */
#define MPI_TAG_SOLUTION                          17027

/*!
 * \}
 */

/*!
 * \name MPI jobs
 * \{
 */

/*!
 * \brief   The definition of an MPI tag used to send out jobs.
 */
#define MPI_TAG_JOB                               17030

/*!
 * \brief   The job code used to indicate the node should stop processing jobs.
 */
#define MPI_JOB_STOP                              17031

/*!
 * \brief   The job code used to indicate the node should process a slice.
 */
#define MPI_JOB_PROCESS_SLICE                     17132

/*!
 * \brief   The job code used to indicate the node should estimate tau.
 */
#define MPI_JOB_ESTIMATE_TAU                      17232

/*!
 * \brief   The job code used to indicate the node should solve a system.
 */
#define MPI_JOB_SOLVE_SYSTEM                      17332

/*!
 * \brief   The job code used to indicate the node should sleep for a fixed
 *          time period and then report back to the server with a job request.
 */
#define MPI_JOB_SLEEP                             17432

/*!
 * \}
 */

/*!
 * \name Slice flags
 * \{
 */

/*!
 * \brief   A flag that indicates that the error bound was not respected.
 */
#define SLICE_FLAGS_ERROR_BOUND_WARNING           0x00000001

/*!
 * \brief   A flag that indicates that the slice was created by mirroring
 *          another slice.
 */
#define SLICE_FLAGS_MIRRORED                      0x00000100

/*!
 * \brief   A flag that indicates that Simpson's method of numerical
 *          integration was used.
 */
#define SLICE_FLAGS_METHOD_SIMPSON                0x00020000

/*!
 * \brief   A flag that indicates that Richardson extrapolation was used.
 *
 * \remark  This flag may be combined with flags indicating which numerical
 *          integration method was used such as #SLICE_FLAGS_METHOD_SIMPSON.
 */
#define SLICE_FLAGS_METHOD_RICHARDSON             0x00080000

/*!
 * \brief   A mask for the integration method in the slice flags.
 */
#define SLICE_FLAGS_MASK_METHOD                   0x000f0000

/*!
 * \brief   A flag that indicates that the slice was created by scaling another
 *          slice to a reduced dimension.
 */
#define SLICE_FLAGS_SCALED                        0x00100000

/*!
 * \}
 */

#endif /* COMMON_H */
