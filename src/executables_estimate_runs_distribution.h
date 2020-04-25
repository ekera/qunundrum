/*!
 * \file    executables_estimate_runs_distribution.h
 * \ingroup volume_executable
 *
 * \brief   The declaration of constants shared by the executables for
 *          estimating volume quotients and the minimum number of runs required
 *          to keep the quotient below a given bound.
 */

/*!
 * \defgroup volume_executable Executables for estimating volume quotients
 * \ingroup  executable
 *
 * \brief    A module for executables for estimating volume quotients.
 */

#ifndef EXECUTABLES_ESTIMATE_RUNS_DISTRIBUTION_H
#define EXECUTABLES_ESTIMATE_RUNS_DISTRIBUTION_H

/*!
 * \brief   The number of problem instances to create per chunk.
 */
#define TAU_CHUNK_SIZE                         1000

/*!
 * \brief   The number of chunks.
 */
#define TAU_CHUNK_COUNT                        1000

/*!
 * \brief   An MPI tag used to send and receive the chunk index.
 */
#define MPI_TAG_CHUNK                         19001

/*!
 * \brief   The maximum number of quantum algorithm outputs that may be
 *          included in a problem instance.
 */
#define MAX_N                                   512

#endif /* EXECUTABLES_ESTIMATE_RUNS_DISTRIBUTION_H */
