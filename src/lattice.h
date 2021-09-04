/*!
 * \file    lattice.h
 * \ingroup lattice
 *
 * \brief   The definition of constants and enumerations used by the various
 *          solvers for lattice problem.
 */

/*!
 * \defgroup lattice Lattices
 * \ingroup  math
 *
 * \brief    A module for functions for solving various lattice problems.
 */

#ifndef LATTICE_H
#define LATTICE_H

#include "parameters.h"
#include "random.h"

/*!
 * \brief   The default precision to use when performing Gram-Schmidt
 *          orthogonalization and executing Babai's algorithm.
 */
#define LATTICE_DEFAULT_PRECISION                 4096

/*!
 * \brief   The maximum block size to use in the block Korkin-Zolotarev (BKZ)
 *          lattice basis reduction algorithm.
 */
#define LATTICE_MAX_BLOCK_SIZE_BKZ                10

/*!
 * \brief   An enumeration of recovery status codes.
 */
typedef enum {
  /*!
   * \brief   Indicates that the quantity was not recovered.
   */
  LATTICE_STATUS_NOT_RECOVERED = 0,

  /*!
   * \brief   Indicates that the quantity was immediately recovered by either
   *          extracting the shortest vector in the reduced lattice basis or
   *          by mapping a target vector to the closest vector in the lattice
   *          using the reduced basis and Babai's nearest plane algorithm.
   */
  LATTICE_STATUS_RECOVERED_IMMEDIATE = 1,

  /*!
   * \brief   Indicates that the quantity was recovered by either extracting the
   *          shortest vector in the reduced lattice basis, or by mapping a
   *          target vector to the closest vector in the lattice using the
   *          reduced basis and Babai's nearest plane algorithm, and performing
   *          a small linear search in the last lattice component.
   */
  LATTICE_STATUS_RECOVERED_SEARCH = 2,

  /*!
   * \brief   Indicates that the quantity was recovered by enumerating the 
   *          lattice.
   */
  LATTICE_STATUS_RECOVERED_ENUMERATION = 3,

  /*!
   * \brief   Indicates that a timeout has occurred.
   */
  LATTICE_STATUS_TIMEOUT = 4,

  /*!
   * \brief   Indicates that the quantity was recovered, except for a cm-smooth
   *          factor that may be efficiently found in an additional classical 
   *          post-processing step. */ 
  LATTICE_STATUS_RECOVERED_IMMEDIATE_SMOOTH = 5,

  /*!
   * \brief   Indicates that the discrete logarithm, modularly reduced by r/z
   *          for some smooth z > 1, was recovered by enumerating the lattice.
   */
  LATTICE_STATUS_RECOVERED_ENUMERATION_SMOOTH = 6
} Lattice_Status_Recovery;

/*!
 * \brief   An enumeration of lattice basis reduction algorithms.
 */
typedef enum {
  /*!
   * \brief   Specifies that the Lenstra-Lenstra-Lovász (LLL) reduction
   *          algorithm should be used to reduce the lattice basis.
   */
  REDUCTION_ALGORITHM_LLL = 1,

  /*!
   * \brief  Specifies that the block Korkin-Zolotarev (BKZ) reduction
   *         algorithm should be used to reduce the lattice basis.
   *
   * A reduced BKZ basis is more strongly reduced than a Lenstra-Lenstra-Lovász-
   * reduced (LLL) basis. However a BKZ-reduced basis takes longer time to
   * compute than an LLL-reduced basis.
   */
  REDUCTION_ALGORITHM_BKZ = 2,

  /*!
   * \brief   Specifies that the Lenstra-Lenstra-Lovász (LLL) reduction
   *          algorithm should be used to reduce the lattice basis. If recovery
   *          of either d or r fails for the LLL basis, then a second attempt
   *          should be made with a block Korkin-Zolotarev (BKZ) reduced basis.
   *
   * A reduced BKZ basis is more strongly reduced than a Lenstra-Lenstra-Lovász-
   * reduced (LLL) basis. However a BKZ-reduced basis takes longer time to
   * compute than an LLL-reduced basis. When computing the BKZ-reduced basis,
   * the LLL-reduced basis is used as the starting point to speed up the
   * computation.
   *
   * The BKZ block size is set to ten if the lattice dimension is greater than
   * ten, otherwise it is set to the dimension of the lattice.
   */
  REDUCTION_ALGORITHM_LLL_BKZ = 3,

  /*!
   * \brief  Specifies that the Hermite Korkin-Zolotarev (BKZ) reduction
   *         algorithm should be used to reduce the lattice basis.
   */
  REDUCTION_ALGORITHM_HKZ = 4,

  /*!
   * \brief   Indicates that no reduction algorithm has been specified.
   */
  REDUCTION_ALGORITHM_DEFAULT = 0
} Lattice_Reduction_Algorithm;

#endif /* LATTICE_H */
