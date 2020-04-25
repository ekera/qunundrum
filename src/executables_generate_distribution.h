/*!
 * \file    executables_generate_distribution.h
 * \ingroup generate_executable
 *
 * \brief   The declaration of enumerations shared by the executables for
 *          generating probability distributions.
 */

/*!
 * \defgroup generate_executable Executables for generating distributions
 * \ingroup  executable
 *
 * \brief    A module for executables for generating probability distributions.
 */

#ifndef EXECUTABLES_GENERATE_DISTRIBUTION_H
#define EXECUTABLES_GENERATE_DISTRIBUTION_H

/*!
 * \brief    An enumeration of linear distribution types.
 */
typedef enum {
  /*!
   * \brief   Indicates that this distribution is for a discrete logartihm.
   */
  DISTRIBUTION_TYPE_LOGARITHM = 1,

  /*!
   * \brief   Indicates that this distribution is for an order.
   */
  DISTRIBUTION_TYPE_ORDER = 2,

  /*!
   * \brief   Indicates that the distribution type is undefined.
   */
  DISTRIBUTION_TYPE_UNDEFINED = 0
} Distribution_Type;

/*!
 * \brief   An enumeration of selection methods for a short discrete logarithm d
 *          or an order r.
 */
typedef enum {
  /*!
   * \brief   Indicates that the value is explicitly specified.
   */
  SELECTION_METHOD_EXPLICIT = 1,

  /*!
   * \brief   Indicates that the value is deterministically selected from
   *          Catalan's constant.
   */
  SELECTION_METHOD_DETERMINISTIC = 2,

  /*!
   * \brief   Indicates that the value is selected uniformly at random on the
   *           interval [2^(m-1), 2^m).
   */
  SELECTION_METHOD_RANDOM = 3,

  /*!
   * \brief   Indicates that the value is set to 2^(m-1) + 1.
   *
   * This is the smallest odd value on the interval [2^(m-1), 2^m).
   */
  SELECTION_METHOD_MINIMAL = 4,

  /*!
   * \brief   Indicates that the value is set to 2^m - 1.
   *
   * This is the greatest odd value on the interval [2^(m-1), 2^m).
   */
  SELECTION_METHOD_MAXIMAL = 5,

  /*!
   * \brief   Indicates that the selection method is undefined.
   */
  SELECTION_METHOD_UNDEFINED = 0
} Selection_Method;

/*!
 * \brief   An enumeration of tradeoff methods used to select l.
 */
typedef enum {
  /*!
   * \brief   Indicates that l = ceil(m / s) for s a tradeoff factor.
   */
  TRADEOFF_METHOD_FACTOR = 1,

  /*!
   * \brief   Indicates that l is explicitly specified.
   */
  TRADEOFF_METHOD_EXPLICIT = 2,

  /*!
   * \brief   Indicates that the tradeoff method is undefined.
   */
  TRADEOFF_METHOD_UNDEFINED = 0
} Tradeoff_Method;

/*!
 * \brief   An enumeration of sigma selection methods.
 */
typedef enum {
  /*!
   * \brief   Indicates that sigma is heuristically selected to minimize the
   *          error bound.
   *
   * For further details, please see the paper.
   */
  SIGMA_METHOD_HEURISTIC = 1,

  /*!
   * \brief   Indicates that sigma is optimized to minimize the error bound in
   *          each slice using an heuristic exhaustive search.
   *
   * This selection method is considerably slower than #SIGMA_METHOD_HEURISTIC.
   *
   * \remark  This selection method is not guaranteed to find the optimum. The
   *          function seeks to minimumize the error bound; not the error.
   */
  SIGMA_METHOD_OPTIMAL = 2,

  /*!
   * \brief   Indicates that the sigma selection method is undefined.
   */
  SIGMA_METHOD_UNDEFINED = 0
} Sigma_Selection_Method;

typedef enum {
  /*!
   * \brief   Indicates that the error-bounded probability approximation of 
   *          Ekerå is used to compute the distribution.
   *
   * The error-bounded approximation in the paper by Ekerå [1] on computing 
   * general discrete logarithms and orders with tradeoffs is used to compute 
   * the distribution.
   *
   * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms
   * and orders with tradeoffs. In: IACR ePrint Archive, 2018/797.
   */
  PROBABILITY_ESTIMATE_BOUNDED_ERROR = 1,

  /*!
   * \brief   Indicates that the quick and dirty heuristic approximation of 
   *          Ekerå is used to compute the distribution.
   *
   * The quick and dirty approximation in the introduction to the paper [1] by 
   * Ekerå on computing general discrete logarithms and orders with tradeoffs is
   * used to compute the distribution.
   * 
   * [1] Ekerå, M.: Quantum algorithms for computing general discrete logarithms
   * and orders with tradeoffs. In: IACR ePrint Archive, 2018/797.
   */
  PROBABILITY_ESTIMATE_QUICK = 2,

  /*!
   * \brief   Indicates that the probability estimation method is undefined.
   */
  PROBABILITY_ESTIMATE_UNDEFINED = 0
} Probability_Estimate;

#endif /* EXECUTABLES_GENERATE_DISTRIBUTION_H */
