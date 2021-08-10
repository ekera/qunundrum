/*!
 * \file    linear_distribution.h
 * \ingroup linear_distribution
 *
 * \brief   The declaration of data structures representing linear probability
 *          distributions, and of functions for manipulating such distributions.
 *
 * Linear probability distributions are used to represent the probability
 * distributions induced by Shor's order-finding algorithm [1], by Seifert's
 * algorithm for computing order with tradeoffs [2], by Ekerå's algorithm for
 * computing short discrete logarithms [3], and by Ekerå-Håstad's algorithm for
 * computing short discrete logarithms with tradeoffs [4, 5]. The algorithm
 * of Ekerå and Håstad also factors RSA integers by classically reducing the RSA
 * integer factoring problem to a short discrete logarithm problem in a group of
 * unknown order and solving this problem quantumly. 
 * 
 * [1] Shor, P.W.: Polynomial-time algorithms for prime factorization and 
 * discrete logarithms on a quantum computer. In: SIAM Journal on Scientific 
 * Computing (SISC), volume 26(5), pp. 1484 (1997). 
 * 
 * [2] Seifert, J.-P.: Using fewer Qubits in Shor's factorization algorithm via 
 * simultaneous Diophantine approximation. In: CT-RSA 2001, Spring LNCS 2020, 
 * pp. 319-327 (2001). DOI: https://doi.org/10.1007/3-540-45353-9_24.

 * [3] Ekerå, M.: Modifying Shor's algorithm to compute short discrete 
 * logarithms. In: IACR ePrint Archive, 2016/1128.
 * 
 * [4] Ekerå, M. and Håstad, J.: Quantum algorithms for computing short discrete
 * logarithms and factor RSA integers. In: PQCrypto 2017, Springer LNCS 10346,
 * pp. 347-363 (2017). DOI: https://doi.org/10.1007/978-3-319-59879-6_20.
 * 
 * [5] Ekerå, M.: On post-processing in the quantum algorithm for computing 
 * short discrete logarithms. In: IACR ePrint Archive, 2017/1122.
 * 
 * Linear probability distributions may also be used to represent the marginal
 * distributions of two-dimensional probability distributions, collapsed 
 * diagonal distributions, and traces in diagonal distributions for fixed 
 * offsets from the optimal alpha_d for given alpha_r.
 * 
 * Note that Shor's algorithm for factoring composite integers (that are not
 * pure powers) reduces the factoring problem to an order-finding problem.
 * Ekerå-Håstad's algorithm for factoring RSA integers reduces the factoring
 * problem to a short discrete logarithm problem. Both the general integer 
 * factoring problem and the RSA integer factoring problems are hence covered.
 */

/*!
 * \defgroup  linear_distribution Linear probability distributions
 * \ingroup   distribution
 *
 * \brief     A module for linear probability distributions.
 *
 * Linear probability distributions are used to represent the probability
 * distributions induced by Shor's order-finding algorithm [1], by Seifert's
 * algorithm for computing order with tradeoffs [2], by Ekerå's algorithm for
 * computing short discrete logarithms [3], and by Ekerå-Håstad's algorithm for
 * computing short discrete logarithms with tradeoffs [4, 5]. The algorithm
 * of Ekerå and Håstad also factors RSA integers by classically reducing the RSA
 * integer factoring problem to a short discrete logarithm problem in a group of
 * unknown order and solving this problem quantumly. 
 * 
 * [1] Shor, P.W.: Polynomial-time algorithms for prime factorization and 
 * discrete logarithms on a quantum computer. In: SIAM Journal on Scientific 
 * Computing (SISC), volume 26(5), pp. 1484 (1997). 
 * 
 * [2] Seifert, J.-P.: Using fewer Qubits in Shor's factorization algorithm via 
 * simultaneous Diophantine approximation. In: CT-RSA 2001, Spring LNCS 2020, 
 * pp. 319-327 (2001). DOI: https://doi.org/10.1007/3-540-45353-9_24.

 * [3] Ekerå, M.: Modifying Shor's algorithm to compute short discrete 
 * logarithms. In: IACR ePrint Archive, 2016/1128.
 * 
 * [4] Ekerå, M. and Håstad, J.: Quantum algorithms for computing short discrete
 * logarithms and factor RSA integers. In: PQCrypto 2017, Springer LNCS 10346,
 * pp. 347-363 (2017). DOI: https://doi.org/10.1007/978-3-319-59879-6_20.
 * 
 * [5] Ekerå, M.: On post-processing in the quantum algorithm for computing 
 * short discrete logarithms. In: IACR ePrint Archive, 2017/1122.
 * 
 * Linear probability distributions may also be used to represent the marginal
 * distributions of two-dimensional probability distributions, collapsed 
 * diagonal distributions, and traces in diagonal distributions for fixed 
 * offsets from the optimal alpha_d for given alpha_r.
 * 
 * Note that Shor's algorithm for factoring composite integers (that are not
 * pure powers) reduces the factoring problem to an order-finding problem.
 * Ekerå-Håstad's algorithm for factoring RSA integers reduces the factoring
 * problem to a short discrete logarithm problem. Both the general integer 
 * factoring problem and the RSA integer factoring problems are hence covered.
 */

#ifndef LINEAR_DISTRIBUTION_H
#define LINEAR_DISTRIBUTION_H

#include "linear_distribution_slice.h"
#include "distribution.h"
#include "random.h"
#include "parameters.h"
#include "common.h"

#include <mpfr.h>
#include <stdint.h>

/*!
 * \brief   An enumeration of flags indicating if the distribution if for a
 *          discrete logarithm d or an order r.
 */
enum {
  /*!
   * \brief   A flag indicating that the linear distribution is over alpha_d.
   */
  LINEAR_DISTRIBUTION_FLAG_D = 0,

  /*!
   * \brief   A flag indicating that the linear distribution is over alpha_r.
   */
  LINEAR_DISTRIBUTION_FLAG_R = 1
};

/*!
 * \brief   An enumeration of flags indicating if the distribution was computed
 *          by means of numerical integration or collapsed from a
 *          two-dimensional distribution.
 */
enum {
  /*!
   * \brief   A flag indicating that the distribution was constructed by means
   *          of numerical integration.
   */
  LINEAR_DISTRIBUTION_FLAG_COMPUTED = 0,

  /*!
   * \brief   A flag indicating that the distribution was constructed by
   *          collapsing a two-dimensional distribution to a marginal
   *          distribution.
   */
  LINEAR_DISTRIBUTION_FLAG_COLLAPSED = 2
};

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
  LINEAR_DISTRIBUTION_FLAG_DETERMINISTIC = 0,

  /*!
   * \brief   A flag indicating that the distribution was constructed by
   *          setting a minimal d or r = 2^(m - 1) + 1.
   */
  LINEAR_DISTRIBUTION_FLAG_MINIMAL = 4,

  /*!
   * \brief   A flag indicating that the distribution was constructed by
   *          setting a maximal d or r = 2^m - 1.
   */
  LINEAR_DISTRIBUTION_FLAG_MAXIMAL = 8,

  /*!
   * \brief   A flag indicating that the distribution was constructed by
   *          randomly selecting d or r on (2^(m-1), 2^m).
   *
   * For more information, see parameters_selection_random_d_or_r().
   */
  LINEAR_DISTRIBUTION_FLAG_RANDOM = 16,

  /*!
   * \brief   A flag indicating that the distribution was constructed by
   *          explicitly setting d or r to a value provided by the caller.
   */
  LINEAR_DISTRIBUTION_FLAG_EXPLICIT = 32,

  /*!
   * \brief   A flag indicating that the linear distribution was constructed by
   *          setting maximal d = 2^m - 1 where m = n / 2 - 1 for n the length
   *          in bits of an RSA integer.
   */
  LINEAR_DISTRIBUTION_FLAG_RSA = 64
};

/*!
 * \brief   A data structure representing a linear probability distribution.
 *
 * \ingroup linear_distribution
 */
typedef struct {
  /*! \brief  A vector of pointers to slices. */
  Linear_Distribution_Slice ** slices;

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

  /*! \brief The distribution flags. */
  uint32_t flags;

  /*!
   * \brief   The precision using which the distribution was computed.
   *
   * \remark  Incorrect precisions may be reported for distributions imported
   *          and exported using binaries compiled to use different precision
   *          levels, especially if the distribution is modified by computing
   *          or recomputing slices inbetween import and export operations. */
  uint32_t precision;
} Linear_Distribution;

/*!
 * \name Memory allocation
 * \{
 */

/*!
 * \brief   Allocates memory for a distribution.
 *
 * \return  The pointer to the memory allocated for the distribution.
 */
Linear_Distribution * linear_distribution_alloc();

/*!
 * \brief   Deallocates the memory previously allocated for a distribution.
 *
 * \param[in, out] distribution   A pointer to a pointer to the memory allocated
 *                                for the distribution. This function assigns
 *                                NULL to the pointer as it deallocates the
 *                                memory to which it refers.
 */
void linear_distribution_dealloc(
  Linear_Distribution ** const distribution);

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
void linear_distribution_init(
  Linear_Distribution * const distribution,
  const Parameters * const parameters,
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
void linear_distribution_init_copy_scale(
  Linear_Distribution * const dst_distribution,
  const Linear_Distribution * const src_distribution,
  const uint32_t dimension);

/*!
 * \brief   Clears a distribution.
 *
 * This function deallocates memory for slices in the distribution.
 *
 * \param[in, out] distribution   The distribution to clear.
 */
void linear_distribution_clear(
  Linear_Distribution * const distribution);

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
void linear_distribution_init_import(
  Linear_Distribution * const distribution,
  FILE * const file);

/*!
 * \brief   Exports a probability distribution to file.
 *
 * \param[in, out] distribution   The distribution to export to file.
 * \param[in, out] file           The file to which to export the distribution.
 */
void linear_distribution_export(
  const Linear_Distribution * const distribution,
  FILE * const file);

/*!
 * \brief   Exports information about a probability distribution to file.
 *
 * \param[in] distribution        The distribution to export to file.
 * \param[in, out] file           The file to which to export information.
 * \param[in] export_slices       Export verbose information on slices.
 * \param[in] export_regions      Export verbose information on regions.
 */
void linear_distribution_export_info(
  const Linear_Distribution * const distribution,
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
void linear_distribution_init_bcast_recv(
  Linear_Distribution * const distribution,
  const int root);

/*!
 * \brief   Broadcasts a distribution to all other nodes.
 *
 * \param[in] distribution        The distribution to broadcast.
 * \param[in] root                The rank of the root node.
 */
void linear_distribution_bcast_send(
  const Linear_Distribution * const distribution,
  const int root);

/*!
 * \}
 */

/*!
 * \name Collapsing
 * \{
 */

/*!
 * \brief   Initializes a distribution by collapsing a two-dimensional
 *          distribution to a marginal distribution in alpha_d.
 *
 * \param[in, out] dst    The distribution to initialize.
 * \param[in] src         The distribution to collapse.
 */
void linear_distribution_init_collapse_d(
  Linear_Distribution * const dst,
  const Distribution * const src);

/*!
 * \brief   Initializes a distribution by collapsing a two-dimensional
 *          distribution to a marginal distribution in alpha_r.
 *
 * \param[in, out] dst    The distribution to initialize.
 * \param[in] src         The distribution to collapse.
 */
void linear_distribution_init_collapse_r(
  Linear_Distribution * const dst,
  const Distribution * const src);

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
void linear_distribution_insert_slice(
  Linear_Distribution * const distribution,
  Linear_Distribution_Slice * const slice);

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
void linear_distribution_sort_slices(
  Linear_Distribution * const distribution);

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
const Linear_Distribution_Slice * linear_distribution_sample_slice(
  const Linear_Distribution * const distribution,
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
 * \param[in, out] min_log_alpha  The minimum signed logarithmic alpha.
 * \param[in, out] max_log_alpha  The maximum signed logarithmic alpha.
 *
 * \return  Returns #TRUE if a region was successfully sampled, #FALSE if the
 *          region sampled is outside the range of the distribution.
 */
bool linear_distribution_sample_region(
  const Linear_Distribution * const distribution,
  Random_State * const random_state,
  double * const min_log_alpha,
  double * const max_log_alpha);

/*!
 * \brief   Samples an approximate argument alpha from the probability
 *          distribution.
 *
 * The argument alpha is coarsely sampled, and is not guaranteed to be
 * admissible.
 *
 * \param[in] distribution        The distribution to sample.
 * \param[in, out] random_state   The random state to use when sampling.
 * \param[in, out] alpha          The argument alpha.
 *
 * \return  Returns #TRUE if the argument alpha was successfully sampled, #FALSE
 *          if the argument sampled is outside the range of the distribution.
 */
bool linear_distribution_sample_approximate_alpha(
  const Linear_Distribution * const distribution,
  Random_State * const random_state,
  mpfr_t alpha);

/*!
 * \brief   Samples an argument alpha from the probability distribution.
 *
 * The argument alpha is sampled with high resolution, and is guaranteed to be
 * admissible.
 *
 * To map output from this function to integers (j, k) for logarithms, or to
 * integers j for orders, depending on the distribution type, see the functions
 * sample_j_from_alpha_r() and sample_j_k_from_alpha_d(), respectively.
 *
 * \param[in] distribution        The distribution to sample.
 * \param[in, out] random_state   The random state to use when sampling.
 * \param[in, out] alpha          The argument alpha.
 *
 * \return  Returns #TRUE if the argument alpha was successfully sampled, #FALSE
 *          if the argument sampled is outside the range of the distribution.
 */
bool linear_distribution_sample_alpha(
  const Linear_Distribution * const distribution,
  Random_State * const random_state,
  mpz_t alpha);

/*!
 * \}
 */

#endif /* LINEAR_DISTRIBUTION_H */
