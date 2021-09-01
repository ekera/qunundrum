/*!
 * \file    test/test_probability.h
 * \ingroup unit_tests_two_dimensional_distribution
 *
 * \brief   The declaration of unit tests for the probability functions in the
 *          \ref two_dimensional_distribution_probability module.
 */

/*!
 * \defgroup unit_tests_two_dimensional_distribution_probability \
 *              Unit tests for two-dimensional probability functions
 * \ingroup  unit_tests
 * \ingroup  two_dimensional_distribution_probability
 *
 * \brief    A module for unit tests for two-dimensional probability
 *           functions.
 */

#ifndef TEST_PROBABILITY_H
#define TEST_PROBABILITY_H

/*!
 * \brief   Executes unit tests for probability_approx() with heuristically
 *          selected sigma via KAT vectors.
 */
void test_probability_approx_heuristic_sigma_kat();

/*!
 * \brief   Execute all unit tests for the
 *          \ref two_dimensional_distribution_probability module.
 */
void test_probability();

#endif /* TEST_PROBABILITY_H */
