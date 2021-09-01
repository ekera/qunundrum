/*!
 * \file    test/test_linear_probability.h
 * \ingroup unit_tests_linear_distribution_probability
 *
 * \brief   The declaration of unit tests for the probability functions in the
 *          \ref linear_distribution_probability module.
 */

/*!
 * \defgroup unit_tests_linear_distribution_probability \
 *              Unit tests for linear probability functions
 * \ingroup  unit_tests
 * \ingroup  linear_distribution_probability
 *
 * \brief    A module for unit tests for linear probability functions.
 */

#ifndef TEST_LINEAR_PROBABILITY_H
#define TEST_LINEAR_PROBABILITY_H

/*!
 * \brief   Executes unit tests for linear_probability_d() via KAT vectors.
 */
void test_linear_probability_d_kat();

/*!
 * \brief   Executes unit tests for linear_probability_r() via KAT vectors.
 */
void test_linear_probability_r_kat();

/*!
 * \brief   Execute all unit tests for the \ref linear_distribution_probability
 *          module.
 */
void test_linear_probability();

#endif /* TEST_LINEAR_PROBABILITY_H */
