/*!
 * \file    test/test_diagonal_probability.h
 * \ingroup unit_tests_diagonal_distribution_probability
 * 
 * \brief   The declaration of unit tests for the probability functions in the 
 *          \ref diagonal_distribution_probability module.
 */

/*!
 * \defgroup unit_tests_diagonal_distribution_probability \
 *              Unit tests for diagonal probability functions
 * \ingroup  unit_tests
 * \ingroup  diagonal_distribution_probability
 * 
 * \brief    A module for unit tests for diagonal probability functions.
 */

#ifndef TEST_DIAGONAL_PROBABILITY_H
#define TEST_DIAGONAL_PROBABILITY_H

/*!
 * \brief   Executes unit tests for diagonal_probability() via KAT 
 *          vectors.
 */
void test_diagonal_probability_approx_kat();

/*!
 * \brief   Execute all unit tests for the 
 *          \ref diagonal_distribution_probability module.
 */
void test_diagonal_probability();

#endif /* TEST_DIAGONAL_PROBABILITY_H */
