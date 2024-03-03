/*!
 * \file    test/test_diagonal_distribution.h
 * \ingroup unit_tests_diagonal_distribution
 *
 * \brief   The declaration of unit tests for diagonal distributions in the
 *          \ref diagonal_distribution module and sub-modules.
 */

/*!
 * \defgroup unit_tests_diagonal_distribution \
 *           Unit tests for diagonal distributions
 * \ingroup  unit_tests
 * \ingroup  diagonal_distribution
 *
 * \brief    A module for unit tests for diagonal distributions.
 */

#ifndef TEST_DIAGONAL_DISTRIBUTION_H
#define TEST_DIAGONAL_DISTRIBUTION_H

/*!
 * \brief   Executes unit tests for diagonal distribution slices when d and r
 *          are deterministically selected.
 */
void test_diagonal_distribution_slice_det();

/*!
 * \brief   Executes all unit tests for the \ref diagonal_distribution module.
 */
void test_diagonal_distribution();

#endif /* TEST_DIAGONAL_DISTRIBUTION_H */
