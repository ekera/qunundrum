/*!
 * \file    test/test_linear_distribution.h
 * \ingroup unit_tests_linear_distribution
 *
 * \brief   The declaration of unit tests for linear distributions in the
 *          \ref linear_distribution module and sub-modules.
 */

/*!
 * \defgroup unit_tests_linear_distribution Unit tests for linear distributions
 * \ingroup  unit_tests
 * \ingroup  linear_distribution
 *
 * \brief    A module for unit tests for linear distributions.
 */

#ifndef TEST_LINEAR_DISTRIBUTION_H
#define TEST_LINEAR_DISTRIBUTION_H

/*!
 * \brief   Executes unit tests for linear distribution slices when d is
 *          deterministically selected.
 */
void test_linear_distribution_slice_det_d();

/*!
 * \brief   Executes unit tests for linear distribution slices when r is
 *          deterministically selected.
 */
void test_linear_distribution_slice_det_r();

/*!
 * \brief   Executes all unit tests for the \ref linear_distribution module.
 */
void test_linear_distribution();

#endif /* TEST_LINEAR_DISTRIBUTION_H */
