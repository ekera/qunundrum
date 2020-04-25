/*!
 * \file    test/test_sample.h
 * \ingroup unit_tests_sample_distribution
 * 
 * \brief   The declaration of unit tests for the \ref sample_distribution 
 *          module.
 */

/*!
 * \defgroup unit_tests_sample_distribution \
 *              Unit tests for sampling distributions
 * \ingroup  unit_tests
 * \ingroup  sample_distribution
 * 
 * \brief    A module for unit tests for sampling probability distributions.
 */

#ifndef TEST_SAMPLE_H
#define TEST_SAMPLE_H

/*!
 * \brief   Executes unit tests for test_sample_approximate_alpha_from_region().
 */
void test_sample_approximate_alpha_from_region();

/*!
 * \brief   Executes unit tests for test_sample_alpha_from_region().
 */
void test_sample_alpha_from_region();

/*!
 * \brief   Executes unit tests for test_sample_j_from_alpha_r().
 */
void test_sample_j_from_alpha_r();

/*!
 * \brief   Executes unit tests for test_sample_j_k_from_alpha_d().
 */
void test_sample_j_k_from_alpha_d();

/*!
 * \brief   Executes unit tests for test_sample_j_k_from_alpha_d_r().
 */
void test_sample_j_k_from_alpha_d_r();

/*!
 * \brief   Executes unit tests for sample_j_k_from_alpha_d_r().
 */
void test_sample_j_k_from_diagonal_alpha_d_r();

/*!
 * \brief   Executes all unit tests for the \ref sample_distribution module.
 */
void test_sample();

#endif /* TEST_SAMPLE_H */
