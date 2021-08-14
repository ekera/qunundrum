/*!
 * \file    test/test_math.h
 * \ingroup unit_tests_math
 * 
 * \brief   The declaration of unit tests for the basic mathematical functions 
 *          in the \ref math module.
 */

/*!
 * \defgroup unit_tests_math Unit tests for mathematical functions
 * \ingroup  unit_tests
 * \ingroup  math
 * 
 * \brief    A module for unit tests for mathematical functions.
 */

#ifndef TEST_MATH_H
#define TEST_MATH_H

/*!
 * \brief   Executes unit tests for the abs_i(), abs_d() and abs_ld() functions.
 *
 * \note    This test function is limited to rudimentary sanity checks.
 */
void test_abs();

/*!
 * \brief   Executes unit tests for the sgn_i(), sgn_d() and sgn_ld() functions.
 *
 * \note    This test function is limited to rudimentary sanity checks.
 */
void test_sgn();

/*!
 * \brief   Executes unit tests for the min_i(), min_ui(), min_d() and min_ld()
 *          functions.
 *
 * \note    This test function is limited to rudimentary sanity checks.
 */
void test_min();

/*!
 * \brief   Executes unit tests for the max_i(), max_ui(), max_d() and max_ld()
 *          functions.
 *
 * \note    This test function is limited to rudimentary sanity checks.
 */
void test_max();

/*!
 * \brief   Executes unit tests for the mod_reduce() function.
 *
 * \note    This test function is limited to rudimentary sanity checks.
 */
void test_math_mod_reduce();

/*!
 * \brief   Executes unit tests for the kappa() function.
 *
 * \note    This test function is limited to rudimentary sanity checks.
 */
void test_kappa();

/*!
 * \brief   Executes unit tests for the is_pow2() function.
 *
 * \note    This test function is limited to rudimentary sanity checks.
 */
void test_is_pow2();

/*!
 * \brief   Executes all unit tests for the \ref math module.
 */
void test_math();

#endif /* TEST_MATH_H */
