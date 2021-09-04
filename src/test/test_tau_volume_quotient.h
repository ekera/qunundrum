/*!
 * \file    test/test_tau_volume_quotient.h
 * \ingroup unit_tests_estimating_volume_quotients
 *
 * \brief   The declaration of unit tests for the
 *          \ref estimating_volume_quotients module.
 */

/*!
 * \defgroup unit_tests_estimating_volume_quotients \
 *              Unit tests for estimating volume quotients
 * \ingroup  unit_tests
 * \ingroup  estimating_volume_quotients
 *
 * \brief    A module for unit tests for estimating volume quotients.
 */

#ifndef TEST_TAU_VOLUME_QUOTIENT_H
#define TEST_TAU_VOLUME_QUOTIENT_H

/*!
 * \brief   Executes unit tests for tau_volume_quotient() via KAT vectors.
 */
void test_tau_volume_quotient_kat();

/*!
 * \brief   Execute all unit tests for the \ref estimating_volume_quotients
 *          module.
 */
void test_tau_volume_quotient();

#endif /* TEST_TAU_VOLUME_QUOTIENT_H */
