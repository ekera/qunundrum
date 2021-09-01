/*!
 * \file    test/test_random.h
 * \ingroup unit_tests_random
 *
 * \brief   The declaration of unit tests for the \ref random module.
 */

/*!
 * \defgroup unit_tests_random Unit tests for random number generation
 * \ingroup  unit_tests
 * \ingroup  random
 *
 * \brief    A module for unit tests for random number generation.
 */

#ifndef TEST_RANDOM_H
#define TEST_RANDOM_H

/*!
 * \brief   Executes unit tests for reading randomness from a device.
 */
void test_random_device();

/*!
 * \brief   Executes unit tests for reading randomness from a Keccak state.
 */
void test_random_keccak();

/*!
 * \brief   Executes unit tests for reading random pivots.
 */
void test_random_generate_pivot();

/*!
 * \brief   Executes unit tests for reading random multiprecision integers.
 */
void test_random_generate_mpz();

/*!
 * \brief   Executes unit tests for the \ref random module.
 */
void test_random();

#endif /* TEST_RANDOM_H */
