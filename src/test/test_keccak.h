/*!
 * \file    test/test_keccak.h
 * \ingroup unit_tests_keccak
 * 
 * \brief   The declaration of unit tests for \ref keccak module.
 */

/*!
 * \defgroup unit_tests_keccak Unit tests for the Keccak-f permutation
 * \ingroup  unit_tests
 * \ingroup  keccak
 * 
 * \brief    A module for unit tests for the Keccak-f permutation.
 */

#ifndef TEST_KECCAK_H
#define TEST_KECCAK_H

#include "../keccak.h"

#include <stdint.h>

/*!
 * \brief  The state resulting from applying keccak_f() once to the zero state.
 * 
 * This KAT vector is copied from KeccakF-1600-IntermediateValues.txt.
 */
static const uint64_t KECCAK_EXPECTED_FIRST[KECCAK_LANE_COUNT] = {
  0xF1258F7940E1DDE7ULL, 0x84D5CCF933C0478AULL, 0xD598261EA65AA9EEULL,
  0xBD1547306F80494DULL, 0x8B284E056253D057ULL, 0xFF97A42D7F8E6FD4ULL,
  0x90FEE5A0A44647C4ULL, 0x8C5BDA0CD6192E76ULL, 0xAD30A6F71B19059CULL,
  0x30935AB7D08FFC64ULL, 0xEB5AA93F2317D635ULL, 0xA9A6E6260D712103ULL,
  0x81A57C16DBCF555FULL, 0x43B831CD0347C826ULL, 0x01F22F1A11A5569FULL,
  0x05E5635A21D9AE61ULL, 0x64BEFEF28CC970F2ULL, 0x613670957BC46611ULL,
  0xB87C5A554FD00ECBULL, 0x8C3EE88A1CCF32C8ULL, 0x940C7922AE3A2614ULL,
  0x1841F924A2C509E4ULL, 0x16F53526E70465C2ULL, 0x75F644E97F30A13BULL,
  0xEAF1FF7B5CECA249ULL
};

/*!
 * \brief  The state resulting from applying keccak_f() twice to the zero state.
 * 
 * This KAT vector is copied from KeccakF-1600-IntermediateValues.txt.
 */
static const uint64_t KECCAK_EXPECTED_SECOND[KECCAK_LANE_COUNT] = {
  0x2D5C954DF96ECB3CULL, 0x6A332CD07057B56DULL, 0x093D8D1270D76B6CULL,
  0x8A20D9B25569D094ULL, 0x4F9C4F99E5E7F156ULL, 0xF957B9A2DA65FB38ULL,
  0x85773DAE1275AF0DULL, 0xFAF4F247C3D810F7ULL, 0x1F1B9EE6F79A8759ULL,
  0xE4FECC0FEE98B425ULL, 0x68CE61B6B9CE68A1ULL, 0xDEEA66C4BA8F974FULL,
  0x33C43D836EAFB1F5ULL, 0xE00654042719DBD9ULL, 0x7CF8A9F009831265ULL,
  0xFD5449A6BF174743ULL, 0x97DDAD33D8994B40ULL, 0x48EAD5FC5D0BE774ULL,
  0xE3B8C8EE55B7B03CULL, 0x91A0226E649E42E9ULL, 0x900E3129E7BADD7BULL,
  0x202A9EC5FAA3CCE8ULL, 0x5B3402464E1C3DB6ULL, 0x609F4E62A44C1059ULL,
  0x20D06CD26A8FBF5CULL
};

/*!
 * \brief   Executes unit tests for the Keccak-f permutation keccak_f().
 */
void test_keccak_f();

/*!
 * \brief   Executes all unit tests for the \ref keccak module.
 */
void test_keccak();

#endif /* TEST_KECCAK_H */