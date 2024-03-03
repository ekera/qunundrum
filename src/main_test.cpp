/*!
 * \file    main_test.cpp
 * \ingroup test_exe
 *
 * \brief   The definition of the main entry point to the test executable.
 */

/*!
 * \defgroup test_exe The test executable
 * \ingroup  executables
 * \ingroup  unit_tests
 *
 * \brief    A module for the test executable.
 */

#include "test/test_diagonal_distribution.h"
#include "test/test_diagonal_probability.h"
#include "test/test_keccak.h"
#include "test/test_keccak_random.h"
#include "test/test_linear_distribution.h"
#include "test/test_linear_probability.h"
#include "test/test_math.h"
#include "test/test_probability.h"
#include "test/test_random.h"
#include "test/test_sample.h"
#include "test/test_tau_volume_quotient.h"

#include "common.h"

#include <mpfr.h>

/*!
 * \brief The main entry point to the test executable.
 *
 * \param[in, out] argc   The arguments count.
 * \param[in, out] argv   The arguments vector.
 *
 * \return Zero upon successful execution, non-zero otherwise.
 */
int main(int argc, char ** argv) {
  (void)argc; /* Not used. */
  (void)argv; /* Not used. */

  mpfr_set_default_prec(PRECISION);

  test_math();

  test_keccak();
  test_keccak_random();
  test_random();

  test_probability();
  test_linear_probability();
  test_diagonal_probability();

  test_linear_distribution();
  test_diagonal_distribution();

  test_sample();

  test_tau_volume_quotient();
}
