/*!
 * \file    test/test_linear_distribution.cpp
 * \ingroup unit_tests_linear_distribution
 *
 * \brief   The definition of unit tests for linear distributions in the
 *          \ref linear_distribution module and sub-modules.
 */

#include "test_linear_distribution.h"

#include "test_common.h"

#include "../errors.h"
#include "../linear_distribution.h"
#include "../linear_distribution_slice.h"
#include "../math.h"
#include "../parameters.h"
#include "../parameters_selection.h"

#include <gmp.h>

#include <stdint.h>
#include <stdio.h>

/*!
 * \brief   The dimension for which slices are computed for the tests.
 */
#define TEST_SLICE_DIMENSION  2048

/*!
 * \brief   The dimension to which computed slices are scaled for the tests.
 */
#define TEST_SLICE_SCALED_DIMENSION  256

void test_linear_distribution_slice_det_d() {
  printf("Testing linear_distribution_slice_det_d...\n");

  /* Initialize a random state. */
  Random_State random_state;
  random_init(&random_state);

  /* Setup constants. */
  const uint32_t m = 128;
  const uint32_t s = 1;
  const uint32_t t = 30;

  /* Setup the distribution parameters. */
  Parameters parameters;
  parameters_init(&parameters);

  {
    mpz_t d;
    mpz_init(d);

    mpz_t r;
    mpz_init(r);

    parameters_selection_deterministic_d_r(d, r, m);
    parameters_explicit_m_s(&parameters, d, r, m, s, t);

    /* Clear memory. */
    mpz_clear(d);
    mpz_clear(r);
  }

  /* Setup linear slices. */
  Linear_Distribution_Slice slice;
  linear_distribution_slice_init(&slice, TEST_SLICE_DIMENSION);

  Linear_Distribution_Slice scaled_slice;
  linear_distribution_slice_init(&scaled_slice, TEST_SLICE_SCALED_DIMENSION);

  /* Tables of expected probabilities. */
  const int32_t offsets[16] = {
    -5,
    -4,
    -3,
    -2,
    -1,
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10
  };

  const long double expected_probabilities[16] = {
    /* From Mathematica: NIntegrate with flags: MaxRecursion->400,
     *    AccuracyGoal->100, WorkingPrecision->128. With truncation. */
    0.0213302826690752331655749215386265799772932747823383596421465804734290714,
      /* [2^(m-5), 2^(m-4)] */
    0.0420576665146783784417293062314817641731320267519831675756683442638383806,
      /* [2^(m-4), 2^(m-3)] */
    0.0794739598552841691767122216585393553672987347723427928770629180786529151,
      /* [2^(m-3), 2^(m-2)] */
    0.1270683539258403722598134507301931804040640413057477695563250535313967743,
      /* [2^(m-2), 2^(m-1)] */
    0.1123141473052097413758361608634787146846746520408801672596436570997174342,
      /* [2^(m-1), 2^m] */
    0.0474329958267535240321051033248911090290445537190744654773640404055915158,
      /* [2^m, 2^(m+1)] */
    0.0243420515407842744484605366781757888593658571505950262454244329669597813,
      /* [2^(m+1), 2^(m+2)] */
    0.0122669383037106757672237067167467742607490758746441612884936554414410219,
      /* [2^(m+2), 2^(m+3)] */
    0.0061462977686566458801273700461281001171955577179703184271225106629400052,
      /* [2^(m+3), 2^(m+4)] */
    0.0030747829961898171349306827232805979318482966192337349240037997222295410,
      /* [2^(m+4), 2^(m+5)] */
    0.0015375967597999937769230085809629277095264861132239557653786153706396106,
      /* [2^(m+5), 2^(m+6)] */
    0.0007688240691602858847495852270694401085800060442789012562081147438420469,
      /* [2^(m+6), 2^(m+7)] */
    0.0003844152467264679931603625004022327590014413821934008934715288060114952,
      /* [2^(m+7), 2^(m+8)] */

    /* From Mathematica: NIntegrate with flags: MaxRecursion->400,
     *    AccuracyGoal->128, WorkingPrecision->128. With truncation. */
    0.0001922080249124473806002123883887596973554682978806125525283932544047340,
      /* [2^(m+8), 2^(m+9)] */
    0.0000961040626508418800326012912687449444393080220590186190377708994639885,
      /* [2^(m+9), 2^(m+10)] */
    0.0000480520375997784187929535416498227530711211962492982194256159256676091
      /* [2^(m+10), 2^(m+11)] */
  };

  /* Declare temporary variables used below. */
  double tmp_min_log_alpha;
  double tmp_max_log_alpha;

  /* Compute slices, */
  for (uint32_t i = 0; i < 16; i++) {
    int32_t min_log_alpha = m + offsets[i];
    int32_t max_log_alpha = m + offsets[i] + 1;

    /* Positive. */
    printf(" Computing slice for min_log_alpha_d = %d...\n", min_log_alpha);

    linear_distribution_slice_compute_richardson(
      &slice,
      &parameters,
      LINEAR_DISTRIBUTION_SLICE_COMPUTE_TARGET_D,
      min_log_alpha);

    /* Sanity checks. */
    if (TRUE != test_cmp_ld(slice.total_probability,
                            expected_probabilities[i]))
    {
      critical("test_linear_distribution_slice_det_d(): "
        "Incorrect probability.");
    }

    if (0 != slice.total_error) {
      critical("test_linear_distribution_slice_det_d(): "
        "Expected a total error of zero.");
    }

    if (min_log_alpha != slice.min_log_alpha) {
      critical("test_linear_distribution_slice_det_d(): "
        "Incorrect min_log_alpha.");
    }

    linear_distribution_slice_coordinates(
      &slice, &tmp_min_log_alpha, &tmp_max_log_alpha);

    if ((min_log_alpha != tmp_min_log_alpha) ||
        (max_log_alpha != tmp_max_log_alpha))
    {
      critical("test_linear_distribution_slice_det_d(): "
        "Incorrect min_log_alpha and max_log_alpha returned.");
    }

    linear_distribution_slice_sample_region(
      &slice,
      &random_state,
      &tmp_min_log_alpha,
      &tmp_max_log_alpha);

    if ((abs_d(tmp_min_log_alpha) < abs_i(min_log_alpha)) ||
        (abs_d(tmp_max_log_alpha) > abs_i(max_log_alpha)) ||
        (sgn_d(tmp_min_log_alpha) != sgn_i(min_log_alpha)) ||
        (sgn_d(tmp_max_log_alpha) != sgn_i(max_log_alpha)))
    {
      critical("test_linear_distribution_slice_det_d(): "
        "Failed to sample a region from the slice.");
    }


    /* Negative. */
    min_log_alpha *= -1;
    max_log_alpha *= -1;

    printf(" Computing slice for min_log_alpha_d = %d...\n", min_log_alpha);

    linear_distribution_slice_compute_richardson(
      &slice,
      &parameters,
      LINEAR_DISTRIBUTION_SLICE_COMPUTE_TARGET_D,
      min_log_alpha);

    /* Sanity checks. */
    if (TRUE != test_cmp_ld(slice.total_probability,
                            expected_probabilities[i]))
    {
      critical("test_linear_distribution_slice_det_d(): "
        "Incorrect probability.");
    }

    if (0 != slice.total_error) {
      critical("test_linear_distribution_slice_det_d(): "
        "Expected a total error of zero.");
    }

    if (min_log_alpha != slice.min_log_alpha) {
      critical("test_linear_distribution_slice_det_d(): "
        "Incorrect min_log_alpha.");
    }

    linear_distribution_slice_coordinates(
      &slice, &tmp_min_log_alpha, &tmp_max_log_alpha);

    if ((min_log_alpha != tmp_min_log_alpha) ||
        (max_log_alpha != tmp_max_log_alpha))
    {
      critical("test_linear_distribution_slice_det_d(): "
        "Incorrect min_log_alpha and max_log_alpha returned.");
    }

    linear_distribution_slice_sample_region(
      &slice,
      &random_state,
      &tmp_min_log_alpha,
      &tmp_max_log_alpha);

    if ((abs_d(tmp_min_log_alpha) < abs_i(min_log_alpha)) ||
        (abs_d(tmp_max_log_alpha) > abs_i(max_log_alpha)) ||
        (sgn_d(tmp_min_log_alpha) != sgn_i(min_log_alpha)) ||
        (sgn_d(tmp_max_log_alpha) != sgn_i(max_log_alpha)))

    {
      critical("test_linear_distribution_slice_det_d(): "
        "Failed to sample a region from the slice.");
    }

    /* Test scaling. */
    linear_distribution_slice_copy_scale(&scaled_slice, &slice);

     if ((scaled_slice.min_log_alpha != slice.min_log_alpha) ||
        (scaled_slice.total_error) != (slice.total_error) ||
        (scaled_slice.total_probability) != (slice.total_probability))
    {
      critical("test_linear_distribution_slice_det_d(): "
        "Failed to correctly scale the slice.");
    }

    long double tmp_sum = 0;

    for (uint32_t i = 0; i < TEST_SLICE_SCALED_DIMENSION; i++) {
      tmp_sum += scaled_slice.norm_vector[i];

      long double tmp_inner_sum = 0;

      for (uint32_t j = 0;
        j < (TEST_SLICE_DIMENSION / TEST_SLICE_SCALED_DIMENSION); j++)
      {
        tmp_inner_sum +=
          slice.norm_vector[i * (TEST_SLICE_DIMENSION / \
            TEST_SLICE_SCALED_DIMENSION) + j];
      }

      if (tmp_inner_sum != scaled_slice.norm_vector[i]) {
        critical("test_linear_distribution_slice_det_d(): "
        "Failed to correctly scale the slice.");
      }
    }

    if (TRUE != test_cmp_ld(slice.total_probability, tmp_sum)) {
      critical("test_diagonal_distribution_slice_det_d(): "
        "Incorrect total probability.");
    }
  }

  /* Clear memory. */
  random_close(&random_state);
  parameters_clear(&parameters);
  linear_distribution_slice_clear(&slice);
  linear_distribution_slice_clear(&scaled_slice);
}

void test_linear_distribution_slice_det_r() {
  printf("Testing linear_distribution_slice_det_r...\n");

  /* Initialize a random state. */
  Random_State random_state;
  random_init(&random_state);

  /* Setup constants. */
  const uint32_t m = 128;
  const uint32_t s = 1;
  const uint32_t t = 30;

  /* Setup the distribution parameters. */
  Parameters parameters;
  parameters_init(&parameters);

  {
    mpz_t d;
    mpz_init(d);

    mpz_t r;
    mpz_init(r);

    parameters_selection_deterministic_d_r(d, r, m);
    parameters_explicit_m_s(&parameters, d, r, m, s, t);

    /* Clear memory. */
    mpz_clear(d);
    mpz_clear(r);
  }

  /* Setup linear slices. */
  Linear_Distribution_Slice slice;
  linear_distribution_slice_init(&slice, TEST_SLICE_DIMENSION);

  Linear_Distribution_Slice scaled_slice;
  linear_distribution_slice_init(&scaled_slice, TEST_SLICE_SCALED_DIMENSION);

  const int32_t offsets[16] = {
    -5,
    -4,
    -3,
    -2,
    -1,
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10
  };

  const long double expected_probabilities[16] = {
    /* From Mathematica, NIntegrate with flags: MaxRecursion->100,
     *    AccuracyGoal->128, WorkingPrecision->128. With truncation. */
    0.0323551566868607057987128637435863424470758464295070103141839652641850335,
      /* [2^(m-5), 2^(m-4)] */
    0.0631410235334714141067135492736508468495944310739498462161862778164728670,
      /* [2^(m-4), 2^(m-3)] */
    0.1144094675627533466623530034715829169829886308167174798581331710970881274,
      /* [2^(m-3), 2^(m-2)] */
    0.1528605577062922115704915499274528852227941442420105076794329273362429410,
      /* [2^(m-2), 2^(m-1)] */
    0.0560891354848479208165243201617735694796343558387728798428957156861911738,
      /* [2^(m-1), 2^m] */
    0.0235836799389394685265174126269091895831236664504133809659564538405832384,
      /* [2^m, 2^(m+1)] */
    0.0124518973150204485762065258224944096778615407823754086102417542721608297,
      /* [2^(m+1), 2^(m+2)] */
    0.0063644343005867865695820209506418126921739943954736524740486117634231883,
      /* [2^(m+2), 2^(m+3)] */
    0.0031562129969183507784989830740790833483145080277677521194518887834426412,
      /* [2^(m+3), 2^(m+4)] */
    0.0014850434739194293584320684379776103990120010665910384595281601067884592,
      /* [2^(m+4), 2^(m+5)] */
    0.0007641798136142296366270209096268393073949169114151517071530251990664484,
      /* [2^(m+5), 2^(m+6)] */
    0.0003777572838517225193305000923278374142932756014648348117470964043042342,
      /* [2^(m+6), 2^(m+7)] */
    0.0001891702832049387255524309144448721795974838470558931114033484951978106,
      /* [2^(m+7), 2^(m+8)] */
    0.0000948931857172883691464838630069342870920082407077185797947270267400545,
      /* [2^(m+8), 2^(m+9)] */
    0.0000474058469341121956593564853031255002486335477251926032407721754481573,
      /* [2^(m+9), 2^(m+10)] */
    0.0000236951428994529638224975620505152365939244258049075053002091394808050
      /* [2^(m+10), 2^(m+11)] */
  };

  /* Declare temporary variables used below. */
  double tmp_min_log_alpha;
  double tmp_max_log_alpha;

  /* Compute slices, */
  for (uint32_t i = 0; i < 16; i++) {
    int32_t min_log_alpha = m + offsets[i];
    int32_t max_log_alpha = m + offsets[i] + 1;

    /* Positive. */
    printf(" Computing slice for min_log_alpha_r = %d...\n", min_log_alpha);

    linear_distribution_slice_compute_richardson(
      &slice,
      &parameters,
      LINEAR_DISTRIBUTION_SLICE_COMPUTE_TARGET_R,
      min_log_alpha);

    /* Sanity checks. */
    if (TRUE != test_cmp_ld(slice.total_probability,
                            expected_probabilities[i],
                            (offsets[i] >= 10) ? 1e-4 : 1e-6))
    {
      critical("test_linear_distribution_slice_det_r(): "
        "Incorrect probability.");
    }

    if (0 != slice.total_error) {
      critical("test_linear_distribution_slice_det_r(): "
        "Expected a total error of zero.");
    }

    linear_distribution_slice_coordinates(
      &slice, &tmp_min_log_alpha, &tmp_max_log_alpha);

    if ((min_log_alpha != tmp_min_log_alpha) ||
        (max_log_alpha != tmp_max_log_alpha))
    {
      critical("test_linear_distribution_slice_det_r(): "
        "Incorrect min_log_alpha and max_log_alpha returned.");
    }

    linear_distribution_slice_sample_region(
      &slice,
      &random_state,
      &tmp_min_log_alpha,
      &tmp_max_log_alpha);

    if ((abs_d(tmp_min_log_alpha) < abs_i(min_log_alpha)) ||
        (abs_d(tmp_max_log_alpha) > abs_i(max_log_alpha)) ||
        (sgn_d(tmp_min_log_alpha) != sgn_i(min_log_alpha)) ||
        (sgn_d(tmp_max_log_alpha) != sgn_i(max_log_alpha)))

    {
      critical("test_linear_distribution_slice_det_r(): "
        "Failed to sample a region from the slice.");
    }

    /* Negative. */
    min_log_alpha *= -1;
    max_log_alpha *= -1;

    printf(" Computing slice for min_log_alpha_r = %d...\n", min_log_alpha);

    linear_distribution_slice_compute_richardson(
      &slice,
      &parameters,
      LINEAR_DISTRIBUTION_SLICE_COMPUTE_TARGET_R,
      min_log_alpha);

    /* Sanity checks. */
    if (TRUE != test_cmp_ld(slice.total_probability,
                            expected_probabilities[i],
                            (offsets[i] >= 10) ? 1e-4 : 1e-6))
    {
      critical("test_linear_distribution_slice_det_r(): "
        "Incorrect probability.");
    }

    if (0 != slice.total_error) {
      critical("test_linear_distribution_slice_det_r(): "
        "Expected a total error of zero.");
    }

    linear_distribution_slice_coordinates(
      &slice, &tmp_min_log_alpha, &tmp_max_log_alpha);

    if ((min_log_alpha != tmp_min_log_alpha) ||
        (max_log_alpha != tmp_max_log_alpha))
    {
      critical("test_linear_distribution_slice_det_r(): "
        "Incorrect min_log_alpha and max_log_alpha returned.");
    }

    linear_distribution_slice_sample_region(
      &slice,
      &random_state,
      &tmp_min_log_alpha,
      &tmp_max_log_alpha);

    if ((abs_d(tmp_min_log_alpha) < abs_i(min_log_alpha)) ||
        (abs_d(tmp_max_log_alpha) > abs_i(max_log_alpha)) ||
        (sgn_d(tmp_min_log_alpha) != sgn_i(min_log_alpha)) ||
        (sgn_d(tmp_max_log_alpha) != sgn_i(max_log_alpha)))

    {
      critical("test_linear_distribution_slice_det_r(): "
        "Failed to sample a region from the slice.");
    }

    /* Test scaling. */
    linear_distribution_slice_copy_scale(&scaled_slice, &slice);

     if ((scaled_slice.min_log_alpha != slice.min_log_alpha) ||
        (scaled_slice.total_error) != (slice.total_error) ||
        (scaled_slice.total_probability) != (slice.total_probability))
    {
      critical("test_linear_distribution_slice_det_r(): "
        "Failed to correctly scale the slice.");
    }

    long double tmp_sum = 0;

    for (uint32_t i = 0; i < TEST_SLICE_SCALED_DIMENSION; i++) {
      tmp_sum += scaled_slice.norm_vector[i];

      long double tmp_inner_sum = 0;

      for (uint32_t j = 0;
        j < (TEST_SLICE_DIMENSION / TEST_SLICE_SCALED_DIMENSION); j++)
      {
        tmp_inner_sum +=
          slice.norm_vector[i * (TEST_SLICE_DIMENSION / \
            TEST_SLICE_SCALED_DIMENSION) + j];
      }

      if (tmp_inner_sum != scaled_slice.norm_vector[i]) {
        critical("test_linear_distribution_slice_det_r(): "
        "Failed to correctly scale the slice.");
      }
    }

    if (TRUE != test_cmp_ld(slice.total_probability, tmp_sum)) {
      critical("test_diagonal_distribution_slice_det_r(): "
        "Incorrect total probability.");
    }
  }

  /* Clear memory. */
  random_close(&random_state);
  parameters_clear(&parameters);
  linear_distribution_slice_clear(&slice);
  linear_distribution_slice_clear(&scaled_slice);
}

void test_linear_distribution() {
  test_linear_distribution_slice_det_d();
  test_linear_distribution_slice_det_r();
}
