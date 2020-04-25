/*!
 * \file    debug_common.h
 * \ingroup debug_trace_print
 * 
 * \brief   The declaration of functions used to generate debug trace printouts.
 */

/*!
 * \defgroup debug_trace_print Debug printouts
 * \ingroup  unit_tests_debug_tracing
 * \ingroup  utility
 *
 * \brief A module for functions for debugging and tracing.
 */

#ifndef DEBUG_COMMON_H
#define DEBUG_COMMON_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief   Prints the contents of a buffer to the console.
 * 
 * \param[in] buffer    The buffer to print.
 * \param[in] length    The length in bytes of the buffer.
 */
void debug_print_buffer(
  const uint8_t * const buffer,
  const uint32_t length);

#ifdef __cplusplus
}
#endif

#endif /* DEBUG_COMMON_H */