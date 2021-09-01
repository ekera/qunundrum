/*!
 * \file    debug_common.c
 * \ingroup debug_trace_print
 *
 * \brief   The definition of functions used to generate debug trace printouts.
 */

#include "debug_common.h"

#include <stdint.h>
#include <stdio.h>

void debug_print_buffer(
  const uint8_t * const buffer,
  const uint32_t length)
{
  for (uint32_t i = 0; i < length; i++) {
    if ((i > 0) && ((i % 16) == 0)) {
      printf("\n");
    } else if ((i > 0) && ((i % 8) == 0)) {
      printf(" ");
    }

    printf("%.2x ", buffer[i]);
  }

  printf("\n");
}
