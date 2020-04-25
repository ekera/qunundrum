/*!
 * \file    errors.c
 * \ingroup errors
 *
 * \brief The definition of functions for error handling.
 */

#include "errors.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

void critical(const char * msg, ...) {
  va_list args;
  va_start(args, msg);

  fprintf(stderr, "Error: ");
  vfprintf(stderr, msg, args);
  fprintf(stderr, "\n");
  fflush(stderr);

  va_end(args);

  exit(-1);
}
