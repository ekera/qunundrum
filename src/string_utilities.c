/*!
 * \file    string_utilities.c
 * \ingroup string_utility
 *
 * \brief   The definition of utility functions for manipulating strings.
 */

#include "string_utilities.h"

#include "errors.h"

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

void safe_snprintf(char * str, size_t size, const char * format, ...)
{
  va_list args;

  va_start(args, format);
  int written = vsnprintf(str, size, format, args);
  va_end(args);

  if ((written < 0) || (written >= (int)size)) {
    critical("safe_snprintf(): "
      "Overflow detected when constructing the string.");
  }
}

void safe_strlcpy(char * dst, const char * src, size_t dst_size) {
  size_t src_len = strlen(src);

  if (dst_size <= src_len) {
    critical("safe_strlcpy(): Overflow detected when copying the string.");
  }

  memcpy(dst, src, (src_len + 1) * sizeof(char));
}

const char * truncate_path(const char * path) {
  for (const char * path_ptr = path; '\0' != *path_ptr; path_ptr++) {
    if ('/' == *path_ptr) {
      path = path_ptr + 1;
    }
  }

  return path;
}
