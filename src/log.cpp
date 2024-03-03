#include "log.h"

#include <time.h>

#include <stdio.h>

/*!
 * \brief   The length in bytes of the timestamp buffer.
 */
#define LOG_TIMESTAMP_BUFFER_LENGTH 256

void log_timestamp_fprintf(
  FILE * const file)
{
  time_t timestamp;
  struct tm* tm_info;

  timestamp = time(NULL);
  tm_info = localtime(&timestamp);

  char buffer[LOG_TIMESTAMP_BUFFER_LENGTH];
  strftime(
    buffer,
    LOG_TIMESTAMP_BUFFER_LENGTH,
    "%Y-%m-%d %H:%M:%S %Z",
    tm_info);

  fprintf(file, "# Timestamp: %s\n", buffer);
}