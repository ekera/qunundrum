#include "log.h"

#include "common.h"
#include "errors.h"
#include "executables.h"
#include "string_utilities.h"

#include <stdint.h>
#include <stdio.h>

#include <unistd.h>
#include <sys/stat.h>

#include <time.h>

/*!
 * \brief   The length in bytes of the timestamp buffer.
 */
#define LOG_TIMESTAMP_BUFFER_LENGTH 256

/*!
 * \brief   The path to the device from which to read random data for use in
 *          generating unique log file names.
 */
static const char * const LOG_RANDOM_DEVICE = "/dev/urandom";

/*!
 * \brief   The length in bytes of the log file name suffix buffer.
 */
#define LOG_SUFFIX_BUFFER_LENGTH 32

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

/*!
 * \brief   Generates a log file name suffix of the form
 *          "YYYYMMDD-HHmmss±ZZZZ-XXXXXXXX".
 *
 * For further details on the format, see the documentation for log_open().
 *
 * \param[in, out] suffix   The buffer in which to store the suffix.
 * \param[in] length        The capacity of the buffer in characters.
 */
static void log_suffix_generate(
  char * const suffix,
  const size_t length)
{
  time_t timestamp = time(NULL);

  struct tm tm_info;
  localtime_r(&timestamp, &tm_info);

  char datetime[LOG_SUFFIX_BUFFER_LENGTH];
  strftime(
    datetime,
    LOG_SUFFIX_BUFFER_LENGTH,
    "%Y%m%d-%H%M%S",
    &tm_info);

  const long offset_seconds = tm_info.tm_gmtoff;
  const char offset_sign = (offset_seconds < 0) ? '-' : '+';
  const long offset_seconds_abs =
    (offset_seconds < 0) ? -offset_seconds : offset_seconds;
  const int offset_hours = (int)(offset_seconds_abs / 3600);
  const int offset_minutes = (int)((offset_seconds_abs / 60) % 60);

  /* Read a random 32-bit unsigned integer from the random device. */
  FILE * random_device = fopen(LOG_RANDOM_DEVICE, "rb");
  if (NULL == random_device) {
    critical("log_suffix_generate(): Failed to open \"%s\" for reading.",
      LOG_RANDOM_DEVICE);
  }

  uint32_t random_value;
  if (1 != fread(&random_value, sizeof(random_value), 1, random_device)) {
    critical("log_suffix_generate(): Failed to read from \"%s\".",
      LOG_RANDOM_DEVICE);
  }

  fclose(random_device);
  random_device = NULL;

  safe_snprintf(suffix, length, "%s%c%02d%02d-%08x",
    datetime, offset_sign, offset_hours, offset_minutes, random_value);
}

FILE * log_open(
  const char * const prefix)
{
  /* Generate the suffix once and cache it for the remainder of the process, so
   * that repeated calls reopen the same log file. */
  static bool log_setup = FALSE;

  static char suffix[LOG_SUFFIX_BUFFER_LENGTH] = { '\0' };

  if (FALSE == log_setup) {
    log_suffix_generate(suffix, LOG_SUFFIX_BUFFER_LENGTH);
  }

  /* Create the log directory if it does not exist. */
  if (0 != access(LOGS_DIRECTORY, F_OK)) {
    if (0 != mkdir(LOGS_DIRECTORY, DIRECTORY_PERMISSIONS)) {
      critical("log_open(): Failed to create the directory \"%s\".",
        LOGS_DIRECTORY);
    }
  }

  /* Setup the path to the log file. */
  char log_path[MAX_SIZE_PATH_BUFFER];
  safe_snprintf(log_path, MAX_SIZE_PATH_BUFFER,
    "%s/%s-%s.txt", LOGS_DIRECTORY, prefix, suffix);

  FILE * file = fopen(log_path, "a+");
  if (NULL == file) {
    critical("log_open(): Failed to open \"%s\" for appending.", log_path);
  }

  if (FALSE == log_setup) {
    fprintf(file, "# Log file: %s-%s.txt\n", prefix, suffix);

    log_setup = TRUE;
  }

  /* Return the file. */
  return file;
}