/*!
 * \file    log.h
 * \ingroup executable
 *
 * \brief   The declaration of convenience functions for logging in executables.
 */

#ifndef LOG_H
#define LOG_H

#include <stdio.h>

/*!
 * \brief   Logs a timestamp.
 *
 * \param[in, out] file   The file to which to log the timestamp.
 */
void log_timestamp_fprintf(
  FILE * const file);

/*!
 * \brief   Opens a log file unique to this process for appending, creating
 *          #LOGS_DIRECTORY first if it does not already exist.
 *
 * The log file is created in #LOGS_DIRECTORY with a path of the form
 * "prefix-YYYYMMDD-HHmmss±ZZZZ-XXXXXXXX.txt" where YYYYMMDD-HHmmss is the local
 * date and time, ±ZZZZ is the local time zone's signed offset from UTC in hours
 * and minutes (HHmm), and XXXXXXXX is a 32-bit unsigned integer read from
 * /dev/urandom and printed in hexadecimal with leading zeros. The prefix is
 * provided when calling this function.
 *
 * The date, offset and random integer are generated once and cached for the
 * remainder of the process, so repeated calls to this function — e.g. to log
 * further entries later on during the same run — reopen the same log file. This
 * prevents concurrently running executables from writing to the same log file.
 *
 * \param[in] prefix   The log file name prefix.
 *
 * \return  The log file.
 */
FILE * log_open(
  const char * const prefix);

#endif /* LOG_H */