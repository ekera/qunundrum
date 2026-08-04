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
 * The log file is created in the #LOGS_DIRECTORY directory with a file name of
 * the form "prefix-YYYYMMDD-HHmmss±ZZZZ-XXXXXXXX.txt" where YYYYMMDD-HHmmss is
 * the local date and time, ±ZZZZ is the local time zone's signed offset from
 * UTC in hours and minutes (HHmm), and XXXXXXXX is a 32-bit unsigned integer
 * read from /dev/urandom and printed in hexadecimal form with leading zeros.
 * The prefix is provided by the caller to this function.
 *
 * The rationale for appending the suffix to the prefix is to prevent
 * concurrently running executables from writing to the same log file. The
 * suffix is cached so that multiple calls to this function within the same
 * executable all result in the same log file being opened.
 *
 * When the log file is first created by this function, a header that contains
 * the name of the log file is written to the file. No header is written to the
 * file if the file is re-opened in subsequent calls to this function.
 *
 * \param[in] prefix   The log file name prefix.
 *
 * \return  The log file.
 */
FILE * log_open(
  const char * const prefix);

#endif /* LOG_H */