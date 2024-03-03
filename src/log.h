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

#endif /* LOG_H */