/*!
 * \file    executables.h
 * \ingroup executable
 *
 * \brief   The declaration of constants shared by all executables.
 */

/*!
 * \defgroup executable Executables
 *
 * \brief    A module for executables.
 */

#ifndef EXECUTABLES_H
#define EXECUTABLES_H

#include <sys/stat.h>

/*!
 * \brief   The relative path to the distribution directory.
 */
const static char * const DISTRIBUTIONS_DIRECTORY = "distributions";

/*!
 * \brief   The relative path to the log directory.
 */
const static char * const LOGS_DIRECTORY = "logs";

/*!
 * \brief   The relative path to the plots directory.
 */
const static char * const PLOTS_DIRECTORY = "plots";

/*!
 * \brief   The maximum path buffer size.
 *
 * The paths are setup in such a matter that this path size cannot be exceeded.
 */
#define MAX_SIZE_PATH_BUFFER 1024

/*!
 * \brief   The default permission to use when creating directories.
 */
#define DIRECTORY_PERMISSIONS \
  (S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)

#endif /* EXECUTABLES_H */
