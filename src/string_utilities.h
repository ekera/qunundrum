/*!
 * \file    string_utilities.h
 * \ingroup string_utility
 *
 * \brief   The declaration of utility functions for manipulating strings.
 */

/*!
 * \defgroup string_utility String utilities
 * \ingroup  utility
 *
 * \brief    A group for string utility functions.
 */

#ifndef STRING_UTILITIES_H
#define STRING_UTILITIES_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief   Calls snprintf(), terminating the executable with an error message
 *          by calling critical() if the buffer capacity is exceeded.
 *
 * \param[in, out] str    The buffer in which to store the string.
 * \param[in] size        The capacity of the buffer in characters.
 * \param[in] format      The formatting string.
 * \param[in] ...         Arguments for the formatting specifiers.
 */
void safe_snprintf(char * str, size_t size, const char * format, ...);

/*!
 * \brief   Copies a string from a source buffer to a destination buffer,
 *          terminating the executable with an error message by calling
 *          critical() if the destination buffer capacity is exceeded.
 *
 * This function will NULL-terminate the destination string. The source string
 * is assumed to be NULL-terminated.
 *
 * \param[in, out] dst    The destination buffer in which to store the string.
 * \param[in] src         The source buffer from which to read.
 * \param[in] dst_size    The capacity of the destination buffer in characters.
 */
void safe_strlcpy(char * dst, const char * src, size_t dst_size);

/*!
 * \brief   Truncates a path by returning only its last component.
 *
 * \param[in] path        The path.
 *
 * \return The last component of the path.
 */
const char * truncate_path(const char * path);

#ifdef __cplusplus
}
#endif

#endif /* STRING_UTILITIES_H */
