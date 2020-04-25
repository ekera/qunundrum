/*!
 * \file    errors.h
 * \ingroup errors
 *
 * \brief The declaration of functions for error handling.
 */

/*!
 * \defgroup  errors Error handling
 * \ingroup   utility
 *
 * \brief A module for error handling functions.
 */

#ifndef ERRORS_H
#define ERRORS_H

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief   Reports an error to standard error and terminates the executable
 *          with an error code without returning control to the caller.
 *
 * For formatting specifiers, see the man page for printf().
 *
 * \param[in] msg   The error message, with optional formatting specifiers and
 *                  corresponding trailing arguments.
 */
void critical(const char * msg, ...);

#ifdef __cplusplus
}
#endif

#endif /* ERRORS_H */
