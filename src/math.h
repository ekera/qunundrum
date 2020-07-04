/*!
 * \file    math.h
 * \ingroup math
 *
 * \brief   The declaration of basic mathematical functions.
 */

/*!
 * \defgroup math Mathematics
 *
 * \brief A module for mathematical functions.
 */

#ifndef MATH_H
#define MATH_H

#include <gmp.h>

#include <stdint.h>

/*!
 * \name Absolute values
 * \{
 */

/*!
 * \brief Returns the absolute value of a signed 32 bit integer.
 *
 * \param[in] x   The value x.
 *
 * \return The absolute value of x.
 */
uint32_t abs_i(const int32_t x);

/*!
 * \brief Returns the absolute value of a double.
 *
 * \param[in] x   The value x.
 *
 * \return The absolute value of x.
 */
double abs_d(const double x);

/*!
 * \brief Returns the absolute value of a long double.
 *
 * \param[in] x   The value x.
 *
 * \return The absolute value of x.
 */
long double abs_ld(const long double x);

/*!
 * \}
 */

/*!
 * \name Signs
 * \{
 */

/*!
 * \brief Returns the sign value of a signed 32 bit integer.
 *
 * \param[in] x   The value x.
 *
 * \return The sign of x.
 */
int32_t sgn_i(const int32_t x);

/*!
 * \brief Returns the sign value of a double.
 *
 * \param[in] x   The value x.
 *
 * \return The sign of x.
 */
int32_t sgn_d(const double x);

/*!
 * \brief Returns the sign value of a long double.
 *
 * \param[in] x   The value x.
 *
 * \return The sign of x.
 */
int32_t sgn_ld(const long double x);

/*!
 * \}
 */

/*!
 * \name Maximum
 * \{
 */

/*!
 * \brief Returns the maximum of two signed 32 bit integers.
 *
 * \param[in] a   The value a.
 * \param[in] b   The value b.
 *
 * \return The maximum of the values a and b.
 */
int32_t max_i(const int32_t a, const int32_t b);

/*!
 * \brief Returns the maximum of two unsigned 32 bit integers.
 *
 * \param[in] a   The value a.
 * \param[in] b   The value b.
 *
 * \return The maximum of the values a and b.
 */
uint32_t max_ui(const uint32_t a, const uint32_t b);

/*!
 * \}
 */

/*!
 * \name Problem-specific
 * \{
 */

/*!
 * \brief Reduces x modulo n and constrains the result to [n/2, n/2).
 *
 * In the notation of Shor and following authors, this function computes {x}_n.
 *
 * \param[in, out] x   The integer x.
 * \param[in]      n   The integer n.
 */
void mod_reduce(mpz_t x, const mpz_t n);

/*!
 * \brief   Computes the greatest integer kappa such that 2^kappa divides x.
 *
 * \param[in] x   The integer x.
 *
 * \return The greatest integer kappa such that 2^kappa divides x.
 */
uint32_t kappa(const mpz_t x);

/*!
 * \brief   Checks if an unsigned 32 bit integer is a power of two.
 *
 * \param[in] x   The integer x.
 *
 * \return Returns #TRUE if the integer is a power of two, #FALSE otherwise.
 */
bool is_pow2(const uint32_t x);

/*!
 * \}
 */

#endif /* MATH_H */
