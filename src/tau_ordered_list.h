/*!
 * \file    tau_ordered_list.h
 * \ingroup estimating_volume_quotients
 *
 * \brief   The definition of data structures for collecting tau estimates in a
 *          truncated order list, and the declaration of functions for
 *          manipulating such lists.
 */

#ifndef TAU_ORDERED_LIST_H
#define TAU_ORDERED_LIST_H

#include <stdint.h>

/*!
 * \brief   A data structure for representing an order truncated list of tau
 *          estimates.
 */
typedef struct {
  /*!
   * \brief   The capacity of the ordered list of tau estimates.
   */
  uint32_t capacity;

  /*!
   * \brief   The number of tau estimates in the list.
   */
  uint32_t size;

  /*!
   * \brief   The list of tau values.
   */
  long double * tau;
} Tau_Ordered_List;


/*!
 * \name Initialization
 * \{
 */

/*!
 * \brief   Initializes an ordered list of tau estimates with a given capacity.
 *
 * This function allocates memory.
 *
 * \param[in, out] list   The ordered list to initialize.
 * \param[in] capacity    The capacity of the list.
 */
void tau_ordered_list_init(
  Tau_Ordered_List * const list,
  const unsigned int capacity);

/*!
 * \brief   Clears an ordered list of tau estimates.
 *
 * This function deallocates the memory allocated by tau_ordered_list_init().
 *
 * \param[in, out] list   The ordered list to clear.
 */
void tau_ordered_list_clear(
  Tau_Ordered_List * const list);

/*!
 * \}
 */

/*!
 * \name Inserting
 * \{
 */

/*!
 * \brief   Inserts a tau estimates into an ordered list.
 *
 * The estimate will be merged into the list to ensure that it is still sorted.
 *
 * If the list is at full capacity, one tau estimate will be excluded from the
 * list. Note that the estimate inserted may in this case be excluded if it is
 * smaller than all other estimates in the list.
 *
 * \param[in, out] list   The ordered list into which to insert the estimate.
 * \param[in] tau         The tau estimate to insert into the ordered list.
 */
void tau_ordered_list_insert(
  Tau_Ordered_List * const list,
  const long double tau);

/*!
 * \}
 */

/*!
 * \name Merging via message passing
 * \{
 */

 /*!
  * \brief    Sends an ordered list to another node.
  *
  * \param[in, out] list    The ordered list to send.
  * \param[in] rank         The rank of the node to which to send the list.
  */
void tau_ordered_list_send_merge(
  const Tau_Ordered_List * const list,
  const int rank);

/*!
  * \brief    Receives an ordered list to another node and merges it into an
  *           ordered list.
  *
  * \param[in, out] list    The ordered list into which the merge the ordered
  *                         list received from the node.
  * \param[in] rank         The rank of the node from which to receive the list.
  */
void tau_ordered_list_recv_merge(
  Tau_Ordered_List * const list,
  const int rank);

/*!
 * \}
 */

/*!
 * \name Printing
 * \{
 */

/*!
 * \brief   Prints an ordered list of tau estimates.
 *
 * \param[in] list   The ordered list to print.
 */
void tau_ordered_list_print(
  const Tau_Ordered_List * const list);

/*!
 * \}
 */

#endif /* TAU_ORDERED_LIST_H */
