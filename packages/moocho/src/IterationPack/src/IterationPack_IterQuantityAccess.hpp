// /////////////////////////////////////////////////////////////////////////////////////
// IterQuantityAccess.h
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef ITER_QUANTITY_ACCESS_H
#define ITER_QUANTITY_ACCESS_H

#include "IterQuantity.h"

namespace GeneralIterationPack {

///
/** Interface to Iteration Quantities.
 *
 * Quantities are updated, read and queried given the
 * offset to the current iteration k.  For example, to set
 * a quantity for the <tt>k+1</tt> iteration you would call <tt>set_k(+1)</tt>.
 * The functions ending with <tt>prefix_k(offset)</tt> are meant to suggest
 * <tt>prefix k + offset</tt>.  For example:
 \verbatim

 has_storage_k(+1) => has storage k+1
 get_k(-1) => get k-1
 \endverbatim
 * Subclasses can implement this interface in a variety of ways.  But
 * they must follow a few simple rules: <ul>
 * <li>  Only forward transitions are allowed.  This effects the behavior of
 *	      \c IterQuantity::has_storage_k() and \c set_k().
 * </ul>
 *
 * The client should not have to worry about how much memory is
 * available.  Instead, it is for the object that configures the client
 * to provide the appropriate subclass to meet the needs of the client.
 */
template<class T_info>
class IterQuantityAccess : public IterQuantity {
public:

	///
	typedef	IterQuantity::NoStorageAvailable	NoStorageAvailable;
	///
	typedef IterQuantity::QuanityNotSet			QuanityNotSet;

	///
	/** Return a reference for the <tt>k + offset</tt> iteration quanity.
	 *
	 * Clients call this member function to access a quantity for a
	 * given iteration or modify the quantity which has already been
	 * set for that iteration.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->updated_k(offset) == true</tt> (throw QuanityNotSet)
	 * </ul>
	 */
	virtual T_info& get_k(int offset) = 0;

	///
	/** Return a const reference for the <tt>k + offset</tt> iteration quanity.
	  *
	  * Preconditions:<ul>
	  * <li> <tt>this->updated_k(offset) == true</tt> (throw QuanityNotSet)
	  * </ul>
	  *
	  * Clients call this member function to access a const quantity for a
	  * given iteration.
	  */
	virtual const T_info& get_k(int offset) const = 0;

	///
	/** Return a reference to the storage location for the <tt>k + offset</tt> iteration quanity.
	 *
	 * Precondtions:<ul>
	 * <li> <tt>this->has_storage_k(offset) == true</tt> (throw <tt>NoStorageAvailable</tt>)
	 * </ul>
	 *
	 * Postcondtions:<ul>
	 * <li> <tt>this->updated_k(offset) == true</tt>
	 * <li> <tt>this->updated_k(i) == false</tt> for <tt>i>/tt> in the set of
	 *      <tt>{ i : <tt>this->will_loose_mem(i,offset) == true</tt> }</tt> before this call.
	 * </ul> 
	 *
	 * This function will return a reference to the storage for the
	 * <tt>k + offset</tt> iteration quantity.  Calling this function
	 * may cause the loss of memory for a back iteration.  If
	 * <tt>will_loose_mem(back_offset,offset)</tt> returns
	 * <tt>true</tt> then <tt>updated_k(back_offset)</tt> will return
	 * <tt>false</tt> after this function returns (assuming an
	 * exception is not thrown).
	 *
	 * The client should expect nothing more than simple default
	 * initialization of the object who's reference is returned from
	 * this method.  The client is expected to use this reference to
	 * initalize this object appropriately.  If the client does not
	 * sufficiently update the object at <tt>k + offset</tt> before
	 * the reference is let go, the object's reference can be required
	 * with a call to the non-const version of <tt>get_k(offset)</tt>.
	 */
	virtual T_info& set_k(int offset) = 0;

	///
	/** Set the iteration quantity for the <tt>k + set_offset</tt>
	 * iteration to the <tt>k + get_offset</tt> iteration and return
	 * the reference to the <tt>k + set_offset</tt> iteration
	 * quantity.
	 *
	 * @param  set_offset  [in] The iteration offset to be set.
	 * @param  get_offset  [in[ The iteration offset to copy into the
	 *                     <tt>k + set_offset</tt> iteration.
	 *
	 * Precondtions:<ul>
	 * <li> <tt>this->has_storage_k(set_offset) == true</tt> (throw <tt>NoStorageAvailable</tt>)
	 * <li> <tt>this->updated_k(get_offset) == true</tt> (throw QuanityNotSet)
	 * </ul>
	 *
	 * Postcondtions:<ul>
	 * <li> <tt>this->updated_k(set_offset) == true</tt>
	 * <li> <tt>this->updated_k(i) == false</tt> for <tt>i>/tt> in the set of
	 *      <tt>{ i : <tt>this->will_loose_mem(i,offset) == true</tt> }</tt> before this call.
	 * </ul> 
	 *
	 * This method blends the functionality of the <tt>get_k()</tt>
	 * and the other <tt>set_k()</tt> methods.  This method ensures
	 * that a quantity from one iteration (<tt>k + get_offset</tt>)
	 * will be properly and efficienlty copied into the storage
	 * location of another iteration (<tt>k + set_offset</tt>) not
	 * matter how storage is handled by the implementing subclass.
	 */
	virtual T_info& set_k(int set_offset, int set_offset) = 0;

};	// end class IterQuantityAccess 

}	// end namespace GeneralIterationPack

#endif	// ITER_QUANTITY_ACCESS_H
