// ///////////////////////////////////////////////////////////////
// AbstractLinAlgPack_InnerProduct.hpp
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

#ifndef ALAP_INNER_PRODUCT_H
#define ALAP_INNER_PRODUCT_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

///
/** Abstract interface for inner products.
 *
 * ToDo: Finish documentaion
 */
class InnerProduct {
public:

	///
	/** Compute the inner product of two vectors.
	 *
	 * Preconditions:<ul>
	 * <li> ToDo: Spell out
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> ToDo: Spell out
	 * </ul>
	 *
	 * @param  v1  [in] First vector
	 * @param  v2  [in] Second vector
	 *
	 * @return  Returns some inner product of two vectors within a vector space..
	 */
	virtual value_type inner_prod(const Vector& v1, const Vector& v2) const = 0;

}; // end class InnerProduct

} // end namespace AbstractLinAlgPack

#endif  // ALAP_INNER_PRODUCT_H
