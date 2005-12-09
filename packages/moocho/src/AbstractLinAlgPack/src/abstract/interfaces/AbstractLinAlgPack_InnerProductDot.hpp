// ///////////////////////////////////////////////////////////////
// AbstractLinAlgPack_InnerProductDot.hpp
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

#ifndef ALAP_INNER_PRODUCT_DOT_H
#define ALAP_INNER_PRODUCT_DOT_H

#include "AbstractLinAlgPack_InnerProduct.hpp"

namespace AbstractLinAlgPack {

///
/** Implements the inner product as the dot product.
 *
 * ToDo: Finish documentaion
 */
class InnerProductDot : public InnerProduct {
public:

	/** @name Overridden from InnerProduct */
	//@{
	///
	value_type inner_prod(const Vector& v1, const Vector& v2) const;
	//@}

}; // end class InnerProductDot

} // end namespace AbstractLinAlgPack

#endif  // ALAP_INNER_PRODUCT_DOT_H
