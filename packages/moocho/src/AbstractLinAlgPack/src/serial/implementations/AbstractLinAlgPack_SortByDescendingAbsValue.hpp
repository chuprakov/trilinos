// ////////////////////////////////////////////////////////////
// AbstractLinAlgPack_SortByDescendingAbsValue.hpp
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

#ifndef SORT_BY_DESCENDING_ABS_VALUE_H
#define SORT_BY_DESCENDING_ABS_VALUE_H

#include <math.h>

#include "AbstractLinAlgPack_SpVectorClass.hpp"

namespace AbstractLinAlgPack {

///
/** Function object class for sorting a sparse vectors in descending order
  * by abs(v(i)).
  */
class SortByDescendingAbsValue {
public:
	bool operator()( const AbstractLinAlgPack::SpVector::element_type& x
		, const AbstractLinAlgPack::SpVector::element_type& y ) const
	{
		return ::fabs(x.value()) > ::fabs(y.value());
	}
};	// end class AbsMultVal


}	// end namespace AbstractLinAlgPack 

#endif 	// SORT_BY_DESCENDING_ABS_VALUE_H
