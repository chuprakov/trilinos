// /////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_extended_value_type.hpp
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

#ifndef EXTENDED_VALUE_TYPE_H
#define EXTENDED_VALUE_TYPE_H

namespace DenseLinAlgPack {

	typedef long double extended_value_type;
	// ToDo: Use doubledouble on platformes where long double
	// is no larger than double

} // end namespace DenseLinAlgPack

#endif // EXTENDED_VALUE_TYPE_H
