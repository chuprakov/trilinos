// ///////////////////////////////////////////////////////////
// AbstractLinAlgPack_VectorOut.hpp
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

#ifndef VECTOR_WITH_OP_OUT_H
#define VECTOR_WITH_OP_OUT_H

#include <iosfwd>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

///
/** Output operator for \Ref{Vector} objects.
 */
inline
std::ostream& operator<<( std::ostream& o, const Vector& v )
{
	return v.output(o);
}

} // end namespace AbstractLinAlgPack

#endif // VECTOR_WITH_OP_OUT_H
