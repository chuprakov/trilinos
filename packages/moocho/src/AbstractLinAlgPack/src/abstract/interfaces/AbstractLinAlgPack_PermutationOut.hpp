// ///////////////////////////////////////////////////////////
// AbstractLinAlgPack_PermutationOut.hpp
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

#ifndef PERMUTATION_OUT_H
#define PERMUTATION_OUT_H

#include <iosfwd>

#include "AbstractLinAlgPack_Permutation.hpp"

namespace AbstractLinAlgPack {

///
/** Output operator for Permutation objects.
 */
inline
std::ostream& operator<<( std::ostream& o, const Permutation& P )
{
	return P.output(o);
}

} // end namespace AbstractLinAlgPack

#endif // PERMUTATION_OUT_H
