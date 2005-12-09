// ////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_SparseElementOut.hpp
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

#ifndef SPARSE_ELEMENT_OUT_H
#define SPARSE_ELEMENT_OUT_H

#include <ostream>

#include "AbstractLinAlgPack_SparseElement.hpp"

namespace AbstractLinAlgPack {

class SparsePtrVector;

///
/** Outputs a SparseElement<> object.
  *
  * Output format is:
  *
  * #ele.value():ele.indice()#\\
  */
template <class T_Indice, class T_Value>
std::ostream& operator<<(std::ostream& o, const SparseElement<T_Indice,T_Value>& ele);

// /////////////////////////////////////////
// Inline definitions
template <class T_Indice, class T_Value>
inline std::ostream& operator<<(std::ostream& o, const SparseElement<T_Indice,T_Value>& ele) {
	return o << ele.value() << ":" << ele.indice();
}

} // end namespace AbstractLinAlgPack

#endif // SPARSE_ELEMENT_OUT_H
