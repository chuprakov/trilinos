// ////////////////////////////////////////////////////////////////////
// SparseElementOut.h

#ifndef SPARSE_ELEMENT_OUT_H
#define SPARSE_ELEMENT_OUT_H

#include <ostream>

#include "SparseElement.h"

namespace SparseLinAlgPack {

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

} // end namespace SparseLinAlgPack

#endif // SPARSE_ELEMENT_OUT_H