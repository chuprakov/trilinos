// ////////////////////////////////////////////////////////////////////////////////////////
// VectorClass.h

#ifndef VECTOR_CLASS_H
#define VECTOR_CLASS_H

#include <valarray>

#include "VectorClassTmpl.h"

namespace LinAlgPack {

	typedef VectorTmpl<value_type>       Vector;
	typedef VectorSliceTmpl<value_type>  VectorSlice;

} // end namespace LinAlgPack

#endif	// end VECTOR_CLASS_H
