// ////////////////////////////////////////////////////////////////////////////////////////
// VectorClassExt.h

#ifndef VECTOR_CLASS_EXT_H
#define VECTOR_CLASS_EXT_H

#include <valarray>

#include "VectorClassTmpl.h"

namespace LinAlgPack {

	typedef VectorTmpl<extended_value_type>       VectorExt;
	typedef VectorSliceTmpl<extended_value_type>  VectorSliceExt;

} // end namespace LinAlgPack

#endif	// end VECTOR_CLASS_EXT_H
