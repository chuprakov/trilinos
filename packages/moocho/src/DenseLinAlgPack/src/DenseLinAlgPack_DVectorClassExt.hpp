// ////////////////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_DVectorClassExt.hpp
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

#ifndef VECTOR_CLASS_EXT_H
#define VECTOR_CLASS_EXT_H

#include <valarray>

#include "DenseLinAlgPack_DVectorClassTmpl.hpp"

namespace DenseLinAlgPack {

	typedef VectorTmpl<extended_value_type>       VectorExt;
	typedef VectorSliceTmpl<extended_value_type>  VectorSliceExt;

} // end namespace DenseLinAlgPack

#endif	// end VECTOR_CLASS_EXT_H
