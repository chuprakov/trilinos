// ///////////////////////////////////////////////////////////////
// AbstractLinAlgPackTypes.h
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

#ifndef ABSTRACT_LIN_ALG_PACK_TYPES_H
#define ABSTRACT_LIN_ALG_PACK_TYPES_H

#include <memory>

#include "RTOpPack/include/RTOp.h"
#include "BLAS_CppTypes.h"
#include "Range1D.h"

namespace AbstractLinAlgPack {

typedef RTOp_index_type  size_type;
typedef RTOp_value_type  value_type;
typedef RTOp_index_type  index_type;

typedef RangePack::Range1D Range1D; // For some reason doxygen likes typedef more than using?

/** @name Main interface library */
//@{

/// Enumeration for returning the amount of overlap between two objects
enum EOverLap { NO_OVERLAP = 0, SOME_OVERLAP, SAME_MEM };	

// pure abstract classes

class VectorSpaceBase;
class VectorSpace;

class VectorBase;
class VectorBaseMutable;
class VectorWithOp;
class VectorWithOpMutable;

class MatrixSpaceBase;
template<class M>
class MatrixSpace;

class MatrixBase;
class MatrixWithOp;
class MatrixNonsingular;
class MatrixWithOpNonsingular;
class MatrixSymWithOp;
class MatrixSymNonsingular;
class MatrixSymWithOpNonsingular;

class MultiVector;
class MultiVectorMutable;

class BasisSystem;

// template classes

template <class T_Indice, class T_Value>	class SparseElement;
template <class T_Element, class T_Alloc>	class SparseVector;
template <class T_Element>					class SparseVectorSlice;

// concrete classes

class EtaVector;
class GenPermMatrixSlice;
typedef SparseVector<
	SparseElement<index_type,value_type>
	, std::allocator<
		SparseElement<index_type,value_type>
	  >
  >												SpVector;
typedef SparseVectorSlice<
	SparseElement<index_type,value_type> >		SpVectorSlice;

//@}

/** @name Standard tools library */
//@{

// pure abstract classes

class PermVector;

class MatrixSymInitDiagonal;
class MatrixSymDiagonal;

// template classes

template <class M_itfc, class M_impl, class T_PostMod>  class MatrixSpaceStd;

// concrete subclasses

class VectorSpaceCompositeStd;
class VectorWithOpMutableCompositeStd;
class MatrixCompositeStd;
class MatrixSymIdentity;
class MatrixSymDiagonalStd;
class MatrixZero;

// testing classes

class VectorSpaceTester;
class VectorSpaceTesterSetOptions;
class BasisSystemTester;
class BasisSystemTesterSetOptions;

//@}

} // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_TYPES_H
