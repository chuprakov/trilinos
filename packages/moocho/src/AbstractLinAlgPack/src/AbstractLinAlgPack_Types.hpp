// ///////////////////////////////////////////////////////////////
// AbstractLinAlgPackTypes.h

#ifndef ABSTRACT_LIN_ALG_PACK_TYPES_H
#define ABSTRACT_LIN_ALG_PACK_TYPES_H

#include <memory>

#include "RTOpPack/include/RTOp.h"
#include "BLAS_CppTypes.h"

namespace RangePack {
	class Range1D;
}

namespace AbstractLinAlgPack {

typedef RTOp_index_type  size_type;
typedef RTOp_value_type  value_type;
typedef RTOp_index_type  index_type;

using RangePack::Range1D;

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
class MatrixWithOpMutable;
class MatrixFactorized;
class MatrixWithOpFactorized;
class MatrixSymWithOp;
class MatrixSymWithOpMutable;
class MatrixSymFactorized;
class MatrixSymWithOpFactorized;

// template classes

template <class T_Indice, class T_Value>	class SparseElement;
template <class T_Element, class T_Alloc>	class SparseVector;
template <class T_Element>					class SparseVectorSlice;

// concrete classes

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

template <class M_itfc, class M_impl>   class MatrixSpaceStd;

// concrete subclasses

class VectorSpaceCompositeStd;
class VectorWithOpMutableCompositeStd;
class MatrixSymIdentity;
class MatrixSymDiagonalStd;

//@}

} // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_TYPES_H
