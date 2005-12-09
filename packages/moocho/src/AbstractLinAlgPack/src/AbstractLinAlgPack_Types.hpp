// ///////////////////////////////////////////////////////////////
// AbstractLinAlgPack_Types.hpp
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
#include <stdexcept>

#include "RTOp.h"
#include "BLAS_Cpp_Types.hpp"
#include "DenseLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

#include "DenseLinAlgPack_PublicTypes.ud"

typedef RTOp_index_type  size_type;
typedef RTOp_value_type  value_type;
typedef RTOp_index_type  index_type;

#ifdef DOXYGEN_COMPILE // Doxygen needs a little help finding these links
///
typedef DenseLinAlgPack::VectorTmpl<value_type>             DVector;
///
typedef DenseLinAlgPack::VectorSliceTmpl<value_type>        DVectorSlice;
///
typedef DenseLinAlgPack::DMatrix                            DMatrix;
///
typedef DenseLinAlgPack::DMatrixSlice                       DMatrixSlice;
///
typedef DenseLinAlgPack::DMatrixSliceTriEle                 DMatrixSliceTriEle;
///
typedef DenseLinAlgPack::DMatrixSliceTri                    DMatrixSliceTri;
///
typedef DenseLinAlgPack::DMatrixSliceSym                    DMatrixSliceSym;
///
typedef RangePack::Range1D Range1D;
#endif


/** @name Exception classes */
//@{

/// Base class for precondition exceptions
class PreConditionException : public std::logic_error
{public: PreConditionException(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Base class for postcondition exceptions
class PostConditionException : public std::runtime_error
{public: PostConditionException(const std::string& what_arg) : std::runtime_error(what_arg) {}};

/// Base class for input exceptions (Preconditions).
class InputException : public PreConditionException
{public: InputException(const std::string& what_arg) : PreConditionException(what_arg) {}};

/// Base class for invalid setup for a class object when an exception is thrown
class SetupException : public PreConditionException
{public: SetupException(const std::string& what_arg) : PreConditionException(what_arg) {}};


//@}

/** @name Main interface library */
//@{

class InnerProduct;

class VectorSpaceFactory;
class VectorSpace;
class Vector;
class VectorMutable;

class MatrixBase;
class MatrixOp;
class MatrixNonsing;
class MatrixOpNonsing;
class MatrixSymOp;
class MatrixSymNonsing;
class MatrixSymOpNonsing;
class MatrixSymDiag;

class MultiVector;
class MultiVectorMutable;

class MatrixSymSecant;

class BasisSystem;
class BasisSystemPerm;
class BasisSystemFactory;

class Permutation;

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

class MatrixSymInitDiag;
class MatrixSymDiag;

// concrete subclasses

class BasisSystemComposite;
class VectorSpaceBlocked;
class VectorMutableBlocked;
class MatrixOpSubView;
class MatrixComposite;
class MatrixSymIdent;
class MatrixSymDiagStd;
class MatrixZero;
class MatrixPermAggr;
class MatrixOpNonsingAggr;

// testing classes

class VectorSpaceTester;
class VectorSpaceTesterSetOptions;
class MatrixOpNonsingTester;
class BasisSystemTester;
class BasisSystemTesterSetOptions;

//@}

/** @name Serial interface library */
//@{

// pure abstract classes

class MatrixOpSerial;
class MatrixNonsingSerial;
class MatrixSymOpSerial;
class MatrixSymNonsingSerial;
class MatrixOpNonsingSerial;
class MatrixSymOpNonsingSerial;
class MatrixSymDenseInitialize;
class MatrixSymDiagSparse;
class MatrixLoadSparseElements;
class MatrixConvertToSparse;
class MatrixExtractSparseElements;
class MatrixExtractInvCholFactor;
class MatrixSymOpGetGMSSymMutable;
class MatrixSymOpGetGMSSym;
class MatrixSymAddDelUpdateable;

//@}

/** @name Serial implementations library */
//@{

class VectorDenseEncap;
class VectorDenseMutableEncap;
class MatrixDenseEncap;
class MatrixDenseMutableEncap;
class MatrixDenseSymEncap;
class MatrixDenseSymMutableEncap;
class MatrixDenseTriEncap;

class PermutationSerial;
class VectorSpaceSerial;
class VectorMutableDense;
class VectorSparse;
class MatrixSparseCOORSerial;
class MatrixSymPosDefCholFactor;
class MatrixConvertToSparseEncap;
class MultiVectorMutableDense;

class MatrixSymDiagSparseStd;

//@}

/** @name Serial solvers library */
//@{

// Matrix scaling classes

class MatrixScaling_Strategy;

// Sparse linear solver classes

class DirectSparseSolver;        // Abstract interface
class DirectSparseSolverImp;     // Node implementation classs
class DirectSparseSolverMA28;    // Concrete subclass
class DirectSparseSolverMA48;    // ""
class DirectSparseSolverSuperLU; // ""

// BasisSystem classes

class BasisSystemPermDirectSparse;
class BasisSystemFactoryStd;

//@}

} // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_TYPES_H
