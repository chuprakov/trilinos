// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// ///////////////////////////////////////////////////////////////
// TSFCoreTypes.hpp

#ifndef TSF_CORE_TYPES_HPP
#define TSF_CORE_TYPES_HPP

#include "TSFCore_ConfigDefs.hpp"

#include "Range1D.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "RTOpPackTypes.hpp"

namespace TSFCore {

/** \defgroup TSFCoreBasicTypes_grp Basic TSFCore types.
 * \ingroup TSFCore_fundamental_interfaces_code_grp
 */
//@{

//
// Basic types
//

/// Type for the dimension of a vector space
typedef RTOp_index_type  Index;

/// Type for a range of indices
typedef RangePack::Range1D   Range1D;

///
/** Enumeration for determining how an operator is applied.
 */
enum ETransp {
  NOTRANS     ///< Use the transposed operator
  ,TRANS      ///< Use the nontransposed operator
  ,CONJTRANS  ///< Use transpose-conjugate if complex and otherwise just means <tt>TRANS</tt>
};

///
/** Return a string name for a <tt>ETransp</tt> value.
 */
inline
const char* toString(ETransp transp) { return ( transp == NOTRANS ? "NOTRANS" : ( transp == TRANS ? "TRANS" : "CONJTRANS" ) ); }

///
/** Perform a not operation on an ETransp value
 */
inline
ETransp not_trans( ETransp trans )
{
	return ( trans == NOTRANS ? TRANS : NOTRANS );
}

///
/** Combine two transpose arguments
 */
inline
ETransp trans_trans( ETransp trans1, ETransp trans2 )
{
	return ( trans1 == trans2 ? NOTRANS : ( trans1==CONJTRANS || trans2==CONJTRANS ? CONJTRANS : TRANS ) );
}

//@}

namespace Exceptions {

/** \defgroup TSFCoreExceptions_grp Basic TSFCore exception types.
 * \ingroup TSFCore_fundamental_interfaces_code_grp
 */
//@{

/// Thrown if any member functions are called before initialize() has been called.
class UnInitialized : public std::logic_error
{public: UnInitialized(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if vector spaces are incompatible
class IncompatibleVectorSpaces : public std::logic_error
{public:
	IncompatibleVectorSpaces(const std::string& what_arg) : std::logic_error(what_arg) {}
//	IncompatibleVectorSpaces(const IncompatibleVectorSpaces& ivs) : std::logic_error(ivs.what()) {}
};

/// Thrown if the argument <tt>M_trans</tt> is not supported,
class OpNotSupported : public std::logic_error
{public: OpNotSupported(const std::string& what_arg) : std::logic_error(what_arg) {}};

//@}

} // namespace Exceptions

//
// TSFCore
//

// Core abstract interface classes

template<class Scalar> class VectorSpaceFactory;
template<class Scalar> class VectorSpace;
template<class Scalar> class Vector;
template<class Scalar> class LinearOp;
template<class Scalar> class MultiVector;

// Misc interface classes

template<class Scalar> class LinearOpHandle;

// Basic node support subclasses and interfaces

template<class Scalar> class ScalarProd;
template<class Scalar> class ScalarProdVectorSpaceBase;
template<class Scalar> class EuclideanLinearOpBase;
template<class Scalar> class SerialVectorSpaceBase;
template<class Scalar> class SerialVectorBase;

// Basic concrete support subclasses

template<class Scalar> class EuclideanScalarProd;
template<class Scalar> class LinearOpScalarProd;
template<class Scalar> class SerialVectorSpaceFactoryStd;
template<class Scalar> class SerialVectorSpaceStd;
template<class Scalar> class SerialVectorStd;
template<class Scalar> class SerialMultiVectorStd;
template<class Scalar> class MultiVectorCols;
template<class Scalar> class VectorMultiVector;

} // end namespace TSFCore

#endif // TSF_CORE_TYPES_HPP
