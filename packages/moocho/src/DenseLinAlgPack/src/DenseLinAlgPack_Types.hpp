// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef LINALGPACK_TYPES_H
#define LINALGPACK_TYPES_H

#include "DenseLinAlgPack_Options.hpp"
#include "Thyra_Range1D.hpp"
#include "BLAS_Cpp_Types.hpp"

namespace DenseLinAlgPack {

/* * @name {\bf DenseLinAlgPack Type Declarations}.
  *
  * These are forward declarations of the types used with the DenseLinAlgPack
  * package (namespace).  In addition the BLAS_Cpp enumerations
  * \Ref{Transp}, \Ref{Side}, \Ref{Uplo}, and \Ref{Diag} and there values
  * are avalible using the qualifier #BLAS_Cpp#.
  */

// @{

///
using RangePack::Range1D;
#ifdef _INTEL_CXX
using RangePack::full_range;
#endif
///
using BLAS_Cpp::rows;
///
using BLAS_Cpp::cols;
///
using BLAS_Cpp::trans_not;

/// Enumeration for returning the amount of overlap between two objects
enum EOverLap { NO_OVERLAP = 0, SOME_OVERLAP, SAME_MEM };	

///
class IVector;
///
class Range1D;
///
template<class T>
class VectorTmpl;
///
template<class T>
class VectorSliceTmpl;
///
typedef VectorTmpl<value_type>                DVector;
///
typedef VectorSliceTmpl<value_type>           DVectorSlice;
///
typedef VectorTmpl<extended_value_type>       VectorExt;
///
typedef VectorSliceTmpl<extended_value_type>  VectorSliceExt;
///
class TransVectorSlice;
///
class DMatrix;
///
class DMatrixSlice;
///
class TransGenMatrixSlice;
///
class DMatrixSliceTriEle;
///
class DMatrixSliceTri;
///
class DMatrixSliceSym;

// @}

}  // namespace DenseLinAlgPack

#endif // LINALGPACK_TYPES_H
