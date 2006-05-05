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

#ifndef MATRIX_BASIS_NON_BASIS_H
#define MATRIX_BASIS_NON_BASIS_H

#include "MatrixWithOpFactorized.hpp"

namespace AbstractLinAlgPack {

///
/** This is the interface to a basis and an nonbasis matrix.
  *
  * The form of this matrix is:
  *
  * M = [ C' ; N' ]
  *
  * It is also allowed that N = NULL and therefore:
  *
  * M = [ C' ]
  *
  * If this matrix is not initialized then rows() == cols() == 0.
  * If rows() > 0 then all of the MatrixOp functions can be
  * called with success.
  *
  * Note that it is not required that the nonsigular version of C be
  * available.  However if C_nonsingular is available, then it is assumed
  * of course that b == op(C) * inv(op(C_nonsingular)) * b so that they
  * are consistent.
  */
class MatrixBasisNonbasis
  : public MatrixOp
{
public:

    /// Base class for postcondition exceptions
  class MatrixNotAvailable : public std::runtime_error
  {public: MatrixNotAvailable(const std::string& what_arg) : std::runtime_error(what_arg) {}};

  /// Get reference to basis Matrix C ([m,m] = size(C))
  /**
    * Only call if rows() > 0 otherwise throws std::logic_error
    */
  virtual const MatrixOp& C() const = 0;

  /// Get reference to nonsingular version of the basis Matrix C ([m,m] = size(C))
  /**
    * Only call if rows() > 0 otherwise throws std::logic_error/
    * If the factorization is not available, then a MatrixNotAvailable exception
    * will be thrown.
    */
  virtual const MatrixFactorized& C_nonsingular() const = 0;

  /// Get reference to non-basis Matrix N ([n-m,n] = size(N))
  /**
    * Only call if cols() > rows() otherwise throws std::logic_error
    */
  virtual const MatrixOp& N() const = 0;

  // //////////////////////////////////////
  // Overridden from MatrixOp

  /// vs_lhs = alpha * op(M_rhs1) * vs_rhs2 + beta * vs_lhs (BLAS xGEMV)
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const;

  /// vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs (BLAS xGEMV)
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const;

};	// end class MatrixBasisNonbasis

}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_BASIS_NON_BASIS_H
