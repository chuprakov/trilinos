// /////////////////////////////////////////////////////////////////////////////
// DecompositionSystemCoordinate.h
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

#ifndef DECOMPOSITION_SYSTEM_COORDINATE_H
#define DECOMPOSITION_SYSTEM_COORDINATE_H

#include "DecompositionSystemVarReductImpNode.h"
#include "IdentZeroVertConcatMatrixClass.h"
#include "SparseSolverPack/test/BasisSystemTesterSetOptions.h"
#include "Misc/include/StandardCompositionMacros.h"

namespace ConstrainedOptimizationPack {

///
/** Coordiante variable redunction interface class {abstract}.
  *
  * This is the interface for the coordinate variable reduction decomposition
  * where:
  *
  * #Y = [I ; 0]# (type \Ref{IdentZeroVertConcatMatrixSubclass}).
  *
  * The solution of the
  * linear systems #op(A(:,indep)'*Y) * x = y#, can be implemented here since,
  * #A(:,indep)' * Y = [C  N] * [I  0]' = C#.
  *
  * Also this fixes the definitions of of U and V as:
  *
  * #U = A(:,dep)' * Y = E#\\
  * #V = A(:,dep)' * Z = F - E * inv(C) * N#\\
  *
  * Here it is, therefore, assumed that all subclasses will use the passed in type
  * for U as E with the BasisSystem object and delete_access_matrices(...)
  * is defined here to refect that.  This behavior can of course be overrriden
  * by overridding this function.
  *
  * For now the copy constructor and the assignment operator are not defined.
  */
class DecompositionSystemCoordinate : public DecompositionSystemVarReductImpNode {
public:

#	///
	typedef SparseSolverPack::TestingPack::BasisSystemTester
		BasisSystemTester;

	/// «std comp» Members for the BasisSystem tester object
	STANDARD_COMPOSITION_MEMBERS( BasisSystemTester, basis_sys_tester )

	/** @name Public Types */
	//@{

	///
	typedef DecompositionSystemVarReductImpNode	inherited;

	///
	typedef inherited::InvalidMatrixType	InvalidMatrixType;

	//@}

	///
	/** Initializes #this# to use the BasisSystem<..> object given to it.
	  */
	DecompositionSystemCoordinate(BasisSystem* basis_sys, bool owns_basis_sys = false)
		: inherited(basis_sys, owns_basis_sys)
	{}

	/// Overridden
	void update_decomp(MatrixWithOp* A, MatrixWithOp* Z, MatrixWithOp* Y
		, MatrixWithOp* U, MatrixWithOp* V);

	/// Overridden
	void solve_transAtY(const VectorSlice& b, BLAS_Cpp::Transp trans, VectorSlice* x) const;

	/// Overridden
	void solve_transAtY(const GenMatrixSlice& B, BLAS_Cpp::Transp trans, GenMatrixSlice* X) const;

protected:

	// Template primative methods to be defined.

	/// Validate the type of Z.
	virtual void validate_Z(MatrixWithOp* Z) = 0;

	/// Update Z and V
	virtual void update_Z_and_V(MatrixWithOp* Z, GenMatrix* cV) = 0;

	///
	/** Overridden.
	  *
	  * Defined here since U = E and therefore we do not want to delete E since
	  * U is allocated from the outside.
	  */
	virtual void delete_access_matrices();

};	// end class DecompositionSystemCoordinate

}	// end namespace ConstrainedOptimizationPack

#endif	// DECOMPOSITION_SYSTEM_COORDINATE_H
