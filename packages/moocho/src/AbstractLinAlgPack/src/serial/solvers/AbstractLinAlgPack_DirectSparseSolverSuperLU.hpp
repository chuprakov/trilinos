// /////////////////////////////////////////////////////////////
// DirectSparseSolverSuperLU.h
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

#ifdef SPARSE_SOLVER_PACK_USE_SUPERLU

#ifndef	DIRECT_SPARSE_SOLVER_SUPERLU_H
#define DIRECT_SPARSE_SOLVER_SUPERLU_H

#include <valarray>
#include <vector>
#include <string>

#include "DirectSparseSolverImp.h"
#include "SuperLUSolver.h"
#include "LinAlgPack/src/VectorClass.h"
#include "LinAlgPack/src/IVector.h"
#include "StandardMemberCompositionMacros.h"

namespace SparseSolverPack {

///
/** Concreate sparse solver subclass that uses SuperLU.
 *
 * ToDo: Finish documentation!
 */
class DirectSparseSolverSuperLU : public DirectSparseSolverImp {
public:

	/** @name Control parameters */
	//@{

	// ToDo: Fill these in!

	//@}

	/** @name Constructors/initializers */
	//@{

	///
	/** Default constructor */
	DirectSparseSolverSuperLU();

	//@}

	/** @name Overridden from DirectSparseSolver */
	//@{

	///
	const basis_matrix_factory_ptr_t basis_matrix_factory() const;
	///
	void estimated_fillin_ratio( value_type estimated_fillin_ratio );

	//@}

protected:

	/** @name Protected types */
	//@{

	///
	/** Implements the BasisMatrix object for SuperLU.
	 */
	class BasisMatrixSuperLU : public BasisMatrixImp {
	public:

		/** @name Overridden from BasisMatrixImp */
		//@{

		///
		MemMngPack::ref_count_ptr<BasisMatrixImp> create_matrix() const;
		///
		void V_InvMtV(
			VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
			,const VectorWithOp& v_rhs2) const ;
		
		//@}

	}; // end class BasisMatrixSuperLU

	///
	/** Stores the factorization structure for SuperLU
	 */
	class FactorizationStructureSuperLU : public FactorizationStructure {
	public:
		friend class DirectSparseSolverSuperLU;
		friend class BasisMatrixSuperLU;
	private:
		MemMngPack::ref_count_ptr<SuperLUPack::SuperLUSolver>
			superlu_solver_;
		MemMngPack::ref_count_ptr<SuperLUPack::SuperLUSolver::FactorizationStructure>
			fact_struct_;
		FactorizationStructureSuperLU();
	}; // end class FactorizationStructureSuperLU

	///
	/** Stores the factorization nonzeros for SuperLU
	 */
	class FactorizationNonzerosSuperLU : public FactorizationNonzeros {
	public:
		friend class DirectSparseSolverSuperLU;
		friend class BasisMatrixSuperLU;
	private:
		MemMngPack::ref_count_ptr<SuperLUPack::SuperLUSolver::FactorizationNonzeros>
			fact_nonzeros_;
	}; // end class FactorizationNonzerosSuperLU

	//@}

	/** @name Overridden from DirectSparseSolverImp */
	//@{

	///
	const MemMngPack::ref_count_ptr<FactorizationStructure> create_fact_struc() const;
	///
	const MemMngPack::ref_count_ptr<FactorizationNonzeros> create_fact_nonzeros() const;
	///
	void imp_analyze_and_factor(
		const SparseLinAlgPack::MatrixConvertToSparse   &A
		,FactorizationStructure                         *fact_struc
		,FactorizationNonzeros                          *fact_nonzeros
		,LinAlgPack::IVector                            *row_perm
		,LinAlgPack::IVector                            *col_perm
		,size_type                                      *rank
		,std::ostream                                   *out
		);
	///
	void imp_factor(
		const SparseLinAlgPack::MatrixConvertToSparse   &A
		,const FactorizationStructure                   &fact_struc
		,FactorizationNonzeros                          *fact_nonzeros
		,std::ostream                                   *out
		);

	//@}

private:

	// /////////////////////////////////
	// Private data members

	// ////////////////////////////////
	// Private member functions

};	// end class DirectSparseSolverSuperLU 

// ////////////////////////////////////////
// Inline members

}	// end namespace SparseSolverPack 

#endif	// DIRECT_SPARSE_SOLVER_SUPERLU_H

#endif // SPARSE_SOLVER_PACK_USE_SUPERLU
