// /////////////////////////////////////////////////////////////
// DirectSparseSolverMA28.h
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

#ifndef	DIRECT_SPARSE_FORTRAN_COMPATIBLE_SOLVER_MA28_H
#define DIRECT_SPARSE_FORTRAN_COMPATIBLE_SOLVER_MA28_H

#include "DirectSparseSolverImp.h"

namespace SparseSolverPack {

///
/** Concreate sparse solver subclass that uses MA28.
 *
 * ToDo: Finish documentation!
 */
class DirectSparseSolverMA28 : public DirectSparseSolverImp {
public:

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
	/** Stores the factorization structure for MA28
	 */
	class FactorizationStructureMA28 : public FactorizationStructure {
	public:
		// ToDo: Implement this class!
	}; // end class FactorizationStructureMA28

	///
	/** Stores the factorization nonzeros for MA28
	 */
	class FactorizationNonzerosMA28 : public FactorizationNonzeros {
	public:
		// ToDo: Implement this class!
	}; // end class FactorizationNonzerosMA28

	///
	/** Implements the BasisMatrix object for MA28
	 */
	class BasisMatrixMA28 : public BasisMatrixImp {
	public:
		// ToDO: Implement this class!
	}; // end class BasisMatrixMA28

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
		,const BasisMatrix::fact_struc_ptr_t            &fact_struc
		,const BasisMatrixImp::fact_nonzeros_ptr_t      &fact_nonzeros
		,LinAlgPack::IVector                            *row_perm
		,LinAlgPack::IVector                            *col_perm
		,size_type                                      *rank
		,BasisMatrixImp                                 *basis_matrix
		,std::ostream                                   *out            = NULL
		);
	///
	void imp_factor(
		const SparseLinAlgPack::MatrixConvertToSparse   &A
		,const BasisMatrix::fact_struc_ptr_t            &fact_struc
		,const BasisMatrixImp::fact_nonzeros_ptr_t      &fact_nonzeros
		,BasisMatrixImp                                 *basis_matrix
		,std::ostream                                   *out            = NULL
		);

	//@}

};	// end class DirectSparseSolverMA28 

}	// end namespace SparseSolverPack 

#endif	// DIRECT_SPARSE_FORTRAN_COMPATIBLE_SOLVER_MA28_H
