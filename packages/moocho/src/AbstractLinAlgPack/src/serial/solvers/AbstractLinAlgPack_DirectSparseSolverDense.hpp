// /////////////////////////////////////////////////////////////
// DirectSparseSolverDense.h
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

#ifndef	DIRECT_SPARSE_SOLVER_DENSE_H
#define DIRECT_SPARSE_SOLVER_DENSE_H

#include <valarray>
#include <vector>
#include <string>

#include "DirectSparseSolverImp.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/IVector.h"
#include "StandardMemberCompositionMacros.h"

namespace SparseSolverPack {

///
/** Concreate sparse solver subclass that uses the dense LAPACK routines.
 *
 * ToDo: Finish documentation!
 */
class DirectSparseSolverDense : public DirectSparseSolverImp {
public:

	/** @name Constructors/initializers */
	//@{

	///
	/** Default constructor */
	DirectSparseSolverDense();

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
	/** Implements the BasisMatrix object for Dense.
	 */
	class BasisMatrixDense : public BasisMatrixImp {
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

	}; // end class BasisMatrixDense

	///
	/** Stores the factorization structure for Dense
	 */
	class FactorizationStructureDense : public FactorizationStructure {
	public:
		friend class DirectSparseSolverDense;
		friend class BasisMatrixDense;
	private:
		FortranTypes::f_int      m_;         // Number of rows in A
		FortranTypes::f_int      n_;         // Number of columns in A
		FortranTypes::f_int      nz_;        // Number of nonzeros in A
		FortranTypes::f_int      rank_;      // Rank of the basis
		IVector                  col_perm_;  // First rank entries selects the basis of A
		IVector                  inv_col_perm_; // Inverse of col_perm_
		FactorizationStructureDense();
	}; // end class FactorizationStructureDense

	///
	/** Stores the factorization nonzeros for Dense
	 */
	class FactorizationNonzerosDense : public FactorizationNonzeros {
	public:
		typedef FortranTypes::f_int f_int;
		friend class DirectSparseSolverDense;
		friend class BasisMatrixDense;
	private:
		GenMatrix                          LU_;
		bool                               rect_analyze_and_factor_; // true for n > m analyze_and_factor()
		std::valarray<f_int>               ipiv_; // The permutation sent to xGETRS (identity if rect_analyze_and_factor_==true)
		IVector                            basis_perm_; // Only used if rect_analyze_and_factor_==true
	}; // end class FactorizationNonzerosDense

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

};	// end class DirectSparseSolverDense 

}	// end namespace SparseSolverPack 

#endif	// DIRECT_SPARSE_SOLVER_DENSE_H
