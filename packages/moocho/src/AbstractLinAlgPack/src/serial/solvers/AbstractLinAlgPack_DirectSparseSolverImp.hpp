// /////////////////////////////////////////////////////////////
// DirectSparseSolverImp.h
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

#ifndef	DIRECT_SPARSE_FORTRAN_COMPATIBLE_SOLVER_IMP_H
#define DIRECT_SPARSE_FORTRAN_COMPATIBLE_SOLVER_IMP_H

#include "DirectSparseSolver.h"

namespace SparseSolverPack {

///
/** Implementation node class for \c DirectSparseSolver that takes
 * care of the memory management details.
 *
 * ToDo: Finish documentation!
 *
 * Subclasses must override the following methods which are pure
 * virtual and are not implemented here: \c basis_matrix_factor(),
 * \c estimated_fillin_ratio(), \c create_fact_struc(),
 * \c create_fact_nonzeros(), \c imp_analyze_and_factor() and
 * \c imp_factor().
 */
class DirectSparseSolverImp : public DirectSparseSolver {
public:

	/** @name Overridden from DirectSparseSolver */
	//@{

	///
	void analyze_and_factor(
		const SparseLinAlgPack::MatrixConvertToSparse   &A
		,LinAlgPack::IVector                            *row_perm
		,LinAlgPack::IVector                            *col_perm
		,size_type                                      *rank
		,BasisMatrix                                    *basis_matrix
		,std::ostream                                   *out
		);
	///
	void factor(
		const SparseLinAlgPack::MatrixConvertToSparse   &A
		,BasisMatrix                                    *basis_matrix
		,const BasisMatrix::fact_struc_ptr_t            &fact_struc
		,std::ostream                                   *out
		);
	///
	const BasisMatrix::fact_struc_ptr_t& get_fact_struc() const;
	///
	void set_uninitialized();

	//@}

protected:

	/** @name Protected types */
	//@{

	///
	/** Abstract class for objects that represent the factorization
	  * nonzeros of a particular matrix.
	  *
	  * This storage for the nonzeros can be reused over and over again
	  * for factorizing matrices with the same structure.
	  */
	class FactorizationNonzeros {
	public:
		///
		virtual ~FactorizationNonzeros() {}
	};

	///
	/** Implementation node subclass that combines factorization structure and
	 * nonzeros into a single basis matrix object.
	 *
	 * The only method that subclasses must override for a complete
	 * basis matrix object is:
	 * 
	 * <tt>AbstractLinAlgPack::MatrixNonsingular::V_InvMtV()</tt>
	 */
	class BasisMatrixImp : public BasisMatrix {
	public:

		///
		typedef MemMngPack::ref_count_ptr<FactorizationNonzeros> fact_nonzeros_ptr_t;

		///
		/** Calls <tt>this->initialize()</tt>
		 */
		BasisMatrixImp(
			size_type                      dim
			,const fact_struc_ptr_t        &fact_struc
			,const fact_nonzeros_ptr_t     &fact_nonzeros
			);

		///
		/** Initialize given initialized factorization structure and factorization nonzeros objects.
		 */
		virtual void initialize(
			size_type                      dim
			,const fact_struc_ptr_t        &fact_struc
			,const fact_nonzeros_ptr_t     &fact_nonzeros
			);

		///
		/** Make the basis matrix uninitialized.
		 *
		 * Postconditions:<ul>
		 * <li> <tt>this->get_fact_struc().get() == NULL</tt>
		 * <li> <tt>this->get_fact_nonzeros().get() == NULL</tt>
		 * </ul>
		 */
		virtual void set_uninitialized();
		
		/// Return a smart pointer to the object that represents the factorization structure.
		virtual const fact_nonzeros_ptr_t&  get_fact_nonzeros() const;

		/** @name Overridden from MatrixBase */
		//@{

		///
		size_type rows() const;
		///
		size_type cols() const;

		//@}

		/** @name Overridden from BasisMatrix */
		//@{

		///
		virtual const fact_struc_ptr_t&  get_fact_struc() const;

		//@}

	private:

#ifdef DOXYGEN_COMPILE
		FactorizationStructure    *fact_struc;
		FactorizationNonzeros     *fact_nonzeros;
#else
		size_type                 dim_;
		fact_struc_ptr_t          fact_struc_;
		fact_nonzeros_ptr_t       fact_nonzeros_;
#endif

	}; // end class BasisMatrixImp

	//@}

	/** @name Protected pure virtual methods to be overridden by concrete direct solver subclasses */
	//@{

	///
	/** Create a new, uninitialized \c FactorizationStructure object.
	 */
	virtual const MemMngPack::ref_count_ptr<FactorizationStructure> create_fact_struc() const = 0;

	///
	/** Create a new, uninitialized \c FactorizationNonzeros object.
	 */
	virtual const MemMngPack::ref_count_ptr<FactorizationNonzeros> create_fact_nonzeros() const = 0;

	///
	/** Called to implement the \c analyze_and_factor() without having to worry about
	 * memory mangagment details.
	 *
	 * ToDo: Finish documentation!
	 */	
	virtual void imp_analyze_and_factor(
		const SparseLinAlgPack::MatrixConvertToSparse   &A
		,const BasisMatrix::fact_struc_ptr_t            &fact_struc
		,const BasisMatrixImp::fact_nonzeros_ptr_t      &fact_nonzeros
		,LinAlgPack::IVector                            *row_perm
		,LinAlgPack::IVector                            *col_perm
		,size_type                                      *rank
		,BasisMatrixImp                                 *basis_matrix
		,std::ostream                                   *out            = NULL
		) = 0;

	///
	/** Called to implement the \c analyze_and_factor() without having to worry about
	 * memory mangagment details.
	 *
	 * ToDo: Finish documentation!
	 */	
	virtual void imp_factor(
		const SparseLinAlgPack::MatrixConvertToSparse   &A
		,const BasisMatrix::fact_struc_ptr_t            &fact_struc
		,const BasisMatrixImp::fact_nonzeros_ptr_t      &fact_nonzeros
		,BasisMatrixImp                                 *basis_matrix
		,std::ostream                                   *out            = NULL
		) = 0;

	//@}

private:

#ifdef DOXYGEN_COMPILE
		FactorizationStructure           *fact_struc;
#else
		BasisMatrix::fact_struc_ptr_t    fact_struc_;
#endif

};	// end class DirectSparseSolverImp 

}	// end namespace SparseSolverPack 

#endif	// DIRECT_SPARSE_FORTRAN_COMPATIBLE_SOLVER_IMP_H
