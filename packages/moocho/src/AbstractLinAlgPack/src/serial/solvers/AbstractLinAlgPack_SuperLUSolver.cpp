// ///////////////////////////////////////////////////////////////////
// SuperLUSolver.cpp
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

#include <assert.h>
#include <valarray>

#include "SparseSolverPack/include/SuperLUSolver.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"
#include "dsp_defs.h"
#include "util.h"

namespace SuperLUPack {

class SuperLUSolverImpl;

///
/** Implementation of SuperLUSolver.
 *
 * ToDo: Finish documentation!
 */
class SuperLUSolverImpl : public SuperLUSolver {
public:

	/** @name Public Types */
	//@{

	///
	class FactorizationStructureImpl : public FactorizationStructure {
	public:
		friend class SuperLUSolverImpl;
	private:
		int                   m_;
		int                   n_;
		int                   nz_;
		std::valarray<int>    perm_r_;
		std::valarray<int>    perm_c_;
		std::valarray<int>    etree_;
	};

	///
	class FactorizationNonzerosImpl : public FactorizationNonzeros {
	public:
		friend class SuperLUSolverImpl;
	private:
		SuperMatrix   L_;
		SuperMatrix   U_;
	};

	//@}

	/** @name Overridden from SuperLUSolver */
	//@{

	///
	void analyze_and_factor(
		int                         m
		,int                        n
		,int                        nz
		,const double               a_val[]
		,const int                  a_row_i[]
		,const int                  a_col_ptr[]
		,FactorizationStructure     *fact_struct
		,FactorizationNonzeros      *fact_nonzeros
		,int                        row_perm[]
		,int                        col_perm[]
		,int                        *rank
		);
	///
	void factor(
		int                             m
		,int                            n
		,int                            nz
		,const double                   a_val[]
		,const int                      a_row_i[]
		,const int                      a_col_ptr[]
		,const FactorizationStructure   &fact_struct
		,FactorizationNonzeros          *fact_nonzeros
		);
	///
	void solve(
		const FactorizationStructure    &fact_struct
		,const FactorizationNonzeros    &fact_nonzeros
		,bool                           transp
		,int                            n
		,int                            nrhs
		,double                         rhs[]
		,int                            ldrhs
		) const;

	//@}

}; // end class SuperLUSolver

//
// SuperLUSolver
//

MemMngPack::ref_count_ptr<SuperLUSolver>
SuperLUSolver::create_solver()
{
	return MemMngPack::rcp(new SuperLUSolverImpl());
}

MemMngPack::ref_count_ptr<SuperLUSolver::FactorizationStructure>
SuperLUSolver::create_fact_struct()
{
	return MemMngPack::rcp(new SuperLUSolverImpl::FactorizationStructureImpl());
}

MemMngPack::ref_count_ptr<SuperLUSolver::FactorizationNonzeros>
SuperLUSolver::create_fact_nonzeros()
{
	return MemMngPack::rcp(new SuperLUSolverImpl::FactorizationNonzerosImpl());
}

//
// SuperLUSolverImp
//

// Overridden from SuperLUSolver

void SuperLUSolverImpl::analyze_and_factor(
	int                         m
	,int                        n
	,int                        nz
	,const double               a_val[]
	,const int                  a_row_i[]
	,const int                  a_col_ptr[]
	,FactorizationStructure     *fact_struct
	,FactorizationNonzeros      *fact_nonzeros
	,int                        perm_r[]
	,int                        perm_c[]
	,int                        *rank
	)
{
	using DynamicCastHelperPack::dyn_cast;

	FactorizationStructureImpl
		&fs = dyn_cast<FactorizationStructureImpl>(*fact_struct);
	FactorizationNonzerosImpl
		&fn = dyn_cast<FactorizationNonzerosImpl>(*fact_nonzeros);

	char refact[] = "N";

	// Setup SuperLU stuff
    const int panel_size = sp_ienv(1);
    const int relax      = sp_ienv(2);
	StatInit(panel_size,relax);

	// Resize storage
	fs.m_ = m;
	fs.n_ = n;
	fs.nz_ = nz;
	fs.perm_r_.resize(m);
	fs.perm_c_.resize(n);
	fs.etree_.resize(n);

    // Create matrix A in the format expected by SuperLU
	SuperMatrix A;
	dCreate_CompCol_Matrix(
		&A, m, n, nz
		,const_cast<double*>(a_val)
		,const_cast<int*>(a_row_i)
		,const_cast<int*>(a_col_ptr)
		,NC, D_, GE
		);

	// Get the columm permutations
	int permc_spec = 0; // ToDo: Make this an external parameter
	get_perm_c(permc_spec, &A, &fs.perm_c_[0]);

	// Permute the columns of the matrix
	SuperMatrix AC;
	sp_preorder(refact,&A,&fs.perm_c_[0],&fs.etree_[0],&AC);

	int info = -1;
	dgstrf(
		refact
		,&AC  
		,1.0    /* diag_pivot_thresh */
		,0.0    /* drop_tol */
		,relax
		,panel_size
		,&fs.etree_[0]
		,NULL   /* work */
		,0      /* lwork */
		,&fs.perm_r_[0]
		,&fs.perm_c_[0]
		,&fn.L_
		,&fn.U_
		,&info
		);

	THROW_EXCEPTION(
		info != 0, std::runtime_error
		,"SuperLUSolverImpl::analyze_and_factor(...): Error, dgstrf(...) returned info = " << info
		);

	std::copy( &fs.perm_r_[0], &fs.perm_r_[0] + m, perm_r );
	std::copy( &fs.perm_c_[0], &fs.perm_c_[0] + n, perm_c );
	*rank = n; // We must assume this until I can figure out a way to do better!

	StatFree();

}

void SuperLUSolverImpl::factor(
	int                             m
	,int                            n
	,int                            nz
	,const double                   a_val[]
	,const int                      a_row_i[]
	,const int                      a_col_ptr[]
	,const FactorizationStructure   &fact_struct
	,FactorizationNonzeros          *fact_nonzeros
	)
{
	using DynamicCastHelperPack::dyn_cast;

	const FactorizationStructureImpl
		&fs = dyn_cast<const FactorizationStructureImpl>(fact_struct);
	FactorizationNonzerosImpl
		&fn = dyn_cast<FactorizationNonzerosImpl>(*fact_nonzeros);

	char refact[] = "Y";

	// Setup SuperLU stuff
    const int panel_size = sp_ienv(1);
    const int relax      = sp_ienv(2);
	StatInit(panel_size,relax);

    // Create matrix A in the format expected by SuperLU
	SuperMatrix A;
	dCreate_CompCol_Matrix(
		&A, m, n, nz
		,const_cast<double*>(a_val)
		,const_cast<int*>(a_row_i)
		,const_cast<int*>(a_col_ptr)
		,NC, D_, GE
		);

	// Permute the columns
	SuperMatrix AC;
	sp_preorder(
		refact,&A
		,const_cast<int*>(&fs.perm_c_[0])
		,const_cast<int*>(&fs.etree_[0])
		,&AC
		);

	int info = -1;
	dgstrf(
		refact
		,&AC  
		,1.0    /* diag_pivot_thresh */
		,0.0    /* drop_tol */
		,relax
		,panel_size
		,const_cast<int*>(&fs.etree_[0])
		,NULL   /* work */
		,0      /* lwork */
		,const_cast<int*>(&fs.perm_r_[0])
		,const_cast<int*>(&fs.perm_c_[0])
		,&fn.L_
		,&fn.U_
		,&info
		);

	THROW_EXCEPTION(
		info != 0, std::runtime_error
		,"SuperLUSolverImpl::factor(...): Error, dgstrf(...) returned info = " << info
		);

	StatFree();

}

void SuperLUSolverImpl::solve(
	const FactorizationStructure    &fact_struct
	,const FactorizationNonzeros    &fact_nonzeros
	,bool                           transp
	,int                            n
	,int                            nrhs
	,double                         rhs[]
	,int                            ldrhs
	) const
{

	using DynamicCastHelperPack::dyn_cast;

	const FactorizationStructureImpl
		&fs = dyn_cast<const FactorizationStructureImpl>(fact_struct);
	const FactorizationNonzerosImpl
		&fn = dyn_cast<const FactorizationNonzerosImpl>(fact_nonzeros);

	SuperMatrix B;
    dCreate_Dense_Matrix(&B, n, nrhs, rhs, ldrhs, DN, D_, GE);

	char transc[1];
	transc[0] = ( transp ? 'T' : 'N' );

	int info = -1;
    dgstrs(
		transc
		,const_cast<SuperMatrix*>(&fn.L_)
		,const_cast<SuperMatrix*>(&fn.U_)
		,const_cast<int*>(&fs.perm_r_[0])
		,const_cast<int*>(&fs.perm_c_[0])
		,&B, &info
		);

	THROW_EXCEPTION(
		info != 0, std::runtime_error
		,"SuperLUSolverImpl::solve(...): Error, dgssv(...) returned info = " << info
		);

}

} // end namespace SuperLUPack

#endif // SPARSE_SOLVER_PACK_USE_SUPERLU
