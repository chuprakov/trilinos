// /////////////////////////////////////////////////////////////
// DirectSparseSolverMA28.cpp
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

#ifdef SPARSE_SOLVER_PACK_USE_MA28

#include <assert.h>

#include <ostream>

#include "SparseSolverPack/include/DirectSparseSolverMA28.h"
#include "SparseSolverPack/include/MatrixScaling_Strategy.h"
#include "AbstractFactoryStd.h"
#include "ThrowException.h"
#include "dynamic_cast_verbose.h"

#ifdef SPARSE_SOLVER_PACK_USE_MC29

// Declarations to call MC29AD to scale matrix.
namespace MC29AD_C_Decl {

typedef FortranTypes::f_int			f_int;
typedef FortranTypes::f_dbl_prec	f_dbl_prec;

extern "C" {
	FORTRAN_FUNC_DECL_UL(void,MC29AD,mc29ad) ( const f_int& M, const f_int& N
		, const f_int& NE, const f_dbl_prec A[], const f_int IRN[]
		, const f_int ICN[], f_dbl_prec R[], f_dbl_prec C[]
		, f_dbl_prec W[], const f_int& LP, f_int* IFAIL );
}

}  // end namespace MC29AD_C_Decl

namespace {

typedef FortranTypes::f_int			f_int;
typedef FortranTypes::f_dbl_prec	f_dbl_prec;

inline 
void mc29ad( const f_int& m, const f_int& n
	, const f_int& ne, const f_dbl_prec a[], const f_int irn[]
	, const f_int icn[], f_dbl_prec r[], f_dbl_prec c[]
	, f_dbl_prec w[], const f_int& lp, f_int* ifail )
{
	MC29AD_C_Decl::FORTRAN_FUNC_CALL_UL(MC29AD,mc29ad)(m,n,ne,a,irn,icn,r,c,w,lp,ifail);
}

}	// end namespace

#endif // SPARSE_SOLVER_PACK_USE_MC29

namespace SparseSolverPack {

// //////////////////////////////////////////////////
// DirectSparseSolverMA28::FactorizationStructureMA28

DirectSparseSolverMA28::FactorizationStructureMA28::FactorizationStructureMA28()
	:m_(0),n_(0),max_n_(0),nz_(0),licn_(0),lirn_(0),iflag_(0)
	 ,u_(0.1),scaling_(NO_SCALING)
	 ,cntr_var_changed_(false)
{}

// //////////////////////////////////////////////////
// DirectSparseSolverMA28::BasisMatrixMA28

// Overridden from BasisMatrixImp

MemMngPack::ref_count_ptr<DirectSparseSolverImp::BasisMatrixImp>
DirectSparseSolverMA28::BasisMatrixMA28::create_matrix() const
{
	return MemMngPack::rcp(new BasisMatrixMA28);
}

void DirectSparseSolverMA28::BasisMatrixMA28::V_InvMtV(
	VectorWithOpMutable* y, BLAS_Cpp::Transp trans_rhs1, const VectorWithOp& x
	) const 
{
	using DynamicCastHelperPack::dyn_cast;

	// Get the concrete factorization and nonzeros objects
	const FactorizationStructureMA28
		&fs = dyn_cast<const FactorizationStructureMA28>(*this->get_fact_struc());
	const FactorizationNonzerosMA28
		&fn = dyn_cast<const FactorizationNonzerosMA28>(*this->get_fact_nonzeros());

	// Validate input
#ifdef _DEBUG
	THROW_EXCEPTION(
		fs.cntr_var_changed_, std::logic_error
		,"DirectSparseSolverMA28::BasisMatrixMA28::V_InvMtV(...), Error, "
		"You can't change control variables between calls of \'analyze_and factorize\' and \'solve\'");	
#endif

	assert(0); // Update the below code!

/*
	// Check for consistancy
	const size_type y_dim = y->dim(); x_dim = x.dim();
	THROW_EXCEPTION(
		fs.m_ != 
	if(m != m_ || n != n_ )
		throw InputException("m, and n must be the same as in the call to \'analyze_and_factorize\'");

	// Allocate workspace memory
	std::valarray<f_dbl_prec>	w(max_n_);

	f_int mtype = ( (trans_matrix == BLAS_Cpp::no_trans) ? 1 : 0 );

	// Copy rhs into x since ma28 operates this way, but only if they are different vectors.
	if(rhs != x)
		std::copy(rhs, rhs+max_n_, x);

	// Scale the rhs
	scale_rhs( trans_matrix, x ); 

	// Solve for the rhs
	ma28_.ma28cd(max_n_, &a_[0], licn_, &icn_[0], &ikeep_[0], x, &w[0], mtype);

	// Scale the lhs
	scale_lhs( trans_matrix, x ); 

*/
	assert(0); // ToDo: Implement!
}

// //////////////////////////////////////////////////
// DirectSparseSolverMA28

// Constructors/initializers

DirectSparseSolverMA28::DirectSparseSolverMA28()
	: estimated_fillin_ratio_(10.0)
{}

// Overridden from DirectSparseSolver

const DirectSparseSolver::basis_matrix_factory_ptr_t
DirectSparseSolverMA28::basis_matrix_factory() const
{
	namespace mmp = MemMngPack;
	return mmp::rcp(new mmp::AbstractFactoryStd<BasisMatrix,BasisMatrixMA28>());
}

void DirectSparseSolverMA28::estimated_fillin_ratio(
	value_type estimated_fillin_ratio
	)
{
	estimated_fillin_ratio_ = estimated_fillin_ratio;
}

// Overridden from DirectSparseSolverImp

const MemMngPack::ref_count_ptr<DirectSparseSolver::FactorizationStructure>
DirectSparseSolverMA28::create_fact_struc() const
{
	return MemMngPack::rcp(new FactorizationStructureMA28);
}

const MemMngPack::ref_count_ptr<DirectSparseSolverImp::FactorizationNonzeros>
DirectSparseSolverMA28::create_fact_nonzeros() const
{
	return MemMngPack::rcp(new FactorizationNonzerosMA28);
}

void DirectSparseSolverMA28::imp_analyze_and_factor(
	const SparseLinAlgPack::MatrixConvertToSparse   &A
	,const BasisMatrix::fact_struc_ptr_t            &fact_struc
	,const BasisMatrixImp::fact_nonzeros_ptr_t      &fact_nonzeros
	,LinAlgPack::IVector                            *row_perm
	,LinAlgPack::IVector                            *col_perm
	,size_type                                      *rank
	,std::ostream                                   *out
	)
{
	using DynamicCastHelperPack::dyn_cast;
	typedef MatrixConvertToSparse MCTS;

	if(out)
		*out << "\nUsing MA28 to analyze and factor a new matrix ...\n";

	// Get the concrete factorization and nonzeros objects
	FactorizationStructureMA28
		&fs = dyn_cast<FactorizationStructureMA28>(*fact_struc);
	FactorizationNonzerosMA28
		&fn = dyn_cast<FactorizationNonzerosMA28>(*fact_nonzeros);

	// Get the dimensions of things.
	const index_type
		m = A.rows(),
		n = A.cols(),
		nz = A.num_nonzeros( MCTS::EXTRACT_FULL_MATRIX ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM );

	// Validate input
	THROW_EXCEPTION(
		n <= 0 || m <= 0 || m > n, std::invalid_argument
		,"DirectSparseSolverMA28::imp_analyze_and_factor(...) : Error!" );

	// Sets control variable change flag to false.
	fs.cntr_var_changed_ = false;

	// Memorize the dimenstions for checks later
	fs.m_ = m; fs.n_ = n; fs.nz_ = nz;
	fs.max_n_ = std::_MAX(fs.m_,fs.n_);

	// By default set licn and ircn equal to fillin_ratio * nz.
	if( fs.licn_ < fs.nz_ ) fs.licn_ = estimated_fillin_ratio_ * fs.nz_;
	if( fs.lirn_ < fs.nz_ ) fs.lirn_ = estimated_fillin_ratio_ * fs.nz_;

	// Initialize matrix factorization storage and temporary storage
	fs.irn_.resize(fs.nz_);   // perminatly stores nz row indexes (needed for matrix scaling)
	fs.icn_.resize(fs.licn_); // first nz entries stores the column indexes
	fn.a_.resize(fs.licn_);
	fs.ikeep_.resize(5*fs.max_n_);
	std::valarray<index_type>   irn_tmp_(fs.lirn_), iw(8*fs.max_n_);  // ToDo: use Workspace<>
	std::valarray<value_type>   w(fs.max_n_);     // ToDo: use Workspace<>

	// Fill in the matrix information
	A.coor_extract_nonzeros(
		MCTS::EXTRACT_FULL_MATRIX
		,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM
		,fs.nz_
		,&fn.a_[0]
		,fs.nz_
		,&fs.irn_[0]
		,&fs.icn_[0]
		);
	std::copy( &fs.irn_[0], &fs.irn_[0] + fs.nz_, &irn_tmp_[0] );

	// Scale the matrix
	if( fs.matrix_scaling_.get() )
		fs.matrix_scaling_->scale_matrix(
			fs.m_, fs.n_, fs.nz_, &fs.irn_[0], &fs.icn_[0], true
			,&fn.a_[0]
			);

	// Factor the matrix
	fs.ma28_.ma28ad(
		fs.max_n_, fs.nz_, &fn.a_[0], fs.licn_, &irn_tmp_[0], fs.lirn_, &fs.icn_[0], fs.u_
		,&fs.ikeep_[0], &iw[0], &w[0], fs.iflag_
		);

	// Todo : deal with resizing for insufficient storage when needed.  Write a loop!

	// Extract the basis matrix selection
	assert(0); // ToDo: Implement!

	// Check for errors and throw exception if you have to.
	ThrowIFlagException(fs.iflag_);

}

void DirectSparseSolverMA28::imp_factor(
	const SparseLinAlgPack::MatrixConvertToSparse   &A
	,const BasisMatrix::fact_struc_ptr_t            &fact_struc
	,const BasisMatrixImp::fact_nonzeros_ptr_t      &fact_nonzeros
	,std::ostream                                   *out
	)
{
	using DynamicCastHelperPack::dyn_cast;
	if(out)
		*out << "\nUsing MA28 to analyze and factor a new matrix ...\n";
	// Get the concrete factorization and nonzeros objects
	const FactorizationStructureMA28
		&fs = dyn_cast<const FactorizationStructureMA28>(*fact_struc);
	FactorizationNonzerosMA28
		&fn = dyn_cast<FactorizationNonzerosMA28>(*fact_nonzeros);
	// Extract the nonzeros out of A
	assert(0); // ToDo: Implement!
	// ...
}

// private

void DirectSparseSolverMA28::ThrowIFlagException(index_type iflag)
{
	E_IFlag e_iflag = static_cast<E_IFlag>(iflag);
	const char msg_err_head[] = "DirectSparseSolverMA28::ThrowIFlagException(iflag) : Error";
	switch(e_iflag) {
		case SLOW_ITER_CONV :
			THROW_EXCEPTION(
				true, std::runtime_error
				,msg_err_head << ", Convergence to slow" );
		case MAXIT_REACHED :
			THROW_EXCEPTION(
				true, std::runtime_error
				,msg_err_head << ", Maximum iterations exceeded");
		case MA28BD_CALLED_WITH_DROPPED :
			THROW_EXCEPTION(
				true, std::logic_error
				,msg_err_head << ", ma28bd called with elements dropped in ma28ad");
		case DUPLICATE_ELEMENTS :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", Duplicate elements have been detected");
		case NEW_NONZERO_ELEMENT :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", A new non-zero element has be passed to ma28bd that was not ot ma28ad");
		case N_OUT_OF_RANGE :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", 1 <=max(n,m) <= 32767 has been violated");
		case NZ_LE_ZERO :
			THROW_EXCEPTION(
				true, std::logic_error
				,msg_err_head << ", nz <= 0 has been violated");
		case LICN_LE_NZ :
			THROW_EXCEPTION(
				true, std::logic_error
				,msg_err_head << ", licn <= nz has been violated");
		case LIRN_LE_NZ :
			THROW_EXCEPTION(
				true, std::logic_error
				,msg_err_head << ", lirn <= nz has been violated");
		case ERROR_DURRING_BLOCK_TRI :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", An error has occured durring block triangularization");
		case LICN_AND_LIRN_TOO_SMALL :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", licn and lirn are to small to hold matrix factorization");
		case LICN_TOO_SMALL :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", licn is to small to hold matrix factorization");
		case LICN_FAR_TOO_SMALL :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", licn is to far small to hold matrix factorization");
		case LIRN_TOO_SMALL :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", lirn is to small to hold matrix factorization");
		case NUMERICALLY_SINGULAR :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", matrix is numerically singular, see \'abort2\'");
		case STRUCTURALLY_SINGULAR :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", matrix is structurally singular, see \'abort1\'");
		default:
			return; // We don't throw exceptions for other values of iflag.
	}
}

}	// end namespace SparseSolverPack 

#endif // SPARSE_SOLVER_PACK_USE_MA28
