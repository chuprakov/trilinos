// //////////////////////////////////////////////////////////////////////////////////
// MatrixVarReductImplicit.cpp

#include "ConstrainedOptimizationPack/include/MatrixVarReductImplicit.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/GenPermMatrixSlice.h"
#include "SparseLinAlgPack/include/dense_Vp_StPtMtV.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/WorkspacePack.h"

namespace {

//
// Implicit matrix vector multiplication:
//
// y = b*y + a*op(inv(C)*N)*x
//
template<class V>
void imp_Vp_StMtV_implicit(
	LinAlgPack::VectorSlice                                             *y
	,LinAlgPack::value_type                                             a
	,const ConstrainedOptimizationPack::DecompositionSystemVarReduct    &decomp_sys
	,BLAS_Cpp::Transp                                                   D_trans
	,const V                                                            &x
	,LinAlgPack::value_type                                             b
	)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	const LinAlgPack::size_type
		r   = decomp_sys.C().rows(),
		dof = decomp_sys.N().cols();

	if( D_trans == no_trans ) {
		//
		// y += a * inv(C) * ( N * x )
		//
		wsp::Workspace<LinAlgPack::value_type> t1_ws(wss,r), t2_ws(wss,r);
		LinAlgPack::VectorSlice t1(&t1_ws[0],t1.size()), t2(&t2_ws[0],t2.size());
		// t1 = N*x
		LinAlgOpPack::V_MtV( &t1, decomp_sys.N(), no_trans, x );
		// t2 = inv(C) * t1
		decomp_sys.solve_C( t1, no_trans, &t2 );
		// y = b*y
		if(b==0.0)       *y = 0.0;
		else if(b!=1.0)  LinAlgPack::Vt_S(y,b);
		// y += a*t2
		LinAlgPack::Vp_StV( y, -a, t2 );
	}
	else {
		//
		// y = b*y + a * N' * ( inv(C') * x )
		//
		wsp::Workspace<LinAlgPack::value_type> t1_ws(wss,r);
		LinAlgPack::VectorSlice t1(&t1_ws[0],t1.size());
		// t1 = inv(C')*x
		decomp_sys.solve_C( x, trans, &t1 );
		// y = b*y + a*N'*t1
	    SparseLinAlgPack::Vp_StMtV( y, a,  decomp_sys.N(), trans, t1, b );
	}
}

//
// Generate a row of inv(C)*N if not already computed.
//
void imp_compute_InvCtN_row(
	LinAlgPack::size_type                                                 r
	,const ConstrainedOptimizationPack::DecompositionSystemVarReduct      &decomp_sys
	,LinAlgPack::size_type                                                j
	,LinAlgPack::VectorSlice                                              *e_j  // Set to all zeros on input and output!
	,ConstrainedOptimizationPack::MatrixVarReductImplicit::InvCtN_rows_t  *InvCtN_rows
	)
{
	typedef  LinAlgPack::value_type  value_type;
	using LinAlgPack::VectorSlice;
	if( (*InvCtN_rows)[j-1] == NULL ) {
		// Generate row j of inv(C)*N
		value_type *vec = (*InvCtN_rows)[j-1] = new value_type[r]; // ToDo: We may want to allocate more vectors at once!
		VectorSlice row_j(vec,r);
		// row_j = N'*inv(C')*e_j
		(*e_j)(j) = 1.0;
		imp_Vp_StMtV_implicit( &row_j, 1.0, decomp_sys, BLAS_Cpp::trans, *e_j, 0.0 );
		(*e_j)(j) = 0.0;
	}
}

//
// Perform the matrix-vector multiplication:
// 
// y = b*y -a * op(P) * [inv(C) * N] * x
//
// by generating rows [inv(C)*N](j,:) for each nonzero entry op(P)(i,j).
//
template<class V>
void imp_Vp_StPtMtV_by_row(
	LinAlgPack::VectorSlice                                               *y
	,LinAlgPack::value_type                                               a
	,const ConstrainedOptimizationPack::GenPermMatrixSlice                &P
	,BLAS_Cpp::Transp                                                     P_trans
	,const ConstrainedOptimizationPack::DecompositionSystemVarReduct      &decomp_sys
	,const V                                                              &x
	,LinAlgPack::value_type                                               b
	,ConstrainedOptimizationPack::MatrixVarReductImplicit::InvCtN_rows_t *InvCtN_rows
	)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	using LinAlgPack::dot;
	using LinAlgPack::VectorSlice;
	using SparseLinAlgPack::dot;
	using SparseLinAlgPack::Vp_StMtV;
	using SparseLinAlgPack::GenPermMatrixSlice;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	
	const LinAlgPack::size_type
		D_rows = decomp_sys.C().rows(),
		D_cols = decomp_sys.N().cols();
	// y = b*y
	if(b==0.0)       *y = 0.0;
	else if(b!=1.0)  Vt_S(y,b);
	// Compute t = N'*inv(C')*e(j) then y(i) += -a*t'*x where op(P)(i,j) = 1.0
	wsp::Workspace<LinAlgPack::value_type>   e_j_ws(wss,D_rows);
	VectorSlice                              e_j(&e_j_ws[0],e_j_ws.size());
	e_j = 0.0;
	for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
		const LinAlgPack::size_type
			i = P_trans == no_trans ? itr->row_i() : itr->col_j(),
			j = P_trans == no_trans ? itr->col_j() : itr->row_i();
		// t = op(M') * e(j)
		imp_compute_InvCtN_row(D_rows,decomp_sys,j,&e_j,InvCtN_rows);
		VectorSlice t((*InvCtN_rows)[j-1],D_cols);
		// y(i) += -a*t'*x
		(*y)(i) += (-a) * dot(t,x);
	}
}

} // end namespace

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

namespace ConstrainedOptimizationPack {

MatrixVarReductImplicit::MatrixVarReductImplicit()
	: decomp_sys_(NULL), use_dense_mat_vec_(true)
{} // Every other member will initialize themselfs

void MatrixVarReductImplicit::initialize(
	const DecompositionSystemVarReduct     *decomp_sys
	,const GenMatrixSlice                  *D_dense
	,const release_resource_ptr_t          &release_resource_ptr
	)
{
	// Validate the inputs
	if( !decomp_sys )
		throw std::invalid_argument(
			"MatrixVarReductImplicit::initialize(...): Error, "
			"decomp_sys must not be NULL" );
	if( D_dense && (D_dense->rows() != decomp_sys->C().rows() || D_dense->cols() != decomp_sys->N().cols() ) )
		throw std::invalid_argument(
			"MatrixVarReductImplicit::initialize(...): Error, "
			"*D_dense does not match the size of the decomposition in *decomp_sys" );
	// Set the members
	decomp_sys_ = decomp_sys;
	if( D_dense )
		D_dense_.bind( const_cast<GenMatrixSlice&>(*D_dense) );
	use_dense_mat_vec_ = true; // ToDo: We should use a timer to determine if this is faster than sparse or not!
	release_resource_ptr_ = release_resource_ptr;
	if(InvCtN_rows_.size()) { // Free previously allocated vectors
		for( InvCtN_rows_t::iterator itr = InvCtN_rows_.begin(); itr != InvCtN_rows_.end(); ++itr ) {
			if( *itr )
				delete [] *itr; // ToDo: We may want to allocate vectors in larger chuncks
			*itr = (value_type*)NULL;
		}
	}
}

void MatrixVarReductImplicit::set_uninitialized()
{
	decomp_sys_            = NULL;
	D_dense_               = GenMatrixSlice();
	use_dense_mat_vec_     = true;
	release_resource_ptr_  = NULL;
	InvCtN_rows_.resize(0);
}

const DecompositionSystemVarReduct& MatrixVarReductImplicit::decomp_sys() const
{
	assert_initialized();
	return *decomp_sys_;
}

const MatrixVarReductImplicit::release_resource_ptr_t&
MatrixVarReductImplicit::release_resource_ptr() const
{
	assert_initialized();
	return release_resource_ptr_;
}

// Overridden from Matrix

size_type MatrixVarReductImplicit::rows() const
{
	return decomp_sys_ ? decomp_sys_->C().rows() : 0;
}

size_type MatrixVarReductImplicit::cols() const
{
	return decomp_sys_ ? decomp_sys_->N().cols() : 0;
}

// Overridden from MatrixWithOp

MatrixWithOp& MatrixVarReductImplicit::operator=(const MatrixWithOp& M)
{
	assert(0); // Finish!
	return *this;
}

void MatrixVarReductImplicit::MatrixVarReductImplicit::Vp_StMtV(
	VectorSlice* y, value_type a, BLAS_Cpp::Transp D_trans
	, const VectorSlice& x, value_type b) const
{
	assert_initialized();
	LinAlgPack::Vp_MtV_assert_sizes(y->size(),rows(),cols(),D_trans,x.size());
	if( use_dense_mat_vec_ && D_dense_.rows() > 0 ) {
		LinAlgOpPack::Vp_StMtV( y, a, D_dense_, D_trans, x, b );
	}
	else {
		imp_Vp_StMtV_implicit( y, a, *decomp_sys_, D_trans, x, b );
	}
}

void MatrixVarReductImplicit::Vp_StMtV(
	VectorSlice* y, value_type a, BLAS_Cpp::Transp D_trans
	, const SpVectorSlice& x, value_type b) const
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	assert_initialized();
	const size_type
		D_rows = this->rows(), D_cols = this->cols(),
		opD_rows = rows(D_rows,D_cols,D_trans), opD_cols = cols(D_rows,D_cols,D_trans);
	LinAlgPack::Vp_MtV_assert_sizes(y->size(),D_rows,D_cols,D_trans,x.size());
	if( use_dense_mat_vec_ && D_dense_.rows() > 0 ) {
		LinAlgOpPack::Vp_StMtV( y, a, D_dense_, D_trans, x, b );
	}
	else {
		if( x.nz() == x.size() ) {  // This is B.S.  Should use MatrixWithOpFactorized object for C!
			VectorSlice dx = SparseLinAlgPack::dense_view(x);
			imp_Vp_StMtV_implicit( y, -a, *decomp_sys_, D_trans, dx, b );
		}
		else if( D_trans == BLAS_Cpp::trans && x.nz() < D_cols ) {
			//
			// y = b*y + (-a)*[N'*inv(C')]*x
			//
			// We can do something crafty here.  We can generate columns of N'*inv(C')
			// and then perform y += -a*[N'*inv(C')](:,j)*x(j) for nonzero x(j)
			//
			wsp::Workspace<LinAlgPack::value_type>   e_j_ws(wss,D_rows);
			VectorSlice                              e_j(&e_j_ws[0],e_j_ws.size());
			e_j = 0.0;
			// y = b*y
			if(b==0.0)       *y = 0.0;
			else if(b!=1.0)  Vt_S(y,b);
			// y += -a*[N'*inv(C')](:,j)*x(j), for j <: { j | x(j) != 0.0 }
			const SpVectorSlice::difference_type o = x.offset();
			for( SpVectorSlice::const_iterator itr = x.begin(); itr != x.end(); ++itr ) {
				const size_type j = itr->indice() + o;
				imp_compute_InvCtN_row(D_rows,*decomp_sys_,j,&e_j,&InvCtN_rows_);
				LinAlgPack::Vp_StV( y, -a * itr->value(), VectorSlice(InvCtN_rows_[j-1],D_cols) );
			}
		}
		else {   // This is B.S.  Should use MatrixWithOpFactorized object for C!
			Vector dx;
			LinAlgOpPack::assign( &dx, x );
			imp_Vp_StMtV_implicit( y, -a, *decomp_sys_, D_trans, dx(), b );
		}
	}
}

void MatrixVarReductImplicit::Vp_StPtMtV(
	VectorSlice* y, value_type a
	, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp D_trans
	, const VectorSlice& x, value_type b) const
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	assert_initialized();
	const size_type
		D_rows = this->rows(), D_cols = this->cols(),
		opD_rows = rows(D_rows,D_cols,D_trans), opD_cols = cols(D_rows,D_cols,D_trans);
	LinAlgPack::Vp_MtV_assert_sizes(y->size(),P.rows(),P.cols(),P_trans,opD_rows);
	LinAlgPack::Vp_MtV_assert_sizes(cols(P.rows(),P.cols(),P_trans),D_rows,D_cols,D_trans,x.size());
	if( D_dense_.rows() > 0 ) {
		dense_Vp_StPtMtV(y,a,P,P_trans,D_dense_,D_trans,x,b);
	}
	else if( P.nz() > D_cols || D_trans == BLAS_Cpp::trans ) {
		// Just use the default implementation
		MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,D_trans,x,b);
	}
	else {
		imp_Vp_StPtMtV_by_row(y,a,P,P_trans,*decomp_sys_,x,b,&InvCtN_rows_);
	}
}

void MatrixVarReductImplicit::Vp_StPtMtV(
	VectorSlice* y, value_type a
	, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp D_trans
	, const SpVectorSlice& x, value_type b) const
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	assert_initialized();
	const size_type
		D_rows = this->rows(), D_cols = this->cols(),
		opD_rows = rows(D_rows,D_cols,D_trans), opD_cols = cols(D_rows,D_cols,D_trans);
	LinAlgPack::Vp_MtV_assert_sizes(y->size(),P.rows(),P.cols(),P_trans,opD_rows);
	LinAlgPack::Vp_MtV_assert_sizes(cols(P.rows(),P.cols(),P_trans),D_rows,D_cols,D_trans,x.size());
	if( D_dense_.rows() > 0 ) {
		dense_Vp_StPtMtV(y,a,P,P_trans,D_dense_,D_trans,x,b);
	}
	else if( P.nz() > D_cols || D_trans == BLAS_Cpp::trans ) {
		// Just use the default implementation
		MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,D_trans,x,b);
	}
	else {
		imp_Vp_StPtMtV_by_row(y,a,P,P_trans,*decomp_sys_,x,b,&InvCtN_rows_);
	}
}

// Private member functions

void MatrixVarReductImplicit::assert_initialized() const
{
	if( !decomp_sys_ )
		throw std::logic_error(
			"MatrixVarReductImplicit::assert_initialized(): Error, "
			"initialize(...) not called" );
}

}	// end namespace ConstrainedOptimizationPack 
