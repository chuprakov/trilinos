// //////////////////////////////////////////////////////////////
// ConstraintsRelaxedStd.cpp
//
// ToDo: 12/29/00: Consider scaling when determining if a
// constraint is violated.  We should consider the size of
// ||d(e[j](x))/d(x)||inf but this is expensive to compute
// given the current interfaces.  We really need to rectify
// this!
//

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include <assert.h>

#include <limits>

#include "ConstrainedOptimizationPack/include/ConstraintsRelaxedStd.h"
#include "SparseLinAlgPack/include/GenPermMatrixSliceOp.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/sparse_bounds_diff.h"
#include "LinAlgPack/include/LinAlgOpPack.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StV;
	using SparseLinAlgPack::Vp_StMtV;
}

namespace {

// Find the maxinum element of a dense vector
// and its indice.
std::pair<LinAlgPack::size_type,LinAlgPack::value_type>
imp_max_element( const LinAlgPack::VectorSlice &v )
{
	const LinAlgPack::VectorSlice::const_iterator
		itr = std::max_element( v.begin(), v.end() );
	typedef std::pair<LinAlgPack::size_type,LinAlgPack::value_type> ele_t;
	return ele_t( itr - v.begin() + 1, *itr );
}

// Update a maxinum violation
//
// If max(v(i)) > max_viol then this function
// sets max_viol = v(i) and max_viol_j = i and returns
// true.  Otherwise it returns false.
//
bool imp_update_max_viol(
	  const LinAlgPack::VectorSlice &v
	, LinAlgPack::value_type		*max_viol
	, LinAlgPack::size_type			*max_viol_j
	)
{
	typedef std::pair<LinAlgPack::size_type,LinAlgPack::value_type> ele_t;
	ele_t ele = imp_max_element(v);
	if( ele.second > *max_viol ) {
		*max_viol_j	= ele.first;
		*max_viol 	= ele.second;
		return true;
	}
	return false;
}

// Get an element from a sparse vector and return zero if it does not exist
LinAlgPack::value_type get_sparse_element( const SparseLinAlgPack::SpVectorSlice& v
	, LinAlgPack::size_type i )
{
	const SparseLinAlgPack::SpVectorSlice::element_type
		*ele_ptr = v.lookup_element(i);
	return ele_ptr ? ele_ptr->value() : 0.0;
}

}	// end namespace

namespace ConstrainedOptimizationPack {
namespace QPSchurPack {

// members for ConstraintsRelaxedStd

ConstraintsRelaxedStd::ConstraintsRelaxedStd()
	:
		inequality_pick_policy_(ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY)
		,etaL_(0.0)
		,dL_(NULL)
		,dU_(NULL)
		,eL_(NULL)
		,eU_(NULL)
		,Ed_(NULL)
		,last_added_j_(0)
		,last_added_bound_(0.0)
		,last_added_bound_type_(FREE)
{}

void ConstraintsRelaxedStd::initialize(
	  size_type nd
	, value_type etaL
	, const SpVectorSlice* dL, const SpVectorSlice* dU
	, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
		, const SpVectorSlice* eL, const SpVectorSlice* eU
	, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f	
	, size_type m_undecomp, const size_type j_f_undecomp[]
	, VectorSlice* Ed
	, bool check_F
	, value_type bounds_tol
	, value_type inequality_tol
	, value_type equality_tol
	)
{
	size_type
		m_in = 0,
		m_eq = 0;

	assert( m_undecomp == (F ? f->size() : 0) ); // ToDo: support undecomposed equalities in future.

	// Validate that the correct sets of constraints are selected
	if( dL && !dU )
		throw std::invalid_argument( "ConstraintsRelaxedStd::initialize(...) : Error, "
			"If dL!=NULL then dU!=NULL must also be true." );
	if( E && ( !b || !eL || !eU ) )
		throw std::invalid_argument( "ConstraintsRelaxedStd::initialize(...) : Error, "
			"If E!=NULL then b!=NULL, eL!=NULL and eU!=NULL must also be true." );
	if( F && !f )
		throw std::invalid_argument( "ConstraintsRelaxedStd::initialize(...) : Error, "
			"If F!=NULL then f!=NULL must also be true." );

	// Validate input argument sizes
	if(dL) {
		if( dL->size() != nd )
			throw std::invalid_argument( "ConstraintsRelaxedStd::initialize(...) : Error, "
				"dL.size() != d->size()." );
		if( dU->size() != nd )
			throw std::invalid_argument( "ConstraintsRelaxedStd::initialize(...) : Error, "
				"dU.size() != d->size()." );
	}
	if(E) {
		m_in = BLAS_Cpp::rows( E->rows(), E->cols(), trans_E );
		if( BLAS_Cpp::cols( E->rows(), E->cols(), trans_E )	!= nd )
			throw std::invalid_argument( "ConstraintsRelaxedStd::initialize(...) : Error, "
				"op(E).cols() != nd." );
		if( b->size() != m_in )
			throw std::invalid_argument( "ConstraintsRelaxedStd::initialize(...) : Error, "
				"b->size() != op(E).rows()." );
		if( eL->size() != m_in )
			throw std::invalid_argument( "ConstraintsRelaxedStd::initialize(...) : Error, "
				"eL->size() != op(E).rows()." );
		if( eU->size() != m_in )
			throw std::invalid_argument( "ConstraintsRelaxedStd::initialize(...) : Error, "
				"eU->size() != op(E).rows()." );
		if( Ed ) {
			if( Ed->size() != m_in )
				throw std::invalid_argument( "ConstraintsRelaxedStd::initialize(...) : Error, "
					"Ed->size() != op(E).rows()." );
		}
	}
	if(F) {
		m_eq = BLAS_Cpp::rows( F->rows(), F->cols(), trans_F );
		if( BLAS_Cpp::cols( F->rows(), F->cols(), trans_F )	!= nd )
			throw std::invalid_argument( "QPSolverRelaxed::solve_qp(...) : Error, "
				"op(F).cols() != nd." );
		if( f->size() != m_eq )
			throw std::invalid_argument( "QPSolverRelaxed::solve_qp(...) : Error, "
				"f->size() != op(F).rows()." );
	}
	
	// Initialize other members
	A_bar_.initialize(nd,m_in,m_eq,E,trans_E,b,F,trans_F,f,m_undecomp,j_f_undecomp);
	etaL_			= etaL;
	dL_				= dL;
	dU_				= dU;
	eL_				= eL;
	eU_				= eU;
	Ed_				= Ed;
	check_F_		= check_F;
	bounds_tol_		= bounds_tol;
	inequality_tol_	= inequality_tol;
	equality_tol_	= equality_tol;
	last_added_j_	= 0;	// No cached value.
}

// Overridden from Constraints

size_type ConstraintsRelaxedStd::n() const
{
	return A_bar_.nd() + 1;
}

size_type ConstraintsRelaxedStd::m_breve() const
{
	return A_bar_.m_in() + A_bar_.m_eq();
}

const MatrixWithOp& ConstraintsRelaxedStd::A_bar() const
{
	return A_bar_;
}

void ConstraintsRelaxedStd::pick_violated_policy( EPickPolicy pick_policy )
{
	switch(pick_policy) {
		case ANY_VIOLATED:
			inequality_pick_policy_ = ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY;
			break;
		case MOST_VIOLATED:
			inequality_pick_policy_ = ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY;
			break;
		default:
			assert(0);
	}
}

Constraints::EPickPolicy
ConstraintsRelaxedStd::pick_violated_policy() const
{
	switch(inequality_pick_policy_) {
		case ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY:
			return ANY_VIOLATED;
		case ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY:
			return ANY_VIOLATED;
		case ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY:
			return MOST_VIOLATED;
		default:
			assert(0);
	}
	return ANY_VIOLATED;	// will never be executed
}

void ConstraintsRelaxedStd::pick_violated(
	  const VectorSlice& x, size_type* j_viol, value_type* constr_val
	, value_type* viol_bnd_val, value_type* norm_2_constr, EBounds* bnd, bool* can_ignore
	) const
{
	using LinAlgPack::norm_inf;
	using SparseLinAlgPack::imp_sparse_bnd_diff;

	if( x.size() != A_bar_.nd()+1 ) {
		throw std::length_error( "ConstraintsRelaxedStd::pick_violated(...) : Error, "
			"x is the wrong length" );
	}

	const size_type
		nd = A_bar_.nd();
	const VectorSlice
		d = x(1,nd);
	const value_type
		eta = x(nd+1);

	Vector r;
	bool Ed_computed = false;

	// //////////////////////////////////////////////
	// Check the equality constraints first
	if( check_F_ && A_bar_.F() ) {
		// ToDo: Finish this!
		assert(0);
	}

	// /////////////////////////////////////////////
	// Find the most violated variable bound.

	size_type	max_bound_viol_j		= 0;
	value_type 	max_bounds_viol 		= 0.0;
	bool		max_bound_viol_upper	= false;
	if( dL_ && ( dL_->nz() || dU_->nz() ) ) {
		const value_type scale = 1.0 / (1.0 + norm_inf(d));
		// r = dL - d
		r.resize(nd);
		imp_sparse_bnd_diff( +1, *dL_, BLAS_Cpp::lower, d, &r() );
		imp_update_max_viol( r(), &max_bounds_viol, &max_bound_viol_j );
		// r = d - dU
		imp_sparse_bnd_diff( -1, *dU_, BLAS_Cpp::upper, d, &r() );
		max_bound_viol_upper
			= imp_update_max_viol( r(), &max_bounds_viol, &max_bound_viol_j );

		if( scale * max_bounds_viol > bounds_tol_ )
		{
			*j_viol			= max_bound_viol_j;
			*constr_val		= d(max_bound_viol_j);
			*viol_bnd_val	= max_bound_viol_upper
				? d(max_bound_viol_j) - max_bounds_viol		// e(i) - (e(i) - eU(i)) = eU(i)
				: max_bounds_viol + d(max_bound_viol_j);	// (eL(i) - e(i)) + e(i) = eL(i)
			*norm_2_constr	= 1.0;
			*bnd			= max_bound_viol_upper ? UPPER : LOWER;
			*can_ignore		= false;
		}
		else {
			max_bound_viol_j = 0;	// No variable bounds sufficiently violated.
		}
	}

	if( (	inequality_pick_policy_ == ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY
			||	inequality_pick_policy_ == ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY )
		&& max_bound_viol_j
			 )
	{
		// A variable bound has been violated so lets just return it!
		last_added_j_			= *j_viol;
		last_added_bound_type_	= *bnd;
		last_added_bound_		= *viol_bnd_val;
		return;	
	}

	// /////////////////////////////////////////////
	// Check the general inequalities

	size_type	max_inequality_viol_j		= 0;
	value_type 	max_inequality_viol			= 0.0;
	bool		max_inequality_viol_upper	= false;

	if( inequality_pick_policy_ == ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY ) {
		// Find the first general inequality violated by more than
		// the defined tolerance.
		throw std::logic_error( "ConstraintsRelaxedStd::pick_violated(...) : Error,\n"
			"The option ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY has not been implemented yet\n" );
	}
	else {
		// Find the most violated inequality constraint
		if( A_bar_.m_in() &&  ( eL_->nz() || eU_->nz() ) ) {
			// e = op(E)*d - b*eta
			Vector e;
			LinAlgOpPack::V_MtV( &e, *A_bar_.E(), A_bar_.trans_E(), d );
			if(Ed_) {
				*Ed_ = e;
				Ed_computed = true;
			}
			const value_type scale = 1.0 / (1.0 + norm_inf(e));
			LinAlgPack::Vp_StV( &e(), -eta, *A_bar_.b() );
			// r = eL - e
			r.resize(A_bar_.m_in());
			imp_sparse_bnd_diff( +1, *eL_, BLAS_Cpp::lower, e(), &r() );
			imp_update_max_viol( r(), &max_inequality_viol, &max_inequality_viol_j );
			// r = e - eU
			imp_sparse_bnd_diff( -1, *eU_, BLAS_Cpp::upper, e(), &r() );
			max_inequality_viol_upper
				= imp_update_max_viol( r(), &max_inequality_viol, &max_inequality_viol_j );
			if( max_inequality_viol > max_bounds_viol
				&& scale * max_inequality_viol > inequality_tol_ )
			{
				*j_viol			= max_inequality_viol_j + nd + 1; // offset into A_bar
				*constr_val		= e(max_inequality_viol_j);
				*viol_bnd_val	= max_inequality_viol_upper
					? e(max_inequality_viol_j) - max_inequality_viol	// e(i) - (e(i) - eU(i)) = eU(i)
					: max_inequality_viol + e(max_inequality_viol_j);	// (eL(i) - e(i)) + e(i) = eL(i)
				*norm_2_constr	= 1.0;	// ToDo: Compute it some how?
				*bnd			= max_inequality_viol_upper ? UPPER : LOWER;
				*can_ignore		= false;
			}
			else {
				max_inequality_viol_j = 0;	// No general inequality constraints sufficiently violated.
			}
		}
	}

	if( max_bound_viol_j || max_inequality_viol_j ) {
		// One of the constraints has been violated so just return it.
		last_added_j_			= *j_viol;
		last_added_bound_type_	= *bnd;
		last_added_bound_		= *viol_bnd_val;
		return;		
	}

	// If we get here then no constraint was found that violated any of the tolerances.
	if(Ed_ && !Ed_computed) {
		// Ed = op(E)*d
		LinAlgOpPack::V_MtV( Ed_, *A_bar_.E(), A_bar_.trans_E(), d );
	}
	*j_viol			= 0;
	*constr_val		= 0.0;
	*viol_bnd_val	= 0.0;
	*norm_2_constr	= 0.0;
	*bnd			= FREE;	// Meaningless
	*can_ignore		= false;
}

void ConstraintsRelaxedStd::ignore( size_type j )
{
	throw std::logic_error(  "ConstraintsRelaxedStd::ignore(...) : Error, "
		"This operation is not supported yet!" );
}

value_type ConstraintsRelaxedStd::get_bnd( size_type j, EBounds bnd ) const
{
	const value_type inf = std::numeric_limits<value_type>::max();

	if( j > A_bar_.cols() ) {
		throw std::range_error( "ConstraintsRelaxedStd::get_bnd(j,bnd) : Error, "
			"j is not in range" );
	}

	// See if this is the last constraint we added to the active set.
	if( j == last_added_j_ && bnd == last_added_bound_type_ ) {
		return last_added_bound_;
	}

	// Lookup the bound! (sparse lookup)
	size_type j_local = j;
	const SpVectorSlice::element_type *ele_ptr = NULL;
	if( j_local <= A_bar_.nd() && dL_ ) {
		switch( bnd ) {
			case EQUALITY:
			case LOWER:
				return ( ele_ptr = dL_->lookup_element(j_local) )
					? ele_ptr->value() : -inf;
			case UPPER:
				return ( ele_ptr = dU_->lookup_element(j_local) )
					? ele_ptr->value() : +inf;
			default:
				assert(0);
		}
	}
	else if( (j_local -= A_bar_.nd()) <= 1 ) {
		switch( bnd ) {
			case EQUALITY:
			case LOWER:
				return etaL_;
			case UPPER:
				return +inf;
			default:
				assert(0);
		}
	}
	else if( (j_local -= 1) <= A_bar_.m_in() ) {
		switch( bnd ) {
			case EQUALITY:
			case LOWER:
				return ( ele_ptr = eL_->lookup_element(j_local) )
					? ele_ptr->value() : -inf;
			case UPPER:
				return ( ele_ptr = eU_->lookup_element(j_local) )
					? ele_ptr->value() : +inf;
			default:
				assert(0);
		}
	}
	else if( (j_local -= A_bar_.m_in()) <= A_bar_.m_eq() ) {
		switch( bnd ) {
			case EQUALITY:
			case LOWER:
			case UPPER:
				return -(*A_bar_.f())(j_local);
			default:
				assert(0);
		}
	}

	return 0.0;	// will never be executed!
}

void ConstraintsRelaxedStd::cache_last_added( size_type last_added_j, value_type last_added_bound
	, EBounds last_added_bound_type ) const
{
	last_added_j_			= last_added_j;
	last_added_bound_		= last_added_bound;
	last_added_bound_type_	= last_added_bound_type;
}

// members for ConstraintsRelaxedStd::MatrixConstraints

ConstraintsRelaxedStd::MatrixConstraints::MatrixConstraints()
	:
		nd_(0)
		,m_in_(0)
		,m_eq_(0)
		,E_(NULL)
		,trans_E_(BLAS_Cpp::no_trans)
		,b_(NULL)
		,F_(NULL)
		,trans_F_(BLAS_Cpp::no_trans)
		,f_(NULL)
{}

void ConstraintsRelaxedStd::MatrixConstraints::initialize(
	  size_type				nd
	, size_type				m_in
	, size_type				m_eq
	, const MatrixWithOp	*E
	, BLAS_Cpp::Transp		trans_E
	, const VectorSlice		*b
	, const MatrixWithOp	*F
	, BLAS_Cpp::Transp		trans_F
	, const VectorSlice		*f
	, size_type             m_undecomp
	, const size_type       j_f_undecomp[]
	)
{
	namespace GPMSIP = SparseLinAlgPack::GenPermMatrixSliceIteratorPack;

	// Setup P_u
	const bool test_setup = true; // Todo: Make this an argument!
	P_u_row_i_.resize(m_undecomp);
	P_u_col_j_.resize(m_undecomp);
	if( m_undecomp > 0 ) {
		const size_type
			*j_f_u = j_f_undecomp;
		row_i_t::iterator
			row_i_itr = P_u_row_i_.begin();
		col_j_t::iterator
			col_j_itr = P_u_col_j_.begin();
		for( size_type i = 1; i <= m_undecomp; ++i, ++j_f_u, ++row_i_itr, ++col_j_itr ) {
			*row_i_itr = *j_f_u;
			*col_j_itr = i;
		}					
		P_u_.initialize_and_sort(nd,m_undecomp,m_undecomp,0,0,GPMSIP::BY_ROW
			,&P_u_row_i_[0],&P_u_col_j_[0],test_setup);
	}
	
	// Set the rest of the members
	nd_			= nd;
	m_in_		= m_in;
	m_eq_		= m_eq;
	E_			= E;
	trans_E_	= trans_E;
	b_			= b;
	F_			= F;
	trans_F_	= trans_F;
	f_			= f;
}

// Overridden from Matrix

size_type ConstraintsRelaxedStd::MatrixConstraints::rows() const
{
	return nd_ + 1;
}

size_type ConstraintsRelaxedStd::MatrixConstraints::cols() const
{
	return nd_ + 1 + m_in_ + m_eq_;
}

// Overridden from MatrixWithOp

MatrixWithOp& ConstraintsRelaxedStd::MatrixConstraints::operator=(
	const MatrixWithOp& m)
{
	// ToDo: Finish me
	assert(0);
	return *this;
}

/* 10/25/00 I don't think we need this function yet!
void ConstraintsRelaxedStd::MatrixConstraints::Mp_StPtMtP(
	  GenMatrixSlice* C, value_type a
	, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
	, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
	) const
{
	using BLAS_Cpp::trans_not;
	using SparseLinAlgPack::dot;
	using SparseLinAlgPack::Vp_StMtV;
	using SparseLinAlgPack::Vp_StPtMtV;
	namespace GPMSIP = SparseLinAlgPack::GenPermMatrixSliceIteratorPack;

	//	
	//	A_bar = [  I   0  op(E')   op(F')  ]
	//	        [  0   1   -b'      -f'    ]
	//

	const size_type
		E_start = nd() + 1 + 1,
		F_start = E_start + m_in(),
		F_end	= F_start + m_eq() - 1;
	const Range1D
		d_rng = Range1D(1,nd()),
		E_rng = m_in() ? Range1D(E_start,F_start-1) : Range1D(),
		F_rng = m_eq() ? Range1D(F_start,F_end) : Range1D();

	// For this to work (as shown below) we need to have P1 sorted by
	// row if op(P1) = P1' or sorted by column if op(P1) = P1.
	// Also, we must have P1 sorted by
	// row if op(P2) = P2 or sorted by column if op(P2) = P2'
	// If P1 or P2 are not sorted properly, we will just use the default
	// implementation of this operation.
	if( 	( P1.ordered_by() == GPMSIP::BY_ROW && P1_trans == BLAS_Cpp::no_trans )
	    || 	( P1.ordered_by() == GPMSIP::BY_COL && P1_trans == BLAS_Cpp::trans )
	    || 	( P2.ordered_by() == GPMSIP::BY_ROW && P2_trans == BLAS_Cpp::trans )
	    || 	( P2.ordered_by() == GPMSIP::BY_COL && P2_trans == BLAS_Cpp::no_trans ) )
	{
		// Call the default implementation
		MatrixWithOp::Vp_StPtMtV(C,a,P1,P1_trans,M_trans,P2,P2_trans);
		return;
	}

	if( M_trans == BLAS_Cpp::no_trans ) {
		//
		// C += a * op(P1) * A_bar * op(P2)
		//
		//    = a * [ op(P11)  op(P12) ] * [ I  0  op(E')  op(F') ] * [ op(P21) ]
		//                                 [ 0  1    -b'    -f'   ]   [ op(P22) ]
		//                                                            [ op(P23) ]
		//                                                            [ op(P24) ]
		//
		// C +=   a*op(P11)*op(P21) + a*op(P21)*op(P22)
		//      + a*op(P11)*op(E')*op(P23) - a*op(P12)*b'*op(P23)
		//      + a*op(P11)*op(F')*op(P24) - a*op(P12)*f'*op(P24)
		//      

		assert(0);	// ToDo: Implement This!

	}
	else {
		assert(0);	// ToDo: Finish This!
	}
}
*/

void ConstraintsRelaxedStd::MatrixConstraints::Vp_StMtV(
	  VectorSlice* y, value_type a, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& x, value_type b) const
{

	assert( !F_ || P_u_.cols() == f_->size() ); // ToDo: Add P_u when needed!

	using BLAS_Cpp::trans_not;
	using LinAlgPack::dot;
	using LinAlgPack::Vt_S;
	using LinAlgPack::Vp_StV;
	using SparseLinAlgPack::Vp_StMtV;

	LinAlgOpPack::Vp_MtV_assert_sizes(y->size(),rows(),cols(),trans_rhs1,x.size());

	//	
	//	A_bar = [  I   0  op(E')   op(F')  ]
	//	        [  0   1   -b'      -f'    ]
	//

	const size_type
		E_start = nd() + 1 + 1,
		F_start = E_start + m_in(),
		F_end	= F_start + m_eq() - 1;
	const Range1D
		d_rng = Range1D(1,nd()),
		E_rng = m_in() ? Range1D(E_start,F_start-1) : Range1D(),
		F_rng = m_eq() ? Range1D(F_start,F_end) : Range1D();

	// y = b * y
	if( b == 0.0 )
		*y = 0.0;
	else
		Vt_S( y, b );
	
	if( trans_rhs1 == BLAS_Cpp::no_trans ) {
		//
		// y += a* A_bar * x
		// 
		//   += a * [ I  0  op(E')  op(F') ] * [ x1 ]
		//          [ 0  1   -b'     -f'   ]   [ x2 ]
		//                                     [ x3 ]
		//                                     [ x4 ]
		//
		// [ y1 ]  += [ a * x1 + a * op(E') * x3 + a * op(F') * x4 ]
		// [ y2 ]     [ a * x2 - a * b' * x3     - a * f' * x4     ]
		//
		VectorSlice
			y1 = (*y)(d_rng);
		value_type
			&y2 = (*y)(nd()+1);
		const VectorSlice
			x1 = x(d_rng);
		const value_type
			x2 = x(nd()+1);
		const VectorSlice
			x3 = m_in() ? x(E_rng) : VectorSlice(),
			x4 = m_eq() ? x(F_rng) : VectorSlice();
		
		// [ y1 ]  += [ a * x1 + a * op(E') * x3 + a * op(F') * x4 ]
		Vp_StV( &y1, a, x1 );
		if( m_in() )
			Vp_StMtV( &y1, a, *E(), trans_not( trans_E() ), x3 );
		if( m_eq() )
			Vp_StMtV( &y1, a, *F(), trans_not( trans_F() ), x4 );
		// [ y2 ]  += [ a * x2 - a * b' * x3     - a * f' * x4     ]
		y2 += a * x2;
		if( m_in() )
			y2 += - a * dot( *this->b(), x3 );
		if( m_eq() )
			y2 += - a * dot( *f(), x4 );
	}
	else if ( trans_rhs1 == BLAS_Cpp::trans ) {
		//
		// y += a* A_bar' * x
		// 
		//   += a * [ I       0 ] * [ x1 ]
		//          [ 0       1 ]   [ x2 ]
		//          [ op(E)  -b ]
		//          [ op(F)  -f ]
		//
		// [ y1 ]    [ a * x1                        ]
		// [ y2 ]    [                + a * x2       ]
		// [ y3 ] += [ a * op(E) * x1 - a * b * x2   ]
		// [ y4 ]    [ a * op(F) * x1 - a * f * x2   ]
		//
		VectorSlice
			y1 = (*y)(d_rng);
		value_type
			&y2 = (*y)(nd()+1);
		VectorSlice
			y3 = m_in() ? (*y)(E_rng) : VectorSlice(),
			y4 = m_eq() ? (*y)(F_rng) : VectorSlice();
		const VectorSlice
			x1 = x(d_rng);
		const value_type
			x2 = x(nd()+1);
		// y1 += a * x1
		Vp_StV( &y1, a, x1 );
		// y2 += a * x2
		y2 += a * x2;
		// y3 += a * op(E) * x1 - (a*x2) * b
		if( m_in() ) {
			Vp_StMtV( &y3, a, *E(), trans_E(), x1 );
			Vp_StV( &y3, - a * x2, *this->b() );
		}
		// y4 += a * op(F) * x1 - (a*x2) * f
		if( m_eq() ) {
			Vp_StMtV( &y4, a, *F(), trans_F(), x1 );
			Vp_StV( &y4, - a * x2, *f() );
		}
	}
	else {
		assert(0);	// Invalid trans value
	}
}

void ConstraintsRelaxedStd::MatrixConstraints::Vp_StPtMtV(
	  VectorSlice* y, value_type a
	, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp M_trans
	, const SpVectorSlice& x, value_type beta) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;

	assert( !F_ || P_u_.cols() == f_->size() ); // ToDo: Add P_u when needed!

	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
 	using SparseLinAlgPack::dot;
	using SparseLinAlgPack::Vp_StMtV;
	using SparseLinAlgPack::Vp_StPtMtV;
	namespace GPMSIP = SparseLinAlgPack::GenPermMatrixSliceIteratorPack;

	LinAlgOpPack::Vp_MtV_assert_sizes(y->size(),P.rows(),P.cols(),P_trans
		, BLAS_Cpp::rows( rows(), cols(), M_trans) );
	LinAlgOpPack::Vp_MtV_assert_sizes( BLAS_Cpp::cols( P.rows(), P.cols(), P_trans)
		,rows(),cols(),M_trans,x.size());

	//	
	//	A_bar = [  I   0  op(E')   op(F')  ]
	//	        [  0   1   -b'      -f'    ]
	//

	const size_type
		E_start = nd() + 1 + 1,
		F_start = E_start + m_in(),
		F_end	= F_start + m_eq() - 1;
	const Range1D
		d_rng = Range1D(1,nd()),
		E_rng = m_in() ? Range1D(E_start,F_start-1) : Range1D(),
		F_rng = m_eq() ? Range1D(F_start,F_end) : Range1D();

	// For this to work (as shown below) we need to have P sorted by
	// row if op(P) = P' or sorted by column if op(P) = P.  If
	// P is not sorted properly, we will just use the default
	// implementation of this operation.
	if( 	( P.ordered_by() == GPMSIP::BY_ROW && P_trans == BLAS_Cpp::no_trans )
	    || 	( P.ordered_by() == GPMSIP::BY_COL && P_trans == BLAS_Cpp::trans ) )
	{
		// Call the default implementation
		MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,M_trans,x,beta);
		return;
	}

	if( M_trans == BLAS_Cpp::no_trans ) {
		//
		// y = beta*y + a * op(P) * A_bar * x
		//
		//   = beta * y
		//   
		//    + a * [op(P1)  op(P2) ] * [ I  0  op(E')  op(F') ] * [ x1 ]
		//                              [ 0  1   -b'     -f'   ]   [ x2 ]
		//                                                         [ x3 ]
		//                                                         [ x4 ]
		//
		// y = beta*y + a*op(P1)*x1 + a*op(P1)*op(E')*x3 + a*op(P1)*op(F')*x4
		//     + a*op(P2)*x2 - a*op(P2)*b'*x3 - a*op(P2)*f'*x4
		//
		// Where:
		//   op(P1) = op(P)(:,1:nd)
		//   op(P2) = op(P)(:,nd+1:nd+1)
		//

		const GenPermMatrixSlice
			P1 = ( P.is_identity() 
				   ? GenPermMatrixSlice( nd(), nd(), GenPermMatrixSlice::IDENTITY_MATRIX )
				   : P.create_submatrix(d_rng,P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL)
				),
			P2 = ( P.is_identity()
				   ? GenPermMatrixSlice(
					   P_trans == no_trans ? nd() : 1
					   , P_trans == no_trans ? 1 : nd()
					   , GenPermMatrixSlice::ZERO_MATRIX )
				   : P.create_submatrix(Range1D(nd()+1,nd()+1),P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL)
				);

		const SpVectorSlice
			x1 = x(d_rng);
		const value_type
			x2 = get_sparse_element(x,nd()+1);
		const SpVectorSlice
			x3 = m_in() ? x(E_rng) : SpVectorSlice(NULL,0,0,0),
			x4 = m_eq() ? x(F_rng) : SpVectorSlice(NULL,0,0,0);

		// y = beta*y + a*op(P1)*x1
		Vp_StMtV( y, a, P1, P_trans, x1, beta );
		// y += a*op(P1)*op(E')*x3
		if( m_in() )
			Vp_StPtMtV( y, a, P1, P_trans, *E(), trans_not(trans_E()), x3 );
		// y += a*op(P1)*op(F')*x4
		if( m_eq() )
			Vp_StPtMtV( y, a, P1, P_trans, *F(), trans_not(trans_F()), x4 );
		//
		// y += a*op(P2)*x2 - a*op(P2)*b'*x3 - a*op(P2)*f'*x4
		//   += a * op(P2) * ( x2 + b'*x3 - f'*x4 )
		//   
		//   ==>
		//   
		// y(i) +=  a * ( x2 - b'*x3 - f'*x4 )
		//   
		if( P2.nz() ){
			assert(P2.nz() == 1);
			const size_type
				i = P_trans == BLAS_Cpp::no_trans
					? P2.begin()->row_i() : P2.begin()->col_j();
			value_type
				&y_i = (*y)(i) += a * x2;
			if(m_in())
				y_i += -a * dot(*this->b(),x3);
			if(m_eq())
				y_i += -a * dot(*this->f(),x4);
		}
	}
	else if ( M_trans == BLAS_Cpp::trans ) {
		//
		// y = beta*y + a*op(P)*A_bar'*x
		// 
		//   = beta*y
		//   
		//    + a * [ P1  P2  P3  P4 ] * [ I       0 ] * [ x1 ]
		//                               [ 0       1 ]   [ x2 ]
		//                               [ op(E)  -b ]
		//                               [ op(F)  -f ]
		//
		// y = beta*y + a*P1*x1 + a*P2*x2 + a*P3*op(E)*x1 - a*P3*b*x2
		//     + a*P4*op(F)*x1 - a*P4*f*x2
		//
		// Where:
		//   P1 = op(P)(:,1:nd)
		//   P2 = op(P)(:,nd+1:nd+1)
		//   P3 = op(P)(:,nd+1+1:nd+1+m_in)
		//   P4 = op(P)(:,nd+1+m_in+1:nd+1+m_in+m_eq)
		//

		assert( !P.is_identity() ); // We should never have this!

		size_type off = 0;
		const GenPermMatrixSlice
			P1 = P.create_submatrix(Range1D(off+1,off+nd())
									,P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL);
		off += nd();
		const GenPermMatrixSlice
			P2 = P.create_submatrix(Range1D(off+1,off+1)
									,P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL);
		off += 1;
		const GenPermMatrixSlice
			P3 = m_in()
				? P.create_submatrix(Range1D(off+1,off+m_in())
									 ,P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL)
				: GenPermMatrixSlice();
		off += m_in();
		const GenPermMatrixSlice
			P4 = m_eq()
				? P.create_submatrix(Range1D(off+1,off+m_eq())
									 ,P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL)
				: GenPermMatrixSlice();

		const SpVectorSlice
			x1 = x(d_rng);
		const value_type
			x2 = get_sparse_element(x,nd()+1);

		// y = beta*y + a*op(P1)*x1
		Vp_StMtV( y, a, P1, P_trans, x1, beta );
		// y += a*op(P2)*x2
		if( P2.nz() ){
			assert(P2.nz() == 1);
			(*y)( P_trans == BLAS_Cpp::no_trans
					? P2.begin()->row_i() : P2.begin()->col_j() )
				+= a * x2;
		}
		if(m_in()) {
			// y += a*P3*op(E)*x1
			Vp_StPtMtV( y, a, P3, P_trans, *E(), trans_E(), x1 );
			// y += (-a*x2)*P3*b
			Vp_StMtV( y, - a * x2, P3, P_trans, *this->b() );
		}
		if(m_eq()) {
			// y += a*P4*op(F)*x1
			Vp_StPtMtV( y, a, P4, P_trans, *F(), trans_F(), x1 );
			// y += (-a*x2)*P4*f
			Vp_StMtV( y, - a * x2, P4, P_trans, *this->f() );
		}
		
	}
	else {
		assert(0);	// Invalid trans value
	}
}

}	// end namespace QPSchurPack 
}	// end namespace ConstrainedOptimizationPack 
