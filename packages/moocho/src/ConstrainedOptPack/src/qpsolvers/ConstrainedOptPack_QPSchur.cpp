// //////////////////////////////////////////////////////////////
// QPSchur.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include <assert.h>

#include <ostream>
#include <iomanip>

#include "ConstrainedOptimizationPack/include/QPSchur.h"
#include "ConstrainedOptimizationPack/include/MatrixSymAddDelUpdateable.h"
#include "ConstrainedOptimizationPack/include/ComputeMinMult.h"
#include "SparseLinAlgPack/include/MatrixWithOpFactorized.h"
#include "SparseLinAlgPack/include/MatrixWithOpOut.h"
#include "SparseLinAlgPack/include/SpVectorOut.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/GenPermMatrixSliceOp.h"
#include "SparseLinAlgPack/include/GenPermMatrixSliceOut.h"
#include "SparseLinAlgPack/include/EtaVector.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/VectorOut.h"
#include "LinAlgPack/include/GenMatrixOut.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StV;
	using SparseLinAlgPack::Vp_StMtV;
}

namespace {

// Some local helper functions.

// Print a bnd as a string
inline
const char* bnd_str( ConstrainedOptimizationPack::QPSchurPack::EBounds bnd ) {
	switch(bnd) {
		case ConstrainedOptimizationPack::QPSchurPack::FREE:
			return "FREE";
		case ConstrainedOptimizationPack::QPSchurPack::UPPER:
			return "UPPER";
		case ConstrainedOptimizationPack::QPSchurPack::LOWER:
			return "LOWER";
		case ConstrainedOptimizationPack::QPSchurPack::EQUALITY:
			return "EQUALITY";
	}
	assert(0);	// should never be executed
	return 0;
}

// print a bool
inline
const char* bool_str( bool b ) {
	return b ? "true" : "false";
}

// Deincrement all indices less that k_remove
void deincrement_indices(
	LinAlgPack::size_type k_remove
	,std::vector<LinAlgPack::size_type> *indice_vector
	,size_t len_vector
	)
{
	typedef LinAlgPack::size_type				size_type;
	typedef std::vector<LinAlgPack::size_type>	vec_t;
	assert( len_vector <= indice_vector->size() );
	for( vec_t::iterator itr = indice_vector->begin(); itr != indice_vector->begin() + len_vector; ++itr ) {
		if( *itr > k_remove )
			--(*itr);
	}
}

// Insert the element (r_v,c_v) into r[] and c[] sorted by r[]
void insert_pair_sorted(
	LinAlgPack::size_type  r_v
	,LinAlgPack::size_type c_v
	,size_t len_vector                       // length of the new vector
	,std::vector<LinAlgPack::size_type> *r
	,std::vector<LinAlgPack::size_type> *c
	)
{
	typedef std::vector<LinAlgPack::size_type> rc_t;
	assert( r->size() >= len_vector && c->size() >= len_vector );
	// find the insertion point in r[]
	rc_t::iterator
		itr = std::lower_bound( r->begin(), r->begin() + len_vector-1, r_v );
	const LinAlgPack::size_type p = itr - r->begin();
	// Shift all of the stuff out of the way to make room for the insert
	{for( rc_t::iterator itr_last = r->begin() + len_vector-1;
			itr_last > r->begin() + p; --itr_last )
	{
		*itr_last = *(itr_last-1);
	}}
	{for( rc_t::iterator itr_last = c->begin() + len_vector-1;
			itr_last > c->begin() + p; --itr_last )
	{
		*itr_last = *(itr_last-1);
	}}
	// Insert the new elements
	(*r)[p] = r_v;
	(*c)[p] = c_v;
}
	
// z_hat = inv(S_hat) * ( d_hat - U_hat'*vo )
void calc_z( const SparseLinAlgPack::MatrixFactorized& S_hat
	, const LinAlgPack::VectorSlice& d_hat, const SparseLinAlgPack::MatrixWithOp& U_hat
	, const LinAlgPack::VectorSlice& vo, LinAlgPack::VectorSlice* z_hat )
{
	using SparseLinAlgPack::Vp_StMtV;
	using SparseLinAlgPack::V_InvMtV;
	LinAlgPack::Vector tmp = d_hat;
	Vp_StMtV( &tmp(), -1.0, U_hat, BLAS_Cpp::trans, vo );
	V_InvMtV( z_hat, S_hat, BLAS_Cpp::no_trans, tmp() );	
}

// v = inv(Ko) * ( fo - U_hat * z_hat )
void calc_v( const SparseLinAlgPack::MatrixFactorized& Ko
	, const LinAlgPack::VectorSlice& fo, const SparseLinAlgPack::MatrixWithOp& U_hat
	, const LinAlgPack::VectorSlice& z_hat, LinAlgPack::VectorSlice* v )
{
	using SparseLinAlgPack::Vp_StMtV;
	using SparseLinAlgPack::V_InvMtV;
	LinAlgPack::Vector tmp = fo; // must use temporary storage.
	Vp_StMtV( &tmp(), -1.0, U_hat, BLAS_Cpp::no_trans, z_hat );
	V_InvMtV( v, Ko, BLAS_Cpp::no_trans, tmp() );	
}

// mu_D_hat =
// 		- Q_XD_hat' * g
// 		- Q_XD_hat' * G * x
// 		- Q_XD_hat' * A * v(n_R+1:n_R+m)
// 		- Q_XD_hat' * A_bar * P_plus_hat * z_hat
void calc_mu_D(
	  const ConstrainedOptimizationPack::QPSchur::ActiveSet& act_set
	, const LinAlgPack::VectorSlice& x
	, const LinAlgPack::VectorSlice& v
	, LinAlgPack::VectorSlice* mu_D )
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using LinAlgOpPack::V_MtV;
	using LinAlgOpPack::V_StMtV;
	using SparseLinAlgPack::V_MtV;
	using SparseLinAlgPack::Vp_MtV;
	using SparseLinAlgPack::Vp_StPtMtV;

	const ConstrainedOptimizationPack::QPSchurPack::QP
		&qp = act_set.qp();
	const LinAlgPack::size_type
		n = qp.n(),
		n_R = qp.n_R(),
		m = qp.m();

	const SparseLinAlgPack::GenPermMatrixSlice &Q_XD_hat = act_set.Q_XD_hat();
	const LinAlgPack::VectorSlice			 g = qp.g();
	const SparseLinAlgPack::MatrixSymWithOp &G = qp.G();
	// mu_D_hat = - Q_XD_hat' * g
	V_StMtV( mu_D, -1.0, Q_XD_hat, trans, g ); 
	// mu_D_hat += - Q_XD_hat' * G * x
	Vp_StPtMtV( mu_D, -1.0, Q_XD_hat, trans, G, no_trans, x );
	// mu_D_hat += - Q_XD_hat' * A * v(n_R+1:n_R+m)
	if( m ) {
		Vp_StPtMtV( mu_D, -1.0, Q_XD_hat, trans, qp.A(), no_trans, v(n_R+1,n_R+m) );
	}
	// p_mu_D_hat += - Q_XD_hat' * A_bar * P_plus_hat * z_hat
	if( act_set.q_plus_hat() && act_set.q_hat() ) {
		const LinAlgPack::VectorSlice z_hat = act_set.z_hat();
		SparseLinAlgPack::SpVector P_plus_hat_z_hat;
		V_MtV( &P_plus_hat_z_hat, act_set.P_plus_hat(), no_trans, z_hat ); 
		Vp_StPtMtV( mu_D, -1.0, Q_XD_hat, trans
			, qp.constraints().A_bar(), no_trans, P_plus_hat_z_hat() );
	}
}

// p_mu_D_hat =
// 		- Q_XD_hat' * G * Q_R * p_v(1:n_R)
// 		- Q_XD_hat' * G * P_XF_hat * p_z_hat
// 		- Q_XD_hat' * A * p_v(n_R+1:n_R+m)
// 		- Q_XD_hat' * A_bar * P_plus_hat * p_z_hat
void calc_p_mu_D(
	  const ConstrainedOptimizationPack::QPSchur::ActiveSet& act_set
	, const LinAlgPack::VectorSlice& p_v
	, const LinAlgPack::VectorSlice& p_z_hat
	, LinAlgPack::VectorSlice* p_mu_D )
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using SparseLinAlgPack::V_MtV;
	using SparseLinAlgPack::Vp_StPtMtV;

	const ConstrainedOptimizationPack::QPSchurPack::QP
		&qp = act_set.qp();
	const LinAlgPack::size_type
		n = qp.n(),
		n_R = qp.n_R(),
		m = qp.m();

	const SparseLinAlgPack::GenPermMatrixSlice &Q_XD_hat = act_set.Q_XD_hat();
	const SparseLinAlgPack::MatrixSymWithOp &G = qp.G();
	// p_mu_D_hat = - Q_XD_hat' * G * Q_R * p_v(1:n_R)
	{
		SparseLinAlgPack::SpVector Q_R_p_v1;
		V_MtV( &Q_R_p_v1, qp.Q_R(), no_trans, p_v(1,n_R) ); 
		Vp_StPtMtV( p_mu_D, -1.0, Q_XD_hat, trans, G, no_trans, Q_R_p_v1(), 0.0 );
	}
	// p_mu_D_hat += - Q_XD_hat' * G * P_XF_hat * p_z_hat
	if( act_set.q_F_hat() ) {
		SparseLinAlgPack::SpVector P_XF_hat_p_z_hat;
		V_MtV( &P_XF_hat_p_z_hat, act_set.P_XF_hat(), no_trans, p_z_hat ); 
		Vp_StPtMtV( p_mu_D, -1.0, Q_XD_hat, trans, G, no_trans, P_XF_hat_p_z_hat() );
	}
	// p_mu_D_hat += - Q_XD_hat' * A * p_v(n_R+1:n_R+m)
	if( m ) {
		Vp_StPtMtV( p_mu_D, -1.0, Q_XD_hat, trans, qp.A(), no_trans, p_v(n_R+1,n_R+m) );
	}
	// p_mu_D_hat += - Q_XD_hat' * A_bar * P_plus_hat * p_z_hat
	if( act_set.q_plus_hat() ) {
		SparseLinAlgPack::SpVector P_plus_hat_p_z_hat;
		V_MtV( &P_plus_hat_p_z_hat, act_set.P_plus_hat(), no_trans, p_z_hat ); 
		Vp_StPtMtV( p_mu_D, -1.0, Q_XD_hat, trans
			, qp.constraints().A_bar(), no_trans, P_plus_hat_p_z_hat() );
	}
}

}	// end namespace

namespace ConstrainedOptimizationPack {

// public member functions for QPSchur::U_hat_t

QPSchur::U_hat_t::U_hat_t()
	:
		 G_(NULL)
		,A_(NULL)
		,A_bar_(NULL)
		,Q_R_(NULL)
		,P_XF_hat_(NULL)
		,P_plus_hat_(NULL)
{}

void QPSchur::U_hat_t::initialize( 
	 const MatrixSymWithOp		*G
	,const MatrixWithOp			*A
	,const MatrixWithOp			*A_bar
	,const GenPermMatrixSlice	*Q_R
	,const GenPermMatrixSlice	*P_XF_hat
	,const GenPermMatrixSlice	*P_plus_hat
	)
{
	G_				= G;
	A_				= A;
	A_bar_			= A_bar;
	Q_R_			= Q_R;
	P_XF_hat_		= P_XF_hat;
	P_plus_hat_		= P_plus_hat;
}

size_type QPSchur::U_hat_t::rows() const
{
	return Q_R_->cols() + ( A_ ? A_->cols() : 0 );
}

size_type QPSchur::U_hat_t::cols() const
{
	return P_plus_hat_->cols();
}

/* 10/25/00: I don't think we need this function!
void QPSchur::U_hat_t::Mp_StM(GenMatrixSlice* C, value_type a
	, BLAS_Cpp::Transp M_trans ) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using LinAlgOpPack::Mp_StMtP;
	using LinAlgOpPack::Mp_StPtMtP;

	// C += a * op(U_hat)

	LinAlgOpPack::Mp_M_assert_sizes( C->rows(), C->cols(), no_trans
		, rows(), cols(), M_trans );

	const size_type
		n_R	= Q_R_->cols(),
		m 	= A() ? A()->cols()  : 0;

	if( M_trans == no_trans ) {
		//
		// C += a * op(U_hat)
		// 
		//    = a * [  Q_R' * G * P_XF_hat + Q_R' * A_bar * P_plus_hat ]
		//          [                  A' * P_XF_hat                   ]
		//          
		// C1 += a * Q_R' * G * P_XF_hat + a * Q_R' * A_bar * P_plus_hat
		// 
		// C2 += a * A' * P_XF_hat
		//
		GenMatrixSlice
			C1 = (*C)(1,n_R,1,C->cols()),
			C2 = m ? (*C)(n_R+1,n_R+m,1,C->cols()) : GenMatrixSlice();
		// C1 += a * Q_R' * G * P_XF_hat
		if( P_XF_hat().nz() )
			Mp_StPtMtP( &C1, a, Q_R(), trans, G(), no_trans, P_XF_hat(), no_trans );
		// C1 += a * Q_R' * A_bar * P_plus_hat
		if( P_plus_hat().nz() )
			Mp_StPtMtP( &C1, a, Q_R(), trans, A_bar(), no_trans, P_plus_hat(), no_trans );
		// C2 += a * A' * P_XF_hat
		if(m && P_XF_hat().nz())
			Mp_StMtP( &C2, a, *A(), trans, P_plus_hat(), no_trans );
	}
	else {
		assert(0);	// Implement this!
	}
}
*/

void QPSchur::U_hat_t::Vp_StMtV(VectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
	, const VectorSlice& x, value_type b) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using LinAlgPack::Vt_S;
	using SparseLinAlgPack::V_MtV;
	using SparseLinAlgPack::Vp_StMtV;
	using SparseLinAlgPack::Vp_StPtMtV;

	LinAlgOpPack::Vp_MtV_assert_sizes(y->size(),rows(),cols(),M_trans,x.size());

	//
	// U_hat = [  Q_R' * G * P_XF_hat + Q_R' * A_bar * P_plus_hat ]
	//         [                  A' * P_XF_hat                   ]
	//         

	const size_type
		n_R	= Q_R_->cols(),
		m 	= A() ? A()->cols()  : 0;

	if( M_trans == BLAS_Cpp::no_trans ) {
		//
		// y =  b*y + a * U_hat * x
		// 
		//   = b*y + a * [  Q_R' * G * P_XF_hat + Q_R' * A_bar * P_plus_hat ] * x
		//               [                  A' * P_XF_hat                   ]
		// 
		//  =>
		// 
		// y1 = b * y1 + a * Q_R' * G * P_XF_hat * x + a * Q_R' * A_bar * P_plus_hat * x
		// y2 = b * y2 + a * A' * P_XF_hat * x
		// 
		VectorSlice
			y1 = (*y)(1,n_R),
			y2 = m ? (*y)(n_R+1,n_R+m) : VectorSlice();
		SpVector
			P_XF_hat_P,
			P_plus_hat_P;
		// P_XF_hat_P = P_XF_hat * x
		if( P_XF_hat().nz() )
			V_MtV( &P_XF_hat_P, P_XF_hat(), no_trans, x );
		// P_plus_hat_P = P_plus_hat * x
		if(P_plus_hat().nz())
			V_MtV( &P_plus_hat_P, P_plus_hat(), no_trans, x );
		// y1 = b * y1
		if(b==0.0)      y1=0.0;
		else if(b!=1.0) Vt_S(&y1,b);
		// y1 += a * Q_R' * G * P_XF_hat_P
		if(P_XF_hat().nz())
			Vp_StPtMtV( &y1, a, Q_R(), trans, G(), no_trans, P_XF_hat_P() );
		// y1 += a * Q_R' * A_bar * P_plus_hat_P
		if(P_plus_hat().nz())
			Vp_StPtMtV( &y1, a, Q_R(), trans, A_bar(), no_trans, P_plus_hat_P() );
		if(m) {
			// y2 = b * y2
			if(b==0.0)      y2=0.0;
			else if(b!=1.0) Vt_S(&y2,b);
			// y2 +=  a * A' * P_XF_hat_P
			if( P_XF_hat().nz() )
				Vp_StMtV( &y2, a, *A(), trans, P_XF_hat_P() );
		}
	}
	else if( M_trans == BLAS_Cpp::trans ) {
		//
		// y =  b*y + a * U_hat' * x
		// 
		//   = b*y + a * [  P_XF_hat' * G * Q_R + P_plus_hat' * A_bar' * Q_R, P_XF_hat' * A ] * [ x1 ]
		//                                                                                      [ x2 ]
		//  =>
		// 
		// y = b * y + a * P_XF_hat' * G * Q_R * x1 + a * P_plus_hat' * A_bar' * Q_R * x1
		//     + a * P_XF_hat' * A * x2
		// 
		const VectorSlice
			x1 = x(1,n_R),
			x2 = m ? x(n_R+1,n_R+m) : VectorSlice();
		SpVector
			Q_R_x1;
		// Q_R_x1 = Q_R * x1
		V_MtV( &Q_R_x1, Q_R(), no_trans, x1 );
		// y = b*y
		if(b==0.0)      *y = 0.0;
		else if(b!=1.0) Vt_S( y, b );
		// y += a * P_XF_hat' * G * Q_R_x1
		if(P_XF_hat().nz())
			Vp_StPtMtV( y, a, P_XF_hat(), trans, G(), no_trans, Q_R_x1() );
		// y += a * P_plus_hat' * A_bar' * Q_R_x1
		if(P_plus_hat().nz())
			Vp_StPtMtV( y, a, P_plus_hat(), trans, A_bar(), trans, Q_R_x1() );
		// y += a * P_XF_hat' * A * x2
		if( m && P_XF_hat().nz() )
			Vp_StPtMtV( y, a, P_XF_hat(), trans, *A(), no_trans, x2 );
	}
	else {
		assert(0);	// Invalid value for M_trans
	}
}

void QPSchur::U_hat_t::Vp_StMtV(VectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
	, const SpVectorSlice& x, value_type b) const
{
//	// Uncomment to use the default version
//	MatrixWithOp::Vp_StMtV(y,a,M_trans,x,b); return;

	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using LinAlgPack::Vt_S;
	using LinAlgOpPack::V_MtV;
	using SparseLinAlgPack::V_MtV;
	using SparseLinAlgPack::Vp_StMtV;
	using SparseLinAlgPack::Vp_StPtMtV;

	LinAlgOpPack::Vp_MtV_assert_sizes(y->size(),rows(),cols(),M_trans,x.size());

	//
	// U_hat = [  Q_R' * G * P_XF_hat + Q_R' * A_bar * P_plus_hat ]
	//         [                  A' * P_XF_hat                   ]
	//         

	const size_type
		n_R	= Q_R_->cols(),
		m 	= A() ? A()->cols()  : 0;

	if( M_trans == BLAS_Cpp::no_trans ) {
		//
		// y =  b*y + a * U_hat * x
		// 
		//   = b*y + a * [  Q_R' * G * P_XF_hat + Q_R' * A_bar * P_plus_hat ] * x
		//               [                  A' * P_XF_hat                   ]
		// 
		//  =>
		// 
		// y1 = b * y1 + a * Q_R' * G * P_XF_hat * x + a * Q_R' * A_bar * P_plus_hat * x
		// y2 = b * y2 + a * A' * P_XF_hat * x
		// 
		VectorSlice
			y1 = (*y)(1,n_R),
			y2 = m ? (*y)(n_R+1,n_R+m) : VectorSlice();
		SpVector
			P_XF_hat_P,
			P_plus_hat_P;
		// P_XF_hat_P = P_XF_hat * x
		if( P_XF_hat().nz() )
			V_MtV( &P_XF_hat_P, P_XF_hat(), no_trans, x );
		// P_plus_hat_P = P_plus_hat * x
		if(P_plus_hat().nz())
			V_MtV( &P_plus_hat_P, P_plus_hat(), no_trans, x );
		// y1 = b * y1
		if(b==0.0)      y1=0.0;
		else if(b!=1.0) Vt_S(&y1,b);
		// y1 += a * Q_R' * G * P_XF_hat_P
		if(P_XF_hat().nz())
			Vp_StPtMtV( &y1, a, Q_R(), trans, G(), no_trans, P_XF_hat_P() );
		// y1 += a * Q_R' * A_bar * P_plus_hat_P
		if(P_plus_hat().nz())
			Vp_StPtMtV( &y1, a, Q_R(), trans, A_bar(), no_trans, P_plus_hat_P() );
		if(m) {
			// y2 = b * y2
			if(b==0.0)      y2=0.0;
			else if(b!=1.0) Vt_S(&y2,b);
			// y2 += a * A' * P_XF_hat_P
			if(P_XF_hat().nz())
				Vp_StMtV( &y2, a, *A(), trans, P_XF_hat_P() );
		}
	}
	else if( M_trans == BLAS_Cpp::trans ) {
		//
		// y =  b*y + a * U_hat' * x
		// 
		//   = b*y + a * [  P_XF_hat' * G * Q_R + P_plus_hat' * A_bar' * Q_R, P_XF_hat' * A ] * [ x1 ]
		//                                                                                      [ x2 ]
		//  =>
		// 
		// y = b * y + a * P_XF_hat' * G * Q_R * x1 + a * P_plus_hat' * A_bar' * Q_R * x1
		//     + a * P_XF_hat' * A * x2
		// 
		const SpVectorSlice
			x1 = x(1,n_R),
			x2 = m ? x(n_R+1,n_R+m) : SpVectorSlice(NULL,0,0,0);
		SpVector
			Q_R_x1;
		// Q_R_x1 = Q_R * x1
		V_MtV( &Q_R_x1, Q_R(), no_trans, x1 );
		// y = b*y
		if(b ==0.0)     *y = 0.0;
		else if(b!=1.0) Vt_S( y, b );
		// y += a * P_XF_hat' * G * Q_R_x1
		if(P_XF_hat().nz())
			Vp_StPtMtV( y, a, P_XF_hat(), trans, G(), no_trans, Q_R_x1() );
		// y += a * P_plus_hat' * A_bar' * Q_R_x1
		if(P_plus_hat().nz())
			Vp_StPtMtV( y, a, P_plus_hat(), trans, A_bar(), trans, Q_R_x1() );
		// y += a * P_XF_hat' * A * x2
		if( m && P_XF_hat().nz() )
			Vp_StPtMtV( y, a, P_XF_hat(), trans, *A(), no_trans, x2 );
	}
	else {
		assert(0);	// Invalid value for M_trans
	}
}

// public member functions for QPSchur::ActiveSet

QPSchur::ActiveSet::ActiveSet(const schur_comp_ptr_t& schur_comp)
	:
		schur_comp_(schur_comp)
		,initialized_(false)
		,test_(false)
		,qp_(NULL)
		,x_init_(NULL)
		,n_(0)
		,n_R_(0)
		,m_(0)
		,q_plus_hat_(0)
		,q_F_hat_(0)
		,q_C_hat_(0)
{}

void QPSchur::ActiveSet::initialize(
	  QP& qp, size_type num_act_change, const int ij_act_change[]
	, const QPSchurPack::EBounds bnds[], bool test
	, std::ostream *out, EOutputLevel output_level )
{
	using LinAlgOpPack::V_mV;
	using SparseLinAlgPack::V_MtV;
	using SparseLinAlgPack::V_InvMtV;
	using SparseLinAlgPack::Vp_StPtMtV;
	using SparseLinAlgPack::Mp_StPtMtP;
	using SparseLinAlgPack::M_StMtInvMtM;
	using LinAlgPack::sym;
	namespace GPMSTP = SparseLinAlgPack::GenPermMatrixSliceIteratorPack;
	
	const size_type
		n		= qp.n(),
		n_R		= qp.n_R(),
		n_X		= n - n_R,
		m		= qp.m();
	const QP::x_init_t
		&x_init = qp.x_init();
	const QP::l_x_X_map_t
		&l_x_X_map = qp.l_x_X_map();
	const QP::i_x_X_map_t
		&i_x_X_map = qp.i_x_X_map();
	const VectorSlice
		b_X = qp.b_X();
	const VectorSlice
		g = qp.g();
	const MatrixSymWithOp
		&G = qp.G();
	const QP::Constraints
		&constraints = qp.constraints();
	const size_type
		m_breve	= constraints.m_breve();

	try {

	// Count the number of each type of change
	size_type 
		q_plus_hat		= 0,
		q_F_hat			= 0,
		q_C_hat			= 0;
	if( num_act_change ) {
		for( size_type k = 1; k <= num_act_change; ++k ) {
			const int ij = ij_act_change[k-1];
			const QPSchurPack::EBounds bnd = bnds[k-1];
			if( ij < 0 ) {
				// Initially fixed variable being freed.
				if( x_init(-ij) == QPSchurPack::FREE ) {
					std::ostringstream omsg;
					omsg
						<< "QPSchur::ActiveSet::initialize(...) : Error, "
						<< "The variable x(" << -ij << ") is not initially fixed and can not "
						<< "be freed by ij_act_change[" << k-1 << "]\n";
					throw std::invalid_argument( omsg.str() );
				}
				if( x_init(-ij) == QPSchurPack::EQUALITY ) {
					std::ostringstream omsg;
					omsg
						<< "QPSchur::ActiveSet::initialize(...) : Error, "
						<< "The variable x(" << -ij << ") is equality fixed and therefore can not "
						<< "be freed by ij_act_change[" << k-1 << "]\n";
					throw std::invalid_argument( omsg.str() );
				}
				++q_F_hat;
			}
			else {
				// Constraint being added to the active-set
				if( ij <= n ) {
					// Fixing a variable to a bound
					QPSchurPack::EBounds x_init_bnd = x_init(ij);
					if( x_init_bnd == QPSchurPack::FREE ) {
						// initially free variable being fixed
						++q_plus_hat;
					}
					else if ( x_init_bnd == QPSchurPack::EQUALITY ) {
						// ToDo: Throw exception
						assert(0);
					}
					else if( x_init_bnd == bnd ) {
						// ToDo: Throw exception
						assert(0);
					}
					else {
						// Initially fixed variable being fixed to another bound
						++q_F_hat;	// We must free the variable first
						++q_C_hat;	// Now we fix it to a different bound.
					}
				}
				else {
					// Adding a general inequality (or equality) constraint
					if( ij > n + m_breve ) {
						// ToDo: Throw exception
						assert(0);
					}		
					++q_plus_hat;
				}
			}
		} 
	}

	const size_type
		q_D_hat = (n - n_R) - q_F_hat,
		q_D_hat_max = n_X;

	// Now let's set stuff up: ij_map, constr_norm, bnds and part of d_hat
	const size_type
		q_hat = q_plus_hat + q_F_hat + q_C_hat,
		q_hat_max = n_X + n,	// If all the initially fixed variables where freed
								// Then all the degrees of freedom where used up with other constraints.
		q_F_hat_max = n_X,
		q_C_hat_max = n_X,
		q_plus_hat_max = n;
		
	ij_map_.resize(q_hat_max);
	constr_norm_.resize(q_hat_max);
	bnds_.resize(q_hat_max);
	d_hat_.resize(q_hat_max);	// set the terms involving the bounds first.

	if( num_act_change ) {
		size_type s = 0;
		for( size_type k = 1; k <= num_act_change; ++k ) {
			const int ij = ij_act_change[k-1];
			const QPSchurPack::EBounds bnd = bnds[k-1];
			if( ij < 0 ) {
				// Initially fixed variable being freed.
				ij_map_[s]		= ij;
				constr_norm_[s]	= 1.0;
				bnds_[s]		= QPSchurPack::FREE;
				d_hat_[s]		= - g(-ij);		// - g^X_{l^{(-)}}
				++s;
			}
			else {
				// Constraint being added to the active-set
				if( ij <= n ) {
					// Fixing a variable to a bound
					QPSchurPack::EBounds x_init_bnd = x_init(ij);
					if( x_init_bnd == QPSchurPack::FREE ) {
						// initially free variable being fixed
						ij_map_[s]		= ij;
						constr_norm_[s]	= 1.0;
						bnds_[s]		= bnd;
						d_hat_[s]		= constraints.get_bnd(ij,bnd);
						++s;
					}
					else {
						// Initially fixed variable being fixed to another bound
						// Free the variable first
						ij_map_[s]		= ij;
						constr_norm_[s]	= 1.0;
						bnds_[s]		= QPSchurPack::FREE;
						d_hat_[s]		= - g(ij);		// - g^X_{l^{(-)}}
						++s;
						// Now fix to a different bound
						ij_map_[s]		= ij;
						constr_norm_[s]	= 1.0;
						bnds_[s]		= bnd;
						d_hat_[s]		= constraints.get_bnd(ij,bnd) - b_X(l_x_X_map(ij));
						++s;
					}
				}
				else {
					// Adding a general inequality (or equality) constraint
					ij_map_[s]		= ij;
					constr_norm_[s]	= 1.0;	// ToDo: We need to compute this in an efficient way!
					bnds_[s]		= bnd;
					d_hat_[s]		= constraints.get_bnd(ij,bnd);	// \bar{b}_{j^{(+)}}
					++s;
				}
			}
		}
		assert(s == q_hat);
	}

	// Setup P_XF_hat_ and P_plus_hat_
	P_XF_hat_row_.resize(q_F_hat_max);
	P_XF_hat_col_.resize(q_F_hat_max);
	P_plus_hat_row_.resize(q_plus_hat_max);
	P_plus_hat_col_.resize(q_plus_hat_max);
	if(q_hat) {
		// See QPSchur.h for description of P_XF_hat and P_plus_hat
		size_type
			k_XF_hat = 0,	// zero based
			k_plus_hat = 0;	// zero based
		ij_map_t::const_iterator
			ij_itr 		= ij_map_.begin(),
			ij_itr_end	= ij_itr + q_hat;
		for( size_type s = 1; ij_itr != ij_itr_end; ++ij_itr, ++s ) {
			const int ij = *ij_itr;
			if( ij < 0 ) {
				const size_type i = -ij;
				assert( i <= n );
				// [P_XF_hat](:,s) = e(i)
				P_XF_hat_row_[k_XF_hat] = i;
				P_XF_hat_col_[k_XF_hat] = s;
				++k_XF_hat;
			}
			else if( !(ij <= n && x_init(ij) != QPSchurPack::FREE ) ) {
				const size_type j = ij;
				assert( 0 < j && j <= n + m_breve );
				// [P_plus_hat](:,s) = e(j)
				P_plus_hat_row_[k_plus_hat] = j;
				P_plus_hat_col_[k_plus_hat] = s;
				++k_plus_hat;
			}
		}
		assert( k_XF_hat == q_F_hat );
		assert( k_plus_hat == q_plus_hat );
	}
	P_XF_hat_.initialize_and_sort(
		  n,q_hat,q_F_hat,0,0,GPMSTP::BY_ROW
		, q_F_hat ? &P_XF_hat_row_[0] : NULL
		, q_F_hat ? &P_XF_hat_col_[0] : NULL
		,test
		);
	P_plus_hat_.initialize_and_sort(
		  n+m_breve,q_hat,q_plus_hat,0,0,GPMSTP::BY_ROW
		, q_plus_hat ? &P_plus_hat_row_[0] : NULL
		, q_plus_hat ? &P_plus_hat_col_[0] : NULL
		,test
		);

	// Setup Q_XD_hat_
	Q_XD_hat_row_.resize(q_D_hat_max);
	Q_XD_hat_col_.resize(q_D_hat_max);
	if(q_D_hat) {
		// See QPSchur.h for description of Q_XD_hat
		size_type
			k_XD_hat = 0;	// zero based
		GenPermMatrixSlice::const_iterator
			Q_X_itr = qp.Q_X().begin();	// This is sorted by row already!
		P_row_t::const_iterator
			XF_search 		= P_XF_hat_row_.begin(),	// These are already sorted by row!
			XF_search_end 	= XF_search + q_F_hat;
		for( size_type l = 1; l <= n_X; ++l, ++Q_X_itr ) {
			const size_type i = Q_X_itr->row_i();	// Already sorted by row
			// Look for i in XF
			for( ; XF_search != XF_search_end && *XF_search < i; ++XF_search ) ;
			if( XF_search == XF_search_end || (XF_search < XF_search_end && *XF_search > i) ) {
				// We went right past i and did not find it so
				// this variable has not been freed so lets add it!
				Q_XD_hat_row_[k_XD_hat] = i;
				Q_XD_hat_col_[k_XD_hat] = k_XD_hat + 1;
				++k_XD_hat;
			}
		}
		assert( k_XD_hat == q_D_hat );
	}
	Q_XD_hat_.initialize(
		  n,q_D_hat,q_D_hat,0,0,GPMSTP::BY_ROW	// Should already be sorted by row!
		, q_D_hat ? &Q_XD_hat_row_[0] : NULL
		, q_D_hat ? &Q_XD_hat_col_[0] : NULL
		,test
		);

	// Setup l_fxfx
	l_fxfx_.resize(q_D_hat_max);
	if(q_D_hat) {
		for( size_type k = 0; k < q_D_hat; ++k ) {
			l_fxfx_[k] = l_x_X_map(Q_XD_hat_row_[k]);
		}
	}

	// Set the rest of the terms in d_hat involving matrices
	//
	// d_hat += - P_XF_hat' * G * b_XX - P_plus_hat' * A_bar' * b_XX
	// 
	// where: b_XX = Q_X * b_X
	// 
	if( q_hat ) {
		SpVector b_XX;
		V_MtV( &b_XX, qp.Q_X(), BLAS_Cpp::no_trans, b_X );
		Vp_StPtMtV( &d_hat_(1,q_hat), -1.0, P_XF_hat_, BLAS_Cpp::trans
			, G, BLAS_Cpp::no_trans, b_XX() );
		Vp_StPtMtV( &d_hat_(1,q_hat), -1.0, P_plus_hat_, BLAS_Cpp::trans
			, constraints.A_bar(), BLAS_Cpp::trans, b_XX() );
	}

	// Setup U_hat
	U_hat_.initialize( &G, m ? &qp.A() : NULL, &constraints.A_bar()
		, &qp.Q_R(), &P_XF_hat_, &P_plus_hat_ );

	// Set the rest of the members
	test_		= test;
	qp_			= &qp;
	x_init_		= &x_init;
	n_			= n;
	n_R_		= n_R;
	m_			= m;
	m_breve_	= m_breve;
	q_plus_hat_	= q_plus_hat;
	q_F_hat_	= q_F_hat;
	q_C_hat_	= q_C_hat;

	// Resize storage for z_hat, p_z_hat, mu_D_hat, and p_mu_D_hat and set to zero by default
	z_hat_.resize(q_hat_max);
	p_z_hat_.resize(q_hat_max);
	mu_D_hat_.resize(n_X);
	p_mu_D_hat_.resize(n_X);

	initialized_ = true;	// set to true tenatively so that we can
							// print this stuff.

	if( out && (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
		*out
			<< "\nPrint definition of Active-Set before the Schur complement is formed...\n";
		dump_act_set_quantities( *this, *out, false );
	}

	// Initialize and factorize the schur complement
	if( q_hat ) {
		// Temporary storage for S (dense)
		GenMatrix S_store(q_hat,q_hat);
		sym_gms S( S_store, BLAS_Cpp::lower );
		// S = -1.0 * U_hat' * inv(Ko) * U_hat
		M_StMtInvMtM( &S, -1.0, U_hat_, BLAS_Cpp::trans, qp.Ko()
			, MatrixSymFactorized::DUMMY_ARG );
		// Now add parts of V_hat
		if( q_F_hat ) {
			// S += P_XF_hat' * G * P_XF_hat
			Mp_StPtMtP( &S, 1.0, MatrixSymWithOp::DUMMY_ARG, qp_->G(), P_XF_hat_, BLAS_Cpp::no_trans );
		}
		if( q_F_hat && q_plus_hat ) {
			// S += P_XF_hat' * A_bar * P_plus_hat + P_plus_hat' * A_bar' * P_XF_hat
			qp_->constraints().A_bar().syr2k(
				BLAS_Cpp::no_trans, 1.0
				,P_XF_hat_, BLAS_Cpp::no_trans
				,P_plus_hat_, BLAS_Cpp::no_trans
				,1.0, &S );
		}
		if( q_F_hat && q_C_hat ) {
			// S += P_F_tilde' * P_C_hat + P_C_hat' * P_F_tilde
			assert(0);	// ToDo: Implement this!
		}

		if( out && (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
			*out
				<< "\nIninitial Schur Complement before it is factorized:\n"
				<< "\nS_hat =\nLower triangular part (ignore nonzeros above diagonal)\n"
				<< S_store;
		}
		// Initialize and factorize the schur complement!
		schur_comp().update_interface().initialize(
			S, q_hat_max, true
			, MatrixSymAddDelUpdateable::Inertia( q_plus_hat + q_C_hat, 0, q_F_hat ) );
		// ToDo: Think about how to deal with the case where we may want to
		// selectively remove some rows/columns of S in order to
		// get a nonsingular schur complement.  This may be complicated though.
		if( out && (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
			*out
				<< "\nSchur Complement after factorization:\n"
				<< "\nS_hat =\n"
				 << schur_comp().op_interface();
		}
	}
	else {
		schur_comp().update_interface().set_uninitialized();
	}

	// Success, we are initialized!
	initialized_ = true;
	return;

	}	// end try
	catch(...) {
		initialized_ = false;
		throw;
	}
}

void QPSchur::ActiveSet::refactorize_schur_comp()
{
	// ToDo: Finish Me
	assert(0);
}

void QPSchur::ActiveSet::add_constraint(
	  size_type ja, QPSchurPack::EBounds bnd_ja
	, bool update_steps, bool force_refactorization )
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using LinAlgPack::dot;
	using LinAlgOpPack::V_StMtV;
	using SparseLinAlgPack::dot;
	using SparseLinAlgPack::Vp_StPtMtV;
	using SparseLinAlgPack::V_InvMtV;
	
	typedef SparseLinAlgPack::EtaVector eta_t;

	assert_initialized();

	const QPSchurPack::QP::Constraints
		&constraints = qp_->constraints();

	if( is_init_fixed(ja) && (*x_init_)(ja) == bnd_ja ) {
		//
		// This is a variable that was initially fixed, then freed and now
		// is being fixed back to its original bound.
		//
		assert(0);	// ToDo: implement
	}
	else {
		//
		// Expand the size of the schur complement to add the new constraint
		//
		 
		// Compute the terms for the update
		
		value_type			d_p = 0.0;
		const size_type		q_hat = this->q_hat();
		Vector				t_hat(q_hat);
		value_type			alpha_hat = 0.0;
		bool				changed_bounds = false;
				
		if( ja <= n_ && !is_init_fixed(ja) ) {
			//
			// Fix an initially free variable is being fixed
			// 
			// u_p = e(ja) <: R^(n_R+m)
			const eta_t u_p = eta_t(ja,n_R_+m_);
			// r = inv(Ko)*u_p
			Vector r;	// ToDo: Make this sparse!
			V_InvMtV( &r, qp_->Ko(), no_trans, u_p() );
			// t_hat = - U_hat' * r
			if(q_hat)
				V_StMtV( &t_hat(), -1.0, U_hat_, trans, r() );
			// alpha_hat = - u_p ' * r
			alpha_hat = - dot( u_p(), r() );
			// d_p = \bar{b}_{j^{(+)}}
			d_p = constraints.get_bnd( ja, bnd_ja );

			changed_bounds = false;
		}
		else if ( is_init_fixed(ja) ) {
			// An intially fixed variable was freed and
			// is now being fixed to the other bound.

			assert(0);	// ToDo: Finish this!

			changed_bounds = true;
		}
		else {	// ja > n
			//
			// Add an extra equality or inequality constraint.
			//
			// u_p = [ Q_R' * A_bar * e(ja) ] n_R
			//       [        0             ] m
			const eta_t e_ja = eta_t(ja,n_+m_breve_);
			const MatrixWithOp &A_bar = constraints.A_bar();
			Vector u_p( n_R_ + m_ );	// ToDo: make this sparse
			Vp_StPtMtV( &u_p(1,n_R_), 1.0, qp_->Q_R(), trans, A_bar, no_trans, e_ja(), 0.0 );
			if( m_ )
				u_p(n_R_+1,n_R_+m_) = 0.0;
			// r = inv(Ko) * u_p
			Vector r;	// ToDo: Make this sparse!
			V_InvMtV( &r, qp_->Ko(), no_trans, u_p() );
			if(q_hat) {
				// t_hat = v_p - U_hat' * r
				//    where: v_p = P_XF_hat' * A_bar * e(ja)
				V_StMtV( &t_hat(), -1.0, U_hat_, trans, r() );
				Vp_StPtMtV( &t_hat(), 1.0, P_XF_hat_, trans, A_bar, no_trans, e_ja() );
			}
			// alpha_hat = - u_p ' * r
			alpha_hat = - dot( u_p(), r() );
			// d_p = \bar{b}_{j^{(+)}} - b_X' * Q_X' * A_bar * e(ja)
			// 
			// d_p = \bar{b}_{j^{(+)}}
			d_p = constraints.get_bnd( ja, bnd_ja );
			if(n_ > n_R_) {
				// d_p += - b_X' * Q_X' * A_bar * e(ja)
				r.resize( n_ - n_R_ );	// reuse storage
				Vp_StPtMtV( &r(), 1.0, qp_->Q_X(), trans, A_bar, no_trans, e_ja(), 0.0 );			
				d_p += - dot( qp_->b_X(), r() );
			}
			
			changed_bounds = false;
		}
		
		// Update the schur complement (if nonsingular)
		try {
			if(q_hat) {
				schur_comp().update_interface().augment_update(
					&t_hat(), alpha_hat, force_refactorization
					,MatrixSymAddDelUpdateable::EIGEN_VAL_NEG  );
			}
			else {
				schur_comp().update_interface().initialize(
					alpha_hat, (n_-n_R_) + n_-m_ );
			}
		}
		catch(const MatrixSymAddDelUpdateable::SingularUpdateException& excpt) {
			throw LDConstraintException( std::string( "QPSchur::ActiveSet::add_constraint(...) : "
				"Error, constraint appears to be linearly dependent:\n" )
				+ std::string( excpt.what() ) );
		}
		catch(const MatrixSymAddDelUpdateable::WrongInertiaUpdateException& excpt) {
			throw BadUpdateException( std::string( "QPSchur::ActiveSet::add_constraint(...) : "
				"Error, The updated schur complement appears to have the wrong inertia:\n" )
				+ std::string( excpt.what() ) );
		}

		// Update the rest of the augmented KKT system
		if( changed_bounds )
			++q_C_hat_;
		else
			++q_plus_hat_;
		const size_type q_hat_new = q_F_hat_ + q_C_hat_ + q_plus_hat_;
		// Add ij_map(q_hat) = ja to ij_map(...)
		ij_map_[q_hat_new - 1]	= ja;
		// Add constr_norm(q_hat) to constr_norm(...)
		constr_norm_[q_hat_new - 1]	= 1.0;	// ToDo: Compute this for real!
		// Add bnds(q_hat)_ = bnd_ja to bnds(...)
		bnds_[q_hat_new - 1] 	= bnd_ja;
		// Augment d_hat = [ d_hat; d_p ]
		d_hat_(q_hat_new) = d_p;
		// Augment z_hat with new (zeroed) multiplier value, z_hat = [ z_hat; 0 ]
		z_hat_(q_hat_new) = 0.0;
		if( update_steps ) {
			// Set the step for this multiplier to 1.0, p_z_hat = [ p_z_hat; 1 ]
			p_z_hat_(q_hat_new) = 1.0;
		}
		if( !changed_bounds ) {
			// Insert (ja, q_hat_new) into P_plus_hat, sorted by row
			insert_pair_sorted(ja,q_hat_new,q_plus_hat_,&P_plus_hat_row_,&P_plus_hat_col_);
		}
	}
	// Update the permutation matrices and U_hat
	reinitialize_matrices(test_);
}

void QPSchur::ActiveSet::drop_constraint(
	 int jd , bool force_refactorization )
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using LinAlgPack::dot;
	using LinAlgOpPack::V_StMtV;
	using LinAlgOpPack::V_MtV;
	using LinAlgOpPack::Vp_MtV;
	using SparseLinAlgPack::dot;
	using SparseLinAlgPack::transVtMtV;
	using SparseLinAlgPack::Vp_StPtMtV;
	using SparseLinAlgPack::V_InvMtV;
	
	typedef SparseLinAlgPack::EtaVector eta_t;

	assert_initialized();

	if( jd < 0 ) {
		//
		// A variable initially fixed is being freed.
		// Increase the dimension of the augmented the KKT system!
		//
		size_type
			q_hat      = this->q_hat(),
			q_F_hat    = this->q_F_hat(),
			q_plus_hat = this->q_plus_hat(),
			q_D_hat    = this->q_D_hat();
		// Get indexes
		const size_type id = -jd;
		assert( 1 <= id && id <= n_ );
		const size_type ld = qp_->l_x_X_map()(-jd);
		assert( 1 <= ld && ld <= n_ - n_R_ );
		size_type kd; // Find kd
		{for( kd = 1; kd <= q_D_hat; ++kd ) {
			if( l_fxfx_[kd-1] == ld ) break;
		}}
		assert( kd < q_D_hat + 1 );
		// Get references
		const MatrixSymWithOp
			&G           = qp_->G();
		const VectorSlice
			g            = qp_->g();
		const MatrixWithOp
			&A_bar       = qp_->constraints().A_bar();
		const MatrixSymWithOpFactorized
			&Ko          = qp_->Ko();
		const MatrixWithOp
			&U_hat       = this->U_hat();
		const GenPermMatrixSlice
			&Q_R         = qp_->Q_R(),
			&Q_X         = qp_->Q_X(),
			&P_XF_hat    = this->P_XF_hat(),
			&P_plus_hat  = this->P_plus_hat();
		const VectorSlice
			b_X          = qp_->b_X();
		//
		// Compute the update quantities to augmented KKT system
		//
		// e_id
		eta_t e_id(id,n_);
		// u_p = [ Q_R'*G*e_id ; A'*e_id ] <: R^(n_R+m)
		Vector u_p(n_R_+m_);
		Vp_StPtMtV( &u_p(1,n_R_), 1.0, Q_R, trans, G, no_trans, e_id(), 0.0 );
		if( m_ )
			V_MtV( &u_p(n_R_+1,n_R_+m_), qp_->A(), trans, e_id() );
		const value_type
			nrm_u_p = LinAlgPack::norm_inf( u_p() );
		// sigma = e_id'*G*e_id <: R
		const value_type
			sigma = transVtMtV( e_id(), G, no_trans, e_id() );
		// d_p = - g(id) - b_X'*(Q_X'*G*e_id) <: R
		Vector Q_X_G_e_id(Q_X.cols());
		Vp_StPtMtV( &Q_X_G_e_id(), 1.0, Q_X, trans, G, no_trans, e_id(), 0.0 );
		const value_type
			d_p = -g(id) - dot( b_X, Q_X_G_e_id() );
		// r = inv(Ko)*u_p <: R^(n_R+m)
		Vector r;
		if( nrm_u_p > 0.0 )
			V_InvMtV( &r, Ko, no_trans, u_p() );
		// t_hat = v_p - U_hat'*r
		// where: v_p = P_XF_hat'*G*e_id + P_plus_hat'*A_bar'*e_id <: R^(q_hat)
		Vector t_hat(q_hat);
		if(q_hat) {
			t_hat = 0.0;
			// t_hat += v_p
			if(q_F_hat_)
				Vp_StPtMtV( &t_hat(), 1.0, P_XF_hat, trans, G, no_trans, e_id() );
			if(q_plus_hat_)
				Vp_StPtMtV( &t_hat(), 1.0, P_plus_hat, trans, A_bar, trans, e_id() );
			// t_hat += U_hat'*r
			if( nrm_u_p > 0.0 )
				Vp_MtV( &t_hat(), U_hat, trans, r() );
		}
		// alpha_hat = sigma - u_p'*r
		const value_type
			alpha_hat = sigma - ( nrm_u_p > 0.0 ? dot(u_p(),r()) : 0.0 );
		//
		// Update the schur complement (if nonsingular)
		//
		try {
			if(q_hat) {
				schur_comp().update_interface().augment_update(
					&t_hat(), alpha_hat, force_refactorization
					,MatrixSymAddDelUpdateable::EIGEN_VAL_POS );
			}
			else {
				schur_comp().update_interface().initialize(
					alpha_hat, (n_-n_R_) + n_-m_ );
			}
		}
		catch(const MatrixSymAddDelUpdateable::SingularUpdateException& excpt) {
			throw LDConstraintException( std::string( "QPSchur::ActiveSet::drop_constraint(...) : "
				"Error, the new KKT system with the freed variable appears to be singular:\n" )
				+ std::string( excpt.what() ) );
		}
		catch(const MatrixSymAddDelUpdateable::WrongInertiaUpdateException& excpt) {
			throw BadUpdateException( std::string( "QPSchur::ActiveSet::drop_constraint(...) : "
				"Error, The updated schur complement appears to have the wrong inertia:\n" )
				+ std::string( excpt.what() ) );
		}
		//
		// Remove multiplier from mu_D_hat(...)
		//
		// remove l_fxfx(kd) == ld from l_fxfx(...)
		std::copy( l_fxfx_.begin() + kd, l_fxfx_.begin() + q_D_hat
			, l_fxfx_.begin() + (kd-1) );
		// remove mu_D_hat(kd) from mu_D_hat(...)
		std::copy( mu_D_hat_.begin() + kd, mu_D_hat_.begin() + q_D_hat
			, mu_D_hat_.begin() + (kd-1) );
		// remove Q_XD_hat(:,kd) = e(id) from Q_XD_hat
		std::copy( Q_XD_hat_row_.begin() + kd, Q_XD_hat_row_.begin() + q_D_hat
			, Q_XD_hat_row_.begin() + (kd-1) );
		std::copy( Q_XD_hat_col_.begin() + kd, Q_XD_hat_col_.begin() + q_D_hat
			, Q_XD_hat_col_.begin() + (kd-1) );
		deincrement_indices( kd, &Q_XD_hat_col_, q_D_hat-1 );
		//
		// Update the counts
		//
		++q_F_hat_;
		q_hat = this->q_hat();
		//
		// Add the elements for newly freed variable
		//
		// add ij_map(q_hat) == -id to ij_map(...)
		ij_map_[q_hat-1] = -id;
		// add s_map(-id) == q_hat to s_map(...)
		// ToDo: implement s_map(...)
		// add bnds(q_hat) == FREE to bnds(...)
		bnds_[q_hat-1] = QPSchurPack::FREE;
		// add d_hat(q_hat) == d_p to d_hat(...)
		d_hat_[q_hat-1] = d_p;
		// add p_X(ld) == 0 to the end of z_hat(...)
		z_hat_[q_hat-1] = 0.0; // This is needed so that (z_hat + beta*t_D*p_z_hat)(q_hat) == 0
		// Insert (id,q_hat) into P_XF_hat sorted by row
		insert_pair_sorted(id,q_hat,q_F_hat_,&P_XF_hat_row_,&P_XF_hat_col_);
	}
	else {
		//
		// Shrink the dimension of the augmented KKT system to remove the constraint!
		//
		const size_type q_hat = this->q_hat();
		const size_type sd = s_map(jd);
		assert(sd);
		// Delete the sd row and column for S_hat
		schur_comp().update_interface().delete_update(
			sd,force_refactorization
			,MatrixSymAddDelUpdateable::EIGEN_VAL_NEG );
		// Remove the ij_map(s) = jd element from ij_map(...)
		std::copy( ij_map_.begin() + sd, ij_map_.begin() + q_hat
			, ij_map_.begin() + (sd-1) );
		// Remove the constr_norm(s) elment from constr_norm(...)
		std::copy( constr_norm_.begin() + sd, constr_norm_.begin() + q_hat
			, constr_norm_.begin() + (sd-1) );
		// Remove the bnds(s) element from bnds(...)
		std::copy( bnds_.begin() + sd, bnds_.begin() + q_hat
			, bnds_.begin() + (sd-1) );
		// Remove the d_hat(s) element from d_hat(...)
		std::copy( d_hat_.begin() + sd, d_hat_.begin() + q_hat
			, d_hat_.begin() + (sd-1) );
		// Remove the z_hat(s) element from z_hat(...)
		std::copy( z_hat_.begin() + sd, z_hat_.begin() + q_hat
			, z_hat_.begin() + (sd-1) );
		// Remove the p_z_hat(s) element from p_z_hat(...)
		std::copy( p_z_hat_.begin() + sd, p_z_hat_.begin() + q_hat
			, p_z_hat_.begin() + (sd-1) );
		if( is_init_fixed( jd ) ) {
			// This must be an intially fixed variable, currently fixed at a different
			// bound.  In this case nothing else has to be modifed.
			--q_C_hat_;
		}
		else {
			// We must remove (jd,sd) from P_plus_hat
			P_row_t::iterator
				itr = std::lower_bound( P_plus_hat_row_.begin(), P_plus_hat_row_.begin()+q_plus_hat_, jd );
			assert( itr != P_plus_hat_row_.end() );
			const size_type p = itr - P_plus_hat_row_.begin();
			std::copy( P_plus_hat_row_.begin() + p + 1, P_plus_hat_row_.begin()+q_plus_hat_,
				P_plus_hat_row_.begin() + p );
			std::copy( P_plus_hat_col_.begin() + p + 1, P_plus_hat_col_.begin()+q_plus_hat_,
				P_plus_hat_col_.begin() + p );
			--q_plus_hat_;
		}
		// Deincrement all counters in permutation matrices for removed element
		deincrement_indices( sd, &P_XF_hat_col_, q_F_hat_ );
		deincrement_indices( sd, &P_plus_hat_col_, q_plus_hat_ );
	}
	// Update the permutation matrices and U_hat
	reinitialize_matrices(test_);
}

void QPSchur::ActiveSet::drop_add_constraints(
	int jd, size_type ja, QPSchurPack::EBounds bnd_ja, bool update_steps )
{
	drop_constraint( jd, false );
	add_constraint( ja, bnd_ja, update_steps, true );
}

QPSchur::ActiveSet::QP&
QPSchur::ActiveSet::qp()
{
	assert_initialized();
	return *qp_;
}

const QPSchur::ActiveSet::QP&
QPSchur::ActiveSet::qp() const
{
	assert_initialized();
	return *qp_;
}

size_type QPSchur::ActiveSet::q_hat() const
{
	assert_initialized();
	return q_plus_hat_ + q_F_hat_ + q_C_hat_;
}

size_type QPSchur::ActiveSet::q_plus_hat() const
{
	assert_initialized();
	return q_plus_hat_;
}

size_type QPSchur::ActiveSet::q_F_hat() const
{
	assert_initialized();
	return q_F_hat_;
}

size_type QPSchur::ActiveSet::q_C_hat() const
{
	assert_initialized();
	return q_C_hat_;
}

size_type QPSchur::ActiveSet::q_D_hat() const
{
	assert_initialized();
	return (n_ - n_R_) - q_F_hat_;  // n^{X} - \hat{q}^{F}
}

int QPSchur::ActiveSet::ij_map( size_type s ) const
{
	assert( 1 <= s && s <= this->q_hat() );
	return ij_map_[s-1];
}

size_type QPSchur::ActiveSet::s_map( int ij ) const
{
	ij_map_t::const_iterator
		begin	= ij_map_.begin(),
		end		= begin + q_hat(),
		itr = std::find( begin, end, ij );
	return ( itr != end ? (itr - begin) + 1 : 0 );
}

value_type QPSchur::ActiveSet::constr_norm( size_type s ) const
{
	assert( 1 <= s && s <= this->q_hat() );
	return constr_norm_(s);
}

QPSchurPack::EBounds QPSchur::ActiveSet::bnd( size_type s ) const
{
	assert( 1 <= s && s <= this->q_hat() );
	return bnds_[s-1];
}

size_type QPSchur::ActiveSet::l_fxfx( size_type k ) const
{
	assert( 1 <= k && k <= this->q_D_hat() );
	return l_fxfx_[k-1];
}

const QPSchur::U_hat_t& QPSchur::ActiveSet::U_hat() const
{
	assert_initialized();
	return U_hat_;
}

const MatrixSymWithOpFactorized& QPSchur::ActiveSet::S_hat() const
{
	assert_initialized();
	return schur_comp().op_interface();
}

const GenPermMatrixSlice& QPSchur::ActiveSet::P_XF_hat() const
{
	assert_initialized();
	return P_XF_hat_;
}

const GenPermMatrixSlice& QPSchur::ActiveSet::P_plus_hat() const
{
	assert_initialized();
	return P_plus_hat_;
}

const GenPermMatrixSlice& QPSchur::ActiveSet::Q_XD_hat() const
{
	assert_initialized();
	return Q_XD_hat_;
}

const VectorSlice QPSchur::ActiveSet::d_hat() const
{
	assert_initialized();
	return d_hat_(1,q_hat());
}

VectorSlice QPSchur::ActiveSet::z_hat()
{
	assert_initialized();
	return z_hat_(1,q_hat());
}

const VectorSlice QPSchur::ActiveSet::z_hat() const
{
	assert_initialized();
	return z_hat_(1,q_hat());
}

VectorSlice QPSchur::ActiveSet::p_z_hat()
{
	assert_initialized();
	return p_z_hat_(1,q_hat());
}

const VectorSlice QPSchur::ActiveSet::p_z_hat() const
{
	assert_initialized();
	return p_z_hat_(1,q_hat());
}

VectorSlice QPSchur::ActiveSet::mu_D_hat()
{
	assert_initialized();
	return mu_D_hat_(1,q_D_hat());
}

const VectorSlice QPSchur::ActiveSet::mu_D_hat() const
{
	assert_initialized();
	return mu_D_hat_(1,q_D_hat());
}

VectorSlice QPSchur::ActiveSet::p_mu_D_hat()
{
	assert_initialized();
	return p_mu_D_hat_(1,q_D_hat());
}

const VectorSlice QPSchur::ActiveSet::p_mu_D_hat() const
{
	assert_initialized();
	return p_mu_D_hat_(1,q_D_hat());
}

bool QPSchur::ActiveSet::is_init_fixed( size_type j ) const
{
	assert_initialized();
	return j <= n_ && (*x_init_)(j) != QPSchurPack::FREE;
}

bool QPSchur::ActiveSet::all_dof_used_up() const
{
	return n_ - m_ == (n_ - n_R_) - q_F_hat_ + q_C_hat_ + q_plus_hat_;
}

// private member functions for QPSchur::ActiveSet

void QPSchur::ActiveSet::assert_initialized() const
{
	if( !initialized_ )
		throw std::logic_error( "QPSchur::ActiveSet::assert_initialized() : Error, "
			"The active set has not been initialized yet!" );
}

void QPSchur::ActiveSet::assert_s( size_type s) const
{
	assert( s <= q_hat() );	// ToDo: Throw an exception
}

void QPSchur::ActiveSet::reinitialize_matrices(bool test)
{
	namespace GPMSTP = SparseLinAlgPack::GenPermMatrixSliceIteratorPack;

	const size_type q_hat = this->q_hat();
	const size_type q_D_hat = this->q_D_hat();

	P_XF_hat_.initialize(
		  n_,q_hat,q_F_hat_,0,0,GPMSTP::BY_ROW
		, q_F_hat_ ? &P_XF_hat_row_[0] : NULL
		, q_F_hat_ ? &P_XF_hat_col_[0] : NULL
		,test
		);
	P_plus_hat_.initialize(
		  n_+m_breve_,q_hat,q_plus_hat_,0,0,GPMSTP::BY_ROW
		, q_plus_hat_ ? &P_plus_hat_row_[0] : NULL
		, q_plus_hat_ ? &P_plus_hat_col_[0] : NULL
		,test
		);
	
	Q_XD_hat_.initialize(
		  n_,q_D_hat,q_D_hat,0,0,GPMSTP::BY_ROW
		, q_D_hat ? &Q_XD_hat_row_[0] : NULL
		, q_D_hat ? &Q_XD_hat_col_[0] : NULL
		,test
		);
	U_hat_.initialize( &qp_->G(), m_ ? &qp_->A() : NULL, &qp_->constraints().A_bar()
		, &qp_->Q_R(), &P_XF_hat_, &P_plus_hat_);
}

// public member functions for QPSchur

value_type QPSchur::DEGENERATE_MULT = std::numeric_limits<value_type>::min();

QPSchur::QPSchur(
		  const schur_comp_ptr_t& 	schur_comp
		, size_type					max_iter
		, value_type				feas_tol
		, value_type				loose_feas_tol
		, value_type				dual_infeas_tol
		, value_type				huge_primal_step
		, value_type				huge_dual_step
		)
	:
		schur_comp_(schur_comp)
		,max_iter_(max_iter)
		,feas_tol_(feas_tol)
		,loose_feas_tol_(loose_feas_tol)
		,dual_infeas_tol_(dual_infeas_tol)
		,huge_primal_step_(huge_primal_step)
		,huge_dual_step_(huge_dual_step)
		,act_set_(schur_comp)
{}

QPSchur::ESolveReturn QPSchur::solve_qp(
	  QP& qp
	, size_type num_act_change, const int ij_act_change[]
		, const QPSchurPack::EBounds bnds[]
	, std::ostream *out, EOutputLevel output_level, ERunTests test_what
	, VectorSlice* x, SpVector* mu, VectorSlice* lambda, SpVector* lambda_breve
	, size_type* iter, size_type* num_adds, size_type* num_drops
	)
{
	using std::setw;
	using std::endl;
	using std::right;
	using LinAlgPack::norm_inf;
	using SparseLinAlgPack::norm_inf;
	using SparseLinAlgPack::V_InvMtV;

	if( !out )
		output_level = NO_OUTPUT;

	const int	dbl_prec = 6;
	int			prec_saved = out ? out->precision() : 0;

	ESolveReturn
		solve_return = SUBOPTIMAL_POINT;

	try {

	// Print QPSchur output header
	if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
		*out
			<< "\n*** Entering QPSchur::solve_qp(...) ***\n";
	}

	// Print the definition of the QP to be solved.
	if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
		*out
			<< "\n*** Dump the definition of the QP to be solved ***\n";
		qp.dump_qp(*out);
	}
	
	// Print warm start info.
	if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
		*out
			<< "\n*** Warm start info\n"
			<< "\nnum_act_change = " << num_act_change << endl;
	}

	if( num_act_change > 0 && (int)output_level >= (int)OUTPUT_ACT_SET ) {
		*out << std::setprecision(dbl_prec);
		*out
			<< endl
			<< right << setw(5) << "s"
			<< right << setw(20) << "ij_act_change"
			<< right << setw(10) << "bnds" << endl
			<< right << setw(5) << "---"
			<< right << setw(20) << "-------------"
			<< right << setw(10) << "--------" << endl;
		for( size_type s = 1; s <= num_act_change; ++s )
			*out
				<< right << setw(5) << s
				<< right << setw(20) << ij_act_change[s-1]
				<< right << setw(10) << bnd_str(bnds[s-1]) << endl;
		*out << std::setprecision(prec_saved);
	}

	// Initialize the active set.
	try {
		act_set_.initialize( qp, num_act_change, ij_act_change, bnds
			, test_what == RUN_TESTS, out, output_level );
	}
	catch( const ActiveSet::LDConstraintException& excpt ) {
		if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
			*out
				<< "\n*** Error in initializing schur complement\n"
				<< excpt.what() << std::endl
				<< "\nSetting num_act_change = 0 and proceeding with a cold start...\n";
		}
		act_set_.initialize( qp, num_act_change = 0, ij_act_change, bnds
			, test_what == RUN_TESTS, out, output_level );
	}

	// Compute vo =  inv(Ko) * fo
	V_InvMtV( &vo_, qp.Ko(), BLAS_Cpp::no_trans, qp.fo() );

	if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
		*out
			<< "\nSolution to the initial KKT system, vo = inv(Ko)*fo:\n\n||vo||inf = " << norm_inf(vo_()) << std::endl;
	}
	if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
		*out
			<< "\nvo =\n" << vo_();
	}

	// ////////////////////////////////////////////////
	// Remove constraints until we are dual feasible.
	// 
	// Here we are essentially performing a primal algorithm where we are only
	// removing constraints.  If the dual variables are not dual feasible then
	// we will remove the one with the largest scaled dual feasibility
	// violation then compute the dual variables again.  Eventually we
	// will come to a point where we have a dual feasible point.  If
	// we have to, we will remove all of the inequality constraints and then
	// this will by default be a dual feasible point (i.e. we picked all the
	// wrong inequality constraints).
	// 
	// The difficulty here is in dealing with near degenerate constraints.
	// If a constraint is near degenerate then we would like to not drop
	// it since we may have to add it again later.
	// 
	*iter = 0;
	*num_adds = 0;
	*num_drops = 0;
	// Print header for removing constraints
	if( (int)output_level >= (int)OUTPUT_ITER_SUMMARY && num_act_change > 0 ) {
		*out
			<< "\n***"
			<< "\n*** Removing constriants until we are dual feasible"
			<< "\n***\n"
			<< "\n*** Start by removing constraints within the Schur complement first\n";
	}
	// Print summary header for max_viol and jd.
	if( (int)OUTPUT_ITER_SUMMARY <= (int)output_level 
		&& (int)output_level < (int)OUTPUT_ITER_QUANTITIES
		&& num_act_change > 0  )
	{
		*out << std::setprecision(dbl_prec);
		*out     
			<< "\nIf max_viol > 0 and jd != 0 then constraint jd will be dropped from the active set\n\n"
			<< setw(20)	<< "max_viol"
			<< setw(5)	<< "sd"
			<< setw(5)	<< "jd"		<< endl
			<< setw(20)	<< "--------------"
			<< setw(5)	<< "----"
			<< setw(5)	<< "----"	<< endl;
	}
	for( int k = num_act_change; k > 0; --k, ++(*iter) ) {
		// Compute z_hat (z_hat = inv(S_hat)*(d_hat - U_hat'*vo))
		VectorSlice z_hat = act_set_.z_hat();
		calc_z( act_set_.S_hat(), act_set_.d_hat(), act_set_.U_hat(), vo_()
			, &z_hat );
		// Determine if we are dual feasible.
		value_type	max_viol = 0.0;	// max scaled violation of dual feasability.
		size_type	jd = 0;			// indice of constraint with max scaled violation.
		VectorSlice::const_iterator
			z_itr = const_cast<const VectorSlice&>(z_hat).begin();
		const size_type q_hat = act_set_.q_hat();	// Size of schur complement system.
		// Print header for s, z_hat(s), bnd(s), viol, max_viol and jd
		if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
			*out
				<< "\nLooking for a constraint with the maximum dual infeasiblity to drop...\n\n"
				<< right << setw(5)		<< "s"
				<< right << setw(20)	<< "z_hat"
				<< right << setw(20)	<< "bnd"
				<< right << setw(20)	<< "viol"
				<< right << setw(20)	<< "max_viol"
				<< right << setw(5)		<< "jd"	<< endl
				<< right << setw(5)		<< "----"
				<< right << setw(20)	<< "--------------"
				<< right << setw(20)	<< "--------------"
				<< right << setw(20)	<< "--------------"
				<< right << setw(20)	<< "--------------"
				<< right << setw(5)		<< "----"	<< endl;
		}
		for( int s = 1; s <= q_hat; ++s, ++z_itr ) {
			int j = act_set_.ij_map(s);
			if( j > 0 ) {
				// This is for an active constraint not initially in Ko so
				// z_hat(s) = mu(j)
				QPSchurPack::EBounds bnd = act_set_.bnd(s);
				const value_type inf = std::numeric_limits<value_type>::max();
				value_type viol = -inf;
				if( bnd == QPSchurPack::LOWER ) {
					viol = (*z_itr) / act_set_.constr_norm(s);	// + scaled violation.
					if( viol > 0 ) {
						if( viol < dual_infeas_tol() ) {
							// We need to fix the sign of this near degenerate multiplier
							assert(0);
						}
						else if( viol > max_viol ) {
							max_viol = viol;
							jd = j;
						}
					}
				}
				else if( bnd == QPSchurPack::UPPER ) {
					viol = -(*z_itr) / act_set_.constr_norm(s);	// + scaled violation.
					if( viol > 0 ) {
						if( viol < dual_infeas_tol() ) {
							// We need to fix the sign of this near degenerate multiplier
							assert(0);
						}
						else if( viol > max_viol ) {
							max_viol = viol;
							jd = j;
						}
					}
				}
				// Print row for s, z_hat(s), bnd(s), viol, max_viol and jd
				if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
					*out << std::setprecision(dbl_prec);
					*out
						<< right << setw(5)		<< s
						<< right << setw(20)	<< *z_itr
						<< right << setw(20)	<< bnd_str(bnd)
						<< right << setw(20)	<< viol
						<< right << setw(20)	<< max_viol
						<< right << setw(5)		<< jd	<< endl;
				}
			}
		}
		// Print row of max_viol and jd
		if( (int)OUTPUT_ITER_SUMMARY <= (int)output_level 
			&& (int)output_level < (int)OUTPUT_ITER_QUANTITIES )
		{
			*out
				<< setw(20)	<< max_viol
				<< setw(5)	<< act_set_.s_map(jd)
				<< setw(5)	<< jd					<< endl;
		}
		if( jd == 0 ) break;	// We have a dual feasible point w.r.t. these constraints
		// Remove the active constraint with the largest scaled violation.
		act_set_.drop_constraint( jd );
		++(*iter);
		++(*num_drops);
		// Print U_hat, S_hat and d_hat.
		if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
			*out << std::setprecision(prec_saved);
			*out
				<< "\nPrinting active set after dropping constraint jd = " << jd << " ...\n";
			dump_act_set_quantities( act_set_, *out );
		}
	}

	if(out) *out << std::setprecision(prec_saved);

	// Print how many constraints where removed from the schur complement
	if( num_act_change > 0 && (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
		*out
			<< "\nThere where " << (*num_drops)
			<< " constraints dropped from the schur complement from the initial guess of the active set.\n";
	}

	// Compute v
	if( act_set_.q_hat() > 0 ) {
		v_.resize( qp.n_R() + qp.m() );
		calc_v( qp.Ko(), qp.fo(), act_set_.U_hat(), act_set_.z_hat(), &v_() );
		if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
			*out
				<< "\nSolution to the system; v = inv(Ko)*(fo - U_hat*z_hat):\n"
				<< "\n||v||inf = " << norm_inf(v_()) << std::endl;
		}
		if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
			*out
				<< "\nv =\n" << v_();
		}
	}
	else {
		v_ = vo_;
	}

	// Set x
	set_x( act_set_, v_(), x );
	if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
		*out
			<< "\nCurrent guess for unknowns x:\n\n||x||inf = " << norm_inf(*x) << std::endl;
	}
	if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
		*out
			<< "\nx =\n" << *x;
	}

	//
	// Determine if any initially fixed variables need to be freed by checking mu_D_hat.
	//
	if( act_set_.q_D_hat() ) {
		if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
			*out << "\n*** Second, free initially fixed variables not in Ko\n\n";
		}
		const QPSchurPack::QP::i_x_X_map_t&  i_x_X_map = act_set_.qp().i_x_X_map();
		const QPSchurPack::QP::x_init_t&     x_init    = act_set_.qp().x_init();
		// Print summary header for max_viol and jd.
		if( (int)OUTPUT_ITER_SUMMARY <= (int)output_level 
			&& (int)output_level < (int)OUTPUT_ITER_QUANTITIES
			&& num_act_change > 0  )
		{
			*out << std::setprecision(dbl_prec);
			*out     
				<< "\nIf max_viol > 0 and id != 0 then the variable x(id) will be freed from its initial bound\n\n"
				<< setw(20)	<< "max_viol"
				<< setw(5)	<< "kd"
				<< setw(5)	<< "id"		<< endl
				<< setw(20)	<< "--------------"
				<< setw(5)	<< "----"
				<< setw(5)	<< "----"	<< endl;
		}
		size_type q_D_hat = act_set_.q_D_hat(); // This will be deincremented
		while( q_D_hat > 0 ) {
			// mu_D_hat = ???
			VectorSlice mu_D_hat = act_set_.mu_D_hat();
			calc_mu_D( act_set_, *x, v_(), &mu_D_hat );
			// Determine if we are dual feasible.
			value_type	max_viol = 0.0;	// max scaled violation of dual feasability.
			int			id = 0;			// indice of variable with max scaled violation.
			size_type   kd = 0;
			VectorSlice::const_iterator
				mu_D_itr = const_cast<const VectorSlice&>(mu_D_hat).begin();
			// Print header for k, mu_D_hat(k), bnd, viol, max_viol and id
			if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
				*out
					<< "\nLooking for a variable bound with the max dual infeasibility to drop...\n\n"
					<< right << setw(5)		<< "k"
					<< right << setw(20)	<< "mu_D_hat"
					<< right << setw(20)	<< "bnd"
					<< right << setw(20)	<< "viol"
					<< right << setw(20)	<< "max_viol"
					<< right << setw(5)		<< "id"	<< endl
					<< right << setw(5)		<< "----"
					<< right << setw(20)	<< "--------------"
					<< right << setw(20)	<< "--------------"
					<< right << setw(20)	<< "--------------"
					<< right << setw(20)	<< "--------------"
					<< right << setw(5)		<< "----"	<< endl;
			}
			for( int k = 1; k <= q_D_hat; ++k, ++mu_D_itr ) {
				int
					i = i_x_X_map(act_set_.l_fxfx(k));
				QPSchurPack::EBounds
					bnd = x_init(i);
				const value_type inf = std::numeric_limits<value_type>::max();
				value_type viol = -inf;
				if( bnd == QPSchurPack::LOWER ) {
					viol = (*mu_D_itr);
					if( viol > 0 ) {
						if( viol < dual_infeas_tol() ) {
							// We need to fix the sign of this near degenerate multiplier
							assert(0);
						}
						else if( viol > max_viol ) {
							max_viol = viol;
							kd = k;
							id = i;
						}
					}
				}
				else if( bnd == QPSchurPack::UPPER ) {
					viol = -(*mu_D_itr);
					if( viol > 0 ) {
						if( viol < dual_infeas_tol() ) {
							// We need to fix the sign of this near degenerate multiplier
							assert(0);
						}
						else if( viol > max_viol ) {
							max_viol = viol;
							kd = k;
							id = i;
						}
					}
				}
				// Print row for k, mu_D_hat(k), bnd, viol, max_viol and jd
				if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
					*out << std::setprecision(dbl_prec);
					*out
						<< right << setw(5)		<< k
						<< right << setw(20)	<< *mu_D_itr
						<< right << setw(20)	<< bnd_str(bnd)
						<< right << setw(20)	<< viol
						<< right << setw(20)	<< max_viol
						<< right << setw(5)		<< id	<< endl;
				}
			}
			// Print row of max_viol and id
			if( (int)OUTPUT_ITER_SUMMARY <= (int)output_level 
				&& (int)output_level < (int)OUTPUT_ITER_QUANTITIES )
			{
				*out
					<< setw(20)	<< max_viol
					<< setw(5)	<< kd
					<< setw(5)	<< id         << endl;
			}
			if( id == 0 ) break;	// We have a dual feasible point w.r.t. these variable bounds
			// Remove the active constraint with the largest scaled violation.
			act_set_.drop_constraint( -id );
			++(*iter);
			++(*num_adds);
			--q_D_hat;
			// Print U_hat, S_hat and d_hat.
			if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
				*out << std::setprecision(prec_saved);
				*out
					<< "\nPrinting active set after freeing initially fixed variable id = " << id << " ...\n";
				dump_act_set_quantities( act_set_, *out );
			}
			if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
				*out
					<< "\nSolution to the new KKT system; z_hat = inv(S_hat)*(d_hat - U_hat'*vo), v = inv(Ko)*(fo - U_hat*z_hat):\n";
			}
			// Compute z_hat (z_hat = inv(S_hat)*(d_hat - U_hat'*vo))
			calc_z( act_set_.S_hat(), act_set_.d_hat(), act_set_.U_hat(), vo_()
					, &act_set_.z_hat() );
			if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
				*out
					<< "\n||z_hat||inf = " << norm_inf(act_set_.z_hat()) << std::endl;
			}
			if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
				*out
					<< "\nz_hat =\n" << act_set_.z_hat();
			}
			// Compute v
			calc_v( qp.Ko(), qp.fo(), act_set_.U_hat(), act_set_.z_hat(), &v_() );
			if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
				*out
					<< "\n||v||inf = " << norm_inf(v_()) << std::endl;
			}
			if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
				*out
					<< "\nv =\n" << v_();
			}
			// Set x
			set_x( act_set_, v_(), x );
			if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
				*out
					<< "\nCurrent guess for unknowns x:\n\n||x||inf = " << norm_inf(*x) << std::endl;
			}
			if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
				*out
					<< "\nx =\n" << *x;
			}
		}
	}

	// Print how many initially fixed variables where freed
	if( *num_adds > 0 && (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
		*out
			<< "\nThere where " << (*num_adds)
			<< " initially fixed variables not in Ko that where freed and added to the schur complement.\n";
	}

	// Run the primal dual algorithm
	solve_return = qp_algo(
		  PICK_VIOLATED_CONSTRAINT
		, out, output_level, test_what
		, vo_(), &act_set_, &v_()
		, x, iter, num_adds, num_drops
		);

	if(out) *out << std::setprecision(prec_saved);

	if( solve_return != OPTIMAL_SOLUTION )
		set_x( act_set_, v_(), x );

	set_multipliers( act_set_, v_(), mu, lambda, lambda_breve );

	// Print Solution x, lambda and mu
	if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
		switch(solve_return) {
			case OPTIMAL_SOLUTION:
				*out
					<< "\n*** Solution found!\n";
				break;
			case MAX_ITER_EXCEEDED:
				*out
					<< "\n*** Maximum iterations exceeded!\n";
				break;
			case MAX_ALLOWED_STORAGE_EXCEEDED:
				*out
					<< "\n*** The maxinum size of the schur complement has been exceeded!\n";
				break;
			case INFEASIBLE_CONSTRAINTS:
				*out
					<< "\n*** The constraints are infeasible!\n";
				break;
			case DUAL_INFEASIBILITY:
				*out
					<< "\n*** The dual variables are infeasible (numerical roundoff?)!\n";
				break;
			case SUBOPTIMAL_POINT:
				*out
					<< "\n*** The current point is suboptimal but we will return it anyway!\n";
				break;
			default:
				assert(0);
		}
	}
	if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
		*out	<< "\n||x||inf                = "	<< norm_inf(*x);
		*out	<< "\nmax(|mu(i)|)            = " 	<< norm_inf((*mu)());
		*out	<< "\nmin(|mu(i)|)            = " 	<< min_abs((*mu)());
		if(lambda)
			*out
				<< "\nmax(|lambda(i)|)        = "	<< norm_inf(*lambda)
				<< "\nmin(|lambda(i)|)        = "	<< min_abs(*lambda);
		if(lambda_breve)
			*out
				<< "\nmax(|lambda_breve(i)|)  = "	<< norm_inf((*lambda_breve)())
				<< "\nmin(|lambda_breve(i)|)  = "	<< min_abs((*lambda_breve)());
		*out << std::endl;
	}
	if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
		*out	<< "\nx =\n" 				<< *x;
		*out	<< "\nmu =\n" 				<< *mu;
		if(lambda)
			*out << "\nlambda =\n"			<< *lambda;
		if(lambda_breve)
			*out << "\nlambda_breve =\n"	<< *lambda_breve;
	}
	// Print 'goodby' header.
	if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
		*out
			<< "\n\n*** Leaving QPSchur::solve_qp(...) ***\n";
	}

	}	// end try
	catch(...) {
		if(out) *out << std::setprecision(prec_saved);
		throw;
	}

	return solve_return;
}

const QPSchur::ActiveSet& QPSchur::act_set() const
{
	return act_set_;
}

// protected member functions for QPSchur

QPSchur::ESolveReturn QPSchur::qp_algo(
	  EPDSteps next_step
	, std::ostream *out, EOutputLevel output_level, ERunTests test_what
	, const VectorSlice& vo, ActiveSet* act_set, VectorSlice* v
	, VectorSlice* x, size_type* iter, size_type* num_adds, size_type* num_drops
	)
{
	// ToDo:
	//	* Generalize for initially fixed variables left out of Ko that are later freed.

	using std::setw;
	using std::endl;
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using LinAlgPack::dot;
	using LinAlgPack::norm_inf;
	using LinAlgPack::Vt_S;
	using LinAlgPack::V_mV;
	using LinAlgPack::Vp_StV;
	using LinAlgPack::V_VmV;
	using LinAlgOpPack::Vp_V;
	using LinAlgOpPack::V_StV;
	using LinAlgOpPack::V_StMtV;
	using SparseLinAlgPack::dot;
	using SparseLinAlgPack::norm_inf;
	using SparseLinAlgPack::EtaVector;
	using SparseLinAlgPack::V_InvMtV;
	using SparseLinAlgPack::Vp_StMtV;
	using SparseLinAlgPack::Vp_StPtMtV;
	using LinAlgOpPack::V_MtV;

	// Print header for "Starting Primal-Dual Iterations"
	if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
		*out
			<< "\n*** Starting Primal-Dual Iterations ***\n";
	}

	const int	dbl_prec = 6;
	int			prec_saved = out ? out->precision() : 0;

	try {

	QPSchurPack::QP
		&qp = act_set->qp();
	const size_type
		n		= qp.n(),
		n_R		= qp.n_R(),
		m		= qp.m(),
		m_breve	= qp.constraints().m_breve();

	Vector
		v_plus(v->size()),
		z_hat_plus,
		p_v(v->size());

	const value_type
		inf = std::numeric_limits<value_type>::max();
	size_type itr;	// move to somewhere else?

	// Put these here because they need to be remembered between iterations if a linearly
	// dependent constriant is dropped.
	size_type				ja = 0;		// + indice of violated constraint to add to active set
	value_type				con_ja_val;	// value of violated constraint.
	value_type				b_a; // value of the violated bound
	value_type				norm_2_constr;	// norm of violated constraint
	QPSchurPack::EBounds	bnd_ja;	// bound of constraint ja which is violated.
	bool					can_ignore_ja;	// true if we can ignore a constraint if it is LD.
	bool					assume_lin_dep_ja;
	value_type				gamma_plus;	// used to store the new multipler value for the added
										// constraint.
	const int				summary_lines_counter_max = 15;
	int						summary_lines_counter = 0;
	int						jd = 0;	// + indice of constraint to delete from active set.
									// - indice of intially fixed variable to be freed
	int						last_jd = 0; // Last jd change to the active set
	value_type				t_P;	// Primal step length (constraint ja made active)
	value_type				t_D;	// Dual step length ( longest step without violating dual
									// feasibility of currently active constraints ).
	value_type				beta;	// +1 if multiplier of constraint being added is positive
									// -1 if multiplier of constraint being added is negative.
	bool					warned_degeneracy = false; // Will warn the user if degeneracy
									// is detected.
	value_type				dual_infeas_scale = 1.0;	// Scaling for determining if a
									// Lagrange multiplier is near degenerate.
	bool					return_to_init_fixed = false;	// True if the constraint being added
									// to the active set is a varible returning to its orginally
	                                // fixed variable bound.
	bool                    all_dof_used_up; // Used to keep track of when we don't need to compute p_v

	for( itr = 0; itr <= max_iter_; ++itr, ++(*iter) ) {

		// Print header for itr, nact, change (ADD, DROP), indice (ja, jb)
		//, bound (LOWER, UPPER, EQUALITY), violation (primal or dual), rank (LD,LI)
		if( (int)output_level == (int)OUTPUT_ITER_SUMMARY ) {
			if(out) *out << std::setprecision(dbl_prec);
			if( summary_lines_counter == 0 ) {
				summary_lines_counter = 1;
				*out
					<< endl
					<< setw(6)	<< "itr"
					<< setw(6)	<< "qhat"
					<< setw(8)	<< "change"
					<< setw(6)	<< "j"
					<< setw(10)	<< "bnd"
					<< setw(20)	<< "viol, p_z(jd)"
					<< setw(6)	<< "rank" << endl
					<< setw(6)	<< "----"
					<< setw(6)	<< "----"
					<< setw(8)	<< "------"
					<< setw(6)	<< "----"
					<< setw(10)	<< "--------"
					<< setw(20) << "--------------"
					<< setw(6)	<< "----"	<< endl;
			}
		}
		if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
			*out
				<< "\n************************************************"
				<< "\n*** qp_iter = " << itr
				<< "\n*** q_hat   = " << act_set->q_hat() << std::endl;
		}
		switch( next_step ) {	// no break; statements in this switch statement.
			case PICK_VIOLATED_CONSTRAINT: {
				if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
					*out
						<< "\n*** PICK_VIOLATED_CONSTRAINT\n";
				}
				// Set parts of x that are not currently fixed and may have changed.
				// Also, we want set specifially set those variables that where
				// initially free and then latter fixed to their bounds so that
				// they will not be seen as violated.
				set_x( *act_set, *v, x );
				if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
					*out
						<< "\n||x||inf = " << norm_inf(*x) << std::endl;
				}
				if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
					*out
						<< "\nx =\n" << *x;
				}
				act_set->qp().constraints().pick_violated( *x, &ja, &con_ja_val
					, &b_a, &norm_2_constr, &bnd_ja, &can_ignore_ja );
				assume_lin_dep_ja = false;	// Assume this initially.
				if( ja > 0 && act_set->is_init_fixed(ja) && qp.x_init()(ja) == bnd_ja )
					return_to_init_fixed = true;
				else
					return_to_init_fixed = false;
				// Print ja, bnd_ja, can_ignore_ja
				if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
					*out
						<< "\nja = " << ja	<< endl;
					if(ja) {
						*out
							<< "\ncon_ja_val           = "		<< con_ja_val
							<< "\nb_a                  = "		<< b_a
							<< "\nnorm_2_constr        = "		<< norm_2_constr
							<< "\nbnd_ja               = "		<< bnd_str(bnd_ja)
							<< "\ncan_ignore_ja        = "		<< bool_str(can_ignore_ja)
							<< "\nreturn_to_init_fixed = "		<< bool_str(return_to_init_fixed)
							<< endl
							;
					}
				}
				// Print first part of row for itr, change, indice, bound, violation
				if( (int)output_level == (int)OUTPUT_ITER_SUMMARY ) {
					*out
						<< setw(6)	<< itr
						<< setw(6)	<< act_set->q_hat()
						<< setw(8)	<< ( ja ? "ADD" : "-" )
						<< setw(6)	<< ja
						<< setw(10)	<< ( ja ? bnd_str(bnd_ja) : "-" )
						<< setw(20);
						if(ja)
							*out << (con_ja_val - b_a); 
						else
							*out << "-";
				}
				if( ja == 0 ) {
					// Todo: Implement iterative refinement if needed.
					return OPTIMAL_SOLUTION;	// current point is optimal.
				}
				const size_type sa = act_set->s_map(ja);
				if( sa != 0 || ( act_set->is_init_fixed(ja) && act_set->s_map(-ja) == 0 ) )
				{
					// Ohps! This constraint is already in the active set.
					std::ostringstream omsg;
					omsg
						<< "\nQPSchur::qp_algo(...) : Error, we have picked the constriant "
						<< "a(" << ja << ") with violation\n"
						<< "(a(ja)'*x - b_a) = (" << con_ja_val
						<< " - " << b_a << ") = " << (con_ja_val - b_a) << "\n"
						<< "to add to the active set but it is already part of the active set.\n"
						<< "This is an indication of instability in the calculations.\n"
						<< "The QP algorithm is terminated!\n";
					if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
						*out << omsg.str() << endl;
					}
					return SUBOPTIMAL_POINT;
				}
				else if(	act_set->all_dof_used_up()
						 && ::fabs( con_ja_val - b_a )
								/ std::_MAX(::fabs(con_ja_val),1.0) < feas_tol()	)
				{
					// Here all of the degrees of freedom are used up and the violated
					// constraint is not part of the active set.
					// In addition, it violation is very small.  Therefore we can assume
					// that this is a degenerate variable bound.
					if( act_set->qp().constraints().pick_violated_policy()
						== QP::Constraints::MOST_VIOLATED )
					{
						// This is the most violated constraint so just ignore it and we
						// are done.
						if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
							*out
								<< "\nWith all the dof used up, the inequality constriant a("
								<< ja
								<< ")' * x is the most violated "
								<< "but the violation = " << con_ja_val - b_a
								<< " is smaller than feas_tol = " << feas_tol_ << " so just assume "
								<< "this is a degenerate constraint and we are done.\n";
						}
						return OPTIMAL_SOLUTION;
					}
					else {
						if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
							*out
								<< "\nThe inequality constriant a("
								<< ja
								<< ")' * x is violated "
								<< "but the violation = " << con_ja_val - b_a
								<< " is smaller than feas_tol = " << feas_tol_ << " but it may not "
								<< "be the most violated constraints so set pick_violated_policy = "
								<< "MOST_VIOLATED and pick another violated constraint.\n";
						}
						act_set->qp().constraints().pick_violated_policy(
							QP::Constraints::MOST_VIOLATED );
						continue;	// Go back and pick the most violated constraint
					}
				} 
			}
			case UPDATE_ACTIVE_SET: {
				if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
					*out
						<< "\n*** UPDATE_ACTIVE_SET\n";
				}
				++(*num_adds);
				if( act_set->all_dof_used_up() || act_set->is_init_fixed(ja) ) {
					// All of the degrees of freedom are currently used up so we must
					// assume that we must remove one of these currently active
					// constraints and replace it with the violated constraint.
					// In this case we know that this will make the schur
					// complement singular so let's just skip the update
					// and set this now.  We also may be here if we are fixing
					// an initially fixed variable to some bound.
					assume_lin_dep_ja = true;
					if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
						if(act_set->all_dof_used_up()) {
							*out
								<< "\nAll of the degrees of freedom are used up so "
									"the constraint ja must be linearly dependent\n";
						}
						else {
							*out
								<< "\nThis is an initially fixed variable that was freed and "
									"now is being fixed again\n";
						}
					}
				}
				else {
					assume_lin_dep_ja = false;
					try {
						act_set->add_constraint( ja, bnd_ja, false, true );
						if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
							*out << "\nNew KKT system is nonsingular (LI constraints)\n";
						}
						if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
							dump_act_set_quantities( *act_set, *out );
						}
					}
					catch( const ActiveSet::BadUpdateException& excpt ) {
						// Constraint really is linearly dependent.
						if( (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
							*out
								<< "\nSchur complement update failed.  Constraint ja = "
								<< ja << " appears to be linearly dependent\n"
								<< "(" << excpt.what() << ")\n";
						}
						summary_lines_counter = summary_lines_counter_max;
						if( !(act_set->q_D_hat() + act_set->q_plus_hat()) ) {
							std::ostringstream omsg;
							omsg
								<< "\nQPSchur::qp_algo(...) : "
								<< "Error, constraint j = "<< ja << " is linearly dependent "
								<< "and there are no other constraints to drop.\n"
								<< "The QP must be infeasible\n";
							// Print omsg.
							if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
								*out << omsg.str() << endl;
							}
							return INFEASIBLE_CONSTRAINTS;
						}
						assume_lin_dep_ja = true;
					}
				}
				if( assume_lin_dep_ja && can_ignore_ja ) {
					act_set->qp().constraints().ignore( ja );
					next_step = PICK_VIOLATED_CONSTRAINT;
					continue;
				}
				// Infer the sign of the multiplier for the new constraint being added
				beta = ( con_ja_val > b_a ? +1.0 : -1.0 );
				if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
					*out
						<< "\nbeta = " << beta << endl;
				}
			}
			case COMPUTE_SEARCH_DIRECTION: {
				if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
					*out
						<< "\n*** COMPUTE_SEARCH_DIRECTION\n";
				}
				// All degrees of freedom used up?
				all_dof_used_up = act_set->all_dof_used_up();
				if( all_dof_used_up && (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
					*out
						<< "\nAll the degrees of freedom are used up.\n"
						<< "Therefore we don't need to compute p_v (set p_v = 0) and gamma_plus = beta*inf ...\n";
				}
				// Print end of row for rank
				if( (int)output_level == (int)OUTPUT_ITER_SUMMARY ) {
					*out
						<< setw(6)	<< ( assume_lin_dep_ja ? "LD" : "LI" ) << endl;
					out->flush();
					// Increment counter for summary header here.
					summary_lines_counter++;
				}
				if( assume_lin_dep_ja ) {
					//
					// The schur complement is not updated so we must compute
					// p_z_hat and p_v explicitly.
					//
					// If all the degrees of freedom
					// are used up then we know that the step of the primal variables
					// will be zero.  However, if m > 0 then v and p_v also contain
					// terms for the Lagrange multipliers for the equality constriants
					// but we don't need to compute these during the algorithm.
					// Therefore we can just set p_v = 0 and save a solve with Ko.
					// If the QP is feasible then a constraint will be dropped, the
					// KKT system will be updated and then v_plus will be computed
					// at the next iteration and be used to compute v so all is good.
					//
					if( act_set->is_init_fixed(ja) ) {
						// Fix a varaible that was fixed and then freed
						assert(0);	// ToDo: Finish this!
					}
					else {
						// Add a constraint that is not an initially fixed
						// variable bound.
						// 
						// p_z_hat = inv(S_hat) * ( - v_a + U_hat' * inv(Ko) * u_a )
						// 
						// p_v = inv(Ko) * ( -u_a - U_hat * p_z_hat )
						// 
						// gamma_plus = ( d_a - u_a'*v - v_a'*z_hat ) / ( u_a'*p_v + v_a'*p_z_hat )
						// 
						// ToDo: (9/25/00): Make u_a and v_a both sparse and combine the following code.
						// 
						if( ja <= n ) {
							// Fix an initially free variable
							//
							// u_a = e(ja) <: R^(n_R + m)
							// 
							// v_a = 0     <: R^(q_hat)
							//
							// d_a = b_a   <: R
							// 
							const EtaVector  u_a = EtaVector( ja, n_R + m );
							const value_type d_a = b_a;
							Vector t1;
							// t1 = inv(Ko) * u_a
							V_InvMtV( &t1, qp.Ko(), no_trans, u_a() );
							if( act_set->q_hat() ) {
								// t2 = U_hat'*t1
								Vector t2;
								V_MtV( &t2, act_set->U_hat(), trans, t1() );
								// p_z_hat = inv(S_hat) * t2
								V_InvMtV( &act_set->p_z_hat(), act_set->S_hat(), no_trans, t2() );
								// t1 = - u_a
								V_StV( &t1, -1.0, u_a() );
								// t1 += - U_hat * p_z_hat
								Vp_StMtV( &t1(), -1.0, act_set->U_hat(), no_trans, act_set->p_z_hat() );
								// p_v = inv(Ko) * t1
								if(!all_dof_used_up)
									V_InvMtV( &p_v, qp.Ko(), no_trans, t1() );
								else
									p_v = 0.0;
							}
							else {
								// p_v = -t1
								V_mV( &p_v, t1() );
							}
							// gamma_plus = ( d_a - u_a'*v) / ( u_a'*p_v )
							if(!all_dof_used_up)
								gamma_plus = ( d_a - dot(u_a(),*v) ) / dot(u_a(),p_v());
							else
								gamma_plus = beta * inf;
						}
						else {
							// Add a general inequality (or equality) constraint
							//
							// u_a = [ Q_R' * A_bar * e(ja) ] <: R^(n_R + m)
							//       [          0           ]
							// 
							// v_a = P_XF_hat' * A_bar * e_ja <: R^(q_hat)
							//
							// d_a = b_a - b_X' * (Q_X' * A_bar * e_ja) <: R
							//
							const EtaVector e_ja = EtaVector( ja, n + m_breve );
							Vector u_a( n_R + m );  // ToDo: Use workspace!
							// u_a(1:n_R) =  Q_R' * A_bar * e(ja)
							Vp_StPtMtV( &u_a(1,n_R), 1.0, qp.Q_R(), trans
										, qp.constraints().A_bar(), no_trans, e_ja(), 0.0 );
							// u_a(n_R+1:n_R+m) = 0.0
							if(m)
								u_a(n_R+1,n_R+m) = 0.0;
							// t0 = Q_X' * A_bar * e_ja
							Vector t0(n-n_R);
							if( n > n_R )
								Vp_StPtMtV( &t0(), 1.0, qp.Q_X(), trans
											, qp.constraints().A_bar(), no_trans, e_ja(), 0.0 );
							// d_a = b_a - b_X'*t0
							const value_type
								d_a = b_a - ( n > n_R ? dot( qp.b_X(), t0() ) : 0.0 );
							// t1 = inv(Ko) * u_a
							Vector t1;
							V_InvMtV( &t1, qp.Ko(), no_trans, u_a );
							if( act_set->q_hat() ) {
								// t2 = U_hat'*t1
								Vector t2;
								V_MtV( &t2, act_set->U_hat(), trans, t1() );
								// v_a = P_XF_hat' * A_bar * e_ja
								Vector v_a(act_set->q_hat()); // ToDo: Use workspace for this!
								Vp_StPtMtV( &v_a(), 1.0, act_set->P_XF_hat(), trans
											, qp.constraints().A_bar(), no_trans, e_ja(), 0.0 );
								// t2 += -v_a
								Vp_StV( &t2(), -1.0, v_a() );
								// p_z_hat = inv(S_hat) * t2
								V_InvMtV( &act_set->p_z_hat(), act_set->S_hat(), no_trans, t2() );
								if(!all_dof_used_up) {
									// t1 = - u_a
									V_StV( &t1, -1.0, u_a() );
								    // t1 += - U_hat * p_z_hat
									Vp_StMtV( &t1(), -1.0, act_set->U_hat(), no_trans, act_set->p_z_hat() );
								    // p_v = inv(Ko) * t1
									V_InvMtV( &p_v, qp.Ko(), no_trans, t1() );
								}
								else {
									p_v = 0.0;
								}
								// gamma_plus = ( d_a - u_a'*v - v_a'*z_hat ) / ( u_a'*p_v + v_a'*p_z_hat )
								if(!all_dof_used_up)
									gamma_plus = ( ( d_a - dot(u_a,*v) - dot(v_a(),act_set->z_hat()) )
												   / ( dot(u_a,p_v()) + dot(v_a(),act_set->p_z_hat()) ) );
								else
									gamma_plus = beta * inf;
							}
							else {
								// p_v = -t1
								if(!all_dof_used_up)
									V_mV( &p_v, t1() );
								else
									p_v = 0.0;
								// gamma_plus = ( d_a - u_a'*v) / ( u_a'*p_v )
								if(!all_dof_used_up)
									gamma_plus = ( d_a - dot(u_a,*v) ) / dot(u_a,p_v());
								else
									gamma_plus = beta * inf;
							}
						}
					}
					// Print the steps p_v and p_z_hat
					if(act_set->q_hat()) {
						if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
							*out
								<< "\n||p_z_hat||inf = " << norm_inf(act_set->p_z_hat()) << endl;
						}
						if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES )
						{
							*out << "\np_z_hat =\n" << act_set->p_z_hat();
						}
					}
					if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
						*out
							<< "\n||p_v||inf = " << norm_inf(p_v()) << endl;
					}
					if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES )
					{
						*out << "\np_v =\n" << p_v();
					}
					if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
						*out
							<< "\ngamma_plus = " << gamma_plus << endl;
					}
					// Compute step for mu_D_hat
					if( act_set->q_D_hat() ) {
						calc_p_mu_D( *act_set, p_v(), act_set->p_z_hat(), &act_set->p_mu_D_hat() );
						if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
							*out
								<< "\n||p_mu_D_hat||inf = " << norm_inf(act_set->p_mu_D_hat()) << std::endl;
						}
						if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
							*out
								<< "\np_mu_D_hat =\n" << act_set->p_mu_D_hat();
						}
					}
				}				
				else {
					// The new schur complement is already updated so compute
					// the solution outright.
				
					// Compute z_hat_plus, v_plus

					// z_hat_plus = inv(S_hat) * ( d_hat - U_hat' * vo  )
					const size_type q_hat = act_set->q_hat();
					z_hat_plus.resize( q_hat );
					calc_z( act_set->S_hat(), act_set->d_hat(), act_set->U_hat(), vo()
						, &z_hat_plus() );
					if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
						*out
							<< "\n||z_hat_plus||inf = " << norm_inf(z_hat_plus()) << std::endl;
					}
					if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
						*out
							<< "\nz_hat_plus =\n" << z_hat_plus();
					}
					// v_plus = inv(Ko) * (fo - U_hat * z_hat_plus)
					calc_v( act_set->qp().Ko(), act_set->qp().fo(), act_set->U_hat(), z_hat_plus()
						, &v_plus() );
					if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
						*out
							<< "\n||v_plus||inf = " << norm_inf(v_plus()) << std::endl;
					}
					if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
						*out
							<< "\nv_plus =\n" << v_plus();
					}
					// Compute p_z_hat (change in z_hat w.r.t newly added constriant multiplier)
					VectorSlice p_z_hat = act_set->p_z_hat();
					// p_z_hat = z_hat_plus - z_hat
					V_VmV( &p_z_hat(), z_hat_plus(), act_set->z_hat() );
					// p_v = v_plus - v
					V_VmV( &p_v(), v_plus(), *v );
					// p_mu_D_hat
					if( act_set->q_D_hat() )
						calc_p_mu_D( *act_set, p_v(), p_z_hat(), &act_set->p_mu_D_hat() );
					// gamma_plus
					if( act_set->is_init_fixed(ja) && act_set->qp().x_init()(ja) == bnd_ja ) {
						gamma_plus = act_set->p_mu_D_hat()( act_set->q_D_hat() );
					}
					else {
						gamma_plus = z_hat_plus(q_hat);
					}
					// p_z_hat = p_z_hat / gamma_plus
					Vt_S( &p_z_hat(), 1.0 / gamma_plus );
					// p_v = p_v / gamma_plus
					Vt_S( &p_v(), 1.0 / gamma_plus );
					// Print gama_plus, p_z_hat, p_v and p_mu_D_hat
					if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
						*out
							<< "\ngamma_plus = " << gamma_plus << std::endl;
					}
					if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
						*out
							<< "\n||p_z_hat||inf = " << norm_inf(p_z_hat()) << std::endl;
					}
					if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
						*out
							<< "\np_z_hat =\n" << p_z_hat();
					}
					if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
						*out
							<< "\n||p_v||inf = " << norm_inf(p_v()) << std::endl;
					}
					if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
						*out
							<< "\np_v =\n" << p_v();
					}
					if( act_set->q_D_hat() ) {
						// p_mu_D_hat = p_mu_D_hat / gamma_plus
						Vt_S( &act_set->p_mu_D_hat(), 1.0 / gamma_plus ); 
						if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
							*out
								<< "\n||p_mu_D_hat||inf =\n" << norm_inf(act_set->p_mu_D_hat()) << std::endl;
						}
						if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
							*out
								<< "\np_mu_D_hat =\n" << act_set->p_mu_D_hat();
						}
					}
				}
			}
			case COMPUTE_STEP_LENGTHS: {

				if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
					*out
						<< "\n*** COMPUTE_STEP_LENGTHS\n";
				}
				// Compute the dual infeasibility scaling
				const size_type q_hat = act_set->q_hat();
				dual_infeas_scale = 1.0;
				if( q_hat )
					dual_infeas_scale = std::_MAX( dual_infeas_scale, norm_inf( act_set->z_hat() ) );
				if( m )
					dual_infeas_scale = std::_MAX( dual_infeas_scale, norm_inf( (*v)(n_R+1,n_R+m) ) );
				if( act_set->q_D_hat() )
					dual_infeas_scale = std::_MAX( dual_infeas_scale, norm_inf( act_set->mu_D_hat() ) );
				
				// Primal step length, t_P = beta * gamma_plus, z_plus = [ z_hat_plus; gama_plus ].
				// Or constraint ja is linearly dependent in which case p_x is zero so
				// t_P is infinite.
				t_P = beta * gamma_plus;	// Could be < 0
				if( t_P < 0.0 ) {
					if( ::fabs( t_P / (norm_2_constr*dual_infeas_scale) )  <= dual_infeas_tol() ) {
						if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
							*out
								<< "\nWarning, A near degenerate inequality constraint ja = " << ja
									<< " is being added that has the wrong sign with:\n"
								<< "    t_P                     = " << t_P 				<< std::endl
								<< "    dual_infeas_scale       = " << dual_infeas_scale	<< std::endl
								<< "    norm_2_constr            = " << norm_2_constr   << std::endl
								<< "    |t_P/(norm_2_constr*dual_infeas_scale)| = "
									<< ::fabs(t_P/(norm_2_constr*dual_infeas_scale))
									<< " <= dual_infeas_tol = " << dual_infeas_tol() << std::endl
								<< "therefore we will adjust things and keep going.\n";
						}
						assert(0);	// ToDo: Finish this!
					}
					else {
						std::ostringstream omsg;
						omsg
							<< "QPSchur::qp_algo(...) :\n"
							<< "Error, an inequality constraint ja = " << ja
							<< " is being added that has the wrong sign and is not near degenerate with:\n"
							<< "    t_P                     = " << t_P 				<< std::endl
							<< "    dual_infeas_scale       = " << dual_infeas_scale	<< std::endl
							<< "    norm_2_constr            = " << norm_2_constr   << std::endl
							<< "    |t_P/(norm_2_constr*dual_infeas_scale)| = "
								<< ::fabs(t_P/(norm_2_constr*dual_infeas_scale))
								<< " < -dual_infeas_tol = " << dual_infeas_tol() << std::endl
							<< "There may be serious illconditioning in the problem.\n";
						if( out && (int)output_level >= (int)OUTPUT_BASIC_INFO )
							*out << omsg.str();
						return DUAL_INFEASIBILITY;
					}
				}
				t_P = beta * gamma_plus;	// Now guaranteed to be > 0
				if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
					*out
						<< "\nt_P = " << t_P << endl;
				}

				/////////////////////////////////////////////////////////////////////////
				// Dual step length.  Largest step t that does not cause violation in
				// dual feasibility (i.e. lagrange multipliers for inequalities are
				// dual feasible, or primal optimal ).
				// lambda_hat_new = lambda_hat + beta * t_D * p_lambda_hat must be dual feasible.
				t_D = inf;
				jd = 0;
				value_type max_feas_viol = 0.0; // Remember the amount of violation.
				int j_degen = 0;	// remember which (if any) constraint was near
									// degenerate and had an incorrect sign.
				QPSchurPack::EBounds bnd_jd;	// The bound of the constraint to be dropped.

				// Search through Lagrange multipliers in z_hat
				if( act_set->q_hat() ) {

					VectorSlice z_hat = act_set->z_hat();
					VectorSlice p_z_hat = act_set->p_z_hat();
					VectorSlice::const_iterator
						z_itr		= const_cast<const VectorSlice&>(z_hat).begin(),
						p_z_itr		= const_cast<const VectorSlice&>(p_z_hat).begin();
					const size_type
						qq = assume_lin_dep_ja || (!assume_lin_dep_ja && return_to_init_fixed)
							? q_hat : q_hat - 1;
					// Print header for s, j, z_hat(s), p_z_hat(s), bnds(s), t, t_D, jd
					if( qq > 0 && (int)output_level >= (int)OUTPUT_ACT_SET ) {
						if(out) *out << std::setprecision(dbl_prec);
						*out
							<< "\nComputing the maximum step for multiplers for dual feasibility\n\n"
							<< setw(5)	<< "s"
							<< setw(5)	<< "j"
							<< setw(20)	<< "z_hat"
							<< setw(20)	<< "p_z_hat"
							<< setw(20)	<< "bnd"
							<< setw(20)	<< "t"
							<< setw(20)	<< "t_D"
							<< setw(5)	<< "jd"	<< endl
							<< setw(5)	<< "----"
							<< setw(5)	<< "----"
							<< setw(20)	<< "--------------"
							<< setw(20)	<< "--------------"
							<< setw(20)	<< "--------------"
							<< setw(20)	<< "--------------"
							<< setw(20)	<< "--------------"
							<< setw(5)	<< "----"	<< endl;
					}
					for( int s = 1; s <= qq; ++s, ++z_itr, ++p_z_itr) {
						int j = act_set->ij_map(s);
						if( j > 0 ) {
							namespace ns = QPSchurPack;
							ns::EBounds bnd = act_set->bnd(s);
							// Print first part of row for s, j, z_hat(s), p_z_hat(s), bnds(s) ....
							if( (int)output_level >= (int)OUTPUT_ACT_SET ) {
								*out
									<< setw(5)	<< s
									<< setw(5)	<< j
									<< setw(20)	<< *z_itr
									<< setw(20)	<< *p_z_itr
									<< setw(20)	<< bnd_str(bnd);
							}
							value_type t = inf;
							bool j_is_degen = false;
							if( (bnd == ns::LOWER && *z_itr > 0.0) || (bnd == ns::UPPER && *z_itr < 0.0) ) {
								if( (int)output_level >= (int)OUTPUT_BASIC_INFO && !warned_degeneracy ) {
									*out
										<< "\nWarning, possible degeneracy and numerical instability detected.\n";
									warned_degeneracy = true;
								}
								assert(0);	// ToDo: Finish this!
							}
							const value_type feas_viol = beta*(*p_z_itr);
							if( bnd == ns::LOWER && feas_viol <= 0.0 )
								;	// dual feasible for all t > 0
							else if( bnd == ns::UPPER && feas_viol >= 0.0 )
								;	// dual feasible for all t > 0
							else {
								// finite t.
								t = -beta*(*z_itr)/(*p_z_itr);
								if( t < t_D ) {	// remember minimum step length
									t_D = t;
									jd = j;
									if(j_is_degen) j_degen = j;
									max_feas_viol = feas_viol;
									bnd_jd = bnd;
								}
							}
							// Print rest of row for ... t, t_D, jd
							if( (int)output_level >= (int)OUTPUT_ACT_SET ) {
								*out
									<< setw(20)	<< t
									<< setw(20)	<< t_D
									<< setw(5)	<< jd	<< endl;
							}
						}
					}
					if(out) *out << std::setprecision(prec_saved);
				}
					
				// Search through Lagrange multipliers in mu_D_hat
				if( act_set->q_D_hat() ) {
					const size_type q_D_hat = act_set->q_D_hat();
					VectorSlice mu_D_hat = act_set->mu_D_hat();
					VectorSlice p_mu_D_hat = act_set->p_mu_D_hat();
					const GenPermMatrixSlice &Q_XD_hat = act_set->Q_XD_hat();
					VectorSlice::const_iterator
						mu_D_itr		= const_cast<const VectorSlice&>(mu_D_hat).begin(),
						p_mu_D_itr		= const_cast<const VectorSlice&>(p_mu_D_hat).begin();
					GenPermMatrixSlice::const_iterator
						Q_XD_itr		= Q_XD_hat.begin();
					const size_type
						qD = assume_lin_dep_ja && return_to_init_fixed ? q_D_hat-1 : q_D_hat;
					// Print header for k, i, mu_D_hat(k), p_mu_D_hat(k), x_init(k), t, t_D, jd
					if( qD > 0 && (int)output_level >= (int)OUTPUT_ACT_SET ) {
						if(out) *out << std::setprecision(dbl_prec);
						*out
							<< "\nComputing the maximum step for multiplers for dual feasibility\n\n"
							<< setw(5)	<< "k"
							<< setw(5)	<< "i"
							<< setw(20)	<< "mu_D_hat"
							<< setw(20)	<< "p_mu_D_hat"
							<< setw(20)	<< "x_init"
							<< setw(20)	<< "t"
							<< setw(20)	<< "t_D"
							<< setw(5)	<< "jd"	<< endl
							<< setw(5)	<< "----"
							<< setw(5)	<< "----"
							<< setw(20)	<< "--------------"
							<< setw(20)	<< "--------------"
							<< setw(20)	<< "--------------"
							<< setw(20)	<< "--------------"
							<< setw(20)	<< "--------------"
							<< setw(5)	<< "----"	<< endl;
					}
					for( int k = 1; k <= qD; ++k, ++mu_D_itr, ++p_mu_D_itr, ++Q_XD_itr )
					{
						int i = Q_XD_itr->row_i();	// ith fixed variable
						{
							namespace ns = QPSchurPack;
							ns::EBounds bnd = qp.x_init()(i);
							// Print first part of row for s, j, z_hat(s), p_z_hat(s), bnds(s) ....
							if( (int)output_level >= (int)OUTPUT_ACT_SET ) {
								*out
									<< setw(5)	<< k
									<< setw(5)	<< i
									<< setw(20)	<< *mu_D_itr
									<< setw(20)	<< *p_mu_D_itr
									<< setw(20)	<< bnd_str(bnd);
							}
							value_type t = inf;
							bool j_is_degen = false;
							if( (bnd == ns::LOWER && *mu_D_itr > 0.0) || (bnd == ns::UPPER && *mu_D_itr < 0.0) ) {
								if( (int)output_level >= (int)OUTPUT_BASIC_INFO && !warned_degeneracy ) {
									*out
										<< "\nWarning, possible degeneracy and numerical instability detected.\n";
									warned_degeneracy = true;
								}
								assert(0);	// ToDo: Finish this!
							}
							const value_type feas_viol = beta*(*p_mu_D_itr);
							if( bnd == ns::LOWER && feas_viol <= 0.0 )
								;	// dual feasible for all t > 0
							else if( bnd == ns::UPPER && feas_viol >= 0.0 )
								;	// dual feasible for all t > 0
							else {
								// finite t.
								t = -beta*(*mu_D_itr)/(*p_mu_D_itr);
								if( t < t_D ) {	// remember minimum step length
									t_D = t;
									jd = -i;
									if(j_is_degen) j_degen = jd;
									max_feas_viol = feas_viol;
									bnd_jd = bnd;
								}
							}
							// Print rest of row for ... t, t_D, jd
							if( (int)output_level >= (int)OUTPUT_ACT_SET ) {
								*out
									<< setw(20)	<< t
									<< setw(20)	<< t_D
									<< setw(5)	<< jd	<< endl;
							}
						}
					}
					if(out) *out << std::setprecision(prec_saved);
				}

				// Print t_D, jd
				if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
					*out
						<< "\nt_D = "	<< t_D	<< endl
						<< "jd = "		<< jd	<< endl;
				}
				if( jd == j_degen && jd != 0 && t_D < t_P ) {
					if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
						*out
							<< "\nWarning, the near degenerate constraint j = "
							<< jd << " which had the incorrect sign\nand was adjusted "
							<< "was selected to be dropped from the active set.\n";
					}
				}

				// Print start of row for itr, change, indice, bound, violation
				if( t_D < t_P && (int)output_level == (int)OUTPUT_ITER_SUMMARY ) {
					*out
						<< setw(6)	<< itr
						<< setw(6)	<< act_set->q_hat()
						<< setw(8)	<< "DROP"
						<< setw(6)	<< jd
						<< setw(10)	<< bnd_jd
						<< setw(20)	<< max_feas_viol;
				}
			}
			case TAKE_STEP: {
				if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
					*out
						<< "\n*** TAKE_STEP\n";
				}
				if( t_P >= huge_primal_step() && t_D >= huge_dual_step() ) {
					if( 	act_set->all_dof_used_up()
						&&	::fabs( con_ja_val - b_a )
								/ std::_MAX(::fabs(con_ja_val),1.0) < loose_feas_tol_	)
					{
						if( act_set->qp().constraints().pick_violated_policy()
							== QP::Constraints::MOST_VIOLATED )
						{
							// This is the most violated constraint so  perform a step of
							// iterative refinement and quit.
							if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
								*out
									<< "\nWith all the dof used up, the inequality constriant a("
									<< ja
									<< ")' * x is the most violated "
									<< "but the violation = " << con_ja_val - b_a
									<< " is smaller than loose_feas_tol = " << loose_feas_tol_
									<< " but larger than feas_tol = " << feas_tol_
									<< " and the iteration failed so just assume "
									<< "this is a degenerate constraint and we are done.\n";
							}
							return OPTIMAL_SOLUTION;
						}
						else {
							act_set->qp().constraints().pick_violated_policy(
								QP::Constraints::MOST_VIOLATED );
							continue;	// Go back and pick the most violated constraint
						}
					}
					// This is an infeasible QP
					if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
						std::ostringstream omsg;
						omsg
							<< "QPSchur::qp_algo(...) : "
							<< "Error, QP is infeasible, inconsistent constraint "
							<< ja << " detected";
						// Print omsg.
						*out << "\n*** (a) " << omsg.str() << endl;
					}
					return INFEASIBLE_CONSTRAINTS;
				}
				else if( t_P > t_D ) {
					if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
						if( t_P >= huge_primal_step() ) {
							*out
								<< "\n*** (b) Dual Step (t_P = " << t_P << " >= huge_primal_step = "
									<< huge_primal_step() << endl;
						}
						else {
							*out
								<< "\n*** (b) Partial Primal-Dual Step\n";
						}
					}
					if( assume_lin_dep_ja ) {
						act_set->drop_add_constraints( jd, ja, bnd_ja, true );
					}
					else {
						act_set->drop_constraint( jd );
					}
					// z_hat = z_hat + beta * t_D * p_z_hat
					if(act_set->q_hat())
						Vp_StV( &act_set->z_hat(), beta * t_D, act_set->p_z_hat() );
					// v = v + beta * t_D * p_v
					Vp_StV( v, beta * t_D, p_v() );
					// mu_D_hat = mu_D_hat + beta * t_D * p_mu_D_hat
					if(act_set->q_D_hat())
						Vp_StV( &act_set->mu_D_hat(), beta * t_D, act_set->p_mu_D_hat() );

					++(*num_drops);
					last_jd = jd;

					if( (int)output_level >= (int)OUTPUT_ITER_STEPS )
					{
						*out
							<< "\nUpdated primal and dual variables:\n"
							<< "\n||v||inf           = " << norm_inf(*v)
							<< "\n||z_hat||inf       = " << norm_inf(act_set->z_hat())
							<< "\nmax(|mu_D_hat(i)|) = " << norm_inf(act_set->mu_D_hat())
							<< "\nmin(|mu_D_hat(i)|) = " << min_abs(act_set->mu_D_hat())
							<< endl;
					}
					if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES )
					{
						*out << "\nv = \n" << *v << endl;
						dump_act_set_quantities( *act_set, *out );
					}
					assume_lin_dep_ja = false;
					next_step = COMPUTE_SEARCH_DIRECTION;
					continue;
				}
				else {	// t_P < t_D
					if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
						*out
							<< "\n*** (c) Full Primal-Dual Step\n";
					}
					++(*num_adds);
					if( !assume_lin_dep_ja ) {
						act_set->z_hat() 	= z_hat_plus;
						*v 					= v_plus;
					}
					else {
						act_set->add_constraint( ja, bnd_ja, true, true );
						// z_hat = z_hat + beta * t_P * p_z_hat
						if(act_set->q_hat())
							Vp_StV( &act_set->z_hat(), beta * t_P, act_set->p_z_hat() );
						// v = v + beta * t_P * p_v
						Vp_StV( v, beta * t_P, p_v() );
					}
					// mu_D_hat = mu_D_hat + beta * t_P * p_mu_D_hat
					if(act_set->q_D_hat())
						Vp_StV( &act_set->mu_D_hat(), beta * t_P, act_set->p_mu_D_hat() );

					if( (int)output_level >= (int)OUTPUT_ITER_STEPS )
					{
						*out
							<< "\n||v||inf           = " << norm_inf(*v)
							<< "\n||z_hat||inf       = " << norm_inf(act_set->z_hat())
							<< "\nmax(|mu_D_hat(i)|) = " << norm_inf(act_set->mu_D_hat())
							<< "\nmin(|mu_D_hat(i)|) = " << min_abs(act_set->mu_D_hat())
							<< endl;
					}
					if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES )
					{
						*out << "\nv = \n" << *v << endl;
						if( assume_lin_dep_ja )
							dump_act_set_quantities( *act_set, *out );
						else
							*out
								<< "\nz_hat =\n" << act_set->z_hat()
								<< "\nmu_D_hat =\n" << act_set->mu_D_hat() << endl;
					}
					next_step = PICK_VIOLATED_CONSTRAINT;
					continue;
				}
			}
			default:
				assert(0);	// only a local programming error
		}
	}

	} // end try
	catch( std::exception& excpt ) {
		if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
			*out
				<< "\n*** Caught a standard exception :\n"
				<< excpt.what() << endl
				<< "*** Rethrowing the exception ...\n";
		}
		if(out) *out << std::setprecision(prec_saved);
		throw;
	}
	catch(...) {
		if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
			*out
				<< "\n*** Caught an unknown exception.  Rethrowing the exception ...\n";
		}
		if(out) *out << std::setprecision(prec_saved);
		throw;
	}

	if(out) *out << std::setprecision(prec_saved);

	// If you get here then the maximum number of QP iterations has been exceeded
	return MAX_ITER_EXCEEDED;
}

void QPSchur::set_x( const ActiveSet& act_set, const VectorSlice& v, VectorSlice* x )
{
	using BLAS_Cpp::no_trans;
	using LinAlgOpPack::V_MtV;
	using LinAlgOpPack::Vp_MtV;
	
	// x = Q_R * v(1:n_R) + Q_X * b_X + P_XF_hat * z_hat
	if( act_set.qp().n_R() )
		V_MtV( x, act_set.qp().Q_R(), no_trans, v(1,act_set.qp().n_R()) );
	if( act_set.qp().n() > act_set.qp().n_R() )
		Vp_MtV( x, act_set.qp().Q_X(), no_trans, act_set.qp().b_X() );
	if( act_set.q_hat() )
		Vp_MtV( x, act_set.P_XF_hat(), no_trans, act_set.z_hat() );
}

void QPSchur::set_multipliers( const ActiveSet& act_set, const VectorSlice& v
	, SpVector* mu, VectorSlice* lambda, SpVector* lambda_breve )
{
	using BLAS_Cpp::no_trans;
	using LinAlgOpPack::V_MtV;
	using LinAlgOpPack::Vp_MtV;
	using LinAlgOpPack::Vp_StMtV;
	using SparseLinAlgPack::V_MtV;
	using SparseLinAlgPack::Vp_MtV;
	namespace GPMSTP = SparseLinAlgPack::GenPermMatrixSliceIteratorPack;

	const size_type
		n = act_set.qp().n(),
		n_R = act_set.qp().n_R(),
		m = act_set.qp().m(),
		m_breve = act_set.qp().constraints().m_breve();

	// mu = P_plus_hat(1:n,:) * z_hat + Q_XD_hat * mu_D + (steps for initially fixed
	// 		variables fixed to the other bounds).

	typedef SpVector::element_type ele_t;
	mu->resize( n, n-m );	// Resize for the maxinum number of fixed variables possible

	// mu += P_plus_hat(1:n,:) * z_hat
	if( act_set.q_plus_hat() )
		Vp_MtV( mu, act_set.P_plus_hat().create_submatrix(Range1D(1,n),GPMSTP::BY_ROW)
			, no_trans, act_set.z_hat() );
	// mu += Q_XD_hat * mu_D_hat
	if( act_set.q_D_hat() )
		Vp_MtV( mu, act_set.Q_XD_hat(), no_trans, act_set.mu_D_hat() );
	// Add multipliers for initially fixed variables fixed to the other bounds.
	if( act_set.q_C_hat() ) {
		assert(0);	// ToDo: Finish This!
	}

	mu->sort();
	
	// lambda = v(n_R+1,n_R+m)
	if( m ) {
		*lambda = v(n_R+1,n_R+m);
	}

	// lambda_breve = P_plus_hat(n+1:n+m_breve,:) * z_hat
	if( m_breve && act_set.q_plus_hat() )
		V_MtV( lambda_breve
			, act_set.P_plus_hat().create_submatrix(Range1D(n+1,n+m_breve),GPMSTP::BY_ROW)
			, no_trans, act_set.z_hat() );
}

// private static member functions for QPSchur

void QPSchur::dump_act_set_quantities(
	   const ActiveSet& act_set, std::ostream& out
	 , bool print_S_hat )
{
	using std::endl;
	using std::setw;
	using std::left;
	using std::right;

	const QPSchurPack::QP
		&qp = act_set.qp();
	const QPSchurPack::Constraints
		&constraints = qp.constraints();

	const int int_w = 12;
	const char int_ul[] = "----------"; 
	const int dbl_prec = 6;
	const int dbl_w = 20;
	const char dbl_ul[] = "------------------";

	const int prec_saved = out.precision();

	try {
    out << "\n*** Dumping the current active set ***\n"
		<< "\nDimensions of the current active set:\n"
		<< "\nn           = " << right << setw(int_w) << qp.n()					<< " (Number of unknowns)"
		<< "\nn_R         = " << right << setw(int_w) << qp.n_R()				<< " (Number of initially free variables in Ko)"
		<< "\nm           = " << right << setw(int_w) << qp.m()					<< " (Number of initially fixed variables not in Ko)"
		<< "\nm_breve     = " << right << setw(int_w) << constraints.m_breve()	<< " (Number of extra general equality/inequality constriants)"
		<< "\nq_hat       = " << right << setw(int_w) << act_set.q_hat()		<< " (Number of augmentations to the initial KKT system Ko)"
		<< "\nq_plus_hat  = " << right << setw(int_w) << act_set.q_plus_hat()	<< " (Number of added variable bounds and general constraints)"
		<< "\nq_F_hat     = " << right << setw(int_w) << act_set.q_F_hat()		<< " (Number of initially fixed variables not at their initial bound)"
		<< "\nq_C_hat     = " << right << setw(int_w) << act_set.q_C_hat()		<< " (Number of initially fixed variables at the other bound)"
		<< "\nq_D_hat     = " << right << setw(int_w) << act_set.q_D_hat()		<< " (Number of initially fixed variables still fixed at initial bound)"
		<< endl;

	// Print table of quantities in augmented KKT system
	out	<< "\nQuantities for augmentations to the initial KKT system:\n";
	out << std::setprecision(dbl_prec);
	const size_type q_hat = act_set.q_hat();
	out	<< endl
		<< right << setw(int_w) << "s"
		<< right << setw(int_w) << "ij_map(s)"
		<< right << setw(int_w) << "bnd(s)"
		<< right << setw(dbl_w) << "constr_norm(s)"
		<< right << setw(dbl_w) << "d_hat(s)"
		<< right << setw(dbl_w) << "z_hat(s)"
		<< right << setw(dbl_w) << "p_z_hat(s)"
		<< endl;
	out	<< right << setw(int_w) << int_ul
		<< right << setw(int_w) << int_ul
		<< right << setw(int_w) << int_ul
		<< right << setw(dbl_w) << dbl_ul
		<< right << setw(dbl_w) << dbl_ul
		<< right << setw(dbl_w) << dbl_ul
		<< right << setw(dbl_w) << dbl_ul
		<< endl;
	{for( size_type s = 1; s <= q_hat; ++s ) {
		out	<< right << setw(int_w) << s
			<< right << setw(int_w) << act_set.ij_map(s)
			<< right << setw(int_w) << bnd_str(act_set.bnd(s))
			<< right << setw(dbl_w) << act_set.constr_norm(s)
			<< right << setw(dbl_w) << act_set.d_hat()(s)
			<< right << setw(dbl_w) << act_set.z_hat()(s)
			<< right << setw(dbl_w) << act_set.p_z_hat()(s)
			<< endl;
	}}
	out << std::setprecision(prec_saved);
	
	// Print P_XF_hat, P_plus_hat, U_hat and S_hat
	out	<< "\nP_XF_hat =\n" 	<< act_set.P_XF_hat();
	out	<< "\nP_plus_hat =\n" 	<< act_set.P_plus_hat();
	out	<< "\nU_hat =\n" 		<< act_set.U_hat();
	if(print_S_hat)
		out	<< "\nS_hat =\n" 		<< act_set.S_hat();
	
	// Print table of multipliers for q_D_hat
	out	<< "\nQuantities for initially fixed variables which are still fixed at their initial bound:\n";
	out << std::setprecision(dbl_prec);
	const size_type q_D_hat = act_set.q_D_hat();
	out	<< endl
		<< right << setw(int_w) << "k"
		<< right << setw(int_w) << "l_fxfx(k)"
		<< right << setw(dbl_w) << "mu_D_hat(k)"
		<< right << setw(dbl_w) << "p_mu_D_hat(s)"
		<< endl;
	out	<< right << setw(int_w) << int_ul
		<< right << setw(int_w) << int_ul
		<< right << setw(dbl_w) << dbl_ul
		<< right << setw(dbl_w) << dbl_ul
		<< endl;
	{for( size_type k = 1; k <= q_D_hat; ++k ) {
		out	<< right << setw(int_w) << k
			<< right << setw(int_w) << act_set.l_fxfx(k)
			<< right << setw(dbl_w) << act_set.mu_D_hat()(k)
			<< right << setw(dbl_w) << act_set.p_mu_D_hat()(k)
			<< endl;
	}}
	out << std::setprecision(prec_saved);
	
	// Print Q_XD_hat
	out	<< "\nQ_XD_hat =\n" << act_set.Q_XD_hat();

	out << "\n*** End dump of current active set ***\n";

	}	// end try
	catch(...) {
		out << std::setprecision(prec_saved);
		throw;
	}

	out << std::setprecision(prec_saved);
}

// QPSchurPack::QP

void QPSchurPack::QP::dump_qp( std::ostream& out )
{
	using std::endl;
	using std::setw;
	using std::left;
	using std::right;

	const Constraints
		&constraints = this->constraints();

	const size_type
		n = this->n(),
		n_R = this->n_R(),
		m = this->m(),
		m_breve = constraints.m_breve();

	out	<< "\n*** Original QP ***\n"
		<< "\nn       = " << n
		<< "\nm       = " << m
		<< "\nm_breve = " << m_breve
		<< endl;
	out	<< "\ng =\n" << g();
	out	<< "\nG =\n" << G();
	if(m) {
		out	<< "\nA =\n" << A();
		// Le'ts recover c from fo(n_R+1:n_R+m) = c - A' * Q_X * b_x
		assert(0);	// Finish this!
	}
	out	<< "\nA_bar =\n" << constraints.A_bar();
	// Get c_L_bar and c_U_bar
	Vector c_L_bar(n+m_breve), c_U_bar(n+m_breve);
	{for( size_type j = 1; j <= n+m_breve; ++j ){
		c_L_bar(j) = constraints.get_bnd(j,LOWER);
		c_U_bar(j) = constraints.get_bnd(j,UPPER);
	}}
	out	<< "\nc_L_bar =\n" << c_L_bar();
	out	<< "\nc_U_bar =\n" << c_U_bar();
	
	out	<< "\n*** Initial KKT system (fixed and free variables) ***\n"
		<< "\nn_R = " << n_R
		<< endl;
	out	<< "\nb_X =\n" << b_X();
	out	<< "\nQ_R =\n" << Q_R();
	out	<< "\nQ_X =\n" << Q_X();
	out	<< "\nKo =\n" << Ko();
	out	<< "\nfo =\n" << fo();

}

}	// end namespace ConstrainedOptimizationPack
