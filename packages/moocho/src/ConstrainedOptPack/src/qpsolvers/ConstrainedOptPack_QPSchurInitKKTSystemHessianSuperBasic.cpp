// /////////////////////////////////////////////////////////
// QPSchurInitKKTSystemHessianSuperBasic.cpp

#include "ConstrainedOptimizationPack/include/QPSchurInitKKTSystemHessianSuperBasic.h"
#include "ConstrainedOptimizationPack/include/MatrixHessianSuperBasic.h"
#include "SparseLinAlgPack/include/GenPermMatrixSlice.h"
#include "SparseLinAlgPack/include/GenPermMatrixSliceOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace ConstrainedOptimizationPack {

void QPSchurInitKKTSystemHessianSuperBasic::initialize_kkt_system(
	const VectorSlice&    g
	,const MatrixWithOp&  G
	,value_type           etaL
	,const SpVectorSlice& dL
	,const SpVectorSlice& dU
	,const MatrixWithOp*  F
	,BLAS_Cpp::Transp     trans_F
	,const VectorSlice*   f
	,i_x_free_t*          i_x_free
	,i_x_fixed_t*         i_x_fixed
	,bnd_fixed_t*         bnd_fixed
	,j_f_decomp_t*        j_f_decomp
	,Vector*              b_X
	,Ko_ptr_t*            Ko
	,Vector*              fo
	) const
{
	using BLAS_Cpp::trans;
	using DynamicCastHelperPack::dyn_cast;
	using LinAlgOpPack::V_mV;
	using LinAlgOpPack::V_StMtV;
	using SparseLinAlgPack::Vp_StMtV;

	// Validate type of and convert G
#ifdef _WINDOWS
	const MatrixHessianSuperBasic&
		G_super = dynamic_cast<const MatrixHessianSuperBasic&>(G);
#else
	const MatrixHessianSuperBasic&
		G_super = dyn_cast<const MatrixHessianSuperBasic>(G);
#endif

	// get some stuff
	const GenPermMatrixSlice
		&Q_R = G_super.Q_R(),
		&Q_X = G_super.Q_X();
	const size_type
		nd   = G_super.rows(),
		nd_R = Q_R.cols(),
		nd_X = Q_X.cols();
	assert( nd_R + nd_X == nd );

	// Setup output arguments

	// i_x_free[l-1] = (G.Q_R.begin()+l-1)->row_i(), l = 1...nd_R
	i_x_free->resize(nd_R);
	if(nd_R) {
		GenPermMatrixSlice::const_iterator
			Q_itr = Q_R.begin();
		i_x_free_t::iterator
			i_itr = i_x_free->begin(); 
		for( ; Q_itr != Q_R.end(); ++Q_itr, ++i_itr ) {
			const size_type i = Q_itr->row_i();
			assert( 0 < i && i <= nd );
			*i_itr = i;
		}
	}
	// i_x_fixed[]
	i_x_fixed->resize(nd_X+1);
	if(nd_X) {
		// i_x_fixed[l-1] = (G.Q_X.begin()+l-1)->row_i(), l = 1...nd_X
		GenPermMatrixSlice::const_iterator
			Q_itr = Q_X.begin();
		i_x_fixed_t::iterator
			i_itr = i_x_fixed->begin(); 
		for( ; Q_itr != Q_X.end(); ++Q_itr, ++i_itr ) {
			const size_type i = Q_itr->row_i();
			assert( 0 < i && i <= nd );
			*i_itr = i;
		}
	}
	(*i_x_fixed)[nd_X] = nd+1; // relaxation is always initially active
	// bnd_fixed[]
	bnd_fixed->resize(nd_X+1);
	if(nd_X) {
		// bnd_fixed[l-1] = G.bnd_fixed[l-1], l = 1...nd_X
		typedef MatrixHessianSuperBasic MHSB;
		const MHSB::bnd_fixed_t
			&bnd_fixed_from = G_super.bnd_fixed();
		assert(bnd_fixed_from.size() == nd_X);
		MHSB::bnd_fixed_t::const_iterator
			bnd_from_itr = bnd_fixed_from.begin();
		bnd_fixed_t::iterator
			bnd_to_itr = bnd_fixed->begin();
		for( ; bnd_from_itr != bnd_fixed_from.end(); ++bnd_from_itr, ++ bnd_to_itr ) {
			switch( *bnd_from_itr ) {
			    case MHSB::LOWER:
					*bnd_to_itr = QPSchurPack::LOWER;
					break;
			    case MHSB::UPPER:
					*bnd_to_itr = QPSchurPack::UPPER;
					break;
			    case MHSB::EQUALITY:
					*bnd_to_itr = QPSchurPack::EQUALITY;
					break;
			    default:
					assert(0);
			}
		}
	}
	(*bnd_fixed)[nd_X] = QPSchurPack::LOWER; // relaxation is always initially active
	// j_f_decomp[]
	j_f_decomp->resize(0);
	// b_X
	b_X->resize(nd_X+1);
	if(nd_X) {
		// b_X[l-1] = { dL(i) if bnd_fixed[l-1] == LOWER or EQUALITY
		//              dU(i) if bnd_fixed[l-1] == UPPER }
		//             , l = 1...nd_X
		//             (where i = i_x_fixed[l-1])
		bnd_fixed_t::const_iterator
			bnd_itr     = const_cast<const bnd_fixed_t&>(*bnd_fixed).begin(),
			bnd_itr_end = const_cast<const bnd_fixed_t&>(*bnd_fixed).begin();
		i_x_fixed_t::const_iterator
			i_x_itr     = const_cast<const i_x_fixed_t&>(*i_x_fixed).begin();
		Vector::iterator
			b_X_itr     = b_X->begin();
		const SpVectorSlice::element_type
			*ele = NULL;
		for( ; bnd_itr != bnd_itr_end; ++bnd_itr, ++i_x_itr, ++b_X_itr ) {
			const size_type i = *i_x_itr;
			switch(*bnd_itr) {
			    case QPSchurPack::LOWER:
			    case QPSchurPack::EQUALITY:
					*b_X_itr = (ele = dL.lookup_element(i))->value(); // Should not be null!
					break;
			    case QPSchurPack::UPPER:
					*b_X_itr = (ele = dU.lookup_element(i))->value(); // Should not be null!
					break;
			    default:
					assert(0);
			}
		}
	}
	(*b_X)[nd_X] = etaL; // relaxation is always initially active
	// Ko = G.B_RR
	*Ko = G_super.B_RR_ptr(); // now B_RR is a shared object
	// fo = - G.Q_R'*g - op(G.B_RX)*b_X(1:nd_X)
	V_StMtV( fo, -1.0, Q_R, trans, g );
	if( nd_X && G_super.B_RX_ptr().get() )
		Vp_StMtV( &(*fo)(), -1.0, *G_super.B_RX_ptr(), G_super.B_RX_trans(), (*b_X)(1,nd_X) );

}

} // end namesapce ConstrainedOptimizationPack
