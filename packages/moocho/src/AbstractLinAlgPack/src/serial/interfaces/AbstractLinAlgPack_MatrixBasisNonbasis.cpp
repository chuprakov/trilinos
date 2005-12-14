// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <assert.h>

#include "AbstractLinAlgPack_MatrixBasisNonbasis.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_DMatrixOp.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

namespace AbstractLinAlgPack {

// Overridden from MatrixOp

void MatrixBasisNonbasis::Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const
{
	using AbstractLinAlgPack::Vp_StMtV;

	size_type n = rows(), m = cols();

	DenseLinAlgPack::Vp_MtV_assert_sizes( vs_lhs->size()
		, n, m, trans_rhs1, vs_rhs2.size() );

	switch(trans_rhs1) {
		case BLAS_Cpp::no_trans:
			//
			// M = [ C' ; N' ]
			//
			// y = b*y + a * M * x
			//
			// y = b*y + a * [ C' ] * x
			//               [ N' ]
			//
			// y(1:m)   = b*y(1:m)   + a * C' * x
			// y(m+1:n) = b*y(m+1:n) + a * N' * x
			
			// y(1:m)   = b*y(1:m)   + a * C' * x
			Vp_StMtV( &(*vs_lhs)(1,m), alpha, C(), BLAS_Cpp::trans, vs_rhs2, beta );

			// y(m+1:n) = b*y(m+1:n) + a * N' * x
			if( rows() > cols() )
				Vp_StMtV( &(*vs_lhs)(m+1,n), alpha, N(), BLAS_Cpp::trans, vs_rhs2, beta );

			break;

		case BLAS_Cpp::trans:
			//
			// M' = [ C , N ]
			//
			// y = b*y + a * M' * x
			//
			// y = b*y + a * [ C , N ] * [ x(1:m)   ]
			//                           [ x(m+1:n) ]
			//
			// y = b*y + a * C * x(1:m) + a * N * x(m+1:n)

			// y = b*y + a * C * x(1:m)
			Vp_StMtV( vs_lhs, alpha, C(), BLAS_Cpp::no_trans, vs_rhs2(1,m), beta );
			
			// y = y + a * N * x(m+1:n)
			if( rows() > cols() )
				Vp_StMtV( vs_lhs, alpha, N(), BLAS_Cpp::no_trans, vs_rhs2(m+1,n) );

			break;

		default:
			assert(0);
	}
}

void MatrixBasisNonbasis::Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const
{
	using AbstractLinAlgPack::Vp_StMtV;

	size_type n = rows(), m = cols();
	// Assert that the dimensions of the aruments match up and if not
	// then throw an excption.
	DenseLinAlgPack::Vp_MtV_assert_sizes( vs_lhs->size()
		, n, m, trans_rhs1, sv_rhs2.size() );

	switch(trans_rhs1) {
		case BLAS_Cpp::no_trans:
			//
			// M = [ C' ; N' ]
			//
			// y = b*y + a * M * x
			//
			// y = b*y + a * [ C' ] * x
			//               [ N' ]
			//
			// y(1:m)   = b*y(1:m)   + a * C' * x
			// y(m+1:n) = b*y(m+1:n) + a * N' * x
			
			// y(1:m)   = b*y(1:m)   + a * C' * x
			Vp_StMtV( &(*vs_lhs)(1,m), alpha, C(), BLAS_Cpp::trans, sv_rhs2, beta );

			// y(m+1:n) = b*y(m+1:n) + a * N' * x
			if( rows() > cols() )
				Vp_StMtV( &(*vs_lhs)(m+1,n), alpha, N(), BLAS_Cpp::trans, sv_rhs2, beta );

			break;

		case BLAS_Cpp::trans:
			//
			// M' = [ C , N ]
			//
			// y = b*y + a * M' * x
			//
			// y = b*y + a * [ C , N ] * [ x(1:m)   ]
			//                           [ x(m+1:n) ]
			//
			// y = b*y + a * C * x(1:m) + a * N * x(m+1:n)

			// y = b*y + a * C * x(1:m)
			Vp_StMtV( vs_lhs, alpha, C(), BLAS_Cpp::no_trans, sv_rhs2(1,m), beta );
			
			// y = y + a * N * x(m+1:n)
			if( rows() > cols() )
				Vp_StMtV( vs_lhs, alpha, N(), BLAS_Cpp::no_trans, sv_rhs2(m+1,n) );

			break;

		default:
			assert(0);
	}
}

}	// end namespace AbstractLinAlgPack 
