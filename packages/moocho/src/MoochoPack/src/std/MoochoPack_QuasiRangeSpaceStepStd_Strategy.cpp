// //////////////////////////////////////////////////////////////////////////
// QuasiRangeSpaceStepStd_Strategy.cpp
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

#include "ReducedSpaceSQPPack/include/std/QuasiRangeSpaceStepStd_Strategy.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"

namespace ReducedSpaceSQPPack {

bool QuasiRangeSpaceStepStd_Strategy::solve_quasi_range_space_step(
	std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
	,const VectorWithOp& xo, const VectorWithOp& c_xo, VectorWithOpMutable* v
  	)
{
	using LinAlgOpPack::V_InvMtV;
	using LinAlgOpPack::V_StMtV;
	const MatrixWithOpNonsingular
		&R_k = s->R().get_k(0);
	VectorSpace::vec_mut_ptr_t
		vy = R_k.space_rows().create_member();
	// vy = inv(R_k) * c_xo
	V_InvMtV( vy.get(), R_k, BLAS_Cpp::no_trans, c_xo );
	// v = -Y_k*vy
	V_StMtV( v, -1.0, s->Y().get_k(0), BLAS_Cpp::no_trans, *vy );

	return true;
}

void QuasiRangeSpaceStepStd_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
	out << L << "*** Compute the approximate range space step:\n"
		<< L << "vy = inv(R_k) * c_xo\n"
		<< L << "v = -Y_k*vy\n";
}

} // end namespace ReducedSpaceSQPPack
