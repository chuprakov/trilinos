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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include "ReducedSpaceSQPPack/include/std/QuasiRangeSpaceStepStd_Strategy.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/WorkspacePack.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

namespace ReducedSpaceSQPPack {

bool QuasiRangeSpaceStepStd_Strategy::solve_quasi_range_space_step(
	std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
	,const VectorSlice& xo, const VectorSlice& c_xo, VectorSlice* v
  	)
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	wsp::Workspace<value_type>  vy_ws(wss,c_xo.size());
	VectorSlice                 vy(&vy_ws[0],vy_ws.size());
	// vy = inv(Gc_k'*Y_k) * c_xo
	s->decomp_sys().solve_transAtY( c_xo, BLAS_Cpp::no_trans, &vy );
	// v = -Y_k*vy
	LinAlgOpPack::V_StMtV( v, -1.0, s->Y().get_k(0), BLAS_Cpp::no_trans, vy );

	return true;
}

void QuasiRangeSpaceStepStd_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
	out << L << "*** Compute the approximate range space step:\n"
		<< L << "vy = inv(Gc_k'*Y_k) * c_xo\n"
		<< L << "v = -Y_k*vy\n";
}

} // end namespace ReducedSpaceSQPPack
