// //////////////////////////////////////////////////////////////////////////
// QuasiRangeSpaceStepStd_Strategy.cpp

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
