// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxed.cpp

#include <assert.h>

#include "ConstrainedOptimizationPack/include/QPSolverRelaxed.h"

namespace ConstrainedOptimizationPack {

// public

QPSolverStats::ESolutionType
QPSolverRelaxed::solve_qp(
	  const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& xL, const SpVectorSlice& xU
	, const MatrixWithOp& E, BLAS_Cpp::Transp trans_E, const VectorSlice& b
		, const SpVectorSlice& eL, const SpVectorSlice& eU
	, const MatrixWithOp& F, BLAS_Cpp::Transp trans_F, const VectorSlice& f
	, std::ostream* out, EOutputLevel olevel, ERunTests test_what
	, value_type* eta, Vector* x, SpVector* nu, SpVector* mu, Vector* lambda
	)
{
	return solve_qp(g,G,etaL,xL,xU,&E,trans_E,&b,&eL,&eU,&F,trans_F,&f
		,out,olevel,test_what,eta,x,nu,mu,lambda);
}

QPSolverStats::ESolutionType
QPSolverRelaxed::solve_qp(
	  const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& xL, const SpVectorSlice& xU
	, const MatrixWithOp& E, BLAS_Cpp::Transp trans_E, const VectorSlice& b
		, const SpVectorSlice& eL, const SpVectorSlice& eU
	, std::ostream* out, EOutputLevel olevel, ERunTests test_what
	, value_type* eta, Vector* x, SpVector* nu, SpVector* mu
	)
{
	return solve_qp(g,G,etaL,xL,xU,&E,trans_E,&b,&eL,&eU,NULL,BLAS_Cpp::no_trans,NULL
		,out,olevel,test_what,eta,x,nu,mu,NULL);
}

QPSolverStats::ESolutionType
QPSolverRelaxed::solve_qp(
	  const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& xL, const SpVectorSlice& xU
	, const MatrixWithOp& F, BLAS_Cpp::Transp trans_F, const VectorSlice& f
	, std::ostream* out, EOutputLevel olevel, ERunTests test_what
	, value_type* eta, Vector* x, SpVector* nu, Vector* lambda
	)
{
	return solve_qp(g,G,etaL,xL,xU,NULL,BLAS_Cpp::no_trans,NULL,NULL,NULL,&F,trans_F,&f
		,out,olevel,test_what,eta,x,nu,NULL,lambda);
}


QPSolverStats::ESolutionType
QPSolverRelaxed::solve_qp(
	  const VectorSlice& g, const MatrixWithOp& G
	, const SpVectorSlice& xL, const SpVectorSlice& xU
	, std::ostream* out, EOutputLevel olevel, ERunTests test_what
	, Vector* x, SpVector* nu
	)
{
	return solve_qp(g,G,0,xL,xU,NULL,BLAS_Cpp::no_trans,NULL,NULL,NULL
		,NULL,BLAS_Cpp::no_trans,NULL
		,out,olevel,test_what,NULL,x,nu,NULL,NULL);
}

// protected

QPSolverStats::ESolutionType
QPSolverRelaxed::solve_qp(
	  const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& xL, const SpVectorSlice& xU
	, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
		, const SpVectorSlice* eL, const SpVectorSlice* eU
	, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
	, std::ostream* out, EOutputLevel olevel, ERunTests test_what
	, value_type* eta, Vector* x, SpVector* nu, SpVector* mu, Vector* lambda
	)
{

	// ToDo: Implement checking sets of present constraints and sizes
	assert(0);

	return solve_qp(g,G,etaL,xL,xU,E,trans_E,b,eL,eU,F,trans_F,f
		,out,olevel,test_what,eta,x,nu,mu,lambda);
}
	
}	// end namespace ConstrainedOptimizationPack
