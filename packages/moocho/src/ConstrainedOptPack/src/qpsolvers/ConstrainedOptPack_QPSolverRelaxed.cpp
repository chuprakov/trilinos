// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxed.cpp

#include <assert.h>

#include "ConstrainedOptimizationPack/include/QPSolverRelaxed.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "LinAlgPack/include/VectorClass.h"

namespace ConstrainedOptimizationPack {

// public

QPSolverStats::ESolutionType
QPSolverRelaxed::solve_qp(
	  std::ostream* out, EOutputLevel olevel, ERunTests test_what
	, const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& dL, const SpVectorSlice& dU
	, const MatrixWithOp& E, BLAS_Cpp::Transp trans_E, const VectorSlice& b
		, const SpVectorSlice& eL, const SpVectorSlice& eU
	, const MatrixWithOp& F, BLAS_Cpp::Transp trans_F, const VectorSlice& f
	, value_type* obj_d
	, value_type* eta, VectorSlice* d
	, SpVector* nu
	, SpVector* mu, VectorSlice* Ed
	, VectorSlice* lambda, VectorSlice* Fd
	)
{
	return solve_qp(out,olevel,test_what,g,G,etaL,dL,dU
		,&E,trans_E,&b,&eL,&eU,&F,trans_F,&f
		,obj_d,eta,d,nu,mu,Ed,lambda,Fd);
}

QPSolverStats::ESolutionType
QPSolverRelaxed::solve_qp(
	  std::ostream* out, EOutputLevel olevel, ERunTests test_what
	, const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& dL, const SpVectorSlice& dU
	, const MatrixWithOp& E, BLAS_Cpp::Transp trans_E, const VectorSlice& b
		, const SpVectorSlice& eL, const SpVectorSlice& eU
	, value_type* obj_d
	, value_type* eta, VectorSlice* d
	, SpVector* nu
	, SpVector* mu, VectorSlice* Ed
	)
{
	return solve_qp(out,olevel,test_what,g,G,etaL,dL,dU
		,&E,trans_E,&b,&eL,&eU,NULL,BLAS_Cpp::no_trans,NULL
		,obj_d,eta,d,nu,mu,Ed,NULL,NULL);
}

QPSolverStats::ESolutionType
QPSolverRelaxed::solve_qp(
	  std::ostream* out, EOutputLevel olevel, ERunTests test_what
	, const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& dL, const SpVectorSlice& dU
	, const MatrixWithOp& F, BLAS_Cpp::Transp trans_F, const VectorSlice& f
	, value_type* obj_d
	, value_type* eta, VectorSlice* d
	, SpVector* nu
	, VectorSlice* lambda, VectorSlice* Fd
	)
{
	return solve_qp(out,olevel,test_what,g,G,etaL,dL,dU
		,NULL,BLAS_Cpp::no_trans,NULL,NULL,NULL,&F,trans_F,&f
		,obj_d,eta,d,nu,NULL,NULL,lambda,Fd);
}


QPSolverStats::ESolutionType
QPSolverRelaxed::solve_qp(
	  std::ostream* out, EOutputLevel olevel, ERunTests test_what
	, const VectorSlice& g, const MatrixWithOp& G
	, const SpVectorSlice& dL, const SpVectorSlice& dU
	, value_type* obj_d
	, VectorSlice* d
	, SpVector* nu
	)
{
	return solve_qp(out,olevel,test_what,g,G,0,dL,dU
		,NULL,BLAS_Cpp::no_trans,NULL,NULL,NULL
		,NULL,BLAS_Cpp::no_trans,NULL
		,obj_d,NULL,d,nu,NULL,NULL,NULL,NULL);
}

QPSolverStats::ESolutionType
QPSolverRelaxed::solve_qp(
	  std::ostream* out, EOutputLevel olevel, ERunTests test_what
	, const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& dL, const SpVectorSlice& dU
	, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
		, const SpVectorSlice* eL, const SpVectorSlice* eU
	, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
	, value_type* obj_d
	, value_type* eta, VectorSlice* d
	, SpVector* nu
	, SpVector* mu, VectorSlice* Ed
	, VectorSlice* lambda, VectorSlice* Fd
	)
{

	validate_input(g,G,etaL,dL,dU
		,E,trans_E,b,eL,eU,F,trans_F,f
		,obj_d,eta,d,nu,mu,Ed,lambda,Fd);
	
	return imp_solve_qp(out,olevel,test_what,g,G,etaL,dL,dU
		,E,trans_E,b,eL,eU,F,trans_F,f
		,obj_d,eta,d,nu,mu,Ed,lambda,Fd);
}

void QPSolverRelaxed::validate_input(
	  const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& dL, const SpVectorSlice& dU
	, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
		, const SpVectorSlice* eL, const SpVectorSlice* eU
	, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
	, const value_type* obj_d
	, const value_type* eta, const VectorSlice* d
	, const SpVector* nu
	, const SpVector* mu, const VectorSlice* Ed
	, const VectorSlice* lambda, const VectorSlice* Fd
	)
{
	// Validate output arguments
	if( !d )
		throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
			"If d!=NULL is not allowed." );
	if( ( E || F ) && !eta )
		throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
			"If eta!=NULL is not allowed if E!=NULL or F!=NULL." );

	// Validate the sets of constraints arguments
	if( E && ( !b || !eL || !eU || !mu ) )
		throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
			"If E!=NULL then b!=NULL, eL!=NULL, eU!=NULL and mu!=NULL must also "
			"be true." );
	if( F && ( !f || !lambda ) )
		throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
			"If F!=NULL then f!=NULL and lambda!=NULL must also "
			"be true." );

	// Validate the sizes of the arguments
	const size_type
		n = d->size();
	if( g.size() != n )
		throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
			"g.size() != d->size()." );
	if( G.rows() != n || G.cols() != n )
		throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
			"G.rows() != d->size() or G.cols() != d->size()." );
	if( dL.size() != n )
		throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
			"dL.size() != d->size()." );
	if( dU.size() != n )
		throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
			"dU.size() != d->size()." );
	if( E ) {
		const size_type
			m_in = BLAS_Cpp::rows( E->rows(), E->cols(), trans_E );
		if( BLAS_Cpp::cols( E->rows(), E->cols(), trans_E )	!= n )
			throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
				"op(E).cols() != d->size()." );
		if( b->size() != m_in )
			throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
				"b->size() != op(E).rows()." );
		if( eL->size() != m_in )
			throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
				"eL->size() != op(E).rows()." );
		if( eU->size() != m_in )
			throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
				"eU->size() != op(E).rows()." );
		if( Ed && Ed->size() != m_in )
			throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
				"Ed->size() != op(E).rows()." );
	}
	if( F ) {
		const size_type
			m_eq = BLAS_Cpp::rows( F->rows(), F->cols(), trans_F );
		if( BLAS_Cpp::cols( F->rows(), F->cols(), trans_F )	!= n )
			throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
				"op(F).cols() != d->size()." );
		if( f->size() != m_eq )
			throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
				"f->size() != op(F).rows()." );
		if( lambda->size() != m_eq )
			throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
				"lambda->size() != op(F).rows()." );
		if( Fd && Fd->size() != m_eq )
			throw std::invalid_argument( "QPSolverRelaxed::validate_input(...) : Error, "
				"Fd->size() != op(F).rows()." );
	}
}

}	// end namespace ConstrainedOptimizationPack
