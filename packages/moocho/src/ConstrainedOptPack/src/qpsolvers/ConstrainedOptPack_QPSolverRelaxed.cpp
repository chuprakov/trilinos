// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxed.cpp

#include <assert.h>

#include "ConstrainedOptimizationPack/include/QPSolverRelaxed.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/MatrixWithOpOut.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/SpVectorOut.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "Misc/include/profile_hack.h"

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
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "QPSolverRelaxed::solve_qp(...)" );
#endif

	validate_input(g,G,etaL,dL,dU
		,E,trans_E,b,eL,eU,F,trans_F,f
		,obj_d,eta,d,nu,mu,Ed,lambda,Fd);

	print_qp_input(
		out,olevel,g,G,etaL,dL,dU,E,trans_E,b,eL,eU
		,F,trans_F,f,eta,d,nu,mu,lambda	);
	
	QPSolverStats::ESolutionType
		solve_return = imp_solve_qp(
			out,olevel,test_what,g,G,etaL,dL,dU
			,E,trans_E,b,eL,eU,F,trans_F,f
			,obj_d,eta,d,nu,mu,Ed,lambda,Fd);

	print_qp_output(out,olevel,obj_d,eta,d,nu,mu,Ed,lambda,Fd);

	return solve_return;

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

void QPSolverRelaxed::print_qp_input( 
	std::ostream* out, EOutputLevel olevel
	, const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& dL, const SpVectorSlice& dU
	, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
	, const SpVectorSlice* eL, const SpVectorSlice* eU
	, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
	, value_type* eta, VectorSlice* d
	, SpVector* nu
	, SpVector* mu
	, VectorSlice* lambda
	)
{
	using LinAlgPack::norm_inf;
	using SparseLinAlgPack::norm_inf;

	if( out && (int)olevel >= (int)PRINT_ITER_STEPS ) {
		*out<< "\n*** Printing input to QPSolverRelaxed::solve_qp(...) ...\n";
		// g
		*out << "\n||g||inf = " << norm_inf(g) << std::endl;
		if( (int)olevel >= (int)PRINT_ITER_VECTORS )
			*out<< "g =\n" << g;
		// G
		if( (int)olevel >= (int)PRINT_EVERY_THING )
			*out<< "\nG =\n" << G;
		// etaL
		*out << "\netaL = " << etaL << std::endl;
		// eta
		*out << "\neta = " << *eta << std::endl;
		// dL
		*out<< "\ndL.nz()   = " << dL.nz() << std::endl
			<< "||dL||inf = " << norm_inf(dL) << std::endl;
		if( (int)olevel >= (int)PRINT_ITER_VECTORS )
			*out<< "dL =\n" << dL;
		// dU
		*out<< "\ndU.nz()   = " << dU.nz() << std::endl
			<< "||dU||inf = " << norm_inf(dU) << std::endl;
		if( (int)olevel >= (int)PRINT_ITER_VECTORS )
			*out<< "dU =\n" << dU;
		// d
		*out << "\n||d||inf = " << norm_inf(*d) << std::endl;
		if( (int)olevel >= (int)PRINT_ITER_VECTORS )
			*out<< "d =\n" << *d;
		// nu
		if(nu) {
			*out<< "\nnu.nz()   = " << nu->nz() << std::endl
				<< "||nu||inf = " << norm_inf((*nu)()) << std::endl;
			if( (int)olevel >= (int)PRINT_ITER_ACT_SET )
				*out<< "nu =\n" << *nu;
		}
		if(E) {
			// op(E)
			if( (int)olevel >= (int)PRINT_EVERY_THING )
				*out<< "\nE" << std::endl << *E
					<< "trans_E = " << BLAS_Cpp::trans_to_string(trans_E) << std::endl;
			// b
			*out << "\n||b||inf = " << norm_inf(*b) << std::endl;
			if( (int)olevel >= (int)PRINT_ITER_VECTORS )
			*out<< "b =\n" << *b;
			// eL
			*out<< "\neL.nz()   = " << eL->nz() << std::endl;
			if( (int)olevel >= (int)PRINT_ITER_VECTORS )
				*out<< "eL =\n" << *eL;
			// eU
			*out<< "\neU.nz()   = " << eU->nz() << std::endl;
			if( (int)olevel >= (int)PRINT_ITER_VECTORS )
				*out<< "eU =\n" << *eU;
			// mu
			*out<< "\nmu.nz()   = " << mu->nz() << std::endl
				<< "||mu||inf = " << norm_inf((*mu)()) << std::endl;
			if( (int)olevel >= (int)PRINT_ITER_ACT_SET )
				*out<< "mu =\n" << *mu;
		}
		if(F) {
			// op(F)
			if( (int)olevel >= (int)PRINT_EVERY_THING )
				*out<< "\nF" << std::endl << *F
					<< "trans_F = " << BLAS_Cpp::trans_to_string(trans_F) << std::endl;
			// f
			*out<< "\n||f||inf = " << norm_inf(*f) << std::endl;
			if( (int)olevel >= (int)PRINT_ITER_VECTORS )
				*out<< "f =\n" << *f;
			// lambda
			*out<< "\n||lambda||inf = " << norm_inf(*lambda) << std::endl;
			if( (int)olevel >= (int)PRINT_ITER_ACT_SET )
				*out<< "lambda =\n" << *lambda;
		}
	}
}

void QPSolverRelaxed::print_qp_output(
	std::ostream* out, EOutputLevel olevel
	, const value_type* obj_d
	, const value_type* eta, const VectorSlice* d
	, const SpVector* nu
	, const SpVector* mu, const VectorSlice* Ed
	, const VectorSlice* lambda, const VectorSlice* Fd
	)
{
	using LinAlgPack::norm_inf;
	using SparseLinAlgPack::norm_inf;

	if( out && (int)olevel > (int)PRINT_ITER_STEPS ) {
		*out<< "\n*** Printing output from QPSolverRelaxed::solve_qp(...) ...\n";
		// obj_d
		if(obj_d)
			*out << "\nobj_d = " << *obj_d << std::endl;
		// eta
		*out << "\neta = " << *eta << std::endl;
		// d
		*out << "\n||d||inf = " << norm_inf(*d) << std::endl;
		if( (int)olevel >= (int)PRINT_ITER_VECTORS )
			*out<< "d =\n" << *d;
		// nu
		if(nu) {
			*out<< "\nnu.nz()   = " << nu->nz() << std::endl
				<< "||nu||inf = " << norm_inf((*nu)()) << std::endl;
			if( (int)olevel >= (int)PRINT_ITER_ACT_SET )
				*out<< "nu =\n" << *nu;
		}
		// Ed
		if(Ed) {
			*out << "\n||Ed||inf = " << norm_inf(*Ed) << std::endl;
			if( (int)olevel >= (int)PRINT_ITER_VECTORS )
				*out<< "Ed =\n" << *Ed;
		}
		// mu
		if(mu) {
			*out<< "\nmu.nz()   = " << mu->nz() << std::endl
				<< "||mu||inf = " << norm_inf((*mu)()) << std::endl;
			if( (int)olevel >= (int)PRINT_ITER_ACT_SET )
				*out<< "mu =\n" << *mu;
		}
		// lambda
		if(lambda) {
			*out<< "\n||lambda||inf = " << norm_inf(*lambda) << std::endl;
			if( (int)olevel >= (int)PRINT_ITER_ACT_SET )
				*out<< "lambda =\n" << *lambda;
		}
		// Fd
		if(Fd) {
			*out << "\n||Fd||inf = " << norm_inf(*Fd) << std::endl;
			if( (int)olevel >= (int)PRINT_ITER_VECTORS )
				*out<< "Fd =\n" << *Fd;
		}
	}
}

}	// end namespace ConstrainedOptimizationPack
