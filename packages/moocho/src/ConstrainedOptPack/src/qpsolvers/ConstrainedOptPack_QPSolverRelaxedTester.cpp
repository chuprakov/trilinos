// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxedTester.cpp

#include "ConstrainedOptimizationPack/include/QPSolverRelaxedTester.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/test/CompareDenseVectors.h"
#include "SparseLinAlgPack/include/sparse_bounds_diff.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "LinAlgPack/include/LinAlgOpPack.h"

namespace {

//
const char* solution_type_str( ConstrainedOptimizationPack::QPSolverStats::ESolutionType solution_type )
{
	typedef ConstrainedOptimizationPack::QPSolverStats qpst;
	switch( solution_type ) {
	
	case qpst::OPTIMAL_SOLUTION:
		return "OPTIMAL_SOLUTION";
	case qpst::PRIMAL_FEASIBLE_POINT:
		return "PRIMAL_FEASIBLE_POINT";
	case qpst::DUAL_FEASIBLE_POINT:
		return "DUAL_FEASIBLE_POINT";
	case qpst::SUBOPTIMAL_POINT:
		return "SUBOPTIMAL_POINT";
	default:
		assert(0);
	}
	return "";	// will never be executed.
}

// Compute complementarity conditions.
// 
// If uplo == upper then:
// 
//                / gamma(i) * constr_resid(i), for gamma(i) > 0
// comp_err(i) = |
//                \ 0 otherwise
// 
// If uplo == lower then:
// 
//                / gamma(i) * constr_resid(i), for gamma(i) < 0
// comp_err(i) = |
//                \ 0 otherwise
// 
// 
void set_complementarity(
	  const SparseLinAlgPack::SpVector	&gamma
	, const LinAlgPack::VectorSlice		&constr_resid
	, BLAS_Cpp::Uplo					uplo
	, LinAlgPack::Vector			 	*comp_err
	)
{
	assert( gamma.size() == constr_resid.size() );
	comp_err->resize( gamma.size() );
	*comp_err = 0.0;
	const SparseLinAlgPack::SpVector::difference_type o = gamma.offset();
	if( gamma.nz() ) {
		for( SparseLinAlgPack::SpVector::const_iterator itr = gamma.begin(); itr != gamma.end(); ++itr ) {
			const LinAlgPack::size_type i = itr->indice() + o;
			if( itr->value() > 0 && uplo == BLAS_Cpp::upper )
				(*comp_err)(i) = itr->value() * constr_resid(i);
			else if( itr->value() < 0 && uplo == BLAS_Cpp::lower )
				(*comp_err)(i) = itr->value() * constr_resid(i);
		}
	}
}

}	// end namespace

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StV;
	using SparseLinAlgPack::Vp_StMtV;
}

namespace ConstrainedOptimizationPack {

// public

QPSolverRelaxedTester::QPSolverRelaxedTester(
		  value_type	opt_warning_tol
		, value_type	opt_error_tol
		, value_type	feas_warning_tol
		, value_type	feas_error_tol
		, value_type	comp_warning_tol
		, value_type	comp_error_tol
		)
	:
		opt_warning_tol_(opt_warning_tol)
		,opt_error_tol_(opt_error_tol)
		,feas_warning_tol_(feas_warning_tol)
		,feas_error_tol_(feas_error_tol)
		,comp_warning_tol_(comp_warning_tol)
		,comp_error_tol_(comp_error_tol)
{}

bool QPSolverRelaxedTester::check_optimality_conditions(
	  QPSolverStats::ESolutionType solution_type
	, std::ostream* out, bool print_all_warnings, bool print_vectors
	, const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& dL, const SpVectorSlice& dU
	, const MatrixWithOp& E, BLAS_Cpp::Transp trans_E, const VectorSlice& b
		, const SpVectorSlice& eL, const SpVectorSlice& eU
	, const MatrixWithOp& F, BLAS_Cpp::Transp trans_F, const VectorSlice& f
	, const value_type* obj_d
	, const value_type* eta, const VectorSlice* d
	, const SpVector* nu
	, const SpVector* mu, const VectorSlice* Ed
	, const VectorSlice* lambda, const VectorSlice* Fd
	)
{
	check_optimality_conditions(solution_type,out,print_all_warnings,print_vectors,g,G
		,etaL,dL,dU,&E,trans_E,&b,&eL,&eU,&F,trans_F,&f
		,obj_d,eta,d,nu,mu,Ed,lambda,Fd);
}

bool QPSolverRelaxedTester::check_optimality_conditions(
	  QPSolverStats::ESolutionType solution_type
	, std::ostream* out, bool print_all_warnings, bool print_vectors
	, const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& dL, const SpVectorSlice& dU
	, const MatrixWithOp& E, BLAS_Cpp::Transp trans_E, const VectorSlice& b
		, const SpVectorSlice& eL, const SpVectorSlice& eU
	, const value_type* obj_d
	, const value_type* eta, const VectorSlice* d
	, const SpVector* nu
	, const SpVector* mu, const VectorSlice* Ed
	)
{
	check_optimality_conditions(solution_type,out,print_all_warnings,print_vectors,g,G
		,etaL,dL,dU,&E,trans_E,&b,&eL,&eU,NULL,BLAS_Cpp::no_trans,NULL
		,obj_d,eta,d,nu,mu,Ed,NULL,NULL);
}

bool QPSolverRelaxedTester::check_optimality_conditions(
	  QPSolverStats::ESolutionType solution_type
	, std::ostream* out, bool print_all_warnings, bool print_vectors
	, const VectorSlice& g, const MatrixWithOp& G
	, value_type etaL
	, const SpVectorSlice& dL, const SpVectorSlice& dU
	, const MatrixWithOp& F, BLAS_Cpp::Transp trans_F, const VectorSlice& f
	, const value_type* obj_d
	, const value_type* eta, const VectorSlice* d
	, const SpVector* nu
	, const VectorSlice* lambda, const VectorSlice* Fd
	)
{
	check_optimality_conditions(solution_type,out,print_all_warnings,print_vectors,g,G
		,etaL,dL,dU,NULL,BLAS_Cpp::no_trans,NULL,NULL,NULL,&F,trans_F,&f
		,obj_d,eta,d,nu,NULL,NULL,lambda,Fd );
}

bool QPSolverRelaxedTester::check_optimality_conditions(
	  QPSolverStats::ESolutionType solution_type
	, std::ostream* out, bool print_all_warnings, bool print_vectors
	, const VectorSlice& g, const MatrixWithOp& G
	, const SpVectorSlice& dL, const SpVectorSlice& dU
	, const value_type* obj_d
	, const VectorSlice* d
	, const SpVector* nu
	)
{
	check_optimality_conditions(solution_type,out,print_all_warnings,print_vectors,g,G
		,0.0,dL,dU,NULL,BLAS_Cpp::no_trans,NULL,NULL,NULL,NULL,BLAS_Cpp::no_trans,NULL
		,obj_d,NULL,d,nu,NULL,NULL,NULL,NULL);
}

bool QPSolverRelaxedTester::check_optimality_conditions(
	  QPSolverStats::ESolutionType solution_type
	, std::ostream* out, bool print_all_warnings, bool print_vectors
	, const VectorSlice& g, const MatrixWithOp& G
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
	QPSolverRelaxed::validate_input(g,G,etaL,dL,dU
		,E,trans_E,b,eL,eU,F,trans_F,f
		,obj_d,eta,d,nu,mu,Ed,lambda,Fd);

	return imp_check_optimality_conditions(solution_type
		,out,print_all_warnings,print_vectors,g,G,etaL,dL,dU
		,E,trans_E,b,eL,eU,F,trans_F,f
		,obj_d,eta,d,nu,mu,Ed,lambda,Fd);
}

// protected

bool QPSolverRelaxedTester::imp_check_optimality_conditions(
	  QPSolverStats::ESolutionType solution_type
	, std::ostream* out, bool print_all_warnings, bool print_vectors
	, const VectorSlice& g, const MatrixWithOp& G
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
	using std::endl;
	using BLAS_Cpp::trans_not;
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::upper;
	using BLAS_Cpp::lower;
	using LinAlgPack::norm_inf;
	using LinAlgPack::Vt_S;
	using LinAlgPack::Vp_StV;
	using LinAlgOpPack::V_MtV;
	using LinAlgOpPack::Vp_MtV;
	using LinAlgOpPack::V_StMtV;
	using LinAlgOpPack::Vp_V;
	using SparseLinAlgPack::norm_inf;
	using SparseLinAlgPack::imp_sparse_bnd_diff;
	typedef QPSolverStats qps_t;

	const value_type really_big_error_tol = std::numeric_limits<value_type>::max();

	value_type scale = 0.0;
	Vector
		u,	// hold the result to pass to comparison function
		c,	// hold complementarity conditions
		t,	// hold temporary vectors
		r,	// - op(F)*d + eta * f
		e;	// op(E)*d + b*eta
	value_type
		d_norm_inf,	// holds ||d||inf
		e_norm_inf;	// holds ||e||inf

	SparseLinAlgPack::TestingPack::CompareDenseVectors comp_v;

	if(out)
		*out
			<< "\n*** Begin checking QP optimality conditions\n"
			<< "\nThe solution type is " << solution_type_str(solution_type) << endl;

	bool force_opt_error_check
		= solution_type==qps_t::OPTIMAL_SOLUTION || solution_type==qps_t::DUAL_FEASIBLE_POINT;
	const bool force_inequality_error_check
		= solution_type==qps_t::OPTIMAL_SOLUTION || solution_type==qps_t::PRIMAL_FEASIBLE_POINT;
	const bool force_equality_error_check
		= solution_type!=qps_t::SUBOPTIMAL_POINT;
	const bool force_complementarity_error_check
		= solution_type!=qps_t::SUBOPTIMAL_POINT;

	// Check some postconditions
	if(nu) nu->assert_valid_and_sorted();
	if(mu) mu->assert_valid_and_sorted();

	////////////////////////////
	// Checking d(L)/d(d) = 0
	if(out)
		*out
			<< "\nChecking d(L)/d(d) = g + G*d + nu + op(E)'*mu - op(F)'*lambda == 0 ...\n";

	if(out && !force_opt_error_check)
		*out
			<< "The optimality error tolerance will not be enforced ...\n";

	scale = 1.0;

	u = g;
	scale += norm_inf(g);
	
	V_MtV( &t, G, no_trans, *d );
	Vp_V( &u(), t() );
	scale += norm_inf(t);

	if(out) {
		*out << "||g + G*d||inf = " << norm_inf(u()) << endl;
		if(print_vectors)
			*out << "g + G*d =\n" << u();
	}

	if( nu ) {
		Vp_V( &u(), (*nu)() );
		scale += norm_inf((*nu)());
	}

	if(E) {
		V_MtV( &t, *E, trans_not(trans_E), (*mu)() );
		Vp_V( &u(), t() );
		scale += norm_inf(t);
		if(out) {
			*out << "||op(E)'*mu||inf = " << norm_inf(t()) << endl;
			if(print_vectors)
				*out << "op(E)'*mu =\n" << t();
		}
	}			

	if(F) {
		V_MtV( &t, *F, trans_not(trans_F), *lambda );
		Vp_V( &u(), t() );
		scale += norm_inf(t);
		if(out && print_vectors)
			*out
				<< "\nop(F)'*lambda =\n" << t();
		if(out) {
			*out << "||op(F)'*lambda||inf = " << norm_inf(t()) << endl;
			if(print_vectors)
				*out << "op(F)'*lambda =\n" << t();
		}
	}			

	if(out && print_vectors)
		*out
			<< "g + G*d + nu + op(E)'*mu - op(F)'*lambda =\n" << u();

	Vt_S( &u(), 1.0/scale );

	if(out) {
		*out
			<< "Comparing:\n"
			<< "u = | g + G*d + nu + op(E)'*mu - op(F)'*lambda | / scale\n"
			<< "scale = " << scale << endl
			<< "v = 0 ...\n";
	}

	if(!comp_v.comp( u(), 0.0, opt_warning_tol()
		, force_opt_error_check ? opt_error_tol() : really_big_error_tol
		, print_all_warnings, out ))	return false;

	if(out) {
		*out
			<< "\nTesting feasibility of the constraints and the complementarity conditions ...\n";
		if(!force_inequality_error_check)
			*out
				<< "The inequality feasibility error tolerance will not be enforced ...\n";
		if(!force_equality_error_check)
			*out
				<< "The equality feasibility error tolerance will not be enforced ...\n";
		if(!force_complementarity_error_check)
			*out
				<< "The complementarity conditions error tolerance will not be enforced ...\n";
	}

	/////////////////////
	// etaL - eta
	if(out)
		*out
			<< "\nChecking etaL - eta <= 0 ...\n";
	if(out && print_vectors)
		*out
			<< "etaL - eta = " << (etaL - (*eta)) << endl;
	if( etaL - (*eta) > feas_warning_tol() ) {
		if(out)
			*out
				<< "Warning, etaL - eta = " << etaL << " - " << (*eta)
				<< " = " << (etaL - (*eta)) << " >  feas_warning_tol = "
				<<  feas_warning_tol() << endl;
	}
	if( force_inequality_error_check && etaL - (*eta) > feas_error_tol() ) {
		if(out)
			*out
				<< "Error, etaL - eta = " << etaL << " - " << (*eta)
				<< " = " << (etaL - (*eta)) << " >  feas_error_tol = "
				<<  feas_error_tol() << endl;
		return false;
	} 

	d_norm_inf = norm_inf(*d);

	///////////////////////////////////
	// dL - d <= 0
	if(out)
		*out
			<< "\nChecking dL - d <= 0 ...\n";
	u.resize(d->size());
	imp_sparse_bnd_diff( +1, dL, lower, *d, &u() );
	if(out && print_vectors)
		*out
			<< "dL - d =\n" << u();
	Vt_S( &u(), 1.0/(1.0+d_norm_inf) );
	if(out) {
		*out
			<< "Comparing: u - v <= 0\n"
			<< "u = (dL - d) | / (1 + ||d||inf ), v = 0 ...\n";
	}
	if(!comp_v.comp_less( u(), 0.0, opt_warning_tol()
		, force_inequality_error_check ? feas_error_tol() : really_big_error_tol
		, print_all_warnings, out )) return false;
	if(nu) {
		if(out)
			*out
				<< "\nChecking nuL(i) * (dL - d)(i) = 0  ...\n";
		set_complementarity( *nu, u(), lower, &c );
		if(out && print_vectors)
			*out
				<< "nuL(i) * (dL - d)(i) =\n" << c();
		if(out) {
			*out
				<< "Comparing:\n"
				<< "u(i) = nuL(i) * (dL - d)(i), v = 0 ...\n";
		}
		if(!comp_v.comp( c(), 0.0, opt_warning_tol()
			, force_complementarity_error_check ? comp_error_tol() : really_big_error_tol
			, print_all_warnings, out )) return false;
	}

	///////////////////////////////////
	// d - dU <= 0
	if(out)
		*out
			<< "\nChecking d - dU <= 0 ...\n";
	u.resize(d->size());
	imp_sparse_bnd_diff( -1, dU, upper, *d, &u() );
	if(out && print_vectors)
		*out
			<< "d - dU =\n" << u();
	Vt_S( &u(), 1.0/(1.0+d_norm_inf) );
	if(out) {
		*out
			<< "Comparing: u - v <= 0\n"
			<< "u = (d - dU) | / (1 + ||d||inf ), v = 0 ...\n";
	}
	if(!comp_v.comp_less( u(), 0.0, opt_warning_tol()
		, force_inequality_error_check ? feas_error_tol() : really_big_error_tol
		, print_all_warnings, out )) return false;
	if(nu) {
		if(out)
			*out
				<< "\nChecking nuU(i) * (d - dU)(i) = 0  ...\n";
		set_complementarity( *nu, u(), upper, &c );
		if(out && print_vectors)
			*out
				<< "nuU(i) * (d - dU)(i) =\n" << c();
		if(out) {
			*out
				<< "Comparing:\n"
				<< "u(i) = nuU(i) * (dL - d)(i), v = 0 ...\n";
		}
		if(!comp_v.comp( c(), 0.0, opt_warning_tol()
			, force_complementarity_error_check ? comp_error_tol() : really_big_error_tol
			, print_all_warnings, out )) return false;
	}

	if( E ) {

		///////////////////////////////////
		// e = op(E)*d + b*eta
		if(out)
			*out
				<< "\nComputing e = op(E)*d + b*eta ...\n";
		V_MtV( &e, *E, trans_E, *d );
		Vp_StV( &e(), *eta, *b );
		e_norm_inf = norm_inf(e());	
		if(out && print_vectors)
			*out
				<< "e = op(E)*d + b*eta  =\n" << e();

		///////////////////////////////////
		// eL - e <= 0
		if(out)
			*out
				<< "\nChecking eL - e <= 0 ...\n";
		u.resize(e.size());
		imp_sparse_bnd_diff( +1, *eL, lower, e(), &u() );
		if(out && print_vectors)
			*out
				<< "eL - e =\n" << u();
		Vt_S( &u(), 1.0/(1.0+e_norm_inf) );
		if(out) {
			*out
				<< "Comparing: u - v <= 0\n"
				<< "u = (eL - e) | / (1 + ||e||inf ), v = 0 ...\n";
		}
		if(!comp_v.comp_less( u(), 0.0, opt_warning_tol()
			, force_inequality_error_check ? feas_error_tol() : really_big_error_tol
			, print_all_warnings, out )) return false;
		if(out)
			*out
				<< "\nChecking muL(i) * (eL - e)(i) = 0  ...\n";
		set_complementarity( *mu, u(), lower, &c );
		if(out && print_vectors)
			*out
				<< "muL(i) * (eL - e)(i) =\n" << c();
			if(out) {
				*out
					<< "Comparing:\n"
					<< "u(i) = muL(i) * (eL - e)(i), v = 0 ...\n";
			}
			if(!comp_v.comp( c(), 0.0, opt_warning_tol()
				, force_complementarity_error_check ? comp_error_tol() : really_big_error_tol
				, print_all_warnings, out )) return false;

		///////////////////////////////////
		// e - eU <= 0
		if(out)
			*out
				<< "\nChecking e - eU <= 0 ...\n";
		u.resize(e.size());
		imp_sparse_bnd_diff( -1, *eU, upper, e(), &u() );
		if(out && print_vectors)
			*out
				<< "\ne - eU =\n" << u();
		Vt_S( &u(), 1.0/(1.0+e_norm_inf) );
		if(out) {
			*out
				<< "\nComparing: u - v <= 0\n"
				<< "u = (e - eU) | / (1 + ||e||inf )\n"
				<< "v = 0 ...\n";
		}
		if(!comp_v.comp_less( u(), 0.0, opt_warning_tol()
			, force_inequality_error_check ? feas_error_tol() : really_big_error_tol
			, print_all_warnings, out )) return false;
		if(out)
			*out
				<< "\nChecking muU(i) * (e - eU)(i) = 0  ...\n";
		set_complementarity( *mu, u(), upper, &c );
		if(out && print_vectors)
			*out
				<< "\nmuU(i) * (e - eU)(i) =\n" << c();
			if(out) {
				*out
					<< "\nComparing:\n"
					<< "u(i) = muU(i) * (e - eU)(i)\n"
					<< "v = 0 ...\n";
			}
			if(!comp_v.comp( c(), 0.0, opt_warning_tol()
				, force_complementarity_error_check ? comp_error_tol() : really_big_error_tol
				, print_all_warnings, out )) return false;

	}

	if( F ) {

		///////////////////////////////////
		// r = - op(F)*d + eta * f 
		if(out)
			*out
				<< "\nComputing r = - op(F)*d + eta * f ...\n";
		V_StMtV( &r, -1.0, *F, trans_F, *d );
		Vp_StV( &r(), *eta, *f );
		if(out && print_vectors)
			*out
				<< "\nr = - op(F)*d + eta * f =\n" << r();

		if(out) {
			*out
				<< "\nChecking r == f:\n"
				<< "u = r, v = f ...\n";
		}
		if(!comp_v.comp( r(), *f, opt_warning_tol()
			, force_equality_error_check ? feas_error_tol() : really_big_error_tol
			, print_all_warnings, out )) return false;

	}

	if(out)
		*out
			<< "\n*** Congradulations, all of the enforced QP optimality conditions were within the specified"
			<< "\n*** error tolerances!\n";

	return true;	// If we get here then success

}

}	// end namespace ConstrainedOptimizationPack
