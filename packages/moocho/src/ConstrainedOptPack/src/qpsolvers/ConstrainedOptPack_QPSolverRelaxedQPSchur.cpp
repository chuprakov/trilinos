// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxedQPSchur.cpp

#include <assert.h>

#include <vector>

#include "ConstrainedOptimizationPack/include/QPSolverRelaxedQPSchur.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "SparseLinAlgPack/include/SortByDescendingAbsValue.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

namespace ConstrainedOptimizationPack {

QPSolverRelaxedQPSchur::QPSolverRelaxedQPSchur(
		  const init_kkt_sys_ptr_t&  init_kkt_sys
		, value_type		max_qp_iter_frac
		, QPSchurPack::ConstraintsRelaxedStd::EInequalityPickPolicy
							inequality_pick_policy
		, ELocalOutputLevel	print_level
		, value_type 		bounds_tol
		, value_type 		inequality_tol
		, value_type 		equality_tol
		, value_type		loose_feas_tol
		, value_type		dual_infeas_tol
		, value_type        pivot_tol
		, value_type		huge_primal_step
		, value_type		huge_dual_step
		, value_type		bigM
		, value_type		warning_tol
		, value_type		error_tol
		)
	:
	init_kkt_sys_(init_kkt_sys)
	,max_qp_iter_frac_(max_qp_iter_frac)
	,inequality_pick_policy_(inequality_pick_policy)
	,print_level_(print_level)
	,bounds_tol_(bounds_tol)
	,inequality_tol_(inequality_tol)
	,equality_tol_(equality_tol)
	,loose_feas_tol_(loose_feas_tol)
	,dual_infeas_tol_(dual_infeas_tol)
	,pivot_tol_(pivot_tol)
	,huge_primal_step_(huge_primal_step)
	,huge_dual_step_(huge_dual_step)
	,bigM_(bigM)
	,warning_tol_(warning_tol)
	,error_tol_(error_tol)
{}

QPSolverRelaxedQPSchur::~QPSolverRelaxedQPSchur()
{
	this->release_memory();
}

// Overridden from QPSolverRelaxed

QPSolverStats
QPSolverRelaxedQPSchur::get_qp_stats() const
{
	return qp_stats_;
}

void QPSolverRelaxedQPSchur::release_memory()
{
	// ToDo: Implement
	assert(0);
}

QPSolverStats::ESolutionType
QPSolverRelaxedQPSchur::imp_solve_qp(
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
	using DynamicCastHelperPack::dyn_cast;
	using LinAlgOpPack::V_mV;
	typedef QPSchurPack::ConstraintsRelaxedStd constr_t;

	const size_type
		nd = g.size(),
		m_in = E ? BLAS_Cpp::rows(E->rows(),E->cols(),trans_E) : 0,
		m_eq = F ? BLAS_Cpp::rows(F->rows(),F->cols(),trans_F) : 0;

	// Validate that G is symmetric
	const MatrixSymWithOp&
#ifdef _WINDOWS
		G_sym = dynamic_cast<const MatrixSymWithOp&>(G);
#else
		G_sym = dyn_cast<const MatrixSymWithOp>(G);
#endif

	// ///////////////////////////////
	// Setup the initial KKT system

	InitKKTSystem::i_x_free_t     i_x_free;
	InitKKTSystem::i_x_fixed_t    i_x_fixed;
	InitKKTSystem::bnd_fixed_t    bnd_fixed;
	InitKKTSystem::j_f_decomp_t   j_f_decomp;
	size_type n_R_tmp;
	init_kkt_sys().initialize_kkt_system(
		g,G,etaL,dL,dU,F,trans_F,f
		,&n_R_tmp,&i_x_free,&i_x_fixed,&bnd_fixed,&j_f_decomp
		,&b_X_,&Ko_,&fo_ );
	const size_type
		n_R = n_R_tmp,
		n_X = nd + 1 - n_R; // fixed variables in d and eta
	assert( i_x_free.size() == 0 || i_x_free.size()  >= n_R );  // Todo: Make an exception!
	assert( i_x_fixed.size() >= n_X );  // Todo: Make an exception!
	assert( bnd_fixed.size() >= n_X ); // Todo: Make and exception!

	// //////////////////////////////
	// Initialize constraints object

	// Setup j_f_undecomp
	const size_type
		m_undecomp = F ? f->size() - j_f_decomp.size() : 0;
	typedef std::vector<size_type> j_f_undecomp_t;
	j_f_undecomp_t j_f_undecomp(m_undecomp);
	if( m_undecomp ) {
		assert(0); // ToDo: Implement this when needed!
	}

	// initialize constraints object
	constraints_.initialize(
		nd,etaL,&dL,&dU,E,trans_E,b,eL,eU,F,trans_F,f
		, m_undecomp, m_undecomp ? &j_f_undecomp[0] : NULL
		, Ed
		, true	// Check the equality constraints since they will not be
		        // added to the initial active set in case they are linearly
		        // dependent!
		);
	// ToDo: Add computation of Fd to above constraints class!
	// ToDo: Add j_d_decomp to the above constraints class!

	// ///////////////////////////
	// Initialize the QP object

	// g_relaxed_ = [ g; bigM ]
	g_relaxed_.resize(nd+1);
	g_relaxed_(1,nd) = g;
	g_relaxed_(nd+1) = bigM();

	// G_relaxed_ = [ G, zeros(...); zeros(...), bigM ]
	G_relaxed_.initialize( G_sym, bigM() );
	
	qp_.initialize(
		g_relaxed_(),G_relaxed_,NULL
		,n_R, i_x_free.size() ? &i_x_free[0] : NULL
		,&i_x_fixed[0],&bnd_fixed[0]
		,b_X_(),*Ko_,fo_(),&constraints_
		,out,test_what==RUN_TESTS,warning_tol(),error_tol()
		,int(olevel)>=int(PRINT_ITER_VECTORS)
		);

	// ///////////////////////////////////////////////////////
	// Setup for a warm start (changes to initial KKT system)

	typedef std::vector<int> 					ij_act_change_t;
	typedef std::vector<EBounds>				bnds_t;
	size_type			num_act_change = 0; // The default is a cold start
	const size_type     max_num_act_change = (nu ? nu->nz() : 0) + (mu ? mu->nz() : 0) + n_X;
	ij_act_change_t		ij_act_change(max_num_act_change);
	bnds_t				bnds(max_num_act_change);

	// Note that we do not add the general equality constraints to the initial
	// guess of the active set since they may be linearly dependent!

	if( ( nu && nu->nz() ) || ( m_in && mu->nz() ) ) {
		//
		// Setup num_act_change, ij_act_change, bnds for a warm start!
		//
		const size_type
			nu_nz = nu ? nu->nz() : 0,
			mu_nz = mu ? mu->nz() : 0;
		// Combine all the multipliers for the bound and general inequality
		// constraints and sort them from the largest to the smallest.  Hopefully
		// the constraints with the larger multiplier values will not be dropped
		// from the active set.
		SpVector gamma( nd + 1 + m_in , nu_nz + mu_nz );
		typedef SpVector::element_type ele_t;
		if(nu && nu->nz()) {
			const SpVector::difference_type o = nu->offset();
			for( SpVector::const_iterator itr = nu->begin(); itr != nu->end(); ++itr ) {
				gamma.add_element( ele_t( itr->indice() + o, itr->value() ) );
			}
		}
		if(mu && mu->nz()) {
			const SpVector::difference_type o = mu->offset() + nd + 1;
			for( SpVector::const_iterator itr = mu->begin(); itr != mu->end(); ++itr ) {
				gamma.add_element( ele_t( itr->indice() + o, itr->value() ) );
			}
		}
		std::sort( gamma.begin(), gamma.end()
			, SparseLinAlgPack::SortByDescendingAbsValue() );
		// Now add the inequality constraints in decreasing order (if they are
		// not already initially fixed variables)
		const QPSchurPack::QP::x_init_t &x_init = qp_.x_init();
		if(gamma.nz()) {
			const SpVector::difference_type o = gamma.offset();
			for( SpVector::const_iterator itr = gamma.begin(); itr != gamma.end(); ++itr ) {
				const size_type i =  itr->indice() + o;
				if( i <= nd && x_init(i) != FREE )
					continue; // This variable is already initially fixed
				// This is not an initially fixed variable so add it
				ij_act_change[num_act_change] = i;
				bnds[num_act_change]
					= itr->value() > 0.0 ? UPPER : LOWER;
				++num_act_change;
			}
		}
	}
	// We need to loop through x_init() and nu() in order and look for variables
	// that are initially fixed in x_init() but are not present in nu().  For these
	// variables we need to free them in ij_act_change[].
	{
		QPSchurPack::QP::x_init_t::const_iterator
			x_init_itr = qp_.x_init().begin();
		const SpVector::difference_type o = nu->offset();
		SpVector::const_iterator
			nu_itr = const_cast<const SpVector*>(nu)->begin(), // Okay even if nu.nz() == 0
			nu_end = const_cast<const SpVector*>(nu)->end();
		for( size_type i = 1; i <= nd; ++i, ++x_init_itr ) {
			if( *x_init_itr != FREE
				&& *x_init_itr != EQUALITY )
			{
				// This is an initially fixed upper or lower bound
				// Look for lagrange multiplier stating that it is
				// still fixed.
				if( nu_itr != nu_end && (nu_itr->indice() + o) < i ) {
					++nu_itr;
				}
				else if( nu_itr != nu_end && (nu_itr->indice() + o) == i ) {
					// This active bound is present but lets make sure
					// that it is still the same bound
					if( ( *x_init_itr == LOWER && nu_itr->value() > 0 )
						|| ( *x_init_itr == UPPER && nu_itr->value() < 0 ) )
					{
						// The bound has changed from upper to lower or visa-versa!
						ij_act_change[num_act_change] = i;
						bnds[num_act_change]
							= nu_itr->value() > 0.0 ? UPPER : LOWER;
						++num_act_change;
					}
				}
				else {
					// This initially fixed variable is not fixed in nu so lets free it!
					ij_act_change[num_act_change] = -i;
					bnds[num_act_change]          = FREE;
					++num_act_change;
				}
			}
		}
	}
	// Set the output level
	QPSchur::EOutputLevel qpschur_olevel;
	switch( print_level() ) {
		case USE_INPUT_ARG: {
			// Use the input print level
			switch( olevel ) {
				case PRINT_NONE:
					qpschur_olevel = QPSchur::NO_OUTPUT;
					break;
				case PRINT_BASIC_INFO:
					qpschur_olevel = QPSchur::OUTPUT_BASIC_INFO;
					break;
				case PRINT_ITER_SUMMARY:
					qpschur_olevel = QPSchur::OUTPUT_ITER_SUMMARY;
					break;
				case PRINT_ITER_STEPS:
					qpschur_olevel = QPSchur::OUTPUT_ITER_STEPS;
					break;
				case PRINT_ITER_ACT_SET:
				case PRINT_ITER_VECTORS:
					qpschur_olevel = QPSchur::OUTPUT_ACT_SET;
					break;
				case PRINT_EVERY_THING:
					qpschur_olevel = QPSchur::OUTPUT_ITER_QUANTITIES;
					break;
				default:
					assert(0);
			}
			break;
		}
		case NO_OUTPUT:
			qpschur_olevel = QPSchur::NO_OUTPUT;
			break;
		case OUTPUT_BASIC_INFO:
			qpschur_olevel = QPSchur::OUTPUT_BASIC_INFO;
			break;
		case OUTPUT_ITER_SUMMARY:
			qpschur_olevel = QPSchur::OUTPUT_ITER_SUMMARY;
			break;
		case OUTPUT_ITER_STEPS:
			qpschur_olevel = QPSchur::OUTPUT_ITER_STEPS;
			break;
		case OUTPUT_ACT_SET:
			qpschur_olevel = QPSchur::OUTPUT_ACT_SET;
			break;
		case OUTPUT_ITER_QUANTITIES:
			qpschur_olevel = QPSchur::OUTPUT_ITER_QUANTITIES;
			break;
		default:
			assert(0);
	}

	//
	// Set options for ConstraintsRelaxedStd.
	// 
	if( bounds_tol() > 0.0 )
		constraints_.bounds_tol(bounds_tol());
	if( inequality_tol() > 0.0 )
		constraints_.inequality_tol(inequality_tol());
	if( equality_tol() > 0.0 )
		constraints_.equality_tol(equality_tol());
	constraints_.inequality_pick_policy(inequality_pick_policy());

	//
	// Set options for QPSchur.
	// 
	qp_solver_.max_iter( max_qp_iter_frac() * nd );
	qp_solver_.feas_tol( constraints_.bounds_tol() );	// Let's assume the bound tolerance is the tightest
	if(loose_feas_tol() > 0.0)
		qp_solver_.loose_feas_tol( loose_feas_tol() );
	else
		qp_solver_.loose_feas_tol( 10.0 * qp_solver_.feas_tol() );
	if(dual_infeas_tol() > 0.0)
		qp_solver_.dual_infeas_tol( dual_infeas_tol() );
	if(huge_primal_step() > 0.0)
		qp_solver_.huge_primal_step( huge_primal_step() );
	if(huge_dual_step() > 0.0)
		qp_solver_.huge_dual_step( huge_dual_step() );
	if(pivot_tol() > 0.0)
		schur_comp_.pivot_tol(pivot_tol());
	qp_solver_.set_schur_comp( QPSchur::schur_comp_ptr_t( &schur_comp_, false ) );
	qp_solver_.warning_tol( warning_tol() );
	qp_solver_.error_tol( error_tol() );
	
	//
	// Solve the QP with QPSchur
	// 
	Vector _x(nd+1);		// solution vector [ d; eta ]
	SpVector _mu;			// lagrange multipliers for variable bounds [ nu; kappa ]
	SpVector _lambda_breve;	// solution for extra general constraints [ mu; lambda ]
	size_type qp_iter = 0, num_adds = 0, num_drops = 0;
	QPSchur::ESolveReturn
		solve_returned
			= qp_solver_.solve_qp(
				  qp_
				, num_act_change, num_act_change ? &ij_act_change[0] : NULL
					, num_act_change ? &bnds[0] : NULL
				, out, qpschur_olevel
					, test_what==RUN_TESTS ? QPSchur::RUN_TESTS : QPSchur::NO_TESTS
				, &_x(), &_mu, NULL, &_lambda_breve
			, &qp_iter, &num_adds, &num_drops
			);

	// Set the solution

	// d
	*d = _x(1,nd);
	// nu
	if( nu ) {
		nu->resize( nd, std::_MIN( nd, _mu.nz() ) );
		const SpVector::difference_type o = _mu.offset();
		if( _mu.nz() ) {
			for(SpVector::const_iterator itr = _mu.begin(); itr != _mu.end(); ++itr)
			{
				typedef SpVector::element_type ele_t;
				if( itr->indice() + o <= nd ) // don't add multiplier for eta <= etaL
					nu->add_element( ele_t( itr->indice() + o, itr->value() ) );
			}
		}
		nu->assume_sorted(true);
	}
	// mu, lambda	
	if( m_in || m_eq ) {
		*eta = _x(nd+1);	// must be non-null
		if(mu) mu->resize(m_in,_lambda_breve.nz());	// resize for storage for all inequalities active?
		const SpVector::difference_type o = _lambda_breve.offset();
		if(_lambda_breve.nz()) {
			for(SpVector::const_iterator itr = _lambda_breve.begin();
				itr != _lambda_breve.end();
				++itr)
			{
				typedef SpVector::element_type ele_t;
				if(itr->indice() <= m_in) {
					assert(mu);	// Local error only (if validated properly)
					mu->add_element(ele_t(itr->indice(),itr->value()));
				}
				else {
					assert(lambda);	// Local error only (if validated properly)
					(*lambda)(itr->indice() - m_in) = itr->value();
				}
			}
		}
		mu->assume_sorted(true);
	}
	// obj_d (This could be updated within QPSchur in the future)
	if(obj_d) {
		// obj_d = g'*d + 1/2 * d' * G * g
		*obj_d = LinAlgPack::dot(g,*d)
			+ 0.5 * SparseLinAlgPack::transVtMtV(*d,G,BLAS_Cpp::no_trans,*d);
	}
	// Ed
	if(Ed && E) {
		switch(constraints_.inequality_pick_policy()) {
			case constr_t::ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY:
				if(solve_returned == QPSchur::OPTIMAL_SOLUTION)
					break; // Ed already computed (see ConstraintsRelaxedStd::pick_violated())
			case constr_t::ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY:
				break; // Ed already computed (see ConstraintsRelaxedStd::pick_violated())
			default:
				// We need to compute Ed
				LinAlgOpPack::V_MtV( Ed, *E, trans_E, *d );
		}
	}
	// Fd (This could be updated within ConstraintsRelaxedStd in the future)
	if(Fd) {
		LinAlgOpPack::V_MtV( Fd, *F, trans_F, *d );
	}
	// Set the QP statistics
	QPSolverStats::ESolutionType solution_type;
	QPSolverStats::EConvexity convexity = QPSolverStats::CONVEX;
	switch( solve_returned ) {
		case QPSchur::OPTIMAL_SOLUTION:
			solution_type = QPSolverStats::OPTIMAL_SOLUTION;
			break;
		case QPSchur::MAX_ITER_EXCEEDED:
			solution_type = QPSolverStats::DUAL_FEASIBLE_POINT;
			break;
		case QPSchur::MAX_ALLOWED_STORAGE_EXCEEDED:
			solution_type = QPSolverStats::DUAL_FEASIBLE_POINT;
			break;
		case QPSchur::INFEASIBLE_CONSTRAINTS:
		case QPSchur::NONCONVEX_QP:
			convexity = QPSolverStats::NONCONVEX;
		case QPSchur::DUAL_INFEASIBILITY:
		case QPSchur::SUBOPTIMAL_POINT:
			solution_type = QPSolverStats::SUBOPTIMAL_POINT;
			break;
		default:
			assert(0);
	}
	qp_stats_.set_stats(
		solution_type,convexity,qp_iter,num_adds,num_drops
		, num_act_change > 0 || n_X > 1, *eta > 0.0 );

	return qp_stats_.solution_type();
}

}	// end namespace ConstrainedOptimizationPack
