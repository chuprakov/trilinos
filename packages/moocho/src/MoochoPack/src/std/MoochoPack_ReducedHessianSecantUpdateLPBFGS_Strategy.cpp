// ///////////////////////////////////////////////////////
// ReducedHessianSecantUpdateLPBFGS_Strategy.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include "ReducedSpaceSQPPack/include/std/ReducedHessianSecantUpdateLPBFGS_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/get_init_fixed_free_indep.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgo.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"
#include "ConstrainedOptimizationPack/include/MatrixSymPosDefLBFGS.h"
#include "ConstrainedOptimizationPack/include/MatrixSymPosDefCholFactor.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/MatrixWithOpOut.h"
#include "SparseLinAlgPack/include/GenPermMatrixSlice.h"
#include "SparseLinAlgPack/include/GenPermMatrixSliceOp.h"
#include "SparseLinAlgPack/include/MatrixSymInitDiagonal.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

namespace ReducedSpaceSQPPack {

ReducedHessianSecantUpdateLPBFGS_Strategy::ReducedHessianSecantUpdateLPBFGS_Strategy(
	const dense_proj_update_ptr_t&  dense_proj_update
	,const bfgs_update_ptr_t&       bfgs_update
	,value_type                     act_set_frac_proj_start
	,value_type                     super_basic_mult_drop_tol
	,size_type                      num_superbasics_switch_dense
	,value_type                     act_set_frac_switch_dense
	,size_type                      min_num_updates_switch_dense
	,size_type                      max_num_updates_switch_dense
	,size_type                      num_add_recent_updates
	)
	: dense_proj_update_(dense_proj_update)
	, bfgs_update_(bfgs_update)
	, act_set_frac_proj_start_(act_set_frac_proj_start)
	, super_basic_mult_drop_tol_(super_basic_mult_drop_tol)
	, num_superbasics_switch_dense_(num_superbasics_switch_dense)
	, act_set_frac_switch_dense_(act_set_frac_switch_dense)
	, min_num_updates_switch_dense_(min_num_updates_switch_dense)
	, max_num_updates_switch_dense_(max_num_updates_switch_dense)
	, num_add_recent_updates_(num_add_recent_updates)
{}

bool ReducedHessianSecantUpdateLPBFGS_Strategy::perform_update(
	VectorSlice* s_bfgs, VectorSlice* y_bfgs, bool first_update
	,std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
	,MatrixWithOp *rHL_k
	)
{
	using std::setw;
	using std::endl;
	using std::right;
	using DynamicCastHelperPack::dyn_cast;
	namespace rcp = ReferenceCountingPack;
	using rcp::ref_count_ptr;
	using LinAlgOpPack::V_MtV;
	using SparseLinAlgPack::norm_inf;
	typedef ConstrainedOptimizationPack::MatrixHessianSuperBasic MHSB_t;
	

	const size_type
		n    = algo->nlp().n(),
		r    = algo->nlp().r(),
		n_pz = n-r;

#ifdef _WINDOWS
	MHSB_t &rHL_super = dynamic_cast<MHSB_t&>(*rHL_k);
#else
	MHSB_t &rHL_super = dyn_cast<MHSB_t>(*rHL_k);
#endif

	// See if we still have a limited memory BFGS update matrix
	ref_count_ptr<MatrixSymPosDefLBFGS> // We don't want this to be deleted until we are done with it
		lbfgs_rHL_RR = rcp::rcp_const_cast<MatrixSymPosDefLBFGS>(
			rcp::rcp_dynamic_cast<const MatrixSymPosDefLBFGS>(rHL_super.B_RR_ptr()) );

	if( lbfgs_rHL_RR.get() ) {
		//
		// We have a limited memory BFGS matrix
		//

		// Determine when it is time to switch to dense BFGS
		const bool consider_switch = lbfgs_rHL_RR->m_bar() >= min_num_updates_switch_dense();
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << "\n*** (LPBFGS) Limited memory BFGS to dense projected BFGS\n" 
				<< "\nnum_previous_updates = " << lbfgs_rHL_RR->m_bar()
				<< ( consider_switch ? " >= " : " < " )
				<< "min_num_updates_switch_dense = " << min_num_updates_switch_dense()
				<< ( consider_switch
					 ? "\nConsidering if we should switch to dense projected BFGS updating of superbasics ...\n"
					 : "\nNot time to consider switching to dense projected BFGS updating of superbasics yet" );
		}
		bool switched_to_dense = false;
		if( consider_switch ) {
			// 
			// Our job here is to determine if it is time to switch to dense projected
			// BFGS updating.
			//
			if( act_set_stats_(*s).updated_k(-1) ) {
				if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
					out	<< "\nDetermining if projected BFGS updating of superbasics should begin ...\n";
				}
				ActSetStats &stats = act_set_stats_(*s).get_k(-1);
				const SpVector
					&nu_km1 = s->nu().get_k(-1);
				const SpVectorSlice
					nu_indep = nu_km1(s->var_indep());
				const size_type
					num_active_indep = nu_indep.nz(),
					num_adds         = stats.num_adds(),
					num_drops        = stats.num_drops();
				const value_type
					frac_same
					= ( num_adds == ActSetStats::NOT_KNOWN || num_active_indep == 0
						? 0.0
						: std::_MAX(((double)(num_active_indep)-num_adds-num_drops) / num_active_indep, 0.0 ) );
				const bool low_num_super_basics = n_pz - num_active_indep <= num_superbasics_switch_dense();
				const bool act_set_calmed_down = ( num_active_indep > 0 && frac_same >= act_set_frac_switch_dense() );
				const bool max_num_updates_exceeded = lbfgs_rHL_RR->m_bar() >= max_num_updates_switch_dense();
				bool switch_to_dense = low_num_super_basics && ( act_set_calmed_down || max_num_updates_exceeded );
				if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
					out << "\n(n-r) - num_act_indep = " << n_pz << " - " << num_active_indep << " = " << n_pz - num_active_indep
						<< ( low_num_super_basics ? " <= " : " > " )
						<< " num_superbasics_switch_dense = " << num_superbasics_switch_dense();
					if( num_active_indep ) {
						out	<< "\nmax(num_active_indep-num_adds-num_drops,0)/(num_active_indep) = "
							<< "max("<<num_active_indep<<"-"<<num_adds<<"-"<<num_drops<<",0)/("<<num_active_indep<<") = "
							<< frac_same;
						if( act_set_calmed_down )
							out << " >= ";
						else
							out << " < ";
						out << "act_set_frac_switch_dense = " << act_set_frac_switch_dense();
					}
					if( low_num_super_basics && act_set_calmed_down ) {
						out << "\nThe active set has calmed down enough and their are not too many super basic variables so lets"
							<< "\nfurther consider switching to dense projected BFGS updating of superbasic variables ...\n";
					}
					else if( low_num_super_basics && max_num_updates_exceeded ) {
						out << "\nThe active set has not calmed down enough but their are not too many super basic variables"
							<< "\nand num_previous_updates = " << lbfgs_rHL_RR->m_bar() << " >= max_num_updates_switch_dense = "
							<< max_num_updates_switch_dense()
							<< "\nso we will further consider switching to dense projected BFGS updating of superbasic variables ...\n";
					}
					else {
						out << "\nIt is not time to switch so just keep performing limited memory BFGS for now ...\n";
					}
				}
				if( switch_to_dense ) {
					// Determine the set of initially fixed and free independent variables
					std::vector<size_type>                            // Todo: use temp workspace
						i_x_free(n_pz), 
						i_x_fixed(n_pz);
					std::vector<ConstrainedOptimizationPack::EBounds> // Todo: use temp workspace
						bnd_fixed(n_pz);
					size_type
						n_pz_X = 0,
						n_pz_R = 0;
					get_init_fixed_free_indep(
						n,r,nu_indep,super_basic_mult_drop_tol(),olevel,out
						,&n_pz_X,&n_pz_R,&i_x_free[0],&i_x_fixed[0],&bnd_fixed[0] );
					if( n_pz_R > num_superbasics_switch_dense() ) {
						if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
							out << "\nOops! the actual number of super basic variables that we would keep n_pz_R = " << n_pz_R
								<< " > num_superbasics_switch_dense = " << num_superbasics_switch_dense()
								<< "\nWe will have to keep preforming limited memory BFGS for now!"
								<< "(consider decreasing the value of super_basic_mult_drop_tol = " << super_basic_mult_drop_tol()
								<< std::endl;
						}
					}
					else {
						// Create new dense rHL_RR matrix and reinitialize rHL for set of superbasics
						if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
							out	<< "\nCreate new dense BFGS matrix rHL_RR for superbasics only and initialize to rHL_RR = (1/gamma_k)*I ("
								<< "where gamma_k = " << lbfgs_rHL_RR->gamma_k() << ") ...\n";
						}
						ref_count_ptr<MatrixSymPosDefCholFactor>
							rHL_RR = new MatrixSymPosDefCholFactor(
								NULL    // Let it allocate its own memory
								,NULL   // ...
								,lbfgs_rHL_RR->maintain_original()
								,lbfgs_rHL_RR->maintain_inverse()
								);
						rHL_RR->init_identity( n_pz_R, 1.0/lbfgs_rHL_RR->gamma_k() ); // rHL_RR = (1/gamma_k)*I
						if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
							out << "\nDense rHL_RR after rHL_RR = (1/gamma_k)*I but before adding previous BFGS updates\nrHL_RR =\n" << *rHL_RR;
						}
						// Initialize rHL_XX = I for now
#ifdef _WINDOWS
						MatrixSymInitDiagonal
							&rHL_XX = dynamic_cast<MatrixSymInitDiagonal&>(
								const_cast<MatrixSymWithOp&>(*rHL_super.B_XX_ptr()));
#else
						MatrixSymInitDiagonal
							&rHL_XX = dyn_cast<MatrixSymInitDiagonal>(
								const_cast<MatrixSymWithOp&>(*rHL_super.B_XX_ptr()));
#endif
						rHL_XX.init_identity( n_pz_X ); // ToDo: We may want scaled element?
						rHL_super.initialize(
							n_pz, n_pz_R, &i_x_free[0], &i_x_fixed[0], &bnd_fixed[0]
							,rcp::rcp_static_cast<const MatrixSymWithOpFactorized>(rHL_RR)
							,NULL,BLAS_Cpp::no_trans,rHL_super.B_XX_ptr()
							);
						if( lbfgs_rHL_RR->m_bar() ) {
							const size_type num_add_updates = std::_MIN(num_add_recent_updates(),lbfgs_rHL_RR->m_bar());
							if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
								out	<< "\nAdd the min(num_previous_updates,num_add_recent_updates) = min(" << lbfgs_rHL_RR->m_bar()
									<< "," << num_add_recent_updates() << ") = " << num_add_updates << " most recent BFGS updates ...\n";
							}
							// Now add previous update vectors
							const GenPermMatrixSlice
								&Q_R = rHL_super.Q_R(),
								&Q_X = rHL_super.Q_X();
							const size_type
								n_pz_R = Q_R.cols(),
								n_pz_X = Q_X.cols();
							assert( n_pz_R + n_pz_X == n_pz );
							const GenMatrixSlice
								S = lbfgs_rHL_RR->S(),
								Y = lbfgs_rHL_RR->Y();
							Vector s_bfgs_R(n_pz_R);  // ToDo: use temp workspace for this
							Vector y_bfgs_R(n_pz_R);  // ToDo: use temp workspace for this
							size_type k = lbfgs_rHL_RR->k_bar();  // Location in S and Y of must recent update vectors
							for( size_type l = 1; l <= num_add_updates; ++l, --k ) {
								if(k == 0) k = lbfgs_rHL_RR->m_bar();  // see MatrixSymPosDefLBFGS
								// s_bfgs_R = Q_R'*s_bfgs
								V_MtV( &s_bfgs_R(), Q_R, BLAS_Cpp::trans, S.col(k) );
								// y_bfgs_R = Q_R'*y_bfgs
								V_MtV( &y_bfgs_R(), Q_R, BLAS_Cpp::trans, Y.col(k) );
								// ToDo: We could print these updates for debugging if needed!
								// ( rHL_RR, s_bfgs_R, y_bfgs_R ) -> rHL_RR (this should not throw an exception!)
								try {
									rHL_RR->secant_update( &s_bfgs_R(), &y_bfgs_R(), NULL );
								    // ToDo: We could test the secant condition if we have to.
									// ToDo: We could print the updated BFGS matrix if we had to.
								}
								catch( const MatrixSymSecantUpdateable::UpdateSkippedException& excpt ) {
									if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
										out	<< "\nOops!  The " << l << "th most recent BFGS update was rejected?:\n"
											<< excpt.what() << std::endl;
									}
								}
							}
							if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
								out << "\nDense rHL_RR after adding previous BFGS updates\nrHL_BRR =\n" << *rHL_RR;
							}
						}
						else {
							if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
								out	<< "\nThere were no previous update vectors to add!\n";
							}
						}
						if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
							out << "\nFull rHL after reinitialization but before adding previous BFGS updates\nrHL =\n" << *rHL_k;
						}
						switched_to_dense = true; // Our work is done now.
					}
				}
			}
		}
		// If we have not switched to dense then just update the limited memory BFGS matrix
		if(!switched_to_dense) {
			bfgs_update().perform_update(
				s_bfgs,y_bfgs,first_update,out,olevel,algo->algo_cntr().check_results()
				,lbfgs_rHL_RR.get()
				,&quasi_newton_stats_(*s).set_k(0)
				);
			return true;
		}
	}
	else {
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out	<< "\nThe matrix rHL_RR is not limited memory so we must have already switched to dense ...\n";
		}
	}
	//
	// If we get here then we must have switched to dense
	// projected updating so lets just pass it on!
	//
	return dense_proj_update().perform_update(
		s_bfgs,y_bfgs,first_update,out,olevel,algo,s,rHL_k);
}

void ReducedHessianSecantUpdateLPBFGS_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Perform limited memory LBFGS updating initially then switch to dense\n"
		<< L << "*** projected BFGS updating when appropriate.\n"
		<< L << "ToDo: Finish implementation!\n";
}

}  // end namespace ReducedSpaceSQPPack
