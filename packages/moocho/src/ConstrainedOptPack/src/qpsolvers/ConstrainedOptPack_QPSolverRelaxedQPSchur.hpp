// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxedQPSchur.h

#ifndef QP_SOLVER_RELAXED_QP_SCHUR_H
#define QP_SOLVER_RELAXED_QP_SCHUR_H

#include "QPSolverRelaxed.h"
#include "QPSchur.h"
#include "QPInitFixedFreeStd.h"
#include "MatrixHessianRelaxed.h"
#include "ConstraintsRelaxedStd.h"
#include "MatrixSymAddDelBunchKaufman.h"
#include "Misc/include/StandardCompositionMacros.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ConstrainedOptimizationPack {

///
/** Solves Quadratic Programming (QP) problems using QPSchur.
  *
  * This is the only subclass needed for QPSchur.  All of the specifics of how the
  * initial KKT system is formed is delegated to a strategy object of type
  * \Ref{InitKKTSystem} (see below).  Note that the matrix #G# must support the
  * #MatrixSymWithOp# interface (which should be no problem!)
  */
class QPSolverRelaxedQPSchur
	: public QPSolverRelaxed
{
public:

	///
	/** Interface for the object that forms the initial KKT system {abstract}.
	 *
	 * Note that this interface is set up such that the relaxation variable
	 * must always be initially fixed (and rightly so to avoid illconditioning).
	 */
	class InitKKTSystem {
	public:
		///
		typedef std::vector<size_type> i_x_free_t;
		///
		typedef std::vector<size_type> i_x_fixed_t;
		///
		typedef std::vector<EBounds>   bnd_fixed_t;
		///
		typedef std::vector<size_type> j_f_decomp_t;
		///
		typedef ReferenceCountingPack::ref_count_ptr<const MatrixSymWithOpFactorized>
			Ko_ptr_t;
		///
		virtual ~InitKKTSystem() {}
		///
		/** Initializes the KKT system.
		 *
		 * Let the following permutation matrices define the selection of the
		 * initial KKT system:
		 *
		 * #Q = [ Q_R, Q_X ]# : Initially fixed #Q_R# and free #Q_X# variables
		 *
		 * #P = [ P_d, P_u ]# : Decomposed #P_d# and undecomposed #P_u# constraints
		 *
		 * Given the definitions of #Q# and #P# above, this function will return
		 * the initial KKT system:
		 \begin{verbatim}
		 Ko = [ Q_R'*G*Q_X     op(F')*P_d ]
		      [ P_d'*op(F)          0     ]

		 fo = [ -Q_R'*g - Q_R'*G*Q_X*b_X    ]
		      [ -P_d'f - P_d'*op(F)*Q_X*b_X ]

		 b_X = ??? (see below)
		 \end{verbatim}
		 *
		 * @param  g    [in] See #QPSolverRelaxed::solve_qp(...)
		 * @param  G    [in] See #QPSolverRelaxed::solve_qp(...)
		 * @param  dL   [in] See #QPSolverRelaxed::solve_qp(...)
		 * @param  dU   [in] See #QPSolverRelaxed::solve_qp(...)
		 * @param  F    [in] See #QPSolverRelaxed::solve_qp(...)
		 * @param  trans_f
		 *              [in] See #QPSolverRelaxed::solve_qp(...)
		 * @param  f    [in] See #QPSolverRelaxed::solve_qp(...)
		 * @param  n_R  [out] Number of initially free variables.
		 * @param  i_x_free
		 *              [out] array (size #n_R# or #0#):
		 *              If #i_x_free.size() > 0# then #i_x_free[l-1], l = 1...n_R#
		 *              defines the matrix #Q_R# as:\\
		 *              #Q_R(:,l) = e(i_x_free[l-1]), l = 1...n_R#\\
		 *              If #i_x_free.size() == 0# then #i_x_free# is implicitly
		 *              identity and #Q_R# is defiend as:\\
		 *              #Q_R(:,l) = e(l), l = 1...n_R#\\
		 *              The ordering of these indices is significant.
		 *@param  i_x_fixed
		 *              [out] array (size #n_X#):
		 *              #i_x_fixed[l-1], l = 1...n_X# defines the matrix #Q_X# as:\\
		 *              #Q_X(:,l) = e(i_x_fixed[l-1]), l = 1...n_X#\\
		 *              The ordering of these indices is significant.
		 * @param  bnd_fixed
		 *             [out] array (size #n_X#):
		 *             #bnd_fixed[l-1], l = 1...n_X# defines the initial active set as:\\
		 *\begin{verbatim}
                           / LOWER : b_X(l) = dL(i_x_fixed[l-1])
		 bnd_fixed[l-1] = |  UPPER : b_X(l) = dU(i_x_fixed[l-1])
		                   \ EQUALITY : b_X(l) = dL(i) = dU(i) (i = i_x_fixed[l-1])
		 \end{verbatim}
		 * @param  j_f_decomp
		 *             [out] array (size #m#):
		 *             #j_f_decomp[p-1], p = 1...m# defines the decomposed equalities included
		 *             in #Ko# as:\\
		 *             #P_d(:,p) = e(j_f_decomp[p-1]), p = 1...m#\\
		 *             The ordering of these indices is significant and are not necessarily
		 *             sorted in assending or decending order.
		 * @param  b_X [out] vector (size #n_X#):
		 *             Initial varaible bounds (see #bnd_fixed# above).  Note that
		 *             the relaxation variable is always one of the initially fixed
		 *             variables.
		 * @param  Ko  [in/out] Initial KKT matrix (size #(n_R+m) x (n_R+m)#).
		 *             On output, Ko will contain a possibly dynamically allocated nonsingular
		 *             matrix object that represents Ko.  In input, if #Ko->get() != NULL#,
		 *             and no other objects have a reference to this object (based on
		 *             #Ko->count()#, and it is of the  proper type, then this matrix may be reused.
		 * @param  fo  [out] vector (size #n_R + m#) of the rhs for the initial KKT system.
		 */
		virtual void initialize_kkt_system(
			const VectorSlice&    g
			,const MatrixWithOp&  G
			,value_type           etaL
			,const SpVectorSlice& dL
			,const SpVectorSlice& dU
			,const MatrixWithOp*  F
			,BLAS_Cpp::Transp     trans_F
			,const VectorSlice*   f
			,size_type*           n_R
			,i_x_free_t*          i_x_free
			,i_x_fixed_t*         i_x_fixed
			,bnd_fixed_t*         bnd_fixed
			,j_f_decomp_t*        j_f_decomp
			,Vector*              b_X
			,Ko_ptr_t*            Ko
			,Vector*              fo
			) const = 0;

	}; // end class InitKKTSystem

	///
	/** Interface for the object that can reform an initial KKT system
	 * dynamically {abstract}.
	 *
	 * This interface allows the definition of the initial KKT system to
	 * be changed on the fly.  This may not be possible for may different
	 * QPs so this is an optional interface.  Allowing a redefinition of the
	 * initial KKT system may allow QPs with more degrees of freedom and lots
	 * of changes to the initial active set to be efficiently solved.
	 */
	class ReinitKKTSystem : public InitKKTSystem {
	public:
		// ToDo: Create method reinitailze_kkt_system(...)
	}; // end class ReinitKKTSystem


	///
	/** Strategy object that sets up the initial KKT system.
	  */
	STANDARD_COMPOSITION_MEMBERS( InitKKTSystem, init_kkt_sys )

	///
	/** Set the maximum number of QP iterations as max_itr = max_qp_iter_frac * n.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_qp_iter_frac )

	///
	/** Set the maximum real runtime in minutes.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_real_runtime )

	///
	/** <<std member comp>> members policy used to select a violated constraint.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( QPSchurPack::ConstraintsRelaxedStd::EInequalityPickPolicy
													, inequality_pick_policy )

	/// Output level
	enum ELocalOutputLevel {
		 USE_INPUT_ARG				= -1	// Use the value input to solve_qp(...)
		,NO_OUTPUT					= 0		//
		,OUTPUT_BASIC_INFO			= 1		// values sent to QPSchur::solve_qp(...)
		,OUTPUT_ITER_SUMMARY		= 2		// ...
		,OUTPUT_ITER_STEPS			= 3
		,OUTPUT_ACT_SET				= 4
		,OUTPUT_ITER_QUANTITIES		= 5
	};

	///
	/** Set the output level for QPSchur.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( ELocalOutputLevel, print_level )

	///
	/** Set the feasibility tolerance for the bound constriants.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, bounds_tol )

	///
	/** Set the feasibility tolerance for the general inequality constraints.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, inequality_tol )

	///
	/** Set the feasibility tolerance for the general equality constriants.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, equality_tol )

	///
	/** Set a looser feasibility tolerance ( > feas_tol )
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, loose_feas_tol )

	///
	/** Set the tolerence where a scaled Langrange multiplier is considered
	  * degenerate.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, dual_infeas_tol )

	///
	/** Set the tolerence for the size of the step in the primal space that is considered
	  * to be a near infinite step.  This is used to determine if the KKT
	  * system is near singular.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, huge_primal_step )

	///
	/** Set the tolerence for the size of the step in the dual space that is considered
	  * to be a near infinite step.  This is used to determine if the constriants
	  * are infeasible.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, huge_dual_step )

	///
	/** <<std member comp>> members for the Big M parameter used in the objective.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, bigM )

	///
	/** <<std member comp>> members for the warning tolerance for tests.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol )

	///
	/** <<std member comp>> members for the error tolerance for tests.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol )

	///
	/** Set the minimum number of refinement iterations to perform
	 * when using iterative refinement.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, iter_refine_min_iter )
		
	///
	/** Set the maximum number of refinement iterations to perform
	 * when using iterative refinement.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, iter_refine_max_iter )

	///
	/** Set the maxinum scaled tolerance the residual of the optimality conditions
	 * must be before terminating iterative refinement.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, iter_refine_opt_tol )

	///
	/** Set the maxinum scaled tolerance the residual of the feasibility conditions
	 * must be before terminating iterative refinement.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, iter_refine_feas_tol )

	///
	/** Set whether iterative refinement is automatically used once the solution
	 * is found.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, iter_refine_at_solution )

	///
	/** Set the relative tolerance for pivots in the schur complement under
	 * which a waning will be printed (see MatrixSymAddDelUpdateable) for
	 * near singular updates.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, pivot_warning_tol )

	///
	/** Set the relative tolerance for pivots in the schur complement under
	 * which a singularity exception will be thrown (see MatrixSymAddDelUpdateable)
	 * for singular updates.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, pivot_singular_tol )

	///
	/** Set the relative tolerance for pivots in the schur complement over
	 * which a wrong inertia exception will be throw (see MatrixSymAddDelUpdateable)
	 * for updates with the wrong inertia.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, pivot_wrong_inertia_tol )

	///
	/** Set whether equality constriants are to be added to the active set
	 * initialy to the schur complement or not.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, add_equalities_initially )


	///
	QPSolverRelaxedQPSchur(
		const init_kkt_sys_ptr_t&    init_kkt_sys       = NULL
		,value_type                  max_qp_iter_frac   = 10.0
		,value_type                  max_real_runtime   = 1e+20
		,QPSchurPack::ConstraintsRelaxedStd::EInequalityPickPolicy
		                             inequality_pick_policy
		                                 = QPSchurPack::ConstraintsRelaxedStd::ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY
		,ELocalOutputLevel           print_level             = USE_INPUT_ARG // Deduce from input arguments
		,value_type                  bounds_tol              = -1.0           // use default
		,value_type                  inequality_tol          = -1.0	          // use default
		,value_type                  equality_tol            = -1.0	          // use default
		,value_type                  loose_feas_tol          = -1.0	          // use default
		,value_type                  dual_infeas_tol         = -1.0	          // use default
		,value_type                  huge_primal_step        = -1.0	          // use defalut
		,value_type                  huge_dual_step          = -1.0	          // use default
		,value_type                  bigM                    = 1e+10
		,value_type                  warning_tol             = 1e-10
		,value_type                  error_tol               = 1e-5
		,size_type                   iter_refine_min_iter    = 1
		,size_type                   iter_refine_max_iter    = 3
		,value_type                  iter_refine_opt_tol     = 1e-12
		,value_type                  iter_refine_feas_tol    = 1e-12
		,bool                        iter_refine_at_solution = true
		,value_type                  pivot_warning_tol       = 1e-8
		,value_type                  pivot_singular_tol      = 1e-11
		,value_type                  pivot_wrong_inertia_tol = 1e-11
		,bool                        add_equalities_initially= true
		);

	///
	~QPSolverRelaxedQPSchur();

	// /////////////////////////////////
	// Overridden from QPSolverRelaxed

	///
	QPSolverStats get_qp_stats() const;

	///
	void release_memory();

protected:

	// /////////////////////////////////
	// Overridden from QPSolverRelaxed

	///
	QPSolverStats::ESolutionType imp_solve_qp(
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
		);

private:

	// ////////////////////////////
	// Private data members

	QPSolverStats                    qp_stats_;
	QPSchur                          qp_solver_;
	QPSchurPack::QPInitFixedFreeStd  qp_;
	MatrixHessianRelaxed	         G_relaxed_;
	QPSchurPack::ConstraintsRelaxedStd
							         constraints_;
	MatrixSymAddDelBunchKaufman      schur_comp_;
	Vector					         g_relaxed_;
	Vector					         b_X_;
	InitKKTSystem::Ko_ptr_t          Ko_;
	Vector					         fo_;
	
	// ////////////////////////////
	// Private member functions

};	// end class QPSolverRelaxedQPSchur

}	// end namespace ConstrainedOptimizationPackTypes

#endif	// QP_SOLVER_RELAXED_QP_SCHUR_H
