// //////////////////////////////////////////////////////////////
// QPSchur.h

#ifndef QPSCHUR_H
#define QPSCHUR_H

#include <ostream>
#include <map>
#include <vector>

#include "ConstrainedOptimizationPackTypes.h"
#include "MatrixSymAddDelUpdateableWithOpFactorized.h"
#include "SparseLinAlgPack/include/MatrixSymWithOpFactorized.h"
#include "SparseLinAlgPack/include/MatrixSymWithOp.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/GenPermMatrixSlice.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "Misc/include/StandardCompositionMacros.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ConstrainedOptimizationPack {

namespace QPSchurPack {

/// Utility class for a ranged check vector
template < class T >
class vector_one_based_checked : public std::vector<T>
{
public:
	/// one based indexing
	T& operator()( size_type i )
	{
#ifdef LINALGPACK_CHECK_RANGE
			return this->at(i-1);
#else
			return this->operator[](i-1);
#endif
	}
	/// one based indexing
	T operator()( size_type i ) const
	{
#ifdef LINALGPACK_CHECK_RANGE
			return this->at(i-1);
#else
			return this->operator[](i-1);
#endif
	}
}; // end class vector_one_based

class Constraints;
class QP;

///
enum EBounds { FREE, UPPER, LOWER, EQUALITY };

///
/** Represents the QP to be solved by QPSchur {abstract}.
  *
  * In order to solve a QP, clients must define subclasses
  * for this interface and the \Ref{Constraints} interface
  * defined later.  This is where the specialized properties
  * of the QP are exploited.  For certain types of QPs, standard
  * implementation classes can be defined.
  *
  * Here the QP is of the form:
  *
  \begin{verbatim}

	(1.a)	min		g'*x + 1/2*x'*G*x
	(1.b)	s.t.	A'*x = c
	(1.c)			cl_bar <= A_bar'*x <= cu_bar

	where:
		x <: R^n
		g <: R^n
		G <: R^(n x n)
		A <: R^(n x m)
		c <: R^m
		A_bar <: R^(n x m_bar)
		cl_bar, cu_bar <: R^m_bar

  \end{verbatim}
  *
  * Above, #cl_bar <= A_bar'*x <= cu_bar# may represent variable bounds, general
  * inequalities and equality constraints and these constraints are represented
  * by the class \Ref{Constraints}.
  *
  * When solving the above QP using the schur complement QP solver we start out with
  * a KKT system with a subset of variables initially fixed at a bound:
  \begin{verbatim}

  [ G_RR    G_RX     A_F     0  ] [ x_R    ]   [ - g_R ]
  [ G_RX'   G_XX     A_X     I  ] [ x_X    ]   [ - g_X ]
  [ A_R'    A_X'      0      0  ]*[ lambda ] = [   c   ]
  [           I       0      0  ] [ mu_X   ]   [   b_X ]

  \end{verbatim}
  * We can simplify the above system by solving for the initially
  * fixed variables and removing them from the initial KKT system to give:
  \begin{verbatim}

  x_X = b_X

  [ G_RR     A_R  ]   [ x_R    ]   [ - g_R - G_RX*b_X ]
  [ A_R'      0   ] * [ lambda ] = [    c - A_X'*b_X  ]
  \_______________/   \________/   \__________________/
          Ko              vo                fo

  mu_X = - g_X - G_RX'*x_R - G_X*b_X - A_X*lambda

  where:
  		n_X = n - n_R
		x_R  = Q_R'*x        <: R^n_R
		g_R  = Q_R'*g        <: R^n_R
		G_RR = Q_R'*G*Q_R    <: R^(n_R x n_R)
		G_RX = Q_R'*G*Q_X    <: R^(n_R x n_X)
		G_XX = Q_X'*G*Q_X    <: R^(n_X x n_X)
		A_R  = Q_R'*A        <: R^(n_R x m) 
		A_X  = Q_X'*A        <: R^(n_X x m)
		Q_R                  <: R^(n x n_R)
		Q_X                  <: R^(n x n_X)

  \end{verbatim}
  * This class is an interface for encapsulating this QP.  Operations are available
  * for accessing #g#, #G#, #A#, #Ko#, #vo#, and #fo# as well as the mapping
  * matrices #Q_R# and #Q_X# (both ordered by row). Also, operations are available
  * for accessing data structures that describe the set of initially fixed and free
  * variables.  See the \Ref{Constraints} interface for how to access the constraints
  * in (1.c) and the matrix #A_bar#.
  */
class QP {
public:

	// /////////////////
	// Public Types

	///
	typedef vector_one_based_checked<EBounds>		x_init_t;
	///
	typedef vector_one_based_checked<size_type>		l_x_X_map_t;
	///
	typedef vector_one_based_checked<size_type>		i_x_X_map_t;

	///
	typedef QPSchurPack::Constraints				Constraints;

	// /////////////////
	// Public Interface

	///
	virtual ~QP()
	{}

	// ///////////////////////////////////////
	// Initial active set independent members 

	///
	virtual size_type n() const = 0;
	///
	virtual size_type m() const = 0;
	///
	virtual const VectorSlice g() const = 0;
	///
	virtual const MatrixSymWithOp& G() const = 0;
	/// If m == 0 then don't call this
	virtual const MatrixWithOp& A() const = 0;

	// /////////////////////////////////////
	// Initial active set specific members

	///
	virtual size_type n_R() const = 0;

	///
	/** Return the status of a variable initially.
	  *
	  * For 1 <= i <= n:
	  *
	\begin{verbatim}
	             / FREE      : x(i) is initially free
	             | LOWER     : x(i) is initially fixed at xl(i)
	 x_init(i) = | UPPER     : x(i) is initially fixed at xu(i)
	             \ EQUALITY  : x(i) fixed at xl(i) = xu(i) always
	\end{verbatim}
	*
	*/
	virtual const x_init_t& x_init() const = 0;

	///
	/** Map from full x(i) to initially fixed x_X(l).
	  *
	  * For 1 <= i <= n:
	  * 
	\begin{verbatim}
	                 / l : x(i) = x_X(l) = b_X(l) initially (1 <= l <= n_X)
	 l_x_X_map(i) =  |
	                 \ 0 : otherwise
	\end{verbatim}
	*
	*/
	virtual const l_x_X_map_t& l_x_X_map() const = 0;

	///
	/** Map from initially fixed x_X(l) to full x(i).
	  *
	  * For 1 <= l <= n_X:
	  * 
	\begin{verbatim}
	 i_x_X_map(l) = i : x(i) = x_X(l) = b_X(l) initially (1 <= i <= n)
	\end{verbatim}
	*
	*/
	virtual const i_x_X_map_t& i_x_X_map() const = 0;

	///
	/** The bounds of the initially fixed variables.
	  *
	  * For 1 <= l <= n_X:
	  *
	\begin{verbatim}
	           / xl(i_x_X_map(l))                     : if x_init(i_x_X_map(l)) == LOWER
	 b_X(l) =  | xu(i_x_X_map(l))                     : if x_init(i_x_X_map(l)) == UPPER
	           \ xl(i_x_X_map(l)) = xu(i_x_X_map(l))  : if x_init(i_x_X_map(l)) == EQUALITY
	\end{verbatim}
	*
	*/
	virtual const VectorSlice b_X() const = 0;

	/// (Q_R().ordered_by() == BY_ROW)
	virtual const GenPermMatrixSlice& Q_R() const = 0;

	/// (Q_X().ordered_by() == BY_ROW)
	virtual const GenPermMatrixSlice& Q_X() const = 0;

	///
	virtual const MatrixSymWithOpFactorized& Ko() const = 0;

	///
	virtual const VectorSlice fo() const = 0;

	// //////////////////////////////////////////////////////////
	// Additional constaints for cl_bar <= A_bar'*x <= cu_bar

	///
	virtual Constraints& constraints() = 0;

	///
	virtual const Constraints& constraints() const = 0;

	///
	/** Dump the definition of the QP to a stream.
	  *
	  * This function is only to be used for debugging small problems.
	  */
	virtual void dump_qp( std::ostream& out );

};	// end class QP

///
/** Represents the extra constraints in the QP to be satisfied
  * by the schur complement QP solver QPSchur {abstract}.
  *
  * This class is only ment to be used in conjunction with the class \Ref{QP}
  * and \Ref{QPSchur}.  Its interface is designed to be minimal with respect to
  * the needs of the #QPSchur# solver.  However, this interface may be useful
  * for any primal-dual QP solver.
  *
  * These constraints are:
  \begin{verbatim}

	(1.c)	cl_bar <= A_bar'*x <= cu_bar

	where:
		A_bar <: R^(n x m_bar)
		cl_bar, cu_bar <: R^m_bar

  \end{verbatim}
  *
  * These constraints are also partitioned as:
  \begin{verbatim}

	s.t.
		[     xl     ]    [   I       ]       [     xu     ]
		[  cl_breve  ] <= [  A_breve' ]*x  <= [  cu_breve  ]

	where:
		I <: R^(n x n)
		xl, xu <: R^n, are the variable bounds for variables that have bounds (sparse)
		A_breve <: R^(n x m_breve), is the Jacobian for the general constraints
		cl_breve, cu_breve <: R^m_breve, are bounds for general constraints

  \end{verbatim}
  *
  * Here #m_bar = n + m_breve#
  *
  * Above, some of the bounds in #xl#, #xu#, #cl_breve#, and #cu_breve# may be #-inf#
  * or #+inf# and will therefore never be violated and never be added to the active set.
  * Also, some of the lower and upper bounds may be equal which turns those
  * inequality constraints into equality constraints (or fixed variables).
  */
class Constraints {
public:

	///
	enum EPickPolicy { ANY_VIOLATED, MOST_VIOLATED };

	///
	virtual ~Constraints()
	{}

	///
	virtual size_type n() const = 0;
	
	///
	virtual size_type m_breve() const = 0;

	///
	virtual const MatrixWithOp& A_bar() const = 0;
	
	/// Set the policy used to pick a violated constraint.
	virtual void pick_violated_policy( EPickPolicy pick_policy ) = 0;
	///
	virtual EPickPolicy pick_violated_policy() const = 0;

	///
	/** Pick a violated constraint.
	  *
	  * @param	x			 [in] Trial point to pick a violated constraint at.
	  * @param	j_viol		 [out] Indice of violated constraint.  j_viol = 0 if
	  *							 no constraint is violated by more that some tolerance.
	  * @param	constr_val	 [out] The value if the violated constraint a_bar(j)'*x.
	  * @param	viol_bnd_val [out] The value if the violated bound.
	  * @param	norm_2_constr[out] The 2 norm of the violated constraint ||a_bar(j)||2
	  * @param	bnd			 [out] Classification of the bound being violated.
	  * @param	can_ignore	 [out] True if the constraint can be ignored if it is linearly
	  *							 dependent.
	  */
	virtual void pick_violated( const VectorSlice& x, size_type* j_viol, value_type* constr_val
		, value_type* viol_bnd_val, value_type* norm_2_constr, EBounds* bnd, bool* can_ignore
		) const = 0;

	///
	/** Inform to ignore the jth constraint the next time pick_violated(...) is called.
	  */
	virtual void ignore( size_type j ) = 0;

	///
	/** Return the bound for a constraint.
	  *
	  * @param	j	[in] Indice of the constraint of the bound to obtain.
	  * @param	bnd	[in] Which bound to obtain (UPPER or LOWER).
	  * @return
	  *		xl(j) [ 0 < j < n, bnd == LOWER ]\\ 
	  *		xu(j) [ 0 < j < n, bnd == UPPER ]\\ 
	  *		cl_breve(j - n) [ n + 1 < j < n + m_breve, bnd == LOWER ]\\ 
	  *		cu_breve(j - n) [ n + 1 < j < n + m_breve, bnd == UPPER ]\\ 
	  */
	virtual value_type get_bnd( size_type j, EBounds bnd ) const = 0;

};	// end class Constraints

}	// end namespace QPSchurPack 

///
/** Solves a Quadratic Program with a primal-dual QP method
  * using a schur complement.
  *
  * See the paper "QPSchur: A Primal-Dual Active-Set Quadratic Programming
  * Algorithm Using a Schur Complement Factorization Method" for a description
  * of what this class does.
  */
class QPSchur {
public:

	// ////////////////////////
	/** @name Public Types */
	//@{

	///
	typedef QPSchurPack::QP		QP;
	/// Thrown if a test failed
	class TestFailed : public std::logic_error
	{public: TestFailed(const std::string& what_arg) : std::logic_error(what_arg) {}};
	/// Thrown if constraints are inconsistant (no feasible region)
	class InconsistantConstraintsException : public std::logic_error
	{public: InconsistantConstraintsException(const std::string& what_arg) : std::logic_error(what_arg) {}};
	/// Thrown if there is some numerical instability
	class NumericalInstabilityException : public std::runtime_error
	{public: NumericalInstabilityException(const std::string& what_arg) : std::runtime_error(what_arg) {}};
	/// Thrown if during the course of the primal-dual iteration a non-dual feasible point if found.
	class DualInfeasibleException : public NumericalInstabilityException
	{public: DualInfeasibleException(const std::string& what_arg)
		: NumericalInstabilityException(what_arg) {}};
	/// Enumeration for if to run internal tests or not.
	enum ERunTests { RUN_TESTS, NO_TESTS };
	/// solve_qp return values
	enum ESolveReturn {
		 OPTIMAL_SOLUTION
		,MAX_ITER_EXCEEDED
		,MAX_ALLOWED_STORAGE_EXCEEDED
		,INFEASIBLE_CONSTRAINTS
		,DUAL_INFEASIBILITY
		,SUBOPTIMAL_POINT
	};
	/// Output level
	enum EOutputLevel {
		 NO_OUTPUT					= 0
		,OUTPUT_BASIC_INFO			= 1
		,OUTPUT_ITER_SUMMARY		= 2
		,OUTPUT_ITER_STEPS			= 3
		,OUTPUT_ACT_SET				= 4
		,OUTPUT_ITER_QUANTITIES		= 5
	};
	/// Value for near degenerate lagrange multipliers
	static value_type DEGENERATE_MULT;

	//		end Public Types
	//@}

	// //////////////////////////////////////
	/** @name Public Member functions */
	//@{

	/// «std comp» members for schur complement matrix object S_hat
	STANDARD_COMPOSITION_MEMBERS( MatrixSymAddDelUpdateableWithOpFactorized, schur_comp )

	///
	/** Set the maximum number of primal-dual QP iterations to take.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, max_iter )

	///
	/** Set the feasibility tolerance for the constriants.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, feas_tol )

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
	QPSchur(
		  const schur_comp_ptr_t& 	schur_comp 			= NULL
		, size_type					max_iter			= 100
		, value_type				feas_tol			= 1e-8
		, value_type				loose_feas_tol		= 1e-6
		, value_type				dual_infeas_tol		= 1e-12
		, value_type				huge_primal_step	= 1e+20
		, value_type				huge_dual_step		= 1e+20
		);

	///
	/** Solve a QP.
	  *
	  * @param	qp	[in] The abstraction for the QP being solved
	  * @param	num_act_change
	  * 			[in] The number of changes to the
	  *					active set before the primal-dual QP algorithm
	  *					starts.
	  * @param	ij_act_change
	  * 			[in] Array (size num_act_change): specifying
	  *					how to initialize the active set.  If i = ij_act_change(s)
	  *					< 0 then the initially fixed variable x(-i) is to be
	  *					freed.  If j = ij_act_change(s) > 0 then the constraint
	  *					a_bar(j)'*x is to be added to the active set.  The order
	  *					of these changes can significantly effect the performance
	  *					of the algorithm if these change.
	  *	@param	bnd [in] Array (size num_act_change):  bnd(s) gives which bound to
	  *					make active.  If ij_act_change(s) < 0 then this is ignored.
	  * @param	out [out] output stream.  Iteration information is printed according
	  *					to output_level.  If #output_level == NO_OUTPUT# then #out# may
	  *					be #NULL#.  If #out==NULL#, then output_level is forced to #NO_OUTPUT#
	  * @param	output_level
	  * 			[in] Specifies the level of output (see \Ref{EOutputLevel}).
	  *	@param	test_what
	  *				[in] Determines if internal validation tests are performed.
	  *					The optimality conditions for the QP are not checked
	  *					internally, since this is something that client can
	  *					(and should) do independently.
	  *					RUN_TESTS : As many validation/consistency tests
	  *						are performed internally as possible.  If a test
	  *						fails then a TestFailed execption will be thrown.
	  *					NO_TEST : No tests are performed internally.  This is
	  *						to allow the fastest possible execution.
	  * @param	x	[out] vector (size qp.n()): Solution or current iteration value
	  * @param	mu 	[out] sparse vector (size qp.n()): Optimal lagrange multipliers for
	  * 				bound constraints.  On output mu->is_sorted() == true.
	  * @param	lambda
	  * 			[out] vector (size q.m()): Optimal lagrange multipliers for
	  *					equality constraints.
	  * @param	lambda_breve
	  * 			[out] sparse vector (size qp.constraints().m_breve()) for the active
	  *					constraints in A_breve.  On output lambda_breve->is_sorted() == true.
	  *	@param	iter [out] The number of warm start drops and primal-dual iterations.
	  *	@param	num_adds [out] The number of updates to the active set where a constraint
	  *					was added.  These do not include initially fixed variables.
	  *	@param	num_drops [out] The number of updates to the active set where a constraint
	  *					was dropped.  These include constraints dropped during a warm start
	  *					as well as during the primal-dual QP iterations.
	  */
	virtual
	ESolveReturn solve_qp(
		  QP& qp
		, size_type num_act_change, const int ij_act_change[]
			, const QPSchurPack::EBounds bnds[]
		, std::ostream *out, EOutputLevel output_level, ERunTests test_what
		, VectorSlice* x, SpVector* mu, VectorSlice* lambda, SpVector* lambda_breve
		, size_type* iter, size_type* num_adds, size_type* num_drops
		);

	//		end Public Member functions
	//@}

	///
	/** Represents the matrix U_hat. 
	  *
	  * This matrix is only ment to be an aggregate of and ActiveSet
	  * object and is only managed by the ActiveSet object.  It is made
	  * public so that clients can developed specialized implementations
	  * if needed.
	  */
	class U_hat_t : public MatrixWithOp {
	public:
		/// Construct uninitialized
		U_hat_t();
		/// Initialize.
		void initialize( 
			 const MatrixSymWithOp		*G
			,const MatrixWithOp			*A
			,const MatrixWithOp			*A_bar
			,const GenPermMatrixSlice	*Q_R
			,const GenPermMatrixSlice	*P_XF_hat
			,const GenPermMatrixSlice	*P_plus_hat
			);
		///
		const MatrixSymWithOp& G() const
		{	return *G_;	}
		///
		const MatrixWithOp* A() const
		{	return A_;	}
		///
		const MatrixWithOp& A_bar() const
		{	return *A_bar_;	}
		///
		const GenPermMatrixSlice& Q_R() const
		{	return *Q_R_; }
		///
		const GenPermMatrixSlice& P_XF_hat() const
		{	return *P_XF_hat_;	}
		///
		const GenPermMatrixSlice& P_plus_hat() const
		{	return *P_plus_hat_;	}

		// /////////////////////////////////////////////////////
		// Overridden from Matrix

		///
		size_type rows() const;
		///
		size_type cols() const;

		// /////////////////////////////////////////////////////////
		// Overridden from MatrixWithOp

//		///
//		void Mp_StM(GenMatrixSlice* gms_lhs, value_type alpha
//			, BLAS_Cpp::Transp trans_rhs) const ;

		///
		void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
			, const VectorSlice& vs_rhs2, value_type beta) const;

		///
		void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
			, const SpVectorSlice& sv_rhs2, value_type beta) const;
		
	private:
		const MatrixSymWithOp		*G_;
		const MatrixWithOp			*A_;
		const MatrixWithOp			*A_bar_;
		const GenPermMatrixSlice	*Q_R_;
		const GenPermMatrixSlice	*P_XF_hat_;
		const GenPermMatrixSlice	*P_plus_hat_;

	};	// end class U_hat_t

	///
	/** Represents and manages the active set for the QPSchur algorithm.
	  *
	  * Concrete type that encapsulates the maintaince of the active set and
	  * abstracts quantities associated with it.
	  *
	  * At each iteration the algorithm must solve the system:
	  \begin{verbatim}

		[ Ko       U_hat ]   [    v  ]   [  fo   ]
		[ U_hat'   V_hat ] * [ z_hat ] = [ d_hat ]

	  \end{verbatim}
	  *
	  * Above, U_hat contains the updates to the KKT system for adding constraints
	  * to the active set and freeing variables that where initially fixed
	  * and therefore left out of Ko.
	  * 
	  * This object maintains references to objects that represent the current
	  * augmented KKT system:
	  *
	  * MatrixWithOp               : U_hat        ( \hat{U} )
	  * MatrixSymWithOP            : V_hat        ( \hat{V} )
	  * MatrixSymWithOpFactorized  : S_hat        ( \hat{S} )
	  * GenPermMatrixSlice         : P_XF_hat     ( \hat{P}^{XF} )
	  * GenPermMatrixSlice         : P_plus_hat   ( \hat{P}^{(+)} )
	  * GenPermMatrixSlice         : Q_XD_hat     ( \hat{Q}^{XD}  )
	  * Vector                     : d_hat        ( \hat{d} )
	  * Vector                     : z_hat        ( \hat{z} )
	  * 
	  */
	class ActiveSet {
	public:

		// /////////////////////
		// Public types

		///
		typedef QPSchurPack::QP																	QP;
		///
		typedef QPSchurPack::EBounds																EBounds;

		/// Thrown if the update failed
		class BadUpdateException : public std::logic_error
		{public: BadUpdateException(const std::string& what_arg) : std::logic_error(what_arg) {}};

		/// Thrown if you attempt to add a constraint to the active set that is linearly dependent.
		class LDConstraintException : public BadUpdateException
		{public: LDConstraintException(const std::string& what_arg) : BadUpdateException(what_arg) {}};

		// /////////////////////
		// Public interface

		///
		/** «std comp» members for schur complement matrix S_hat.
		  *
		  * Warning: Resetting schur_comp will cause a reinitialization to
		  * an empty active set.
		  */
		STANDARD_COMPOSITION_MEMBERS( MatrixSymAddDelUpdateableWithOpFactorized, schur_comp )

		///
		ActiveSet(const schur_comp_ptr_t& schur_comp);

		// ///////////////////////////////////////
		/** @name Update the active set. */
		//@{

		///
		/** Initialize with an additional active set.
		  *
		  * If the initial schur complement is not full rank
		  * then an #LDConstraintException# exception will be thrown.
		  * The active set will contain all of the constraints it
		  * can such that the schur complement is nonsingular.
		  */
		void initialize( QP& qp, size_type num_act_change, const int ij_act_change[]
		, const QPSchurPack::EBounds bnds[], bool test
		, std::ostream *out, EOutputLevel output_level );

		///
		/** Reinitialize the schur complement factorization for the current active set
		  *
		  * If the freashly computed schur complement is still nonsingular (or indefinite)
		  * then an #LDConstraintException# exception will be thrown.
		  */
		void refactorize_schur_comp();

		///
		/** Add a constraint to the active set then refactorize the schur complemnt
		  * (if forced).
		  *
		  * ToDo: Finish documentation
		  *
		  * If the constraint is linearly dependent this function will throw
		  * the exception #LDConstraintException# and the active set will not
		  * be updated.
		  */
		void add_constraint( size_type ja, QPSchurPack::EBounds bnd_ja
			, bool update_steps, bool force_refactorization = true );

		///
		/** Drop a constraint from the active set then refactorize the schur
		  * complement (if forced).
		  *
		  * ToDo: Finish documentation
		  */
		void drop_constraint( size_type jd , bool force_refactorization = true );

		///
		/** Drop a constraint from, then add a constraint to the active set
		  * and refactorize the schur complement.
		  *
		  * ToDo: Finish documentation
		  */
		void drop_add_constraints( size_type jd, size_type ja, QPSchurPack::EBounds bnd_ja
			, bool update_steps );

		//@}

		///
		QP& qp();
		///
		const QP& qp() const;

		// ///////////////////////////////////////////////
		/** @name Access the active sets quantities. */
		//@{

		///
		/** Return the total size of the schur complement.
		  *
		  * q_hat = q_plus_hat + q_F_hat + q_C_hat.
		  */
		size_type q_hat() const;

		///
		/** Return the number of constraints from A_bar added
		  * to the active set.
		  */
		size_type q_plus_hat() const;

		///
		/** Return the number of variables that where
		  * initially fixed but are currently free or
		  * fixed to another bound.
		  */
		size_type q_F_hat() const;

		///
		/** Return the number of variables that where
		  * initially fixed but are currently
		  * fixed to another bound.
		  */
		size_type q_C_hat() const;

		///
		/** Return the number of variables that where
		  * initially fixed and are still currently
		  * fixed to their intial bounds.
		  */
		size_type q_D_hat() const;

		///
		/** Returns -i for row & column of S_bar for an initially
		  * fixed variable left out of Ko that became free and returns
		  * j for the constraint a(j)'*x that was added to the active
		  * set.
		  *
		  * 1 <= s <= q_hat
		  */
		int ij_map( size_type s ) const;

		///
		/** Map from a constraint or initially fixed variable
		  * to a row and column in the schur complement S_bar.
		  *
		  * To determine if an initially fixed varible x(i) is now
		  * free call s_map(-i).  If s_map(-i) returns zero then
		  * x(i) is still fixed.  Otherwise s_map(-i) returns the
		  * row and column in S_bar for this change in the
		  * active set.
		  *
		  * To determine if a constraint a(j)'*x is part of the
		  * active set call s_map(j).  If s_map(j) returns zero
		  * then a(j)'*x is not part of the active set.
		  * Otherwise s_map(j) returns the row and column
		  * in S_bar for this change in the active set.
		  */
		size_type s_map( int ij ) const;

		///
		/** Returns ||a(j)||2 where j = ij_map(s).
		  * 
		  * If ij_map(s) < 0, the this function returns zero.
		  * 
		  * 1 <= s <= q_hat
		  */
		value_type constr_norm( size_type s ) const;

		///
		/** Return which bound is active for the active constraint.
		  */
		QPSchurPack::EBounds bnd( size_type s ) const;

		///
		const U_hat_t& U_hat() const;
		///
		const MatrixSymWithOpFactorized& S_hat() const;
		///
		const GenPermMatrixSlice& P_XF_hat() const;
		///
		const GenPermMatrixSlice& P_plus_hat() const;
		///
		const GenPermMatrixSlice& Q_XD_hat() const;
		///
		const VectorSlice d_hat() const;
		///
		VectorSlice z_hat();
		///
		const VectorSlice z_hat() const;
		///
		VectorSlice p_z_hat();
		///
		const VectorSlice p_z_hat() const;
		///
		VectorSlice mu_D_hat();
		///
		const VectorSlice mu_D_hat() const;
		///
		VectorSlice p_mu_D_hat();
		///
		const VectorSlice p_mu_D_hat() const;

		///
		/** Determine if a constriant was an initially fixed variable.
		  *
		  * This function will return true if:
		  * 
		  * j <= n && x_init(j) != FREE
		  * 
		  * This is just a function of convienience
		  * 
		  */
		bool is_init_fixed( size_type j ) const;

		/// Returns true if all the degrees of freedom of the QP are used up
		bool all_dof_used_up() const;

		//@}

	private:

		// ///////////////////////////
		// Private types

		///
		typedef std::vector<int>			ij_map_t;
		///
		typedef std::map<int,size_type>		s_map_t;
		///
		typedef std::vector<EBounds>		bnds_t;
		///
		typedef std::vector<size_type>		P_row_t;
		///
		typedef std::vector<size_type>		P_col_t;

		// ///////////////////////////
		// Private data members

		bool				initialized_;
		bool				test_;
		QP*					qp_;	// QP being solved.
		const QP::x_init_t	*x_init_;
		size_type			n_;
		size_type			n_R_;
		size_type			m_;
		size_type			m_breve_;
		size_type			q_plus_hat_;
		size_type			q_F_hat_;
		size_type			q_C_hat_;
		ij_map_t			ij_map_;
//		s_map_t				s_map_;
		Vector				constr_norm_;
		bnds_t				bnds_;
		U_hat_t				U_hat_;
		//
		// for s = 1...q_hat
		//
		//                     /  e(i)    if i > 0 (where: i = -ij_map(s))
		// [P_XF_hat](:,s)   = |
		//                     \  0       otherwise
		//                     
		GenPermMatrixSlice	P_XF_hat_;		// \hat{P}^{XF}
		P_row_t				P_XF_hat_row_;	// i
		P_row_t				P_XF_hat_col_;	// s
		//
		// for s = 1...q_hat
		//
		//                     /  e(j)    if j > 0 && !is_init_fixed(j) (where: j = ij_map(s))
		// [P_plus_hat](:,s) = |
		//                     \  0       otherwise
		//
		GenPermMatrixSlice	P_plus_hat_;	// \hat{P}^{(+)}
		P_row_t				P_plus_hat_row_;	// j
		P_row_t				P_plus_hat_col_;	// s
		//
		// for k = 1...q_D_hat
		//
		// [Q_XD_hat](:,k) = e(i)  (where is_init_fixed(i) && s_map(-i) == 0)
		//
		GenPermMatrixSlice	Q_XD_hat_;		// \hat{Q}^{XD}
		P_row_t				Q_XD_hat_row_;	// i
		P_row_t				Q_XD_hat_col_;	// k
		Vector				d_hat_;			// \hat{d}		
		Vector				z_hat_;			// \hat{z}
		Vector				p_z_hat_;
		Vector				mu_D_hat_;		// \hat{\mu}^{D}
		Vector				p_mu_D_hat_;	// p^{\hat{\mu}^{D}}

		// ///////////////////////////
		// Private member functions

		//
		void assert_initialized() const;

		// Assert in range.
		void assert_s( size_type s) const;

		// Reinitialize P_XF_hat, P_plus_hat, Q_XD_hat, and U_hat
		void reinitialize_matrices(bool test);

		// not defined and not to be called.
		ActiveSet();

	};	// end class ActiveSet

	/// Return a reference to the active set object
	const ActiveSet& act_set() const;

	/// Dump all the active set quantities for debugging
	static void dump_act_set_quantities( const ActiveSet& act_set, std::ostream& out
		, bool print_S_hat = true );

protected:

	// /////////////////////////
	// Protected types

	///
	enum EPDSteps { PICK_VIOLATED_CONSTRAINT, UPDATE_ACTIVE_SET, COMPUTE_SEARCH_DIRECTION
		, COMPUTE_STEP_LENGTHS, TAKE_STEP };

	// ///////////////////////////
	// Protected Member functions

	///
	/** Run the algorithm from a dual feasible point.
	  *
	  * By default, the algorithm should start with
	  * first_step = PICK_VIOLATED_CONSTRAINT if we are starting
	  * with a dual feasible point.
	  */
	virtual
	ESolveReturn qp_algo(
		  EPDSteps first_step
		, std::ostream *out, EOutputLevel output_level, ERunTests test_what
		, const VectorSlice& vo, ActiveSet* act_set, VectorSlice* v
		, VectorSlice* x, size_type* iter, size_type* num_adds, size_type* num_drops
		);

	///
	/** Set the values in x for all the variables.
	  *
	  * @param	set_fixed	[in] If true then those variables that where initially
	  *							free are specifically fixed to their bounds.
	  */
	virtual void set_x( const ActiveSet& act_set, const VectorSlice& v, VectorSlice* x );

	/// Map from the active set to the sparse multipliers for the inequality constraints
	virtual void set_multipliers( const ActiveSet& act_set, const VectorSlice& v
		, SpVector* mu, VectorSlice* lambda, SpVector* lambda_breve );

private:

	// /////////////////////////
	// Private data members

	ActiveSet		act_set_;	// The active set.

	// The rest of these are just workspace so as to reduce allocations.
	Vector			vo_,			// vo = inv(Ko) * fo
					v_,				// v = [ x_R; lambda ]
					z_hat_plus_,	// z_hat for if ja is added to active set.
					v_plus_,		// v = [ x; lambda ] if ja is added to active set.
					p_z_hat_;		// Search direction for z_hat for adding ja to active set.

	// /////////////////////////
	// Private Member functions.

};	// end class QPSchur

}	// end namespace ConstrainedOptimizationPack 

#endif	// QPSCHUR_H
