// //////////////////////////////////////////////////////////////
// ConstraintsRelaxedStd.h

#ifndef QP_SCHUR_CONSTRAINTS_RELAXED_STD_H
#define QP_SCHUR_CONSTRAINTS_RELAXED_STD_H

#include "QPSchur.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"

namespace ConstrainedOptimizationPack {

namespace QPSchurPack {

///
/** Constraints subclass that is used to represent generic
  * varaible bounds, and general inequality and equality constraints.
  *
  * The generic constraints represented by this class are those
  * of the QPSolverRelaxed interface which are:
  \begin{verbatim}
  
  (1.2)               etaL <=  eta
  (1.3)               dL   <=  d                       <= dU
  (1.4)               eL   <=  op(E)*d + b*eta         <= eU
  (1.5)                        op(F)*d + (1 - eta) * f  = 0
  
  \end{verbatim}
  *
  * These constraints are converted into the form:
  \begin{verbatim}

       [ dL   ]     [ I           ]              [ dU   ]
       [ etaL ] <=  [          1  ] * [  d  ] <= [ inf  ]
  (2)  [ eL   ]     [ op(E)    b  ]   [ eta ]    [ eU   ]
       [ -f   ]     [ op(F)    -f ]              [ -f   ]
       \______/     \_____________/   \_____/    \______/
        cL_bar           A_bar'          x        cU_bar

       =>

  (3)   [     xl     ]    [   I       ]       [     xu     ]
        [  cl_breve  ] <= [  A_breve' ]*x  <= [  cu_breve  ]

       =>

  (4)	cl_bar <= A_bar'*x <= cu_bar

  \end{verbatim}
  * 
  * The main responcibility of this class is to expose a
  * MatrixWithOp object for A_bar shown in (2) and to compute
  * constraint violations.
  */
class ConstraintsRelaxedStd : public Constraints {
public:

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
	enum EInequalityPickPolicy {
		ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY
		,ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY
		,ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY
	};

	///
	/** <<std member comp>> members policy used to select a violated constraint.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EInequalityPickPolicy, inequality_pick_policy )


	/// Constructs to uninitialized
	ConstraintsRelaxedStd();

	///
	/** Initialize constriants.
	  *
	  * If there are no variable bounds then set:
	  * void(dL) == void(dU) == NULL
	  * 
	  * If there are no general inequality constraints
	  * then set:
	  * void(E) == void(b) == void(eL) == void(eU) == NULL
	  * 
	  * If there are no general equality constraints then
	  * set:
	  * void(F) = void(f) == NULL
	  * 
	  * If check_F == false, then the equality constriants
	  * in op(F) will not be checked as violated constriants.
	  * This is to facilitate the addition of the equality
	  * constraints to the initial schur complement and therefore
	  * these constraints should never be violated (except for
	  * illconditioning).
	  * The tolerances below which a constriant will not be considered
	  * violated are given by bounds_tol, inequality_tol and equality_tol.
	  * 
	  * Here, Ed is updated within the function pick_violated(...) when the
	  * options pick_policy == ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY
	  * or pick_policy == ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY
	  * are used.  For these polices it is guarrentteed that
	  * Ed is set to the value op(E)*d on the last call the pick_violated(...)
	  * in which none of the constriants is violated by mor than the
	  * set tolerances and therefore the Primal-Dual algorithm will
	  * terminate immediately.  This is to
	  * save some computational work.  To skip computing this value, just
	  * set Ed == NULL.  On return, if the current value for this->pick_policy()
	  * does not equal any of the two options given above then it should
	  * be assumed that Ed was not update on the last call to
	  * pick_violated(...) and therefore Ed should not be considered current.
	  * 
	  * ToDo: Specify more concretely exactly what the criteria is for
	  * considering that a constraint is violated or in picking the most
	  * violated constraint.
	  */
	void initialize(
		  size_type nd
		, value_type etaL
		, const SpVectorSlice* dL, const SpVectorSlice* dU
		, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
			, const SpVectorSlice* eL, const SpVectorSlice* eU
		, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
		, VectorSlice* Ed
		, bool check_F = true
		, value_type bounds_tol			= 1e-10
		, value_type inequality_tol		= 1e-10
		, value_type equality_tol		= 1e-10
		);

	// /////////////////////////////////////
	// Overridden from Constraints

	///
	size_type n() const;
	///
	size_type m_breve() const;
	///
	/** Represents the constraints matrix.
	  *
	  \begin{verbatim}

		A_bar = [  I   0  op(E')   op(F')  ]
		        [  0   1    b'      -f'    ]

	  \end{verbatim}
	  *
	  */
	const MatrixWithOp& A_bar() const;
	///
	void pick_violated_policy( EPickPolicy pick_policy );
	///
	EPickPolicy pick_violated_policy() const;
	///
	/** Here the next violated constraint to add to the active set is selected.
	  *
	  * Violated constraints are selected to to add to the active set in the following
	  * order.
	  * \begin{itemize}
	  * \item The equality constraints are added first, one at a time (if not already added
	  *		as part of the warm start).
	  * \item Add inequality constraints according according to the following options:
	  *		\begin{itemize}   
	  *		\item ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY
	  *			Check the variable bounds first and add the most violated.  If no
	  *			variable bounds are violated by more than bounds_tol then check for
	  *			the most violated inequality constraint by computing r = op(E)*d+b*eta then
	  *			add the most violated bound (eL, eU).
	  *		\item ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY
	  *			Check the variable bounds first and add the most violated.  If no
	  *			variable bounds are violated by more than bounds_tol then check for
	  *			the first violated inequality constraint by computing e(j)'*(op(E)*d+b*eta)
	  *			one constraint at a time (or a few at a time).  This option may be
	  *			better if the cost of computing op(E)*d is significant.
	  *		\item ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY
	  *			Select the most violated constraint from the variable bounds and the
	  *			general inequality constraints by computing  r = op(E)*d+b*eta then
	  *			add the most violated variable bound.  This option is always the most
	  *			expensive but may result in less QP iterations.
	  * \end{itemize}
	  */
	void pick_violated( const VectorSlice& x, size_type* j_viol, value_type* constr_val
		, value_type* viol_bnd_val, value_type* norm_2_constr, EBounds* bnd, bool* can_ignore
		) const;
	///
	void ignore( size_type j );
	///
	value_type get_bnd( size_type j, EBounds bnd ) const;

protected:

	// /////////////////////////////////////////////
	// Protected types

	///
	/** Matrix type for A_bar.
	  *
	  \begin{verbatim}

		A_bar = [  I   0  op(E')   op(F')  ]
		        [  0   1    b'      -f'    ]

	  \end{verbatim}
	  *
	  */
	class MatrixConstraints : public MatrixWithOp {
	public:

		///
		/** Construct to unitinitialized.
		  *
		  * this->nd() == 0 after construction.
		  *
		  */
		MatrixConstraints();

		///
		/** Initialize.
		  *
		  * The sizes of the arguments are
		  * not checked.
		  * 
		  * It is expected that the objects being
		  * pointed to will not be resized or invalidated
		  * since copies of data are not made!
		  * 
		  */
		void initialize(
			  size_type				nd
			, size_type				m_in
			, size_type				m_eq
			, const MatrixWithOp	*E
			, BLAS_Cpp::Transp		trans_E
			, const VectorSlice		*b
			, const MatrixWithOp	*F
			, BLAS_Cpp::Transp		trans_F
			, const VectorSlice		*f
			);

		// ///////////////////////
		// Access constituents

		///
		size_type			nd() const
		{	return nd_;	}
		///
		size_type			m_in() const
		{	return m_in_;	}
		///
		size_type			m_eq() const
		{	return m_eq_;	}
		///
		const MatrixWithOp*	E() const
		{	return E_;	}
		///
		BLAS_Cpp::Transp	trans_E() const
		{	return trans_E_;	}
		///
		const VectorSlice*	b() const
		{	return b_;	}
		///
		const MatrixWithOp*	F() const
		{	return F_;	}
		///
		BLAS_Cpp::Transp	trans_F() const
		{	return trans_F_;	}
		///
		const VectorSlice*	f() const
		{	return f_;	}


		// ///////////////////////////////
		// Overridden from Matrix

		///
		size_type rows() const;
		///
		size_type cols() const;

		// //////////////////////////////
		// Overridden from MatrixWithOp

		///
		MatrixWithOp& operator=(const MatrixWithOp& m);
//		///
//		void Mp_StPtMtP(GenMatrixSlice* gms_lhs, value_type alpha
//			, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
//			, BLAS_Cpp::Transp M_trans
//			, const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
//			) const ;
		///
		void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
			, const VectorSlice& vs_rhs2, value_type beta) const;
		///
		void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
			, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
			, BLAS_Cpp::Transp M_rhs2_trans
			, const SpVectorSlice& sv_rhs3, value_type beta) const;
		
	private:
		size_type			nd_;	// # unknowns d
		size_type			m_in_;	// # op(E)*d inequality constraints
		size_type			m_eq_;	// # op(F)*d equality constraints
		const MatrixWithOp	*E_;	// If NULL then no general inequalities
		BLAS_Cpp::Transp	trans_E_;
		const VectorSlice	*b_;
		const MatrixWithOp	*F_;	// If NULL then no general equalities
		BLAS_Cpp::Transp	trans_F_;
		const VectorSlice	*f_;
	};	// end class MatrixConstraints

private:
		MatrixConstraints	A_bar_;
		value_type			etaL_;
		const SpVectorSlice	*dL_;	// If NULL then no simple bounds
		const SpVectorSlice *dU_;
		const SpVectorSlice	*eL_;
		const SpVectorSlice	*eU_;
		const VectorSlice	*f_;
		VectorSlice			*Ed_;
		bool				check_F_;
		mutable size_type	last_added_j_;			// Remember the last bound added so that
		mutable value_type	last_added_bound_;		// we can save our selfs some work.
		mutable EBounds		last_added_bound_type_;	//

		///
		void cache_last_added( size_type last_added_j, value_type last_added_bound
			, EBounds last_added_bound_type ) const;


};	// end class ConstraintsRelaxedStd

}	// end namespace QPSchurPack 

}	// end namespace ConstrainedOptimizationPack 

#endif	// QP_SCHUR_CONSTRAINTS_RELAXED_STD_H
