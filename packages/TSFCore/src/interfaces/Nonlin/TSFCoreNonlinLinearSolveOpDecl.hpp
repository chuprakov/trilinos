// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreNonlinLinearSolveOpDecl.hpp

#ifndef TSFCORE_NONLIN_LINEAR_SOLVE_OP_DECL_HPP
#define TSFCORE_NONLIN_LINEAR_SOLVE_OP_DECL_HPP

#include "TSFCoreNonlinTypes.hpp"
#include "TSFCoreOpBase.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Base class the implements solves with a linear operator.
 *
 * Note that this is not a mathematical abstraction in itself but is
 * an implementation artifact.
 *
 * This interface allows the solution of the following linear systems:
 *
 * <ul>
 * <li>Compute <tt>x</tt> such that <tt>op(M)*x = y</tt> is sufficiently solved
 * <li>Compute <tt>X</tt> such that <tt>1/alpha*op(M)*X(j) - Y(j)</tt>, for <tt>j=1..Y.domain()->dim()</tt>
 *     is sufficiently solved.
 * </ul>
 *
 * through the <tt>solve()</tt> methods where <tt>x</tt> and
 * <tt>y</tt> are <tt>Vector</tt> objects while <tt>X</tt> and
 * <tt>Y</tt> are <tt>MultiVector</tt> objects.  Above, the definition
 * of "sufficiently solved" is determined by a <tt>Solvers::ConvergenceTester</tt> object.
 * The reason for the exact form of the above conditions
 * is that there are direct BLAS and LAPACK equivalent versions
 * of these operations (for direct solvers).
 *
 * Note that it is not required that <tt>*this</tt> operator be
 * nonsingular, only that the right-hand-sides <tt>y</tt> and
 * <tt>Y</tt> be consistent with <tt>*this</tt> operator object.
 * The exact nonsingular status of <tt>*this</tt> operator is
 * returned by the method <tt>nonsingStatus()</tt>.  Note that
 * if <tt>this->nonsingStatus()==OP_SINGULAR</tt> then the
 * right-hand-sides <tt>y</tt> and <tt>Y</tt> had better be
 * consistent with <tt>*this</tt> or a <tt>Solvers::Exceptions::SolverBreakdown</tt>
 * exception will be thown.
 *
 * Note that it is strictly forbidden to alias the input/output
 * objects <tt>x</tt> and <tt>X</tt> with the input objects <tt>y</tt>
 * and <tt>Y</tt>.
 *
 * If a <tt>%LinearSolveOp</tt> subclass can not support a
 * particular value of <tt>M_tans</tt> in the <tt>solve()</tt>
 * methods, then the method <tt>opSupported()</tt> returns
 * <tt>false</tt> for that particular value of <tt>M_trans</tt>.
 *
 * <b>Notes for subclass developers</tt>
 *
 * A subclass must only override the single-vector version of
 * <tt>solve()</tt>.
 *
 * If a <tt>LinearSolveOp</tt> subclass can not support the
 * transposed argument of <tt>solve()</tt>
 * (<tt>M_trans==TRANS</tt>) then the method <tt>opSupported()</tt>
 * must be overriden to return <tt>false</tt>.
 *
 * If possible, the subclass should also override the <tt>clone_ilo()</tt>
 * method with allows clients to create copies of a <tt>LinearSolveOp</tt>
 * object.  This functionality is very important in some
 * circumstances.  However, this functionality is not required and
 * <tt>clone_ilo()</tt> returns a null smart pointer object.
 *
 * If multi-vectors are supported in general by the application and
 * linear algebra library then, if possible, the subclass should also
 * override the multi-vector version of <tt>solve()</tt>.  In many
 * cases, a specialized multi-vector version will outperform the
 * default implementation (which is based on the single vector
 * version) in this class.
 */
template<class Scalar>
class LinearSolveOp : virtual public OpBase<Scalar> {
public:

	/** @name Pure virtual methods (must be overridden by subclass) */
	//@{

	///
	/** Solve a linear system with this linear operator (or its transpose):
	 * Compute <tt>x</tt> such that <tt>norm(op(M)*x - y) <= tol</tt>
	 *
	 * @param  M_trans
	 *                [in] Determines whether the transposed or non-trnasposed
	 *                operator is applied as:
	 *                <ul>
	 *                <li> <tt>op(M) = M</tt>, for <tt>M_trans==NOTRANS</tt>
	 *                <li> <tt>op(M) = M'</tt>, for <tt>M_trans==TRANS</tt>
	 *                </ul>
	 *                where <tt>M == *this</tt>
	 * @param  y      [in] The right hand side vector.
	 * @param  x      [in/out] In input, <tt>x</tt> should contain a guess for the solution
	 *                of the linear system (by default you should use 0).  In output, contains
	 *                the approximate solution to the linear system as defined above.
	 * @param  convTester
	 *                [in] Determines the convergence criteria for the linear solver.
	 * 
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>this->opSupported(M_trans) (throw <tt>Exceptions::OpNotSupported</tt>)
	 * <li> <tt>y.space()->isCompatible(M_trans==NOTRANS ? *this->range() : *this->domain()) == true</tt>
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>x->space()->isCompatible(M_trans==NOTRANS ? *this->domain() : *this->range()) == true</tt>
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>x</tt> can not alias <tt>y</tt>.  It is up to the client to ensure that <tt>x</tt>
	 *      and <tt>y</tt> are distinct since in general this can not be verified by the implementation until,
	 *      perhaps, it is too late.  If possible, an exception will be thrown if aliasing is detected.
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>convTester==NULL</tt>] If the solver implementation has determined that a suitable
	 *      solution has been found, this method will return silently.
	 * <li> [<tt>convTester!=NULL</tt>] A solution will have been considered to have been found and this
	 *      method will return silently if <tt>convTester->convStatus(...,isConverged[])</tt> returned
	 *      <tt>isConverged[0] == true</tt> from at least one call the method <tt>convStatus(...)</tt>.
	 * <li> If <tt>this->nonsingStatus()==OP_SINGULAR</tt> and the right-hand-side <tt>y</tt> is not
	 *      consistent with <tt>*this</tt>, then a <tt>Solvers::Exceptions::SolverBreakdown</tt> exception will be thown.
	 * <li> If a solution could not be found bacause the linear operator was nearly singular,
	 *      then a <tt>Solvers::Exceptions::SolverBreakdown</tt> exception will be thrown.
	 * <li> If a solution could not be found for some other reason (i.e. the <tt>*convTester</tt> object
	 *      did not state that the system was solved), then a <tt>Solvers::Exceptions::FailureToConverge</tt> exception will be thrown.
	 * </ul>
	 */
	virtual void solve(
		const ETransp                        M_trans
		,const Vector<Scalar>                &y
		,Vector<Scalar>                      *x
		,Solvers::ConvergenceTester<Scalar>  *convTester = NULL
		) const = 0;

	//@}

	/** @name Virtual functions with default implemenations */
	//@{

	///
	/** Returns the status of the opeator.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * </ul>
	 *
	 * The default implementation returns <tt>OP_SINGULARITY_UNKNOWN</tt>
	 */
	virtual ENonsingStatus nonsingStatus() const;

	///
	/** Clone this object (if supported).
	 *
	 * The primary purpose for this method is to allow a client to
	 * capture the current state of this object and be guaranteed that
	 * some other client will not alter its behavior.  A smart
	 * implementation will use reference counting and lazy evaluation
	 * internally and will not actually copy any large amount of data
	 * unless it has to.
	 *
	 * The default implementation returns <tt>return.get()==NULL</tt>
	 * which is allowable by this specification.  An object is not
	 * required to return a non-NULL value but almost every good
	 * implementation should and will.
	 */
	virtual MemMngPack::ref_count_ptr<const LinearSolveOp<Scalar> > clone_lso() const;

	///
	/** Solve a set of linear systems with this linear operator (or its transpose):
	 * Compute <tt>X</tt> such that <tt>norm(1/alpha*op(M)*X(j) - Y(j)) <= tols(j)</tt>, for <tt>j=1..Y.domain()->dim()</tt>
	 *
	 * @param  alpha  [in] Scalar that multiplies <tt>M</tt> where <tt>M == *this</tt>
	 * @param  M_trans
	 *                [in] Determines whether the transposed or non-trnasposed
	 *                operator is applied as:
	 *                <ul>
	 *                <li> <tt>op(M) = M</tt>, for <tt>M_trans==NOTRANS</tt>
	 *                <li> <tt>op(M) = M'</tt>, for <tt>M_trans==TRANS</tt>
	 *                </ul>
	 *                where <tt>M == *this</tt>
	 * @param  Y      [in] The right hand side multi-vector.
	 * @param  X      [in/out] In input, <tt>X</tt> should contain a guess for the solution
	 *                of the linear systems (by default you should use 0).  In output, contains
	 *                the approximate solution to the linear systems as defined above.
	 * @param  convTester
	 *                [in] Determines the convergence criteria for the linear solver.
	 * 
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>this->opSupported(M_trans) (throw <tt>Exceptions::OpNotSupported</tt>)
	 * <li> <tt>Y.range()->isCompatible(M_trans==NOTRANS ? *this->range() : *this->domain()) == true</tt>
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>Y.domain()->isCompatible(*X->domain()) == true</tt>
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>X->range()->isCompatible(M_trans==NOTRANS ? *this->domain() : *this->range()) == true</tt>
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>X</tt> can not alias <tt>Y</tt>.  It is up to the client to ensure that <tt>X</tt>
	 *      and <tt>Y</tt> are distinct since in general this can not be verified by the implementation until,
	 *      perhaps, it is too late.  If possible, an exception will be thrown if aliasing is detected.
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>convTester==NULL</tt>] If the solver implementation has determined that suitable
	 *      solutions has been found, this method will return silently.
	 * <li> [<tt>convTester!=NULL</tt>] A solution will have been considered to have been found and this
	 *      method will return silently if <tt>convTester->convStatus(...,isConverged[])</tt> returned
	 *      <tt>isConverged[k] == true</tt> for every linear system being solved  from at least one call
	 *      the method <tt>convStatus(...)</tt>.
	 * <li> If <tt>this->nonsingStatus()==OP_SINGULAR</tt> and the right-hand-sides <tt>Y</tt> are not
	 *      consistent with <tt>*this</tt>, then a <tt>Solvers::Exceptions::SolverBreakdown</tt> exception will be thown.
	 * <li> If a solution could not be found bacause the linear operator was nearly singular,
	 *      then a <tt>Solvers::Exceptions::SolverBreakdown</tt> exception will be thrown.
	 * <li> If a solution could not be found for some other reason (i.e. the <tt>*convTester</tt> object
	 *      did not state that the system was solved), then a <tt>Solvers::Exceptions::FailureToConverge</tt> exception will be thrown.
	 * </ul>
	 *
	 * This method has a default implementation in terms of the
	 * <tt>solve()>/tt> method for vectors.
	 */
	virtual void solve(
		const ETransp                           M_trans
		,const MultiVector<Scalar>              &Y
		,MultiVector<Scalar>                    *X
		,const Scalar                           alpha       = 1.0
		,Solvers::ConvergenceTester<Scalar>     *convTester = NULL
		) const;

	//@}

};	// class LinearSolveOp

} // namespace Nonlin
} // namespace TSFCore

#endif	// TSFCORE_NONLIN_LINEAR_SOLVE_OP_DECL_HPP
