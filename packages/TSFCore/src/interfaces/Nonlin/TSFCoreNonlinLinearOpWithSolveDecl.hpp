// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreNonlinLinearOpWithSolveDecl.hpp

#ifndef TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_DECL_HPP

#include "TSFCoreNonlinLinearSolveOp.hpp"
#include "TSFCoreLinearOp.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Base class for all linear operators that support <tt>apply()</tt> and <tt>solve()</tt>.
 */
template<class Scalar>
class LinearOpWithSolve : virtual public LinearOp<Scalar>, virtual public LinearSolveOp<Scalar> {
public:

	/** @name Virtual functions with default implemenations */
	//@{

	///
	/** Clone the nonsingular linear operator object (if supported).
	 *
	 * The primary purpose for this method is to allow a client to
	 * capture the current state of an nonsingular linear operator
	 * object and be guaranteed that some other client will not alter
	 * its behavior.  A smart implementation will use reference
	 * counting and lazy evaluation internally and will not actually
	 * copy any large amount of data unless it has to.
	 *
	 * The default implementation returns <tt>return.get()==NULL</tt>
	 * which is allowable by this specification.  An nonsingular linear
	 * operator object is not required to return a non-NULL value but
	 * almost every good implementation should and will.
	 */
	virtual MemMngPack::ref_count_ptr<const LinearOpWithSolve> clone_lows() const;

	///
	/** Get access to the preconditioner to this nonsingular linear operator.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * </ul>
	 *
	 * The default implementation returns <tt>return.get()==NULL</tt>.
	 */
	MemMngPack::ref_count_ptr<const LinearOp<Scalar> > preconditioner() const;

	//@}

	/** @name Overridden methods from LinearOp */
	//@{

	/// This method is simply overridden to return <tt>this->clone_lons()</tt>.
	MemMngPack::ref_count_ptr<const LinearOp<Scalar> > clone() const;

	//@}

	/** @name Overridden methods from LinearSolveOp */
	//@{

	/// This method is simply overridden to return <tt>this->clone_lons()</tt>.
	MemMngPack::ref_count_ptr<const LinearSolveOp<Scalar> > clone_lso() const;

	//@}

private:

#ifdef DOXYGEN_COMPILE
	LinearOp<Scalar>*  preconditioner; // doxygen only!
#endif

};	// class LinearOpWithSolve

} // namespace Nonlin
} // namespace TSFCore

#endif	// TSFCORE_NONLIN_LINEAR_OP_WITH_SOLVE_DECL_HPP
