// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraLinearOp.hpp

#ifndef TSFCORE_EPETRA_LINEAR_OP_HPP
#define TSFCORE_EPETRA_LINEAR_OP_HPP

#include "TSFCoreLinearOp.hpp"
#include "ReleaseResource.hpp"

class Epetra_Operator;

namespace TSFCore {

class EpetraVectorSpace;

///
/** Implementation of <tt>LinearOp</tt> using an <tt>Epetra_Operator</tt> object.
 *
 * This subclass can be used to represent the non-transposed operator
 * or transposed operator defined by an <tt>Epetra_Operator</tt>
 * object.  This class assumes that both the non-transposed and
 * transposed applications of the operator can be performed.
 */
class EpetraLinearOp : public LinearOp<RTOp_value_type> {
public:

	/** @name Public types */
	//@{

	///
	typedef RTOp_value_type Scalar;

	/// Defines the state of a <tt>EpetraLinearOp</tt> object.
	struct State {
		///
		State() : opTrans(NOTRANS) {}
		///
		MemMngPack::ref_count_ptr<Epetra_Operator>                  op;
		///
		ETransp                                                     opTrans;
		///
		MemMngPack::ref_count_ptr<MemMngPack::ReleaseResource>      extra;
	};

	//@}

	/** @name Constructors / initializers / accessors */
	//@{

	///
	/** Construct to uninitialized.
	 *
	 * See the postconditions for <tt>setUninitialized()</tt>
	 */
	EpetraLinearOp();

	/// Calls <tt>initialize()</tt>.
	EpetraLinearOp(
		const MemMngPack::ref_count_ptr<Epetra_Operator>                &op
		,ETransp                                                        opTrans = NOTRANS
		,const MemMngPack::ref_count_ptr<MemMngPack::ReleaseResource>   &extra  = MemMngPack::null
		);

	///
	/** Initialize
	 *
	 * @param  op       [in] The <tt>Epetra_Operator</tt> this <tt>*this</tt> will wrap.
	 * @param  opTrans  [in] If <tt>opTrans==NOTRANS</tt> then <tt>op</tt> will be viewed as <tt>op</tt>
	 *                  and if <tt>opTrans==TRANS</tt> then <tt>op</tt> will be viewed as its transpose
	 *                  <tt>op'</tt> for the behavior of <tt>apply()</tt>.
	 * @param  extra    [in] Any extra data that is associated with the use of <tt>op</tt> but who's
	 *                  memory is not managed in <tt>op</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>op.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->state().op.get() == op.get()</tt>
	 * <li> <tt>this->state().opTrans == opTrans</tt>
	 * <li> <tt>this->state().extra.get() == extra.get()</tt>
	 * <li> <tt>this->domain().get() != NULL</tt>
	 * <li> <tt>this->range().get() != NULL</tt>
	 * </ul>
	 */
	void initialize(
		const MemMngPack::ref_count_ptr<Epetra_Operator>                &op
		,ETransp                                                        opTrans = NOTRANS
		,const MemMngPack::ref_count_ptr<MemMngPack::ReleaseResource>   &extra  = MemMngPack::null
		);
	
	///
	/** Set to uninitialized and return the current state.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->state().op.get() == NULL</tt>
	 * <li> <tt>this->state().opTrans == NOTRANS</tt>
	 * <li> <tt>this->state().extra.get() == NULL</tt>
	 * <li> <tt>this->domain().get() == NULL</tt>
	 * <li> <tt>this->range().get() == NULL</tt>
	 * </ul>
	 */
	State setUninitialized();
	
	///
	/** Return the current state.
	 *
	 * <b>Warning!</b> Do not attempt to change the state of any of
	 * the objects exposed by <tt>return</tt>.  As long as a client
	 * does not use any dynamic, const or static casting, then the
	 * compiler will not allow this so in a sense this method is safe.
	 */
	const State& state() const;

	//@}
	
	/** @name Overridden from OpBase */
	//@{
	
	///
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> > range() const;
	///
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> > domain() const;
	
	//@}
	
	/** @name Overridden from LinearOp */
	//@{
	
	///
	MemMngPack::ref_count_ptr<const LinearOp<Scalar> > clone() const;
	///
	void apply(
		const ETransp            M_trans
		,const Vector<Scalar>    &x
		,Vector<Scalar>          *y
		,const Scalar            alpha
		,const Scalar            beta
		) const;

	//@}

private:

	// ////////////////////////////////////
	// Private data members

	State                                               state_;
	MemMngPack::ref_count_ptr<const EpetraVectorSpace>  domain_;
	MemMngPack::ref_count_ptr<const EpetraVectorSpace>  range_;

};	// end class EpetraLinearOp

// ////////////////////////////////////////
// Inline members

inline
const EpetraLinearOp::State&
EpetraLinearOp::state() const
{
	return state_;
}

}	// end namespace TSFCore

#endif	// TSFCORE_EPETRA_LINEAR_OP_HPP
