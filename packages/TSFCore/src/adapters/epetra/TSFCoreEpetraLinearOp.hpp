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
		const Teuchos::RefCountPtr<Epetra_Operator>   &op
		,ETransp                                      opTrans = NOTRANS
		);

	///
	/** Initialize
	 *
	 * @param  op       [in] The <tt>Epetra_Operator</tt> this <tt>*this</tt> will wrap.
	 * @param  opTrans  [in] If <tt>opTrans==NOTRANS</tt> then <tt>op</tt> will be viewed as <tt>op</tt>
	 *                  and if <tt>opTrans==TRANS</tt> then <tt>op</tt> will be viewed as its transpose
	 *                  <tt>op'</tt> for the behavior of <tt>apply()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>op.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->domain().get() != NULL</tt>
	 * <li> <tt>this->range().get() != NULL</tt>
	 * </ul>
	 */
	void initialize(
		const Teuchos::RefCountPtr<Epetra_Operator>   &op
		,ETransp                                      opTrans = NOTRANS
		);
	
	///
	/** Set to uninitialized and return the current state.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->domain().get() == NULL</tt>
	 * <li> <tt>this->range().get() == NULL</tt>
	 * </ul>
	 */
	void setUninitialized(
		Teuchos::RefCountPtr<Epetra_Operator>    *op      = NULL
		,ETransp                                 *opTrans = NULL
		);
	
	//@}
	
	/** @name Overridden from OpBase */
	//@{
	
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> > range() const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> > domain() const;
	
	//@}
	
	/** @name Overridden from LinearOp */
	//@{
	
	///
	Teuchos::RefCountPtr<const LinearOp<Scalar> > clone() const;
	///
	void apply(
		const ETransp            M_trans
		,const Vector<Scalar>    &x
		,Vector<Scalar>          *y
		,const Scalar            alpha
		,const Scalar            beta
		) const;
	///
	void apply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha
		,const Scalar                 beta
		) const;

	//@}

private:

	// ////////////////////////////////////
	// Private data members

	Teuchos::RefCountPtr<Epetra_Operator>          op_;
	ETransp                                        opTrans_;
	Teuchos::RefCountPtr<const EpetraVectorSpace>  domain_;
	Teuchos::RefCountPtr<const EpetraVectorSpace>  range_;

};	// end class EpetraLinearOp

}	// end namespace TSFCore

#endif	// TSFCORE_EPETRA_LINEAR_OP_HPP
