// /////////////////////////////////////////////////////////////////
// Epetra_NonlinearProblemFirstOrder.hpp

#ifndef EPETRA_NONLINEAR_PROBLEM_FIRST_ORDER_HPP
#define EPETRA_NONLINEAR_PROBLEM_FIRST_ORDER_HPP

#include "Epetra_NonlinearProblem.hpp"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"

namespace Epetra {

///
/** Simple class for passing and Epetra_Operator or Epetra_MultiVector object.
 *
 * <b>Assertions:</b>
 * <ul>
 * <li> <tt>( op==NULL && mv==NULL )|| ( op!=NULL != mv!=NULL )</tt>
 * </ul>
 */
class EpetraOp_or_EpetraMV {
public:
  EpetraOp_or_EpetraMV() : op_(NULL), mv_(NULL) {}
  EpetraOp_or_EpetraMV( Epetra_Operator *op ) : op_(op), mv_(NULL) {}
  EpetraOp_or_EpetraMV( Epetra_MultiVector  *mv ) : op_(NULL), mv_(mv) {}
  Epetra_Operator*     op() const { return op_; }
  Epetra_MultiVector*  mv() const { return mv_; }
private:
  Epetra_Operator      *op_;
  Epetra_MultiVector   *mv_;
};

///
/** Enumeration for determining how an operator is applied.
 */
enum ETransp {
	NOTRANS  ///< Use the transposed operator
	,TRANS   ///< Use the nontransposed operator
};

///
/** Epetra namespace for a nonlinear problem.
 *
 * ToDo: Finish documentation!
 */
class NonlinearProblemFirstOrder : virtual public NonlinearProblem {
public:

	/** @name Adjoints supported? */
	//@{

	///
	virtual bool adjointSupported() const;

	//@}

	/** @name Factories for linear operators */
	//@{

	///
	virtual Teuchos::RefCountPtr<Epetra_Operator> create_DcDy_op() const = 0;
  
  ///
  /** Determines if <tt>this</tt> can create and define a specialized preconditioner for DcDy.
   *
   * The default implementation returns <tt>false</tt>.
   */
  virtual bool specialized_DcDy_prec() const;

	///
  /** Create a (possibly uninitialized) specialized preconditioner for DcDy.
   *
   * <b>Postconditions:</b>
   * <ul>
   * <li> [<tt>this->specialized_DcDy_prec()==true</tt>] <tt>this->create_DcDy_prec().get()!=NULL</tt>
   * <li> [<tt>this->specialized_DcDy_prec()==false</tt>] <tt>this->create_DcDy_prec().get()==NULL</tt>
   * </ul>
   *
   * The default implementation returns <tt>return.get()==NULL</tt>.
   */
	virtual Teuchos::RefCountPtr<Epetra_Operator> create_DcDy_prec() const;

  /// Return if a <tt>Epetra_Operator</tt> is used for <tt>DcDu(l)</tt>.
  /**
   * @return Returns <tt>true</tt> if <tt>this->create_DcDu(l).get()!=NULL</tt>.
   * otherwise returns <tt>false</tt>.
   *
   * If this function returns returns <tt>false</tt> then
   * <tt>Epetra_MultiVector(*this->map_c(),this->map_u(l)->NumGlobalElements())</tt>
   * should be used to create the object for <tt>DcDu(l)</tt> that is passed
   * to <tt>calc_Dc(..)</tt>.
   *
   * The default implementation returns <tt>false</tt>.
   */
  virtual bool use_EO_DcDu(int l) const;

	///
  /** Returns an (uninitialized) <tt>Epetra_Operator</tt> that will store
   * <tt>D(c)/D(u(l))</tt>.
   *
   * <b>Postconditions:</b>
   * <ul>
   * <li> [<tt>this->use_EO_DcDu(l)==true</tt>] <tt>this->create_DcDu(l).get()!=NULL</tt>
   * <li> [<tt>this->use_EO_DcDu(l)==false</tt>] <tt>this->create_DcDu(l).get()==NULL</tt>
   * </ul>
   *
   * The default implementation returns <tt>return.get()==NULL</tt>.
   */
 	virtual Teuchos::RefCountPtr<Epetra_Operator> create_DcDu_op(int l) const;

	///
  /** Returns a (potentially uninitialized) <tt>Epetra_MultiVector</tt> that will store
   * <tt>D(c)/D(u(l))</tt>.
   *
   * <b>Postconditions:</b>
   * <ul>
   * <li> [<tt>this->use_EO_DcDu(l)==false</tt>] <tt>this->create_DcDu_mv(l).get()!=NULL</tt>
   * <li> [<tt>this->use_EO_DcDu(l)==true</tt>] <tt>this->create_DcDu_mv(l).get()==NULL</tt>
   * </ul>
   *
   * The default implementation returns
   *\verbatim
    if(this->use_EO_DcDu(l)) {
      return Teuchos::null;
    }
    else {
      return Teuchos::rcp(
        new Epetra_MultiVector(
          *this->map_c()
          ,map_u(l)->NumGlobalElements()
          )
        );
    }
   *\endverbatim

   */
 	virtual Teuchos::RefCountPtr<Epetra_MultiVector> create_DcDu_mv(int l) const;

	//@}

	/** @name Transpose arguments */
	//@{

	///
	virtual ETransp opDcDy() const = 0;

	///
	virtual ETransp opDcDu(int l) const;

	//@}

	/** @name Calculation methods */
	//@{

	///
	virtual void calc_Dc(
		const Epetra_Vector           &y
		,const Epetra_Vector*         u[]
    ,Epetra_Vector                *c
    ,Epetra_Operator              *DcDy_op
    ,Epetra_Operator              *DcDy_prec
    ,const EpetraOp_or_EpetraMV   DcDu[]
    ) const = 0;

	///
	virtual void calc_Dg(
		const Epetra_Vector           &y
		,const Epetra_Vector*         u[]
    ,Epetra_Vector                *g
    ,Epetra_MultiVector           *DgDy
    ,Epetra_MultiVector*          DgDu[]
		) const;

	//@}

}; // class NonlinearProblemFirstOrder

} // namespace Epetra

#endif // EPETRA_NONLINEAR_PROBLEM_FIRST_ORDER_HPP

