// /////////////////////////////////////////////////////////////////
// Epetra_NonlinearProblemFirstOrder.hpp

#ifndef EPETRA_NONLINEAR_PROBLEM_FIRST_ORDER_HPP
#define EPETRA_NONLINEAR_PROBLEM_FIRST_ORDER_HPP

#include "Epetra_NonlinearProblem.hpp"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"

namespace Epetra {

///
/** Struct for passing Epetra_Operator or Epetra_MultiVector
 *
 * <b>Assertions:</b>
 * <ul>
 * <li> <tt>!( op!=NULLL && mv!=NULL )</tt>
 * </ul>
 */
struct EpetraOp_or_EpetraMV {
  Epetra_Operator      *op;
  Epetra_MultiVector   *mv;
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
	virtual Teuchos::RefCountPtr<Epetra_Operator> create_DcDy() const = 0;

	///
  /** Returns an (uninitialized) <tt>Epetra_Operator</tt> that will store
   * <tt>D(c)/D(u(l))</tt>.
   *
   * If this function return <tt>return.get()==NULL</tt> then
   * a <tt>Epetra_MultiVector</tt> must be used instead.
   */
 	virtual Teuchos::RefCountPtr<Epetra_Operator> create_DcDu(int l) const;

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
    ,Epetra_Operator              *DcDy
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

