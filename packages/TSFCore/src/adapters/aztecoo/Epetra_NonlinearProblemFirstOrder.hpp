// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

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
	/** Returns true if the operator <tt>DcDy_op</tt> computed by <tt>calc_Dc(...)</tt> is a constant object.
	 *
	 * The default implementation returns false.
	 */
	virtual bool DcDy_op_is_const() const;
  
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

	///
	/** Returns true if the preconditioner operator <tt>DcDy_prec</tt> computed by <tt>calc_Dc(...)</tt> is a constant object.
	 *
	 * The default implementation returns false.
	 */
	virtual bool DcDy_prec_is_const() const;

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
  virtual bool use_DcDu_op(int l) const;

	///
  /** Returns an (uninitialized) <tt>Epetra_Operator</tt> that will store
   * <tt>D(c)/D(u(l))</tt>.
   *
   * <b>Postconditions:</b>
   * <ul>
   * <li> [<tt>this->use_DcDu_op(l)==true</tt>] <tt>this->create_DcDu(l).get()!=NULL</tt>
   * <li> [<tt>this->use_DcDu_op(l)==false</tt>] <tt>this->create_DcDu(l).get()==NULL</tt>
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
   * <li> [<tt>this->use_DcDu_op(l)==false</tt>] <tt>this->create_DcDu_mv(l).get()!=NULL</tt>
   * <li> [<tt>this->use_DcDu_op(l)==true</tt>] <tt>this->create_DcDu_mv(l).get()==NULL</tt>
   * </ul>
   *
   * The default implementation returns
   *\verbatim
    if(this->use_DcDu_op(l)) {
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

	///
	/** Returns true if the object (Epetra_MultiVector or Epetra_Operator) for
	 * <tt>DcDu(l)</tt> computed by <tt>calc_Dc(...)</tt> is a constant object.
	 *
	 * The default implementation returns false.
	 */
	virtual bool DcDu_is_const(int l) const;

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

	/** @name Overridden from NonlinearProblem */
	//@{
	
	/// Calls <tt>calc_Dc(y,u,c,NULL,NULL,NULL)
	void calc_c(
		const Epetra_Vector     &y
		,const Epetra_Vector*   u[]
    ,Epetra_Vector          *c
		) const;
	/// Calls <tt>calc_Dg(y,u,g,NULL,NULL,NULL)
	void calc_g(
		const Epetra_Vector     &y
		,const Epetra_Vector*   u[]
    ,Epetra_Vector          *g
		) const;

	//@}


}; // class NonlinearProblemFirstOrder

} // namespace Epetra

#endif // EPETRA_NONLINEAR_PROBLEM_FIRST_ORDER_HPP

