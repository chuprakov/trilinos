/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFTRANSPOSEOPERATOR_HPP
#define TSFTRANSPOSEOPERATOR_HPP

#include "TSFConfigDefs.hpp"

namespace TSFExtended
{
  /** 
   * Class TransposeOperator represents the unformed transpose \f$A^T\f$
   * of another operator \f$A\f$. The transpose itself is never formed,
   * rather, any operator applications are carried out with the 
   * ETransp flag set to TRANS in the apply() method. 
   *
   * Should one ever want to form an explicit transpose, the 
   * form() method should be called.
   * 
   */
  template <class Scalar>
  class TransposeOperator : public TSFCore::LinearOp<Scalar>,
                            public Formable<Scalar>,
                            public ExplicitlyTransposeableOp<Scalar>
  {
  public:
    /** Create with an operator. */
    TransposeOperator(const LinearOperator& op) : op_(op) {;}

    /** Virtual dtor */
    virtual ~TransposeOperator(){;}

    /** Return the domain, which is the range of the operator being
     * transposed. */
    RefCountPtr<const TSFCore::VectorSpace<Scalar> > domain() const 
    {return op_.range();}

    /** Return the range, which is the domain of the operator being
     * transposed. */
    RefCountPtr<const TSFCore::VectorSpace<Scalar> > range() const 
    {return op_.domain();}

    /** Apply the transpose of the underlying operator */
    void apply(
               const ETransp            M_trans
               ,const Vector<Scalar>    &x
               ,Vector<Scalar>          *y
               ,const Scalar            alpha
               ,const Scalar            beta
               ) const 
    {
      op_.ptr()->apply(not_trans(M_trans), x, y, alpha, beta);
    }

    /** Form an explicit representation if supported by the underlying
     * operator. */
    LinearOperator<Scalar> form() const 
    {
      const ExplicitlyTransposeable* et = 
        dynamic_cast<const ExplicitlyTransposeable*>(op_.ptr().get());


      TEST_FOR_EXCEPTION(et==0, runtime_error,
                         "TransposeOperator<Scalar>::form() called where "
                         "the operator is unable to form an explicit "
                         "transpose. The operator is " << op_.describe());

      return et->formTranspose();
    }

    /** Form the transpose of this operator, which is a copy of the
     * underlying operator. */
    LinearOperator<Scalar> formTranspose() const
    {
      return op_.clone();
    }
    
  private:
    
    LinearOperator<Scalar> op_;
  };
}
