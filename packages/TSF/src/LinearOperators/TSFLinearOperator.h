#ifndef TSFLINEAROPERATOR_H
#define TSFLINEAROPERATOR_H

#include "TSFVector.h"
#include "TSFSmartPtr.h"
#include "TSFTimeMonitor.h"
#include "TSFMatrixOperator.h"

namespace TSF
{

  class TSFLinearOperatorBase;
  class TSFMatrixOperator;
  class TSFLinearSolver;
  
  /** \ingroup LinearOperator 
   * LinearOperator is the user-level object representing linear mappings
   * from a domain TSFVectorSpace to a range TSFVectorSpace.
   * 
   */
  
  class TSFLinearOperator
    {
    public:
      /**\name constructors */
      //@{
      /** empty ctor constructs a null linear operator. This is primarily for
       * use with templated container classes. */
      TSFLinearOperator();
      /** create a TSFLinearOperator from a pointer to a subtype. */
      TSFLinearOperator(TSFLinearOperatorBase* ptr);
      //@}
      
      /** \name introspection methods */
      //@{
      /** return domain space */
      const TSFVectorSpace& domain() const ;
      /** return range space */
      const TSFVectorSpace& range() const ;
      
      /** am I a BlockOperator?  */
      bool isBlockOperator() const;
      
      /** am I a zero operator? */
      bool isZeroOperator() const ;
      //@}
      
      /** \name Block structure information */
      //@{
      /** get the number of block rows */
      int numBlockRows() const ;
      
      /** get the number of block columns */
      int numBlockCols() const ;
      
      /** get the (i,j)-th submatrix of a block operator */
      TSFLinearOperator getBlock(int i, int j) const ;
      
      /** set the (i,j)-th submatrix of a block operator */
      TSFLinearOperator& setBlock(int i, int j, 
				  const TSFLinearOperator& sub);
      //@}
      
      
      /** \name Application of forward, inverse, adjoint, 
       * and adjoint image operators */
      //@{
      /** apply the operator to a vector. The input vector is tested to ensure 
       * that it is in the domain of the operator */
      void apply(const TSFVector& arg, TSFVector& out) const ;
      
      /** apply the operator's inverse to a vector */
      void applyInverse(const TSFVector& arg, TSFVector& out) const ;
      /** apply the operator's inverse to a vector, using the specified 
       * solver. */
      void applyInverse(const TSFLinearSolver& solver, 
			const TSFVector& arg, TSFVector& out) const ;
      
      /** apply the operator's adjoint to a vector */
      void applyAdjoint(const TSFVector& arg, TSFVector& out) const ;
      
      /** apply the operator's inverse to a vector */
      void applyInverseAdjoint(const TSFVector& arg, TSFVector& out) const ;
      /** apply the operator's inverse to a vector, using the specified 
       * solver. */
      void applyInverseAdjoint(const TSFLinearSolver& solver, 
			       const TSFVector& arg, TSFVector& out) const ;
      
      /** overloaded matrix-vector multiply.  */
      TSFVector operator*(const TSFVector& x) const ;
      //@}
      
      /** \name Creation of auxiliary operators */
      //@{
      /** get an operator object representing the inverse. */
      TSFLinearOperator inverse() const ;
      
      /** get an operator object representing the inverse. The inverse
       * operator will use the specified solver to solve systems */
      TSFLinearOperator inverse(const TSFLinearSolver& solver) const ;
      
      /** get an operator object representing the adjoint */
      TSFLinearOperator adjoint() const ;
      
      /** get an operator object representing the inverse adjoint */
      TSFLinearOperator inverseAdjoint() const ;
      
      /** get an operator object representing the inverse adjoint */
      TSFLinearOperator inverseAdjoint(const TSFLinearSolver& solver) const ;
      //@}
      
      /** \name */
      //@{
      /** scale by a constant (deferred) */
      TSFLinearOperator operator*(const double& scale) const ;
      
      /** addition of two operators (deferred) */
      TSFLinearOperator operator+(const TSFLinearOperator& op) const ;
      
      /** subtraction of two operators (deferred) */
      TSFLinearOperator operator-(const TSFLinearOperator& op) const ;
      
      /** multiplication (composition) of two operators (deferred) */
      TSFLinearOperator operator*(const TSFLinearOperator& op) const ;
      
      /** negation of an operator (deferred) */
      TSFLinearOperator operator-() const ;
      
      //@}
      
      /** \name low-level access */
      //@{
      /** access to the pointer to the linear operator implementation. 
       *  For developer use only. */
      TSFSmartPtr<TSFLinearOperatorBase>& getPtr() {return ptr_;}
      
      /** access to the pointer to the linear operator implementation. 
       *  For developer use only. */
      const TSFSmartPtr<TSFLinearOperatorBase>& getPtr() const {return ptr_;}
      
      /** determine whether this operator is a matrix */
      bool isMatrixOperator() const ;
      
      /** access to the pointer to the matrix implementation. If the operator
       * is not a TSFMatrixOperator, throw an error */
      const TSFSmartPtr<const TSFMatrixOperator> getMatrix() const ;
      //@}
      
      /** \name output */
      //@{
      /** write to a stream */
      friend ostream& operator<<(ostream& os, const TSFLinearOperator& op);		
      
      /** write to a string */
      string toString() const ;
      //@
    private:
      
      /* pointer to a concrete type */
      TSFSmartPtr<TSFLinearOperatorBase> ptr_;
      static TSFTimer opTimer_;
    };
  
  /** \relates TSFLinearOperator left multiplication by a scalar */
  TSFLinearOperator operator*(const TSFReal& scale, 
			      const TSFLinearOperator& op);
  
}

#endif
