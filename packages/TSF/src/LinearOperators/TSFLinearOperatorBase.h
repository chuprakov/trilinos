#ifndef TSFLINEAROPERATORBASE_H
#define TSFLINEAROPERATORBASE_H

#include "TSFConfig.h"
#include "TSFVectorSpace.h"

namespace TSF
{
	class TSFLinearSolver;
	class TSFLinearOperator;
	

	/** \ingroup LinearOperatorSubtypes
	 * Base class for all linear operators. The apply() method is pure virtual, 
	 * and so must be overridden by any derived class. The applyInverse() and
	 * applyAdjoint() methods have default implementations that throw exceptions.
	 */

	class TSFLinearOperatorBase
		{
		public:
			/** Empty ctor, sets domain and range to empty spaces */
			TSFLinearOperatorBase();
			/** Ctor sets domain and range spaces to specified values */
			TSFLinearOperatorBase(const TSFVectorSpace& domain,
														const TSFVectorSpace& range);
			/* the usual virtual dtor */
			virtual ~TSFLinearOperatorBase(){;}


			/** get the number of block rows */
			virtual int numBlockRows() const {return 1;}

			/** get the number of block columns */
			virtual int numBlockCols() const {return 1;}

			/** get the (i,j)-th submatrix */
			virtual void getBlock(int i, int j, TSFLinearOperator& sub) const ;

			/** set the (i,j)-th submatrix */
			virtual void setBlock(int i, int j, 
														const TSFLinearOperator& sub);

			/** say whether we are a block operator */
			virtual bool isBlockOperator() const {return false;}

			/** return domain space */
			const TSFVectorSpace& domain() const {return domain_;}
			/** return range space */
			const TSFVectorSpace& range() const {return range_;}

			/** say whether we are a zero operator */
			virtual	bool isZeroOperator() const {return false;}

			/** apply operator to a vector in the domain space, returning a vector
			 * in the range space */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const = 0 ;

			/** apply inverse operator to a vector in the range space, returning
			 * its preimage as a vector in the domain space. The default implementation
			 * throws an exception */
			virtual void applyInverse(const TSFVector& in,
																TSFVector& out) const ;	

			/** apply inverse operator using the specified solver */
			virtual void applyInverse(const TSFLinearSolver& solver,
																const TSFVector& in,
																TSFVector& out) const ;

			/** apply adjoint operator to a vector in the domain space, returning
			 * a vector in the range space. The default implementation throws an
			 * exception */
			virtual void applyAdjoint(const TSFVector& in,
																TSFVector& out) const ;

			/** apply inverse adjoint operator */
			virtual void applyInverseAdjoint(const TSFVector& in,
																			 TSFVector& out) const ;

			/** apply inverse adjoint operator */
			virtual void applyInverseAdjoint(const TSFLinearSolver& solver,
																			 const TSFVector& in,
																			 TSFVector& out) const ;

			/** */
			virtual void getInverse(const TSFLinearSolver& solver,
															const TSFLinearOperator& self,
															TSFLinearOperator& inv) const ;

			/** */
			virtual void getInverse(const TSFLinearOperator& self,
															TSFLinearOperator& inv) const ;

			/** */
			virtual void getAdjoint(const TSFLinearOperator& self,
															TSFLinearOperator& adj) const ;

			/** */
			virtual void getInverseAdjoint(const TSFLinearOperator& self,
															TSFLinearOperator& invAdj) const ;

			/** */
			virtual void getInverseAdjoint(const TSFLinearSolver& solver,
																		 const TSFLinearOperator& self,
																		 TSFLinearOperator& invAdj) const ;

			/** */
			virtual bool isMatrixOperator() const {return false;}

			/**
			 * Write to a stream 
			 */
			virtual void print(ostream& os) const ;

		protected:
			TSFVectorSpace domain_;
			TSFVectorSpace range_;
		};
}

#endif
