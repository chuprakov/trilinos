#ifndef BVP1D_H
#define BVP1D_H

#include "TSFConfig.h"
#include "TSFArray.h"
#include "TSFDefaultMatrixProblem.h"

namespace TSF
{
	
	
	/** \ingroup TestProblems
	 * Sets up a boundary value problem of the form
	 *
	 * D[alpha(x) D[u(x)]] + beta(x) u(x) = gamma(x)
	 *
	 * with dirichlet boundary conditions
	 *
	 * u(a) = uLeft
	 * u(b) = uRight
	 */
	class BVP1D : public TSFDefaultMatrixProblem 
		{
		public:
			BVP1D(int n, TSFReal a, TSFReal b, TSFMatrixOperator* matrix);
			
			virtual ~BVP1D(){;}

			virtual TSFReal a() const {return a_;}

			virtual TSFReal b() const {return b_;}

			virtual int n() const {return n_;}

			virtual TSFReal h() const {return h_;}

			virtual TSFReal uLeft() const {return 0.0;}

			virtual TSFReal uRight() const {return 0.0;}

			virtual TSFReal alpha(TSFReal x) const {return 1.0;}

			virtual TSFReal beta(TSFReal x) const {return 0.0;}

			virtual TSFReal gamma(TSFReal x) const {return 1.0;}

			virtual TSFReal solution(TSFReal x) const = 0 ;
			
		private:
			virtual int nGlobalRows() const {return n_;}

			virtual int nLocalRows() const {return nLocal_;}

			virtual int lowestLocalRow() const {return lowestLocalRow_;}

			TSFReal getX(int i) const {return a_ + ((TSFReal) i)*h_;}

			virtual void getRowValues(int row, TSFArray<int>& indices, 
																TSFArray<TSFReal>& values) const ;

			virtual int getRowBandwidth(int row) const ;

			/* */
			virtual TSFReal getRHSValue(int row) const ;

			/* */
			virtual TSFReal getSolutionValue(int row) const ;

			/** create the update list */
			virtual TSFSmartPtr<TSFArray<int> > formUpdateList() const ;

			int n_;

			double a_;

			double b_;

			TSFReal h_;

			int nLocal_;

			int lowestLocalRow_;
		};
}

#endif
