#ifndef TSFTESTPROBLEM1D_H
#define TSFTESTPROBLEM1D_H

#include "TSFConfig.h"
#include "TSFDefaultMatrixProblem.h"

namespace TSF
{
	using std::ostream;

	/** \ingroup CoreSubtypes
	 * TSFMatrixProblem
	 * 
	 */

	class TSFTestProblem1D : public TSFDefaultMatrixProblem
		{
		public:
			/** empty ctor only */
			TSFTestProblem1D(TSFReal a, TSFReal b, int n, 
												 TSFMatrixOperator* matrix);

			/** virtual dtor */
			virtual ~TSFTestProblem1D(){;}

			/** */
			virtual TSFReal secondOrderTerm(const TSFReal& x) const {return 1.0;}

			/** */
			virtual TSFReal firstOrderTerm(const TSFReal& x) const {return 0.0;}

			/** */
			virtual TSFReal zerothOrderTerm(const TSFReal& x) const {return 0.0;}

			/** */
			virtual TSFReal forcingTerm(const TSFReal& x) const {return 1.0;}

			/** */
			virtual TSFReal leftBCForcingTerm() const {return 0.0;}

			/** */
			virtual TSFReal rightBCForcingTerm() const {return 0.0;}

			/** */
			virtual TSFReal leftBCZerothOrderTerm() const {return 0.0;}
			
			/** */
			virtual TSFReal rightBCZerothOrderTerm() const {return 0.0;}


			/** */
			virtual TSFReal leftBCFirstOrderTerm() const {return 0.0;}
			
			/** */
			virtual TSFReal rightBCFirstOrderTerm() const {return 0.0;}

			/** */
			virtual TSFReal solution(const TSFReal& x) const {return 0.0;}
		protected:

			/* */
			virtual void getRowValues(int row, TSFArray<int>& indices, 
																TSFArray<TSFReal>& values) const;

			/* */
			virtual void getRowStruct(int row, 
																TSFNonDupTSFArray<int>& indices) const;

			/* */
			virtual TSFReal getRHSValue(int row) const ;


			/* */
			virtual TSFReal getSolutionValue(int row) const ;

			/** */
			virtual int nGlobalRows() const {return nGlobal_;}

			/** */
			virtual int nLocalRows() const {return nLocal_;}

			/** */
			virtual int lowestLocalRow() const {return myLowestRow_;}

			/** */
			double getX(int i) const ;

			TSFReal h_;

			TSFReal a_;

			TSFReal b_;
			
			int nGlobal_;

			int nLocal_;

			int myLowestRow_;
		};
}

#endif
