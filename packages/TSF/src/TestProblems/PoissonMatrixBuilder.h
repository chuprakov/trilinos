#ifndef POISSONMATRIXBUILDER_H
#define POISSONMATRIXBUILDER_H

#include "TSFConfig.h"
#include "TSFLinearOperator.h"
#include "TSFVector.h"
#include "MatrixOperatorBase.h"
#include "PComm.h"

namespace TSF
{
	

	/** \ingroup TestProblems
	 * Builds the matrix for a 1D poisson operator 
	 */
	class PoissonMatrixBuilder
		{
		public:
			
			PoissonMatrixBuilder(int nGlobal, const PComm& pComm);
			
			void buildSystem(TSFLinearOperator& op, 
											 TSFVector& rhs,
											 TSFVector& answer) const ;
			
			
		private:
			void buildMatrix(MatrixOperatorBase* A) const ;
			void buildRHS(TSFVector& rhs) const ;
			void buildAnswer(TSFVector& rhs) const ;
			
			double getRHS(int row) const ;
			double getSoln(int row) const ;
			int getRowSize(int row) const ;
			void getRowStruct(int row, TSFNonDupArray<Int>& indices) const ;
			void getRow(int row, TSFArray<int>& indices, TSFArray<TSFReal>& values) const ;
			
			int nGlobal_;
			TSFSmartPtr<TSFArray<int> > updates_;
			TSFArray<TSFNonDupArray<Int> > columns_;
			PComm comm_;
			double h_;
			int nProc_;
			int myPid_;
			int myNRows_;
		};
}

#endif
