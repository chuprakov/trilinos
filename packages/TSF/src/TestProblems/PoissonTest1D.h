#ifndef POISSONTEST1D_H
#define POISSONTEST1D_H

#include "TSFConfig.h"
#include "TSFLinearOperator.h"
#include "TSFVector.h"
#include "MatrixOperatorBase.h"


namespace TSF
{
	using std::string;
	/** \ingroup TestProblems
	 * Builds the matrix for a 1D poisson operator 
	 */
	class MatrixTestProblem : public TSFTestProblemBase
		{
		public:
			/** */
 			MatrixTestProblem(int nGlobal, 
												MatrixOperatorBase* matrix,
												const PComm& pComm);
			
			/** */
			virtual ~MatrixTestProblem(){;}

			/** */
			virtual TSFLinearOperator getOperator();

			/** */
			virtual TSFVector getRHS();

			/** */
			virtual TSFVector getKnownSolution();
			
			
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
