#ifndef TSFADJOINTMATRIXOPERATOR_H
#define TSFADJOINTMATRIXOPERATOR_H

#include "TSFMatrixOperator.h"

namespace TSF
{

	/** \ingroup LinearOperatorSubtypes
	 * Wrapper for adjoint matrix operator.
	 */
	class TSFAdjointMatrixOperator : public TSFMatrixOperator {
	public:
		///
		TSFAdjointMatrixOperator( const TSFSmartPtr<TSFMatrixOperator>& mat_op );
		///
		void getILUKPreconditioner(
			int fillLevels, int overlapFill, TSFPreconditioner& rtn) const;
		///
		bool isFactored() const;
		///
		void factor();
		///
		void apply(const TSFVector& in, TSFVector& out) const;
		///
		void addToRow(int globalRowIndex,
					  int nCols,
					  const int* globalColumnIndices,
					  const TSFReal* a);
		///
		void zero();
				
	private:
		const TSFSmartPtr<TSFMatrixOperator>& mat_op_;
		TSFAdjointMatrixOperator(); // not defined and not to be called!

	};

}

#endif
