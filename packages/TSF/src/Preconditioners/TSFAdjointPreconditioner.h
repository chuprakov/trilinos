#ifndef TSFADJOINTPRECONDITIONER_H
#define TSFADJOINTPRECONDITIONER_H

#include "TSFConfig.h"
#include "TSFPreconditionerBase.h"
#include "TSFPreconditioner.h"
#include "TSFLinearOperator.h"

namespace TSF
{
	
	
	/** \ingroup Core
	 * Wrapper for adjoint preconditioners.
	 */
	class TSFAdjointPreconditioner : public TSFPreconditionerBase {
	public:
		///
		TSFAdjointPreconditioner(const TSFPreconditioner& prec);
		///
		TSFLinearOperator left();
		///
		TSFLinearOperator right();
		///
		bool isTwoSided() const;
		///
		bool hasLeft() const;
		///
		bool hasRight() const;
		///
		bool isIdentity() const;
		///
		string toString() const;
		
		
	private:
		TSFPreconditioner prec_;
		TSFAdjointPreconditioner(); // not defined and not to be called!
	};
	
}

#endif // TSFPRECONDITIONER_H
