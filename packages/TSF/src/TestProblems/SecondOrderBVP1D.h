#ifndef SECONDORDERBVP1D_H
#define SECONDORDERBVP1D_H

#include "TSFConfig.h"
#include "TSFDefaultMatrixProblem.h"

namespace TSF
{
	using std::ostream;

	/** \ingroup CoreSubtypes
	 * TSFMatrixProblem
	 * 
	 */

	class SecondOrderBVP1D : public TSFOneDTestProblem 
		{
		public:
			/** empty ctor only */
			TSFOneDTestProblem(TSFReal a, TSFReal b, int n, 
												 TSFMatrixOperator* matrix);

			/** virtual dtor */
			virtual ~TSFOneDTestProblem(){;}

			/** */
			virtual TSFReal secondOrderTerm(const TSFReal& x) const {return 1.0;}

			/** */
			virtual TSFReal firstOrderTerm(const TSFReal& x) const {return 0.0;}

			/** */
			virtual TSFReal zerothOrderTerm(const TSFReal& x) const {return 0.0;}

			/** */
			virtual TSFReal forcingTerm(const TSFReal& x) const {return 1.0;}
