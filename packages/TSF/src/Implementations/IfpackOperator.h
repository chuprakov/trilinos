#ifndef IFPACKOPERATOR_H
#define IFPACKOPERATOR_H

#include "TSFConfig.h"
#include "TSFVector.h"
#include "TSFLinearOperatorBase.h"
#include "TSFTimeMonitor.h"

#if HAVE_PETRA


#define PETRA_BOOL_SUPPORTED

#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"
#include "PetraVector.h"
#include "PetraMatrix.h"


namespace TSF
{
	using std::ostream;


	/** \ingroup Petra
	 * 
	 */
	
	class IfpackOperator : public TSFLinearOperatorBase 
		{
		public:
			/** empty ctor */
			IfpackOperator(const TSFVectorSpace& domain,
				       const TSFVectorSpace& range,
				       Ifpack_CrsRiluk* prec,
				       Ifpack_IlukGraph* graph);
			/** TUVD */
			virtual ~IfpackOperator();

			/** apply operator to a vector in the domain space, returning a vector
			 * in the range space */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;

				/** append my timings to a list of timings */
			static void collectTimings(TSFArray<TSFTimer>& timers) ;

			

			
		private:


			Ifpack_IlukGraph* precondGraph_;

			Ifpack_CrsRiluk* precond_;



		

			static TSFTimer opTimer_;
		};


}

#endif /* HAVE_PETRA */

#endif
