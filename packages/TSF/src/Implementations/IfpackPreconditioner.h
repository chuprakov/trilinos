#ifndef IFPACKPRECONDITIONER_H
#define IFPACKPRECONDITIONER_H

#include "TSFConfig.h"
#include "TSFVector.h"
#include "TSFLinearOperatorBase.h"

#if HAVE_PETRA
#include "Ifpack_CrsRiluk.h"
#include "Ifpack_IlukGraph.h"

namespace TSF
{
	using std::ostream;
	
	/** \ingroup Petra
	 * 
	 */
	
	class IfpackLeftPreconditioner : public TSFLinearOperatorBase 
		{
		public:
			/** empty ctor */
			IfpackLeftPreconditioner(const TSFVectorSpace& domain,
															 const TSFVectorSpace& range,
															 const TSFSmartPtr<Ifpack_CrsRiluk>& prec,
															 const TSFSmartPtr<Ifpack_IlukGraph>& graph);
			/** TUVC */
			virtual ~IfpackLeftPreconditioner(){;}

			/** apply operator to a vector in the domain space, returning a vector
			 * in the range space */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;

			

			
		private:
			TSFSmartPtr<Ifpack_CrsRiluk> precond_;

			TSFSmartPtr<Ifpack_IlukGraph> precondGraph_;
				
		};


}

#endif

#endif
