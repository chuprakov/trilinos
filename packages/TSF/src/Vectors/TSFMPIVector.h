#ifndef TSFMPIVECTOR_H
#define TSFMPIVECTOR_H

#include "TSFConfig.h"
#include "TSFInCoreVector.h"
#include <typeinfo>

#if HAVE_MPI
#include "mpi.h"
#endif

namespace TSF
{
	
	using std::string;

	/** \ingroup VectorSubtypes
	 *  TSFMPIVector is a base class for vectors that use MPI to do 
	 * reduction operations */

	class TSFMPIVector : public TSFInCoreVector
		{
		public: 
			/** \name Constructor and Destructors */
			//@{
			/** construct with a given space */
			TSFMPIVector(const TSFVectorSpace& space);

			/** the usual virtual dtor */
			virtual ~TSFMPIVector(){;}
			//@}

		protected:
			/** Reduce a calculation across processors */
			TSFReal reduce(const TSFReal& localValue, TSFReductionOp op) const ;

#if HAVE_MPI
			/** communicator */
			MPI_Comm comm_;
#endif

			/** */
			static TSFTimer& commTimer();
		};
};

#endif
