#ifndef TSFMPIVECTORSPACE_H
#define TSFMPIVECTORSPACE_H

#include "TSFConfig.h"
#include "TSFVectorSpaceBase.h"
#include <string>

#if HAVE_MPI
#include "mpi.h"
#endif


namespace TSF
{
	using std::string;

	class TSFVectorSpace;
	class TSFVectorBase;
	


	/** 
	 * Base class for vector spaces that use MPI
	 */

	class TSFMPIVectorSpace : public TSFVectorSpaceBase
		{
		public: 
			/** empty ctor */
			TSFMPIVectorSpace();

#if HAVE_MPI
			/** construct with an MPI communicator */
			TSFMPIVectorSpace(const MPI_Comm& comm);
#endif
			
			/** the usual virtual dtor */
			virtual ~TSFMPIVectorSpace(){;}
	
#if HAVE_MPI
			/** access to MPI communicator */
			virtual const MPI_Comm& getMPIComm() const {return comm_;}

			/** access to MPI communicator */
			virtual void setMPIComm(const MPI_Comm& comm) {comm_ = comm;}
#endif
				protected:
#if HAVE_MPI		
			MPI_Comm comm_;
#endif
		};
}

#endif
