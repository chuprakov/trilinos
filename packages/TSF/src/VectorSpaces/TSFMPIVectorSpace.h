#ifndef TSFMPIVECTORSPACE_H
#define TSFMPIVECTORSPACE_H

#include "TSFDefs.h"
#include "TSFVectorSpaceBase.h"
#include <string>

#ifdef HAVE_MPI
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

#ifdef HAVE_MPI
      /** construct with an MPI communicator */
      TSFMPIVectorSpace(const MPI_Comm& comm);
#endif

      /** the usual virtual dtor */
      virtual ~TSFMPIVectorSpace(){;}

#ifdef HAVE_MPI
      /** access to MPI communicator */
      virtual const MPI_Comm& getMPIComm() const {return comm_;}

      /** access to MPI communicator */
      virtual void setMPIComm(const MPI_Comm& comm) {comm_ = comm;}
#endif
    protected:
#ifdef HAVE_MPI
      MPI_Comm comm_;
#endif
    };
}

#endif
