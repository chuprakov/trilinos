// @HEADER

// @HEADER

#ifndef PLAYA_MPICOMM_H
#define PLAYA_MPICOMM_H

/*! \file PlayaMPIComm.hpp
    \brief Object representation of a MPI communicator
*/

#include "PlayaDefs.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#endif


namespace Playa
{
  /**
   * \brief Object representation of an MPI communicator.
   *
   * At present, groups are not implemented so the only communicator
   * is MPI_COMM_WORLD.
   */
  class TEUCHOS_LIB_DLL_EXPORT MPIComm
    {
    public:

      //! Empty constructor builds an object for MPI_COMM_WORLD
      MPIComm();

#ifdef HAVE_MPI
      //! Construct a MPIComm for a given MPI communicator
      MPIComm(MPI_Comm comm);
#endif

      //! Get an object representing MPI_COMM_WORLD 
      static MPIComm& world();
      //! Get an object representing MPI_COMM_SELF
      static MPIComm& self();

      //! Return process rank
      int getRank() const {return myRank_;}

      //! Return number of processors in the communicator
      int getNProc() const {return nProc_;}

      //! Synchronize all the processors in the communicator
      void synchronize() const ;

      //! @name Collective communications 
      //@{

      //! All-to-all gather-scatter
      void allToAll(void* sendBuf, int sendCount, int sendType,
                    void* recvBuf, int recvCount, int recvType) const ;

      //! Variable-length gather-scatter
      void allToAllv(void* sendBuf, int* sendCount, int* sendDisplacements,
                     int sendType,
                     void* recvBuf, int* recvCount,
                     int* recvDisplacements,
                     int recvType) const ;

      //! Do a collective operation, scattering the results to all processors
      void allReduce(void* input, void* result, int inputCount, int type,
                     int op) const ;


      //! Gather to root 
      void gather(void* sendBuf, int sendCount, int sendType,
                  void* recvBuf, int recvCount, int recvType,
                  int root) const ;

      //! Gather variable-sized arrays to root 
      void gatherv(void* sendBuf, int sendCount, int sendType,
                   void* recvBuf, int* recvCount, int* displacements, 
                   int recvType, int root) const ;

      //! Gather to all processors
      void allGather(void* sendBuf, int sendCount, int sendType,
                     void* recvBuf, int recvCount, int recvType) const ;

      //! Variable-length gather to all processors
      void allGatherv(void* sendBuf, int sendCount, int sendType,
                      void* recvBuf, int* recvCount, int* recvDisplacements,
                      int recvType) const ;

      //! Broadcast 
      void bcast(void* msg, int length, int type, int src) const ;

      //@}

#ifdef HAVE_MPI
      //! Get the MPI_Comm communicator handle 
      MPI_Comm getComm() const {return comm_;}
#endif

      //! @name Data types
      //@{ 
      //! Integer data type
      static const int INT;
      //! Float data type
      static const int FLOAT;
      //! Double data type
      static const int DOUBLE;
      //! Double/int structdata type
      static const int DOUBLE_INT;
      //! Character data type
      static const int CHAR;
      //@}

      //! @name Operations
      //@{ 
      //! Summation operation
      static const int SUM;
      //! Minimize operation
      static const int MIN;
      //! Maximize operation
      static const int MAX;
      //! Minimize operation
      static const int MINLOC;
      //! Maximize operation
      static const int MAXLOC;
      //! Dot-product (Multiplication) operation
      static const int PROD;
      //@}

      // errCheck() checks the return value of an MPI call and throws
      // a ParallelException upon failure.
      static void errCheck(int errCode, const std::string& methodName);

#ifdef HAVE_MPI
      //! Converts a PMachine data type code to a MPI_Datatype
      static MPI_Datatype getDataType(int type);

      //! Converts a PMachine operator code to a MPI_Op operator code.
      static MPI_Op getOp(int op);
#endif
    private:
#ifdef HAVE_MPI
      MPI_Comm comm_;
#endif

      int nProc_;
      int myRank_;

      /** common initialization function, called by all ctors */
      void init();

      /** Indicate whether MPI is currently running */
      int mpiIsRunning() const ;
    };
}
#endif

