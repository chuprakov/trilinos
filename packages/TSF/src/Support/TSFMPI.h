#ifndef TSFMPI_H
#define TSFMPI_H

#include "TSFConfig.h"

#if HAVE_MPI
#include "mpi.h"
#endif

namespace TSF
{
	/**
	 *
	 */
	class TSFMPI
		{
		public:
			/** initializer, calls MPI_Init() if necessary */
			static void init(int argc, void** argv);

			/** returns the process rank relative to MPI_COMM_WORLD */
			static int getRank() {return rank_;}

			/** returns the number of processors in MPI_COMM_WORLD */
			static int getNProc() {return nProc_;}

			/** finalizer, calls MPI_Finalize() if necessary */
			static void finalize();
		private:
			static int rank_;
			static int nProc_;
		};
}

#endif
