
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_Time.h"

#ifdef EPETRA_MPI
#include <time.h>
#endif
#ifdef Epetra_HAVE_OMP
#include <omp.h>
#endif

//=============================================================================
Epetra_Time::Epetra_Time(const Epetra_Comm& Comm) 
  : StartTime_(0.0),
    Comm_(&Comm)
{
  StartTime_ = WallTime();
}
//=============================================================================
Epetra_Time::Epetra_Time(const Epetra_Time& Time) 
  : StartTime_(Time.StartTime_),
    Comm_(Time.Comm_)
{
}
//=============================================================================
Epetra_Time::~Epetra_Time(void)  
{
}
//=============================================================================
double Epetra_Time::WallTime(void) const
{
#ifdef EPETRA_MPI

	int mpiInitialized;
	MPI_Initialized(&mpiInitialized);

	if( mpiInitialized ) {

		return(MPI_Wtime());

	}
	else {

		clock_t start;

		start = clock();
		return( (double)( start ) / CLOCKS_PER_SEC );

	}

#else

#ifdef Epetra_HAVE_OMP
       return(omp_get_wtime());
#else
#if ICL || defined(_WIN32)

   clock_t start;
   //double duration;

   start = clock();
  return (double)( start ) / CLOCKS_PER_SEC;

#else

#ifndef MINGW
   struct timeval tp;
   static long start=0, startu;
   if (!start)
   {
      gettimeofday(&tp, NULL);
      start = tp.tv_sec;
      startu = tp.tv_usec;
      return(0.0);
   }
   gettimeofday(&tp, NULL);
   return( ((double) (tp.tv_sec - start)) + (tp.tv_usec-startu)/1000000.0 );
#else
   return( (double) clock() / CLOCKS_PER_SEC );
#endif // MINGW

#endif // ICL || WIN32

#endif // Epetra_HAVE_OMP
#endif // EPETRA_MPI

}
//=============================================================================
void Epetra_Time::ResetStartTime(void)
{
  StartTime_ = WallTime();
  return;
}
//=============================================================================
double Epetra_Time::ElapsedTime(void) const
{
  return(WallTime()-StartTime_);
}
