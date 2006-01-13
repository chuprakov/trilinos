#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "Teuchos_TestForException.hpp"

#include <iostream>

int meshreader(const Epetra_Comm & Comm,
               Epetra_IntSerialDenseVector & ipindx,
               Epetra_SerialDenseMatrix & ipcoords,
               Epetra_IntSerialDenseVector & pindx,
               Epetra_SerialDenseMatrix & pcoords,
               Epetra_IntSerialDenseMatrix & t,
               Epetra_IntSerialDenseMatrix & e,
               const char geomFileBase[],
               const bool trace,
               const bool dumpAll
               )
{
  int MyPID = Comm.MyPID();

  const int FileNameSize = 120;
  char FileName[FileNameSize];
  TEST_FOR_EXCEPT(static_cast<int>(std::strlen(geomFileBase) + 5) > FileNameSize);
  sprintf(FileName, "%s.%03d", geomFileBase, MyPID);

  FILE* fid;
  fid = fopen(FileName, "r");

  if(trace) printf("\nReading node info from %s ...\n", FileName);
  int numip = 0, numcp = 0; // # owned nodes, # shared nodes
  fscanf(fid, "%d %d", &numip, &numcp);
  const int nump = numip + numcp; // # total nodes
  if(trace) printf( "\nnumip = %d, numcp = %d, nump = %d\n", numip, numcp, nump );
  ipindx.Size(numip);
  ipcoords.Shape(numip, 2);
  pindx.Size(nump);
  pcoords.Shape(nump, 2);
  for (int i=0; i<numip; i++) {
    fscanf(fid,"%d %lf %lf %*d",&ipindx(i),&ipcoords(i,0),&ipcoords(i,1));
    if(trace&&dumpAll) printf("%d %lf %lf\n",ipindx(i),ipcoords(i,0),ipcoords(i,1));
    pindx(i) = ipindx(i);
    pcoords(i,0) = ipcoords(i,0); pcoords(i,1) = ipcoords(i,1);
  }
  for (int i=numip; i<nump; i++) {
    fscanf(fid,"%d %lf %lf %*d",&pindx(i),&pcoords(i,0),&pcoords(i,1));
  }

  fscanf(fid,"%*[^\n]");   // Skip to the End of the Line
  fscanf(fid,"%*1[\n]");   // Skip One Newline

  fscanf(fid,"%*[^\n]");   // Skip to the End of the Line
  fscanf(fid,"%*1[\n]");   // Skip One Newline

  for (int i=0; i<nump; i++) {
    fscanf(fid,"%*[^\n]"); // Skip to the End of the Line
    fscanf(fid,"%*1[\n]"); // Skip One Newline
  }

  if(trace) printf("\nReading element info from %s ...\n", FileName);
  int numelems = 0; // # elements on this processor
  fscanf(fid, "%d", &numelems);
  if(trace) printf( "\nnumelems = %d\n", numelems );
  t.Shape(numelems,3);
  for (int i=0; i<numelems; i++) {
    fscanf(fid, "%d %d %d", &t(i,0), &t(i,1), &t(i,2));
    if(trace&&dumpAll) printf("%d %d %d\n", t(i,0), t(i,1), t(i,2));
  }

  if(trace) printf("\nReading boundry edge info from %s ...\n", FileName);
  int numedges = 0; // # boundry edges on this processor
  fscanf(fid,"%d",&numedges);
  if(trace) printf( "\nnumedges = %d\n", numedges );
  e.Shape(numedges,3);
  for (int i=0; i<numedges; i++) {
    fscanf(fid, "%d %d %d", &e(i,0), &e(i,1), &e(i,2));
    if(trace&&dumpAll) printf("%d %d %d\n", e(i,0), e(i,1), e(i,2));
  }

  fclose(fid);
  if(trace) printf("Done reading mesh.\n");

  return(0);

}
