#include "CTrilinos_config.h"


#ifdef HAVE_MPI



/* True C header file! */


#ifndef CEPETRA_MPICOMM_H
#define CEPETRA_MPICOMM_H


#include "mpi.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Distributor.h"
#include "CEpetra_Directory.h"
#include "CEpetra_BlockMap.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Cast ( 
  CTrilinos_Object_ID_t id );

/* Original C++ prototype:
   Epetra_MpiComm(MPI_Comm comm);
*/
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create ( MPI_Comm comm );

/* Original C++ prototype:
   Epetra_MpiComm(const Epetra_MpiComm & Comm);
*/
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( 
  CT_Epetra_MpiComm_ID_t CommID );

/* Original C++ prototype:
   Epetra_Comm * Clone() const;
*/
CT_Epetra_Comm_ID_t Epetra_MpiComm_Clone ( 
  CT_Epetra_MpiComm_ID_t selfID );

/* Original C++ prototype:
   virtual ~Epetra_MpiComm();
*/
void Epetra_MpiComm_Destroy ( CT_Epetra_MpiComm_ID_t * selfID );

/* Original C++ prototype:
   void Barrier() const;
*/
void Epetra_MpiComm_Barrier ( CT_Epetra_MpiComm_ID_t selfID );

/* Original C++ prototype:
   int Broadcast(double * MyVals, int Count, int Root) const;
*/
int Epetra_MpiComm_Broadcast_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, int Count, 
  int Root );

/* Original C++ prototype:
   int Broadcast(int * MyVals, int Count, int Root) const;
*/
int Epetra_MpiComm_Broadcast_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int Count, int Root );

/* Original C++ prototype:
   int Broadcast(long * MyVals, int Count, int Root) const;
*/
int Epetra_MpiComm_Broadcast_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, int Count, int Root );

/* Original C++ prototype:
   int Broadcast(char * MyVals, int Count, int Root) const;
*/
int Epetra_MpiComm_Broadcast_Char ( 
  CT_Epetra_MpiComm_ID_t selfID, char * MyVals, int Count, int Root );

/* Original C++ prototype:
   int GatherAll(double * MyVals, double * AllVals, int Count) const;
*/
int Epetra_MpiComm_GatherAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count );

/* Original C++ prototype:
   int GatherAll(int * MyVals, int * AllVals, int Count) const;
*/
int Epetra_MpiComm_GatherAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count );

/* Original C++ prototype:
   int GatherAll(long * MyVals, long * AllVals, int Count) const;
*/
int Epetra_MpiComm_GatherAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count );

/* Original C++ prototype:
   int SumAll(double * PartialSums, double * GlobalSums, int Count) const;
*/
int Epetra_MpiComm_SumAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count );

/* Original C++ prototype:
   int SumAll(int * PartialSums, int * GlobalSums, int Count) const;
*/
int Epetra_MpiComm_SumAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialSums, int * GlobalSums, 
  int Count );

/* Original C++ prototype:
   int SumAll(long * PartialSums, long * GlobalSums, int Count) const;
*/
int Epetra_MpiComm_SumAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count );

/* Original C++ prototype:
   int MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const;
*/
int Epetra_MpiComm_MaxAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count );

/* Original C++ prototype:
   int MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const;
*/
int Epetra_MpiComm_MaxAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMaxs, int * GlobalMaxs, 
  int Count );

/* Original C++ prototype:
   int MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const;
*/
int Epetra_MpiComm_MaxAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count );

/* Original C++ prototype:
   int MinAll(double * PartialMins, double * GlobalMins, int Count) const;
*/
int Epetra_MpiComm_MinAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count );

/* Original C++ prototype:
   int MinAll(int * PartialMins, int * GlobalMins, int Count) const;
*/
int Epetra_MpiComm_MinAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMins, int * GlobalMins, 
  int Count );

/* Original C++ prototype:
   int MinAll(long * PartialMins, long * GlobalMins, int Count) const;
*/
int Epetra_MpiComm_MinAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count );

/* Original C++ prototype:
   int ScanSum(double * MyVals, double * ScanSums, int Count) const;
*/
int Epetra_MpiComm_ScanSum_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * ScanSums, 
  int Count );

/* Original C++ prototype:
   int ScanSum(int * MyVals, int * ScanSums, int Count) const;
*/
int Epetra_MpiComm_ScanSum_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count );

/* Original C++ prototype:
   int ScanSum(long * MyVals, long * ScanSums, int Count) const;
*/
int Epetra_MpiComm_ScanSum_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count );

/* Original C++ prototype:
   MPI_Comm Comm() const;
*/
MPI_Comm Epetra_MpiComm_Comm ( CT_Epetra_MpiComm_ID_t selfID );

/* Original C++ prototype:
   int MyPID() const;
*/
int Epetra_MpiComm_MyPID ( CT_Epetra_MpiComm_ID_t selfID );

/* Original C++ prototype:
   int NumProc() const;
*/
int Epetra_MpiComm_NumProc ( CT_Epetra_MpiComm_ID_t selfID );

/* Original C++ prototype:
   Epetra_Distributor * CreateDistributor() const;
*/
CT_Epetra_Distributor_ID_t Epetra_MpiComm_CreateDistributor ( 
  CT_Epetra_MpiComm_ID_t selfID );

/* Original C++ prototype:
   Epetra_Directory * CreateDirectory(const Epetra_BlockMap & Map) const;
*/
CT_Epetra_Directory_ID_t Epetra_MpiComm_CreateDirectory ( 
  CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );

/* Original C++ prototype:
   int GetMpiTag() const;
*/
int Epetra_MpiComm_GetMpiTag ( CT_Epetra_MpiComm_ID_t selfID );

/* Original C++ prototype:
   MPI_Comm GetMpiComm() const;
*/
MPI_Comm Epetra_MpiComm_GetMpiComm ( CT_Epetra_MpiComm_ID_t selfID );

/* Original C++ prototype:
   Epetra_MpiComm & operator=(const Epetra_MpiComm & Comm);
*/
void Epetra_MpiComm_Assign ( 
  CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_MpiComm_ID_t CommID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_MPICOMM_H */

#endif /* HAVE_MPI */


