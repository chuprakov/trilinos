#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_COMM_H
#define CEPETRA_COMM_H


#include "CEpetra_Distributor.h"
#include "CEpetra_Directory.h"
#include "CEpetra_BlockMap.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_Comm_ID_t Epetra_Comm_Cast ( CTrilinos_Object_ID_t id );

/* Original C++ prototype:
   virtual Epetra_Comm * Clone() const = 0;
*/
CT_Epetra_Comm_ID_t Epetra_Comm_Clone ( CT_Epetra_Comm_ID_t selfID );

/* Original C++ prototype:
   virtual ~Epetra_Comm();
*/
void Epetra_Comm_Destroy ( CT_Epetra_Comm_ID_t * selfID );

/* Original C++ prototype:
   virtual void Barrier() const = 0;
*/
void Epetra_Comm_Barrier ( CT_Epetra_Comm_ID_t selfID );

/* Original C++ prototype:
   virtual int Broadcast(double * MyVals, int Count, int Root) const = 0;
*/
int Epetra_Comm_Broadcast_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, int Count, int Root );

/* Original C++ prototype:
   virtual int Broadcast(int * MyVals, int Count, int Root) const = 0;
*/
int Epetra_Comm_Broadcast_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int Count, int Root );

/* Original C++ prototype:
   virtual int Broadcast(long * MyVals, int Count, int Root) const = 0;
*/
int Epetra_Comm_Broadcast_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, int Count, int Root );

/* Original C++ prototype:
   virtual int Broadcast(char * MyVals, int Count, int Root) const = 0;
*/
int Epetra_Comm_Broadcast_Char ( 
  CT_Epetra_Comm_ID_t selfID, char * MyVals, int Count, int Root );

/* Original C++ prototype:
   virtual int GatherAll(double * MyVals, double * AllVals, int Count) const = 0;
*/
int Epetra_Comm_GatherAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count );

/* Original C++ prototype:
   virtual int GatherAll(int * MyVals, int * AllVals, int Count) const = 0;
*/
int Epetra_Comm_GatherAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int * AllVals, int Count );

/* Original C++ prototype:
   virtual int GatherAll(long * MyVals, long * AllVals, int Count) const = 0;
*/
int Epetra_Comm_GatherAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count );

/* Original C++ prototype:
   virtual int SumAll(double * PartialSums, double * GlobalSums, int Count) const = 0;
*/
int Epetra_Comm_SumAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count );

/* Original C++ prototype:
   virtual int SumAll(int * PartialSums, int * GlobalSums, int Count) const = 0;
*/
int Epetra_Comm_SumAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialSums, int * GlobalSums, 
  int Count );

/* Original C++ prototype:
   virtual int SumAll(long * PartialSums, long * GlobalSums, int Count) const = 0;
*/
int Epetra_Comm_SumAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialSums, long * GlobalSums, 
  int Count );

/* Original C++ prototype:
   virtual int MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const = 0;
*/
int Epetra_Comm_MaxAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count );

/* Original C++ prototype:
   virtual int MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const = 0;
*/
int Epetra_Comm_MaxAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialMaxs, int * GlobalMaxs, 
  int Count );

/* Original C++ prototype:
   virtual int MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const = 0;
*/
int Epetra_Comm_MaxAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialMaxs, long * GlobalMaxs, 
  int Count );

/* Original C++ prototype:
   virtual int MinAll(double * PartialMins, double * GlobalMins, int Count) const = 0;
*/
int Epetra_Comm_MinAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count );

/* Original C++ prototype:
   virtual int MinAll(int * PartialMins, int * GlobalMins, int Count) const = 0;
*/
int Epetra_Comm_MinAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialMins, int * GlobalMins, 
  int Count );

/* Original C++ prototype:
   virtual int MinAll(long * PartialMins, long * GlobalMins, int Count) const = 0;
*/
int Epetra_Comm_MinAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialMins, long * GlobalMins, 
  int Count );

/* Original C++ prototype:
   virtual int ScanSum(double * MyVals, double * ScanSums, int Count) const = 0;
*/
int Epetra_Comm_ScanSum_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, double * ScanSums, 
  int Count );

/* Original C++ prototype:
   virtual int ScanSum(int * MyVals, int * ScanSums, int Count) const = 0;
*/
int Epetra_Comm_ScanSum_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count );

/* Original C++ prototype:
   virtual int ScanSum(long * MyVals, long * ScanSums, int Count) const = 0;
*/
int Epetra_Comm_ScanSum_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count );

/* Original C++ prototype:
   virtual int MyPID() const = 0;
*/
int Epetra_Comm_MyPID ( CT_Epetra_Comm_ID_t selfID );

/* Original C++ prototype:
   virtual int NumProc() const = 0;
*/
int Epetra_Comm_NumProc ( CT_Epetra_Comm_ID_t selfID );

/* Original C++ prototype:
   virtual Epetra_Distributor * CreateDistributor() const = 0;
*/
CT_Epetra_Distributor_ID_t Epetra_Comm_CreateDistributor ( 
  CT_Epetra_Comm_ID_t selfID );

/* Original C++ prototype:
   virtual Epetra_Directory * CreateDirectory(const Epetra_BlockMap & Map) const = 0;
*/
CT_Epetra_Directory_ID_t Epetra_Comm_CreateDirectory ( 
  CT_Epetra_Comm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_COMM_H */

