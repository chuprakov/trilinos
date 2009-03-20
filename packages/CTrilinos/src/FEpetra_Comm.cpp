#include "FEpetra_Comm.h"
#include "FEpetra_Comm_Cpp.hpp"
#include "Epetra_CombineMode.h"
#include "Epetra_Comm.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

namespace {
  using Teuchos::Array;
  using Teuchos::RCP;

  Array<RCP<const Epetra_Comm> >& tableOfComms(){
    static Array< RCP<const Epetra_Comm> > loc_tableOfComms;
    return loc_tableOfComms;
  }
} // namespace

//
// Definitions from Epetra_Comm.h
//

#ifdef __cplusplus
extern "C" {
#endif
  int Epetra_Comm_SumAll(double *PartialSums, double *GlobalSums,int Count,CommID comm){
    Epetra_Comm &epcomm = tableOfComms().getComm(comm);
    return epcomm.SumAll(PartialSums,GlobalSums,Count);
  }
  
#ifdef __cplusplus
} // extern "C"
#endif

//
// Definitions from Epetra_Comm_Cpp.hpp
//


const Teuchos::RCP<const Epetra_Comm>
FEpetra::getComm( CommID id ){
  return tableOfComms()[id];
}
