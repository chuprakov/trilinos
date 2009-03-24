#include "FEpetra_MpiComm.h"
#include "FEpetra_MpiComm_Cpp.hpp"
#include "Epetra_CombineMode.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

namespace {
  using Teuchos::Array;
  using Teuchos::RCP;

  Array<RCP<const Epetra_MpiComm> >& tableOfMpiComms(){
    static Array< RCP<const Epetra_MpiComm> > loc_tableOfMpiComms;
    return loc_tableOfMpiComms;
  }
} // namespace

//
// Definitions from Epetra_MpiComm.h
//

#ifdef __cplusplus
extern "C" {
#endif
  FEpetra_MpiCommID Epetra_MpiComm_Create(MPI_Comm comm){
    using Teuchos::rcp;
    tableOfMpiComms().push_back(rcp(new Epetra_MpiComm(comm)));
    return tableOfMpiComms().size() - 1;
  }
  
  void Epetra_Comm_Destroy( FEpetra_MpiCommID id ){
    tableOfMpiComms()[id] = Teuchos::null;
  }
#ifdef __cplusplus
} // extern "C"
#endif

//
// Definitions from Epetra_MpiComm_Cpp.hpp
//


const Teuchos::RCP<const Epetra_MpiComm>
FEpetra::getMpiComm( FEpetra_MpiCommID id ){
  return tableOfMpiComms()[id];
}
