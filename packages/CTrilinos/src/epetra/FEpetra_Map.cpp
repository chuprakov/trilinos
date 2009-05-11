#include "FEpetra_Map.h"
#include "FEpetra_Map_Cpp.hpp"
#include "Epetra_Map.h"
#include "Epetra_SerialComm.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

namespace {
  using Teuchos::Array;
  using Teuchos::RCP;

  // Note that by storing RCP objects that when the program goes down, any Maps
  // not destoryed explicitly will be destroyed here so valgrind will not
  // complain!  I don't know if this is a good thing or a bad thing?

  Array<RCP<const Epetra_Map> >& tableOfMaps(){
    static Array<RCP<const Epetra_Map> > loc_tableOfMaps;
    return loc_tableOfMaps;
  }
  // 2007/12/18: rabartl: ToDo: Validate that the above table is empty when the
  // program finishes!  
} // namespace

//
// Definitions from Epetra_Map.h
//

#ifdef __cplusplus
extern "C" {
#endif
  MapID Epetra_Map_Create( int numGlobalElements ){
    using Teuchos::rcp;
    Epetra_SerialComm comm;
    tableOfMaps().push_back(rcp(new Epetra_Map(numGlobalElements,0,comm)));
    return tableOfMaps().size() - 1;
  }

  void Epetra_Map_Destroy( MapID mapID ){
    tableOfMaps()[mapID] = Teuchos::null;
    // 2007/12/18: rabartl: ToDo: Maintain a free list of mapIDs that can be
    // reused.  This will allow for many many allocations and deallocations!
  }

  int Epetra_Map_NumGlobalElements( MapID mapID ){
    return FEpetra::getMap(mapID)->NumGlobalElements();
  }


  int Epetra_Map_NumMyElements( MapID mapID ){
    return FEpetra::getMap(mapID)->NumMyElements();
  }
#ifdef __cplusplus
} // extern "C"
#endif

//
// Definitions from Epetra_Map_Cpp.hpp
//

const Teuchos::RCP<const Epetra_Map>
FEpetra::getMap( MapID mapID ){
  return tableOfMaps()[mapID];
}
