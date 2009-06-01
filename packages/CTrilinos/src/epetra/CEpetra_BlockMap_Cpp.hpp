#include "CEpetra_BlockMap.h"
#include "Epetra_BlockMap.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_BLOCKMAP_CPP_HPP
#define CEPETRA_BLOCKMAP_CPP_HPP


class Epetra_BlockMap;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_BlockMap>
getBlockMap( CT_Epetra_BlockMap_ID_t id );

void
purgeBlockMapTable(  );


} // namespace CEpetra


#endif // CEPETRA_BLOCKMAP_CPP_HPP


