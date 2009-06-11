#include "CEpetra_BlockMapData.h"
#include "Epetra_BlockMapData.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_BLOCKMAPDATA_CPP_HPP
#define CEPETRA_BLOCKMAPDATA_CPP_HPP


class Epetra_BlockMapData;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_BlockMapData>
getBlockMapData( CT_Epetra_BlockMapData_ID_t id );

CT_Epetra_BlockMapData_ID_t
storeBlockMapData( const Epetra_BlockMapData *pobj );

void
purgeBlockMapDataTable(  );


} // namespace CEpetra


#endif // CEPETRA_BLOCKMAPDATA_CPP_HPP


