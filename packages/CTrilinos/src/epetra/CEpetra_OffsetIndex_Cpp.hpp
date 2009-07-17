#include "CTrilinos_config.h"

#include "CEpetra_OffsetIndex.h"
#include "Epetra_OffsetIndex.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_OFFSETINDEX_CPP_HPP
#define CEPETRA_OFFSETINDEX_CPP_HPP


class Epetra_OffsetIndex;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_OffsetIndex>
getOffsetIndex( CT_Epetra_OffsetIndex_ID_t id );

const RCP<Epetra_OffsetIndex>
getOffsetIndex( CTrilinos_Object_ID_t id );

CT_Epetra_OffsetIndex_ID_t
storeOffsetIndex( const Epetra_OffsetIndex *pobj );

void
purgeOffsetIndexTable(  );


} // namespace CEpetra


#endif // CEPETRA_OFFSETINDEX_CPP_HPP


