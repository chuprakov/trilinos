#include "CEpetra_CrsGraphData.h"
#include "Epetra_CrsGraphData.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_CRSGRAPHDATA_CPP_HPP
#define CEPETRA_CRSGRAPHDATA_CPP_HPP


class Epetra_CrsGraphData;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_CrsGraphData>
getCrsGraphData( CT_Epetra_CrsGraphData_ID_t id );

CT_Epetra_CrsGraphData_ID_t
storeCrsGraphData( const Epetra_CrsGraphData *pobj );

void
purgeCrsGraphDataTable(  );


} // namespace CEpetra


#endif // CEPETRA_CRSGRAPHDATA_CPP_HPP


