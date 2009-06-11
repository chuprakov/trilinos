#include "CEpetra_CrsGraph.h"
#include "Epetra_CrsGraph.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_CRSGRAPH_CPP_HPP
#define CEPETRA_CRSGRAPH_CPP_HPP


class Epetra_CrsGraph;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_CrsGraph>
getCrsGraph( CT_Epetra_CrsGraph_ID_t id );

CT_Epetra_CrsGraph_ID_t
storeCrsGraph( const Epetra_CrsGraph *pobj );

void
purgeCrsGraphTable(  );


} // namespace CEpetra


#endif // CEPETRA_CRSGRAPH_CPP_HPP


