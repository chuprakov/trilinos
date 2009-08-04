#include "CTrilinos_config.h"

#ifndef CEPETRA_CRSGRAPH_CPP_HPP
#define CEPETRA_CRSGRAPH_CPP_HPP


#include "CEpetra_CrsGraph.h"
#include "Epetra_CrsGraph.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


class Epetra_CrsGraph;


namespace CEpetra {


using Teuchos::RCP;


/* get Epetra_CrsGraph from non-const table using CT_Epetra_CrsGraph_ID */
const RCP<Epetra_CrsGraph>
getCrsGraph( CT_Epetra_CrsGraph_ID_t id );

/* get Epetra_CrsGraph from non-const table using CTrilinos_Object_ID_t */
const RCP<Epetra_CrsGraph>
getCrsGraph( CTrilinos_Object_ID_t id );

/* get const Epetra_CrsGraph from either the const or non-const table
 * using CT_Epetra_CrsGraph_ID */
const RCP<const Epetra_CrsGraph>
getConstCrsGraph( CT_Epetra_CrsGraph_ID_t id );

/* get const Epetra_CrsGraph from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Epetra_CrsGraph>
getConstCrsGraph( CTrilinos_Object_ID_t id );

/* store Epetra_CrsGraph in non-const table */
CT_Epetra_CrsGraph_ID_t
storeCrsGraph( Epetra_CrsGraph *pobj );

/* store const Epetra_CrsGraph in const table */
CT_Epetra_CrsGraph_ID_t
storeConstCrsGraph( const Epetra_CrsGraph *pobj );

/* dump contents of Epetra_CrsGraph and const Epetra_CrsGraph tables */
void
purgeCrsGraphTables(  );


} // namespace CEpetra


#endif // CEPETRA_CRSGRAPH_CPP_HPP


