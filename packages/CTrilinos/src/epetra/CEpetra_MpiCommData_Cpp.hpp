#ifdef HAVE_MPI


#include "CEpetra_MpiCommData.h"
#include "Epetra_MpiCommData.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_MPICOMMDATA_CPP_HPP
#define CEPETRA_MPICOMMDATA_CPP_HPP


class Epetra_MpiCommData;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_MpiCommData>
getMpiCommData( CT_Epetra_MpiCommData_ID_t id );

CT_Epetra_MpiCommData_ID_t
storeMpiCommData( const Epetra_MpiCommData *pobj );

void
purgeMpiCommDataTable(  );


} // namespace CEpetra


#endif // CEPETRA_MPICOMMDATA_CPP_HPP


#endif /* HAVE_MPI */


