#include "CTrilinos_config.h"

#include "CEpetra_Vector.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_VECTOR_CPP_HPP
#define CEPETRA_VECTOR_CPP_HPP


class Epetra_Vector;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_Vector>
getVector( CT_Epetra_Vector_ID_t id );

CT_Epetra_Vector_ID_t
storeVector( const Epetra_Vector *pobj );

void
purgeVectorTable(  );


} // namespace CEpetra


#endif // CEPETRA_VECTOR_CPP_HPP


