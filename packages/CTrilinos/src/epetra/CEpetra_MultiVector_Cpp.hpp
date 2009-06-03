#include "CEpetra_MultiVector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


#ifndef CEPETRA_MULTIVECTOR_CPP_HPP
#define CEPETRA_MULTIVECTOR_CPP_HPP


class Epetra_MultiVector;


namespace CEpetra {


using Teuchos::RCP;


const RCP<Epetra_MultiVector>
getMultiVector( CT_Epetra_MultiVector_ID_t id );

CT_Epetra_MultiVector_ID_t
storeMultiVector( const Epetra_MultiVector *pobj );

void
purgeMultiVectorTable(  );


} // namespace CEpetra


#endif // CEPETRA_MULTIVECTOR_CPP_HPP


