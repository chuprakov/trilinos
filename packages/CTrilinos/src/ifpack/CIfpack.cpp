
/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"


#ifdef HAVE_CTRILINOS_IFPACK


#include "CIfpack_Preconditioner_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CIfpack_Cpp.hpp"
#include "CIfpack.h"
#include "Ifpack.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CIfpack.h
//


extern "C" {


CT_Ifpack_ID_t Ifpack_Create (  )
{
    return CTrilinos::tableRepos().store<Ifpack, CT_Ifpack_ID_t>(
        new Ifpack());
}

void Ifpack_Destroy ( CT_Ifpack_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

CT_Ifpack_Preconditioner_ID_t Ifpack_CreatePreconditioner_UsingName ( 
  CT_Ifpack_ID_t selfID, const char PrecType[], 
  CT_Epetra_RowMatrix_ID_t MatrixID, const int overlap )
{
    return CIfpack::storePreconditioner(CIfpack::getIfpack(
        selfID)->Create(std::string(PrecType), CEpetra::getRowMatrix(
        MatrixID).getRawPtr(), overlap));
}

int Ifpack_SetParameters ( 
  CT_Ifpack_ID_t selfID, int argc, char * argv[], 
  CT_Teuchos_ParameterList_ID_t ListID, char * PrecType[], 
  int * Overlap )
{
    std::string *tmp_PrecType = NULL;
    CTrilinos::pass_string_in(PrecType, tmp_PrecType);
    int ret = CIfpack::getIfpack(selfID)->SetParameters(argc, argv, 
        *CTeuchos::getParameterList(ListID), *tmp_PrecType, 
        *Overlap);
    CTrilinos::pass_string_out(tmp_PrecType, PrecType);
    delete tmp_PrecType;

    return ret;
}

const char * Ifpack_toString ( const CT_EPrecType_E_t precType )
{
    return Ifpack::toString(CTrilinos::convert_to_difficult_enum(
        precType));
}

CT_Ifpack_Preconditioner_ID_t Ifpack_CreatePreconditioner_UsingType ( 
  CT_EPrecType_E_t PrecType, CT_Epetra_RowMatrix_ID_t MatrixID, 
  const int overlap )
{
    return CIfpack::storePreconditioner(Ifpack::Create(
        CTrilinos::convert_to_difficult_enum(PrecType), 
        CEpetra::getRowMatrix(MatrixID).getRawPtr(), overlap));
}


} // extern "C"


//
// Definitions from CIfpack_Cpp.hpp
//


/* get Ifpack from non-const table using CT_Ifpack_ID */
const Teuchos::RCP<Ifpack>
CIfpack::getIfpack( CT_Ifpack_ID_t id )
{
    return CTrilinos::tableRepos().get<Ifpack, CT_Ifpack_ID_t>(id);
}

/* get const Ifpack from either the const or non-const table
 * using CT_Ifpack_ID */
const Teuchos::RCP<const Ifpack>
CIfpack::getConstIfpack( CT_Ifpack_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Ifpack, CT_Ifpack_ID_t>(id);
}

/* store Ifpack in non-const table */
CT_Ifpack_ID_t
CIfpack::storeIfpack( Ifpack *pobj )
{
    return CTrilinos::tableRepos().store<Ifpack, CT_Ifpack_ID_t>(pobj, false);
}

/* store const Ifpack in const table */
CT_Ifpack_ID_t
CIfpack::storeConstIfpack( const Ifpack *pobj )
{
    return CTrilinos::tableRepos().store<Ifpack, CT_Ifpack_ID_t>(pobj, false);
}



#endif /* HAVE_CTRILINOS_IFPACK */


