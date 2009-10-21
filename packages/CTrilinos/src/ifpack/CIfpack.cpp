
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
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Ifpack */
Table<Ifpack>& tableOfIfpacks()
{
    static Table<Ifpack>
        loc_tableOfIfpacks(CT_Ifpack_ID, "CT_Ifpack_ID", FALSE);
    return loc_tableOfIfpacks;
}

/* table to hold objects of type const Ifpack */
Table<const Ifpack>& tableOfConstIfpacks()
{
    static Table<const Ifpack>
        loc_tableOfConstIfpacks(CT_Ifpack_ID, "CT_Ifpack_ID", TRUE);
    return loc_tableOfConstIfpacks;
}


} // namespace


//
// Definitions from CIfpack.h
//


extern "C" {


CT_Ifpack_ID_t Ifpack_Create (  )
{
    return CTrilinos::concreteType<CT_Ifpack_ID_t>(
        tableOfIfpacks().store(new Ifpack()));
}

void Ifpack_Destroy ( CT_Ifpack_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Ifpack_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstIfpacks().remove(&aid);
    } else {
        tableOfIfpacks().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Ifpack_ID_t>(aid);
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

const char * Ifpack_toString ( const  CT_EPrecType_E_t precType )
{
    return Ifpack::toString((Ifpack::EPrecType) precType);
}

CT_Ifpack_Preconditioner_ID_t Ifpack_CreatePreconditioner_UsingType ( 
  CT_EPrecType_E_t PrecType, CT_Epetra_RowMatrix_ID_t MatrixID, 
  const int overlap )
{
    return CIfpack::storePreconditioner(Ifpack::Create(
        (Ifpack::EPrecType) PrecType, CEpetra::getRowMatrix(
        MatrixID).getRawPtr(), overlap));
}


} // extern "C"


//
// Definitions from CIfpack_Cpp.hpp
//


/* get Ifpack from non-const table using CT_Ifpack_ID */
const Teuchos::RCP<Ifpack>
CIfpack::getIfpack( CT_Ifpack_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Ifpack_ID_t>(id);
    return tableOfIfpacks().get(aid);
}

/* get Ifpack from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Ifpack>
CIfpack::getIfpack( CTrilinos_Object_ID_t id )
{
    return tableOfIfpacks().get(id);
}

/* get const Ifpack from either the const or non-const table
 * using CT_Ifpack_ID */
const Teuchos::RCP<const Ifpack>
CIfpack::getConstIfpack( CT_Ifpack_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Ifpack_ID_t>(id);
    if (id.is_const) {
        return tableOfConstIfpacks().get(aid);
    } else {
        return tableOfIfpacks().get(aid);
    }
}

/* get const Ifpack from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Ifpack>
CIfpack::getConstIfpack( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstIfpacks().get(id);
    } else {
        return tableOfIfpacks().get(id);
    }
}

/* store Ifpack in non-const table */
CT_Ifpack_ID_t
CIfpack::storeIfpack( Ifpack *pobj )
{
    return CTrilinos::concreteType<CT_Ifpack_ID_t>(
            tableOfIfpacks().storeShared(pobj));
}

/* store const Ifpack in const table */
CT_Ifpack_ID_t
CIfpack::storeConstIfpack( const Ifpack *pobj )
{
    return CTrilinos::concreteType<CT_Ifpack_ID_t>(
            tableOfConstIfpacks().storeShared(pobj));
}

/* dump contents of Ifpack and const Ifpack tables */
void
CIfpack::purgeIfpackTables(  )
{
    tableOfIfpacks().purge();
    tableOfConstIfpacks().purge();
}



#endif /* HAVE_CTRILINOS_IFPACK */


