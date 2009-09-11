
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

#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_DistObject_Cpp.hpp"
#include "CEpetra_DistObject.h"
#include "Epetra_DistObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_DistObject */
Table<Epetra_DistObject>& tableOfDistObjects()
{
    static Table<Epetra_DistObject>
        loc_tableOfDistObjects(CT_Epetra_DistObject_ID, "CT_Epetra_DistObject_ID", false);
    return loc_tableOfDistObjects;
}

/* table to hold objects of type const Epetra_DistObject */
Table<const Epetra_DistObject>& tableOfConstDistObjects()
{
    static Table<const Epetra_DistObject>
        loc_tableOfConstDistObjects(CT_Epetra_DistObject_ID, "CT_Epetra_DistObject_ID", true);
    return loc_tableOfConstDistObjects;
}


} // namespace


//
// Definitions from CEpetra_DistObject.h
//


extern "C" {


CT_Epetra_DistObject_ID_t Epetra_DistObject_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstDistObjects(), id);
    } else {
        newid = CTrilinos::cast(tableOfDistObjects(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_DistObject_Abstract ( 
  CT_Epetra_DistObject_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id);
}

void Epetra_DistObject_Destroy ( CT_Epetra_DistObject_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstDistObjects().remove(&aid);
    } else {
        tableOfDistObjects().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(aid);
}

int Epetra_DistObject_Import ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    return CEpetra::getDistObject(selfID)->Import(
        *CEpetra::getConstSrcDistObject(AID), 
        *CEpetra::getConstImport(ImporterID), 
        (Epetra_CombineMode) CombineMode, 
        CEpetra::getConstOffsetIndex(IndexorID).getRawPtr());
}

int Epetra_DistObject_Import_UsingExporter ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    return CEpetra::getDistObject(selfID)->Import(
        *CEpetra::getConstSrcDistObject(AID), 
        *CEpetra::getConstExport(ExporterID), 
        (Epetra_CombineMode) CombineMode, 
        CEpetra::getConstOffsetIndex(IndexorID).getRawPtr());
}

int Epetra_DistObject_Export_UsingImporter ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    return CEpetra::getDistObject(selfID)->Export(
        *CEpetra::getConstSrcDistObject(AID), 
        *CEpetra::getConstImport(ImporterID), 
        (Epetra_CombineMode) CombineMode, 
        CEpetra::getConstOffsetIndex(IndexorID).getRawPtr());
}

int Epetra_DistObject_Export ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    return CEpetra::getDistObject(selfID)->Export(
        *CEpetra::getConstSrcDistObject(AID), 
        *CEpetra::getConstExport(ExporterID), 
        (Epetra_CombineMode) CombineMode, 
        CEpetra::getConstOffsetIndex(IndexorID).getRawPtr());
}

CT_Epetra_BlockMap_ID_t Epetra_DistObject_Map ( 
  CT_Epetra_DistObject_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstDistObject(
        selfID)->Map() ));
}

CT_Epetra_Comm_ID_t Epetra_DistObject_Comm ( 
  CT_Epetra_DistObject_ID_t selfID )
{
    return CEpetra::storeConstComm(&( CEpetra::getConstDistObject(
        selfID)->Comm() ));
}

boolean Epetra_DistObject_DistributedGlobal ( 
  CT_Epetra_DistObject_ID_t selfID )
{
    return CEpetra::getConstDistObject(selfID)->DistributedGlobal();
}


} // extern "C"


//
// Definitions from CEpetra_DistObject_Cpp.hpp
//


/* get Epetra_DistObject from non-const table using CT_Epetra_DistObject_ID */
const Teuchos::RCP<Epetra_DistObject>
CEpetra::getDistObject( CT_Epetra_DistObject_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id);
    return tableOfDistObjects().get(aid);
}

/* get Epetra_DistObject from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_DistObject>
CEpetra::getDistObject( CTrilinos_Object_ID_t id )
{
    return tableOfDistObjects().get(id);
}

/* get const Epetra_DistObject from either the const or non-const table
 * using CT_Epetra_DistObject_ID */
const Teuchos::RCP<const Epetra_DistObject>
CEpetra::getConstDistObject( CT_Epetra_DistObject_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id);
    if (id.is_const) {
        return tableOfConstDistObjects().get(aid);
    } else {
        return tableOfDistObjects().get(aid);
    }
}

/* get const Epetra_DistObject from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_DistObject>
CEpetra::getConstDistObject( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstDistObjects().get(id);
    } else {
        return tableOfDistObjects().get(id);
    }
}

/* store Epetra_DistObject in non-const table */
CT_Epetra_DistObject_ID_t
CEpetra::storeDistObject( Epetra_DistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
            tableOfDistObjects().storeShared(pobj));
}

/* store const Epetra_DistObject in const table */
CT_Epetra_DistObject_ID_t
CEpetra::storeConstDistObject( const Epetra_DistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
            tableOfConstDistObjects().storeShared(pobj));
}

/* dump contents of Epetra_DistObject and const Epetra_DistObject tables */
void
CEpetra::purgeDistObjectTables(  )
{
    tableOfDistObjects().purge();
    tableOfConstDistObjects().purge();
}



