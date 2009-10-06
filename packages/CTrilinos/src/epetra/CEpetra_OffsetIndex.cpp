
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

#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_OffsetIndex.h"
#include "Epetra_OffsetIndex.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_OffsetIndex */
Table<Epetra_OffsetIndex>& tableOfOffsetIndexs()
{
    static Table<Epetra_OffsetIndex>
        loc_tableOfOffsetIndexs(CT_Epetra_OffsetIndex_ID, "CT_Epetra_OffsetIndex_ID", FALSE);
    return loc_tableOfOffsetIndexs;
}

/* table to hold objects of type const Epetra_OffsetIndex */
Table<const Epetra_OffsetIndex>& tableOfConstOffsetIndexs()
{
    static Table<const Epetra_OffsetIndex>
        loc_tableOfConstOffsetIndexs(CT_Epetra_OffsetIndex_ID, "CT_Epetra_OffsetIndex_ID", TRUE);
    return loc_tableOfConstOffsetIndexs;
}


} // namespace


//
// Definitions from CEpetra_OffsetIndex.h
//


extern "C" {


CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstOffsetIndexs(), id);
    } else {
        newid = CTrilinos::cast(tableOfOffsetIndexs(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_OffsetIndex_Abstract ( 
  CT_Epetra_OffsetIndex_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id);
}

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromImporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Import_ID_t ImporterID )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
        tableOfOffsetIndexs().store(new Epetra_OffsetIndex(
        *CEpetra::getConstCrsGraph(SourceGraphID), 
        *CEpetra::getConstCrsGraph(TargetGraphID), 
        *CEpetra::getImport(ImporterID))));
}

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromExporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Export_ID_t ExporterID )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
        tableOfOffsetIndexs().store(new Epetra_OffsetIndex(
        *CEpetra::getConstCrsGraph(SourceGraphID), 
        *CEpetra::getConstCrsGraph(TargetGraphID), 
        *CEpetra::getExport(ExporterID))));
}

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Duplicate ( 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
        tableOfOffsetIndexs().store(new Epetra_OffsetIndex(
        *CEpetra::getConstOffsetIndex(IndexorID))));
}

void Epetra_OffsetIndex_Destroy ( 
  CT_Epetra_OffsetIndex_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstOffsetIndexs().remove(&aid);
    } else {
        tableOfOffsetIndexs().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(aid);
}

int ** Epetra_OffsetIndex_SameOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID )
{
    return CEpetra::getConstOffsetIndex(selfID)->SameOffsets();
}

int ** Epetra_OffsetIndex_PermuteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID )
{
    return CEpetra::getConstOffsetIndex(selfID)->PermuteOffsets();
}

int ** Epetra_OffsetIndex_RemoteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID )
{
    return CEpetra::getConstOffsetIndex(selfID)->RemoteOffsets();
}


} // extern "C"


//
// Definitions from CEpetra_OffsetIndex_Cpp.hpp
//


/* get Epetra_OffsetIndex from non-const table using CT_Epetra_OffsetIndex_ID */
const Teuchos::RCP<Epetra_OffsetIndex>
CEpetra::getOffsetIndex( CT_Epetra_OffsetIndex_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id);
    return tableOfOffsetIndexs().get(aid);
}

/* get Epetra_OffsetIndex from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_OffsetIndex>
CEpetra::getOffsetIndex( CTrilinos_Object_ID_t id )
{
    return tableOfOffsetIndexs().get(id);
}

/* get const Epetra_OffsetIndex from either the const or non-const table
 * using CT_Epetra_OffsetIndex_ID */
const Teuchos::RCP<const Epetra_OffsetIndex>
CEpetra::getConstOffsetIndex( CT_Epetra_OffsetIndex_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id);
    if (id.is_const) {
        return tableOfConstOffsetIndexs().get(aid);
    } else {
        return tableOfOffsetIndexs().get(aid);
    }
}

/* get const Epetra_OffsetIndex from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_OffsetIndex>
CEpetra::getConstOffsetIndex( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstOffsetIndexs().get(id);
    } else {
        return tableOfOffsetIndexs().get(id);
    }
}

/* store Epetra_OffsetIndex in non-const table */
CT_Epetra_OffsetIndex_ID_t
CEpetra::storeOffsetIndex( Epetra_OffsetIndex *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
            tableOfOffsetIndexs().storeShared(pobj));
}

/* store const Epetra_OffsetIndex in const table */
CT_Epetra_OffsetIndex_ID_t
CEpetra::storeConstOffsetIndex( const Epetra_OffsetIndex *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
            tableOfConstOffsetIndexs().storeShared(pobj));
}

/* dump contents of Epetra_OffsetIndex and const Epetra_OffsetIndex tables */
void
CEpetra::purgeOffsetIndexTables(  )
{
    tableOfOffsetIndexs().purge();
    tableOfConstOffsetIndexs().purge();
}



