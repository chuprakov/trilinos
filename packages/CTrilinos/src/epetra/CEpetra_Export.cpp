
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

#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Export.h"
#include "Epetra_Export.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Export */
Table<Epetra_Export>& tableOfExports()
{
    static Table<Epetra_Export>
        loc_tableOfExports(CT_Epetra_Export_ID, "CT_Epetra_Export_ID", false);
    return loc_tableOfExports;
}

/* table to hold objects of type const Epetra_Export */
Table<const Epetra_Export>& tableOfConstExports()
{
    static Table<const Epetra_Export>
        loc_tableOfConstExports(CT_Epetra_Export_ID, "CT_Epetra_Export_ID", true);
    return loc_tableOfConstExports;
}


} // namespace


//
// Definitions from CEpetra_Export.h
//


extern "C" {


CT_Epetra_Export_ID_t Epetra_Export_Cast ( CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstExports(), id);
    } else {
        newid = CTrilinos::cast(tableOfExports(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_Export_Abstract ( 
  CT_Epetra_Export_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id);
}

CT_Epetra_Export_ID_t Epetra_Export_Create ( 
  CT_Epetra_BlockMap_ID_t SourceMapID, 
  CT_Epetra_BlockMap_ID_t TargetMapID )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        tableOfExports().store(new Epetra_Export(
        *CEpetra::getConstBlockMap(SourceMapID), 
        *CEpetra::getConstBlockMap(TargetMapID))));
}

CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( 
  CT_Epetra_Export_ID_t ExporterID )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        tableOfExports().store(new Epetra_Export(
        *CEpetra::getConstExport(ExporterID))));
}

void Epetra_Export_Destroy ( CT_Epetra_Export_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_Export_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstExports().remove(&aid);
    } else {
        tableOfExports().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_Export_ID_t>(aid);
}

int Epetra_Export_NumSameIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->NumSameIDs();
}

int Epetra_Export_NumPermuteIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->NumPermuteIDs();
}

int * Epetra_Export_PermuteFromLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->PermuteFromLIDs();
}

int * Epetra_Export_PermuteToLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->PermuteToLIDs();
}

int Epetra_Export_NumRemoteIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->NumRemoteIDs();
}

int * Epetra_Export_RemoteLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->RemoteLIDs();
}

int Epetra_Export_NumExportIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->NumExportIDs();
}

int * Epetra_Export_ExportLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->ExportLIDs();
}

int * Epetra_Export_ExportPIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->ExportPIDs();
}

int Epetra_Export_NumSend ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->NumSend();
}

int Epetra_Export_NumRecv ( CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::getConstExport(selfID)->NumRecv();
}

CT_Epetra_BlockMap_ID_t Epetra_Export_SourceMap ( 
  CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstExport(
        selfID)->SourceMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_Export_TargetMap ( 
  CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstExport(
        selfID)->TargetMap() ));
}

CT_Epetra_Distributor_ID_t Epetra_Export_Distributor ( 
  CT_Epetra_Export_ID_t selfID )
{
    return CEpetra::storeDistributor(&( CEpetra::getConstExport(
        selfID)->Distributor() ));
}


} // extern "C"


//
// Definitions from CEpetra_Export_Cpp.hpp
//


/* get Epetra_Export from non-const table using CT_Epetra_Export_ID */
const Teuchos::RCP<Epetra_Export>
CEpetra::getExport( CT_Epetra_Export_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id);
    return tableOfExports().get(aid);
}

/* get Epetra_Export from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_Export>
CEpetra::getExport( CTrilinos_Object_ID_t id )
{
    return tableOfExports().get(id);
}

/* get const Epetra_Export from either the const or non-const table
 * using CT_Epetra_Export_ID */
const Teuchos::RCP<const Epetra_Export>
CEpetra::getConstExport( CT_Epetra_Export_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id);
    if (id.is_const) {
        return tableOfConstExports().get(aid);
    } else {
        return tableOfExports().get(aid);
    }
}

/* get const Epetra_Export from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_Export>
CEpetra::getConstExport( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstExports().get(id);
    } else {
        return tableOfExports().get(id);
    }
}

/* store Epetra_Export in non-const table */
CT_Epetra_Export_ID_t
CEpetra::storeExport( Epetra_Export *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
            tableOfExports().storeShared(pobj));
}

/* store const Epetra_Export in const table */
CT_Epetra_Export_ID_t
CEpetra::storeConstExport( const Epetra_Export *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
            tableOfConstExports().storeShared(pobj));
}

/* dump contents of Epetra_Export and const Epetra_Export tables */
void
CEpetra::purgeExportTables(  )
{
    tableOfExports().purge();
    tableOfConstExports().purge();
}



