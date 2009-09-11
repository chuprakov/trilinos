
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

#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Distributor.h"
#include "Epetra_Distributor.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Distributor */
Table<Epetra_Distributor>& tableOfDistributors()
{
    static Table<Epetra_Distributor>
        loc_tableOfDistributors(CT_Epetra_Distributor_ID, "CT_Epetra_Distributor_ID", false);
    return loc_tableOfDistributors;
}

/* table to hold objects of type const Epetra_Distributor */
Table<const Epetra_Distributor>& tableOfConstDistributors()
{
    static Table<const Epetra_Distributor>
        loc_tableOfConstDistributors(CT_Epetra_Distributor_ID, "CT_Epetra_Distributor_ID", true);
    return loc_tableOfConstDistributors;
}


} // namespace


//
// Definitions from CEpetra_Distributor.h
//


extern "C" {


CT_Epetra_Distributor_ID_t Epetra_Distributor_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstDistributors(), id);
    } else {
        newid = CTrilinos::cast(tableOfDistributors(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_Distributor_Abstract ( 
  CT_Epetra_Distributor_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id);
}

void Epetra_Distributor_Destroy ( 
  CT_Epetra_Distributor_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstDistributors().remove(&aid);
    } else {
        tableOfDistributors().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(aid);
}

CT_Epetra_Distributor_ID_t Epetra_Distributor_Clone ( 
  CT_Epetra_Distributor_ID_t selfID )
{
    return CEpetra::storeDistributor(CEpetra::getDistributor(
        selfID)->Clone());
}

int Epetra_Distributor_CreateFromSends ( 
  CT_Epetra_Distributor_ID_t selfID, int NumExportIDs, 
  const int * ExportPIDs, boolean Deterministic, 
  int * NumRemoteIDs )
{
    return CEpetra::getDistributor(selfID)->CreateFromSends(
        NumExportIDs, ExportPIDs, Deterministic, *NumRemoteIDs);
}

int Epetra_Distributor_CreateFromRecvs ( 
  CT_Epetra_Distributor_ID_t selfID, int NumRemoteIDs, 
  const int * RemoteGIDs, const int * RemotePIDs, 
  boolean Deterministic, int * NumExportIDs, int ** ExportGIDs, 
  int ** ExportPIDs )
{
    return CEpetra::getDistributor(selfID)->CreateFromRecvs(
        NumRemoteIDs, RemoteGIDs, RemotePIDs, Deterministic, 
        *NumExportIDs, *ExportGIDs, *ExportPIDs);
}

int Epetra_Distributor_Do ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->Do(export_objs, obj_size, 
        *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReverse ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReverse(export_objs, 
        obj_size, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoPosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoPosts(export_objs, 
        obj_size, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoWaits ( CT_Epetra_Distributor_ID_t selfID )
{
    return CEpetra::getDistributor(selfID)->DoWaits();
}

int Epetra_Distributor_DoReversePosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReversePosts(
        export_objs, obj_size, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReverseWaits ( 
  CT_Epetra_Distributor_ID_t selfID )
{
    return CEpetra::getDistributor(selfID)->DoReverseWaits();
}

int Epetra_Distributor_Do_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->Do(export_objs, obj_size, 
        *sizes, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReverse_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReverse(export_objs, 
        obj_size, *sizes, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoPosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoPosts(export_objs, 
        obj_size, *sizes, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReversePosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReversePosts(
        export_objs, obj_size, *sizes, *len_import_objs, 
        *import_objs);
}


} // extern "C"


//
// Definitions from CEpetra_Distributor_Cpp.hpp
//


/* get Epetra_Distributor from non-const table using CT_Epetra_Distributor_ID */
const Teuchos::RCP<Epetra_Distributor>
CEpetra::getDistributor( CT_Epetra_Distributor_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id);
    return tableOfDistributors().get(aid);
}

/* get Epetra_Distributor from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_Distributor>
CEpetra::getDistributor( CTrilinos_Object_ID_t id )
{
    return tableOfDistributors().get(id);
}

/* get const Epetra_Distributor from either the const or non-const table
 * using CT_Epetra_Distributor_ID */
const Teuchos::RCP<const Epetra_Distributor>
CEpetra::getConstDistributor( CT_Epetra_Distributor_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id);
    if (id.is_const) {
        return tableOfConstDistributors().get(aid);
    } else {
        return tableOfDistributors().get(aid);
    }
}

/* get const Epetra_Distributor from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_Distributor>
CEpetra::getConstDistributor( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstDistributors().get(id);
    } else {
        return tableOfDistributors().get(id);
    }
}

/* store Epetra_Distributor in non-const table */
CT_Epetra_Distributor_ID_t
CEpetra::storeDistributor( Epetra_Distributor *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(
            tableOfDistributors().storeShared(pobj));
}

/* store const Epetra_Distributor in const table */
CT_Epetra_Distributor_ID_t
CEpetra::storeConstDistributor( const Epetra_Distributor *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(
            tableOfConstDistributors().storeShared(pobj));
}

/* dump contents of Epetra_Distributor and const Epetra_Distributor tables */
void
CEpetra::purgeDistributorTables(  )
{
    tableOfDistributors().purge();
    tableOfConstDistributors().purge();
}



