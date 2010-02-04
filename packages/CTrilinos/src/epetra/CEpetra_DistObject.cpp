
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
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_DistObject.h
//


extern "C" {


void Epetra_DistObject_Destroy ( CT_Epetra_DistObject_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
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
    return ((CEpetra::getConstDistObject(
        selfID)->DistributedGlobal()) ? TRUE : FALSE);
}


} // extern "C"


//
// Definitions from CEpetra_DistObject_Cpp.hpp
//


/* get Epetra_DistObject from non-const table using CT_Epetra_DistObject_ID */
const Teuchos::RCP<Epetra_DistObject>
CEpetra::getDistObject( CT_Epetra_DistObject_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_DistObject, CT_Epetra_DistObject_ID_t>(id);
}

/* get Epetra_DistObject from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_DistObject>
CEpetra::getDistObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_DistObject, CTrilinos_Universal_ID_t>(id);
}

/* get const Epetra_DistObject from either the const or non-const table
 * using CT_Epetra_DistObject_ID */
const Teuchos::RCP<const Epetra_DistObject>
CEpetra::getConstDistObject( CT_Epetra_DistObject_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_DistObject, CT_Epetra_DistObject_ID_t>(id);
}

/* get const Epetra_DistObject from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_DistObject>
CEpetra::getConstDistObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_DistObject, CTrilinos_Universal_ID_t>(id);
}

/* store Epetra_DistObject in non-const table */
CT_Epetra_DistObject_ID_t
CEpetra::storeDistObject( Epetra_DistObject *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_DistObject, CT_Epetra_DistObject_ID_t>(pobj, false);
}

/* store const Epetra_DistObject in const table */
CT_Epetra_DistObject_ID_t
CEpetra::storeConstDistObject( const Epetra_DistObject *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_DistObject, CT_Epetra_DistObject_ID_t>(pobj, false);
}



