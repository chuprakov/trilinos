
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
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
//
// Definitions from CEpetra_Export.h
//


extern "C" {


CT_Epetra_Export_ID_t Epetra_Export_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Export_Generalize ( 
  CT_Epetra_Export_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id);
}

CT_Epetra_Export_ID_t Epetra_Export_Create ( 
  CT_Epetra_BlockMap_ID_t SourceMapID, 
  CT_Epetra_BlockMap_ID_t TargetMapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> SourceMap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(SourceMapID);
    const Teuchos::RCP<const Epetra_BlockMap> TargetMap = 
        CTrilinos::tableRepos().getConst<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(TargetMapID);
    return CTrilinos::tableRepos().store<Epetra_Export, CT_Epetra_Export_ID_t>(
        new Epetra_Export(*SourceMap, *TargetMap));
}

CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( 
  CT_Epetra_Export_ID_t ExporterID )
{
    const Teuchos::RCP<const Epetra_Export> Exporter = 
        CTrilinos::tableRepos().getConst<Epetra_Export, CT_Epetra_Export_ID_t>(
        ExporterID);
    return CTrilinos::tableRepos().store<Epetra_Export, CT_Epetra_Export_ID_t>(
        new Epetra_Export(*Exporter));
}

void Epetra_Export_Destroy ( CT_Epetra_Export_ID_t * selfID )
{
    CTrilinos::tableRepos().remove(selfID);
}

int Epetra_Export_NumSameIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->NumSameIDs();
}

int Epetra_Export_NumPermuteIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->NumPermuteIDs();
}

int * Epetra_Export_PermuteFromLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->PermuteFromLIDs();
}

int * Epetra_Export_PermuteToLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->PermuteToLIDs();
}

int Epetra_Export_NumRemoteIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->NumRemoteIDs();
}

int * Epetra_Export_RemoteLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->RemoteLIDs();
}

int Epetra_Export_NumExportIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->NumExportIDs();
}

int * Epetra_Export_ExportLIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->ExportLIDs();
}

int * Epetra_Export_ExportPIDs ( CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->ExportPIDs();
}

int Epetra_Export_NumSend ( CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->NumSend();
}

int Epetra_Export_NumRecv ( CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->NumRecv();
}

CT_Epetra_BlockMap_ID_t Epetra_Export_SourceMap ( 
  CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->SourceMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_Export_TargetMap ( 
  CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_BlockMap, 
        CT_Epetra_BlockMap_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->TargetMap() ));
}

CT_Epetra_Distributor_ID_t Epetra_Export_Distributor ( 
  CT_Epetra_Export_ID_t selfID )
{
    return CTrilinos::tableRepos().store<Epetra_Distributor, 
        CT_Epetra_Distributor_ID_t>(
        &( CTrilinos::tableRepos().getConst<Epetra_Export, 
        CT_Epetra_Export_ID_t>(selfID)->Distributor() ));
}


} // extern "C"


//
// Definitions from CEpetra_Export_Cpp.hpp
//


/* get Epetra_Export from non-const table using CT_Epetra_Export_ID */
const Teuchos::RCP<Epetra_Export>
CEpetra::getExport( CT_Epetra_Export_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Export, CT_Epetra_Export_ID_t>(id);
}

/* get Epetra_Export from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Export>
CEpetra::getExport( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Export, CTrilinos_Universal_ID_t>(id);
}

/* get const Epetra_Export from either the const or non-const table
 * using CT_Epetra_Export_ID */
const Teuchos::RCP<const Epetra_Export>
CEpetra::getConstExport( CT_Epetra_Export_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, CT_Epetra_Export_ID_t>(id);
}

/* get const Epetra_Export from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Export>
CEpetra::getConstExport( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export, CTrilinos_Universal_ID_t>(id);
}

/* store Epetra_Export in non-const table */
CT_Epetra_Export_ID_t
CEpetra::storeExport( Epetra_Export *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_Export, CT_Epetra_Export_ID_t>(pobj, false);
}

/* store const Epetra_Export in const table */
CT_Epetra_Export_ID_t
CEpetra::storeConstExport( const Epetra_Export *pobj )
{
    return CTrilinos::tableRepos().store<Epetra_Export, CT_Epetra_Export_ID_t>(pobj, false);
}



