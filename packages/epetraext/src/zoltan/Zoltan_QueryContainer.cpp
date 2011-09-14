/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

//-------------------------------------------------------------------------
// Filename       : $Zoltan_QueryContainer.C$
//
// Purpose        : Static Container object to allow Static (C-style)
//                  functions to access methods of dynamic objects.
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra
//
// Creation Date  : 08/04/2000
//
// Revision Information:
// ---------------------
//
// Revision Number: $$
//
// Revision Date  : $$
//
// Current Owner  : $$
//-------------------------------------------------------------------------

#include <Zoltan_QueryContainer.h>

#include <Zoltan_QueryObject.h>


int Zoltan::QueryContainer::CurrentObject = 0;

std::map< int, Zoltan::QueryObject * > Zoltan::QueryContainer::StaticMap;

EPETRAEXT_DEPRECATED
void Zoltan::QueryContainer::setQueryID( const int & id )
{
  CurrentObject = id;
}

EPETRAEXT_DEPRECATED
const int & Zoltan::QueryContainer::getQueryID()
{
  return CurrentObject;
}

EPETRAEXT_DEPRECATED
bool Zoltan::QueryContainer::registerQueryObject( const int & id, 
		Zoltan::QueryObject * obj_ptr )
{
  if( StaticMap.find( id ) == StaticMap.end() )
  {
    StaticMap[ id ] = obj_ptr;
    return true;
  }
  else
  {
    // Redundant id
    return false;
  }
}

EPETRAEXT_DEPRECATED
Zoltan::QueryObject * Zoltan::QueryContainer::getQueryObject(
	const int & id )
{
  if( StaticMap.find( id ) != StaticMap.end() )
  {
    return StaticMap[ id ];
  }
  else
  {
    // Not found
    return 0;
  }
}

