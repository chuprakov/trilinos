/*! \@HEADER */
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
Questions? Contact M. Nicole Lemaster (mnlemas\@sandia.gov)

************************************************************************
*/
/*! \@HEADER */


#include "CTrilinos_config.h"
#include "CEpetra_Comm.h"
#include "CEpetra_BlockMap.h"
#include "Epetra_Map.h"
#include "CEpetra_Map.h"
#include "CEpetra_Map_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"

#include "CTrilinos_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Map_Cast ( CTrilinos_Object_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Cast )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 13);
  ECHO(int ElementSize = 4);
  ECHO(int IndexBase = 0);

  /* This cast should be allowed */
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create(
       NumGlobalElements, IndexBase, CommID));
  ECHO(CT_Epetra_Map_ID_t dupID = Epetra_Map_Cast(Epetra_Map_Abstract(selfID)));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), true);

  /* These casts should be allowed */
  ECHO(CT_Epetra_BlockMap_ID_t bmapID = Epetra_BlockMap_Cast(Epetra_Map_Abstract(selfID)));
  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(bmapID), NumGlobalElements);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, bmapID), true);
  ECHO(CT_Epetra_Map_ID_t mapID = Epetra_Map_Cast(Epetra_BlockMap_Abstract(bmapID)));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(bmapID, mapID), true);

  /* This cast should not be allowed */
  ECHO(CT_Epetra_BlockMap_ID_t bmapID2 = Epetra_BlockMap_Create(
       NumGlobalElements, ElementSize, IndexBase, CommID));
  TEST_THROW(Epetra_Map_Cast(Epetra_BlockMap_Abstract(bmapID2)), Teuchos::m_bad_cast);
}

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Map_Create ( 
  int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumGlobalElements = 9);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create(
       NumGlobalElements, IndexBase, CommID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Map_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Create_Linear )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(int NumMyElements = 5);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create_Linear(
       NumGlobalElements, NumMyElements, IndexBase, CommID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Map_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Create_Arbitrary )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());

  ECHO(const int NumMyElements = 4);
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Epetra_Map_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 4);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  ECHO(CT_Epetra_Map_ID_t dupID = Epetra_Map_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.type, CT_Epetra_Map_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
void Epetra_Map_Destroy ( CT_Epetra_Map_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 6);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  ECHO(Epetra_Map_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
void Epetra_Map_Assign ( 
  CT_Epetra_Map_ID_t selfID, CT_Epetra_Map_ID_t mapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Map , Assign )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create a map to duplicate */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int IndexBase = 0);
  ECHO(int NumGlobalElements1 = 4);
  ECHO(CT_Epetra_Map_ID_t  mapID = Epetra_Map_Create(NumGlobalElements1, IndexBase, CommID));

  /* Create the one to operate on */
  ECHO(int NumGlobalElements2 = 5);
  ECHO(CT_Epetra_Map_ID_t selfID = Epetra_Map_Create(NumGlobalElements2, IndexBase, CommID));

  /* Check the initial state */
  ECHO(CT_Epetra_BlockMap_ID_t bselfID = Epetra_BlockMap_Cast(Epetra_Map_Abstract(selfID)));
  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(bselfID), NumGlobalElements2);

  /* Test out the wrapper and check that it worked */
  ECHO(Epetra_Map_Assign(selfID, mapID));
  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(bselfID), NumGlobalElements1);
}

/**********************************************************************/

//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( T ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  DEBUG_UNIT_TEST_GROUP( T )


} // namespace

