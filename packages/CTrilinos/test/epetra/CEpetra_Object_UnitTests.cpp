#include <cstring>

#include "Epetra_Object.h"
#include "CEpetra_Object.h"
#include "CEpetra_Object_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"

#include "CEpetra_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_Object_ID_t Epetra_Object_Cast ( 
  CTrilinos_Object_ID_t id );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , Cast )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(int TracebackModeIn = -1);
  ECHO(bool set_label = true);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create(TracebackModeIn, set_label));

  /* This cast should be allowed */
  ECHO(CT_Epetra_Object_ID_t dupID = Epetra_Object_Cast(Epetra_Object_Abstract(selfID)));
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), true);
}

/**********************************************************************
CT_Epetra_Object_ID_t Epetra_Object_Create ( 
  int TracebackModeIn, boolean set_label );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(int TracebackModeIn = -1);
  ECHO(bool set_label = true);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create(TracebackModeIn, set_label));

  TEST_EQUALITY(selfID.type, CT_Epetra_Object_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( 
  const char * const Label, int TracebackModeIn );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , Create_WithLabel )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(const char * Label = "blah");
  ECHO(int TracebackModeIn = -1);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create_WithLabel(Label, TracebackModeIn));

  TEST_EQUALITY(selfID.type, CT_Epetra_Object_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( 
  CT_Epetra_Object_ID_t ObjectID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(int TracebackModeIn = -1);
  ECHO(bool set_label = true);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create(TracebackModeIn, set_label));

  ECHO(CT_Epetra_Object_ID_t dupID = Epetra_Object_Duplicate(selfID));

  TEST_EQUALITY(dupID.type, CT_Epetra_Object_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(int TracebackModeIn = -1);
  ECHO(bool set_label = true);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create(TracebackModeIn, set_label));

  ECHO(Epetra_Object_Destroy(&selfID));

  TEST_EQUALITY(selfID.type, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
void Epetra_Object_SetLabel ( 
  CT_Epetra_Object_ID_t selfID, const char * const Label );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , SetLabel )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(int TracebackModeIn = -1);
  ECHO(bool set_label = true);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create(TracebackModeIn, set_label));

  ECHO(const char * Label = "blah");
  ECHO(Epetra_Object_SetLabel(selfID, Label));

  ECHO(const char * Label2 = Epetra_Object_Label(selfID));
  TEST_EQUALITY_CONST(strcmp(Label, Label2), 0);
}

/**********************************************************************
const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Object , Label )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(const char * Label = "blah");
  ECHO(int TracebackModeIn = -1);
  ECHO(CT_Epetra_Object_ID_t selfID = Epetra_Object_Create_WithLabel(Label, TracebackModeIn));

  ECHO(const char * Label2 = Epetra_Object_Label(selfID));
  TEST_EQUALITY_CONST(strcmp(Label, Label2), 0);
}

/**********************************************************************
int Epetra_Object_ReportError ( 
  CT_Epetra_Object_ID_t selfID, const char * Message, int ErrorCode );
 **********************************************************************/



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

