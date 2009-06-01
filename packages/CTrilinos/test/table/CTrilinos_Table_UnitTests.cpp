#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"
#include "Teuchos_RCP.hpp"

#include "Teuchos_UnitTestHarness.hpp"


#define JOIN_SET_0(A, B, C) A ## B ## C
#define JOIN_SET(A, B, C)   JOIN_SET_0(A, B, C)

#define BUILD_CALL(A, F) JOIN_SET( A , _ , F )
#define CLASS_TYPE(A)    JOIN_SET( CT_ , A , _ID_t )
#define CLASS_ENUM(A)    JOIN_SET( CT_ , A , _ID )
#define CLASS_ESTR(A)    XSTRFY(CLASS_ENUM(A))
#define STRFY(A)         #A
#define XSTRFY(A)        STRFY(A)
#define CONSTRUCTOR(A)   A, XSTRFY(A)


#define T Epetra_SerialComm
#define T1 Epetra_SerialComm
#define T2 Epetra_Comm
#define T3 Epetra_SerialComm
#define T4 Epetra_Vector


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::RangeError;
using Teuchos::NullReferenceError;
using Teuchos::m_bad_cast;
using CTrilinos::CTrilinosTypeMismatchError;
using CTrilinos::Table;


/* Table::store() */

TEUCHOS_UNIT_TEST( Table, store )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Object_ID_t id = table.store(new T));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.type, CLASS_ENUM(T));
  TEST_EQUALITY_CONST(nonnull(table.get(id)), true);
  TEST_EQUALITY_CONST(is_null(table.get(id)), false);
}

TEUCHOS_UNIT_TEST( Table, storeWrong )
{
  ECHO(Table<T2> table(CONSTRUCTOR(CLASS_ENUM(T2))));
  TEST_THROW(table.store(new T1), CTrilinosTypeMismatchError);
}

TEUCHOS_UNIT_TEST( Table, storeNull )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  TEST_THROW(table.store(NULL), NullReferenceError); 
}


/* Table::remove() */

TEUCHOS_UNIT_TEST( Table, remove )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Object_ID_t id = table.store(new T));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.type, CLASS_ENUM(T));
  ECHO(table.remove(&id));
  TEST_EQUALITY_CONST(id.index, -1);
  TEST_EQUALITY(id.type, CLASS_ENUM(Invalid));
}

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( Table, removeInvalid )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Object_ID_t id);
  ECHO(id.index = -1);
  ECHO(id.type = CLASS_ENUM(T));
  TEST_THROW(table.remove(&id), RangeError);
}

#endif /* TEUCHOS_DEBUG */

TEUCHOS_UNIT_TEST( Table, removeWrong )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Object_ID_t id1 = table.store(new T));
  ECHO(CTrilinos_Object_ID_t id);
  ECHO(id.index = id1.index);
  ECHO(id.type = CLASS_ENUM(T4));
  TEST_THROW(table.remove(&id), CTrilinosTypeMismatchError);
}


/* Table::get() */

TEUCHOS_UNIT_TEST( Table, get )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Object_ID_t id = table.store(new T));
  ECHO(RCP<T> rcpT = table.get(id));
  TEST_EQUALITY_CONST(nonnull(rcpT), true);
  TEST_EQUALITY_CONST(is_null(rcpT), false);
}

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( Table, getInvalid )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Object_ID_t id);
  ECHO(id.index = 0);
  ECHO(id.type = CLASS_ENUM(T));
  TEST_THROW(table.get(id), RangeError);
}

#endif /* TEUCHOS_DEBUG */

TEUCHOS_UNIT_TEST( Table, getWrong )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Object_ID_t id1 = table.store(new T));
  ECHO(CTrilinos_Object_ID_t id);
  ECHO(id.index = id1.index);
  ECHO(id.type = CLASS_ENUM(T4));
  TEST_THROW(table.get(id), CTrilinosTypeMismatchError);
}


/* Table::cast() */

TEUCHOS_UNIT_TEST( Table, cast )
{
  ECHO(Table<T1> table1(CONSTRUCTOR(CLASS_ENUM(T1))));
  ECHO(Table<T2> table2(CONSTRUCTOR(CLASS_ENUM(T2))));

  ECHO(CTrilinos_Object_ID_t id1 = table1.store(new T1));
  ECHO(CTrilinos_Object_ID_t id2 = table2.cast(table1.get(id1)));

  TEST_EQUALITY(id2.type, CLASS_ENUM(T2));
  TEST_EQUALITY_CONST(id2.index, 0);
}

TEUCHOS_UNIT_TEST( Table, castBad )
{
  ECHO(Table<T3> table3(CONSTRUCTOR(CLASS_ENUM(T3))));
  ECHO(Table<T4> table4(CONSTRUCTOR(CLASS_ENUM(T4))));

  ECHO(CTrilinos_Object_ID_t id3 = table3.store(new T3));
  TEST_THROW(table4.cast(table3.get(id3)), m_bad_cast);
}


/* Table::purge() */

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( Table, purge )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Object_ID_t id = table.store(new T));
  TEST_EQUALITY_CONST(nonnull(table.get(id)), true);
  ECHO(table.purge());
  TEST_THROW(table.get(id), RangeError);
}

#endif /* TEUCHOS_DEBUG */


//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( TT ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( TT )

#endif


#define UNIT_TEST_GROUP( TT ) \
  DEBUG_UNIT_TEST_GROUP( TT )


} // namespace
