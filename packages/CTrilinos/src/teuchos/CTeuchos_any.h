#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CTEUCHOS_ANY_H
#define CTEUCHOS_ANY_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Teuchos_any_Cast() and Teuchos_any_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Teuchos_any table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

/* Original C++ prototype:
   any();
*/
CT_Teuchos_any_ID_t Teuchos_any_Create (  );

/* Original C++ prototype:
   any(const any & other);
*/
CT_Teuchos_any_ID_t Teuchos_any_Duplicate ( 
  CT_Teuchos_any_ID_t otherID );

/* Original C++ prototype:
   ~any();
*/
void Teuchos_any_Destroy ( CT_Teuchos_any_ID_t * selfID );

/* Original C++ prototype:
   any & swap(any & rhs);
*/
CT_Teuchos_any_ID_t Teuchos_any_swap ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID );

/* Original C++ prototype:
   any & operator=(const any & rhs);
*/
void Teuchos_any_Assign ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID );

/* Original C++ prototype:
   bool empty() const;
*/
boolean Teuchos_any_empty ( CT_Teuchos_any_ID_t selfID );

/* Original C++ prototype:
   std::string typeName() const;
*/
const char * Teuchos_any_typeName ( CT_Teuchos_any_ID_t selfID );

/* Original C++ prototype:
   bool same( const any &other ) const;
*/
boolean Teuchos_any_same ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t otherID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CTEUCHOS_ANY_H */

