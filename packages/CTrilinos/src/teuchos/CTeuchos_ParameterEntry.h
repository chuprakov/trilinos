#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CTEUCHOS_PARAMETERENTRY_H
#define CTEUCHOS_PARAMETERENTRY_H


#include "CTeuchos_any.h"
#include "CTeuchos_ParameterList.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Teuchos_ParameterEntry_Cast() and Teuchos_ParameterEntry_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Teuchos_ParameterEntry table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

/* Original C++ prototype:
   ParameterEntry();
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Create (  );

/* Original C++ prototype:
   ParameterEntry(const ParameterEntry& source);
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Duplicate ( 
  CT_Teuchos_ParameterEntry_ID_t sourceID );

/* Original C++ prototype:
   ~ParameterEntry();
*/
void Teuchos_ParameterEntry_Destroy ( 
  CT_Teuchos_ParameterEntry_ID_t * selfID );

/* Original C++ prototype:
   ParameterEntry& operator=(const ParameterEntry& source);
*/
void Teuchos_ParameterEntry_Assign ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, 
  CT_Teuchos_ParameterEntry_ID_t sourceID );

/* Original C++ prototype:
   void setAnyValue( const any &value, bool isDefault = false );
*/
void Teuchos_ParameterEntry_setAnyValue ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, 
  CT_Teuchos_any_ID_t valueID, boolean isDefault );

/* Original C++ prototype:
   void setDocString(const std::string &docString);
*/
void Teuchos_ParameterEntry_setDocString ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, char docString[] );

/* Original C++ prototype:
   ParameterList& setList( bool isDefault = false, const std::string &docString = "" );
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterEntry_setList ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean isDefault, 
  char docString[] );

/* Original C++ prototype:
   template<typename T> inline T& getValue(T *ptr) const;
*/
double Teuchos_ParameterEntry_getValue_double ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, double * ptr );

/* Original C++ prototype:
   template<typename T> inline T& getValue(T *ptr) const;
*/
int Teuchos_ParameterEntry_getValue_int ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, int * ptr );

/* Original C++ prototype:
   inline any& getAny(bool activeQry = true);
*/
CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean activeQry );

/* Original C++ prototype:
   inline const any& getAny(bool activeQry = true) const;
*/
CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny_const ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean activeQry );

/* Original C++ prototype:
   inline bool isUsed() const;
*/
boolean Teuchos_ParameterEntry_isUsed ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/* Original C++ prototype:
   bool isList() const;
*/
boolean Teuchos_ParameterEntry_isList ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/* Original C++ prototype:
   template <typename T> inline bool isType() const;
*/
boolean Teuchos_ParameterEntry_isType_double ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/* Original C++ prototype:
   template <typename T> inline bool isType() const;
*/
boolean Teuchos_ParameterEntry_isType_int ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/* Original C++ prototype:
   inline bool isDefault() const;
*/
boolean Teuchos_ParameterEntry_isDefault ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/* Original C++ prototype:
   inline std::string docString() const;
*/
const char * Teuchos_ParameterEntry_docString ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CTEUCHOS_PARAMETERENTRY_H */

