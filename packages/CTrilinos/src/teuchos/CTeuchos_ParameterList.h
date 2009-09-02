#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CTEUCHOS_PARAMETERLIST_H
#define CTEUCHOS_PARAMETERLIST_H


#include "CTeuchos_ParameterEntry.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Teuchos_ParameterList_Cast() and Teuchos_ParameterList_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Teuchos_ParameterList table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Teuchos_ParameterList_Abstract ( 
  CT_Teuchos_ParameterList_ID_t id );

/* Original C++ prototype:
   ParameterList();
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create (  );

/* Original C++ prototype:
   ParameterList(const std::string &name);
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_WithName ( 
  const char name[] );

/* Original C++ prototype:
   ParameterList(const ParameterList& source);
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_FromSource ( 
  CT_Teuchos_ParameterList_ID_t sourceID );

/* Original C++ prototype:
   virtual ~ParameterList();
*/
void Teuchos_ParameterList_Destroy ( 
  CT_Teuchos_ParameterList_ID_t * selfID );

/* Original C++ prototype:
   ParameterList& setName( const std::string &name );
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setName ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   ParameterList& operator=(const ParameterList& source);
*/
void Teuchos_ParameterList_Assign ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID );

/* Original C++ prototype:
   ParameterList& setParameters(const ParameterList& source);
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID );

/* Original C++ prototype:
   ParameterList& setParametersNotAlreadySet(const ParameterList& source);
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParametersNotAlreadySet ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID );

/* Original C++ prototype:
   ParameterList& disableRecursiveValidation();
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_disableRecursiveValidation ( 
  CT_Teuchos_ParameterList_ID_t selfID );

/* Original C++ prototype:
   ParameterList& set( std::string const& name, ParameterList const& value, std::string const& docString = "" );
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  CT_Teuchos_ParameterList_ID_t valueID, char const docString[] );

/* Original C++ prototype:
   ParameterList& setEntry(const std::string& name, const ParameterEntry& entry);
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  CT_Teuchos_ParameterEntry_ID_t entryID );

/* Original C++ prototype:
   template<typename T> T& get(const std::string& name, T def_value);
*/
double Teuchos_ParameterList_get_double_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double def_value );

/* Original C++ prototype:
   template<typename T> T& get(const std::string& name, T def_value);
*/
int Teuchos_ParameterList_get_int_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int def_value );

/* Original C++ prototype:
   std::string& get(const std::string& name, char def_value[]);
*/
const char * Teuchos_ParameterList_get_char_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  char def_value[] );

/* Original C++ prototype:
   std::string& get(const std::string& name, const char def_value[]);
*/
const char * Teuchos_ParameterList_get_const_char_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  const char def_value[] );

/* Original C++ prototype:
   template<typename T> T& get(const std::string& name);
*/
double Teuchos_ParameterList_get_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   template<typename T> T& get(const std::string& name);
*/
int Teuchos_ParameterList_get_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   template<typename T> const T& get(const std::string& name) const;
*/
double Teuchos_ParameterList_get_double_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   template<typename T> const T& get(const std::string& name) const;
*/
int Teuchos_ParameterList_get_int_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   template<typename T> inline T* getPtr(const std::string& name);
*/
double * Teuchos_ParameterList_getPtr_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   template<typename T> inline T* getPtr(const std::string& name);
*/
int * Teuchos_ParameterList_getPtr_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   template<typename T> inline const T* getPtr(const std::string& name) const;
*/
const double * Teuchos_ParameterList_getPtr_double_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   template<typename T> inline const T* getPtr(const std::string& name) const;
*/
const int * Teuchos_ParameterList_getPtr_int_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   ParameterEntry& getEntry(const std::string& name);
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   inline const ParameterEntry& getEntry(const std::string& name) const;
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   inline ParameterEntry* getEntryPtr(const std::string& name);
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   inline const ParameterEntry* getEntryPtr(const std::string& name) const;
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   bool remove( std::string const& name, bool throwIfNotExists = true );
*/
boolean Teuchos_ParameterList_remove ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  boolean throwIfNotExists );

/* Original C++ prototype:
   ParameterList& sublist( const std::string& name, bool mustAlreadyExist = false ,const std::string& docString = "" );
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  boolean mustAlreadyExist, const char docString[] );

/* Original C++ prototype:
   const ParameterList& sublist(const std::string& name) const;
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist_existing ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   const std::string& name() const;
*/
const char * Teuchos_ParameterList_name_it ( 
  CT_Teuchos_ParameterList_ID_t selfID );

/* Original C++ prototype:
   bool isParameter(const std::string& name) const;
*/
boolean Teuchos_ParameterList_isParameter ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   bool isSublist(const std::string& name) const;
*/
boolean Teuchos_ParameterList_isSublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   template<typename T> bool isType(const std::string& name) const;
*/
boolean Teuchos_ParameterList_isType_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   template<typename T> bool isType(const std::string& name) const;
*/
boolean Teuchos_ParameterList_isType_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/* Original C++ prototype:
   template<typename T> bool isType(const std::string& name, T* ptr) const;
*/
boolean Teuchos_ParameterList_isType_double_type ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double * ptr );

/* Original C++ prototype:
   template<typename T> bool isType(const std::string& name, T* ptr) const;
*/
boolean Teuchos_ParameterList_isType_int_type ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int * ptr );

/* Original C++ prototype:
   std::string currentParametersString() const;
*/
const char * Teuchos_ParameterList_currentParametersString ( 
  CT_Teuchos_ParameterList_ID_t selfID );

/* Original C++ prototype:
   void validateParameters( ParameterList const& validParamList, int const depth = 1000, EValidateUsed const validateUsed = VALIDATE_USED_ENABLED, EValidateDefaults const validateDefaults = VALIDATE_DEFAULTS_ENABLED ) const;
*/
void Teuchos_ParameterList_validateParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth, 
  const CT_EValidateUsed_E_t validateUsed, 
  const CT_EValidateDefaults_E_t validateDefaults );

/* Original C++ prototype:
   void validateParametersAndSetDefaults( ParameterList const& validParamList, int const depth = 1000 );
*/
void Teuchos_ParameterList_validateParametersAndSetDefaults ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CTEUCHOS_PARAMETERLIST_H */

