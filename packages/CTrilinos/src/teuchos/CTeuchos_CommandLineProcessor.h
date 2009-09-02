#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CTEUCHOS_COMMANDLINEPROCESSOR_H
#define CTEUCHOS_COMMANDLINEPROCESSOR_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Teuchos_CommandLineProcessor_Cast() and Teuchos_CommandLineProcessor_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Teuchos_CommandLineProcessor table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

/* Original C++ prototype:
   ~CommandLineProcessor();
*/
void Teuchos_CommandLineProcessor_Destroy ( 
  CT_Teuchos_CommandLineProcessor_ID_t * selfID );

/* Original C++ prototype:
   CommandLineProcessor( bool throwExceptions = true ,bool recogniseAllOptions = true ,bool addOutputSetupOptions = false );
*/
CT_Teuchos_CommandLineProcessor_ID_t Teuchos_CommandLineProcessor_Create ( 
  boolean throwExceptions, boolean recogniseAllOptions, 
  boolean addOutputSetupOptions );

/* Original C++ prototype:
   void throwExceptions( const bool & throwExceptions );
*/
void Teuchos_CommandLineProcessor_throwExceptions_set ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  boolean throwExceptions );

/* Original C++ prototype:
   bool throwExceptions() const;
*/
boolean Teuchos_CommandLineProcessor_throwExceptions_get ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID );

/* Original C++ prototype:
   void recogniseAllOptions( const bool & recogniseAllOptions );
*/
void Teuchos_CommandLineProcessor_recogniseAllOptions_set ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  boolean recogniseAllOptions );

/* Original C++ prototype:
   bool recogniseAllOptions() const;
*/
boolean Teuchos_CommandLineProcessor_recogniseAllOptions_get ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID );

/* Original C++ prototype:
   void addOutputSetupOptions( const bool &addOutputSetupOptions );
*/
void Teuchos_CommandLineProcessor_addOutputSetupOptions_set ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  boolean addOutputSetupOptions );

/* Original C++ prototype:
   bool addOutputSetupOptions() const;
*/
boolean Teuchos_CommandLineProcessor_addOutputSetupOptions_get ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID );

/* Original C++ prototype:
   void setDocString( const char doc_string[] );
*/
void Teuchos_CommandLineProcessor_setDocString ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char doc_string[] );

/* Original C++ prototype:
   void setOption( const char option_true[] ,const char option_false[] ,bool *option_val ,const char documentation[] = NULL );
*/
void Teuchos_CommandLineProcessor_setOption_bool ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char option_true[], const char option_false[], 
  boolean * option_val, const char documentation[] );

/* Original C++ prototype:
   void setOption( const char option_name[] ,int *option_val ,const char documentation[] = NULL ,const bool required = false );
*/
void Teuchos_CommandLineProcessor_setOption_int ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char option_name[], int * option_val, 
  const char documentation[], const boolean required );

/* Original C++ prototype:
   void setOption( const char option_name[] ,double *option_val ,const char documentation[] = NULL ,const bool required = false );
*/
void Teuchos_CommandLineProcessor_setOption_double ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char option_name[], double * option_val, 
  const char documentation[], const boolean required );

/* Original C++ prototype:
   void setOption( const char option_name[] ,std::string *option_val ,const char documentation[] = NULL ,const bool required = false );
*/
void Teuchos_CommandLineProcessor_setOption_str ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char option_name[], char * option_val[], 
  const char documentation[], const boolean required );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CTEUCHOS_COMMANDLINEPROCESSOR_H */

