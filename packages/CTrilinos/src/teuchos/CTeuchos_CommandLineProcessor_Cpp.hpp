#include "CTrilinos_config.h"

#ifndef CTEUCHOS_COMMANDLINEPROCESSOR_CPP_HPP
#define CTEUCHOS_COMMANDLINEPROCESSOR_CPP_HPP


#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"


namespace CTeuchos {


using Teuchos::RCP;


/* get Teuchos::CommandLineProcessor from non-const table using CT_Teuchos_CommandLineProcessor_ID */
const RCP<Teuchos::CommandLineProcessor>
getCommandLineProcessor( CT_Teuchos_CommandLineProcessor_ID_t id );

/* get Teuchos::CommandLineProcessor from non-const table using CTrilinos_Object_ID_t */
const RCP<Teuchos::CommandLineProcessor>
getCommandLineProcessor( CTrilinos_Object_ID_t id );

/* get const Teuchos::CommandLineProcessor from either the const or non-const table
 * using CT_Teuchos_CommandLineProcessor_ID */
const RCP<const Teuchos::CommandLineProcessor>
getConstCommandLineProcessor( CT_Teuchos_CommandLineProcessor_ID_t id );

/* get const Teuchos::CommandLineProcessor from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const RCP<const Teuchos::CommandLineProcessor>
getConstCommandLineProcessor( CTrilinos_Object_ID_t id );

/* store Teuchos::CommandLineProcessor in non-const table */
CT_Teuchos_CommandLineProcessor_ID_t
storeCommandLineProcessor( Teuchos::CommandLineProcessor *pobj );

/* store const Teuchos::CommandLineProcessor in const table */
CT_Teuchos_CommandLineProcessor_ID_t
storeConstCommandLineProcessor( const Teuchos::CommandLineProcessor *pobj );

/* dump contents of Teuchos::CommandLineProcessor and const Teuchos::CommandLineProcessor tables */
void
purgeCommandLineProcessorTables(  );


} // namespace CTeuchos


#endif // CTEUCHOS_COMMANDLINEPROCESSOR_CPP_HPP


