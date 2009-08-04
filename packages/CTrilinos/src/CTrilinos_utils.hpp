#ifndef CTRILINOS_UTILS_HPP
#define CTRILINOS_UTILS_HPP


#include <string>

#include "CTrilinos_config.h"
#include "CTrilinos_enums.h"
#include "CTrilinos_Table.hpp"


namespace CTrilinos {


/* convert string to const char * since C does not support former */
const char *str2cptr( std::string s );

/* convert const char * to string since C does not support latter */
std::string cptr2str( const char *pc );

/* stringify the enum name */
std::string enum2str( CTrilinos_Type_ID_t ty );

/* isSameObject(generic_id, generic_id) */
bool isSameObject( CTrilinos_Object_ID_t id1, CTrilinos_Object_ID_t id2 );


} // namespace CTrilinos


#endif

