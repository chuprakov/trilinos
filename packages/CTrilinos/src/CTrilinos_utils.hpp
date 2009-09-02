#ifndef CTRILINOS_UTILS_HPP
#define CTRILINOS_UTILS_HPP


#include <string>

#include "CTrilinos_config.h"
#include "CTrilinos_enums.h"
#include "CTrilinos_Table.hpp"


namespace CTrilinos {


void pass_string_out( const std::string * const s, char *c[] );

void pass_string_in( const char * const c[], std::string *s );

/* stringify the enum name */
std::string enum2str( CTrilinos_Type_ID_t ty );

/* isSameObject(generic_id, generic_id) */
bool isSameObject( CTrilinos_Object_ID_t id1, CTrilinos_Object_ID_t id2 );


} // namespace CTrilinos


#endif

