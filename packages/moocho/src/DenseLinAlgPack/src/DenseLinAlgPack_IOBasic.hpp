// ///////////////////////////////////////////////////////////////////////////////////
// LinAlgPackIOBasic.h
//
// Basic declarations for I/O

#ifndef LINALGPACK_IO_BASIC_H
#define LINALGPACK_IO_BASIC_H

#include <stdexcept>

#include "LinAlgPackTypes.h"

namespace LinAlgPack {
namespace LinAlgPackIO {

/// Exception throw on input error
class InputException : public std::logic_error
{public: InputException(const std::string& what_arg) : std::logic_error(what_arg) {}};

///
typedef int fmtflags;	
	
/// Format flags
enum { ignore_dim_bit = 0x0001, no_insert_newlines_bit = 0x0002 };

}	// end namespace LinAlgPackIO
}	// end namespace LinAlgPack



#endif	// LINALGPACK_IO_BASIC_H