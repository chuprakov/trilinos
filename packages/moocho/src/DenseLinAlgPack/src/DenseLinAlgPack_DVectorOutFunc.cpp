// /////////////////////////////////////////////////////////////////////////////////
// VectorOutFunc.cpp

#include <ostream>
#include <iomanip>

#include "../include/VectorOutFunc.h"
#include "../include/VectorClass.h"

std::ostream& LinAlgPack::output(std::ostream& os, const VectorSlice& vs
	, LinAlgPackIO::fmtflags extra_flags)
{
	int w = os.width(0) - 1; // get the set width (minus 1 since a space is inserted)

	if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) )
		os << std::setw(0) << std::left << vs.size() << std::endl << std::right;

	VectorSlice::const_iterator itr = vs.begin();
	for( size_type i = 1; itr != vs.end(); ++i, ++itr ) {
		os << " " << std::setw(w) << (*itr) << ":" << i; // insert a space to be sure there is white space
		                                                 // inbetween adjacent elements.
	}

	if( !(extra_flags & LinAlgPackIO::no_insert_newlines_bit) )
		os << std::endl;

	return os;
}
