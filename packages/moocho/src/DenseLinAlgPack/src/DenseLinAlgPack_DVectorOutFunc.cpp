// /////////////////////////////////////////////////////////////////////////////////
// DVectorOutFunc.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#include <ostream>
#include <iomanip>

#include "DenseLinAlgPack_DVectorOutFunc.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"

std::ostream& DenseLinAlgPack::output(std::ostream& os, const DVectorSlice& vs
	, LinAlgPackIO::fmtflags extra_flags)
{
	int w = os.width(0) - 1; // get the set width (minus 1 since a space is inserted)

	if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) )
		os << std::setw(0) << std::left << vs.dim() << std::endl << std::right;

	DVectorSlice::const_iterator itr = vs.begin();
	for( size_type i = 1; itr != vs.end(); ++i, ++itr ) {
		os << " " << std::setw(w) << (*itr) << ":" << i; // insert a space to be sure there is white space
		                                                 // inbetween adjacent elements.
	}

	if( !(extra_flags & LinAlgPackIO::no_insert_newlines_bit) )
		os << std::endl;

	return os;
}
