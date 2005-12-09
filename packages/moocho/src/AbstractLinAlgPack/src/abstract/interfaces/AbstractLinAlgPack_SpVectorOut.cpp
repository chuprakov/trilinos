// //////////////////////////////////////////////////////////////////////////////////////
// SpVectorOut.cpp
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

#include <iomanip>
#include <ostream>

#include "AbstractLinAlgPack_SpVectorOut.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"

std::ostream&
AbstractLinAlgPack::operator<<(std::ostream& os, const SpVectorSlice& svs)
{
	os << std::left << std::setw(0) << svs.dim() << "  " << svs.nz() << std::endl << std::right;
	if( !svs.dim() ) return os;
	const SpVectorSlice::difference_type offset = svs.offset();
	for(SpVectorSlice::const_iterator itr = svs.begin(); itr != svs.end(); ++itr )
		os << "  " << itr->value() << ":" << (itr->index() + offset);
	return os << std::endl;
}
