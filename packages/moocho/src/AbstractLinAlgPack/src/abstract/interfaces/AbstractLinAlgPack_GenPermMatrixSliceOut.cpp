// ///////////////////////////////////////////
// GenPermMatrixSliceOut.cpp
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

#include "AbstractLinAlgPack_GenPermMatrixSliceOut.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSlice.hpp"

std::ostream& AbstractLinAlgPack::operator<<(
	std::ostream& out, const GenPermMatrixSlice& P
	)
{
	out	<< P.rows() << " " << P.cols() << " " << P.nz() << std::endl;
	if( P.is_identity() ) {
		out << "identity matrix\n";
	}
	else if( P.nz() ) {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr )
			out << " " << itr->row_i() << ":" << itr->col_j();
		out << std::endl;
	}
	return out;
}
