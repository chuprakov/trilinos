// /////////////////////////////////////////////////////////////////////////////////
// DMatrixOutFunc.cpp
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

#include "DenseLinAlgPack_DMatrixOutFunc.hpp"
#include "DenseLinAlgPack_DVectorOutFunc.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"

std::ostream& DenseLinAlgPack::output(std::ostream& os, const DMatrixSlice& gms
	, LinAlgPackIO::fmtflags extra_flags )
{
	int w = os.width(0) - 1; // get the set width (minus 1 since a space is inserted)

	if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
		std::ios_base::fmtflags old = os.flags();
		os	<< std::setw(0) << std::left << gms.rows() << ' ' << gms.cols()
			<< std::endl;
		os.flags(old);
	}

	if( gms.rows() && gms.cols() ) {
		for(size_type i = 1; i <= gms.rows();++i) {
			const DVectorSlice& vs =gms.row(i); 
			DVectorSlice::const_iterator itr = vs.begin();
			for( size_type j = 1; itr != vs.end(); ++j, ++itr ) {
				os << " " << std::setw(w) << (*itr) << ":" << i << ":" << j;
			}
			os << std::endl;
		}
	}
	
	return os;
}
