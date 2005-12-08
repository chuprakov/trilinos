// /////////////////////////////////////////////////////
// MatlabPack.cpp
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

#include <limits>
#include <iomanip>
#include <ostream>

#include "DenseLinAlgPack_MatlabPack.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"

namespace DenseLinAlgPack {

std::ostream& MatlabPack::out( std::ostream& o, const char* name, const DVectorSlice& vs
	, BLAS_Cpp::Transp trans )
{
	int p = o.precision();
	o.precision( std::numeric_limits<value_type>::digits10 + 3 );
	try {
		o << name << " =  [ ";
		for( DVectorSlice::const_iterator itr = vs.begin(); itr != vs.end(); ++itr )
			o << ' ' << *itr << ";";
		o << "];" << (trans == BLAS_Cpp::no_trans ? ' ' : '\'' ) << std::endl;
	}
	catch(...) {
		o.precision(p);
		throw;
	}
	o.precision(p);
	return o;
}

std::ostream& MatlabPack::out( std::ostream& o, const char* name, const DMatrixSlice& gms
	, BLAS_Cpp::Transp trans )
{
	int p = o.precision();
	o.precision( std::numeric_limits<value_type>::digits10 + 3 );
	try {
		o << name << " =  [\n";
		for( size_type i = 1; i <= gms.rows(); ++i ) {
			const DVectorSlice vs = gms.row(i);
			for( DVectorSlice::const_iterator itr = vs.begin(); itr != vs.end(); ++itr )
				o << *itr << ", ";
			o << ";\n";
		}
		o << "];" << (trans == BLAS_Cpp::no_trans ? ' ' : '\'' ) << std::endl;
	}
	catch(...) {
		o.precision(p);
		throw;
	}
	o.precision(p);
	return o;
}

}	// end namespace DenseLinAlgPack
