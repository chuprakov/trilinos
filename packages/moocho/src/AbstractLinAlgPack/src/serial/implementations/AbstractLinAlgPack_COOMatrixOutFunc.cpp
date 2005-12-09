// //////////////////////////////////////////////////////////////////////////////
// COOMatrixOutFunc.cpp
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

#include "AbstractLinAlgPack_COOMatrixOutFunc.hpp"
#include "AbstractLinAlgPack_COOMatrixClass.hpp"

std::ostream& AbstractLinAlgPack::output(std::ostream& o, const COOMatrix& coom) {

	o	<< coom.rows() << " " << coom.cols() << " " << coom.nz() << "\n";

	const COOMatrix::value_type
		*itr_val		= coom.const_val(),
		*itr_val_end	= coom.const_val() + coom.nz();
	const COOMatrix::indice_type
		*itr_ivect	= coom.const_ivect(),
		*itr_jvect	= coom.const_jvect();

	for(; itr_val != itr_val_end; ++itr_val, ++itr_ivect, ++itr_jvect)
		o << *itr_val << ":" << *itr_ivect << ":" << *itr_jvect << " ";
	o << "\n";

	return o;
}
