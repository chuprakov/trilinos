// //////////////////////////////////////////////////////////////////////////////////////
// COOMatrixClass.cpp
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

#include <sstream>

#include "AbstractLinAlgPack_COOMatrixClass.hpp"
#include "AbstractLinAlgPack_SparseCOOReadMatrix.hpp"
#include "DenseLinAlgPack_IVector.hpp"
#include "DenseLinAlgPack/src/Thyra_Range1D.hpp"

// Junk, test compilation
//#include "MemMngPackDef.h"
//MemMngPack::RefCount<double> ref1;

//#include "SequentialAllocatorPack.hpp"
//template SequentialAllocatorPack::SequentialAllocator<double>;

// ///////////////////////////////////////////////////////////////////////////////////
// COOMatrix

AbstractLinAlgPack::COOMatrix& AbstractLinAlgPack::COOMatrix::operator=(const COOMatrix& coom)
{
	if(this == &coom) return *this;	// assignment to self

	val_.resize(coom.nz_);	// must resize, this is why you can't use default assignment.

	// Now assign the members(this is what the default would have done).
	rows_				= coom.rows_;
	cols_				= coom.cols_;
	nz_					= coom.nz_;
	val_				= coom.val_;
	ivect_ref_			= coom.ivect_ref_;
	jvect_ref_			= coom.jvect_ref_;

	return *this;
}

void AbstractLinAlgPack::COOMatrix::resize(size_type rows, size_type cols, size_type nz)
{
	 // don't resize if you don't have to to presearve row and column access just in case.
	if(rows == rows_ && cols == cols_ && nz == nz_) return;
	
	rows_ = rows;
	cols_ = cols;
	nz_ = nz;
	val_.resize(nz);
	ivect_ref_.obj().resize(nz);
	jvect_ref_.obj().resize(nz);
}

void AbstractLinAlgPack::COOMatrix::initialize(std::istream& istrm) {
	// Read COO matrix into val, ivect, jvect
	read_coo_into_valarrays(istrm,rows_,cols_,nz_,val_,ivect_ref_.obj()
		,jvect_ref_.obj());
}
