// //////////////////////////////////////////////////////////////////////////////////////
// COOMatrixClass.cpp

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include <sstream>

#include "../include/COOMatrixClass.h"
#include "../include/SparseCOOReadMatrix.h"
#include "LinAlgPack/include/IVector.h"
#include "LinAlgPack/include/Range1D.h"

// Junk, test compilation
//#include "ReferenceCountingPackDef.h"
//ReferenceCountingPack::RefCount<double> ref1;

//#include "SequentialAllocatorPack.h"
//template SequentialAllocatorPack::SequentialAllocator<double>;

// ///////////////////////////////////////////////////////////////////////////////////
// COOMatrix

SparseLinAlgPack::COOMatrix& SparseLinAlgPack::COOMatrix::operator=(const COOMatrix& coom)
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

void SparseLinAlgPack::COOMatrix::resize(size_type rows, size_type cols, size_type nz)
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

void SparseLinAlgPack::COOMatrix::initialize(std::istream& istrm) {
	// Read COO matrix into val, ivect, jvect
	read_coo_into_valarrays(istrm,rows_,cols_,nz_,val_,ivect_ref_.obj()
		,jvect_ref_.obj());
}
