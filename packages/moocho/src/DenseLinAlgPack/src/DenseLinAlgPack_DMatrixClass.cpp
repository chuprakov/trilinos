// //////////////////////////////////////////////////////////////////////////////////
// GenMatrixClass.cpp
//

#include <iomanip>

#include "../include/GenMatrixClass.h"

// ////////////////////////////////////////////////////////////////////////////////
// GenMatrixSlice

LinAlgPack::VectorSlice LinAlgPack::GenMatrixSlice::p_diag(difference_type k) const {
	if(k > 0) {
		validate_col_subscript(k+1);
		// upper diagonal (k > 0)
		return VectorSlice( const_cast<value_type*>(col_ptr(1)) + k*max_rows()
			, cols()-k > rows() ? rows() : cols()-k, max_rows()+1 );
	}
	// lower diagonal (k < 0) or center diagonal (k = 0)
	validate_row_subscript(-k+1);
	return VectorSlice( const_cast<value_type*>(col_ptr(1)) - k
		, rows()+k > cols() ? cols() : rows()+k, max_rows()+1 );
}

LinAlgPack::EOverLap LinAlgPack::GenMatrixSlice::overlap(const GenMatrixSlice& gms) const
{
	typedef GenMatrixSlice::size_type size_type;
	
	const VectorSlice::value_type
		*raw_ptr1 = col_ptr(1),
		*raw_ptr2 = gms.col_ptr(1);

	if( !raw_ptr1 || !raw_ptr2 )
		return NO_OVERLAP;		// If one of the views is unbound there can't be any overlap

	VectorSlice::size_type
		max_rows1 = max_rows(),
		max_rows2 = gms.max_rows(),
		rows1 = rows(),
		rows2 = gms.rows(),
		cols1 = cols(),
		cols2 = gms.cols();

	// Establish a frame of reference where raw_ptr1 < raw_ptr2
	if(raw_ptr1 > raw_ptr2) {
		std::swap(raw_ptr1,raw_ptr2);
		std::swap(max_rows1,max_rows2);
		std::swap(rows1,rows2);
		std::swap(cols1,cols2);
	}

	if( raw_ptr1 + (cols1 - 1) * max_rows1 + (rows1 - 1) < raw_ptr2 ) {
		return NO_OVERLAP; // can't be any overlap
	}

	VectorSlice::size_type
		start1 = 0,
		start2 = raw_ptr2 - raw_ptr1;

	if(start1 == start2 && max_rows1 == max_rows2 && rows1 == rows2 && cols1 == cols2)
		return SAME_MEM;
	if(start1 + (rows1 - 1) + (cols1 - 1) * max_rows1 < start2)
		return NO_OVERLAP;	// start2 is past the last element in m1 so no overlap.
	// There may be some overlap so determine if start2 lays in the region of m1.
	// Determine the number of elements in v that start2 is past the start of a
	// column of m1.  If start2 was the first element in one of m1's cols
	// row_i would be 1, and if it was just before for first element of the next
	// column of m1 then row_i would equal to max_rows1.
	size_type row_i = (start2 - start1 + 1) % max_rows1;
	if(row_i <= rows1)
		return SOME_OVERLAP; // start2 is in a column of m1 so there is some overlap
	// Determine how many rows in the original matrix are below the last row in m1
	size_type lower_rows = max_rows1 - (start1 % max_rows1 + rows1);
	if(row_i < rows1 + lower_rows)
		return NO_OVERLAP; // m2 lays below m1 in the original matrix
	// If you get here start2 lays above m1 in original matix so if m2 has enough
	// rows then the lower rows of m2 will overlap the upper rows of m1.
	if(row_i + rows2 - 1 <= max_rows1)
		return NO_OVERLAP; // m2 completely lays above m1
	return SOME_OVERLAP; // Some lower rows of m2 extend into m1 
}

// /////////////////////////////////////////////////////////////////////////////////
// GenMatrix

LinAlgPack::VectorSlice LinAlgPack::GenMatrix::p_diag(difference_type k) const {	
	if(k > 0) {
		validate_col_subscript(k+1);
		// upper diagonal (k > 0)
		return VectorSlice( const_cast<value_type*>(col_ptr(1)) + k * rows()
			, cols()-k > rows() ? rows() : cols()-k, rows()+1 );
	}
	// lower diagonal (k < 0) or center diagonal (k = 0)
	validate_row_subscript(-k+1);
	return VectorSlice( const_cast<value_type*>(col_ptr(1)) - k
		, rows()+k > cols() ? cols() : rows()+k, rows()+1 );
}

LinAlgPack::EOverLap LinAlgPack::GenMatrix::overlap(const GenMatrixSlice& gms) const {
	return (*this)().overlap(gms);
}

// ///////////////////////////////////////////////////////////////////////////////
// Non-member funcitons

void LinAlgPack::assert_gms_sizes(const GenMatrixSlice& gms1, BLAS_Cpp::Transp trans1
	, const GenMatrixSlice& gms2, BLAS_Cpp::Transp trans2)
{
	if	(
			(trans1 == trans2) ? 
				gms1.rows() == gms2.rows() && gms1.cols() == gms2.cols() 
				: gms1.rows() == gms2.cols() && gms1.cols() == gms2.rows()
		)
		return; // compatible sizes so exit
	// not compatible sizes.
	throw std::length_error("Matrix sizes are not the compatible");
}