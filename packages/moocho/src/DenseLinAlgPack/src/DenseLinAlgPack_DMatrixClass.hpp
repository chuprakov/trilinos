// ///////////////////////////////////////////////////////////////////////////////////
// GenMatrixClass.h
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
//
// General matrix and matrix region (slice) classes

#ifndef GEN_MATRIX_CLASS_H
#define GEN_MATRIX_CLASS_H

#include "VectorClass.h"
#include "GenMatrixAssign.h"

/** @name {\bf Dense 2-D Rectangular Matrix Absractions}.
  *
  * The class GenMatrix is a storage class for 2-D matrices while the class GenMatrixSlice
  * is used to represent rectangular regions of a GenMatrix object.
  */
//@{
//		begin General Rectangular 2-D Matrices scope

namespace LinAlgPack {

class GenMatrix;

/** @name {\bf General Matrix Classes}. */
//@{

// ////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenMatrixClass
//

///
/** 2-D General Rectangular Matrix Subregion Class (Slice) (column major).
  *
  * This class is used to represent a rectangular matrix.  It uses a BLAS-like
  * slice of a raw C++ array.  Objects of this class can represent
  * an entire matrix or any rectangular subregion of a matrix.
  */
class GenMatrixSlice {
public:
	/** @name {\bf Nested Member Types (STL)}.
	  *
	  * These nested types give the types used in the interface to the class.
	  *
	  * \begin{description}
	  *	<li>[#value_type#]				- type being stored in the underlying valarray<>			
	  *	<li>[#size_type#]				- type for the number rows and coluns
	  *	<li>[#difference_type#]		- type for the stride between elements (in a row or column)
	  *	<li>[#reference#]				- value_type&
	  *	<li>[#const_reference#]		- const value_type&
	  * \end{description}
	  */
	//@{
	//@}
	typedef LinAlgPack::value_type					value_type;
	typedef LinAlgPack::size_type					size_type;
	typedef ptrdiff_t								difference_type;
	typedef value_type&								reference;
	typedef const value_type&						const_reference;

	/** @name {\bf Constructors}.
	  *
	  * These constructors are used by the other entities in the library
	  * to create GenMatrixSlices.  In general user need not use these
	  * constructors directly.  Instead, the general user should use the 
	  * subscript operators to create subregions of GenMatrix and GenMatrixSlice
	  * objects.
	  * 
	  * The default copy constructor is used and is therefore not shown here.
	  */

	//@{

	/// 
	/** Construct an empty view.
	  *
	  * The client can then call bind(...) to bind the view.
	  */
	GenMatrixSlice();

	/// 
	/** Construct a veiw of a matrix from a raw C++ array.
	  *
	  * The GenMatrixSlice constructed represents a 2-D matrix whos elements are stored
	  * in the array starting at #ptr#.  This is how the BLAS represent general rectangular
	  * matrices.
	  * The class can be used to provide a non-constant view the elements (#GenMatrix#)
	  * or a constant view (#const GenMatrixSlice#).  Here is an example of how to
	  * create a constant view:
	  *
	  \verbatim
		const GenMatrixSlice::size_type m = 4, n = 4;
		const GenMatrixSlice::value_type ptr[m*n] = { ... };
		const GenMatrixslice mat(cosnt_cast<GenMatrixSlice::value_type*>(ptr),m*n,m,m,n);
	  \endverbatim
	  *
	  * The #const_cast<...># such as in the example above is perfectly safe to use
	  * when constructing  #const# veiw of #const# data.  On the other hand casting
	  * away #const# and then non-#const# use is not safe in general since the original
	  * #const# data may reside in ROM for example.  By using a non-#const# pointer in the
	  * constructor you avoid accidentally making a non-#const# view of #const# data.
	  *
	  * Preconditions: <ul>
	  *		<li> #rows <= max_rows# (throw out_of_range)
	  *		<li> #start + rows + max_rows * (cols - 1) <= v.size()# (throw out_of_range)
	  *		</ul>
	  *
	  * @param	ptr			pointer to the storage of the elements of the matrix (column oriented).
	  * @param	size		total size of the storage pointed to by #ptr# (for size checking)
	  * @param	max_rows	number of rows in the full matrix this sub-matrix was taken from.
	  *						This is equivalent to the leading dimension argument (LDA) in
	  *						the BLAS.
	  * @param	rows		number of rows this matrix represents
	  *	@param	cols		number of columns this matrix represents
	  */
	GenMatrixSlice( value_type* ptr, size_type size
		, size_type max_rows, size_type rows, size_type cols );

	/// 
	/** Construct a submatrix region of another GenMatrixSlice.
	  *
	  * This constructor simplifies the creation of subregions using the subscript
	  * operators.
	  *
	  * I and J must be bounded ranges (full_range() == false).
	  *
	  * Preconditions: <ul>
	  *		<li> #I.full_range() == false# (throw out_of_range)
	  *		<li> #J.full_range() == false# (throw out_of_range)
	  *		<li> #I.ubound() <= gms.rows()# (throw out_of_range)
	  *		<li> #J.ubound() <= gms.cols()# (throw out_of_range)
	  *		</ul>
	  */
	GenMatrixSlice( GenMatrixSlice& gms, const Range1D& I
		, const Range1D& J );
	//@}

	///
	/** Set to the view of the input GenMatrixSlice.
	  *
	  *
	  */
	void bind( GenMatrixSlice gms );

	/** @name {\bf Dimensionality, Misc}. */
	//@{

	/// Return the number of rows
	size_type		rows() const;
	/// Return the number of columns
	size_type		cols() const;

	/// 
	/** Returns the degree of memory overlap of #this# and the GenMatrixSlice object #gms#.
	  *
	  * @return 
	  *		\begin{description}
	  *		<li>[NO_OVERLAP]	There is no memory overlap between #this# and #gms#.
	  *		<li>[SOME_OVERLAP]	There is some memory locations that #this# and #gms# share.
	  *		<li>[SAME_MEM]		The GenMatrixSlice objects #this# and #gms# share the exact same memory locations.
	  *		\end{description}
	  */
	EOverLap overlap(const GenMatrixSlice& gms) const;

	//@}

	/** @name {\bf Individual Element Access Subscripting (lvalue)}. */
	//@{

	/// Return element at row i, col j (i,j) (1-based) (throws std::out_of_range if i, j are out of bounds)
	reference		operator()(size_type i, size_type j);
	/// Return element at row i, col j (i,j) (1-based) (throws std::out_of_range if i, j are out of bounds)
	const_reference	operator()(size_type i, size_type j) const;

	//@}

	/** @name {\bf Subregion Access (1-based)}.
	  *
	  * These member functions allow access to subregions of the GenMatrixSlice object.
	  * The functions, row(i), col(j), and diag(k) return VectorSlice objects while
	  * the subscripting operators opeator()(I,J) return GenMatrixSlice objects for
	  * rectangular subregions.
	  */
	//@{

	/// Return VectorSlice object representing the ith row (1-based; 1,2,..,#this->rows()#, or throw std::out_of_range)
	VectorSlice			row(size_type i);
	/// Same as above
	const VectorSlice	row(size_type i) const;
	/// Return VectorSlice object representing the jth column (1-based; 1,2,..,#this->cols()#, or throw std::out_of_range)
	VectorSlice			col(size_type j);
	/// Same as above
	const VectorSlice	col(size_type j) const;
	/// 
	/** Return VectorSlice object representing a diagonal.
	  *
	  * Passing k == 0 returns the center diagonal.  Values of k < 0 are the lower diagonals
	  * (k = -1, -2, ..., -#this->rows()# + 1).  Values of k > 0 are the upper diagonals
	  * (k = 1, 2, ..., #this->cols()# - 1).
	  *
	  * Preconditions: <ul>
	  *		<li> #[k < 0] k <= this->rows() + 1# (throw out_of_range)
	  *		<li> #[k > 0] k <= this->cols() + 1# (throw out_of_range)
	  *		</ul>
	  */
	VectorSlice			diag(difference_type k = 0);
	/// Same as above.
	const VectorSlice	diag(difference_type k = 0) const;
	/// 
	/** Extract a rectangular subregion containing rows I, and columns J.
	  *
	  * This operator function returns a GenMatrixSlice that represents a
	  * rectangular region of this GenMatrixSlice.  This submatrix region
	  * represents the rows I.lbound() to I.ubound() and columns J.lbound()
	  * to J.lbound().  If I or J is unbounded (full_range() == true, constructed
	  * with Range1D()), the all of the rows or columns respectively will be
	  * selected.  For example. To select all the rows and the first 5 columns of
	  * a matrix #A# you would use #A(Range1D(),Range1D(1,5))#.
	  *
	  * Preconditions: <ul>
	  *		<li> #[I.full_range() == false] I.ubound() <= this->rows()# (throw out_of_range)
	  *		<li> #[J.full_range() == false] J.ubound() <= this->cols()# (throw out_of_range)
	  *		</ul>
	  */
	GenMatrixSlice operator()(const Range1D& I, const Range1D& J);
	/// Same as above.
	const GenMatrixSlice operator()(const Range1D& I, const Range1D& J) const;
	/// 
	/** Extract a rectangular subregion containing rows i1 to i2, and columns j1 to j2.
	  *
	  * This operator function returns a GenMatrixSlice that represents a
	  * rectangular region of this GenMatrixSlice.  This submatrix region
	  * represents the rows i1 to 12 and colunms j1 to j2.
	  * 
	  * Preconditions: <ul>
	  *		<li> #i1 <= i2# (throw out_of_range)
	  *		<li> #i2 <= this->rows()# (throw out_of_range)
	  *		<li> #j1 <= j2# (throw out_of_range)
	  *		<li> #j2 <= this->cols()# (throw out_of_range)
	  *		</ul>
	  */
	GenMatrixSlice operator()(size_type i1, size_type i2, size_type j1
		, size_type j2);
	/// Same as above.
	const GenMatrixSlice operator()(size_type i1, size_type i2, size_type j1
		, size_type j2) const;
	/// Allow the address to be taken of an rvalue of this object.
	GenMatrixSlice* operator&() {
	  return this;
	}
	///
	const GenMatrixSlice* operator&() const {
	  return this;
	}
	/// Return reference of this.  Included for iniformity with GenMatrix
	GenMatrixSlice& operator()();
	/// Same as above
	const GenMatrixSlice& operator()() const;

	//@}

	/** @name {\bf Assignment operators}. */
	//@{

	///
	/** Sets all elements = alpha
	  *
	  * If the underlying valarray is unsized (#this->v().size() == 0#) the matrix is sized to 1 x 1
	  * and the single element is set to alpha.
	  *
	  * Postcondtions: <ul>
	  *		<li> #this->operator()(i,j) == alpha#, i = 1,2,...,#this->rows()#, j = 1,2,...,#this->cols()#
	  *		</ul>
	  */
	GenMatrixSlice& operator=(value_type alpha);
	///
	/**  Copies all of the elements of the GenMatrixSlice, #rhs#, into the elements of #this#.
	  *
	  * If the underlying valarray is unsized (#this->v().size() == 0#) the matrix is sized to
	  * the size of the rhs matrix.
	  *
	  * Precondtions: <ul>
	  *		<li> #this->rows() == gms_rhs.rows()# (throw length_error)
	  *		<li> #this->cols() == gms_rhs.cols()# (throw length_error)
	  *		</ul>
	  *
	  * Postcondtions: <ul>
	  *		<li> #this->operator()(i,j) == gms_rhs(i,j)#, i = 1,2,...,#this->rows()#, j = 1,2,...,#this->cols()#
	  *		</ul>
	  */
	GenMatrixSlice& operator=(const GenMatrixSlice& gms_rhs);

	//@}

	/** @name {\bf Raw data access}.
	  */
	//@{

	/// Return the number of rows in the full matrix. Equivalent to BLAS LDA argument.
	size_type			max_rows() const;
	///
	/** Return pointer to the first element in the underlying array the jth
	   * col (j is 1-based here [1,cols]).  If unsized col_ptr(1) returns zero if unsized.
	   */
	value_type*			col_ptr(size_type j);
	/// Same as above.
	const value_type*	col_ptr(size_type j) const;

	//@}

private:
	value_type		*ptr_;		// contains the data for the matrix region
	size_type		max_rows_,	// the number of rows in the full matrix
					rows_,		// the number of rows in this matrix region
					cols_;		// the number of cols in this matrix region

	// Assert the row subscript is in bounds (1-based), (throw std::out_of_range)
	void validate_row_subscript(size_type i) const;
	// Assert the column subscript is in bounds (1-based), (throw std::out_of_range)
	void validate_col_subscript(size_type j) const;
	// Assert that a constructed GenMatrixSlice has a valid range, (throw std::out_of_range)
	void validate_setup(size_type size) const;
	
	// Get a diagonal
	VectorSlice p_diag(difference_type k) const;

};	// end class GenMatrixSlice

///
/** 2-D General Rectangular Matrix (column major) Storage Class.
  *
  * This class provides the storage for 2-D rectangular matrices.
  */
class GenMatrix {
public:
	/** @name {\bf Nested Member Types (STL)}.
	  *
	  * These nested types give the types used in the interface to the class.
	  *
	  * \begin{description}
	  *	<li>[#value_type#]				type being stored in the underlying valarray<>			
	  *	<li>[#size_type#]				type for the number rows and coluns
	  *	<li>[#difference_type#]		type for the stride between elements (in a row or column)
	  *	<li>[#reference#]				value_type&
	  *	<li>[#const_reference#]		const value_type&
	  * \end{description}
	  */
	//@{
	//@}
	
	typedef LinAlgPack::value_type					value_type;
	typedef LinAlgPack::size_type					size_type;
	typedef ptrdiff_t								difference_type;
	typedef value_type&								reference;
	typedef const value_type&						const_reference;
	typedef std::valarray<value_type>				valarray;

	/** @name {\bf Constructors}.
	  *
	  * The general user uses these constructors to create a matrix.  
	  *
	  * The default constructor is used and is therefore not shown here.
	  */
	//@{

	/// Construct a matrix with rows = cols = 0
	GenMatrix();
	/// Construct an uninitialied rectangular matrix (rows x cols) 
	explicit GenMatrix(size_type rows, size_type cols);
	///
	/** Construct rectangular matrix (rows x cols) with elements initialized to val.
	  *
	  * Postconditions: <ul>
	  *		<li> #this->operator()(i,j) == val#, i = 1,2,...,#rows#, j = 1,2,...,#cols#  
	  *		</ul>
	  */
	explicit GenMatrix(value_type val, size_type rows, size_type cols);
	/// 
	/** Construct rectangular matrix (rows x cols) initialized to elements of p (by column).
	  *
	  * Postconditions: <ul>
	  *		<li> #this->operator()(i,j) == p[i-1 + rows * (j - 1)]#, i = 1,2,...,#rows#, j = 1,2,...,#cols#  
	  *		</ul>
	  */
	explicit GenMatrix(const value_type* p, size_type rows, size_type cols);
	///
	/** Construct a matrix from the elements in another GenMatrixSlice, #gms#.
	  *
	  * Postconditions: <ul>
	  *		<li> #this->operator()(i,j) == gms(i,j)#, i = 1,2,...,#rows#, j = 1,2,...,#cols#  
	  *		</ul>
	  */
	GenMatrix(const GenMatrixSlice& gms);

	//@}

	/** @name {\bf Memory Management, Dimensionality, Misc}. */
	//@{
	
	/// Resize matrix to a (rows x cols) matrix and initializes any added elements by val
	void resize(size_type rows, size_type cols, value_type val = value_type());

	/// frees memory and leaves a (0 x 0) matrix
	void free();
	
	/// Return the number of rows
	size_type		rows() const;

	/// Return the number of columns
	size_type		cols() const;

	/// 
	/** Returns the degree of memory overlap of #this# and the GenMatrixSlice object #gms#.
	  *
	  * @return 
	  *		\begin{description}
	  *		<li>[NO_OVERLAP]	There is no memory overlap between #this# and #gms#.
	  *		<li>[SOME_OVERLAP]	There is some memory locations that #this# and #gms# share.
	  *		<li>[SAME_MEM]		The GenMatrixSlice objects #this# and #gms# share the exact same memory locations.
	  *		\end{description}
	  */
	EOverLap overlap(const GenMatrixSlice& gms) const;

	//@}

	/** @name {\bf Individual Element Access Subscripting (lvalue)}. */
	//@{

	/// Return element at row i, col j (i,j) (1-based)
	reference		operator()(size_type i, size_type j);

	/// Return element at row i, col j (i,j) (1-based)
	const_reference	operator()(size_type i, size_type j) const;

	//@}

	/** @name {\bf Subregion Access (1-based)}.
	  *
	  * These member functions allow access to subregions of the GenMatrix object.
	  * The functions, row(i), col(j), and diag(k) return VectorSlice objects while
	  * the subscripting operators opeator()(I,J) return GenMatrixSlice objects for
	  * rectangular subregions.
	  */
	//@{

	/// Return VectorSlice object representing the ith row (1-based; 1,2,..,#this->rows()#)
	VectorSlice			row(size_type i);

	///
	const VectorSlice	row(size_type i) const;

	/// Return VectorSlice object representing the jth column (1-based; 1,2,..,#this->cols()#)
	VectorSlice			col(size_type j);

	///
	const VectorSlice	col(size_type j) const;

	/// 
	/** Return VectorSlice object representing a diagonal.
	  *
	  * Passing k == 0 returns the center diagonal.  Values of k < 0 are the lower diagonals
	  * (k = -1, -2, ..., #this->rows()# - 1).  Values of k > 0 are the upper diagonals
	  * (k = 1, 2, ..., #this->cols()# - 1).
	  *
	  * Preconditions: <ul>
	  *		<li> #[k < 0] k <= this->rows() + 1# (throw out_of_range)
	  *		<li> #[k > 0] k <= this->cols() + 1# (throw out_of_range)
	  *		</ul>
	  */
	VectorSlice			diag(difference_type k = 0);

	///
	const VectorSlice	diag(difference_type k = 0) const;

	/// 
	/** Extract a rectangular subregion containing rows I, and columns J.
	  *
	  * This operator function returns a GenMatrixSlice that represents a
	  * rectangular region of this GenMatrixSlice.  This submatrix region
	  * represents the rows I.lbound() to I.ubound() and columns J.lbound()
	  * to J.lbound().  If I or J is unbounded (full_range() == true, constructed
	  * with Range1D()), the all of the rows or columns respectively will be
	  * selected.  For example. To select all the rows and the first 5 columns of
	  * a matrix #A# you would use #A(Range1D(),Range1D(1,5))#.
	  *
	  * Preconditions: <ul>
	  *		<li> #[I.full_range() == false] I.ubound() <= this->rows()# (throw out_of_range)
	  *		<li> #[J.full_range() == false] J.ubound() <= this->cols()# (throw out_of_range)
	  *		</ul>
	  */
	GenMatrixSlice operator()(const Range1D& I, const Range1D& J);

	///
	const GenMatrixSlice operator()(const Range1D& I, const Range1D& J) const;

	/// 
	/** Extract a rectangular subregion containing rows i1 to i2, and columns j1 to j2.
	  *
	  * This operator function returns a GenMatrixSlice that represents a
	  * rectangular region of this GenMatrixSlice.  This submatrix region
	  * represents the rows i1 to 12 and colunms j1 to j2.
	  * 
	  * Preconditions: <ul>
	  *		<li> #i1 <= i2# (throw out_of_range)
	  *		<li> #i2 <= this->rows()# (throw out_of_range)
	  *		<li> #j1 <= j2# (throw out_of_range)
	  *		<li> #j2 <= this->cols()# (throw out_of_range)
	  *		</ul>
	  */
	GenMatrixSlice operator()(size_type i1, size_type i2, size_type j1
		, size_type j2);

	///
	const GenMatrixSlice operator()(size_type i1, size_type i2, size_type j1
		, size_type j2) const;

	/// Return a GenMatrixSlice that represents this entire matrix.
	GenMatrixSlice operator()();

	///
	const GenMatrixSlice operator()() const;

	//@}

	/** @name {\bf Implicit conversion operators}.
	  *
	  * These functions allow for the implicit converstion from a GenMatrix to a GenMatrixSlice.
	  * This implicit converstion is important for the proper usage of much of the
	  * libraries functionality.
	  */
	//@{

	///
	operator GenMatrixSlice();
	///
	operator const GenMatrixSlice() const;

	//@}

	/** @name {\bf Assignment Operators}. */
	//@{
	
	///
	/** Sets all elements = alpha
	  *
	  * If the underlying valarray is unsized (#this->v().size() == 0#) the matrix is sized to 1 x 1
	  * and the single element is set to alpha.
	  *
	  * Postcondtions: <ul>
	  *		<li> #this->operator()(i,j) == alpha#, i = 1,2,...,#this->rows()#, j = 1,2,...,#this->cols()#
	  */
	GenMatrix& operator=(value_type rhs);
	///
	/** Copies all of the elements of the GenMatrixSlice, #rhs#, into the elements of #this#.
	  *
	  * If #this# is not the same size as gms_rhs the #this# is resized.
	  *
	  * Postcondtions: <ul>
	  *		<li> #this->operator()(i,j) == gms_rhs(i,j)#, i = 1,2,...,#this->rows()#, j = 1,2,...,#this->cols()#
	  */
	GenMatrix& operator=(const GenMatrixSlice& gms_rhs);
	/// Same as above.  Needed to override the default assignment operator.
	GenMatrix& operator=(const GenMatrix& rhs);

	//@}

	/** @name {\bf Raw data access}.
	  */
	//@{

	/// Return the number of rows in the full matrix. Equivalent to BLAS LDA argument.
	size_type			max_rows() const;
	///
	/** Return pointer to the first element in the underlying array the jth
	   * col (j is 1-based here [1,cols]).  If unsized col_ptr(1) returns zero if unsized.
	   */
	value_type*			col_ptr(size_type j);
	/// Same as above.
	const value_type*	col_ptr(size_type j) const;

	//@}

private:
	std::valarray<value_type>	v_;
	size_type					rows_;

	// Assert the row subscript is in bounds (1-based), (throw std::out_of_range)
	void validate_row_subscript(size_type i) const;
	// Assert the column subscript is in bounds (1-based), (throw std::out_of_range)
	void validate_col_subscript(size_type j) const;

	// Get a diagonal, (throw std::out_of_range)
	VectorSlice p_diag(difference_type k) const;

};	// end class GenMatrix

//		end General Matix Classes scope
//@}

// ///////////////////////////////////////////////////////////////////////////////
// Non-member function declarations												//
// ///////////////////////////////////////////////////////////////////////////////

/** @name {\bf GenMatrix / GenMatrixSlice Associated Non-Member Functions}. */
//@{
//		begin non-member functions scope

inline 
///  
/** Explicit conversion function from GenMatrix to GenMatrixSlice.
  *
  * This is needed to allow a defered evaluation class (TCOL) to be evaluated using its
  * implicit conversion operator temp_type() (which returns GenMatrix for GenMatrixSlice
  * resulting expressions).
  */
GenMatrixSlice EvaluateToGenMatrixSlice(const GenMatrix& gm)
{	return GenMatrixSlice(gm); }

/// Assert two matrices are the same size and throws length_error if they are not (LINALGPACK_CHECK_RHS_SIZES).
void assert_gms_sizes(const GenMatrixSlice& gms1, BLAS_Cpp::Transp trans1, const GenMatrixSlice& gms2
	, BLAS_Cpp::Transp trans2);

inline 
/// Assert a matrix is square and throws length_error if it is not (LINALGPACK_CHECK_SLICE_SETUP).
void assert_gms_square(const GenMatrixSlice& gms) {
#ifdef LINALGPACK_CHECK_SLICE_SETUP
	if(gms.rows() != gms.cols())
		throw std::length_error("Matrix must be square");
#endif
} 

inline 
///
/** Utility to check if a lhs matrix slice is the same size as a rhs matrix slice.
  *
  * A GenMatrixSlice can not be resized since the rows_ property of the
  * GenMatrix it came from will not be updated.  Allowing a GenMatrixSlice
  * to resize from unsized would require that the GenMatrixSlice carry
  * a reference to the GenMatrix it was created from.  If this is needed
  * then it will be added.
  */
void assert_gms_lhs(const GenMatrixSlice& gms_lhs, size_type rows, size_type cols
	, BLAS_Cpp::Transp trans_rhs = BLAS_Cpp::no_trans)
{
	if(trans_rhs == BLAS_Cpp::trans) std::swap(rows,cols);
	if(gms_lhs.rows() == rows && gms_lhs.cols() == cols) return; // same size
	// not the same size so is an error
	throw std::length_error("assert_gms_lhs(...):  lhs GenMatrixSlice dim does not match rhs dim");
}

/** @name Return rows or columns from a possiblly transposed GenMatrix or GenMatrixSlice. */
//@{

inline 
///
VectorSlice row(GenMatrixSlice& gms, BLAS_Cpp::Transp trans, size_type i) {
	return (trans ==  BLAS_Cpp::no_trans) ? gms.row(i) : gms.col(i);
} 

inline 
///
VectorSlice col(GenMatrixSlice& gms, BLAS_Cpp::Transp trans, size_type j) {
	return (trans ==  BLAS_Cpp::no_trans) ? gms.col(j) : gms.row(j);
} 

inline 
///
const VectorSlice row(const GenMatrixSlice& gms, BLAS_Cpp::Transp trans, size_type i) {
	return (trans ==  BLAS_Cpp::no_trans) ? gms.row(i) : gms.col(i);
} 

inline 
///
const VectorSlice col(const GenMatrixSlice& gms, BLAS_Cpp::Transp trans, size_type j) {
	return (trans ==  BLAS_Cpp::no_trans) ? gms.col(j) : gms.row(j);
} 

inline 
///
VectorSlice row(GenMatrix& gm, BLAS_Cpp::Transp trans, size_type i) {
	return (trans ==  BLAS_Cpp::no_trans) ? gm.row(i) : gm.col(i);
} 

inline 
///
VectorSlice col(GenMatrix& gm, BLAS_Cpp::Transp trans, size_type j) {
	return (trans ==  BLAS_Cpp::no_trans) ? gm.col(j) : gm.row(j);
} 

inline 
///
const VectorSlice row(const GenMatrix& gm, BLAS_Cpp::Transp trans, size_type i) {
	return (trans ==  BLAS_Cpp::no_trans) ? gm.row(i) : gm.col(i);
} 

inline 
///
const VectorSlice col(const GenMatrix& gm, BLAS_Cpp::Transp trans, size_type j) {
	return (trans ==  BLAS_Cpp::no_trans) ? gm.col(j) : gm.row(j);
} 

//@}

inline 
/// Utility to resize a GenMatrix to the size of a rhs matrix.
void resize_gm_lhs(GenMatrix* gm_rhs, size_type rows, size_type cols
	, BLAS_Cpp::Transp trans_rhs)
{
	if(trans_rhs == BLAS_Cpp::trans) std::swap(rows,cols);
	gm_rhs->resize(rows,cols);
}

//		end non-member functions scope
//@}

//		end General Rectangular 2-D Matrices scope
//@}

// ////////////////////////////////////////////////////////////////////////////////
// Inline definitions of computationally independent member function definitions //
// ////////////////////////////////////////////////////////////////////////////////

// /////////////////////////////////////////////////////////////////////////////
// GenMatrixSlice inline member function definitions

// Private utilities

#ifndef LINALGPACK_CHECK_RANGE
inline
void GenMatrixSlice::validate_row_subscript(size_type i) const
{}
#endif

#ifndef LINALGPACK_CHECK_RANGE
inline
void GenMatrixSlice::validate_col_subscript(size_type j) const
{}
#endif

#ifndef LINALGPACK_CHECK_SLICE_SETUP
inline
void GenMatrixSlice::validate_setup(size_type size) const
{}
#endif

// Constructors

inline
GenMatrixSlice::GenMatrixSlice()
	: ptr_(0), max_rows_(0), rows_(0), cols_(0)
{}

inline
GenMatrixSlice::GenMatrixSlice( value_type* ptr, size_type size
		, size_type max_rows, size_type rows, size_type cols )
	: ptr_(ptr), max_rows_(max_rows), rows_(rows), cols_(cols)
{	
	validate_setup(size);
}

inline
GenMatrixSlice::GenMatrixSlice( GenMatrixSlice& gms, const Range1D& I
		, const Range1D& J)
	: ptr_( gms.col_ptr(1) + (I.lbound() - 1) + (J.lbound() - 1) * gms.max_rows() )
	, max_rows_(gms.max_rows())
	, rows_(I.size())
	, cols_(J.size())
{	
	gms.validate_row_subscript(I.ubound());
	gms.validate_col_subscript(J.ubound());
}

inline
void GenMatrixSlice::bind(GenMatrixSlice gms) {
	ptr_		= gms.ptr_;
	max_rows_	= gms.max_rows_;
	rows_		= gms.rows_;
	cols_		= gms.cols_;
}

// Size / Dimensionality

inline
GenMatrixSlice::size_type GenMatrixSlice::rows() const {
	return rows_;
}

inline
GenMatrixSlice::size_type GenMatrixSlice::cols() const {
	return cols_;
}

// Misc

// Element access

inline
GenMatrixSlice::reference GenMatrixSlice::operator()(size_type i, size_type j)
{	
	validate_row_subscript(i);
	validate_col_subscript(j);
	return ptr_[(i-1) + (j-1) * max_rows_];
}

inline
GenMatrixSlice::const_reference	GenMatrixSlice::operator()(size_type i, size_type j) const
{
	validate_row_subscript(i);
	validate_col_subscript(j);
	return ptr_[(i-1) + (j-1) * max_rows_];
}

// Subregion access (validated by constructor for GenMatrixSlice)

inline
VectorSlice  GenMatrixSlice::row(size_type i) {
	validate_row_subscript(i);
	return VectorSlice( ptr_ + (i-1), cols(), max_rows() );
} 

inline
const VectorSlice GenMatrixSlice::row(size_type i) const {
	validate_row_subscript(i);
	return VectorSlice( const_cast<value_type*>(ptr_) + (i-1), cols(), max_rows() );
} 

inline
VectorSlice	GenMatrixSlice::col(size_type j) {
	validate_col_subscript(j);
	return VectorSlice( ptr_ + (j-1)*max_rows(), rows(), 1 );
} 

inline
const VectorSlice GenMatrixSlice::col(size_type j) const {
	validate_col_subscript(j);
	return VectorSlice( const_cast<value_type*>(ptr_) + (j-1)*max_rows(), rows(), 1 );
} 

inline
VectorSlice GenMatrixSlice::diag(difference_type k) {
	return p_diag(k);
}

inline
const VectorSlice GenMatrixSlice::diag(difference_type k) const {
	return p_diag(k);
}

inline
GenMatrixSlice GenMatrixSlice::operator()(const Range1D& I, const Range1D& J) {
	return GenMatrixSlice(*this, full_range(I, 1, rows()), full_range(J,1,cols()));
}

inline
const GenMatrixSlice GenMatrixSlice::operator()(const Range1D& I, const Range1D& J) const {
	return GenMatrixSlice( const_cast<GenMatrixSlice&>(*this)
		, full_range(I, 1, rows()), full_range(J,1,cols()) );
}

inline
GenMatrixSlice GenMatrixSlice::operator()(size_type i1, size_type i2, size_type j1
	, size_type j2)
{
	return GenMatrixSlice(*this, Range1D(i1,i2), Range1D(j1,j2));
}

inline
const GenMatrixSlice GenMatrixSlice::operator()(size_type i1, size_type i2, size_type j1
	, size_type j2) const
{
	return GenMatrixSlice( const_cast<GenMatrixSlice&>(*this), Range1D(i1,i2)
		, Range1D(j1,j2) );
}

inline
GenMatrixSlice& GenMatrixSlice::operator()() {
	return *this;
}

inline
const GenMatrixSlice& GenMatrixSlice::operator()() const {
	return *this;
}

// Assignment operators

inline
GenMatrixSlice& GenMatrixSlice::operator=(value_type alpha) {
	assign(this, alpha);
	return *this;
}

inline
GenMatrixSlice& GenMatrixSlice::operator=(const GenMatrixSlice& rhs) {
	assign(this, rhs, BLAS_Cpp::no_trans);
	return *this;
}

// Raw data access

inline
GenMatrixSlice::size_type GenMatrixSlice::max_rows() const
{	return max_rows_; }

inline
GenMatrixSlice::value_type* GenMatrixSlice::col_ptr(size_type j) {
	if( ptr_ )
		validate_col_subscript(j);
	return ptr_ + (j-1) * max_rows();	// will be 0 if not bound to a view.
}

inline
const GenMatrixSlice::value_type* GenMatrixSlice::col_ptr(size_type j) const {
	if( ptr_ )
		validate_col_subscript(j);
	return ptr_ + (j-1) * max_rows();	// will be 0 if not bound to a view.
}

// ////////////////////////////////////////////////////////////////////////////////////////
// GenMatrix inline member function definitions

// Private utilities

#ifndef LINALGPACK_CHECK_RANGE
inline
void GenMatrix::validate_row_subscript(size_type i) const
{}
#endif

#ifndef LINALGPACK_CHECK_RANGE
inline
void GenMatrix::validate_col_subscript(size_type j) const
{}
#endif

// constructors

inline
GenMatrix::GenMatrix() : v_(), rows_(0)
{}

inline
GenMatrix::GenMatrix(size_type rows, size_type cols)
	: v_(rows*cols), rows_(rows)
{}

inline
GenMatrix::GenMatrix(value_type val, size_type rows, size_type cols)
	: v_(val,rows*cols), rows_(rows)
{}

inline
GenMatrix::GenMatrix(const value_type* p, size_type rows, size_type cols)
	: v_(rows*cols), rows_(rows)
{
// 6/7/00: valarray<> in libstdc++-2.90.7 has a bug in v_(p,size) so we do not
// use it.  This is a hack until I can find the time to remove valarray all
// together.
	std::copy( p, p + rows*cols, &v_[0] );
}

inline
GenMatrix::GenMatrix(const GenMatrixSlice& gms)
	: v_(gms.rows() * gms.cols()), rows_(gms.rows())
{	
	assign(this, gms, BLAS_Cpp::no_trans);
}

// Memory management

inline
void GenMatrix::resize(size_type rows, size_type cols, value_type val)
{
	v_.resize(rows*cols,val);
	v_ = val;
	rows_ = rows;
}

inline
void GenMatrix::free() {
	v_.resize(0);
	rows_ = 0;
}

// Size / Dimensionality

inline
GenMatrix::size_type GenMatrix::rows() const {
	return rows_;
}

inline
GenMatrix::size_type GenMatrix::cols() const {
	return rows_ > 0 ? v_.size() / rows_ : 0;
}

// Element access

inline
GenMatrix::reference GenMatrix::operator()(size_type i, size_type j)
{	 
	validate_row_subscript(i); validate_col_subscript(j);
	return v_[(i-1) + (j-1) * rows_];
}

inline
GenMatrix::const_reference GenMatrix::operator()(size_type i, size_type j) const
{
	validate_row_subscript(i); validate_col_subscript(j);
	return (const_cast<std::valarray<value_type>&>(v_))[(i-1) + (j-1) * rows_];
}

// subregion access (range checked by constructors)

inline
VectorSlice GenMatrix::row(size_type i)
{
	validate_row_subscript(i);
	return VectorSlice( col_ptr(1) + (i-1), cols(), rows() );
} 

inline
const VectorSlice GenMatrix::row(size_type i) const
{
	validate_row_subscript(i);
	return VectorSlice( const_cast<value_type*>(col_ptr(1)) + (i-1), cols(), rows() );
} 

inline
VectorSlice	GenMatrix::col(size_type j)
{
	validate_col_subscript(j);
	return VectorSlice( col_ptr(1) + (j-1) * rows(), rows(), 1 );
} 

inline
const VectorSlice GenMatrix::col(size_type j) const
{
	validate_col_subscript(j);
	return VectorSlice( const_cast<value_type*>(col_ptr(1)) + (j-1) * rows(), rows(), 1 ) ;
} 

inline
VectorSlice GenMatrix::diag(difference_type k)
{
	return p_diag(k);
}	

inline
const VectorSlice GenMatrix::diag(difference_type k) const
{
	return p_diag(k);
}	

inline
GenMatrixSlice GenMatrix::operator()(const Range1D& I, const Range1D& J)
{
	Range1D Ix = full_range(I,1,rows()), Jx = full_range(J,1,cols());
	return GenMatrixSlice( col_ptr(1) + (Ix.lbound() - 1) + (Jx.lbound() - 1) * rows()
		, max_rows() * cols(), max_rows(), Ix.size(), Jx.size() );
}

inline
const GenMatrixSlice GenMatrix::operator()(const Range1D& I, const Range1D& J) const
{
	Range1D Ix = full_range(I,1,rows()), Jx = full_range(J,1,cols());
	return GenMatrixSlice( const_cast<value_type*>(col_ptr(1)) + (Ix.lbound() - 1) + (Jx.lbound() - 1) * rows()
		, max_rows() * cols(), max_rows(), Ix.size(), Jx.size() );
}

inline
GenMatrixSlice GenMatrix::operator()(size_type i1, size_type i2, size_type j1
	, size_type j2)
{
	return GenMatrixSlice( col_ptr(1) + (i1 - 1) + (j1 - 1) * rows()
		, max_rows() * cols(), max_rows(), i2 - i1 + 1, j2 - j1 + 1 );
}

inline
const GenMatrixSlice GenMatrix::operator()(size_type i1, size_type i2, size_type j1
	, size_type j2) const
{
	return GenMatrixSlice( const_cast<value_type*>(col_ptr(1)) + (i1 - 1) + (j1 - 1) * rows()
		, max_rows() * cols(), max_rows(), i2 - i1 + 1, j2 - j1 + 1 );
}

inline
GenMatrixSlice GenMatrix::operator()()
{
	return GenMatrixSlice( col_ptr(1), max_rows() * cols(), max_rows(), rows(), cols() );
}

inline
const GenMatrixSlice GenMatrix::operator()() const
{
	return GenMatrixSlice( const_cast<value_type*>(col_ptr(1)), max_rows() * cols(), max_rows()
		, rows(), cols() );
}

// Implicit conversion operators

inline
GenMatrix::operator GenMatrixSlice() {
	return (*this)();
}

inline
GenMatrix::operator const GenMatrixSlice() const
{
	return (*this)();
}

// Assignment operators

inline
GenMatrix& GenMatrix::operator=(value_type alpha)
{
	assign(this, alpha);
	return *this;
}

inline
GenMatrix& GenMatrix::operator=(const GenMatrix& rhs)
{
	assign(this, rhs, BLAS_Cpp::no_trans);
	return *this;
}

inline
GenMatrix& GenMatrix::operator=(const GenMatrixSlice& rhs)
{
	assign(this, rhs, BLAS_Cpp::no_trans);
	return *this;
}

// Raw data access

inline
GenMatrix::size_type GenMatrix::max_rows() const
{	return rows_; }

inline
GenMatrix::value_type* GenMatrix::col_ptr(size_type j)
{
	if( v_.size() ) {
		validate_col_subscript(j);
		return &v_[ (j-1) * max_rows() ];
	}
	else {
		return 0;
	}
}

inline
const GenMatrix::value_type* GenMatrix::col_ptr(size_type j) const 
{
	if( v_.size() ) {
		validate_col_subscript(j);
		return &const_cast<valarray&>(v_)[ (j-1) * max_rows() ];
	}
	else {
		return 0;
	}
}

}	// end namespace LinAlgPack

#endif	// GEN_MATRIX_CLASS_H
