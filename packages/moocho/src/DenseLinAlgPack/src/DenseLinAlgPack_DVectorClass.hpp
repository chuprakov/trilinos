// ////////////////////////////////////////////////////////////////////////////////////////
// VectorClass.h
//
// 5/20/99: Changed the internal representation of VectorSlice from valarray<value_type> to
//		value_type*.  This makes this class more flexible.  valarray<> was not removed from
//		Vector because of laziness but all access to valarray has been removed so in the
//		future the underlying storage class may change.

#ifndef VECTORCLASS_H
#define VECTORCLASS_H

#include <valarray>

#include "Misc/include/StrideIter.h"
#include "Range1D.h"
#include "VectorAssign.h"

namespace LinAlgPack{

/** @name {\bf Dense 1-D Vector Abstractions}.
  *
  * These are classes that abstract 1-D vectors. The class \Ref{Vector} is a storage class
  * for vectors while the class \Ref{VectorSlice} is used to represent regions of vectors
  * , for rows, columns or diagonals of matrices (see \Ref{GenMatrix}
  * and \Ref{GenMatrixSlice}).
  */

//@{

/** @name {\bf Vector Classes}. */
//@{
	
class Vector;
	
// //////////////////////////////////////////////////////////////////////////////////////
// VectorSlice
//

///
/** Slice of a 1-D sequential C++ array treated as a vector.
  *
  * Objects of this class represent regions of vectors (continuous), rows of matrices
  * , columns of matrices or diagonals of matrices.
  * The underlying representation is of a continuous C++ array with non unit stride.
  * It uses the same convention that the BLAS use where a vector is represented as
  * the first element of
  * in an array, the stride between elements in that array and the number of elements.
  * Changes to elements through a VectorSlice object result in changes to the elements
  * in the underlying value_type* data.
  *
  * VectorSlice provides many STL compliant features such as typedef type members
  *, iterator returning functions
  * and the size() function.  It also provides access to individual elements (lvalue)
  * through 0-based
  * and 1-based subscripting with operator[](i) and operator()(i) respectively.
  *  In addition and subregions can be created
  * by subscripting (operator()()) with an \Ref{Range1D} object or using lower (>0)
  * and upper bounds of the region.
  */
class VectorSlice {
public:
	
	/** @name {\bf Nested Member Types (STL)}.
	  *
	  * These nested types give the types used in the interface to the class.
	  *
	  * \begin{description}
	  *	\item[#value_type#]				- type being stored in the underlying C++ array			
	  *	\item[#size_type#]				- type used as an indice and for the number of elements
	  *										in the vector
	  *	\item[#difference_type#]		- type for the distance between elements and the stride
	  *	\item[#iterator#]				- type for the forward non-constant iterator
	  *	\item[#const_iterator#]			- type for the forward constant iterator (can't change elements)
	  *	\item[#reverse_iterator#]		- type for the reverse non-constant iterator
	  *	\item[#const_reverse_iterator#]	- type for the reverse constant iterator (can't change elements)
	  *	\item[#reference#]				- type returned from subscripting, iterator deferencing etc.
	  *	\item[#const_reference#]		- "" "" for const vector slice objects
	  * \end{description}
	  */

	//@{
	//@}

	typedef LinAlgPack::value_type							value_type;
	typedef LinAlgPack::size_type							size_type;
	typedef ptrdiff_t										difference_type;
	typedef StrideIterPack::stride_iter<value_type*
		, value_type, value_type&, value_type*
		, difference_type>									iterator;
	typedef StrideIterPack::stride_iter<const value_type*
		, value_type, const value_type&, const value_type*
		, difference_type>									const_iterator;
#ifdef _WINDOWS
	typedef std::reverse_iterator<iterator, value_type
		, value_type&, value_type*, difference_type>		reverse_iterator;
	typedef std::reverse_iterator<const_iterator
		, value_type, const value_type&, const value_type*
		, difference_type>									const_reverse_iterator;
#else
	typedef std::reverse_iterator<iterator>					reverse_iterator;
	typedef std::reverse_iterator<const_iterator>			const_reverse_iterator;
#endif
	typedef value_type&										reference;
	typedef const value_type&								const_reference;

	/** @name {\bf Constructors}.
	  *
	  * The user usually does not need not call any of these constructors
	  * explicitly to create a vector slice.
	  * These
	  * constructors are used by the classes in the library to construct VectorSlice objects.
	  * Instead, users create VectorSlice objects by indexing (\Ref{Range1D}) a \Ref{Vector}
	  * , or \Ref{VectorSlice}
	  * object or calling row(...), col(...) or diag(...) on a \Ref{GenMatrix} or
	  * \Ref{GenMatrixSlice} object.
	  * The default C++ copy constructor is used, and is therefore not show here.
	  *
	  * Constructors are also included for creating views of raw C++ arrays.
	  * These constructors take non-#const# pointers.  However you can savely
	  * create a #const# view of a #const# C++ array by using a constant cast.
	  * For example:
	  *
	  \begin{verbatim}
	  const VectorSlice::size_type n = 5;
	  const VectorSlice::value_type ptr[n] = { 1.0, 2.0, 3.0, 4.0, 5,0 };
	  const VectorSlice vec(const_cast<VectorSlice::value_type*>(ptr),n);
	  \end{verbatim}
	  */

	//@{

	///
	/** Creates an empty view.
	  *
	  * You must use bind(...) to bind to a view to initialize after construction.
	  */
	VectorSlice();
	///
	/** Creates a VectorSice object that represents a non-continous region of a raw C++ array.
	  *
	  * Of course the sequence of elements #ptr[stride * i]# for #i# = 0, 1, ..., #size#-1
	  * must yield valid properly allocated regions of memory.  It is up the the user to insure
	  * that they do.
	  *
	  *	@param	ptr		Pointer to the first element in the raw C++ array
	  * @param	size	number of elements in the vector slice
	  * @param	stride	the distance (may be negative) between each successive element (default = 1)
	  */
	VectorSlice( value_type* ptr, size_type size, difference_type stride = 1 );
	///
	/** Creates a VectorSlice object that represents a continous region of a raw C++ array.
	  *
	  * The VectorSlice Object represents the following elements of raw array:
	  * 
	  *      #ptr[rng.lbound()-1+i]#, for #i# = 0, 1, ..., #rng.ubound()-1#
	  * 
	  * Preconditions: \begin{itemize}
	  *		\item #rng.ubound() + 1 <= n# (throw std::out_of_range)
	  *		\end{itemize}
	  *
	  * @param	v		valarray<> the vector slice region is representing.
	  * @param	rng		index range (1-based) of the region being represented.  
	  *					Here rng.full_range() can be true.
	  */
	VectorSlice( value_type* ptr, size_type size, const Range1D& rng );
	/// 
	/** Create a VectorSlice that represents a continous region of the existing VectorSlice object, vs.
	  *
	  * The index, rng, is relative to the VectorSlice object, vs.
	  * For example rng = [1,3] would create a VectorSlice object
	  * representing the elements 2, 4 and 6.  The following
	  * shows the elements represented by each of the participating objects.
	  \begin{verbatim}

	  vs =   [2, 4, 6, 8, 10]
	  this = [2, 4, 6]

	  \end{verbatim}

	  * Preconditions: \begin{itemize}
	  *		\item rng.full_range() == false (throw #std::out_of_range#)
	  *		\item rng.size() <= vs.size() (throw #std::out_of_range#) 
	  *		\end{itemize}
	  * 
	  * @param	vs		VectorSlice object that this VectorSlice object is being created from
	  * @param	rng		Range1D range of the vector slice being created.
	  */
	VectorSlice( VectorSlice& vs, const Range1D& rng );

	//@}

	/// Bind to the view of another VectorSlice
	void bind(VectorSlice vs);

	/** @name {\bf STL Iterator Access Functions}.
	  *
	  * These member functions return valid STL random access iterators to the elements in the
	  * VectorSlice object.
	  *
	  * The forward iterators returned by begin() and end() iterator sequentialy from the first
	  * element (same element as returned by operator()(1)) to the last
	  * element (same element as returned by operator()(size()).  This goes for reverse 
	  * (stride() < 0) VectorSlice objects as well.  The reverse iterators returned by
	  * rbegin() and rend() iterate in the reverse sequence.
	  *
	  * Warning! Beware of using iterations in a reverse vector slice (stride() < 0). 
	  * In a reverse vector slice end() returns a slice iterator which is the current
	  * element is one before the first allocated element.  Strictly speaking this is
	  * not allowed so use iterators with reversed VectorSlice objects at own risk.
	  */

	//@{

	/// 
	iterator begin();
	///
	iterator end();
	///
	const_iterator begin() const;
	///
	const_iterator end() const;
	///
	reverse_iterator rbegin();
	///
	reverse_iterator rend();
	///
	const_reverse_iterator rbegin() const;
	///
	const_reverse_iterator rend() const;

	//@}

	/** @name {\bf Individual Element Access Subscripting (lvalue)}.
	  *
	  * These operator functions allow access (lvalue) to the individual elements
	  * of the VectorSlice object.
	  * 
	  * The subscript i must be, 1 <= i <= this->size(), for the 1-based element access
	  * operators and, 0 <= i <= this->size() - 1, for the 0-based element access operators.
	  * If they are not then an #std::out_of_range# exception will be thrown.
	  */

	//@{

	/// 1-based element access (lvalue)
	reference operator()(size_type i);
	/// 1-based element access (rvalue)
	const_reference operator()(size_type i) const;
	/// 1-based element access (lvalue)
	reference operator[](size_type i);
	/// 0-based element access (rvalue)
	const_reference operator[](size_type i) const;

	//@}

	/** @name {\bf Subvector Access Operators}.
	  *
	  * These operator functions are used to create views of continous regions of the VectorSlice.
	  * Each of them returns a VectorSlice object for the region.  Constant (const) VectorSlice objects
	  * are returned for a const VectorSlice.  This means that the elements can not be changed
	  * as should be the case.
	  * 
	  * Beware!  VC++ is returning non-const VectorSlice objects for the 
	  * #VectorSlice operator()(...) const;# member functions and therefore a const \Ref{Vector} or
	  * \Ref{VectorSlice} can be modifed my subsetting it.  Hopefully this problem will
	  * be fixed in future versions of the compiler or I when will get another compiler.
	  */

	//@{

	/// 
	/** Returns a VectorSlice object representing the entire vector slice.
	  *
	  * Included for uniformity with vector.
	  */
	/// Allow the address to be taken of an rvalue of this object.
	VectorSlice* operator&() {
	  return this;
	}
	///
	const VectorSlice* operator&() const {
	  return this;
	}
	VectorSlice& operator()();
	/// Same as above.
	const VectorSlice& operator()() const;
	/// 
	VectorSlice operator()(const Range1D& rng);
	/// 
	/** Returns a continous subregion of the VectorSlice object.
	  *
	  * The returned VectorSlice object represents the range of the rng argument.
	  *
	  * Preconditions: \begin{itemize}
	  *		\item #rng.ubound() - 1 <= this->size()# (throw #out_of_range#)
	  *		\end{itemize}
	  *
	  * @param	rng		Indece range [lbound,ubound] of the region being returned.
	  */
	const VectorSlice operator()(const Range1D& rng) const;
	/// 
	/** Returns a VectorSlice object for the continous subregion [ubound, lbound].
	  * 
	  * Preconditions: \begin{itemize}
	  *		\item #lbound > 1# (throw out_of_range)
	  *		\item #lbound < ubound# (throw out_of_range)
	  *		\item #ubound <= this->size()# (throw out_of_range)
	  *		\end{itemize}
	  *
	  * @param	rng		Range [lbound,ubound] of the region being returned.
	  */
	VectorSlice operator()(size_type lbound, size_type ubound);
	/// Same as above.
	const VectorSlice operator()(size_type lbound, size_type ubound) const;
	/// 
	/** Return a const VectorSlice object the reverse of this VectorSlice.
	  *
	  * In the reverse VectorSlice,
	  * the first element becomes the last element and visa-versa.  For example, for 
	  * #VectorSlice r = x.rev()#, #&x(1) == &z(z.size())# and #&x(x.size()) == &z(1)# are both true.
	  * The iterators returned by \Ref{begin()} iterate from the first conceptual element to the last.
	  */
	VectorSlice rev();
	/// Same as above.
	const VectorSlice rev() const;

	//@}
	
	/** @name {\bf Assignment operators}. */

	//@{

	///
	/** vs = alpha (Sets all the elements to the constant alpha).
	  *
	  * Preconditions: \begin{itemize}
	  *		\item #this->size() > 0# (throw #std::length_error#)
	  *		\end{itemize}
	  *
	  * Postconditions: \begin{itemize}
	  *		\item #this->operator()(i) == alpha#, i = 1, 2, ... , #this->size()#
	  *		\end{itemize}
	  */
	VectorSlice& operator=(value_type alpha);
	///
	/** vs = rhs (Copies the elements of rhs into the elements of this).
	  *
	  * Preconditions: \begin{itemize}
	  *		\item #this->size() == rhs.size()# (throw #out_of_range#)
	  *		\item #rhs.size() > 0# (throw #out_of_range#)
	  *		\end{itemize}
	  *
	  * Postconditions: \begin{itemize}
	  *		\item #this->operator()(i) == rhs(i)#, i = 1, 2, ..., #this->size()#
	  *		\end{itemize}
	  */
	VectorSlice& operator=(const VectorSlice& rhs);

	//@}

	/** @name {\bf Misc. Member Functions}. */

	//@{

	/// Returns the number of elements of the VectorSlice.
	size_type size() const;
	/// 
	/** Returns the degree of memory overlap of the two VectorSlice objects this and vs.
	  *
	  * @return 
	  *		\begin{description}
	  *		\item[NO_OVERLAP]	There is no memory overlap between this and vs
	  *		\item[SOME_OVERLAP]	There is some memory locations that this and vs share
	  *		\item[SAME_MEM]		The VectorSlice objects this and vs share the exact same memory locations.
	  *		\end{description}
	  */
	EOverLap overlap(const VectorSlice& vs) const;

	//@}

	/** @name {\bf Raw data access}.
	  *
	  * Provides access to underlying raw data.
	  */

	//@{

	/// Return a pointer to the address of the first memory location of underlying array.
	value_type*			raw_ptr();
	///
	const value_type*	raw_ptr() const;
	/// Return a pointer to the conceptual first element in the underlying array.
	value_type*			start_ptr();
	///
	const value_type*	start_ptr() const;
	/// Return the distance (+,-) (in units of elements) between adjacent elements in the underlying array.
	difference_type		stride() const;

	//@}

private:
	value_type					*ptr_;	// Pointer to first element in array.
	size_type					size_;  // # elements represented in v_
	difference_type				stride_;// # positions to skip between elements. Must be positive

	// Above the sequence represented is:
	//  ptr_[ i * stride_ ], for i = 0, ..., size_ - 1

	// Helper functions
	// Assert a valarray is sized on construction
	size_type validate_sized(size_type size);
	// Assert a range if valid on construction
	void validate_range(size_type ubound, size_type max_ubound);
	// Assest a subscript is in bounds
	void validate_subscript(size_type i) const;

}; // end class VectorSlice

// /////////////////////////////////////////////////////////////////////////////////////////
// Vector
//

///
/** 1-D Vector Abstraction Storage Class.
  *
  * Holds the storage space for a 1-D vector of element type value_type.  The storage space class
  * used in a standard valarray<> private member.  Vector provides much of the
  * same functionaliy of a VectorSlice object accept that Vector object can be resized at any time by
  * either explicitly calling #resize(...)# or to match an assignment to a rhs linear algebra expression.
  */
class Vector {
public:
	/** @name {\bf Nested Member Types (STL)}.
	  *
	  * These nested types give the types used in the interface to the class.
	  *
	  * \begin{description}
	  *	\item[#value_type#]				- type being stored in the underlying valarray<>			
	  *	\item[#size_type#]				- type for the number of elements in the valarray<>
	  *	\item[#difference_type#]		- type for the distance between elements
	  *	\item[#iterator#]				- type for the forward non-constant iterator
	  *	\item[#const_iterator#]			- type for the forward constant iterator (can't change elements)
	  *	\item[#reverse_iterator#]		- type for the reverse non-constant iterator
	  *	\item[#const_reverse_iterator#]	- type for the reverse constant iterator (can't change elements)
	  *	\item[#reference#]				- type returned from subscripting, iterator deferencing etc.
	  *	\item[#const_reference#]		- "" "" for const vector slice objects
	  * \end{description}
	  */

	//@{
	//@}

	typedef LinAlgPack::value_type					value_type;
	typedef LinAlgPack::size_type					size_type;
	typedef ptrdiff_t								difference_type;
	typedef value_type*								iterator;
	typedef const value_type*						const_iterator;
#ifdef _WINDOWS
	typedef std::reverse_iterator<iterator, value_type
		, value_type&, value_type*, difference_type>		reverse_iterator;
	typedef std::reverse_iterator<const_iterator
		, value_type, const value_type&, const value_type*
		, difference_type>									const_reverse_iterator;
#else
	typedef std::reverse_iterator<iterator>					reverse_iterator;
	typedef std::reverse_iterator<const_iterator>			const_reverse_iterator;
#endif
	typedef value_type&								reference;
	typedef const value_type&						const_reference;
	typedef std::valarray<value_type>				valarray;

	/** @name {\bf Constructors}.
	  * 
	  * These constructors allocate and may initialize the elements of a 1-D vector.
	  * The default C++ copy constructor is used and is therefore not show here.
	  */

	//@{

	/// Constructs a vector with 0 elements (this->size()==0).
	Vector();
	/// Constructs a vector with n elements of initialized memory.
	Vector(size_type n);
	/// Constructs a vector with n elements initialized to val.
	Vector(value_type val, size_type n);
	///
	/** Constructs a vector with n elements and intializes elements to those of an array.
	  *
	  * Postconditions: \begin{itemize}
	  *		\item #this->operator[](i) == p[i]#, i = 0, 1, ... n
	  *		\end{itemize}
	  */  
	Vector(const value_type* p, size_type n);
	///
	/** Constructs a Vector object fron a VectorSlice object.
	  *
	  * Postconditions: \begin{itemize}
	  *		\item #this->size() == vs.size()#
	  *		\item #this->operator[](i) == vs[i]#, i = 0, 1, ... n
	  *		\end{itemize}
	  */  
	Vector(const VectorSlice& vs);

	//@}

	/** @name {\bf Memory Management / Misc}. */

	//@{

	///
	/** Resize the vector to hold n elements.
	  *
	  * Any new elements added are initialized to val.
	  *
	  * Postconditions: \begin{itemize}
	  *		\item #this->size() == n#
	  *		\end{itemize}
	  */  
	void resize(size_type n, value_type val = value_type());
	///
	/** Free memory and resize Vector to this->size() == 0.
	  *
	  * Postconditions: \begin{itemize}
	  *		\item #this->size() == 0#
	  *		\end{itemize}
	  */  
	void free();
	/// Returns the number of elements of the Vector.
	size_type size() const;
	/// 
	/** Returns the degree of memory overlap of this and the VectorSlice object vs.
	  *
	  * @return 
	  *		\begin{description}
	  *		\item[NO_OVERLAP]	There is no memory overlap between this and vs
	  *		\item[SOME_OVERLAP]	There is some memory locations that this and vs share
	  *		\item[SAME_MEM]		The VectorSlice objects this and vs share the exact same memory locations.
	  *		\end{description}
	  */
	EOverLap overlap(const VectorSlice& vs) const;
	/// Conversion operator for implicit conversions from Vector to VectorSlice.
	operator VectorSlice();
	/// Conversion operator for implicit conversions from const Vector to const VectorSlice.
	operator const VectorSlice() const;

	//@}
	

	/** @name {\bf STL Iterator Access functions}.
	  *
	  * The iterators returned are valid STL random access iterators.
	  * The forward iterators returned iterate from the first element to the last element.
	  * The reverse iterators returned iterate from the last element to the first element.
	  */

	//@{

	///
	iterator begin();
	///
	iterator end();
	///
	const_iterator begin() const;
	///
	const_iterator end() const;
	///
	reverse_iterator rbegin();
	///
	reverse_iterator rend();
	///
	const_reverse_iterator rbegin() const;
	///
	const_reverse_iterator rend() const;

	//@}

	/** @name {\bf Individual Element Access Subscripting (lvalue)}.
	  *
	  * These operator functions allow access (lvalue) to the individual elements
	  * of the Vector object.
	  * 
	  * The subscript i must be, 1 <= i <= this->size(), for the 1-based element access
	  * operators and, 0 <= i <= this->size() - 1, for the 0-based element access operators.
	  * If they are not then an #std::out_of_range# exception will be thrown.
	  */

	//@{

	/// 1-based element access (lvalue)
	reference operator()(size_type i);
	/// 1-based element access (rvalue)
	const_reference operator()(size_type i) const;
	/// 1-based element access (lvalue)
	reference operator[](size_type i);
	/// 0-based element access (rvalue)
	const_reference operator[](size_type i) const;

	//@}

	/** @name {\bf Subvector Access Operators}.
	  *
	  * These operator functions are used to create views of continous regions of the Vector.
	  * Each of them returns a VectorSlice object for the region.  Constant (const) VectorSlice objects
	  * are returned for a const Vector.  This means that the elements can not be changed
	  * as should be the case.
	  * 
	  * Beware!  VC++ is returning non-const VectorSlice objects for the 
	  * #VectorSlice operator()(...) const;# member functions and therefore a const \Ref{Vector} or
	  * \Ref{VectorSlice} can be modifed my subsetting it.  Hopefully this problem will
	  * be fixed in future versions of the compiler or I when will get another compiler.
	  */

	//@{

	/// 
	/** Returns a VectorSlice object representing the entire Vector. 
	  * 
	  * Call this member function to force a type conversion to VectorSlice.  Using the
	  * VectorSlice of a Vector for algebraic expressions used with the TCOL allows a for simplier
	  * implementaion of those operations by cutting down on the number combinations.  This is
	  * especialy true for longer optimized expression.
	  */
	VectorSlice operator()();
	/// Same as above
	const VectorSlice operator()() const;
	/// 
	/** Returns a continous subregion of the Vector object.
	  *
	  * The returned VectorSlice object represents the range of the rng argument.
	  *
	  * Preconditions: \begin{itemize}
	  *		\item #rng.ubound() - 1 <= this->size()# (throw #out_of_range#)
	  *		\end{itemize}
	  *
	  * @param	rng		Indece range [lbound,ubound] of the region being returned.
	  */
	VectorSlice operator()(const Range1D& rng);
	/// Same as above
	const VectorSlice operator()(const Range1D& rng) const;
	/// 
	/** Returns a VectorSlice object for the continous subregion [ubound, lbound].
	  * 
	  * Preconditions: \begin{itemize}
	  *		\item #lbound > 1# (throw #out_of_range#)
	  *		\item #lbound < ubound# (throw #out_of_range#)
	  *		\item #ubound <= this->size()# (throw #out_of_range#)
	  *		\end{itemize}
	  *
	  * @param	rng		Range [lbound,ubound] of the region being taken.
	  */
	VectorSlice operator()(size_type lbound, size_type ubound);
	/// Same as above.
	const VectorSlice operator()(size_type lbound, size_type ubound) const;
	/// 
	/** Return a VectorSlice object the reverse of this Vector.
	  *
	  * In the reverse VectorSlice,
	  * the first element becomes the last element and visa-versa.  For example, for 
	  * #VectorSlice r = x.rev()#, #&x(1) == &z(z.size())# and #&x(x.size()) == &z(1)# are both true.
	  * The iterators returned by \Ref{begin()} iterate from the first conceptual element to the last.
	  */
	VectorSlice rev();
	/// Same as above.
	const VectorSlice rev() const;

	//@}

	/** @name {\bf Assignment Operators}. */

	//@{

	///
	/** vs = alpha (Sets all the elements to the constant alpha).
	  *
	  * Preconditions: \begin{itemize}
	  *		\item #this->size() > 0# (throw #std::length_error#)
	  *		\end{itemize}
	  *
	  * Postconditions: \begin{itemize}
	  *		\item #this->operator()(i) == alpha#, i = 1, 2, ... , #this->size()#
	  *		\end{itemize}
	  */
	Vector& operator=(value_type alpha);
	///
	/** vs = rhs (Copies the elements of rhs into the elements of this).
	  *
	  * Preconditions: \begin{itemize}
	  *		\item #this->size() == rhs.size()# (throw #out_of_range#)
	  *		\item #rhs.size() > 0# (throw #out_of_range#)
	  *		\end{itemize}
	  *
	  * Postconditions: \begin{itemize}
	  *		\item #this->operator()(i) == rhs(i)#, i = 1, 2, ..., #this->size()#
	  *		\end{itemize}
	  */
	Vector& operator=(const VectorSlice& rhs);
	///
	/** Needed to override the default assignment operator.
	  */
	Vector& operator=(const Vector& rhs);

	//@}

	/** @name {\bf Implementation Access}.
	  *
	  * Provides access to underlying raw data.
	  */

	//@{

	/// Return a pointer to the address of the first memory location of underlying array.
	value_type*			raw_ptr();
	///
	const value_type*	raw_ptr() const;
	/// Return a pointer to the conceptual first element in the underlying array.
	value_type*			start_ptr();
	///
	const value_type*	start_ptr() const;
	/// Return the distance (+,-) (in units of elements) between adjacent elements in the underlying array.
	difference_type		stride() const;

	//@}
	
private:
	valarray v_;

	void validate_subscript(size_type i) const;

}; // end class Vector

//		end Vector Classes scope
//@}

// ///////////////////////////////////////////////////////////////////////////////
// Non-member function declarations												//
// ///////////////////////////////////////////////////////////////////////////////


/** @name {\bf Non-Member Functions}. */

//@{
//		begin non-member functions scope

inline
///
/** Utility for checking the sizes of two VectorSlice objects and throwing an exception
  * if the sizes are not the same.
  */
void assert_vs_sizes(const VectorSlice& v1, const VectorSlice& v2)
{
#ifdef LINALGPACK_CHECK_RHS_SIZES
	if(v1.size() != v2.size())
		throw std::length_error("Sizes of vector regions must be the same.\n");
#endif
}

inline
///
/** Create a general vector slice.
  */
VectorSlice gen_vs( VectorSlice& vs, size_type start, size_type size, ptrdiff_t stride )
{
	return VectorSlice( vs.start_ptr() + vs.stride() * (start-1), size, vs.stride() * stride );
}

inline
///
const VectorSlice gen_vs( const VectorSlice& vs, size_type start, size_type size
	, ptrdiff_t stride )
{
	return VectorSlice( const_cast<VectorSlice::value_type*>(vs.start_ptr()) + vs.stride() * (start-1)
		, size, vs.stride() * stride );
}

//		end non-member functions scope
//@}

//		end Vectors scope
//@}

// ////////////////////////////////////////////////////////////////////////////////
// Inline definitions of member function definitions							 //
// ////////////////////////////////////////////////////////////////////////////////

// /////////////////////////////////////////////////////////////////////////////
// VectorSlice inline member function definitions

#ifndef LINALGPACK_CHECK_SLICE_SETUP
inline
VectorSlice::size_type VectorSlice::validate_sized(size_type size)
{
	return size;
}
#endif

#ifndef LINALGPACK_CHECK_SLICE_SETUP
inline
void VectorSlice::validate_range(size_type ubound, size_type max_ubound)
{}
#endif

#ifndef LINALGPACK_CHECK_SLICE_SETUP
inline
void VectorSlice::validate_subscript(size_type i) const
{}
#endif

inline
VectorSlice::difference_type my_abs(VectorSlice::difference_type i) { return i > 0 ? i : -i; }

// Constructors.  Use default copy constructor

inline
VectorSlice::VectorSlice()
	: ptr_(0)
	, size_(0)
	, stride_(0)
{}

inline
VectorSlice::VectorSlice( value_type* ptr, size_type size, difference_type stride)
	: ptr_(ptr)
	, size_(size)
	, stride_(stride)
{}

inline
VectorSlice::VectorSlice( value_type* ptr, size_type size, const Range1D& rng )
	: ptr_( ptr + rng.lbound() - 1 )
	, size_( rng.full_range() ? validate_sized(size) : rng.size() )
	, stride_(1)
{	validate_range( rng.full_range() ? size : rng.ubound(), size ); }

inline
VectorSlice::VectorSlice( VectorSlice& vs, const Range1D& rng )
	: ptr_( vs.start_ptr() + (rng.lbound() - 1) * vs.stride() )
	, size_( rng.full_range() ? validate_sized(vs.size()) : rng.size() )
	, stride_( vs.stride() )
{	validate_range(  rng.full_range() ? vs.size() : rng.ubound(), vs.size() ); }

inline
void VectorSlice::bind(VectorSlice vs)
{
	ptr_	= vs.ptr_;
	size_	= vs.size_;
	stride_	= vs.stride_;
}

// Iterator functions
inline
VectorSlice::iterator	VectorSlice::begin()
{	return iterator(start_ptr(), stride()); }

inline
VectorSlice::iterator	VectorSlice::end()
{	return iterator(start_ptr() + size() * stride(), stride()); }

inline
VectorSlice::const_iterator VectorSlice::begin() const
{	return const_iterator(start_ptr(), stride()); }

inline
VectorSlice::const_iterator VectorSlice::end() const
{	return const_iterator(start_ptr() + size() * stride(), stride()); }

inline
VectorSlice::reverse_iterator	VectorSlice::rbegin()
{	return reverse_iterator(end()); }

inline
VectorSlice::reverse_iterator	VectorSlice::rend()
{	return reverse_iterator(begin()); }

inline
VectorSlice::const_reverse_iterator VectorSlice::rbegin() const
{	return const_reverse_iterator(end()); }

inline
VectorSlice::const_reverse_iterator VectorSlice::rend() const
{	return const_reverse_iterator(begin()); }

// Element access
inline
VectorSlice::reference VectorSlice::operator()(size_type i) // 1 based
{
	validate_subscript(i);
	return ptr_[(i-1)*stride_];
}

inline
VectorSlice::const_reference VectorSlice::operator()(size_type i) const
{
	validate_subscript(i);
	return ptr_[(i-1)*stride_];
}

inline
VectorSlice::reference VectorSlice::operator[](size_type i) // 0 based		
{
	validate_subscript(i+1);
	return ptr_[(i)*stride_];
}

inline
VectorSlice::const_reference VectorSlice::operator[](size_type i) const
{
	validate_subscript(i+1);
	return ptr_[(i)*stride_];
}

// Subregion Access.  Let the constructors of VectorSlice validate the ranges
inline
VectorSlice& VectorSlice::operator()() 
{	return *this; }

inline
const VectorSlice& VectorSlice::operator()() const 
{	return *this; }

inline
VectorSlice VectorSlice::operator()(const Range1D& rng) 
{	return VectorSlice(*this, full_range(rng,1,size())); }

inline
const VectorSlice VectorSlice::operator()(const Range1D& rng) const
{	return VectorSlice(const_cast<VectorSlice&>(*this), full_range(rng,1,size())); }

inline
VectorSlice VectorSlice::operator()(size_type lbound, size_type ubound)
{	return VectorSlice(*this, Range1D(lbound, ubound)); }

inline
const VectorSlice VectorSlice::operator()(size_type lbound, size_type ubound) const
{	return VectorSlice(const_cast<VectorSlice&>(*this), Range1D(lbound, ubound)); }

inline
VectorSlice VectorSlice::rev()
{	return VectorSlice( start_ptr() + stride() * (size()-1), size(), - stride() ); }

inline
const VectorSlice VectorSlice::rev() const
{	return VectorSlice( const_cast<value_type*>(start_ptr()) + stride() * (size()-1), size(), - stride() ); }

// Assignment Operators
inline
VectorSlice& VectorSlice::operator=(value_type alpha) 
{	assign(this, alpha); return *this; }

inline
VectorSlice& VectorSlice::operator=(const VectorSlice& rhs) 
{	assign(this, rhs); return *this; }

// Misc. member functions

inline
VectorSlice::size_type VectorSlice::size() const
{	return size_; }

// Raw pointer access

inline
VectorSlice::value_type*	VectorSlice::raw_ptr()
{	return stride() > 0 ? start_ptr() : start_ptr() + stride() * (size() - 1); }

inline
const VectorSlice::value_type* VectorSlice::raw_ptr() const
{	return stride() > 0 ? start_ptr() : start_ptr() + stride() * (size() - 1); }

inline
VectorSlice::value_type*	VectorSlice::start_ptr()
{	return ptr_; }

inline
const VectorSlice::value_type* VectorSlice::start_ptr() const
{	return ptr_; }

inline
VectorSlice::difference_type VectorSlice::stride() const
{	return stride_; }

// /////////////////////////////////////////////////////////////////////////////
// Vector inline member function definitions

#ifndef LINALGPACK_CHECK_RANGE
inline
void Vector::validate_subscript(size_type i) const
{}
#endif

// Constructors
inline
Vector::Vector()
{}	// used to shut satisfy compiler

inline
Vector::Vector(size_type n) 
	: v_(n)
{}

inline
Vector::Vector(value_type val, size_type n) 
	: v_(val, n)
{}

inline
Vector::Vector(const value_type* p, size_type n)
	: v_(p, n)
{}

inline
Vector::Vector(const VectorSlice& vs)
	: v_(vs.size())
{	*this = vs;	}

// Memory management
inline
void Vector::resize(size_type n, value_type val)
{	v_.resize(n,val); }

inline
void Vector::free()
{
	v_.resize(0);
}

// Size
inline
Vector::size_type Vector::size() const
{	return v_.size(); }

// Iterator functions
inline
Vector::iterator Vector::begin()
{	return start_ptr(); }

inline
Vector::iterator Vector::end()
{	return start_ptr() + size(); }

inline
Vector::const_iterator Vector::begin() const
{	return start_ptr(); }

inline
Vector::const_iterator Vector::end() const 
{	return start_ptr() + size(); }

inline
Vector::reverse_iterator	Vector::rbegin()
{	return reverse_iterator(end()); }

inline
Vector::reverse_iterator	Vector::rend()
{	return reverse_iterator(begin()); }

inline
Vector::const_reverse_iterator Vector::rbegin() const
{	return const_reverse_iterator(end()); }

inline
Vector::const_reverse_iterator Vector::rend() const 
{	return const_reverse_iterator(begin()); }

// Element access
inline
Vector::reference Vector::operator()(size_type i)
{
	validate_subscript(i);
	return start_ptr()[i-1];
}

inline
Vector::const_reference Vector::operator()(size_type i) const
{
	validate_subscript(i);
	return start_ptr()[i-1];
}

inline
Vector::reference Vector::operator[](size_type i)
{
	validate_subscript(i+1);
	return start_ptr()[i];
}

inline
Vector::const_reference Vector::operator[](size_type i) const
{
	validate_subscript(i+1);
	return start_ptr()[i];
}

// Subregion Access.  Leave validation to the VectorSlice constructors.
inline
VectorSlice Vector::operator()() 
{	return VectorSlice(start_ptr(),size()); }

inline
const VectorSlice Vector::operator()() const
{	return VectorSlice(const_cast<value_type*>(start_ptr()),size()); }

inline
VectorSlice Vector::operator()(const Range1D& rng) 
{	return VectorSlice(start_ptr(),size(),rng); }

inline
const VectorSlice Vector::operator()(const Range1D& rng) const
{	return VectorSlice(const_cast<value_type*>(start_ptr()),size(),rng); }

inline
VectorSlice Vector::operator()(size_type lbound, size_type ubound)
{	return VectorSlice(start_ptr(), size(), Range1D(lbound, ubound)); }

inline
const VectorSlice Vector::operator()(size_type lbound, size_type ubound) const
{	return VectorSlice(const_cast<value_type*>(start_ptr()), size(), Range1D(lbound, ubound)); }

inline
VectorSlice Vector::rev()
{	return VectorSlice( start_ptr() + size() - 1, size(), -1 ); }

inline
const VectorSlice Vector::rev() const
{	return VectorSlice( const_cast<value_type*>(start_ptr()) + size() - 1, size(), -1 ); }

// Conversion operators
inline
Vector::operator VectorSlice()
{	return VectorSlice(start_ptr(), size()); }

inline
Vector::operator const VectorSlice() const
{	return VectorSlice(const_cast<value_type*>(start_ptr()), size()); }

// Assignment Operators
inline
Vector& Vector::operator=(value_type alpha) 
{
	assign(this, alpha);
	return *this;
}

inline
Vector& Vector::operator=(const Vector& rhs) 
{
	assign(this, rhs);
	return *this;
}

inline
Vector& Vector::operator=(const VectorSlice& rhs) 
{
	assign(this, rhs);
	return *this;
}

// Raw pointer access

inline
Vector::value_type*	Vector::raw_ptr()
{	return start_ptr(); }

inline
const Vector::value_type* Vector::raw_ptr() const
{	return start_ptr(); }

inline
Vector::value_type*	Vector::start_ptr()
{	return size() ? &(v_)[0] : 0; }

inline
const Vector::value_type* Vector::start_ptr() const
{	return &const_cast<valarray&>((v_))[0]; }

inline
Vector::difference_type Vector::stride() const
{	return 1; }

} // end namespace LinAlgPack

// //////////////////////////////////////////////////////////////////////////////////////////////////////
// Non-member functions / utilities

#endif	// end VECTORCLASS_H
