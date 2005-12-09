// //////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_SparsePtrElement.hpp
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

#ifndef SPARSE_PTR_ELEMENT_H
#define SPARSE_PTR_ELEMENT_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

///
/** Sparse pointer element type.
  *
  * This class abstracts a sparse element of a templated
  * type.  It is ment to be used in a sparse vector.  It
  * has a pointer to the value of the element.
  *
  * The default assignment operator and copy constructor
  * are allowed.
  */
template <class T_Indice, class T_Value>
class SparsePtrElement {
public:
	/** @name Public Typedefs. */
	//@{

	///
	typedef T_Value							value_type;
	///
	typedef T_Indice						indice_type;

	//@}

	/** @name Constructors */
	//@{

	/// Construct uninitialized (poiner to value set to zero) (#indice() == 0#).
	SparsePtrElement() : indice_(0), pvalue_(0)
	{}

	/// Construct with a pointer to the value and indice set
	SparsePtrElement(indice_type indice, value_type* pvalue) : indice_(indice), pvalue_(pvalue)
	{}
	
	//@}

	/** @name Value and indice access */
	//@{ 

	///
	value_type& value()
	{
		return *pvalue_;
	}
	///
	value_type value() const
	{
		return *pvalue_;
	}
	///
	indice_type indice() const
	{
		return indice_;
	}
	/// Change the indice
	void change_indice(indice_type indice)
	{
		indice_ = indice;
	}
	/// Change the element pointer
	void change_value_ptr(value_type* pvalue)
	{
		pvalue_ = pvalue;
	}

	//@}
private:
	indice_type				indice_;
	value_type*				pvalue_;

};	// end class SparsePtrElement

} // end namespace AbstractLinAlgPack 

#endif // SPARSE_PTR_ELEMENT_H
