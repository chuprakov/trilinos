// ////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_compare_element_indexes.hpp
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

#ifndef COMPARE_ELEMENT_INDEXES_H
#define COMPARE_ELEMENT_INDEXES_H

#include <functional>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** @name Comparision function objects for sparse element objects.
  *
  * These function objects are used to compare element indexes so the
  * templated objects must support the index() operation.
  */
//@{

/// ele1.index() < ele2.index() 
template<class T_Element>
struct compare_element_indexes_less
	: public std::unary_function< T_Element, typename T_Element::index_type >
{
	bool operator()(const T_Element& ele1, const T_Element& ele2) {
		return ele1.index() < ele2.index();
	}
	bool operator()(typename T_Element::index_type i, const T_Element& ele2) {
		return i < ele2.index();
	}
	bool operator()(const T_Element& ele1, typename T_Element::index_type i) {
		return ele1.index() < i;
	}
};

/// ele.index() == i
template<class T_Element>
struct compare_element_indexes_equal_to
	: public std::unary_function< T_Element, typename T_Element::index_type >
{
	compare_element_indexes_equal_to( typename T_Element::index_type i )
		: i_(i)
	{}
	bool operator()(const T_Element& ele) {
		return ele.index() == i_;
	}
private:
	typename T_Element::index_type i_;
	// Not defined and not to be called
	compare_element_indexes_equal_to();
};

//@}

}	// end namespace AbstractLinAlgPack

#endif   // COMPARE_ELEMENT_INDEXES_H
