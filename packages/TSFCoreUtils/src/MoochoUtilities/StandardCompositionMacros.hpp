// @HEADER
// ***********************************************************************
// 
//      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef STANDARD_COMPOSITION_MACROS_H
#define STANDARD_COMPOSITION_MACROS_H

#include "Teuchos_RefCountPtr.hpp"

///
/** Macro that adds <<std comp>> members for a composition association.
  * \ingroup StandardContainmentMacros_grp
  *
  * This form is for when the object being held will have const attributes
  * the same as the <tt>this</tt> object.
  *
  * For example, if you want to include a <<std comp>> association
  * with an non-const object of type MyClass of the name my_object you
  * would include the macro in the public section of YourClass
  * declaration as follows:
  *
  \verbatim
	class YourClass {
	public:
		STANDARD_COMPOSITION_MEMBERS( MyClass, my_object )
	};
  \endverbatim
  *
  * Note that the macro adds the following type and data members
  * to the class declaration:<br>
  \verbatim
	public:
		Teuchos::RefCountPtr<TYPE> NAME_ptr_t;
	private:
		NAME_ptr_t NAME_;		
  \endverbatim
  *
  * Note that since this macro defines a public type it must be
  * included before any other declarations that may need to refere to
  * the type NAME_ptr_t (such as constuctors)
  */
#define STANDARD_COMPOSITION_MEMBERS( TYPE, NAME )					\
	typedef Teuchos::RefCountPtr<TYPE>				\
		NAME ## _ptr_t;												\
	void set_ ## NAME (const NAME ## _ptr_t& NAME )					\
	{	NAME ## _ = NAME ; }										\
	const NAME ## _ptr_t& get_ ## NAME() const						\
	{	return NAME ## _; }											\
	TYPE& NAME()													\
	{	return *NAME ## _; }										\
	const TYPE& NAME() const										\
	{	return *NAME ## _; }										\
private:															\
	NAME ## _ptr_t	NAME ## _;										\
public:


///
/** Macro that adds <<std comp>> members for a composition association.
  * \ingroup StandardContainmentMacros_grp
  *
  * This form is for when the object being held will have non-const attributes
  * irrespective of the const of <tt>this</tt>.
  *
  * For example, if you want to include a <<std comp>> association
  * with an non-const object of type MyClass of the name my_object you
  * would include the macro in the public section of YourClass
  * declaration as follows:
  *
  \verbatim
	class YourClass {
	public:
		STANDARD_NONCONST_COMPOSITION_MEMBERS( MyClass, my_object )
	};
  \endverbatim
  *
  * Note that the macro adds the following type and data members
  * to the class declaration:<br>
  \verbatim
	public:
		Teuchos::RefCountPtr<TYPE> NAME_ptr_t;
	private:
		NAME_ptr_t NAME_;		
  \endverbatim
  *
  * Note that since this macro defines a public type it must be
  * included before any other declarations that may need to refere to
  * the type NAME_ptr_t (such as constuctors)
  */
#define STANDARD_NONCONST_COMPOSITION_MEMBERS( TYPE, NAME ) \
	typedef Teuchos::RefCountPtr<TYPE> \
		NAME ## _ptr_t; \
	void set_ ## NAME (const NAME ## _ptr_t& NAME ) \
	{	NAME ## _ = NAME ; } \
	const NAME ## _ptr_t& get_ ## NAME() const \
	{	return NAME ## _; } \
	TYPE& NAME() const \
	{	return *NAME ## _; } \
private: \
	NAME ## _ptr_t	NAME ## _; \
public:

///
/** Macro that adds <<std comp>> members for a composition association.
  * \ingroup StandardContainmentMacros_grp
  *
  * This form is for when the object being held will have const attributes
  * irrespective of the const of <tt>this</tt>.
  *
  * For example, if you want to include a <<std comp>> association
  * with a const object of type MyClass of the name my_object you
  * would include the macro in the public section of YourClass
  * declaration as follows:
  *
  \verbatim
	class YourClass {
	public:
		STANDARD_CONST_COMPOSITION_MEMBERS( MyClass, my_object )
	};
  \endverbatim
  *
  * Note that the macro adds the following type and data members
  * to the class declaration:<br>
  \verbatim
	public:
		Teuchos::RefCountPtr<const TYPE> NAME_ptr_t;
	private:
		NAME_ptr_t NAME_;		
  \endverbatim
  *
  * Note that since this macro defines a public type it must be
  * included before any other declarations that may need to refere to
  * the type NAME_ptr_t (such as constuctors)
  */
#define STANDARD_CONST_COMPOSITION_MEMBERS( TYPE, NAME )			\
public:																\
	typedef Teuchos::RefCountPtr<const TYPE>		\
		NAME ## _ptr_t;												\
	void set_ ## NAME (const NAME ## _ptr_t& NAME )					\
	{	NAME ## _ = NAME ; }										\
	NAME ## _ptr_t& get_ ## NAME()									\
	{	return NAME ## _ ; }										\
	const NAME ## _ptr_t& get_ ## NAME() const						\
	{	return NAME ## _; }											\
	const TYPE& NAME() const										\
	{	return *NAME ## _; }										\
private:															\
	NAME ## _ptr_t	NAME ## _;										\
public:

#endif	// STANDARD_COMPOSITION_MACROS_H
