// /////////////////////////////////////////////////////////////////
// AbstractFactoryStd.hpp
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

#ifndef ABSTRACT_FACTORY_PACK_ABSTRACT_FACTORY_STD_H
#define ABSTRACT_FACTORY_PACK_ABSTRACT_FACTORY_STD_H

#include "AbstractFactory.hpp"

namespace MemMngPack {

///
/** Default post-modification class for \c AbstractFactorStd which does nothing!
 */
template<class T_impl>
class PostModNothing {
public:
	///
	void initialize(T_impl* p) const {} // required!
};

///
/** Default allocation class for \c AbstractFactoryStd which returns <tt>new T_impl()</tt>.
 */
template<class T_impl>
class AllocatorNew {
public:
	///
	typedef MemMngPack::ref_count_ptr<T_impl>  ptr_t;                       // required!
	///
	const ptr_t allocate() const { return MemMngPack::rcp(new T_impl()); }  // required!
};

///
/** Simple, templated concrete subclass of universal "Abstract Factory"
 * interface for the creation of objects.
 *
 * The only requirements for the derived classes type is that it allow the
 * default constructor <tt>T_impl()</tt> (i.e. dont make T_impl::T_impl() private)
 * and that it allow <tt>T_impl::new()</tt> and <tt>T_impl::delete</tt> (i.e.
 * don't make them private functions, see Meyers, More Effective C++, Item 27).
 *
 * This subclass is also templated with two other types \c T_PostMod and
 * \c T_Alocator which allow for complete customizataion of object initialization
 * and memory management policies.
 *
 * The type <tt>T_PostMod</tt> is responsible for performing any post
 * modifications on a dynamically allocated object before returning it
 * from \c this->create().  The requirements for the type <tt>T_PostMod</tt>
 * are that it has a default constructor, a copy constructor and a method
 * <tt>T_PostMod::initialize(T_itfc2*) const</tt> that will perform any
 * required post modifications (initializations).  The type \c T_itfc2 argument
 * for this function must be a base class of \c T_impl of course.  The default type
 * for <tt>T_PostMod</tt> is <tt>PostModNothing<T_impl></tt> which does nothing.
 *
 * The type \c T_Allocator allows for specialized memory allocation and cleanup.
*  This type must allow the default constructor and copy constructor and have a method
 * <tt>MemMngPack::ref_count_ptr<T_impl> T_Allocator::allocate() const</tt>
 * which creates a smart reference counted pointer to the allocated object.
 * Also, in returning a <tt>ref_count_ptr<></tt> object, the client can set
 * a function object that can specialize the deallocation of the object (see
 * \c ref_count_ptr).  In defining a specialized \c T_Allocator class, the
 * client can all initialize the object using more than just the default
 * constructor.  Therefore, if the client provides a specialized \c T_Allocator
 * class, there are no restrictions on the class \c T_impl (i.e. does not
 * have to have a default constructor or allow \c new or \c delete).
 * The default class for \c T_Allocator is \c AllocatorNew who's \c allocate()
 * method just returns <tt>MemMngPack::rcp(new T_impl())</tt>.  
 * 
 * Since the \c T_Allocator class can specialize both the memory management
 * and can initialize the object using more that the default constructor, the
 * class \c T_PostMod may seem unecessary.  However, it is more likely that
 * the client will want to define an initialization for a set of classes through
 * an abstract interface and can not for a particular concrete subclass.
 * Also the initialization for an object is orthogonal to how it is created and
 * destroyed, thus the two classes \c T_PostMod and \c T_Allocator.
 */
template<class T_itfc, class T_impl, class T_PostMod = PostModNothing<T_impl>, class T_Allocator = AllocatorNew<T_impl> >
class AbstractFactoryStd : public AbstractFactory<T_itfc> {
public:

#ifdef _MIPS_CXX
	typedef typename MemMngPack::AbstractFactory<T_itfc>::obj_ptr_t   obj_ptr_t;
    typedef T_PostMod                                                 post_mod_t;
    typedef T_Allocator                                               allocator_t;
#endif

    ///
    AbstractFactoryStd( const T_PostMod& post_mod = T_PostMod(), const T_Allocator& alloc = T_Allocator() );

	/** @name Overriden from AbstractFactory */
	//@{
	///
    obj_ptr_t create() const;
	//@}

private:
    T_PostMod    post_mod_;
    T_Allocator  alloc_;

}; // end class AbstractFactorStd

///
template<class T_itfc, class T_impl, class T_Allocator >
const MemMngPack::ref_count_ptr<const AbstractFactory<T_itfc> >
abstract_factory_std_alloc(
	const T_Allocator&  alloc = T_Allocator()
	)
{
	namespace rcp = MemMngPack;
	return rcp::rcp(
		new AbstractFactoryStd<T_itfc,T_impl,PostModNothing<T_impl>,T_Allocator>(
			PostModNothing<T_impl>(), alloc )
		);
}

// ///////////////////////////////////////////////////////
// Template member definitions

template<class T_itfc, class T_impl, class T_PostMod, class T_Allocator>
inline
AbstractFactoryStd<T_itfc,T_impl,T_PostMod,T_Allocator>::AbstractFactoryStd(
	const T_PostMod& post_mod, const T_Allocator& alloc
	)
	:post_mod_(post_mod)
	,alloc_(alloc)
{}

template<class T_itfc, class T_impl, class T_PostMod, class T_Allocator>
inline
AbstractFactoryStd<T_itfc,T_impl,T_PostMod,T_Allocator>::obj_ptr_t
AbstractFactoryStd<T_itfc,T_impl,T_PostMod,T_Allocator>::create() const
{
	namespace rcp = MemMngPack;
	typename T_Allocator::ptr_t
		ptr = alloc_.allocate();
	post_mod_.initialize(ptr.get());
	return ptr;
}

} // end MemMngPack

#endif // ABSTRACT_FACTORY_PACK_ABSTRACT_FACTORY_STD_H
