// ///////////////////////////////////////////////////////////////
// AbstractLinAlgPack_sparse_bounds.hpp
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

#ifndef SPARSE_LIN_ALG_PACK_SPARSE_BOUNDS_H
#define SPARSE_LIN_ALG_PACK_SPARSE_BOUNDS_H

#include "DenseLinAlgPack_DVectorClass.hpp"

namespace AbstractLinAlgPack {

///
/** Iterate through a set of sparse bounds.
 *
 * Finish documentation.
 *
 * Allow default copy constructor and assignment operator.
 */
class sparse_bounds_itr {
private:
	enum EBound { LOWER, UPPER, BOTH };
public:

	///
	sparse_bounds_itr(	
		const DVectorSlice::const_iterator   &bl_begin
		,const DVectorSlice::const_iterator  &bl_end
		,const DVectorSlice::const_iterator  &bu_begin
		,const DVectorSlice::const_iterator  &bu_end
		,value_type                         big_bnd
		)
		:bl_itr_(bl_begin), bl_end_(bl_end), bu_itr_(bu_begin), bu_end_(bu_end)
		,big_bnd_(big_bnd), index_(1)
	{
		if( !at_end() && ( *bl_itr_ <= -big_bnd_ && +big_bnd_ <= *bu_itr_ ) )
			this->operator++();
		else
			update();
	}
	///
	const value_type& big_bnd() const
	{    return big_bnd_; }
	///
	bool at_end() const
	{	return bl_itr_ == bl_end_; }
	///
	sparse_bounds_itr& operator++() {
		if(!at_end()) { ++bl_itr_; ++bu_itr_; ++index_; }
		for( ; !at_end() && ( *bl_itr_ <= -big_bnd_ && +big_bnd_ <= *bu_itr_ )
			 ; ++bl_itr_, ++bu_itr_, ++index_ );
		update();
		return *this;
	}
	///
	index_type index() const
	{	return index_; }
	///
	value_type lbound() const
	{	return lbound_; }
	///
	value_type ubound() const
	{	return ubound_; }

private:
	DVectorSlice::const_iterator
		bl_itr_, bl_end_, bu_itr_, bu_end_;
	value_type
		big_bnd_, lbound_, ubound_;
	index_type
		index_;
	EBound
		at_bound_;

	void update() {
		if( bl_itr_ == bl_end_ ) {
			return;
		}
		else if( -big_bnd_ < *bl_itr_ && *bu_itr_ < +big_bnd_ ) {
			lbound_ = *bl_itr_;
			ubound_ = *bu_itr_;
			at_bound_ = BOTH;
		}
		else if( -big_bnd_ < *bl_itr_ ) {
			lbound_ = *bl_itr_;
			ubound_ = +big_bnd_;
			at_bound_ = LOWER;
		}
		else if( *bu_itr_ < +big_bnd_ ) {
			lbound_ = -big_bnd_;
			ubound_ = *bu_itr_;
			at_bound_ = UPPER;
		}
	}

	// not defined and not to be called
	sparse_bounds_itr();

};	// end class sparse_bounds_itr

}	// end namespace AbstractLinAlgPack

#endif // SPARSE_LIN_ALG_PACK_SPARSE_BOUNDS_H
