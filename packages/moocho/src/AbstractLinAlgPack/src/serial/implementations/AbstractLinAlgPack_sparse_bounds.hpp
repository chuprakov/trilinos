// ///////////////////////////////////////////////////////////////
// sparse_bounds.h

#include "SpVectorClass.h"

namespace SparseLinAlgPack {

///
/** Count the number of sparse bounds where at least one bound is
  * finite.
  */
size_type num_bounds( const SpVectorSlice& bl, const SpVectorSlice& bu );

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
	typedef	SpVectorSlice::element_type::indice_type	indice_type;
	///
	typedef	SpVectorSlice::element_type::value_type		value_type;


	///
	sparse_bounds_itr(	
			  const SpVectorSlice::const_iterator& bl_begin
			, const SpVectorSlice::const_iterator& bl_end
			, SpVectorSlice::difference_type bl_offset
			, const SpVectorSlice::const_iterator& bu_begin
			, const SpVectorSlice::const_iterator& bu_end
			, SpVectorSlice::difference_type bu_offset
			, value_type big_bnd = std::numeric_limits<value_type>::max()		)
		: bl_itr_(bl_begin), bl_end_(bl_end), bu_itr_(bu_begin), bu_end_(bu_end)
			, bl_offset_(bl_offset), bu_offset_(bu_offset)
			, big_bnd_(big_bnd)
	{	if(!at_end()) update(); }

	///
	bool at_end() const
	{	return bl_itr_ == bl_end_ && bu_itr_ == bu_end_; }

	///
	sparse_bounds_itr& operator++() {
		switch(at_bound_) {
		case LOWER:
			++bl_itr_;
			break;
		case UPPER:
			++bu_itr_;
			break;
		case BOTH:
			++bl_itr_;
			++bu_itr_; 
			break;
		}
		update();
		return *this;
	}

	///
	indice_type indice() const
	{	return indice_; }

	///
	value_type lbound() const
	{	return lbound_; }

	///
	value_type ubound() const
	{	return ubound_; }

private:
	SpVectorSlice::const_iterator
		bl_itr_, bl_end_, bu_itr_, bu_end_;
	SpVectorSlice::difference_type
		bl_offset_, bu_offset_;
	value_type
		big_bnd_, lbound_, ubound_;
	indice_type
		indice_;
	EBound
		at_bound_;

	void update() {
		if( bl_itr_ == bl_end_ && bu_itr_ == bu_end_ ) {
			return;
		}
		else if( ( bl_itr_ != bl_end_ ) 
			&& ( bu_itr_ == bu_end_ || bl_itr_->indice() + bl_offset_ < bu_itr_->indice() + bu_offset_ ) )
		{
			indice_ = bl_itr_->indice() + bl_offset_;
			lbound_ = bl_itr_->value();
			ubound_ = big_bnd_;
			at_bound_ = LOWER;
		}
		else if( ( bu_itr_ != bu_end_ )
			&& ( bl_itr_ == bl_end_ || bu_itr_->indice() + bu_offset_ < bl_itr_->indice() + bl_offset_ ) ) 
		{
			indice_ = bu_itr_->indice() + bu_offset_;
			lbound_ = - big_bnd_;
			ubound_ = bu_itr_->value();
			at_bound_ = UPPER;
		}
		else if(bl_itr_->indice() + bl_offset_ == bu_itr_->indice() + bu_offset_) {
			indice_ = bl_itr_->indice() + bl_offset_;
			lbound_ = bl_itr_->value();
			ubound_ = bu_itr_->value();
			at_bound_ = BOTH;
		}
	}

	// not defined and not to be called
	sparse_bounds_itr();

};	// end class sparse_bounds_itr

}	// end namespace SparseLinAlgPack