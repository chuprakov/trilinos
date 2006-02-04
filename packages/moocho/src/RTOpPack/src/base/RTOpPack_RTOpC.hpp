// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// ///////////////////////////////
// RTOpPack_RTOpC.hpp

#ifndef RTOPPACK_RTOP_NEW_C_HPP
#define RTOPPACK_RTOP_NEW_C_HPP

#include "RTOpPack_OldTypes.hpp"
#include "RTOpPack_RTOpT.hpp"
#include "RTOp.h"
#include "Teuchos_dyn_cast.hpp"

namespace RTOpPack {

///
/** Adapter subclass that uses a <tt>RTOp_RTOp</tt> object.
 *
 * ToDo: Finish documentation!
 */
class RTOpC : public RTOpT<RTOp_value_type> {
public:

  ///
  typedef RTOp_value_type Scalar;
	///
	RTOpC();
	///
	~RTOpC();
	///
	RTOp_RTOp& op();
	///
	const RTOp_RTOp& op() const;
  ///
  RTOp_ReductTarget& operator()(ReductTarget& reduct_obj) const;
  ///
  const RTOp_ReductTarget& operator()(const ReductTarget& reduct_obj) const;

  /** @name Overridden from RTOpT */
  //@{

  ///
	void get_reduct_type_num_entries(
		int*   num_values
		,int*  num_indexes
		,int*  num_chars
		) const;
	///
	Teuchos::RefCountPtr<ReductTarget> reduct_obj_create() const;
	///
	void reduce_reduct_objs(
		const ReductTarget& _in_reduct_obj, ReductTarget* _inout_reduct_obj
		) const;
	///
	void reduct_obj_reinit( ReductTarget* reduct_obj ) const;
	///
	void extract_reduct_obj_state(
		const ReductTarget     &reduct_obj
		,int                      num_values
		,primitive_value_type     value_data[]
		,int                      num_indexes
		,index_type               index_data[]
		,int                      num_chars
		,char_type                char_data[]
		) const;
	///
	void load_reduct_obj_state(
		int                            num_values
		,const primitive_value_type    value_data[]
		,int                           num_indexes
		,const index_type              index_data[]
		,int                           num_chars
		,const char_type               char_data[]
		,ReductTarget               *reduct_obj
		) const;
  ///
	void get_op_type_num_entries(
		int*  num_values
		,int* num_indexes
		,int* num_chars
		) const;
	///
	void extract_op_state(
		int                             num_values
		,primitive_value_type           value_data[]
		,int                            num_indexes
		,index_type                     index_data[]
		,int                            num_chars
		,char_type                      char_data[]
		) const;
	///
	void load_op_state(
		int                           num_values
		,const primitive_value_type   value_data[]
		,int                          num_indexes
		,const index_type             index_data[]
		,int                          num_chars
		,const char_type              char_data[]
		);
	///
	bool coord_invariant() const;
  ///
  const char* op_name() const;
  ///
	void apply_op(
		const int   num_vecs,       const SubVectorT<Scalar>         sub_vecs[]
		,const int  num_targ_vecs,  const MutableSubVectorT<Scalar>  targ_sub_vecs[]
		,ReductTarget *_reduct_obj
		) const;

  //@}

private:

	RTOp_RTOp op_;

}; // class RTOpC

///
/** Adapter subclass for <tt>RTOp_ReductTarget</tt>
 */
class ReductTargetC : public ReductTarget {
public:
  inline ReductTargetC( const RTOp_RTOp& op, RTOp_ReductTarget obj );
  inline ~ReductTargetC();
  inline RTOp_ReductTarget& obj();
  inline const RTOp_ReductTarget& obj() const;
private:
  const RTOp_RTOp      &op_;
  RTOp_ReductTarget    obj_;
  ReductTargetC(); // Not defined and not to be called
};

// ///////////////////////////////
// Inline member functions

// RTOpC

inline
RTOp_RTOp& RTOpC::op()
{
	return op_;
}

inline
const RTOp_RTOp& RTOpC::op() const
{
	return op_;
}

inline
RTOp_ReductTarget&
RTOpC::operator()(ReductTarget& reduct_obj) const
{
  return Teuchos::dyn_cast<ReductTargetC>(reduct_obj).obj();
}

inline
const RTOp_ReductTarget&
RTOpC::operator()(const ReductTarget& reduct_obj) const
{
  return Teuchos::dyn_cast<const ReductTargetC>(reduct_obj).obj();
}

// ReductTargetC

inline
ReductTargetC::ReductTargetC( const RTOp_RTOp& op, RTOp_ReductTarget obj )
  : op_(op), obj_(obj)
{} 

inline
ReductTargetC::~ReductTargetC()
{
  if( obj() != RTOp_REDUCT_OBJ_NULL ) {
    TEST_FOR_EXCEPTION(
      0!=RTOp_reduct_obj_free(&op_,&obj_)
      ,UnknownError
      ,"RTOpC::reduct_obj_free(...): Error, "
      "RTOp_reduct_obj_free(...) returned != 0"
      );
  }
} 

inline
RTOp_ReductTarget& ReductTargetC::obj()
{
  return obj_;
}

inline
const RTOp_ReductTarget& ReductTargetC::obj() const
{
  return obj_;
}

} // namespace RTOpPack

#endif // RTOPPACK_RTOP_NEW_C_HPP
