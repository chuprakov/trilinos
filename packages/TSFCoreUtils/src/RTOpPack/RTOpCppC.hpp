// /////////////////////////////////////////////////////////////////////
// RTOpCppC.hpp
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

#ifndef RTOP_CPP_C_H
#define RTOP_CPP_C_H

#include "RTOpCpp.hpp"

namespace RTOpPack {

///
/** C++ subclass for RTOp_RTOp objects.
 *
 * ToDo: Finish Documentation
 */
class RTOpC : public RTOp {
public:

	///
	RTOpC( const RTOp_RTOp_vtbl_t* = NULL );
	///
	~RTOpC();
	///
	RTOp_RTOp& op();
	///
	const RTOp_RTOp& op() const;

	/** @name Overridden from RTOp */
	//@{

	///
	op_create_func_t get_op_create_func() const;
	///
	op_free_func_t get_op_free_func() const;
	///
	void get_op_type_num_entries(
		int*  num_values
		,int* num_indexes
		,int* num_chars
		) const;
	///
	void extract_op_state(
		int               num_values
		,RTOp_value_type  value_data[]
		,int              num_indexes
		,RTOp_index_type  index_data[]
		,int              num_chars
		,RTOp_char_type   char_data[]
		) const;
	///
	void load_op_state(
		int                       num_values
		,const RTOp_value_type    value_data[]
		,int                      num_indexes
		,const RTOp_index_type    index_data[]
		,int                      num_chars
		,const RTOp_char_type     char_data[]
		);
	///
	void get_reduct_type_num_entries(
		int*   num_values
		,int*  num_indexes
		,int*  num_chars
		) const;
	///
	void reduct_obj_create_raw( RTOp_ReductTarget* reduct_obj ) const;
	///
	void reduct_obj_reinit( RTOp_ReductTarget reduct_obj ) const;
	///
	void reduct_obj_free( RTOp_ReductTarget* reduct_obj ) const;
	///
	void extract_reduct_obj_state(
		const RTOp_ReductTarget   reduct_obj
		,int                      num_values
		,RTOp_value_type          value_data[]
		,int                      num_indexes
		,RTOp_index_type          index_data[]
		,int                      num_chars
		,RTOp_char_type           char_data[]
		) const;
	///
	void load_reduct_obj_state(
		int                      num_values
		,const RTOp_value_type   value_data[]
		,int                     num_indexes
		,const RTOp_index_type   index_data[]
		,int                     num_chars
		,const RTOp_char_type    char_data[]
		,RTOp_ReductTarget       reduct_obj
		) const;
	///
	void apply_op(
		const int   num_vecs,       const SubVector         sub_vecs[]
		,const int  num_targ_vecs,  const MutableSubVector  targ_sub_vecs[]
		,RTOp_ReductTarget reduct_obj
		) const;
	///
	void reduce_reduct_objs(
		RTOp_ReductTarget in_reduct_obj, RTOp_ReductTarget inout_reduct_obj
		) const;
	///
    void get_reduct_op( RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr ) const;

	//@}

private:

	const RTOp_RTOp_vtbl_t  *vtbl_;
	RTOp_RTOp               op_;

}; // end class RTOpC

// //////////////////////////////////////////
// Inline member functions

// Inline members for RTOpC

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

} // end namespace RTOpPack

#endif // RTOP_CPP_C_H
