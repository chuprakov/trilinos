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

// ////////////////////////////////
// RTOpPack_RTOpTHelpers.hpp

#ifndef RTOPPACK_RTOP_NEW_T_HELPERS_HPP
#define RTOPPACK_RTOP_NEW_T_HELPERS_HPP

#include "RTOpPack_RTOpT.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "dynamic_cast_verbose.hpp"

namespace RTOpPack {

///
/** Simple <tt>ReductTarget</tt> subclass for simple scalar objects
 */
template<class Scalar>
class ReductTargetScalar : public ReductTarget {
public:
  ReductTargetScalar( const Scalar &scalar = Scalar(0.0) ) : scalar_(scalar) {}
  void set( const Scalar &scalar ) { scalar_ = scalar; }
  const Scalar& get() const { return scalar_; }
private:
  Scalar scalar_;
};

///
/** Simple class for all reduction operators that return a simple
 * scalar reduction object.
 *
 * Subclasses have to minimally define <tt>op_name()</tt> and
 * <tt>apply_op()</tt>.  Subclasses should also override
 * <tt>reduce_reduct_objs()</tt> if the reduction is not a simple
 * summation.
 */
template<class Scalar>
class ROpScalarReductionBase : public RTOpT<Scalar> {
public:
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  ///
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, initReductObjValue );
  ///
  ROpScalarReductionBase( const Scalar &initReductObjValue = static_cast<const Scalar&>(0.0) )
    : initReductObjValue_(initReductObjValue)
    {}
  ///
  const Scalar& getRawVal( const ReductTarget &reduct_obj ) const
    {
      using DynamicCastHelperPack::dyn_cast;
      return dyn_cast<const ReductTargetScalar<Scalar> >(reduct_obj).get();
    }
  ///
  void setRawVal( const Scalar &rawVal, ReductTarget *reduct_obj ) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION( reduct_obj==NULL, std::invalid_argument, "Error!" );
#endif
      using DynamicCastHelperPack::dyn_cast;
      dyn_cast<ReductTargetScalar<Scalar> >(*reduct_obj).set(rawVal);
    }
  /** @name Overridden from RTOpT */
  //@{
  ///
	void get_reduct_type_num_entries(
		int*   num_values
		,int*  num_indexes
		,int*  num_chars
		) const
    {
      const Scalar dummy(0.0);
      *num_values = Teuchos::PrimitiveTypeTraits<Scalar>::numPrimitiveObjs(dummy);
      *num_indexes = 0;
      *num_chars = 0;
    }
	///
	Teuchos::RefCountPtr<ReductTarget> reduct_obj_create() const
    {
      return Teuchos::rcp(new ReductTargetScalar<Scalar>(initReductObjValue()));
    }
	/// Default implementation here is for a sum
	void reduce_reduct_objs(
		const ReductTarget& _in_reduct_obj, ReductTarget* _inout_reduct_obj
		) const
    {
      using DynamicCastHelperPack::dyn_cast;
      const ReductTargetScalar<Scalar> &in_reduct_obj    = dyn_cast<const ReductTargetScalar<Scalar> >(_in_reduct_obj); 
      ReductTargetScalar<Scalar>       &inout_reduct_obj = dyn_cast<ReductTargetScalar<Scalar> >(*_inout_reduct_obj); 
      inout_reduct_obj.set( inout_reduct_obj.get() + in_reduct_obj.get() );
    }
	///
	void reduct_obj_reinit( ReductTarget* reduct_obj ) const
    {
      setRawVal( initReductObjValue(), reduct_obj );
    }
	///
	void extract_reduct_obj_state(
		const ReductTarget     &reduct_obj
		,int                      num_values
		,primitive_value_type     value_data[]
		,int                      num_indexes
		,RTOp_index_type          index_data[]
		,int                      num_chars
		,RTOp_char_type           char_data[]
		) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION(
        num_values==0 || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Teuchos::PrimitiveTypeTraits<Scalar>::extractPrimitiveObjs( getRawVal(reduct_obj), num_values, value_data );
    }
	///
	void load_reduct_obj_state(
		int                            num_values
		,const primitive_value_type    value_data[]
		,int                           num_indexes
		,const RTOp_index_type         index_data[]
		,int                           num_chars
		,const RTOp_char_type          char_data[]
		,ReductTarget               *reduct_obj
		) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION(
        num_values==0 || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Scalar val(0.0);
      Teuchos::PrimitiveTypeTraits<Scalar>::loadPrimitiveObjs( num_values, value_data, &val );
      setRawVal( val, reduct_obj );
    }
  //@}
}; // class ROpScalarReductionBase

} // namespace RTOpPack

///
/** Use within an apply_op(...) function implemention here num_vecs==1, num_targ_vecs==0.
 */
#define RTOP_APPLY_OP_1_0( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=1 || (SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==1, sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=0 || (TARG_SUB_VECS)!=NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
    ); \
  const RTOpPack::index_type   subDim  = (SUB_VECS)[0].subDim(); \
  const Scalar                 *v0_val = (SUB_VECS)[0].values(); \
  const ptrdiff_t              v0_s    = (SUB_VECS)[0].stride()

///
/** Use within an apply_op(...) function implemention here num_vecs==0, num_targ_vecs==1.
 */
#define RTOP_APPLY_OP_0_1( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=0 || (SUB_VECS)!=NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==0, sub_vecs==NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=1 || (TARG_SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==1, targ_sub_vecs!=NULL" \
    ); \
  const RTOpPack::index_type   subDim  = (TARG_SUB_VECS)[0].subDim(); \
  Scalar                       *z0_val = (TARG_SUB_VECS)[0].values(); \
  const ptrdiff_t              z0_s    = (TARG_SUB_VECS)[0].stride()

#endif // RTOPPACK_RTOP_NEW_T_HELPERS_HPP
