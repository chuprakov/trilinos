#include "TSFSerialVector.h"
#include "TSFError.h"
#ifdef HAVE_RTOP
#include "RTOpPack/include/RTOp_apply_op_serial.h"
#include <typeinfo>
#endif

using namespace TSF;

TSFSerialVector::TSFSerialVector(const TSFVectorSpace& space)
	: TSFInCoreVector(space), x_(space.dim())
{}

#ifdef HAVE_RTOP

void TSFSerialVector::apply_reduction(
	const RTOp_RTOp &op, int num_vecs, const TSFVectorBase* vecs[]
	,int num_targ_vecs, TSFVectorBase* targ_vecs[], RTOp_ReductTarget reduct_obj
	,const RTOp_index_type first_ele, const RTOp_index_type sub_dim, const RTOp_index_type global_offset
	) const
{
	this->apply_op(
		this,NULL,op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele,sub_dim,global_offset );
}

void TSFSerialVector::apply_transformation(
	const RTOp_RTOp &op, int num_vecs, const TSFVectorBase* vecs[]
	,int num_targ_vecs, TSFVectorBase* targ_vecs[], RTOp_ReductTarget reduct_obj
	,const RTOp_index_type first_ele, const RTOp_index_type sub_dim, const RTOp_index_type global_offset
	)
{
	this->apply_op(
		NULL,this,op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele,sub_dim,global_offset );
}

void TSFSerialVector::apply_op(
	const TSFSerialVector* const_this, TSFSerialVector* nonconst_this
	,const RTOp_RTOp& op
	,const int num_vecs_in,      const TSFVectorBase**   vecs
	,const int num_targ_vecs_in, TSFVectorBase**         targ_vecs
	,RTOp_ReductTarget reduct_obj
	,const int first_ele  , const int sub_dim  , const int global_offset
	) const
{
	typedef const RTOp_value_type*   RTOp_value_type_const_ptr_t;
	typedef  RTOp_value_type*        RTOp_value_type_ptr_t;
	int                           num_vecs           = num_vecs_in      + ( const_this      ? 1 : 0 );
	int                           num_targ_vecs      = num_targ_vecs_in + ( nonconst_this ? 1 : 0 );
	RTOp_value_type_const_ptr_t   *vec_ptrs          = NULL;
	ptrdiff_t                     *vec_strides       = NULL;
	RTOp_value_type_ptr_t         *targ_vec_ptrs     = NULL;
	ptrdiff_t                     *targ_vec_strides  = NULL;
	int k1, k2;
	// Fill up arrays of pointers to vector data
	if( num_vecs ) {
		vec_ptrs    = new RTOp_value_type_const_ptr_t[num_vecs];
		vec_strides = new ptrdiff_t[num_vecs];
		if(const_this) {
			vec_ptrs[0]    = &dynamic_cast<const TSFSerialVector*>(const_this)->x_[0];
			vec_strides[0] = 1;
		}
		for( k1 = 0, k2 = (const_this ? +1: 0); k1 < num_vecs_in; ++k1, ++k2 ) {
			vec_ptrs[k2]    = &dynamic_cast<const TSFSerialVector*>(vecs[k1])->x_[0];
			vec_strides[k2] = 1;
		}
	}
	if( num_targ_vecs ) {
		targ_vec_ptrs    = new RTOp_value_type_ptr_t[num_targ_vecs];
		targ_vec_strides = new ptrdiff_t[num_targ_vecs];
		if(nonconst_this) {
			targ_vec_ptrs[0]    = &dynamic_cast<TSFSerialVector*>(nonconst_this)->x_[0];
			targ_vec_strides[0] = 1;
		}
		for( k1 = 0, k2 = (nonconst_this ? +1: 0); k1 < num_targ_vecs_in; ++k1, ++k2 ) {
			targ_vec_ptrs[k2]    = &dynamic_cast<TSFSerialVector*>(targ_vecs[k1])->x_[0];
			targ_vec_strides[k2] = 1;
		}
	}
	// Call the helper function to setup the sub-vectors for the logical vector
	// specified by (first_ele,sub_dim,global_offset) and then perform the
	// reduction/transformation operation.
	const int err = RTOp_apply_op_serial(
		this->x_.length(), num_vecs, vec_ptrs, vec_strides
		,num_targ_vecs, targ_vec_ptrs, targ_vec_strides
		,first_ele, sub_dim, global_offset, &op, reduct_obj
		);
	assert( err == 0 ); // Should throw an exception or something if this fails!
	// Delete dynamically allocated memory
	if(num_vecs)      { delete [] vec_ptrs;       delete [] vec_strides;      }
	if(num_targ_vecs) {	delete [] targ_vec_ptrs;  delete [] targ_vec_strides; }
}

#endif

void TSFSerialVector::setElements(int n, const int* globalIndices,
																	const TSFReal* values) 
{
	for (int i=0; i<n; i++)
		{
			x_[globalIndices[i]] = values[i];
		}
}

void TSFSerialVector::getElements(int n, const int* globalIndices,
																	TSFReal* values) const 
{
	for (int i=0; i<n; i++)
		{
			values[i] = x_[globalIndices[i]];
		}
}

void TSFSerialVector::addToElements(int n, const int* globalIndices,
																		const TSFReal* values) 
{
	for (int i=0; i<n; i++)
		{
			x_[globalIndices[i]] += values[i];
		}
}

TSFVectorBase* TSFSerialVector::deepCopy() const 
{
	TSFVectorBase* rtn = new TSFSerialVector(*this);
	if (rtn==0) TSFError::raise("TSFSerialVector::deepCopy()");
	return rtn;
}

ostream& TSFSerialVector::print(ostream& os) const 
{
	os << x_ ;
	return os;
}

const DenseSerialVector& TSFSerialVector::getConcrete(const TSFVector& x)
{
	const TSFSerialVector* v = dynamic_cast<const TSFSerialVector*>(x.ptr());
	if (v==0) {
	    std::cerr << "bad vector is ";
	    x.print(cerr);
	    cerr << endl;
#ifdef HAVE_RTOP
		cerr << "With the type \'" << typeid(*x.ptr()).name() << endl;
#endif
	    TSFError::raise("bad cast in TSFSerialVector::getConcrete");
	}
	return v->x_;
}

DenseSerialVector& TSFSerialVector::getConcrete(TSFVector& x)
{
	TSFSerialVector* v = dynamic_cast<TSFSerialVector*>(x.ptr());
	if (v==0) {
	    std::cerr << "bad vector is ";
	    x.print(cerr);
	    cerr << endl;
#ifdef HAVE_RTOP
		cerr << "With the type \'" << typeid(*x.ptr()).name() << endl;
#endif
	    TSFError::raise("bad cast in TSFSerialVector::getConcrete");
	}
	return v->x_;
}










