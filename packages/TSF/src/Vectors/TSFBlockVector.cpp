#include "TSFBlockVector.h"
#include "TSFVectorSpace.h"
#include "TSFVectorSpaceBase.h"

#include "TSFProductSpace.h"
#include "TSFError.h"
#include "TSFUtils.h"

#ifdef HAVE_RTOP
#include "RTOpPack/include/RTOp_parallel_helpers.h"
#endif

using namespace TSF;
using namespace std;

#ifdef HAVE_RTOP

void TSFBlockVector::apply_reduction(
	const RTOp_RTOp &op, int num_vecs, const TSFVectorBase* vecs_in[]
	,int num_targ_vecs, TSFVectorBase* targ_vecs_in[], RTOp_ReductTarget reduct_obj
	,const RTOp_index_type first_ele, const RTOp_index_type sub_dim, const RTOp_index_type global_offset
	) const
{
	// Convert vectors to TSFBlockVector and append this to the beginning of vecs[]
	typedef const TSFBlockVector*   const_vec_base_ptr_t;
	typedef TSFBlockVector*         vec_base_ptr_t;
	const_vec_base_ptr_t    *vecs      = new const_vec_base_ptr_t[num_vecs+1];
	vec_base_ptr_t          *targ_vecs = NULL;
	vecs[0] = this;
	for(int k = 1; k <= num_vecs; ++k) {
		vecs[k] = dynamic_cast<const TSFBlockVector*>(vecs_in[k-1]);
	}
	if(num_targ_vecs) {
		targ_vecs = new vec_base_ptr_t[num_targ_vecs];
		for(int k = 0; k < num_targ_vecs; ++k) {
			targ_vecs[k] = dynamic_cast<TSFBlockVector*>(targ_vecs_in[k]);
		}
	}
	apply_op(op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele,sub_dim,global_offset);
}

void TSFBlockVector::apply_transformation(
	const RTOp_RTOp &op, int num_vecs, const TSFVectorBase* vecs_in[]
	,int num_targ_vecs, TSFVectorBase* targ_vecs_in[], RTOp_ReductTarget reduct_obj
	,const RTOp_index_type first_ele, const RTOp_index_type sub_dim, const RTOp_index_type global_offset
	)
{
	// Convert vectors to TSFBlockVector and append this to the beginning of targ_vecs[]
	typedef const TSFBlockVector*   const_vec_base_ptr_t;
	typedef TSFBlockVector*         vec_base_ptr_t;
	const_vec_base_ptr_t    *vecs      = NULL;
	vec_base_ptr_t          *targ_vecs = new vec_base_ptr_t[num_targ_vecs+1];
	if(num_vecs) {
		for(int k = 0; k < num_vecs; ++k) {
			vecs[k] = dynamic_cast<const TSFBlockVector*>(vecs_in[k]);
		}
	}
	targ_vecs[0] = this;
	for(int k = 1; k <= num_targ_vecs; ++k) {
		targ_vecs[k] = dynamic_cast<TSFBlockVector*>(targ_vecs_in[k-1]);
	}
	apply_op(op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele,sub_dim,global_offset);
}

void TSFBlockVector::apply_op(
	const RTOp_RTOp &op, int num_vecs, const TSFBlockVector* vecs[]
	,int num_targ_vecs, TSFBlockVector* targ_vecs[], RTOp_ReductTarget reduct_obj
	,const RTOp_index_type first_ele, const RTOp_index_type sub_dim, const RTOp_index_type global_offset
	) const
{
	// Storage of pointers to constituent TSFVector subvectors.
	typedef const TSFVector*     const_vec_ptr_t;
	typedef TSFVector*           vec_ptr_t;
	const_vec_ptr_t              *sub_vecs      = NULL;
	vec_ptr_t                    *targ_sub_vecs = NULL;
	if( num_vecs ) {
		sub_vecs = new const_vec_ptr_t[num_vecs];
	}
	if( num_targ_vecs ) {
		targ_sub_vecs = new vec_ptr_t[num_targ_vecs];
	}
	// Loop through the subvectors and apply the operator to regions of overlap
	// Here we can use the utility function for parallel vectors since the idea
	// is very similar.
	RTOp_index_type    global_dim = this->space().dim();
	RTOp_index_type
		local_sub_dim = 0,
		local_offset  = 0,
		overlap_first_local_ele  = 0,
		overalap_local_sub_dim   = 0,
		overlap_global_offset    = 0;
	const int numBlocks = this->numBlocks();
	for( int k = 0; k < numBlocks; ++k ) { // Loop through the subvector blocks
		const TSFVector &curr_sv = subvectors_[k];
		local_sub_dim = curr_sv.space().dim();
		RTOp_parallel_calc_overlap(
			global_dim, local_sub_dim, local_offset, first_ele, sub_dim, global_offset
			,&overlap_first_local_ele, &overalap_local_sub_dim, &overlap_global_offset
			);
		if( overlap_first_local_ele ) {
			// This block subvector overlaps the logical subvector requested by the client.
			// Fill up the subvector arguments
			if( num_vecs ) {
				for( int j = 0; j < num_vecs; ++j )
					sub_vecs[j] = &vecs[j]->subvectors_[k];
			}
			if( num_targ_vecs ) {
				for( int j = 0; j < num_targ_vecs; ++j )
					targ_sub_vecs[j] = &targ_vecs[j]->subvectors_[k];
			}
			// Apply the operator and accumulate the result of the reduction in reduct_obj
			if( num_vecs )
				sub_vecs[0]->apply_reduction(
					op
					,num_vecs-1, num_vecs > 1 ? &sub_vecs[1] : NULL
					,num_targ_vecs, targ_sub_vecs
					,reduct_obj
					,overlap_first_local_ele, overalap_local_sub_dim, overlap_global_offset
					);
			else
				targ_sub_vecs[0]->apply_transformation(
					op
					,num_vecs, sub_vecs
					,num_targ_vecs-1, num_targ_vecs > 1 ? &targ_sub_vecs[1] : NULL
					,reduct_obj
					,overlap_first_local_ele, overalap_local_sub_dim, overlap_global_offset
					);
		}
		//
		local_offset += local_sub_dim;
	}
	// Delete arrays
	if( sub_vecs ) delete [] sub_vecs;
	if( targ_sub_vecs ) delete [] targ_sub_vecs;
}

#endif

TSFBlockVector::TSFBlockVector(const TSFVectorSpace& space)
	: TSFVectorBase(space), subvectors_(space.numBlocks())
{
	for (int i=0; i<space.numBlocks(); i++)
		{
			subvectors_[i] = space.getBlock(i).createMember();
		}
}

int TSFBlockVector::numBlocks() const 
{
	return subvectors_.size();
}

void TSFBlockVector::getBlock(int i, const TSFVector& /* self */, 
																TSFVector& sub) const 
{
	sub = subvectors_[i];
}

void TSFBlockVector::setBlock(int i, const TSFVector& sub)
{
	subvectors_[i] = sub;
}

TSFReal& TSFBlockVector::setElement(int /* g */)
{
	TSFError::raise("TSFBlockVector::setElement should not be called. Set the"
									"elements in each block");
	return dummyElement_;
}

const TSFReal& TSFBlockVector::getElement(int /* g */) const 
{
	TSFError::raise("TSFBlockVector::getElement should not be called. Get the"
									"elements from each block");
	return dummyElement_; // -Wall
}

void TSFBlockVector::setElements(int /* n */, const int* /* globalIndices */,
																 const TSFReal* /* values */) 
{
	TSFError::raise("TSFBlockVector::setElements should not be called.");
}

void TSFBlockVector::getElements(int /* n */, const int* /* globalIndices */,
																 TSFReal* /* values */) const 
{
	TSFError::raise("TSFBlockVector::getElements should not be called. Get the"
									"elements from each block");
}

void TSFBlockVector::addToElements(int n, const int* /* globalIndices */,
																	 const TSFReal* /* values */) 
{	
	TSFError::raise("TSFBlockVector::addToElements should not be called. Get the"
									"elements from each block");
}

void TSFBlockVector::axpy(const TSFReal& a, const TSFVector& other) 
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].selfModifyingAxpy(a, other.getBlock(i));
		}
}

void TSFBlockVector::acceptCopyOf(const TSFVector& other) 
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].acceptCopyOf(other);
		}
}

void TSFBlockVector::scalarMult(const TSFReal& a)
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].selfModifyingScalarMult(a);
		}
}

void TSFBlockVector::dotStar(const TSFVector& other) 
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].selfModifyingDotStar(other.getBlock(i));
		}
}

void TSFBlockVector::dotSlash(const TSFVector& other) 
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].selfModifyingDotSlash(other.getBlock(i));
		}
}


TSFReal TSFBlockVector::dot(const TSFVector& other) const
{
	TSFReal sum = 0.0;

	for (int i=0; i<numBlocks(); i++)
		{
			sum += subvectors_[i].dot(other.getBlock(i));
		}
	return sum;
}

TSFReal TSFBlockVector::norm1() const
{
	TSFReal sum = 0.0;

	for (int i=0; i<numBlocks(); i++)
		{
			sum += subvectors_[i].norm1();
		}
	return sum;
}

TSFReal TSFBlockVector::normInf() const 
{
	TSFReal biggest = TSFUtils::negativeInfinity();
	for (int i=0; i<numBlocks(); i++)
		{
			TSFReal x = subvectors_[i].normInf();
			if (biggest < x) biggest = x;
		}
	return biggest;
}

TSFReal TSFBlockVector::sumElements() const
{
	TSFReal sum = 0.0;

	for (int i=0; i<numBlocks(); i++)
		{
			sum += subvectors_[i].sumElements();
		}
	return sum;
}

void TSFBlockVector::setScalar(const TSFReal& a) 
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].setScalar(a);
		}
}

TSFReal TSFBlockVector::findExtremeValue(MinOrMax type, TSFGeneralizedIndex& location, 
																				 const TSFReal& tol) const
{
	TSFReal current;
	if (type==MIN) current = TSFUtils::infinity();
	else current = TSFUtils::negativeInfinity();

	for (int i=0; i<numBlocks(); i++)
		{
			TSFGeneralizedIndex j;
			TSFReal blockVal = subvectors_[i].findExtremeValue(type, j, tol);
			if (type==MIN && blockVal < current)
				{
					location = TSFGeneralizedIndex(j, i);
					current = blockVal;
				}
			if (type==MAX && blockVal > current)
				{
					location = TSFGeneralizedIndex(j, i);
					current = blockVal;
				}
		}
	return current;
}

void TSFBlockVector::randomize(const TSFRandomNumberGenerator& r)
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].randomize(r);
		}
}

void TSFBlockVector::abs()
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].abs();
		}
}

TSFVectorBase* TSFBlockVector::deepCopy() const 
{
	TSFVectorBase* rtn = new TSFBlockVector(space());
	for (int i=0; i<numBlocks(); i++)
		{
			rtn->setBlock(i, subvectors_[i].copy());
		}
	return rtn;
}

void TSFBlockVector::synchronizeGhostValues() const
{
	for (int i=0; i<numBlocks(); i++) subvectors_[i].synchronizeGhostValues();
}

void TSFBlockVector::invalidateGhostValues() 
{
	for (int i=0; i<numBlocks(); i++) subvectors_[i].invalidateGhostValues();
}

ostream& TSFBlockVector::print(ostream& os) const 
{
	for (int i=0; i<numBlocks(); i++) os << subvectors_[i] << endl;
	return os;
}

