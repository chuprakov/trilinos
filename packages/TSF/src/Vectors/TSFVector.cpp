#include "TSFVector.h"
#include "TSFVectorBase.h"
#include "TSFError.h"
#include "TSFDeferredLinearCombination.h"

#include "TSFOut.h"
#include "DenseSerialVector.h"

using namespace TSF;

TSFVector::TSFVector()
	: ptr_(0)
{;}

TSFVector::TSFVector(TSFVectorBase* ptr)
	: ptr_(ptr)
{;}

TSFVector::TSFVector(const TSFDeferredLinearCombination& op)
	: ptr_(0)
{
	TSFVector tmp = op.evaluateIntoNewVector();
	ptr_ = tmp.ptr_;
}

TSFVector::TSFVector(const TSFDeferredCopy& copy)
	: ptr_(0)
{
	TSFVector tmp;
	copy.createCopy(tmp);
	ptr_ = tmp.ptr_;
}

TSFVector& TSFVector::operator=(const TSFDeferredLinearCombination& op)
{
	if (ptr_.isNull()) 
		{
			TSFVector tmp = op.evaluateIntoNewVector();
			ptr_ = tmp.ptr_;
		}
	else
		{
			op.evaluateIntoExistingVector(*this);
		}
	return *this;
}

TSFVector& TSFVector::operator=(const TSFDeferredCopy& copy)
{
	if (ptr_.isNull()) 
		{
			TSFVector tmp;
			copy.createCopy(tmp);
			ptr_ = tmp.ptr_;
		}
	else
		{
			copy.copyInto(*this);
		}
	return *this;
}

bool TSFVector::isNull() const 
{
	return ptr_.isNull();
}

bool TSFVector::isIdenticalTo(const TSFVector& other) const
{
	return &(*ptr_) == &(*(other.ptr_));
}

const TSFVectorSpace& TSFVector::space() const
{
	return ptr_->space();
}

int TSFVector::numBlocks() const 
{
	return ptr_->numBlocks();
}

TSFVector TSFVector::getBlock(int i) const 
{
	TSFVector rtn;
	if (i<0 || i>=numBlocks()) 
		{
			TSFError::raise("range error in TSFVector::getBlock()");
		}
	ptr_->getBlock(i, *this, rtn);
	return rtn;
}

TSFVector& TSFVector::setBlock(int i, const TSFVector& sub)
{
	if (i<0 || i>=numBlocks()) 
		{
			TSFError::raise("range error in TSFVector::setBlock()");
		}
	ptr_->setBlock(i, sub);
	return *this;
}

#if HAVE_RTOP

void TSFVector::apply_reduction(
	const RTOp_RTOp &op, int num_vecs, const TSFVector* vecs_in[]
	,int num_targ_vecs, TSFVector* targ_vecs_in[], RTOp_ReductTarget reduct_obj
	,const RTOp_index_type first_ele, const RTOp_index_type sub_dim, const RTOp_index_type global_offset
	) const
{
	// Convert from TSFVector arrays to TSFVectorBase arrays
	typedef const TSFVectorBase*        const_vec_base_ptr_t;
	typedef TSFVectorBase*              vec_base_ptr_t;
	const_vec_base_ptr_t      *vecs      = NULL;
	vec_base_ptr_t            *targ_vecs = NULL;
	if(num_vecs) {
		vecs = new const_vec_base_ptr_t[num_vecs];
		for( int k = 0; k < num_vecs; ++k )
			vecs[k] = vecs_in[k]->ptr();
	}
	if(num_targ_vecs) {
		targ_vecs = new vec_base_ptr_t[num_targ_vecs];
		for( int k = 0; k < num_targ_vecs; ++k )
			targ_vecs[k] = targ_vecs_in[k]->ptr();
	}
	// Call the implementation to invoke the operator
	ptr_->apply_reduction(
		op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj
		,first_ele,sub_dim,global_offset
		);
	// Delete the arrays
	if(vecs)      delete [] vecs;
	if(targ_vecs) delete [] targ_vecs;
}

void TSFVector::apply_transformation(
	const RTOp_RTOp &op, int num_vecs, const TSFVector* vecs_in[]
	,int num_targ_vecs, TSFVector* targ_vecs_in[], RTOp_ReductTarget reduct_obj
	,const RTOp_index_type first_ele, const RTOp_index_type sub_dim, const RTOp_index_type global_offset
	)
{
	// Convert from TSFVector arrays to TSFVectorBase arrays
	typedef const TSFVectorBase*        const_vec_base_ptr_t;
	typedef TSFVectorBase*              vec_base_ptr_t;
	const_vec_base_ptr_t      *vecs      = NULL;
	vec_base_ptr_t            *targ_vecs = NULL;
	if(num_vecs) {
		vecs = new const_vec_base_ptr_t[num_vecs];
		for( int k = 0; k < num_vecs; ++k )
			vecs[k] = vecs_in[k]->ptr();
	}
	if(num_targ_vecs) {
		targ_vecs = new vec_base_ptr_t[num_targ_vecs];
		for( int k = 0; k < num_targ_vecs; ++k )
			targ_vecs[k] = targ_vecs_in[k]->ptr();
	}
	// Call the implementation to invoke the operator
	ptr_->apply_transformation(
		op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj
		,first_ele,sub_dim,global_offset
		);
	// Delete the arrays
	if(vecs)      delete [] vecs;
	if(targ_vecs) delete [] targ_vecs;
}

#endif

const TSFReal& TSFVector::operator[](int i) const 
{
	return ptr_->getElement(i);
}

TSFReal& TSFVector::operator[](int i)
{
	return ptr_->setElement(i);
}


void TSFVector::getElements(const TSFArray<int>& globalIndices, 
														DenseSerialVector& sub) const 
{
	sub.resize(globalIndices.size());
	ptr_->getElements(sub.length(), &(globalIndices[0]), &(sub[0]));
}

void TSFVector::setElements(const TSFArray<int>& globalIndices, 
														const DenseSerialVector& sub)
{
	ptr_->setElements(sub.length(), &(globalIndices[0]), &(sub[0]));
}

void TSFVector::getElements(const int* globalIndices, int length,
														DenseSerialVector& sub) const 
{
	sub.resize(length);
	ptr_->getElements(length, globalIndices, &(sub[0]));
}

void TSFVector::setElements(const int* globalIndices, 
														const DenseSerialVector& sub)
{
	ptr_->setElements(sub.length(), globalIndices, &(sub[0]));
}

void TSFVector::addToElements(const TSFArray<int>& globalIndices, 
															const DenseSerialVector& sub)
{
	ptr_->addToElements(sub.length(), &(globalIndices[0]), &(sub[0]));
}

void TSFVector::addToElements(const int* globalIndices, 
															const DenseSerialVector& sub)
{
	ptr_->addToElements(sub.length(), globalIndices, &(sub[0]));
}


TSFDeferredCopy TSFVector::copy() const
{
	return ptr_;
}

TSFReal TSFVector::dot(const TSFVector& other) const
{
	TSFTimeMonitor t(opTimer());
	return ptr_->dot(other);
}

TSFReal TSFVector::norm2() const
{
	TSFTimeMonitor t(opTimer());
	return sqrt(dot(*this));
}

TSFReal TSFVector::norm1() const
{
	TSFTimeMonitor t(opTimer());
	return ptr_->norm1();
}

TSFReal TSFVector::normInf() const
{
	TSFTimeMonitor t(opTimer());
	return ptr_->normInf();
}

TSFReal TSFVector::max() const
{
	TSFGeneralizedIndex i;
	return max(i);
}

TSFReal TSFVector::min() const
{
	TSFGeneralizedIndex i;
	return min(i);
}

TSFReal TSFVector::min(TSFGeneralizedIndex& i) const
{
	return findExtremeValue(TSFVectorBase::MIN, i, TSFUtils::negativeInfinity());
}

TSFReal TSFVector::max(TSFGeneralizedIndex& i) const
{
	return findExtremeValue(TSFVectorBase::MAX, i, TSFUtils::infinity());
}

TSFReal TSFVector::max(const TSFReal& tol, TSFGeneralizedIndex& i) const
{
	return findExtremeValue(TSFVectorBase::MAX, i, tol);
}

TSFReal TSFVector::min(const TSFReal& tol, TSFGeneralizedIndex& i) const
{
	return findExtremeValue(TSFVectorBase::MIN, i, tol);
}

TSFReal TSFVector::findExtremeValue(TSFVectorBase::MinOrMax type,
																		TSFGeneralizedIndex& i,
																		const TSFReal& tol) const
{
	TSFTimeMonitor t(opTimer());
	return ptr_->findExtremeValue(type, i, tol);
}

TSFVector TSFVector::abs() const
{
	TSFVector rtn = copy();
	TSFTimeMonitor t(opTimer());
	rtn.ptr_->abs();
	return rtn;
}

TSFReal TSFVector::sumElements() const
{
	TSFTimeMonitor t(opTimer());
	return ptr_->sumElements();
}

void TSFVector::setScalar(const TSFReal& a)
{
	TSFTimeMonitor t(opTimer());
	return ptr_->setScalar(a);
}

void TSFVector::randomize(const TSFRandomNumberGenerator& r)
{
	ptr_->randomize(r);
}

void TSFVector::zero() 
{
	setScalar(0.0);
}

void TSFVector::print(ostream& os) const
{
	ptr_->print(os);
}

void TSFVector::synchronizeGhostValues() const
{
	ptr_->synchronizeGhostValues();
}

void TSFVector::invalidateGhostValues() 
{
	ptr_->invalidateGhostValues();
}

void TSFVector::axpy(const TSFReal& a, const TSFVector& x, const TSFVector& y)
{
	/* we may copy y into this, so if x==this we need to make a copy of x before
	 * modifying this */
	TSFVector w;
	if (isIdenticalTo(x))
		{
			w = x.copy();
		}
	else
		{
			w = x;
		}
	
	/* save a copy if this==y */
	if (!isIdenticalTo(y))
		{
			acceptCopyOf(y);
		}
	selfModifyingAxpy(a, w);
}

void TSFVector::dotStar(const TSFVector& x, const TSFVector& y)
{
	/* we may copy y into this, so if x==this we need to make a copy of x before
	 * modifying this */
	TSFVector w;
	if (isIdenticalTo(x))
		{
			w = x.copy();
		}
	else
		{
			w = x;
		}
	
	/* save a copy if this==y */
	if (!isIdenticalTo(y))
		{
			acceptCopyOf(y);
		}
	selfModifyingDotStar(w);
}

void TSFVector::dotSlash(const TSFVector& x, const TSFVector& y)
{
	/* we may copy y into this, so if x==this we need to make a copy of x before
	 * modifying this */
	TSFVector w;
	if (isIdenticalTo(x))
		{
			w = x.copy();
		}
	else
		{
			w = x;
		}
	
	/* save a copy if this==y */
	if (!isIdenticalTo(y))
		{
			acceptCopyOf(y);
		}
	selfModifyingDotSlash(w);
}

void TSFVector::scalarMult(const TSFReal& a, const TSFVector& x)
{
	/* save a copy if this==x */
	if (!isIdenticalTo(x))
		{
			acceptCopyOf(x);
		}
	selfModifyingScalarMult(a);
}


void TSFVector::acceptCopyOf(const TSFVector& x)
{
	if (ptr_.isNull())
		{
			TSFVector tmp = x.space().createMember();
			ptr_ = tmp.ptr_;
		}
	ptr_->acceptCopyOf(x);
}

void TSFVector::selfModifyingAxpy(const TSFReal& a, 
																	const TSFVector& x)
{
	TSFTimeMonitor t(opTimer());
	ptr_->axpy(a, x);
}

void TSFVector::selfModifyingScalarMult(const TSFReal& a)
{
	TSFTimeMonitor t(opTimer());
	ptr_->scalarMult(a);
}

void TSFVector::selfModifyingDotStar(const TSFVector& x)
{
	TSFTimeMonitor t(opTimer());
	ptr_->dotStar(x);
}

void TSFVector::selfModifyingDotSlash(const TSFVector& x)
{
	TSFTimeMonitor t(opTimer());
	ptr_->dotSlash(x);
}

TSFTimer& TSFVector::opTimer()
{
	static TSFSmartPtr<TSFTimer> timer= TSFTimer::getNewTimer("TSFVector operations");
	return *timer;
}

TSFTimer& TSFVector::copyTimer()
{
	static TSFSmartPtr<TSFTimer> timer= TSFTimer::getNewTimer("TSFVector copy");
	return *timer;
}

