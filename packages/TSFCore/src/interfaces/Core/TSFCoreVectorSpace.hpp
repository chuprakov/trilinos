// //////////////////////////////////////////////////////////////////////
// TSFCoreVectorSpace.hpp

#ifndef TSFCORE_VECTOR_SPACE_HPP
#define TSFCORE_VECTOR_SPACE_HPP

#include "TSFCoreVectorSpaceDecl.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreSerialVectorSpaceFactory.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "TSFCoreMultiVectorCols.hpp"
#include "RTOp_ROp_dot_prod.h"
#include "RTOpCppC.hpp"

namespace TSFCore {

// Virtual functions with default implementations

template<class Scalar>
bool VectorSpace<Scalar>::isInCore() const
{
	return false;
}

template<class Scalar>
MemMngPack::ref_count_ptr<const VectorSpaceFactory<Scalar> >
VectorSpace<Scalar>::smallVecSpcFcty() const
{
	return MemMngPack::rcp(new SerialVectorSpaceFactory<Scalar>());
}

template<class Scalar>
MemMngPack::ref_count_ptr<MultiVector<Scalar> > 
VectorSpace<Scalar>::createMembers(int numMembers) const
{
	namespace mmp = MemMngPack;
	return mmp::rcp(new MultiVectorCols<Scalar> (mmp::rcp(this,false),this->smallVecSpcFcty()->createVecSpc(numMembers)));
}

template<class Scalar>
Scalar VectorSpace<Scalar>::scalarProd( const Vector<Scalar>& x, const Vector<Scalar>& y ) const
{
	const MultiVectorCols<Scalar>
		X( MemMngPack::rcp( const_cast<Vector<Scalar>*>(&x), false ) ),
		Y( MemMngPack::rcp( const_cast<Vector<Scalar>*>(&y), false ) );
	Scalar scalar_prods[1];
	this->scalarProds(X,Y,scalar_prods);
	return scalar_prods[0];
}

template<class Scalar>
void VectorSpace<Scalar>::scalarProds( const MultiVector<Scalar>& X, const MultiVector<Scalar>& Y, Scalar scalar_prods[] ) const
{
	dot(X,Y,scalar_prods);
}

template<class Scalar>
MemMngPack::ref_count_ptr< const VectorSpace<Scalar> >
VectorSpace<Scalar>::clone() const
{
	return MemMngPack::null;
}


} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SPACE_HPP
