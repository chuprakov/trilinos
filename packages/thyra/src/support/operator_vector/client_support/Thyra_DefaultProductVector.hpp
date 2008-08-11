// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_PRODUCT_VECTOR_HPP
#define THYRA_PRODUCT_VECTOR_HPP


#include "Thyra_DefaultProductVectorDecl.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Teuchos_Workspace.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template <class Scalar>
DefaultProductVector<Scalar>::DefaultProductVector()
{
  uninitialize();
}


template <class Scalar>
DefaultProductVector<Scalar>::DefaultProductVector(
  const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> >  &productSpace
  )
{
  initialize(productSpace);
}


template <class Scalar>
DefaultProductVector<Scalar>::DefaultProductVector(
  const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> >  &productSpace
  ,const Teuchos::RCP<VectorBase<Scalar> >                      vecs[]
  )
{
  initialize(productSpace,vecs);
}


template <class Scalar>
DefaultProductVector<Scalar>::DefaultProductVector(
  const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> >  &productSpace
  ,const Teuchos::RCP<const VectorBase<Scalar> >                vecs[]
  )
{
  initialize(productSpace,vecs);
}


template <class Scalar>
void DefaultProductVector<Scalar>::initialize(
  const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> >  &productSpace
  )
{
  // ToDo: Validate input!
  numBlocks_ = productSpace->numBlocks();
  productSpace_ = productSpace;
  vecs_.resize(numBlocks_);
  for( int k = 0; k < numBlocks_; ++k )
    vecs_[k].initialize(createMember(productSpace->getBlock(k)));
}


template <class Scalar>
void DefaultProductVector<Scalar>::initialize(
  const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> >  &productSpace
  ,const Teuchos::RCP<VectorBase<Scalar> >                      vecs[]
  )
{
  // ToDo: Validate input!
  numBlocks_ = productSpace->numBlocks();
  productSpace_ = productSpace;
  vecs_.resize(numBlocks_);
  for( int k = 0; k < numBlocks_; ++k )
    vecs_[k].initialize(vecs[k]);
}


template <class Scalar>
void DefaultProductVector<Scalar>::initialize(
  const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> >  &productSpace
  ,const Teuchos::RCP<const VectorBase<Scalar> >                vecs[]
  )
{
  // ToDo: Validate input!
  numBlocks_ = productSpace->numBlocks();
  productSpace_ = productSpace;
  vecs_.resize(numBlocks_);
  for( int k = 0; k < numBlocks_; ++k )
    vecs_[k].initialize(vecs[k]);
}


template <class Scalar>
void DefaultProductVector<Scalar>::uninitialize()
{
  productSpace_ = Teuchos::null;
  vecs_.resize(0);
  numBlocks_ = 0;
}


// Overridden from Teuchos::Describable

                                                
template<class Scalar>
std::string DefaultProductVector<Scalar>::description() const
{
  std::ostringstream oss;
  oss
    << Teuchos::Describable::description()
    << "{"
    << "dim="<<this->space()->dim()
    << ",numBlocks = "<<numBlocks_
    << "}";
  return oss.str();
}


template<class Scalar>
void DefaultProductVector<Scalar>::describe(
  Teuchos::FancyOStream &out_arg,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  using Teuchos::RCP;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::describe;
  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      *out
        << Teuchos::Describable::description() << "{"
        << "dim=" << this->space()->dim()
        << "}\n";
      OSTab tab(out);
      *out
        <<  "numBlocks="<< numBlocks_ << std::endl
        <<  "Constituent vector objects v[0], v[1], ... v[numBlocks-1]:\n";
      tab.incrTab();
      for( int k = 0; k < numBlocks_; ++k ) {
        *out << "v["<<k<<"] = " << describe(*vecs_[k].getConstObj(),verbLevel);
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


// Extensions to ProductVectorBase suitable for physically-blocked vectors


template <class Scalar>
void DefaultProductVector<Scalar>::setBlock(
  int i, const Teuchos::RCP<const VectorBase<Scalar> >& b
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(i < 0 || i >= numBlocks_);
  TEST_FOR_EXCEPT(!productSpace_->getBlock(i)->isCompatible(*(b->space())));
#endif
  vecs_[i] = b;
}


template <class Scalar>
void DefaultProductVector<Scalar>::setNonconstBlock(
  int i, const Teuchos::RCP<VectorBase<Scalar> >& b
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(i < 0 || i >= numBlocks_);
  TEST_FOR_EXCEPT(!productSpace_->getBlock(i)->isCompatible(*(b->space())));
#endif
  vecs_[i] = b;
}


// Overridden from ProductVectorBase


template <class Scalar>
Teuchos::RCP<VectorBase<Scalar> >
DefaultProductVector<Scalar>::getNonconstVectorBlock(const int k)
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( k < 0 || numBlocks_-1 < k);
#endif
  return vecs_[k].getNonconstObj();
}


template <class Scalar>
Teuchos::RCP<const VectorBase<Scalar> >
DefaultProductVector<Scalar>::getVectorBlock(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( k < 0 || numBlocks_-1 < k);
#endif
  return vecs_[k].getConstObj();
}


// Overridden from ProductMultiVectorBase


template <class Scalar>
Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >
DefaultProductVector<Scalar>::productSpace() const
{
  return productSpace_;
}


template <class Scalar>
bool DefaultProductVector<Scalar>::blockIsConst(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( k < 0 || numBlocks_-1 < k);
#endif
  return vecs_[k].isConst();
}


template <class Scalar>
Teuchos::RCP<MultiVectorBase<Scalar> >
DefaultProductVector<Scalar>::getNonconstMultiVectorBlock(const int k)
{
  return getNonconstVectorBlock(k);
}


template <class Scalar>
Teuchos::RCP<const MultiVectorBase<Scalar> >
DefaultProductVector<Scalar>::getMultiVectorBlock(const int k) const
{
  return getVectorBlock(k);
}


// Overridden from VectorBase


template <class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultProductVector<Scalar>::space() const
{
  return productSpace_;
}


template <class Scalar>
void DefaultProductVector<Scalar>::applyOpImpl(
  const RTOpPack::RTOpT<Scalar> &op,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
  const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Index first_ele_offset_in,
  const Index sub_dim_in,
  const Index global_offset_in
  ) const
{

  // 2008/02/20: rabartl: ToDo: Upgrade Teuchos::Workspace<T> to implicitly
  // convert to Teuchos::ArrayView<T>.  This will allow the calls to
  // applyOp(...) with sub_vecs and sub_targ_vecs to work without trouble!
  // For now, I just want to get this done.  It is likely that this function
  // is going to change in major ways soon anyway!

  //using Teuchos::Workspace;
  using Teuchos::ptr_dynamic_cast;
  using Teuchos::null;

  //Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  const Index	n = productSpace_->dim();
  const int num_vecs = vecs.size();
  const int num_targ_vecs = targ_vecs.size();

  // Validate the compatibility of the vectors!
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    !(0 <= first_ele_offset_in && first_ele_offset_in < n), std::out_of_range
    ,"DefaultProductVector::applyOp(...): Error, "
    "first_ele_offset_in = " << first_ele_offset_in << " is not in range [0,"<<(n-1)<<"]" );
  TEST_FOR_EXCEPTION(
    global_offset_in < 0, std::invalid_argument
    ,"DefaultProductVector::applyOp(...): Error, "
    "global_offset_in = " << global_offset_in << " is not acceptable" );
  TEST_FOR_EXCEPTION(
    sub_dim_in >= 0 && sub_dim_in - first_ele_offset_in > n, std::length_error
    ,"DefaultProductVector::applyOp(...): Error, "
    "global_offset_in = " << global_offset_in << ", sub_dim_in = " << sub_dim_in
    << "first_ele_offset_in = " << first_ele_offset_in << " and n = " << n
    << " are not compatible" );
  bool test_failed;
  for(int k = 0; k < num_vecs; ++k) {
    test_failed = !this->space()->isCompatible(*vecs[k]->space());
    TEST_FOR_EXCEPTION(
      test_failed, Exceptions::IncompatibleVectorSpaces
      ,"DefaultProductVector::applyOp(...): Error vecs["<<k<<"]->space() "
      <<"of type \'"<<typeName(*vecs[k]->space())<<"\' is not compatible with this "
      <<"\'VectorSpaceBlocked\' vector space!"
      );
  }
  for(int k = 0; k < num_targ_vecs; ++k) {
    test_failed = !this->space()->isCompatible(*targ_vecs[k]->space());
    TEST_FOR_EXCEPTION(
      test_failed, Exceptions::IncompatibleVectorSpaces
      ,"DefaultProductVector::applyOp(...): Error targ_vecs["<<k<<"]->space() "
      <<"of type \'"<<typeName(*vecs[k]->space())<<"\' is not compatible with this "
      <<"\'VectorSpaceBlocked\' vector space!"
      );
  }
#endif

  //
  // The first thing that we do is to see if all of the vectors involved are
  // incore, serial vectors. In this case, the vectors should be compatible.
  // To accomplish this we will pick the continguous vector to implement
  // the applyOp(...) function.
  //
  // Get the index of an incore-only input vector and input/output vector?
  const bool this_isInCore = productSpace_->hasInCoreView(
    Range1D(),VIEW_TYPE_DETACHED,STRIDE_TYPE_NONUNIT
    );
  int incore_vec_k = -1, incore_targ_vec_k = -1;

  // Dynamic cast the pointers for the vector arguments
  Array<Ptr<const ProductVectorBase<Scalar> > >
    vecs_args(num_vecs);
  for(int k = 0; k < num_vecs; ++k) {
    vecs_args[k] = ptr_dynamic_cast<const ProductVectorBase<Scalar> >(vecs[k]);
    if( vecs_args[k] == null ) {
      const bool isInCore_k = vecs[k]->space()->hasInCoreView();
      if( this_isInCore && isInCore_k ) {
        incore_vec_k = k;
        break;
      }
      TEST_FOR_EXCEPTION(
        !this_isInCore || (this_isInCore && !isInCore_k),
        Exceptions::IncompatibleVectorSpaces
        ,"DefaultProductVector::applyOp(...): Error vecs["<<k<<"] "
        <<"of type \'"<<typeName(*vecs[k])<<"\' does not support the "
        <<"\'DefaultProductVector<Scalar>\' interface and is not an incore"
        "vector or this is not an incore vector!"
        );
    }
  }
  Array<Ptr<ProductVectorBase<Scalar> > >
    targ_vecs_args(num_targ_vecs);
  for(int k = 0; k < num_targ_vecs; ++k) {
    targ_vecs_args[k] = ptr_dynamic_cast<ProductVectorBase<Scalar> >(targ_vecs[k]);
    if( targ_vecs_args[k] == null ) {
      const bool isInCore_k = targ_vecs[k]->space()->hasInCoreView();
      if( this_isInCore && isInCore_k ) {
        incore_targ_vec_k = k;
        break;
      }
      TEST_FOR_EXCEPTION(
        !this_isInCore || (this_isInCore && !isInCore_k),
        Exceptions::IncompatibleVectorSpaces
        ,"DefaultProductVector::applyOp(...): Error targ_vecs["<<k<<"] "
        <<"of type \'"<<typeName(*targ_vecs[k])<<"\' does not support the "
        <<"\'DefaultProductVector<Scalar>\' interface and is not an incore"
        " vector or this is not an incore vector!"
        );
    }
  }

  // Let a incore-only vector with a contiguous view handle this through
  // explicit vector access?
  if( incore_vec_k >= 0 ) {
    vecs[incore_vec_k]->applyOp(
      op, vecs, targ_vecs, reduct_obj,
      first_ele_offset_in, sub_dim_in, global_offset_in );
    return;
  }
  else if ( incore_targ_vec_k >= 0 ) {
    targ_vecs[incore_targ_vec_k]->applyOp(
      op, vecs, targ_vecs, reduct_obj,
      first_ele_offset_in, sub_dim_in, global_offset_in
      );
    return;
  }

  //
  // If we get here, then we will implement the applyOpImpl(...) one vector
  // block at a time.
  //
  const Index this_dim = n;
  const Index sub_dim
    = (
      sub_dim_in < 0
      ? this_dim - first_ele_offset_in
      : sub_dim_in
      );
  Index num_elements_remaining = sub_dim;
  const int numBlocks = productSpace_->numBlocks();
  Array<RCP<const VectorBase<Scalar> > >
    sub_vecs_rcps(num_vecs);
  Array<Ptr<const VectorBase<Scalar> > >
    sub_vecs(num_vecs);
  Array<RCP<VectorBase<Scalar> > >
    sub_targ_vecs_rcps(num_targ_vecs);
  Array<Ptr<VectorBase<Scalar> > >
    sub_targ_vecs(num_targ_vecs);
  Index g_off = -first_ele_offset_in;
  for(int k = 0; k < numBlocks; ++k) {
    const Index local_dim = productSpace_->getBlock(k)->dim();
    if( g_off < 0 && -g_off+1 > local_dim ) {
      g_off += local_dim;
      continue;
    }
    const Index
      local_sub_dim
      = ( g_off >= 0
        ? std::min( local_dim, num_elements_remaining )
        : std::min( local_dim + g_off, num_elements_remaining ) );
    if( local_sub_dim <= 0 )
      break;
    // Fill constituent vectors for block k
    for( int i = 0; i < num_vecs; ++i ) {
      sub_vecs_rcps[i] = vecs_args[i]->getVectorBlock(k);
      sub_vecs[i] = sub_vecs_rcps[i].ptr();
    }
    // Fill constituent target vectors for block k
    for( int j = 0; j < num_targ_vecs; ++j ) {
      sub_targ_vecs_rcps[j] = targ_vecs_args[j]->getNonconstVectorBlock(k);
      sub_targ_vecs[j] = sub_targ_vecs_rcps[j].ptr();
    }
    Thyra::applyOp<Scalar>(
      op, sub_vecs(), sub_targ_vecs(),
      reduct_obj,
      g_off < 0 ? -g_off : 0, // first_ele_offset
      local_sub_dim, // sub_dim
      g_off < 0 ? global_offset_in : global_offset_in + g_off // global_offset
      );
    g_off += local_dim;
    num_elements_remaining -= local_sub_dim;
  }
  TEST_FOR_EXCEPT(!(num_elements_remaining==0));

}


// protected


// Overridden protected functions from VectorBase


template <class Scalar>
void DefaultProductVector<Scalar>::acquireDetachedVectorViewImpl(
  const Range1D& rng_in, RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
  const Range1D
    rng = rng_in.full_range() ? Range1D(0,productSpace_->dim()-1) : rng_in;
  int    kth_vector_space  = -1;
  Index  kth_global_offset = 0;
  productSpace_->getVecSpcPoss(rng.lbound(),&kth_vector_space,&kth_global_offset);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
  if(
    rng.lbound() + rng.size()
    <= kth_global_offset + vecs_[kth_vector_space].getConstObj()->space()->dim()
    )
  {
    // This involves only one sub-vector so just return it.
    const_cast<const VectorBase<Scalar>*>(
      &*vecs_[kth_vector_space].getConstObj()
      )->acquireDetachedView( rng - kth_global_offset, sub_vec );
    sub_vec->setGlobalOffset( sub_vec->globalOffset() + kth_global_offset );
  }
  else {
    // Just let the default implementation handle this.  ToDo: In the future
    // we could manually construct an explicit sub-vector that spanned
    // two or more constituent vectors but this would be a lot of work.
    // However, this would require the use of temporary memory but
    // so what.
    VectorDefaultBase<Scalar>::acquireDetachedVectorViewImpl(rng_in,sub_vec);
  }
}


template <class Scalar>
void DefaultProductVector<Scalar>::releaseDetachedVectorViewImpl(
  RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
  if( sub_vec->values().get() == NULL ) return;
  int    kth_vector_space  = -1;
  Index  kth_global_offset = 0;
  productSpace_->getVecSpcPoss(sub_vec->globalOffset(),&kth_vector_space,&kth_global_offset);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
  if(
    sub_vec->globalOffset() + sub_vec->subDim()
    <= kth_global_offset +  vecs_[kth_vector_space].getConstObj()->space()->dim()
    )
  {
    // This sub_vec was extracted from a single constituent vector
    sub_vec->setGlobalOffset( sub_vec->globalOffset() - kth_global_offset );
    vecs_[kth_vector_space].getConstObj()->releaseDetachedView(sub_vec);
  }
  else {
    // This sub_vec was created by the default implementation!
    VectorDefaultBase<Scalar>::releaseDetachedVectorViewImpl(sub_vec);
  }
}


template <class Scalar>
void DefaultProductVector<Scalar>::acquireNonconstDetachedVectorViewImpl(
  const Range1D& rng_in, RTOpPack::SubVectorView<Scalar>* sub_vec
  )
{
  const Range1D
    rng = rng_in.full_range() ? Range1D(0,productSpace_->dim()-1) : rng_in;
  int    kth_vector_space  = -1;
  Index  kth_global_offset = 0;
  productSpace_->getVecSpcPoss(rng.lbound(),&kth_vector_space,&kth_global_offset);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
  if(
    rng.lbound() + rng.size()
    <= kth_global_offset + vecs_[kth_vector_space].getConstObj()->space()->dim()
    )
  {
    // This involves only one sub-vector so just return it.
    vecs_[kth_vector_space].getConstObj()->acquireDetachedView(
      rng - kth_global_offset, sub_vec
      );
    sub_vec->setGlobalOffset( sub_vec->globalOffset() + kth_global_offset );
  }
  else {
    // Just let the default implementation handle this.  ToDo: In the future
    // we could manually construct an explicit sub-vector that spanned
    // two or more constituent vectors but this would be a lot of work.
    // However, this would require the use of temporary memory but
    // so what.
    VectorDefaultBase<Scalar>::acquireNonconstDetachedVectorViewImpl(rng_in,sub_vec);
  }
}


template <class Scalar>
void DefaultProductVector<Scalar>::commitNonconstDetachedVectorViewImpl(
  RTOpPack::SubVectorView<Scalar>* sub_vec
  )
{
  if( sub_vec->values().get() == NULL ) return;
  int    kth_vector_space  = -1;
  Index  kth_global_offset = 0;
  productSpace_->getVecSpcPoss(sub_vec->globalOffset(),&kth_vector_space,&kth_global_offset);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
  if(
    sub_vec->globalOffset() + sub_vec->subDim()
    <= kth_global_offset +  vecs_[kth_vector_space].getConstObj()->space()->dim()
    )
  {
    // This sub_vec was extracted from a single constituent vector
    sub_vec->setGlobalOffset( sub_vec->globalOffset() - kth_global_offset );
    vecs_[kth_vector_space].getNonconstObj()->commitDetachedView(sub_vec);
  }
  else {
    // This sub_vec was created by the default implementation!
    VectorDefaultBase<Scalar>::commitNonconstDetachedVectorViewImpl(sub_vec);
  }
}


template <class Scalar>
void DefaultProductVector<Scalar>::setSubVectorImpl(
  const RTOpPack::SparseSubVectorT<Scalar>& sub_vec
  )
{
  int    kth_vector_space  = -1;
  Index  kth_global_offset = 0;
  productSpace_->getVecSpcPoss(sub_vec.globalOffset(),&kth_vector_space,&kth_global_offset);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
  if(
    sub_vec.globalOffset() + sub_vec.subDim()
    <= kth_global_offset + vecs_[kth_vector_space].getConstObj()->space()->dim()
    )
  {
    // This sub-vector fits into a single constituent vector
    RTOpPack::SparseSubVectorT<Scalar> sub_vec_g = sub_vec;
    sub_vec_g.setGlobalOffset( sub_vec_g.globalOffset() - kth_global_offset );
    vecs_[kth_vector_space].getNonconstObj()->setSubVector(sub_vec_g);
  }
  else {
    // Let the default implementation take care of this.  ToDo: In the future
    // it would be possible to manually set the relevant constituent
    // vectors with no temp memory allocations.
    VectorDefaultBase<Scalar>::setSubVector(sub_vec);
  }
}


} // namespace Thyra


#endif // THYRA_PRODUCT_VECTOR_HPP
