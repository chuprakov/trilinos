// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
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

// /////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraMultiVector.cpp

#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreMPIMultiVectorBase.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "TSFCore_get_Epetra_MultiVector.hpp"
#include "Teuchos_TestForException.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialComm.h"
#include "Epetra_LocalMap.h"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace TSFCore {

// Constructors/Initializers

EpetraMultiVector::EpetraMultiVector()
{}

EpetraMultiVector::EpetraMultiVector(
  const Teuchos::RefCountPtr<Epetra_MultiVector>                        &epetra_multi_vec
  ,const Teuchos::RefCountPtr<const EpetraVectorSpace>                  &epetra_range
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
  ,const Teuchos::RefCountPtr<const EpetraVectorSpace>                  &epetra_domain
#else
  ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> > &domain
#endif
	)
{
	this->initialize(
    epetra_multi_vec,epetra_range
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
    ,epetra_domain
#else
    ,domain
#endif
    );
}

void EpetraMultiVector::initialize(
	const Teuchos::RefCountPtr<Epetra_MultiVector>                        &epetra_multi_vec
	,const Teuchos::RefCountPtr<const EpetraVectorSpace>                  &epetra_range
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
	,const Teuchos::RefCountPtr<const EpetraVectorSpace>                  &epetra_domain
#else
	,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> > &domain
#endif
	)
{
  using Teuchos::dyn_cast;
#ifdef _DEBUG
	const char err_msg[] = "EpetraMultiVector::initialize(...): Error!";
	TEST_FOR_EXCEPTION( epetra_multi_vec.get() == NULL, std::invalid_argument, err_msg ); 
	TEST_FOR_EXCEPTION( !epetra_range.get() || epetra_range->dim() == 0, std::invalid_argument, err_msg ); 
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
	TEST_FOR_EXCEPTION( epetra_domain.get() && epetra_domain->dim() == 0,std::invalid_argument, err_msg );
#else
	TEST_FOR_EXCEPTION( domain.get() && domain->dim() == 0,std::invalid_argument, err_msg );
#endif
	// ToDo: Check the compatibility of the vectors in col_vecs!
#endif // _DEBUG
	// Set multi-vector
	epetra_multi_vec_ = epetra_multi_vec;
	// set range
  epetra_range_  = epetra_range;
	// set domain
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
	if(epetra_domain.get()) {
		epetra_domain_  = epetra_domain;
	}
	else {
		epetra_domain_ = Teuchos::rcp_dynamic_cast<const EpetraVectorSpace>(
			epetra_range_->smallVecSpcFcty()->createVecSpc(epetra_multi_vec->NumVectors())
			,true);
	}
#else
	if(domain.get()) {
		domain_  = domain;
	}
	else {
		domain_ = Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(
			epetra_range_->smallVecSpcFcty()->createVecSpc(epetra_multi_vec->NumVectors())
			,true);
	}
#endif
	// Tell the base class object that the vector space is updated
	updateMpiSpace();
}

void EpetraMultiVector::setUninitialized(
	Teuchos::RefCountPtr<Epetra_MultiVector>                         *epetra_multi_vec
	,Teuchos::RefCountPtr<const EpetraVectorSpace>                   *epetra_range
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
	,Teuchos::RefCountPtr<const EpetraVectorSpace>                   *epetra_domain
#else
	,Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  *domain
#endif
	)
{
	if(epetra_multi_vec) *epetra_multi_vec = epetra_multi_vec_;;
	if(epetra_range) *epetra_range = epetra_range_;
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
	if(epetra_domain) *epetra_domain = epetra_domain_;
#else
	if(domain) *domain = domain_;
#endif
	epetra_multi_vec_ = Teuchos::null;
	epetra_range_ = Teuchos::null;
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
	epetra_domain_ = Teuchos::null;
#else
	domain_ = Teuchos::null;
#endif
	updateMpiSpace();
}

// Overridden from OpBase

Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<EpetraMultiVector::Scalar> >
EpetraMultiVector::domainScalarProdVecSpc() const
{
#ifdef TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE
	return epetra_domain_;
#else
  return domain_;
#endif
}

// Overridden from LinearOp

#if defined(TSFCORE_EPETRA_USE_EPETRA_MULTI_VECTOR_MULTIPLY) && defined(TSFCORE_EPETRA_USE_EPETRA_DOMAIN_VECTOR_SPACE)

void EpetraMultiVector::euclideanApply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
#ifdef _DEBUG
	// ToDo: Check for compatibility of vector spaces
#endif
	//
	// Get Epetra_MultiVector objects for the arguments
	//
	Teuchos::RefCountPtr<const Epetra_MultiVector>
		epetra_X = get_Epetra_MultiVector(
			M_trans==NOTRANS ? *epetra_domain_ : *epetra_range_
			,Teuchos::rcp(&X,false)
			);
	Teuchos::RefCountPtr<Epetra_MultiVector>
		epetra_Y = get_Epetra_MultiVector(
			M_trans==NOTRANS ? *epetra_range_ : *epetra_domain_
			,Teuchos::rcp(Y,false)
			);
	//
	// Do the multiplication
	//
	//   Y = scalarAB * transA(A) * transB(B) + scalarThis*Y
	//
	epetra_Y->Multiply(
		M_trans==NOTRANS ? 'N' : 'T'    // TransA
		,'N'                            // TransB
		,alpha                          // ScalarAB
		,*epetra_multi_vec_             // A
		,*epetra_X                      // B
		,beta                           // ScalarThis
		);
}

#endif

// Overridden from MultiVector

Teuchos::RefCountPtr<Vector<EpetraMultiVector::Scalar> >
EpetraMultiVector::col(Index j)
{
  using Teuchos::dyn_cast;
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( !(  1 <= j  && j <= domain()->dim() ), std::logic_error, "EpetraMultiVector::col(j): Error!" );
#endif
  const Teuchos::RefCountPtr<Epetra_Vector>
    epetra_col_vec = Teuchos::rcp( new Epetra_Vector( ::View, *epetra_range_->epetra_map(), (*epetra_multi_vec_)[j-1] ) );
	return Teuchos::rcp(
#ifdef TSFCORE_EPETRA_USE_TSFCORE_EPETRA_VECTOR
    new EpetraVector(epetra_col_vec,epetra_range_)
#else
    new VectorMultiVector<Scalar>(Teuchos::rcp(new EpetraMultiVector(epetra_col_vec,epetra_range_)))
#endif
		);
}

Teuchos::RefCountPtr<MultiVector<EpetraMultiVector::Scalar> >
EpetraMultiVector::subView( const Range1D& col_rng_in )
{
//	return MultiVector<Scalar>::subView(col_rng_in); // Uncomment to use the default implementation!
	const Range1D colRng = validateColRange(col_rng_in);
	return Teuchos::rcp(
		new EpetraMultiVector(
			Teuchos::rcp(
				new Epetra_MultiVector(
					::View                         // CV
					,*epetra_multi_vec_            // Source
					,colRng.lbound()-1             // StartIndex (zero-based)
					,colRng.size()                 // NumVectors
					)
				)
      ,epetra_range_
			)
		);
}

Teuchos::RefCountPtr<MultiVector<EpetraMultiVector::Scalar> >
EpetraMultiVector::subView( const int numCols, const int cols[] )
{
//	return MultiVector<Scalar>::subView(numCols,cols); // Uncomment to use the default implementation!
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
	// Translate from 1-based indexes to zero-based column indexes
#ifdef _DEBUG
	const int numTotalCols = this->domain()->dim();
#endif
	Workspace<int> zb_cols(wss,numCols,false);
	for( int k = 0; k < numCols; ++k ) {
#ifdef _DEBUG
		TEST_FOR_EXCEPTION(
			cols[k] < 1 || numTotalCols < cols[k], std::invalid_argument
			,"Error, cols["<<k<<"] = " << cols[k] << " does not fall in the range "
			"[1,"<<numTotalCols<<"]!"
			);
#endif
		zb_cols[k] = cols[k] - 1;
	}
	// Create the view
	return Teuchos::rcp(
		new EpetraMultiVector(
			Teuchos::rcp(
				new Epetra_MultiVector(
					::View                         // CV
					,*epetra_multi_vec_            // Source
					,&zb_cols[0]                   // Indices
					,numCols                       // NumVectors
					)
				)
      ,epetra_range_
			)
		);
}

// Overridden from MPIMultiVectorBase

Teuchos::RefCountPtr<const MPIVectorSpaceBase<EpetraMultiVector::Scalar> >
EpetraMultiVector::mpiSpace() const
{
	return epetra_range_;
}

void EpetraMultiVector::getLocalData( const Scalar **values, Index *leadingDim ) const
{
	if( epetra_multi_vec_->ConstantStride() ) {
		Scalar *non_const_values;
		epetra_multi_vec_->ExtractView( &non_const_values, leadingDim );
		*values = non_const_values; 
	}
	else {
		assert(0); // ToDo: Implement!
	}
}

void EpetraMultiVector::freeLocalData( const Scalar *values ) const
{
	if( epetra_multi_vec_->ConstantStride() ) {
		// There is no memory to free!
	}
	else {
		assert(0); // ToDo: Implement!
	}
}

void EpetraMultiVector::getLocalData( Scalar **values, Index *leadingDim )
{
	if( epetra_multi_vec_->ConstantStride() ) {
		epetra_multi_vec_->ExtractView( values, leadingDim );
	}
	else {
		assert(0); // ToDo: Implement!
	}
}

void EpetraMultiVector::commitLocalData( Scalar *values )
{
	if( epetra_multi_vec_->ConstantStride() ) {
		// The data was a direct view so there is
		// no need to commit anything or free anything.
	}
	else {
		assert(0); // ToDo: Implement!
	}
}
	
} // end namespace TSFCore
