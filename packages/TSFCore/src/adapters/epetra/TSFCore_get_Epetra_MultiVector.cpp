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

// ///////////////////////////////////////////////////////
// TSFCoreget_Epetra_MultiVector.cpp

#include "TSFCore_get_Epetra_MultiVector.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#ifdef TSFCORE_EPETRA_USE_TSFCORE_EPETRA_VECTOR
#  include "TSFCoreEpetraVector.hpp"
#endif
#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreVectorMultiVector.hpp"
#include "TSFCoreExplicitMultiVectorView.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

Teuchos::RefCountPtr<Epetra_MultiVector>
TSFCore::get_Epetra_MultiVector(
	const EpetraVectorSpace                              &vs
	,const Teuchos::RefCountPtr<MultiVector<double> >    &mv
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		!vs.isCompatible(*mv->range()), Exceptions::IncompatibleVectorSpaces
		,"TSFCore::get_Epetra_MultiVector(vs,mv): Error, must have compatible vector spaces!"
		);
#endif
	//
	// First, try to dynamic_cast to get the Epetra object directly
	// since this has proven to be the fastest way.  Note, it is not
	// trivial to create an Epetra view of an object and a lot of
	// dynamic memory allocations have to be performed to make this
	// work.  One must be very careful not to create
	// Epetra_MultiVector or Epetra_Vector objects on the fly!
	//
	EpetraMultiVector
		*epetra_mv_adapter = dynamic_cast<EpetraMultiVector*>(&*mv);
	if(epetra_mv_adapter)
		return epetra_mv_adapter->epetra_multi_vec(); // Note, this is dangerous!
#ifdef TSFCORE_EPETRA_USE_TSFCORE_EPETRA_VECTOR
  //
  // Second, try to dynamic_cast to EpetraVector (which derives from
  // MultiVector through Vector).
  //
	EpetraVector
		*epetra_v_adapter = dynamic_cast<EpetraVector*>(&*mv);
	if(epetra_v_adapter)
		return epetra_v_adapter->epetra_vec(); // Note, this is dangerous!
#else
  //
  // Second, try to dynamic_cast to VectorMultiVector which is
  // used to create a Vector implementation from EpetraMultiVector.
  // Note, other implementations may also use VectorMultiVector
  // but we don't want to involve these.
  //
  VectorMultiVector<double>
    *v_mv_adapter = dynamic_cast<VectorMultiVector<double>*>(&*mv);
  if(v_mv_adapter) {
    epetra_mv_adapter = dynamic_cast<EpetraMultiVector*>(&*mv);
    if(epetra_mv_adapter)
      return epetra_mv_adapter->epetra_multi_vec();
  }
#endif
	//
	// The assumption that we (rightly) make here is that if the
	// vector spaces are compatible, that either the multi-vectors are
	// both in-core or the vector spaces are both derived from
	// MPIVectorSpaceBase and have compatible maps.  Remarkably, this
	// code works in the cases serial and or parallel!
	// 
	const Index localOffset = vs.localOffset();
	const Index localSubDim = vs.localSubDim();
	//
	// Here we will extract a view of the local elements in the
	// underlying MultiVector.  Note, that we even do this if mv is a
	// direct subclass of EpetraMultiVector.  The reason that we don't
	// do a dynamic_cast is that the by explicitly calling the
	// SubMultiVector(...) access methods we insure that only the
	// elements of the underlying Epetra MultiVector are changed and
	// when the update is complete the commitSubMultiVector(...)
	// method will give the subclass a chance to update anything else
	// that it needs to update.  In most cases, no data will be
	// allocated or copy and only some small objects will be created
	// and a few virtual functions will be called so the overhead is
	// low and fixed.
	//
	// Create a *mutable* view of the local elements, this view will
	// be set on the RefCountPtr that is returned and then view will
	// be relased when the returned Epetra_MultiVector is returned.
	Teuchos::RefCountPtr<ExplicitMutableMultiVectorView<double> >
		emmvv = Teuchos::rcp(
			new ExplicitMutableMultiVectorView<double>(
				*mv
				,Range1D(localOffset+1,localOffset+localSubDim)
				)
			);
	// Create a temporary Epetra_MultiVector object and give it
	// the above local view
	Teuchos::RefCountPtr<Epetra_MultiVector>
		epetra_mv = Teuchos::rcp(
			new Epetra_MultiVector(
				::View                                  // CV
				,*vs.epetra_map()                       // Map
				,const_cast<double*>(emmvv->values())   // A
				,emmvv->leadingDim()                    // MyLDA
				,emmvv->numSubCols()                    // NumVectors
				)
			);
	// Give the explict view object to the above
	// Epetra_MultiVector smart pointer object.  In this way, When
	// the client is finished with the Epetra_MultiVector view The
	// destructor from the object in emmvv will automatically
	// commit the changes to the elements in the input mv
	// MultiVector object (reguardless of its implementation).
	// This is truly an elegant result!
	Teuchos::set_extra_data( emmvv, "emmvv", &epetra_mv );
	return epetra_mv;
}

Teuchos::RefCountPtr<const Epetra_MultiVector>
TSFCore::get_Epetra_MultiVector(
	const EpetraVectorSpace                                   &vs 
	,const Teuchos::RefCountPtr<const MultiVector<double> >   &mv
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		!vs.isCompatible(*mv->range()), Exceptions::IncompatibleVectorSpaces
		,"TSFCore::get_Epetra_MultiVector(vs,mv): Error, must have compatible vector spaces!"
		);
#endif
	//
	// First try to dynamic_cast to get the Epetra object directly
	// since this has proven to be the fastest way.  Note, it is not
	// trivial to create an Epetra view of an object and a lot of
	// dynamic memory allocations have to be performed to make this
	// work.  One must be very careful not to create
	// Epetra_MultiVector or Epetra_Vector objects on the fly!
	//
	const EpetraMultiVector
		*epetra_mv_adapter = dynamic_cast<const EpetraMultiVector*>(&*mv);
	if(epetra_mv_adapter)
		return epetra_mv_adapter->epetra_multi_vec();
#ifdef TSFCORE_EPETRA_USE_TSFCORE_EPETRA_VECTOR
  //
  // Second, try to dynamic_cast to EpetraVector (which derives from
  // MultiVector through Vector).
  //
	const EpetraVector
		*epetra_v_adapter = dynamic_cast<const EpetraVector*>(&*mv);
	if(epetra_v_adapter)
		return epetra_v_adapter->epetra_vec(); // Note, this is dangerous!
#else
  //
  // Second, try to dynamic_cast to VectorMultiVector which is
  // used to create a Vector implementation from EpetraMultiVector.
  // Note, other implementations may also use VectorMultiVector
  // but we don't want to involve these.
  //
  const VectorMultiVector<double>
    *v_mv_adapter = dynamic_cast<const VectorMultiVector<double>*>(&*mv);
  if(v_mv_adapter) {
    epetra_mv_adapter = dynamic_cast<const EpetraMultiVector*>(&*mv);
    if(epetra_mv_adapter)
      return epetra_mv_adapter->epetra_multi_vec();
  }
#endif
	//
	// Same as above function except as stated below
	//
	const Index localOffset = vs.localOffset();
	const Index localSubDim = vs.localSubDim();
	// Get an explicit *non-mutable* view of all of the elements in
	// the multi vector.
	Teuchos::RefCountPtr<ExplicitMultiVectorView<double> >
		emvv = Teuchos::rcp(
			new ExplicitMultiVectorView<double>(
				*mv
				,Range1D(localOffset+1,localOffset+localSubDim)
				)
			);
	// Create a temporary Epetra_MultiVector object and give it
	// the above view
	Teuchos::RefCountPtr<Epetra_MultiVector>
		epetra_mv = Teuchos::rcp(
			new Epetra_MultiVector(
				::View                                  // CV
				,*vs.epetra_map()                       // Map
				,const_cast<double*>(emvv->values())    // A
				,emvv->leadingDim()                     // MyLDA
				,emvv->numSubCols()                     // NumVectors
				)
			);
	// This will cause the destructor to free the view if needed (see
	// above function).  Since this view is non-mutable, only a
	// freeSubMultiVector(...)  and not a commit will be called.  This
	// is the whole reason there is a seperate implementation for the
	// const and non-const cases.
	Teuchos::set_extra_data( emvv, "emvv", &epetra_mv );
	return epetra_mv;
}
