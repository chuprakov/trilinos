// ///////////////////////////////////////////////////////
// TSFCoreget_Epetra_MultiVector.cpp

#include "TSFCore_get_Epetra_MultiVector.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreExplicitMultiVectorView.hpp"
#include "Epetra_MultiVector.h"
#include "dynamic_cast_verbose.hpp"

Teuchos::RefCountPtr<Epetra_MultiVector>
TSFCore::get_Epetra_MultiVector(
	const EpetraVectorSpace                              &vs
	,const Teuchos::RefCountPtr<MultiVector<double> >    &mv
	)
{
	using DynamicCastHelperPack::dyn_cast;
	const bool isInCore = vs.isInCore();
	if(isInCore) {
		//
		// Here we will extract a view of all of the elements in the
		// underlying MultiVector since this is easy to do when all of
		// the elements are in-core on this processor.  Note, that we
		// even do this if mv is a direct subclass of
		// EpetraMultiVector.  The reason that we don't do a
		// dynamic_cast is that the by explicitly calling the
		// SubMultiVector(...) access methods we insure that only the
		// elements of the underlying Epetra MultiVector are changed
		// and when the update is complete the
		// commitSubMultiVector(...)  method will give the subclass a
		// chance to update anything else that it needs to update.  In
		// most cases, no data will be allocated or copy and only some
		// small objects will be created and a few virtual functions
		// will be called so the overhead is low and fixed.
		//
		// Get an explicit *mutable* view of all of the elements in
		// the multi vector
		Teuchos::RefCountPtr<ExplicitMutableMultiVectorView<double> >
			emmvv = Teuchos::rcp(new ExplicitMutableMultiVectorView<double>(*mv));
		// Create a temporary Epetra_MultiVector object and give it
		// the above view
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
		Teuchos::set_extra_data( emmvv, &epetra_mv );
		return epetra_mv;
	}
	else {
		const MPIVectorSpaceBase<double>
			&mpi_mv_spc = dyn_cast<const MPIVectorSpaceBase<double> >(*mv->range());
		// Create a temporary Epetra_MultiVector object, give it a view
		// of the local elements on each processor, and then wrap in a 
		// EpetraMultiVector object and return!
		assert(0); // ToDo: Implement!
	}
	assert(0); // Should never get here (but some compilers complain if this return statement is not here!)
	return Teuchos::null;
}

Teuchos::RefCountPtr<const Epetra_MultiVector>
TSFCore::get_Epetra_MultiVector(
	const EpetraVectorSpace                                   &vs 
	,const Teuchos::RefCountPtr<const MultiVector<double> >   &mv
	)
{
	using DynamicCastHelperPack::dyn_cast;
	const bool isInCore = vs.isInCore();
	if(isInCore) {
		// Same as above function except ...
		//
		// Get an explicit *non-mutable* view of all of the elements
		// in the multi vector
		Teuchos::RefCountPtr<ExplicitMultiVectorView<double> >
			emvv = Teuchos::rcp(new ExplicitMultiVectorView<double>(*mv));
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
		// This will cause the destructor to free the view if needed (see above function)
		Teuchos::set_extra_data( emvv, &epetra_mv );
		return epetra_mv;
	}
	else {
		const MPIVectorSpaceBase<double>
			&mpi_mv_spc = dyn_cast<const MPIVectorSpaceBase<double> >(*mv->range());
		// Create a temporary Epetra_MultiVector object, give it a view
		// of the local elements on each processor, and then wrap in a 
		// EpetraMultiVector object and return!
		assert(0); // ToDo: Implement!
	}
	assert(0); // Should never get here (but some compilers complain if this return statement is not here!)
	return Teuchos::null;
}
