// ///////////////////////////////////////////////////////
// TSFCoreget_Epetra_MultiVector.hpp

#ifndef TSFCORE_GET_EPETRA_MULTI_VECTOR_HPP
#define TSFCORE_GET_EPETRA_MULTI_VECTOR_HPP

#include "TSFCoreEpetraTypes.hpp"

namespace TSFCore {

///
/** Get a non-<tt>const</tt> <tt>Epetra_MultiVector</tt> object from a
 * non-<tt>const</tt> <tt>MultiVector</tt> object if possible.
 *
 * ToDo: Finish Documentation!
 */
Teuchos::RefCountPtr<Epetra_MultiVector>
get_Epetra_MultiVector(
	const EpetraVectorSpace                              &vs
	,const Teuchos::RefCountPtr<MultiVector<double> >    &mv
	);

///
/** Get a <tt>const</tt> <tt>Epetra_MultiVector</tt> object from a
 * <tt>const</tt> <tt>MultiVector</tt> object if possible.
 *
 * ToDo: Finish Documentation!
 */
Teuchos::RefCountPtr<const Epetra_MultiVector>
get_Epetra_MultiVector(
	const EpetraVectorSpace                                   &vs 
	,const Teuchos::RefCountPtr<const MultiVector<double> >   &mv
	);

} // namespace TSFCore

#endif // TSFCORE_GET_EPETRA_MULTI_VECTOR_HPP
