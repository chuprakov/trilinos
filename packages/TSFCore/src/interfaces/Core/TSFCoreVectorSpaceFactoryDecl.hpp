// ///////////////////////////////////////////////////////////////
// TSFCoreVectorSpaceFactoryDecl.hpp

#ifndef TSFCORE_VECTOR_SPACE_FACTORY_DECL_HPP
#define TSFCORE_VECTOR_SPACE_FACTORY_DECL_HPP

#include "TSFCoreTypes.hpp"

namespace TSFCore {

///
/** Abstract interface for objects that can create vector spaces of a specified dimension.
 *
 * The primary role that a <tt>%VectorSpaceFactory</tt> object takes
 * is defined in the documentation for the class <tt>VectorSpace</tt>
 * and is related to the domain space of <tt>MultiVector</tt> objects.
 *
 * <b>Notes to subclass developers</b>
 * 
 * 
 */
template<class Scalar>
class VectorSpaceFactory {
public:

	/** @name Pure virtual functions that must be overridden */
	//@{

	///
	/** Create a vector space of the given dimension.
	 *
	 * @param  dim  [in] The dimension of the vector space to create.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>dim > 0</tt> (throw <tt>std::invalid_argument</tt>).
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get() != NULL</tt>
	 * <li> <tt>return->dim() == dim</tt>
	 * </ul>
	 *
	 * @return  Returns a smart reference counted pointer to a dynamically
	 * allocated vector space object that can be used to create vector.
	 */
	virtual Teuchos::RefCountPtr< const VectorSpace<Scalar> > createVecSpc(int dim) const = 0;

	//@}

}; // end class VectorSpaceFactory

} // end namespace TSFCore

#endif  // TSFCORE_VECTOR_SPACE_FACTORY_DECL_HPP
