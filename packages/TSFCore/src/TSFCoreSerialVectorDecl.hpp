// /////////////////////////////////////////////////////////////////
// SerialVectorDecl.hpp

#ifndef TSFCORE_VECTOR_SERIAL_DECL_HPP
#define TSFCORE_VECTOR_SERIAL_DECL_HPP

#include "TSFCoreSerialVectorBaseDecl.hpp"
#include "TSFCoreSerialVectorSpaceDecl.hpp"

namespace TSFCore {

///
/** Implementation of serial vectors.
 *
 * This class can be used either as a view of a vector data or as a
 * storage for vector data.
 *
 * To create with storage with the dimension of <tt>dim</tt> just call
 * the constructor <tt>SerialVector(dim)</tt> or after construction
 * you can call <tt>this->initialize(dim)</tt>.
 *
 * To simply create a view of a vector <tt>v</tt> with stride
 * <tt>vs</tt>, without ownership just call
 * <tt>SerialVector(v,vs)</tt> or after construction call
 * <tt>this->initialize(v,vs)</tt>.
 */
template<class Scalar>
class SerialVector : public SerialVectorBase<Scalar> {
public:

	/** @name Constructors/initializers */
	//@{

	///
	/** Frees memory if <tt>this</tt> owns it.
	 */
	~SerialVector();

	///
	/** Calls <tt>this->initialize(dim)</tt>.
	 */
	SerialVector(
		int dim = 0
		);
	///
	/** Calls <tt>this->initialize(v,vs,dim,ownsMem)</tt>.
	 */
	SerialVector(
		Scalar  v[]
		,int    vs
		,int    dim
		,bool   ownsMem = false
		);
	///
	/** Call <tt>this->initialize(v,vs,true)</tt> with allocated data <tt>v</tt>.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->dim() == dim</tt>
	 * </ul>
	 */
	void initialize(
		int dim
		);
	///
	/** Initialize with a dense vector slice.
	 *
	 * 
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->dim() == dim</tt>
	 * <tt>sub_vec.values == v</tt>
	 * <tt>sub_vec.values_stride == vs</tt>
	 * </ul>
	 *  where <tt>sub_vec</tt> above is returned from <tt>this->getSubVector(Range1D(),&sub_vec)</tt>
	 */
	void initialize(
		Scalar  v[]
		,int    vs
		,int    dim
		,bool   ownsMem = false
		);

	//@}

	/** @name Overridden from Vector */
	//@{

	/// Returns 0 if <tt>this</tt> is not initialized yet.
	Index dim() const;
	///
	MemMngPack::ref_count_ptr< const VectorSpace<Scalar> > space() const;
	///
	void getSubVector( const Range1D& rng, RTOpPack::SubVectorT<Scalar>* sub_vec ) const;
	///
	void freeSubVector( RTOpPack::SubVectorT<Scalar>* sub_vec ) const;
	///
	void getSubVector( const Range1D& rng, RTOpPack::MutableSubVectorT<Scalar>* sub_vec );
	///
	void commitSubVector( RTOpPack::MutableSubVectorT<Scalar>* sub_vec );
	///
	void setSubVector( const RTOpPack::SparseSubVectorT<Scalar>& sub_vec );

	//@}

private:

	// ///////////////////////////////////////
	// Private data members
	
	Scalar                                                  *v_;
	int                                                     vs_;
	int                                                     dim_;
	bool                                                    ownsMem_;
	SerialVectorSpace<Scalar>                               space_serial_;

	// ////////////////////////////////
	// Private member functions

	void free_mem();

	// Not defined and not to be called
	SerialVector(const SerialVector&);
	SerialVector& operator=(const SerialVector&);

}; // end class SerialVector

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SERIAL_DECL_HPP
