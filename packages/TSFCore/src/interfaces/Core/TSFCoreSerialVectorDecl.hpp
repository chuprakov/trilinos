// /////////////////////////////////////////////////////////////////
// TSFCoreSerialVectorDecl.hpp

#ifndef TSFCORE_VECTOR_SERIAL_DECL_HPP
#define TSFCORE_VECTOR_SERIAL_DECL_HPP

#include "TSFCoreSerialVectorBase.hpp"
#include "TSFCoreSerialVectorSpace.hpp"

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
	/** Calls <tt>this->initialize(vecSpc)</tt>.
	 */
	SerialVector(
		const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vecSpc
		);
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
		,const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vecSpc = Teuchos::null
		);
	///
	/** Call <tt>this->initialize(v,vs,true,vecSpc)</tt> with internally dynamically allocated data <tt>v</tt>.
	 */
	void initialize(
		const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vecSpc
		);
	///
	/** Call <tt>this->initialize(v,vs,true)</tt> with internally dynamically allocated data <tt>v</tt>.
	 */
	void initialize(
		int dim
		);
	///
	/** Initialize with storage.
	 *
	 * @param  v      [in] Pointer to array of storage that <tt>*this</tt> will represent.
	 * @param  vs     [in] Stride for the storage in <tt>v[]</tt> (see Postconditions).
	 * @param  dim    [in] Number of elements in <tt>v[]</tt> this this will represent (see Postconditions).
	 * @param  ownsMem
	 *                [in] If <tt>true</tt> then <tt>delete [] v</tt> will be called after <tt>*this</tt>
	 *                no longer needs access to this memory.  If <tt>false</tt> then the client is responsible
	 *                for memory management.  The default is <tt>false</tt>.
	 * @param  vecSpc
	 *                [in] Smart pointer to a <tt>VectorSpace</tt> object that will be used to represent the
	 *                vector space for <tt>*this</tt>.  If <tt>vecSpc.get()==NULL</tt> on input, then
	 *                a <tt>SerialVectorSpace</tt> object of dimension <tt>dim</tt> is allocated for this
	 *                role.  The default is <tt>Teuchos::null</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> [<tt>vecSpc.get()!=NULL</tt>] <tt>dim == vecSpc->dim()</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> [<tt>vecSpc.get()!=NULL</tt>] <tt>vecSpc->createMember()</tt> must create vectors that are compatible
	 *      with <tt>*this</tt> (i.e. <tt>getSubVector()</tt>, <tt>commitSubVector()</tt> behave the same as with
	 *      this class).
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>vecSpc.get()!=NULL</tt>] <tt>vecSpc.get() == this->space().get()</tt>
	 * <li> [<tt>vecSpc.get()==NULL</tt>] <tt>dynamic_cast<const SerialVectorSpace<Scalar>*>(this->space().get()) != NULL</tt>
	 * <li> <tt>this->space()->dim() == dim</tt>
	 * <li> <tt>this->getPtr() == v</tt>
	 * <li> <tt>this->getStride() == vs</tt>
	 * <li> <tt>this->getOwnsMem() == ownsMem</tt>
	 * </ul>
	 *
	 * Note that this function is declared virtual so that subclasses
	 * can override it to be informed whenever <tt>*this</tt> vector
	 * is resized.  An override should call this function as
	 * <tt>this->SerialVector<Scalar>::initialize(...)</tt>.
	 */
	virtual void initialize(
		Scalar  v[]
		,int    vs
		,int    dim
		,bool   ownsMem = false
		,const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vecSpc = Teuchos::null
		);

	//@}

	/** @name Accessors (inlined for minimal overhead) */
	//@{

	///
	Scalar* getPtr();
	///
	const Scalar* getPtr() const;
	///
	int getStride() const;
	///
	int getDim() const;
	///
	bool getOwnsMem() const;
	
	//@}

	/** @name Overridden from Vector */
	//@{

	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > space() const;
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
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >   space_serial_;

	// ////////////////////////////////
	// Private member functions

	void free_mem();

	// Not defined and not to be called
	SerialVector(const SerialVector&);
	SerialVector& operator=(const SerialVector&);

}; // end class SerialVector

// /////////////////////////////////////////////////////
// Inline members

template<class Scalar>
inline
Scalar* SerialVector<Scalar>::getPtr()
{
	return v_;
}

template<class Scalar>
inline
const Scalar* SerialVector<Scalar>::getPtr() const
{
	return v_;
}

template<class Scalar>
inline
int SerialVector<Scalar>::getStride() const
{
	return vs_;
}	

template<class Scalar>
inline
int SerialVector<Scalar>::getDim() const
{
	return dim_;
}	

template<class Scalar>
inline
bool SerialVector<Scalar>::getOwnsMem() const
{
	return ownsMem_;
}	

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SERIAL_DECL_HPP
