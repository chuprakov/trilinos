// ////////////////////////////////////////////////////////////////////////////
// TSFCoreMPIVectorSpaceBaseDecl.hpp

#ifndef TSFCORE_MPI_VECTOR_SPACE_BASE_DECL_HPP
#define TSFCORE_MPI_VECTOR_SPACE_BASE_DECL_HPP

#include "TSFCoreVectorSpace.hpp"

namespace TSFCore {

///
/** <tt>%VectorSpace</tt> node subclass for all MPI-based vectors with
 * contiguous local storage.
 *
 * This interface collaborates with the <tt>MPIVetorBase</tt> class
 * to implement MPI-based parallel (or serial of course) vectors
 * that allows different concrete implementations to mix.
 *
 * This interface provides all the information necessary to implement
 * <tt>MPIVectorBase::applyOp()</tt>.  This interface returns an MPI
 * communicator (of which all compatible vector spaces must have the
 * same communicator obviously) through the method <tt>mpiComm()</tt>.
 *
 * <A NAME="MPIVectorSpaceBase_Vector_layout"></A>
 * <b>%Vector data layout:</b>
 *
 * This interface base class assumes that vector data is partitioned
 * to processors in contiguous chunks and dense subvectors.  To spell
 * this out, let <tt>v</tt> be the local vector that is sorted on this
 * processor and let <tt>g</tt> be the global vector.  Then these two
 * vectors are related (using one-based indexing) as:
 *
 * <tt>v(k) == g(k + this->localOffset()), for k = 1...this->localSubDim()</tt>
 *
 * Any type of mapping of vector data to processors that can not be
 * interpreted in this way can not rely on this base class
 * for interoperability.  In the case that it is found at runtime
 * that the mapping of vector elements to processors is not
 * as shown above, then the method <tt>mapCode()</tt> should be
 * overridden to return an invalid code.  It will then be required
 * that the subclass also override the <tt>isCompatible()</tt> method.
 * if this is needed, then vector interoperability between MPI-based
 * vectors will not be possible.
 *
 * <b>Notes to subclass developers:</b>
 *
 * The pure virtual methods <tt>mpiComm()</tt>, <tt>localOffset()</tt>
 * and <tt>localSubDim()</tt> defined in this interface along with the
 * pure virtual methods <tt>dim()</tt> and <tt>createMember()</tt> are
 * the only methods that must be overridden.
 *
 * If it is possible that the mapping of vector elements to processors
 * is not as described above, then the subclass should override the
 * <tt>mapCode()</tt> and <tt>isCompatible()</tt> methods as described
 * above and below.
 *
 * If optimized implementations of multi-vectors can be supported,
 * then the <tt>createMembers()</tt> method should also be overridden.
 */
template<class Scalar>
class MPIVectorSpaceBase : virtual public VectorSpace<Scalar> {
public:

	///
	MPIVectorSpaceBase();

	/** @name Pure virtual methods to be overridden by subclasses */
	//@{

	///
	/** Returns the MPI communicator.
	 */
	virtual MPI_Comm mpiComm() const = 0;
	///
	/** Returns the offset for the local sub-vector stored on this
	 * processor.
	 */
	virtual Index localOffset() const = 0;
	///
	/** Returns the number of local elements stored on this processor.
	 */
 	virtual Index localSubDim() const = 0;

	//@}

	/** @name Virtual methods with default implementations */
	//@{

	///
	/** Invalidate the internal state in case things have changed.
	 *
	 * This method must be called whenever the values
	 * <tt>mpiComm()</tt>, <tt>numProc</tt> (where <tt>numProc</tt> is
	 * returned from <tt>MPI_Comm_size(this->mpiComm(),&numProc)</tt>,
	 * <tt>localOffset()</tt> or <tt>localSubDim()</tt> change.
	 * 
	 * This is method is required since <tt>*this</tt> object
	 * maintains some cached data that is determined by the above
	 * data.  After ths method is called the cached data will
	 * automatically be recomputed when it is needed (i.e. by a call
	 * to <tt>mapCode()</tt> or <tt>isInCore()</tt>).
	 *
	 * Subclass developers: play nice and call this method when any of
	 * the above data changes please!
	 */
	virtual void invalidateState();
	///
	/** Returns the code for the mapping of elements to processors.
	 *
	 * This method takes the data <tt>mpiComm()</tt>, <tt>numProc</tt>
	 * (where <tt>numProc</tt> is returned from
	 * <tt>MPI_Comm_size(this->mpiComm(),&numProc)</tt>,
	 * <tt>localOffset()</tt> or <tt>localSubDim()</tt> on each
	 * processor and then uses it to compute a value for
	 * <tt>mapCode</tt> (using a single global reduction if
	 * <tt>numProc > 1</tt>) which is returned from this function.
	 *
	 * The value returned from this default implementation of this
	 * method must not be changed or this approach breaks down.  The
	 * only reason for overridding this method is for the subclass to
	 * be alerted of <em>when</em> this method is called but not
	 * <em>what</em> is returned from this method.  If a subclass
	 * developer does not understand what this means then <b>don't</b>
	 * override this method!
	 *
	 * The default implementation will always return <tt>return >
	 * 0</tt> so that if this method is overriden to return <tt>return
	 * <= </tt> then this is a flag that the underlying vector map
	 * does not satisfy the assumptions of this vector space interface
	 * and vectors that are in <tt>*this</tt> vector space can not
	 * collaborate with other MPI-based vector implementations.
	 */
	virtual Index mapCode() const;

	//@}

	/** @name Overridden from VectorSpace */
	//@{

	///
	/** Returns true if all of the elements are stored on one processor.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return == (numProc==1)</tt>, where <tt>numProc</tt> is
	 *      returned from <tt>MPI_Comm_size(this->mpiComm(),&numProc)</tt>.
	 * </ul>
	 */
 	bool isInCore() const;
	///
	/** Checks to general compatibility of parallel (or serial on one
	 * processor) MPI-based vector spaces.
	 *
	 * @return Returns true if <tt>*this</tt> and <tt>vecSpace</tt>
	 * are both serial in-core vectors or if <tt>vecSpc</tt> is of
	 * type <tt>MPIVectorSpaceBase<Scalar></tt> and both <tt>*this</tt>
	 * and <tt>vecSpc</tt> have the same MPI communicators and the same
	 * mapping of elements to processors.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return = ( this->isInCore() &&
	 *      vecSpc->isInCore() ) || ( (mpiVecSpc = dynamic_cast<const
	 *      MPIVectorSpaceBase<Scalar>*>(&vecSpc)) && this->mpiComm() ==
	 *      mpiVecSpc->mpiComm() && this->mapCode() ==
	 *      mpiVecSpc->mpiComm())</tt>.
	 *
	 * </ul>
	 *
	 * If the mapping of vector elements to processors is not as
	 * described
	 * <A HREF="classTSFCore_1_1MPIVectorSpaceBase.html#MPIVectorSpaceBase_Vector_layout>above</A>
	 * then this method should be overridden in a way that is specific
	 * to the vector implementation.
	 */
 	bool isCompatible(const VectorSpace<Scalar>& vecSpc) const;
	
	//@}

private:

	// //////////////////////////////////////
	// Private data members

	mutable Index     mapCode_;    // < 0 is a flag that everything needs initialized
	mutable bool      isInCore_;

	// //////////////////////////////////////
	// Private member functions

	void updateState() const;
	
}; // end class MPIVectorSpaceBase

} // end namespace TSFCore

#endif // TSFCORE_MPI_VECTOR_SPACE_BASE_DECL_HPP
