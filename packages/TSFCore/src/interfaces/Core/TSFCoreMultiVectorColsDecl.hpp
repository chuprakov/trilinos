// /////////////////////////////////////////////////////////////////////////////
// TSFCoreMultiVectorColsDecl.hpp

#ifndef TSF_MULTI_VECTOR_COlS_DECL_HPP
#define TSF_MULTI_VECTOR_COlS_DECL_HPP

#include "TSFCore_ConfigDefs.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"

namespace TSFCore {

///
/** Default subclass for <tt>MultiVector</tt> implemented using columns
 * of separate abstract vectors.
 *
 * This is a very bad implementation of a multi-vector but this will
 * work in situations where you need a multi-vector but some
 * underlying linear algebra library does not directly support them.
 *
 * This subclass can be used to represent a <tt>%MultiVector</tt>
 * wrapper around a single <tt>Vector</tt> object so that a single
 * vector can be passed to a method that expects a <tt>%MultiVector</tt>
 * object.
 */
template<class Scalar>
class MultiVectorCols : virtual public MultiVector<Scalar> {
public:

	/** @name Constructors/Initializers */
	//@{

	///
	/** Construct to initialized.
	 *
	 * Postconditions:<ul>
	 * <tt> <tt>this->range().get() == NULL</tt>
	 * <tt> <tt>this->domain().get() == NULL</tt>
	 * </ul>
	 */
	MultiVectorCols();

	/// Calls <tt>initalize()</tt>.
	MultiVectorCols(
		MemMngPack::ref_count_ptr<Vector<Scalar> >  col_vec
		);

	/// Calls <tt>initalize()</tt>.
	MultiVectorCols(
		const MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >          &range
		,const MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >         &domain
		,MemMngPack::ref_count_ptr<Vector<Scalar> >                          col_vecs[] = NULL
		);
	
	///
	/** Initialize given a single vector object.
	 *
	 * @param  col_vec  [in] A single column vector.  It is not allowed for
	 *                  <tt>col_vecs==NULL</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>col_vec.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>col_vec->dim() > 0</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <tt> <tt>this->range().get() == col_vec.space().get()</tt>
	 * <tt> <tt>this->domain()->dim() == 1</tt>
	 * <li> <tt>this->col(1).get() == col_vec.get()</tt>
	 * </ul>
	 */
	void initialize(
		MemMngPack::ref_count_ptr<Vector<Scalar> >  col_vec
		);

	///
	/** Initialize given the spaces for the columns and rows and possibly the column vectors.
	 *
	 * @param  range    [in] The space that the columns must lie in.  The underlying
	 *                  vector space must not be changed while <tt>this</tt> is in use.
	 * @param  domain   [in] The space that the rows must lie in.  The underlying
	 *                  vector space must not be changed while <tt>this</tt> is in use.
	 *                  What this argument really specifies is what vector type
	 *                  will be compatible with the vectors that the client may
	 *                  try to use to interact with the rows of this multivector.
	 * @param  col_vecs [in] Array (size <tt>domain->dim()</tt>) of column
	 *                  vectors to use for the columns of <tt>this</tt>.
	 *                  It is allowed for <tt>col_vecs==NULL</tt> in which case
	 *                  <tt>range->createMember()</tt> will be used to
	 *                  create the colmns of <tt>this</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>range.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>domain.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>range->dim() > 0</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>domain->dim() > 0</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> [<tt>col_vecs != NULL</tt>]
	 *      <tt>col_vecs[j-1].get() != NULL && col_vecs[j-1]->space()->is_compatible(*range) == true</tt>,
	 *      for <tt>j=1..domain->dim()</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <tt> <tt>this->range().get() == range.get()</tt>
	 * <tt> <tt>this->domain().get() == domain.get()</tt>
	 * <li> [<tt>col_vecs != NULL</tt>] <tt>this->col(j).get() == col_vecs[j-1].get()</tt>,
	 *      for <tt>j=1..domain->dim()</tt>
	 * </ul>
	 */
	void initialize(
		const MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >          &range
		,const MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >         &domain
		,MemMngPack::ref_count_ptr<Vector<Scalar> >                          col_vecs[] = NULL
		);

	/// Set uninitalized.
	void set_uninitialized();

	//@}

	/** @name Overridden from LinearOp */
	//@{
	///
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> > range() const;
	///
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> > domain() const;
	//@}

	/** @name Overridden from MultiVector */
	//@{
	///
	MemMngPack::ref_count_ptr<Vector<Scalar> > col(Index j);
	///
	MemMngPack::ref_count_ptr<MultiVector<Scalar> > subView( const Range1D& col_rng );
	//@}

private:
	
#ifdef DOXYGEN_COMPILE
	const VectorSpace<Scalar>                  *range;
	const VectorSpace<Scalar>                  *domain;
	Vector<Scalar>                             *column_vectors;
#else
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >        range_;
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >        domain_;
	std::vector< MemMngPack::ref_count_ptr<Vector<Scalar> > >    col_vecs_;
#endif
	
}; // end class MultiVectorCols

} // end namespace TSFCore

#endif // TSF_MULTI_VECTOR_COlS_DECL_HPP
