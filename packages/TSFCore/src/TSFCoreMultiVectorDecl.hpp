// ///////////////////////////////////////////////////////////
// MultiVectorDecl.hpp

#ifndef TSFCORE_MULTI_VECTOR_DECL_HPP
#define TSFCORE_MULTI_VECTOR_DECL_HPP

#include "TSFCoreLinearOpDecl.hpp"
#include "RTOpCpp.hpp"

namespace TSFCore {

///
/** Interface for a collection of column vectors called a multi-vector.
 *
 * The primary purpose for this interface is to allow for convienent
 * aggregations of column vectors.  Such an orderly arrangement allows
 * for better optimized linear algebra operations such as
 * matrix-matrix multiplication and the solution of linear systems for
 * multiple right hand sides.  Every application area (serial
 * parallel, out-of-core etc.) should be able to define at least one
 * reasonbly efficient implementation of <tt>%MultiVector</tt>.
 *
 * The <tt>%MultiVector</tt> interface is derived from the
 * <tt>LinearOp</tt> interface and therefore a <tt>%MultiVector</tt>
 * can be considered as a linear operator which has some interesting
 * implications.
 *
 * Note that since, this interface is derived from <tt>LinearOp</tt>
 * that it must support the methods <tt>domain()<tt> and
 * <tt>range()</tt>.
 *
 * Another very powerful feature of this interface is the ability to
 * apply reduction/transformation operators over a sub-set of rows and
 * columns in a set of multi-vector objects.  The behavior is
 * identical as if the client extracted the rows or columns in a set
 * of multi-vectors and called <tt>Vector::applyOp()</tt> itself.
 * However, the advantage of using the multi-vector methods is that
 * there may be greater opportunities to exploit parallelism.  Also,
 * the intermediate reduction objects over a set of rows or columns
 * can be reduced by a secondary reduction object.
 *
 * <b>Notes for subclass developers</b>
 *
 * Only one method override is required (in addition to the methods
 * that must be overridden from <tt>LinearOp</tt>, except
 * <tt>apply()</tt>, see below) to create a concrete
 * <tt>MultiVector<Scalar></tt> subclass which is the non-constant
 * version of <tt>col()</tt>.  All of the other methods have default
 * implementations that are quite good in most cases.
 *
 * Note that through the magic of the <tt>applyOp()</tt> methods that
 * this interface is able to implement the pure virtual
 * <tt>apply()</tt> method from the <tt>LinearOp</tt> interface.  This
 * implementation is not optimal, but will be sufficient in many
 * different contexts.
 *
 * The <tt>applyOp()</tt> methods should only be overridden if the
 * subclass can do something more efficient than simply applying the
 * reduction/transformation operators one column at a time.
 *
 * The non-const versions of <tt>subView()</tt> should be overridden
 * if the performance of these methods are important.  The default
 * implementation will not achieve near-optimal performance.
 *
 * Methods that subclasses should almost never need (or want) to
 * override are the const version of <tt>col()</tt>, the const version
 * of <tt>clone_mv()</tt> or the the const versions of
 * <tt>subView()</tt>.
 */
template<class Scalar>
class MultiVector : virtual public LinearOp<Scalar> {
public:
	
	/** @name Provide access to the columns as <tt>Vector</tt> objects */
	//@{

	///
	/** Get a non-mutable view of column vector.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>1 <= j && j <= this->range().dim()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get() != NULL</tt>
	 * <li> <tt>this->range()->isCompatible(*return->space()) == true</tt>
	 * </ul>
	 *
	 * The default implementation of this method (which is the only
	 * implementaion needed by any subclass) is based on the non-const
	 * version <tt>col()</tt>.
	 */
	virtual MemMngPack::ref_count_ptr<const Vector<Scalar> > col(Index j) const;

	///
	/** Get a mutable view of column vector.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>1 <= j && j <= this->range()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get() != NULL</tt>
	 * <li> <tt>this->range()->isCompatible(*return->space()) == true</tt>
	 * </ul>
	 */
	virtual MemMngPack::ref_count_ptr<Vector<Scalar> > col(Index j) = 0;

	//@}

	/** @name Cloning */
	//@{

	///
	/** Clone the multi-vector object (if supported).
	 *
	 * The default implementation returns the object from
	 * the non-const version of <tt>clone_mv()</tt>.
	 */
	virtual MemMngPack::ref_count_ptr<const MultiVector<Scalar> > clone_mv() const;

	///
	/** Clone the multi-vector object (if supported).
	 *
	 * The default implementation uses the vector space to create a
	 * new multi-vector object and then uses a transformation operator
	 * to assign the vector elements.  A subclass should only override
	 * this method if it can do something more sophisticated
	 * (i.e. lazy evaluation) but in general, this is not needed.
	 */
	virtual MemMngPack::ref_count_ptr<MultiVector<Scalar> > clone_mv();

	//@}

	/** @name Sub-view methods */
	//@{

	///
	/** Returns a const sub-view of a contiguous set of columns of the this multi-vector.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>col_rng.ubound() <= this->range()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->range()->isCompatible(*return->range()) == true</tt>
	 * <li> <tt>return->domain()->dim() == RangePack::full_range(col_rng,1,this->domain()->dim()).size()</tt>
	 * <li> <tt>*return->col(1+k)</tt> represents the same column vector as <tt>this->col(col_rng.lbound()+k)</tt>,
	 *      for <tt>k=0...RangePack::full_range(col_rng,1,this->domain()->dim()).ubound()</tt>
	 * </ul>
	 *
	 * The default implementation (which is the only implementation needed by subclasses)
	 * is to return object from the non-const verstion <tt>subView()</tt>.
	 */
	virtual MemMngPack::ref_count_ptr<const MultiVector<Scalar> > subView( const Range1D& col_rng ) const;
	
	///
	/** Returns a non-const sub-view of a contiguous set of columns of the this multi-vector.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>col_rng.ubound() <= this->range()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->range()->isCompatible(*return->range()) == true</tt>
	 * <li> <tt>return->domain()->dim() == RangePack::full_range(col_rng,1,this->domain()->dim()).size()</tt>
	 * <li> <tt>*return->col(1+k)</tt> represents the same column vector as <tt>this->col(col_rng.lbound()+k)</tt>,
	 *      for <tt>k=0...RangePack::full_range(col_rng,1,this->domain()->dim()).ubound()</tt>
	 * </ul>
	 *
	 * The default implementation of this function uses <tt>MultiVectorCols</tt> but this is not
	 * a good default implementation in general.
	 */
	virtual MemMngPack::ref_count_ptr<MultiVector<Scalar> > subView( const Range1D& col_rng );

	///
	/** Returns a const sub-view of a non-contiguous set of columns of the this multi-vector.
	 *
	 * @param  numCols  [in] The number of columns to extract a view for.
	 * @param  cols[]   [in] Array (length <tt>numCols</tt>) of the columns indices to use in the
	 *                  returned view.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>numCols <= this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>1 <= cols[k] <= this->domain()->dim()</tt>, for <tt>k=0...numCols-1</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>col[k1] != col[k2]</tt>, for all <tt>k1 != k2</tt> in the range <tt>[0,numCols-1]</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->range()->isCompatible(*return->range()) == true</tt>
	 * <li> <tt>return->domain()->dim() == numCols</tt>
	 * <li> <tt>*return->col(k+1)</tt> represents the same column vector as <tt>this->col(cols[k])</tt>,
	 *      for <tt>k=0...numCols</tt>
	 * </ul>
	 *
	 * The default implementation (which is the only implementation needed by subclasses)
	 * is to return object from the non-const verstion <tt>subView()</tt>.
	 */
	virtual MemMngPack::ref_count_ptr<const MultiVector<Scalar> > subView( const int numCols, const int cols[] ) const;


	///
	/** Returns a non-const sub-view of a non-contiguous set of columns of the this multi-vector.
	 *
	 * @param  numCols  [in] The number of columns to extract a view for.
	 * @param  cols[]   [in] Array (length <tt>numCols</tt>) of the columns indices to use in the
	 *                  returned view.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>1 <= numCols <= this->domain()->dim()</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>1 <= cols[k] <= this->domain()->dim()</tt>, for <tt>k=0...numCols-1</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>col[k1] != col[k2]</tt>, for all <tt>k1 != k2</tt> in the range <tt>[0,numCols-1]</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->range()->isCompatible(*return->range()) == true</tt>
	 * <li> <tt>return->domain()->dim() == numCols</tt>
	 * <li> <tt>*return->col(k+1)</tt> represents the same column vector as <tt>this->col(cols[k])</tt>,
	 *      for <tt>k=0...numCols</tt>
	 * </ul>
	 *
	 * The default implementation of this function uses <tt>MultiVectorCols</tt> but this is not
	 * a good default implementation in general.
	 */
	virtual MemMngPack::ref_count_ptr<MultiVector<Scalar> > subView( const int numCols, const int cols[] );
	
	//@}

	/** @name Collective applyOp() methods */
	//@{

	///
	/** Apply a reduction/transformation operator column by column and
	 * return an array of the reduction objects.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <tt> See the preconditions for <tt>TSFCore::applyOp()</tt>
	 * </ul>
	 *
	 * See the documentation for the method <tt>TSFCore::applyOp()</tt>
	 * for a description of the arguments.
	 *
	 * This method is not to be called directly by the client but instead
	 * through the nonmember function <tt>TSFCore::applyOp()</tt>.
	 *
	 * It is expected that <tt>this</tt> will be one of the multi-vector
	 * objects in <tt>multi_vecs[]</tt> or <tt>targ_multi_vecs[]</tt>.
	 *
	 * The default implementation calls <tt>Vector::applyOp()</tt> on
	 * each colunn <tt>this->col(j)</tt> for <tt>j = 1
	 * ... this->range()->dim()</tt>.
	 */
	virtual void applyOp(
		const RTOpPack::RTOpT<Scalar>   &primary_op
		,const size_t                   num_multi_vecs
		,const MultiVector<Scalar>*     multi_vecs[]
		,const size_t                   num_targ_multi_vecs
		,MultiVector<Scalar>*           targ_multi_vecs[]
		,RTOp_ReductTarget              reduct_objs[]
		,const Index                    primary_first_ele
		,const Index                    primary_sub_dim
		,const Index                    primary_global_offset
		,const Index                    secondary_first_ele
		,const Index                    secondary_sub_dim
		) const;

	///
	/** Apply a reduction/transformation operator column by column and reduce the intermediate
	 * reduction objects into one reduction object.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <tt> See the preconditions for <tt>TSFCore::applyOp()</tt>
	 * </ul>
	 *
	 * See the documentation for the method <tt>TSFCore::applyOp()</tt>
	 * for a description of the arguments.
	 *
	 * This method is not to be called directly by the client but instead
	 * through the nonmember function <tt>TSFCore::applyOp()</tt>.
	 *
	 * It is expected that <tt>this</tt> will be one of the multi-vector
	 * objects in <tt>multi_vecs[]</tt> or <tt>targ_multi_vecs[]</tt>.
	 *
	 * The default implementation calls <tt>applyOp()</tt> where an
	 * array of reduction objects is taken.
	 */
	virtual void applyOp(
		const RTOpPack::RTOpT<Scalar>   &primary_op
		,const RTOpPack::RTOpT<Scalar>  &secondary_op
		,const size_t                   num_multi_vecs
		,const MultiVector<Scalar>*     multi_vecs[]
		,const size_t                   num_targ_multi_vecs
		,MultiVector<Scalar>*           targ_multi_vecs[]
		,RTOp_ReductTarget              reduct_obj
		,const Index                    primary_first_ele
		,const Index                    primary_sub_dim
		,const Index                    primary_global_offset
		,const Index                    secondary_first_ele
		,const Index                    secondary_sub_dim
		) const;

	//@}

	/** @name Overridden methods from LinearOp */
	//@{

	///
	/** This method is implemented in terms of the mulit-vector <tt>applyOp()</tt> method.
	 *
	 * The implementation takes care of two types of operations.  One
	 * (<tt>M_trans==TRANS</tt>) is the block dot product of two
	 * vectors to form scalar (stored as the vector <tt>y</tt> which
	 * as one component).  The other (<tt>M_trans==NOTRANS</tt>) is
	 * essentially an axpy operation where <tt>x</tt> is a vector with
	 * one element.  Both of these operations are performed using
	 * reduction/transformation operators.
	 *
	 * This implementation is near optimal but the default
	 * implementation of the multi-vector verstion of <tt>apply()</tt>
	 * as implemented in the base class <tt>LinearOp</tt> will not be
	 * a near optimal implementation in its current form do to
	 * multiple, sequential reductions but it could be made to be with
	 * a little work.
	 */
	void apply(
		const ETransp            M_trans
		,const Vector<Scalar>    &x
		,Vector<Scalar>          *y
		,const Scalar            alpha
		,const Scalar            beta
		) const;

	/// This method is simply overridden to return <tt>this->clone_lons()</tt>.
	MemMngPack::ref_count_ptr<const LinearOp<Scalar> > clone() const;

	//@}

private:
	
#ifdef DOXYGEN_COMPILE
	Vector<Scalar>  *columns; // doxygen only
#endif	

}; // end class MultiVector<Scalar>

///
/** Apply a reduction/transformation operator column by column and
 * return an array of the reduction objects.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
inline
void applyOp(
	const RTOpPack::RTOpT<Scalar>   &primary_op
	,const size_t                   num_multi_vecs
	,const MultiVector<Scalar>*     multi_vecs[]
	,const size_t                   num_targ_multi_vecs
	,MultiVector<Scalar>*           targ_multi_vecs[]
	,RTOp_ReductTarget              reduct_objs[]
	,const Index                    primary_first_ele      = 1
	,const Index                    primary_sub_dim        = 0
	,const Index                    primary_global_offset  = 0
	,const Index                    secondary_first_ele    = 1
	,const Index                    secondary_sub_dim      = 0
	)
{
	if(num_multi_vecs)
		multi_vecs[0]->applyOp(
			primary_op
			,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
			,reduct_objs,primary_first_ele,primary_sub_dim,primary_global_offset
			,secondary_first_ele,secondary_sub_dim
			);
	else if(num_targ_multi_vecs)
		targ_multi_vecs[0]->applyOp(
			primary_op
			,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
			,reduct_objs,primary_first_ele,primary_sub_dim,primary_global_offset
			,secondary_first_ele,secondary_sub_dim
			);
}

///
/** Apply a reduction/transformation operator column by column and reduce the intermediate
 * reduction objects into one reduction object.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
inline
void applyOp(
	const RTOpPack::RTOpT<Scalar>   &primary_op
	,const RTOpPack::RTOpT<Scalar>  &secondary_op
	,const size_t                   num_multi_vecs
	,const MultiVector<Scalar>*     multi_vecs[]
	,const size_t                   num_targ_multi_vecs
	,MultiVector<Scalar>*           targ_multi_vecs[]
	,RTOp_ReductTarget              reduct_obj
	,const Index                    primary_first_ele      = 1
	,const Index                    primary_sub_dim        = 0
	,const Index                    primary_global_offset  = 0
	,const Index                    secondary_first_ele    = 1
	,const Index                    secondary_sub_dim      = 0
	)
{
	if(num_multi_vecs)
		multi_vecs[0]->applyOp(
			primary_op,secondary_op
			,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
			,reduct_obj,primary_first_ele,primary_sub_dim,primary_global_offset
			,secondary_first_ele,secondary_sub_dim
			);
	else if(num_targ_multi_vecs)
		targ_multi_vecs[0]->applyOp(
			primary_op,secondary_op
			,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
			,reduct_obj,primary_first_ele,primary_sub_dim,primary_global_offset
			,secondary_first_ele,secondary_sub_dim
			);
}

} // namespace TSFCore

#endif // TSFCORE_MULTI_VECTOR_DECL_HPP
