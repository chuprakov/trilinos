#ifndef TSFVECTOR_H
#define TSFVECTOR_H

#include "TSFDefs.h"
#include <typeinfo>
#include "TSFSmartPtr.h"
#include "TSFTimer.h"
#include "TSFTimeMonitor.h"
#include "TSFRandomNumberGenerator.h"
#include "SystemRand.h"
#include "TSFDeferredCopy.h"
#include "TSFArray.h"
#include "TSFVectorBase.h"
#include <iostream>
#include <string>

namespace TSF
{
  class DenseSerialVector;
  class TSFDeferredLinearCombination;

  using std::string;
  using std::ostream;


  /** \ingroup Vector
   * TSFVector:  User-level handle class for vectors.
   *
   * <b> Constructing a TSFVector </b>
   *
   * Ordinarily, you will never construct a TSFVector directly from a derived type.
   * Rather, the createMember() method of TSFVectorSpace is used to build a vector
   * of the appropriate type, for example,
   * \code
   * TSFVector x = mySpace.createMember();
   * \endcode
   * This hides from you all the ugly details of creating a particular concrete type.
   *
   * You will frequently create an empty vector to be filled in later, for example,
   * \code
   * TSFVector y;
   * \endcode
   * Note that this vector isn't just empty, it's null. Not only does it have no
   * values assigned, it does not have a concrete type. An attempt to set an element
   * of or operate on a null vector will result in an error. What you can do is
   * assign another vector to it,
   * \code
   * TSFVector y;
   * y = x.deepCopy();
   * \endcode
   * or fill it by passing it as a return-by-reference argument to a mathematical
   * operation,
   * \code
   * TSFVector y;
   * x.scalarMult(3.14159, y);
   * \endcode
   *
   *
   * Another way a TSFVector can be created is as the result of an operation between
   * existing vectors, for example
   * \code
   * TSFVector z = a*x + b*y;
   * \endcode
   * See the section on operator overloading for more information on this.
   *
   *
   * <b> Element access </b>
   *
   * There are a number of methods to get or set the value of an element or group
   * of elements at a particular index or group of indices. In a distributed environment,
   * there is a distinction between an element's <i> local </i> index and its
   * <i> global </i> index. All the element access methods of TSFVector use global
   * indices.
   *
   * Element access methods must make virtual function calls, so they
   * should be used as sparingly as possible. In particular, you
   * should <i> never </i> do mathematical operations on the whole
   * vector using the element access methods. Use one of the
   * predefined operations (e.g., eMult, update) or a custom reduction
   * and transformation operator instead. If you are needing to set
   * elements, such as when creating the load vector in a
   * finite-element problem, it is more efficient to set groups of
   * elements rather than single elements.
   *
   * <b> Overloaded operators </b>
   *
   * TSFVector supports the full set of overloaded operators for vector and scalar-vector
   * arithmetic. Overloaded operators are not often used in scientific computing
   * because a naive implementation incurs a large performance penalty due to
   * the creation and copying of temporary vectors. TSFVector operators use a design
   * which avoids the creation of such temporary vectors by flattening the expression
   * tree. Small temporary objects are created in the process, so for small problems
   * (under \f$N \approx 1000 \f$) there is still a noticeable performance penatly.
   * However, for larger problems the performance penalty is less than a few percent
   * compared to FORTRAN BLAS, and decreases quickly with problem size.
   * */

  class TSFVector
    {
    public:
      /** \name Constructors, Destructors, and Assignment Operators */
      //@{
      /** empty ctor. Constructs a null vector */
      TSFVector();

      /** Construct with a pointer to a derived vector type. */
      TSFVector(TSFVectorBase* ptr);

      /** construct with the result of an operation */
      TSFVector(const TSFDeferredLinearCombination& op);

      /** construct with the result of an operation */
      TSFVector(const TSFDeferredCopy& copy);

      /** assign the result of an operation into this vector */
      TSFVector& operator=(const TSFDeferredLinearCombination& op);

      /** assign the result of an operation into this vector */
      TSFVector& operator=(const TSFDeferredCopy& copy);
      //@}


      /** \name Vector space support */
      //@{
      /** return the vector space in which this vector lives. */
      const TSFVectorSpace& space() const ;

      /** indicate whether this vector is a member of a given space */
      bool isMemberOf(const TSFVectorSpace& space) const ;
      //@}

      /** \name Block access */
      //@{

      /** returns true if BlockVector  */
      bool isBlockVector() const;


      /** return the number of subvector blocks */
      int numBlocks() const ;

      /** return the i-th subvector */
      TSFVector getBlock(int i) const ;

      /** set the i-th subvector */
      TSFVector& setBlock(int i, const TSFVector& sub);

      //@}

#if HAVE_RTOP

      /** \name Reduction/Transformation operators. */
      //@{

      ///
        /** Apply a reduction/transformation,operation over a set of vectors:
         * <tt>op(op((*this),v[0]...v[nv-1],z[0]...z[nz-1]),(*reduct_obj)) -> z[0]...z[nz-1],(*reduct_obj)</tt>.
         *
         * The first nonmutable vector in the argument list to <tt>op</tt> will be <tt>this</tt>
         * vector then followed by those in <tt>vecs[k], k = 0...num_vecs-1</tt>.  Therefore, the number
         * of nonmutable vectors passed to <tt>\ref RTOp_apply_op "apply_op"(&op,...)"</tt> will be
         * <tt>num_vecs+1</tt>.
         *
         * The vector to be represented <tt>v</tt> by <tt>this</tt> and passed to the
         * <tt>\ref RTOp_apply_op "RTOp_apply_op"(&op,...)</tt> method is:
         \verbatim

         v(k + global_offset) = this->get_ele(first_ele + k - 1)
         , for k = 1 ... sub_dim
         \endverbatim
         * The other vector arguments are represented similarly.  The situation where
         * <tt>first_ele == 1</tt> and <tt>global_offset > 1</tt> corresponds to the
         * case where the vectors are representing consitituent vectors in a larger
         * aggregrate vector.  The situation where <tt>first_ele > 1</tt> and
         * <tt>global_offset == 0</tt> is for when a sub-view of the vectors are being
         * treated as full vectors.  Other combinations of these arguments is possible.
         *
         * Preconditions:<ul>
         * <li> All of the input vectors is compatible with <tt>this</tt> vector object.
         * <li> <tt>1 <= first_ele <= this->dim()</tt>
         * <li> <tt>global_offset >= 0</tt>
         * <li> <tt>sub_dim - (first_ele - 1) <= this->dim()</tt>
         * </ul>
         *
         * @param  op [in] Reduction/transformation operator to apply over each sub-vector
         *        and assemble the intermediate targets into <tt>reduct_obj</tt> (if
         *              <tt>reduct_obj != RTOp_REDUCT_OBJ_NULL</tt>).
         * @param  num_vecs
         *        [in] Number of nonmutable vectors in <tt>vecs[]</tt>.
         *              If <tt>vecs==NULL</tt> then this argument is ignored but should be set to zero.
         * @param  vecs
         *        [in] Array (length <tt>num_vecs</tt>) of a set of pointers to
         *        nonmutable vectors to include in the operation.
         *        The order of these vectors is significant to <tt>op</tt>.
         *        If <tt>vecs==NULL</tt> then <tt>op</tt> is called with the
         *        single vector represented by <tt>this</tt> object.
         * @param  num_targ_vecs
         *        [in] Number of mutable vectors in <tt>targ_vecs[]</tt>.
         *              If <tt>targ_vecs==NULL</tt> then this argument is ignored but should be set to zero.
         * @param  targ_vecs
         *        [in] Array (length <tt>num_targ_vecs</tt>) of a set of pointers to
         *        mutable vectors to include in the operation.
         *        The order of these vectors is significant to <tt>op</tt>.
         *        If <tt>targ_vecs==NULL</tt> then <tt>op</tt> is called with no mutable vectors.
         * @param  reduct_obj
         *        [in/out] Target object of the reduction operation.
         *        This object must have been created by the <tt>RTOp_reduct_obj_create(&op,&reduct_obj)</tt>
         *              function first.  The reduction operation will be added to <tt>(*reduct_obj)</tt> if
         *              <tt>(*reduct_obj)</tt> has already been through a reduction.  By allowing the info in
         *              <tt>(*reduct_obj)</tt> to be added to the reduction over all of these vectors, the reduction
         *              operation can be accumulated over a set of abstract vectors which can be useful for implementing
         *              composite vectors instance.  If <tt>RTOp_get_reduct_type_num_entries(&op,...)</tt> returns
         *              <tt>num_values == 0</tt>, <tt>num_indexes == 0</tt> and <tt>num_chars == 0</tt> then
         *              <tt>reduct_obj</tt> should be set to <tt>RTOp_REDUCT_OBJ_NULL</tt> and no reduction will be performed.
         * @param  first_ele
         *        [in] (default = 1) The index of the first element in <tt>this</tt> to be included.
         * @param  sub_dim
         *              [in] (default = 0) The number of elements in these vectors to include in the reduction/transformation
         *              operation.  The value of <tt>sub_dim == 0</tt> means to include all available elements.
         * @param  global_offset
         *        [in] (default = 0) The offset applied to the included vector elements.
         */
        void apply_reduction(
                             const RTOp_RTOp &op, int num_vecs, const TSFVector* vecs[]
                             ,int num_targ_vecs, TSFVector* targ_vecs[], RTOp_ReductTarget reduct_obj
                             ,const RTOp_index_type first_ele = 1, const RTOp_index_type sub_dim = 0, const RTOp_index_type global_offset = 0
                             ) const;

        ///
          /** Apply a reduction/transformation,operation over a set of vectors:
           * <tt>op(op(v[0]...v[nv-1],(*this),z[0]...z[nz-1]),(*reduct_obj)) -> (*this),z[0]...z[nz-1],(*reduct_obj)</tt>.
           *
           * The first mutable vector in the argument list to <tt>op</tt> will be
           * <tt>this</tt> vector then followed by those in <tt>targ_vecs[k]</tt>, <tt>k = 0...num_targ_vecs-1</tt>
           * Therefore, the number of mutable vectors passed to <tt>\ref RTOp_apply_op "apply_op"(&op,...)"</tt>
           * will be <tt>num_targ_vecs+1</tt>.
           *
           * See <tt>\ref TSFVector::apply_reduction "apply_reduction(...)" for a discussion of the significance of the
           * arguments <tt>first_ele</tt>, <tt>sub_dim</tt> and <tt>global_offset</tt>.
           *
           * Preconditions:<ul>
           * <li> See <tt>\ref TSFVector::apply_reduction "apply_reduction(...)".
           * </ul>
           *
           * @param  op [in] Reduction/transformation operator to apply over each sub-vector
           *        and assemble the intermediate targets into <tt>reduct_obj</tt> (if
           *              <tt>reduct_obj != RTOp_REDUCT_OBJ_NULL</tt>).
           * @param  num_vecs
           *        [in] Number of nonmutable vectors in <tt>vecs[]</tt>.  If <tt>vecs==NULL</tt>
           *        then this argument is ignored but should be set to zero.
           * @param  vecs
           *        [in] Array (length <tt>num_vecs</tt>) of a set of pointers to
           *        nonmutable vectors to include in the operation.
           *        The order of these vectors is significant to <tt>op</tt>.  if <tt>vecs==NULL</tt>,
           *              then <tt>op.apply_op(...)</tt> is called with no non-mutable sub-vector arguments.
           * @param  num_targ_vecs
           *        [in] Number of mutable vectors in <tt>targ_vecs[]</tt>.  If <tt>targ_vecs==NULL</tt>
           *        then this argument is ignored but should be set to zero.
           * @param  targ_vecs
           *        [in] Array (length <tt>num_targ_vecs</tt>) of a set of pointers to
           *        mutable vectors to include in the operation. The order of these vectors
           *            is significant to <tt>op</tt>.  If <tt>targ_vecs==NULL</tt> then <tt>op</tt> is called with
           *        only one mutable vector (<tt>*this</tt>).
           * @param  reduct_obj
           *        [in/out] Target object of the reduction operation.
           *        This object must have been created by the <tt>op.reduct_obj_create_raw(&reduct_obj)</tt>
           *              function first.  The reduction operation will be added to <tt>(*reduct_obj)</tt> if
           *              <tt>(*reduct_obj)</tt> has already been through a reduction.  By allowing the info in
           *              <tt>(*reduct_obj)</tt> to be added to the reduction over all of these vectors, the reduction
           *              operation can be accumulated over a set of abstract vectors which can be useful for implementing
           *              composite vectors instance.  If <tt>op.get_reduct_type_num_entries<tt>(...)</tt> returns
           *              <tt>num_values == 0</tt>, <tt>num_indexes == 0</tt> and <tt>num_chars == 0</tt> then
           *              <tt>reduct_obj</tt> should be set to #RTOp_REDUCT_OBJ_NULL and no reduction will be performed.
           * @param  first_ele
           *        [in] (default = 1) The index of the first element in <tt>this</tt> to be included.
           * @param  sub_dim
           *              [in] (default = 0) The number of elements in these vectors to include in the reduction/transformation
           *              operation.  The value of <tt>sub_dim == 0</tt> means to include all available elements.
           * @param  global_offset
           *        [in] (default = 0) The offset applied to the included vector elements.
           */
          void apply_transformation(
                                    const RTOp_RTOp &op, int num_vecs, const TSFVector* vecs[]
                                    ,int num_targ_vecs, TSFVector* targ_vecs[], RTOp_ReductTarget reduct_obj
                                    ,const RTOp_index_type first_ele = 1, const RTOp_index_type sub_dim = 0, const RTOp_index_type global_offset = 0
                                    );

          //@}

#endif

          /** \name Setting, getting, and adding to individual elements or groups of elements. */
          //@{
          /**  Read-only access to an element specified by its global index  */
          const TSFReal& operator[](int globalIndex) const ;

          /** Read-write access to an element specified by its global index  */
          TSFReal& operator[](int globalIndex) ;

          /** pack selected elements into a vector */
          void getElements(const TSFArray<int>& globalIndices,
                           DenseSerialVector& sub) const ;

          /** set several elements */
          void setElements(const TSFArray<int>& globalIndices,
                           const DenseSerialVector& sub);

          /** add to several elements */
          void addToElements(const TSFArray<int>& globalIndices,
                             const DenseSerialVector& sub);

          /** pack selected elements into a vector */
          void getElements(const int* globalIndices, const int length,
                           DenseSerialVector& sub) const ;

          /** set several elements */
          void setElements(const int* globalIndices,
                           const DenseSerialVector& sub);

          /** add to several elements */
          void addToElements(const int* globalIndices,
                             const DenseSerialVector& sub);

          /** add to a selected element */
          void addToElement(int globalIndex, const TSFReal& val);
          //@}

          /** \name Math operations */
          //@{
          /** axpy (this = a*x + y) */
          void axpy(const TSFReal& a, const TSFVector& x, const TSFVector& y);

          /** multiplication by a scalar (this = a*x) */
          void scalarMult(const TSFReal& a, const TSFVector& x) ;

          /** addition (this = x + y) */
          inline void add(const TSFVector& x, const TSFVector& y) {axpy(1.0, x, y);}

          /** subtraction (this = x - y) */
          inline void subtract(const TSFVector& x, const TSFVector& y) {axpy(-1.0, y, x);}

          /** element-by-element multiplication (this = x .* y) */
          void dotStar(const TSFVector& x, const TSFVector& y) ;


          /** element-by-element division (this = x ./ y) */
          void dotSlash(const TSFVector& x, const TSFVector& y) ;

          /** dot product with another vector */
          TSFReal dot(const TSFVector& other) const ;

          /** 2-norm */
          TSFReal norm2() const ;

          /** 1-norm */
          TSFReal norm1() const ;

          /** infinity-norm */
          TSFReal normInf() const ;

          /** max value */
          TSFReal max() const ;

          /** min value */
          TSFReal min() const ;

          /** max with location */
          TSFReal max(TSFGeneralizedIndex& location) const ;

          /** min with location */
          TSFReal min(TSFGeneralizedIndex& location) const ;

          /** max less than tol with location */
          TSFReal max(const TSFReal& tol, TSFGeneralizedIndex& location) const ;

          /** min greater than tol with location */
          TSFReal min(const TSFReal& tol, TSFGeneralizedIndex& location) const ;

          /** element-wise absolute value */
          TSFVector abs() const ;

          /** set all elements to zero */
          void zero();

          /** sum all elements */
          TSFReal sumElements() const ;

          /** set all elements to a scalar value */
          void setScalar(const TSFReal& a);

          //@}

          /** \name generating random vectors */
          //@{
          /** Fill a vector with random elements */
          void randomize(const TSFRandomNumberGenerator& r = new SystemRand()) ;
          //@}

          /** \name copying */
          //@{
          /** make a deep copy of this vector. The execution of the copy is
           * deferred until assignment, so that we can test whether we need
           * to allocate space for the copy. */
          TSFDeferredCopy copy() const ;
          //@}

          /** \name output */
          //@{
          /** print, called by stream output */
          void print(ostream& os) const ;

          /** write a representation of this object to a string */
          string toString() const ;
          //@}


          /** \name introspection */
          //@{
          /** determine if a vector has not been initialized */
          bool isNull() const ;

          /** see if two vectors share the same pointer */
          bool isIdenticalTo(const TSFVector& other) const ;
          //@}

          /** \name Hooks for parallel support */
          //@{
          /** gather valid ghost values from other procs */
          void synchronizeGhostValues() const ;

          /** mark ghost values as invalid, meaning that they need to be
           * synchronized */
          void invalidateGhostValues() ;
          //@}

          /** \name dangerous developer-only methods. */
          //@{
          /** read-only pointer access */
          const TSFSmartPtr<TSFVectorBase>& smartPtr() const {return ptr_;}
          /** read-write pointer access */
          TSFSmartPtr<TSFVectorBase>& smartPtr() {return ptr_;}
          /** read-only pointer access */
          const TSFVectorBase* ptr() const {return &(*ptr_);}
          /** read-write pointer access */
          TSFVectorBase* ptr() {return &(*ptr_);}
          //@}

          /** \name internal math operations */
          //@{
          /** copy another's contents into self */
          void acceptCopyOf(const TSFVector& x);
          /** set self = self + a*x */
          void selfModifyingAxpy(const TSFReal& a, const TSFVector& x);
          /** set self = self .* x */
          void selfModifyingDotStar(const TSFVector& x);
          /** set self = self .* x */
          void selfModifyingDotSlash(const TSFVector& x);
          /** set self = a*self */
          void selfModifyingScalarMult(const TSFReal& a);
          /** find a min or max */
          TSFReal findExtremeValue(TSFVectorBase::MinOrMax type, TSFGeneralizedIndex& location,
                                   const TSFReal& tol) const ;

          //@}

          /** Describe the vector.  This gives just the number of
              elements, if the vector is a simple vector.  It gives
              the block structure if the vector is a TSFBlockVector
              if the vector is a block vector.  */
          void describe() const;

          /** The companion to describe that indents for readability  */
          void describe(const int& depth) const;





          //@}
          /** timer for math operations */
          static TSFTimer& opTimer();
          /** timer for deep copies */
          static TSFTimer& copyTimer();

    private:

          TSFSmartPtr<TSFVectorBase> ptr_;

    };

  /** \relates TSFVector
   * write an TSFVector to an output stream
   */
  inline ostream& operator<<(ostream& os, const TSFVector& x)
    {
      x.print(os);
      return os;
    }

}


#endif
