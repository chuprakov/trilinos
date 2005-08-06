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

// //////////////////////////////////////////////////////////////
// TSFProductVectorSpace.hpp

#ifndef TSFPRODUCTVECTORDECL_HPP
#define TSFPRODUCTVECTORDECL_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "TSFDescribableByTypeID.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFProductVectorSpaceDecl.hpp"
#include "Thyra_ProductVectorBase.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

  /** Implementation of a product vector.
   *
   *  This class allows the construction of a product vector.  This
   *  vector can be built in the usual way by constructing it with a
   *  ProductVectorSpace or it can be built on the fly by using an
   *  empty constructor.  In this case, the ProductVectorSpace is
   *  inferred from the blocks and built at the same time.  The only
   *  error is attempting to redefine a block.
   */

  template<class Scalar>
  class ProductVector : public Thyra::ProductVectorBase<Scalar>, 
                        public DescribableByTypeID,
                        public Handleable<Thyra::VectorBase<Scalar> >
  {
  public:

    GET_RCP(Thyra::VectorBase<Scalar>);
    
    /** Constructor for completed BlockVector, i.e., all the spaced
     * are set and the vector can be built.
     *
     * @param space: a ProductVectorSpace for producing the vector
     */
    ProductVector(const VectorSpace<Scalar> &space);



    /** Sets block k with vec.
     *
     * @param k: int specifying the block
     * @param vec: Vector to be put in block k
     */
    void setBlock(int k, const Vector<Scalar>& vec);

    /** Finalizes the Vector  */
    void finalize();
    

    /** Return the number of blocks  */
    int numBlocks(){return numBlocks_;}

    /** Tests that the space is a ProductVectorSpace 
     *
     * @param space: Space to be tested
     * @param method: string with name of method calling it
     */
    void testSpace(const VectorSpace<Scalar> &space, const string &method);
 

    /** \name Math operations */
    //@{
    /** Multiply this vector by a constant scalar factor 
     * \code
     * this = alpha * this;
     * \endcode
     */
    Vector<Scalar>& scale(const Scalar& alpha);

    //       /** 
    //        * Add a scaled vector to this vector:
    //        * \code
    //        * this = this + alpha*x 
    //        * \endcode
    //        */
    //       Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x);

    //       /** 
    //        * Add a scaled vector to this vector times a constant:
    //        * \code
    //        * this = gamma*this + alpha*x 
    //        * \endcode
    //        */
    //       Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
    //                              const Scalar& gamma);
    //       /** 
    //        * Add two scaled vectors to this vector times a constant:
    //        * \code
    //        * this = alpha*x + beta*y + gamma*this
    //        * \endcode
    //        */
    //       Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
    //                              const Scalar& beta, const Vector<Scalar>& y, 
    //                              const Scalar& gamma);

    /** 
     * Copy the values of another vector into this vector
     * \code
     * this = x
     * \endcode
     */
    Vector<Scalar>& acceptCopyOf(const Vector<Scalar>& x);

    /** 
     * Create a new vector that is a copy of this vector 
     */
    Vector<Scalar> copy() const ;

    /** 
     * Element-by-element product (Matlab dot-star operator)
     */
    Vector<Scalar> dotStar(const Vector<Scalar>& other) const ;

    /** 
     * Element-by-element division (Matlab dot-slash operator)
     */
    Vector<Scalar> dotSlash(const Vector<Scalar>& other) const ;

    /** 
     * Return element-by-element reciprocal as a new vector
     */
    Vector<Scalar> reciprocal() const ;

    /** 
     * Return element-by-element absolute value as a new vector
     */
    Vector<Scalar> abs() const ;

    /** 
     * Overwrite self with element-by-element reciprocal
     */
    Vector<Scalar>& reciprocal() ;

    /** 
     * Overwrite self with element-by-element absolute value 
     */
    Vector<Scalar>& abs() ;

    /** 
     * Set all elements to a constant value
     */
    void setToConstant(const Scalar& alpha) ;

      
    /** 
     * Take dot product with another vector
     */
    Scalar dot(const Vector<Scalar>& other) const ;

    /**
     * Compute the 1-norm of this vector
     */
    Scalar norm1() const ;

    /**
     * Compute the 2-norm of this vector
     */
    Scalar norm2() const ;

    /**
     * Compute the infinity-norm of this vector
     */
    Scalar normInf() const ;

    /**
     * Set all elements to zero 
     */
    void zero();

    //@}

    /** \name Element loading interface */
    //@{

    

    /** Describe the product vector */
    string describe(int depth) const;

  protected:

    /** Empty constructor not to be used.  */
    ProductVector();

    /** bool to denote if all spaces are final  */
    bool isFinal_;


    VectorSpace<Scalar> space_;
    
    //Teuchos::RefCountPtr<VectorSpace<Scalar> > space_;
    Teuchos::Array<Vector<Scalar> > vecsE_;
    Teuchos::Array<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > vecsCore_;
    int numBlocks_;

    /** Make sure that the underlying Thyra::ProductVector is correct. */
    void setCore();
    
  };

} // namespace TSFExtended

#endif // TSF_PRODUCT_VECTOR__HPP
