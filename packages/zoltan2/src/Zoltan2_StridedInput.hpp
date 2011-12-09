#ifndef _ZOLTAN2_STRIDEDINPUT_HPP_
#define _ZOLTAN2_STRIDEDINPUT_HPP_

#include <Zoltan2_Standards.hpp>
#include <typeinfo>

/*! \file Zoltan2_StridedInput.hpp
 *  \brief This file defines the StridedInput class.
 */

namespace Zoltan2{

/*! Zoltan2::StridedInput
 *  \brief The StridedInput class manages lists of weights or coordinates.
 *
 * A likely representation for multi-dimensional weights or coordinates is
 *  an array ordered by identifier by dimension. This class is designed
 *  to make access to these arrays efficient. And conversion for TPLs too.
 */

template<typename lno_t, typename scalar_t>
class StridedInput {
private:
  RCP<const Environment> env_;
  ArrayView<const scalar_t> vec_;
  int stride_;

public:

  /*! Constructor
   *    x[0] is the first element of the array.  The subsequent
   *  elements are at x[i*stride].
   */
  StridedInput(RCP<const Environment> env, ArrayView<const scalar_t> x, 
    LNO stride) :  env_(env), vec_(x), stride_(stride) 
  {
  }

  /*! Access an element of the input array.
   */
  scalar_t operator[](lno_t idx) { return vec_[idx*stride_]; }

  /*! Create a contiguous array of the required type, perhaps for a TPL.
   *
   *  TODO: if are there particular conversions we would like to do
   *   for TPLs, we can add methods to do that.  Here we just
   *   essentially cast.  If the cast is not valid (like double to float)
   *   an exception is thrown.
   */

  template <typename T>
    void getInputArray(ArrayRCP<const T> &array)
  {
    size_t n = vec_.size();

    if (n < 1){
      array = ArrayRCP<const T>(Teuchos::Enull); 
    }
    else if (stride_==1 && typeid(T()) == typeid(scalar_t())){
      array = arcpFromArrayView<const T>(vec_);
    }
    else{
      T *tmp = new T [n];
      Z2_LOCAL_MEMORY_ASSERTION(*env_, n, tmp);
      for (lno_t i=0,j=0; i < n; i++,j+=stride_){
        tmp[i] = Teuchos::as<T>(vec_[j]);
      }
      array = arcp(tmp, 0, n, true);
    }
    
    return;
  }

  /*! The raw input information.
   */
  void getStridedList(size_t &len, const scalar_t *&vec, int &stride)
  {
    len = vec_.size();
    vec = vec_.getRawPtr();
    stride = stride_;
  }
};

}  // namespace Zoltan2

#endif
