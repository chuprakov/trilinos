/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_ICCFACTORIZABLEOP_HPP
#define PLAYA_ICCFACTORIZABLEOP_HPP

#include "PlayaDefs.hpp"


namespace Playa
{
template <class Scalar> class Preconditioner;

/** 
 * Base interface for operators for which incomplete Cholesky factorizations
 * can be obtained. 
 */
template <class Scalar>
class ICCFactorizableOp
{
public:
  /** Virtual dtor */
  virtual ~ICCFactorizableOp(){;}


  /** \name incomplete factorization preconditioning interface */
  //@{
  /** create an incomplete factorization. 
   * @param fillLevels number of levels of fill on the local processor
   * @param overlapFill number of levels of fill on remote processors
   * @param relaxationValue fraction of dropped values to be added to the
   * diagonal
   * @param relativeThreshold relative diagonal perutrbation
   * @param absoluteThreshold absolute diagonal perturbation
   * @param rtn newly created preconditioner, returned 
   * by reference argument.
   */
  virtual void getICCPreconditioner(int fillLevels,
    int overlapFill,
    double relaxationValue,
    double relativeThreshold,
    double absoluteThreshold,
    Preconditioner<Scalar>& rtn) const=0;
  //@}
     
      
private:
};
}


#endif
