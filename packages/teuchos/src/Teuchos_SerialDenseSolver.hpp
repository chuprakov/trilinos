//@HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
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

#ifndef _TEUCHOS_SERIALDENSESOLVER_HPP_
#define _TEUCHOS_SERIALDENSESOLVER_HPP_
/*! \file Teuchos_SerialDenseSolver.hpp
  \brief Templated class for solving dense linear problems.
*/

#include "Teuchos_CompObject.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

/*! \class Teuchos::SerialDenseSolver
  \brief A class for solving dense linear problems.

  The Teuchos::SerialDenseSolver class enables the definition, in terms of Teuchos::SerialDenseMatrix 
  and Teuchos::SerialDenseVector objects, of a dense linear problem, followed by the solution of that 
  problem via the most sophisticated techniques available in LAPACK.

  The Teuchos::SerialDenseSolver class is intended to provide full-featured support for solving linear 
  problems for general dense rectangular (or square) matrices.  It is written on top of BLAS and LAPACK 
  and thus has excellent performance and numerical capabilities.  Using this class, one can either perform 
  simple factorizations and solves or apply all the tricks available in LAPACK to get the best possible 
  solution for very ill-conditioned problems.
  
  <b>Teuchos::SerialDenseSolver vs. Teuchos::LAPACK</b>
  
  The Teuchos::LAPACK class provides access to most of the same functionality as Teuchos::SerialDenseSolver.
  The primary difference is that Teuchos::LAPACK is a "thin" layer on top of LAPACK and Teuchos::SerialDenseSolver
  attempts to provide easy access to the more sophisticated aspects of solving dense linear and eigensystems.
  <ul>
  <li> When you should use Teuchos::LAPACK:  If you are simply looking for a convenient wrapper around 
  the Fortran LAPACK routines and you have a well-conditioned problem, you should probably use Teuchos::LAPACK directly.
  <li> When you should use Teuchos::SerialDenseSolver: If you want to (or potentially want to) solve 
  ill-conditioned problems or want to work with a more object-oriented interface, you should probably use 
  Teuchos::SerialDenseSolver.
  </ul>

  <b>Constructing Teuchos::SerialDenseSolver Objects</b>
  
  There is a single Teuchos::SerialDenseSolver constructor.   However, the matrix, right hand side and solution 
  vectors must be set prior to executing most methods in this class.
  
  <b>Setting vectors used for linear solves</b>
  
  The matrix A, the left hand side X and the right hand side B (when solving AX = B, for X), can be set by appropriate set
  methods.  Each of these three objects must be an Teuchos::SerialDenseMatrix or and Teuchos::SerialDenseVector object.  The
  set methods are as follows:
  <ul>
  <li> setMatrix()  - Sets the matrix.
  <li> setVectors() - Sets the left and right hand side vector(s).
  </ul>
  
  <b>Vector and Utility Functions</b>
  
  Once a Teuchos::SerialDenseSolver is constructed, several mathematical functions can be applied to
  the object.  Specifically:
  <ul>
  <li> Factorizations.
  <li> Solves.
  <li> Condition estimates.
  <li> Equilibration.
  <li> Norms.
  </ul>
  
  <b>Strategies for Solving Linear Systems</b>
  In many cases, linear systems can be accurately solved by simply computing the LU factorization
  of the matrix and then performing a forward back solve with a given set of right hand side vectors.  However,
  in some instances, the factorization may be very poorly conditioned and this simple approach may not work.  In
  these situations, equilibration and iterative refinement may improve the accuracy, or prevent a breakdown in
  the factorization. 
  
  Teuchos::SerialDenseSolver will use equilibration with the factorization if, once the object
  is constructed and \e before it is factored, you call the function factorWithEquilibration(true) to force 
  equilibration to be used.  If you are uncertain if equilibration should be used, you may call the function
  shouldEquilibrate() which will return true if equilibration could possibly help.  shouldEquilibrate() uses
  guidelines specified in the LAPACK User Guide, namely if SCOND < 0.1 and AMAX < Underflow or AMAX > Overflow, to 
  determine if equilibration \e might be useful. 
  
  Teuchos::SerialDenseSolver will use iterative refinement after a forward/back solve if you call
  solveToRefinedSolution(true).  It will also compute forward and backward error estimates if you call
  estimateSolutionErrors(true).  Access to the forward (back) error estimates is available via FERR() (BERR()).
  
  Examples using Teuchos::SerialDenseSolver can be found in the Teuchos test directories.
*/

namespace Teuchos {

  template<typename OrdinalType, typename ScalarType>
  class SerialDenseSolver : public CompObject, public Object, public BLAS<OrdinalType, ScalarType>,
			    public LAPACK<OrdinalType, ScalarType>
  {

  public:

    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
  
    //! @name Constructor/Destructor Methods
    //@{ 
    //! Default constructor; matrix should be set using setMatrix(), LHS and RHS set with setVectors().
    SerialDenseSolver();
    
    //! SerialDenseSolver destructor.  
    virtual ~SerialDenseSolver();
    //@}
    
    //! @name Set Methods
    //@{ 
    
    //! Sets the pointers for coefficient matrix
    int setMatrix(const RCP<SerialDenseMatrix<OrdinalType, ScalarType> >& A);

    //! Sets the pointers for left and right hand side vector(s).
    /*! Row dimension of X must match column dimension of matrix A, row dimension of B 
      must match row dimension of A.  X and B must have the same dimensions.
    */
    int setVectors(const RCP<SerialDenseMatrix<OrdinalType, ScalarType> >& X, 
		   const RCP<SerialDenseMatrix<OrdinalType, ScalarType> >& B);
    //@}
    
    //! @name Strategy Modifying Methods
    //@{ 
    
    //! Causes equilibration to be called just before the matrix factorization as part of the call to \c factor.
    /*! This function must be called before the factorization is performed. 
   */
    void factorWithEquilibration(bool flag) {equilibrate_ = flag; return;}
    
    //! If \c flag is true, causes all subsequent function calls to work with the transpose of \e this matrix, otherwise not.
    void solveWithTranspose(bool flag) {transpose_ = flag; if (flag) TRANS_ = Teuchos::TRANS; else TRANS_ = Teuchos::NO_TRANS; return;}
    
    //! Causes all solves to compute solution to best ability using iterative refinement.
    void solveToRefinedSolution(bool flag) {refineSolution_ = flag; return;}
    
    //! Causes all solves to estimate the forward and backward solution error. 
    /*! Error estimates will be in the arrays FERR and BERR, resp, after the solve step is complete.
      These arrays are accessible via the FERR() and BERR() access functions.
    */
    void estimateSolutionErrors(bool flag);
    //@}
    
    //! @name Factor/Solve/Invert Methods
    //@{ 
    
    //! Computes the in-place LU factorization of the matrix using the LAPACK routine \e _GETRF.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int factor();
    
    //! Computes the solution X to AX = B for the \e this matrix and the B provided to SetVectors()..
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int solve();
    
    //! Inverts the \e this matrix.
    /*!
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int invert();
    
    //! Computes the scaling vector S(i) = 1/sqrt(A(i,i)) of the \e this matrix.
    /*! 
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int computeEquilibrateScaling();
    
    //! Equilibrates the \e this matrix.
    /*! 
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int equilibrateMatrix();
    
    //! Equilibrates the current RHS.
    /*! 
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int equilibrateRHS();
        
    //! Apply Iterative Refinement.
    /*! 
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int applyRefinement();
    
    //! Unscales the solution vectors if equilibration was used to solve the system.
    /*! 
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int unequilibrateLHS();
    
    //! Returns the reciprocal of the 1-norm condition number of the \e this matrix.
    /*! 
      \param Value Out
      On return contains the reciprocal of the 1-norm condition number of the \e this matrix.
      
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
     int reciprocalConditionEstimate(MagnitudeType & Value);
    //@}
    
    //! @name Query methods
    //@{ 
    
    //! Returns true if transpose of \e this matrix has and will be used.
    bool transpose() {return(transpose_);}
    
    //! Returns true if matrix is factored (factor available via AF() and LDAF()).
    bool factored() {return(factored_);}
    
    //! Returns true if factor is equilibrated (factor available via AF() and LDAF()).
    bool equilibratedA() {return(equilibratedA_);}
    
    //! Returns true if RHS is equilibrated (RHS available via B() and LDB()).
    bool equilibratedB() {return(equilibratedB_);}
    
    //! Returns true if the LAPACK general rules for equilibration suggest you should equilibrate the system.
     bool shouldEquilibrate() {computeEquilibrateScaling(); return(shouldEquilibrate_);}
    
    //! Returns true if forward and backward error estimated have been computed (available via FERR() and BERR()).
    bool solutionErrorsEstimated() {return(solutionErrorsEstimated_);}
    
    //! Returns true if matrix inverse has been computed (inverse available via AF() and LDAF()).
    bool inverted() {return(inverted_);}
    
    //! Returns true if the condition number of the \e this matrix has been computed (value available via ReciprocalConditionEstimate()).
    bool reciprocalConditionEstimated() {return(reciprocalConditionEstimated_);}
    
    //! Returns true if the current set of vectors has been solved.
    bool solved() {return(solved_);}
    
    //! Returns true if the current set of vectors has been refined.
    bool solutionRefined() {return(solutionRefined_);}
    //@}
    
    //! @name Data Accessor methods
    //@{ 
    
    //! Returns pointer to current matrix.
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getMatrix()  const {return(Matrix_);}
    
    //! Returns pointer to factored matrix (assuming factorization has been performed).
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getFactoredMatrix()  const {return(Factor_);}
    
    //! Returns pointer to current LHS.
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getLHS() const {return(LHS_);}
    
    //! Returns pointer to current RHS.
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getRHS() const {return(RHS_);}
    
    //! Returns row dimension of system.
    OrdinalType numRows()  const {return(M_);}
    
    //! Returns column dimension of system.
    OrdinalType numCols()  const {return(N_);}
    
    //! Returns pointer to pivot vector (if factorization has been computed), zero otherwise.
    std::vector<OrdinalType> IPIV()  const {return(IPIV_);}
    
    //! Returns the 1-Norm of the \e this matrix (returns -1 if not yet computed).
    MagnitudeType ANORM()  const {return(ANORM_);}
    
    //! Returns the reciprocal of the condition number of the \e this matrix (returns -1 if not yet computed).
    MagnitudeType RCOND()  const {return(RCOND_);}
    
    //! Ratio of smallest to largest row scale factors for the \e this matrix (returns -1 if not yet computed).
    /*! If ROWCND() is >= 0.1 and AMAX() is not close to overflow or underflow, then equilibration is not needed.
     */
    MagnitudeType ROWCND()  const {return(ROWCND_);}
    
    //! Ratio of smallest to largest column scale factors for the \e this matrix (returns -1 if not yet computed).
    /*! If COLCND() is >= 0.1 then equilibration is not needed.
     */
    MagnitudeType COLCND()  const {return(COLCND_);}
    
    //! Returns the absolute value of the largest entry of the \e this matrix (returns -1 if not yet computed).
    MagnitudeType AMAX()  const {return(AMAX_);}
    
    //! Returns a pointer to the forward error estimates computed by LAPACK.
    std::vector<MagnitudeType> FERR()  const {return(FERR_);}
    
    //! Returns a pointer to the backward error estimates computed by LAPACK.
    std::vector<MagnitudeType> BERR()  const {return(BERR_);}
    
    //! Returns a pointer to the row scaling vector used for equilibration.
    std::vector<MagnitudeType> R()  const {return(R_);}
    
    //! Returns a pointer to the column scale vector used for equilibration.
    std::vector<MagnitudeType> C()  const {return(C_);}
    //@}
    
    //! @name I/O methods
    //@{ 
    //! Print service methods; defines behavior of ostream << operator.
    void Print(std::ostream& os) const;
    //@}
  protected:
    
    void allocateWORK() { LWORK_ = 4*N_; WORK_.resize( LWORK_ ); return;}
    void allocateIWORK() { IWORK_.resize( N_ ); return;}
    void resetMatrix();
    void resetVectors();
    
    
    bool equilibrate_;
    bool shouldEquilibrate_;
    bool equilibratedA_;
    bool equilibratedB_;
    bool transpose_;
    bool factored_;
    bool estimateSolutionErrors_;
    bool solutionErrorsEstimated_;
    bool solved_;
    bool inverted_;
    bool reciprocalConditionEstimated_;
    bool refineSolution_;
    bool solutionRefined_;
    
    Teuchos::ETransp TRANS_;
    
    OrdinalType M_;
    OrdinalType N_;
    OrdinalType Min_MN_;
    OrdinalType LDA_;
    OrdinalType LDAF_;
    OrdinalType INFO_;
    OrdinalType LWORK_;
    
    std::vector<OrdinalType> IPIV_;
    std::vector<int> IWORK_;
    
    MagnitudeType ANORM_;
    MagnitudeType RCOND_;
    MagnitudeType ROWCND_;
    MagnitudeType COLCND_;
    MagnitudeType AMAX_;
    
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > Matrix_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > LHS_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > RHS_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > Factor_;
    
    ScalarType * A_;
    ScalarType * AF_;
    std::vector<MagnitudeType> FERR_;
    std::vector<MagnitudeType> BERR_;
    std::vector<ScalarType> WORK_;
    std::vector<MagnitudeType> R_;
    std::vector<MagnitudeType> C_;
    
  private:
    // SerialDenseSolver copy constructor (put here because we don't want user access)
    
    SerialDenseSolver(const SerialDenseSolver<OrdinalType, ScalarType>& Source);
    SerialDenseSolver & operator=(const SerialDenseSolver<OrdinalType, ScalarType>& Source);
    
  };  

//=============================================================================

template<typename OrdinalType, typename ScalarType>
SerialDenseSolver<OrdinalType,ScalarType>::SerialDenseSolver()
  : CompObject(),
    Object("Teuchos::SerialDenseSolver"),
    equilibrate_(false),
    shouldEquilibrate_(false),
    equilibratedA_(false),
    equilibratedB_(false),
    transpose_(false),
    factored_(false),
    estimateSolutionErrors_(false),
    solutionErrorsEstimated_(false),
    solved_(false),
    inverted_(false),
    reciprocalConditionEstimated_(false),
    refineSolution_(false),
    solutionRefined_(false),
    TRANS_(Teuchos::NO_TRANS),
    M_(0),
    N_(0),
    Min_MN_(0),
    LDA_(0),
    LDAF_(0),
    INFO_(0),
    LWORK_(0),
    ANORM_(ScalarTraits<ScalarType>::zero()),
    RCOND_(ScalarTraits<ScalarType>::zero()),
    ROWCND_(ScalarTraits<ScalarType>::zero()),
    COLCND_(ScalarTraits<ScalarType>::zero()),
    AMAX_(ScalarTraits<ScalarType>::zero()),
    A_(0),
    AF_(0)
{
  resetMatrix();
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
SerialDenseSolver<OrdinalType,ScalarType>::~SerialDenseSolver()
{}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialDenseSolver<OrdinalType,ScalarType>::resetVectors()
{
  LHS_ = Teuchos::null;
  RHS_ = Teuchos::null;
  reciprocalConditionEstimated_ = false;
  solutionRefined_ = false;
  solved_ = false;
  solutionErrorsEstimated_ = false;
  equilibratedB_ = false;
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialDenseSolver<OrdinalType,ScalarType>::resetMatrix()
{
  resetVectors();
  equilibratedA_ = false;
  factored_ = false;
  inverted_ = false;
  M_ = 0;
  N_ = 0;
  Min_MN_ = 0;
  LDA_ = 0;
  LDAF_ = 0;
  ANORM_ = -ScalarTraits<MagnitudeType>::one();
  RCOND_ = -ScalarTraits<MagnitudeType>::one();
  ROWCND_ = -ScalarTraits<MagnitudeType>::one();
  COLCND_ = -ScalarTraits<MagnitudeType>::one();
  AMAX_ = -ScalarTraits<MagnitudeType>::one();
  A_ = 0;
  AF_ = 0;
  INFO_ = 0;
  LWORK_ = 0;
  R_.resize(0);
  C_.resize(0);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialDenseSolver<OrdinalType,ScalarType>::setMatrix(const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& A) 
{
  resetMatrix();
  Matrix_ = A;
  Factor_ = A;
  M_ = A->numRows();
  N_ = A->numCols();
  Min_MN_ = TEUCHOS_MIN(M_,N_);
  LDA_ = A->stride();
  LDAF_ = LDA_;
  A_ = A->values();
  AF_ = A->values();
  return(0);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialDenseSolver<OrdinalType,ScalarType>::setVectors(const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& X, 
							   const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& B)
{
  // Check that these new vectors are consistent.
  TEST_FOR_EXCEPTION(B->numRows()!=X->numRows() || B->numCols() != X->numCols(), std::invalid_argument,
		     "SerialDenseSolver<T>::setVectors: X and B are not the same size!");
  TEST_FOR_EXCEPTION(B->values()==0, std::invalid_argument,
		     "SerialDenseSolver<T>::setVectors: B is an empty SerialDenseMatrix<T>!");
  TEST_FOR_EXCEPTION(X->values()==0, std::invalid_argument,
		     "SerialDenseSolver<T>::setVectors: X is an empty SerialDenseMatrix<T>!");
  TEST_FOR_EXCEPTION(B->stride()<1, std::invalid_argument,
		     "SerialDenseSolver<T>::setVectors: B has an invalid stride!");
  TEST_FOR_EXCEPTION(X->stride()<1, std::invalid_argument,
		     "SerialDenseSolver<T>::setVectors: X has an invalid stride!");

  resetVectors(); 
  LHS_ = X;
  RHS_ = B;
  return(0);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialDenseSolver<OrdinalType,ScalarType>::estimateSolutionErrors(bool flag) 
{
  estimateSolutionErrors_ = flag;

  // If the errors are estimated, this implies that the solution must be refined
  refineSolution_ = refineSolution_ || flag;
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialDenseSolver<OrdinalType,ScalarType>::factor() {

  if (factored()) return(0); // Already factored

  TEST_FOR_EXCEPTION(inverted(), std::logic_error,
		     "SerialDenseSolver<T>::factor: Cannot factor an inverted matrix!");

  ANORM_ = Matrix_->normOne(); // Compute 1-Norm of A


  // If we want to refine the solution, then the factor must
  // be stored separatedly from the original matrix

  if (A_ == AF_)
    if (refineSolution_ ) {
      Factor_ = rcp( new SerialDenseMatrix<OrdinalType,ScalarType>(*Matrix_) );
      AF_ = Factor_->values();
      LDAF_ = Factor_->stride();
    }
  
  int ierr = 0;
  if (equilibrate_) ierr = equilibrateMatrix();

  if (ierr!=0) return(ierr);
  
  if (IPIV_.size()==0) IPIV_.resize( Min_MN_ ); // Allocated Pivot vector if not already done.

  INFO_ = 0;
  this->GETRF(M_, N_, AF_, LDAF_, &IPIV_[0], &INFO_);
  factored_ = true;

  return(INFO_);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialDenseSolver<OrdinalType,ScalarType>::solve() {

  // We will call one of four routines depending on what services the user wants and 
  // whether or not the matrix has been inverted or factored already.
  //
  // If the matrix has been inverted, use DGEMM to compute solution.
  // Otherwise, if the user want the matrix to be equilibrated or wants a refined solution, we will
  // call the X interface.
  // Otherwise, if the matrix is already factored we will call the TRS interface.
  // Otherwise, if the matrix is unfactored we will call the SV interface.

  int ierr = 0;
  if (equilibrate_) {
    ierr = equilibrateRHS();
    equilibratedB_ = true;
  }
  if (ierr != 0) return(ierr);  // Can't equilibrate B, so return.

  TEST_FOR_EXCEPTION( (equilibratedA_ && !equilibratedB_) || (!equilibratedA_ && equilibratedB_) , 
                     std::logic_error, "SerialDenseSolver<T>::solve: Matrix and vectors must be similarly scaled!");
  TEST_FOR_EXCEPTION( RHS_==Teuchos::null, std::invalid_argument, 
                     "SerialDenseSolver<T>::solve: No right-hand side vector (RHS) has been set for the linear system!");
  TEST_FOR_EXCEPTION( LHS_==Teuchos::null, std::invalid_argument, 
                     "SerialDenseSolver<T>::solve: No solution vector (LHS) has been set for the linear system!");

  if (shouldEquilibrate() && !equilibratedA_)
    std::cout << "WARNING!  SerialDenseSolver<T>::solve: System should be equilibrated!" << std::endl;

  if (inverted()) {

    TEST_FOR_EXCEPTION( RHS_->values() == LHS_->values(), std::invalid_argument, 
                        "SerialDenseSolver<T>::solve: X and B must be different vectors if matrix is inverted.");

    INFO_ = 0;
    this->GEMM(TRANS_, Teuchos::NO_TRANS, N_, RHS_->numCols(), N_, 1.0, AF_, LDAF_,
	       RHS_->values(), RHS_->stride(), 0.0, LHS_->values(), LHS_->stride());
    if (INFO_!=0) return(INFO_);
    solved_ = true;
  }
  else {

    if (!factored()) factor(); // Matrix must be factored
    
    if (RHS_->values()!=LHS_->values()) {
       (*LHS_) = (*RHS_); // Copy B to X if needed
    }
    INFO_ = 0;
    this->GETRS(ETranspChar[TRANS_], N_, RHS_->numCols(), AF_, LDAF_, &IPIV_[0], LHS_->values(), LHS_->stride(), &INFO_);
    if (INFO_!=0) return(INFO_);
    solved_ = true;

  }
  int ierr1=0;
  if (refineSolution_ && !inverted()) ierr1 = applyRefinement();
  if (ierr1!=0) 
    return(ierr1);
  else
    return(ierr);
  
  if (equilibrate_) ierr1 = unequilibrateLHS();
  return(ierr1);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialDenseSolver<OrdinalType,ScalarType>::applyRefinement()
{
  TEST_FOR_EXCEPTION(!solved(), std::logic_error,
		     "SerialDenseSolver<T>::applyRefinement: Must have an existing solution!");
  TEST_FOR_EXCEPTION(A_==AF_, std::logic_error,
		     "SerialDenseSolver<T>::applyRefinement: Cannot apply refinement if no original copy of A!");

  OrdinalType NRHS = RHS_->numCols();
  FERR_.resize( NRHS );
  BERR_.resize( NRHS );
  allocateWORK();
  allocateIWORK();
  
  INFO_ = 0;
  this->GERFS(ETranspChar[TRANS_], N_, NRHS, A_, LDA_, AF_, LDAF_, &IPIV_[0],
	      RHS_->values(), RHS_->stride(), LHS_->values(), LHS_->stride(), 
              &FERR_[0], &BERR_[0], &WORK_[0], &IWORK_[0], &INFO_);
  
  solutionErrorsEstimated_ = true;
  reciprocalConditionEstimated_ = true;
  solutionRefined_ = true;
  
  return(INFO_);

}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialDenseSolver<OrdinalType,ScalarType>::computeEquilibrateScaling() 
{
  if (R_.size()!=0) return(0); // Already computed
 
  R_.resize( M_ );
  C_.resize( N_ );
  
  INFO_ = 0;
  this->GEEQU (M_, N_, AF_, LDAF_, &R_[0], &C_[0], &ROWCND_, &COLCND_, &AMAX_, &INFO_);

  if (COLCND_<0.1*ScalarTraits<MagnitudeType>::one() || 
      ROWCND_<0.1*ScalarTraits<MagnitudeType>::one() || 
      AMAX_ < ScalarTraits<ScalarType>::rmin() || AMAX_ > ScalarTraits<ScalarType>::rmax() ) 
    shouldEquilibrate_ = true;

  return(INFO_);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialDenseSolver<OrdinalType,ScalarType>::equilibrateMatrix()
{
  OrdinalType i, j;
  int ierr = 0;

  if (equilibratedA_) return(0); // Already done.
  if (R_.size()==0) ierr = computeEquilibrateScaling(); // Compute R and C if needed.
  if (ierr!=0) return(ierr);     // If return value is not zero, then can't equilibrate.
  if (A_==AF_) {
    ScalarType * ptr;
    for (j=0; j<N_; j++) {
      ptr = A_ + j*LDA_;
      ScalarType s1 = C_[j];
      for (i=0; i<M_; i++) {
	*ptr = *ptr*s1*R_[i];
	ptr++;
      }
    }
  }
  else {
    ScalarType * ptr;
    ScalarType * ptr1;
    for (j=0; j<N_; j++) {
      ptr = A_ + j*LDA_;
      ptr1 = AF_ + j*LDAF_;
      ScalarType s1 = C_[j];
      for (i=0; i<M_; i++) {
	*ptr = *ptr*s1*R_[i];
	ptr++;
	*ptr1 = *ptr1*s1*R_[i];
	ptr1++;
      }
    }
  }
  
  equilibratedA_ = true;
  
  return(0);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialDenseSolver<OrdinalType,ScalarType>::equilibrateRHS()
{
  OrdinalType i, j;
  int ierr = 0;

  if (equilibratedB_) return(0); // Already done.
  if (R_.size()==0) ierr = computeEquilibrateScaling(); // Compute R and C if needed.
  if (ierr!=0) return(ierr);     // Can't count on R and C being computed.

  ScalarType * R_tmp = &R_[0];
  if (transpose_) R_tmp = &C_[0];

  OrdinalType LDB = RHS_->stride(), NRHS = RHS_->numCols();
  ScalarType * B = RHS_->values();
  ScalarType * ptr;
  for (j=0; j<NRHS; j++) {
    ptr = B + j*LDB;
    for (i=0; i<M_; i++) {
      *ptr = *ptr*R_tmp[i];
      ptr++;
    }
  }
  
  equilibratedB_ = true;

  return(0);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialDenseSolver<OrdinalType,ScalarType>::unequilibrateLHS()
{
  OrdinalType i, j;

  if (!equilibratedB_) return(0); // Nothing to do

  ScalarType * C_tmp = &C_[0];
  if (transpose_) C_tmp = &R_[0];

  OrdinalType LDX = RHS_->stride(), NRHS = RHS_->numCols();
  ScalarType * X = RHS_->values();
  ScalarType * ptr;
  for (j=0; j<NRHS; j++) {
    ptr = X + j*LDX;
    for (i=0; i<N_; i++) {
      *ptr = *ptr*C_tmp[i];
      ptr++;
    }
  }

  return(0);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialDenseSolver<OrdinalType,ScalarType>::invert()
{
  if (!factored()) factor(); // Need matrix factored.

  /* This section work with LAPACK Version 3.0 only 
  // Setting LWORK = -1 and calling GETRI will return optimal work space size in WORK_TMP
  OrdinalType LWORK_TMP = -1;
  ScalarType WORK_TMP;
  GETRI( N_, AF_, LDAF_, &IPIV_[0], &WORK_TMP, &LWORK_TMP, &INFO_);
  LWORK_TMP = WORK_TMP; // Convert to integer
  if (LWORK_TMP>LWORK_) {
  LWORK_ = LWORK_TMP;
  WORK_.resize( LWORK_ );
  }
  */
  // This section will work with any version of LAPACK 
  allocateWORK();

  INFO_ = 0;
  this->GETRI( N_, AF_, LDAF_, &IPIV_[0], &WORK_[0], LWORK_, &INFO_);

  inverted_ = true;
  factored_ = false;
  
  return(INFO_);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialDenseSolver<OrdinalType,ScalarType>::reciprocalConditionEstimate(MagnitudeType & Value)
{
  if (reciprocalConditionEstimated()) {
    Value = RCOND_;
    return(0); // Already computed, just return it.
  }

  if ( ANORM_<ScalarTraits<MagnitudeType>::zero() ) ANORM_ = Matrix_->normOne();

  int ierr = 0;
  if (!factored()) ierr = factor(); // Need matrix factored.
  if (ierr!=0) return(ierr);

  allocateWORK();
  allocateIWORK();

  // We will assume a one-norm condition number
  INFO_ = 0;
  this->GECON( '1', N_, AF_, LDAF_, ANORM_, &RCOND_, &WORK_[0], &IWORK_[0], &INFO_);
  reciprocalConditionEstimated_ = true;
  Value = RCOND_;

  return(INFO_);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialDenseSolver<OrdinalType,ScalarType>::Print(std::ostream& os) const {

  if (Matrix_!=Teuchos::null) os << "Solver Matrix"          << std::endl << *Matrix_ << std::endl;
  if (Factor_!=Teuchos::null) os << "Solver Factored Matrix" << std::endl << *Factor_ << std::endl;
  if (LHS_   !=Teuchos::null) os << "Solver LHS"             << std::endl << *LHS_    << std::endl;
  if (RHS_   !=Teuchos::null) os << "Solver RHS"             << std::endl << *RHS_    << std::endl;

}

} // namespace Teuchos

#endif /* _TEUCHOS_SERIALDENSESOLVER_HPP_ */
