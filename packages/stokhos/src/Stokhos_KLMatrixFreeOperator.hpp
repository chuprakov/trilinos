// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Christopher W. Miller(cmiller@math.umd.edu).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_KLMATRIX_FREE_OPERATOR_HPP
#define STOKHOS_KLMATRIX_FREE_OPERATOR_HPP

#include "Stokhos_SGOperator.hpp"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

namespace Stokhos {
    
  /*! 
   * \brief An Epetra operator representing the block stochastic Galerkin
   * operator.
   */
  class KLMatrixFreeOperator : public Stokhos::SGOperator {
      
  public:

    //! Constructor 
    KLMatrixFreeOperator(
     const Teuchos::RCP<const Epetra_Map>& domain_base_map,
     const Teuchos::RCP<const Epetra_Map>& range_base_map,
     const Teuchos::RCP<const Epetra_Map>& domain_sg_map,
     const Teuchos::RCP<const Epetra_Map>& range_sg_map,
     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
    
    //! Destructor
    virtual ~KLMatrixFreeOperator();

    /** \name Stokhos::SGOperator methods */
    //@{

    //! Setup operator
    virtual void setupOperator(
      const Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> >& poly,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk);

    //! Get SG polynomial
    virtual Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > 
    getSGPolynomial();

    //! Get SG polynomial
    virtual Teuchos::RCP<const Stokhos::VectorOrthogPoly<Epetra_Operator> > 
    getSGPolynomial() const;

    //! Get triple product tensor
    virtual Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > 
    getTripleProduct() const;

    //@}

    /** \name Epetra_Operator methods */
    //@{
    
    //! Set to true if the transpose of the operator is requested
    virtual int SetUseTranspose(bool UseTranspose);
    
    /*! 
     * \brief Returns the result of a Epetra_Operator applied to a 
     * Epetra_MultiVector Input in Result as described above.
     */
    virtual int Apply(const Epetra_MultiVector& Input, 
                      Epetra_MultiVector& Result) const;

    /*! 
     * \brief Returns the result of the inverse of the operator applied to a 
     * Epetra_MultiVector Input in Result as described above.
     */
    virtual int ApplyInverse(const Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const;
    
    //! Returns an approximate infinity norm of the operator matrix.
    virtual double NormInf() const;
    
    //! Returns a character string describing the operator
    virtual const char* Label () const;
  
    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;
    
    /*! 
     * \brief Returns true if the \e this object can provide an 
     * approximate Inf-norm, false otherwise.
     */
    virtual bool HasNormInf() const;

    /*! 
     * \brief Returns a reference to the Epetra_Comm communicator 
     * associated with this operator.
     */
    virtual const Epetra_Comm & Comm() const;

    /*!
     * \brief Returns the Epetra_Map object associated with the 
     * domain of this matrix operator.
     */
    virtual const Epetra_Map& OperatorDomainMap () const;

    /*! 
     * \brief Returns the Epetra_Map object associated with the 
     * range of this matrix operator.
     */
    virtual const Epetra_Map& OperatorRangeMap () const;

    //@}

  private:
    
    //! Private to prohibit copying
    KLMatrixFreeOperator(const KLMatrixFreeOperator&);
    
    //! Private to prohibit copying
    KLMatrixFreeOperator& operator=(const KLMatrixFreeOperator&);
    
  protected:
    
    //! Label for operator
    string label;
    
    //! Stores domain base map
    Teuchos::RCP<const Epetra_Map> domain_base_map;

    //! Stores range base map
    Teuchos::RCP<const Epetra_Map> range_base_map;

    //! Stores domain SG map
    Teuchos::RCP<const Epetra_Map> domain_sg_map;

    //! Stores range SG map
    Teuchos::RCP<const Epetra_Map> range_sg_map;

    //! Stochastic Galerking basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > sg_basis;

    //! Short-hand for Cijk
    typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

    //! Stores triple product tensor
    Teuchos::RCP<const Cijk_type> Cijk;

    //! Stores operators
    Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> > block_ops;

    //! Flag indicating whether operator be scaled with <\psi_i^2>
    bool scale_op;

    //! Flag indicating whether to include mean term
    bool include_mean;

    //! Flag indicating whether transpose was selected
    bool useTranspose;

     //! Number of terms in expansion
    int expansion_size;

    //! Number of Jacobian blocks (not necessarily equal to expansion_size)
    int num_blocks;

    //! MultiVectors for each block for Apply() result
    mutable Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > result_block;

    //! Temporary multivector used in Apply()
    mutable Teuchos::RCP<Epetra_MultiVector> block_products;

    //! Temporary multivector used in Apply() for transpose
    mutable Teuchos::RCP<Epetra_MultiVector> block_products_trans;

  }; // class KLMatrixFreeOperator
  
} // namespace Stokhos

#endif // STOKHOS_KLMATRIX_FREE_OPERATOR_HPP
