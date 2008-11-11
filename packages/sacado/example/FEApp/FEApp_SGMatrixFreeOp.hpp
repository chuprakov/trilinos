// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef FEAPP_SGMATRIXFREEOP_HPP
#define FEAPP_SGMATRIXFREEOP_HPP

#include "FEApp_TemplateTypes.hpp"

#if SG_ACTIVE

#include "Teuchos_RCP.hpp"

#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "EpetraExt_BlockMultiVector.h"

namespace FEApp {
    
  /*! 
   * \brief An Epetra operator representing the identity matrix
   */
  class SGMatrixFreeOp : public Epetra_Operator {
      
  public:

    typedef SGType::expansion_type::tp_type tp_type;

    //! Constructor 
    SGMatrixFreeOp(
     const Teuchos::RCP<const Epetra_Map>& base_map,
     const Teuchos::RCP<const Epetra_Map>& sg_map,
     const Teuchos::RCP<const tp_type >& Cijk,
     const Teuchos::RCP<std::vector< Teuchos::RCP<Epetra_CrsMatrix> > >& jacs);
    
    //! Destructor
    virtual ~SGMatrixFreeOp();

    //! Get Jacobian blocks
    virtual std::vector< Teuchos::RCP<Epetra_CrsMatrix> >&
    getJacobianBlocks();
    
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

  private:
    
    //! Private to prohibit copying
    SGMatrixFreeOp(const SGMatrixFreeOp&);
    
    //! Private to prohibit copying
    SGMatrixFreeOp& operator=(const SGMatrixFreeOp&);
    
  protected:
    
    //! Label for operator
    string label;
    
    //! Stores base map
    Teuchos::RCP<const Epetra_Map> base_map;

    //! Stores SG map
    Teuchos::RCP<const Epetra_Map> sg_map;

    //! Stores triple product tensor
    Teuchos::RCP<const tp_type > Cijk;

    //! Stores Jacobian's
    Teuchos::RCP<std::vector< Teuchos::RCP<Epetra_CrsMatrix> > > jacs;

    //! Flag indicating whether transpose was selected
    bool useTranspose;

    //! Number of blocks
    unsigned int num_blocks;

    //! BlockMultiVector for Apply() input
    mutable Teuchos::RCP<EpetraExt::BlockMultiVector> sg_input;

    //! BlockMultiVector for Apply() result
    mutable Teuchos::RCP<EpetraExt::BlockMultiVector> sg_result;

    //! MultiVectors for each block for Apply() input
    mutable std::vector< Teuchos::RCP<Epetra_MultiVector> > input_block;

    //! MultiVectors for each block for Apply() result
    mutable std::vector< Teuchos::RCP<Epetra_MultiVector> > result_block;

    //! Temporary multivector
    mutable Teuchos::RCP<Epetra_MultiVector> tmp;

  }; // class SGMatrixFreeOp
  
} // namespace FEApp

#endif // SG_ACTIVE

#endif // FEAPP_SGMATRIXFREEOP_HPP
