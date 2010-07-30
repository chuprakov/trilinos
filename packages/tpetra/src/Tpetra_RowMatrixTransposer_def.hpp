//@HEADER
// ************************************************************************
// 
//               Tpetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "Tpetra_Map.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef DOXYGEN_USE_ONLY
  // #include "Tpetra_RowMatrixtransposer_decl.hpp"
#endif

namespace Tpetra {

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::RowMatrixTransposer(const Teuchos::RCP<const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > origMatrix)
  : origMatrix_(origMatrix), comm_(origMatrix->getComm()), indexBase_(origMatrix_->getIndexBase()) {}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::~RowMatrixTransposer() {}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
void RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::createTranspose (const OptimizeOption optimizeTranspose, 
    Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> > &transposeMatrix/*, Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > transposeRowMap*/)
{
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;

  optimizeTranspose_ = optimizeTranspose;
  const size_t LST0 = Teuchos::OrdinalTraits<size_t>::zero();
  const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
  const LocalOrdinal LO0 = Teuchos::OrdinalTraits<LocalOrdinal>::zero();

  // RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > tRowMap;
  // if (transposeRowMap.is_null()) {
  //   tRowMap = origMatrix_->getDomainMap();
  // }
  // else {
  //   tRowMap = transposeRowMap;
  // }

  size_t numMyRows = origMatrix_->getNodeNumRows();
  size_t numMyCols = origMatrix_->getNodeNumCols();
  ArrayRCP<size_t> transNumNz = Teuchos::arcp<size_t>(numMyCols);
  std::fill(transNumNz.begin(), transNumNz.end(), LST0);
  ArrayRCP<ArrayRCP<GlobalOrdinal> > transIndices = Teuchos::arcp<ArrayRCP<LocalOrdinal> >(numMyCols);
  ArrayRCP<ArrayRCP<Scalar>        > transValues  = Teuchos::arcp<ArrayRCP<Scalar>       >(numMyCols);

  // loop over all column indices, counting the number of non-zeros per column
  {
    const RowGraph<LocalOrdinal,GlobalOrdinal,Node> &graph = *origMatrix_->getGraph();
    ArrayRCP<const LocalOrdinal> indices;
    for (LocalOrdinal r=LO0; (size_t)r<numMyRows;++r) {
      indices = graph.getLocalRowView(r);
      for (size_t j=LO0; j<(size_t)indices.size(); ++j) {
        ++transNumNz[indices[j]];
      }
    }
  }

  // now that we know how big the columns are, allocate space.
  for (LocalOrdinal c=LO0; (size_t)c<numMyCols; c++) {
    const size_t numIndices = transNumNz[c];
    if (numIndices>0) {
      transIndices[c] = Teuchos::arcp<GlobalOrdinal>(numIndices);
      transValues[c]  = Teuchos::arcp<Scalar       >(numIndices);
    }
  }


  // pass through the matrix/graph again, transposing the data into transIndices,transValues
  // get local views, add the entries to the transposed arrays. 
  // indices are translated to global indices
  {
    const Map<LocalOrdinal,GlobalOrdinal,Node> &rowMap = *origMatrix_->getRowMap();
    ArrayRCP<const LocalOrdinal> indices;
    ArrayRCP<const Scalar> values;
    // clear these and use them for offsets into the current row
    std::fill( transNumNz.begin(), transNumNz.end(), 0 );
    for (LocalOrdinal lrow=LO0; (size_t)lrow<numMyRows; ++lrow) {
      origMatrix_->getLocalRowView(lrow,indices,values);
      const GlobalOrdinal grow = rowMap.getGlobalElement(lrow);
      for (size_t nnz=LO0; nnz<(size_t)indices.size(); ++nnz) {
        const LocalOrdinal transRow = indices[nnz];
        const size_t loc = transNumNz[transRow];
        transIndices[transRow][loc] = grow;
        transValues [transRow][loc] = values[nnz];
        ++transNumNz[transRow];
      }
    }
  }

  // we already have the matrix data, organized by row. too bad.
  const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > transRowMap = origMatrix_->getColMap();
  transposeMatrix_ = rcp(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(transRowMap, transNumNz, Tpetra::StaticProfile));
  for (size_t c = GST0; c<numMyCols; ++c) {
    transposeMatrix_->insertGlobalValues(transRowMap->getGlobalElement(c), transIndices[c](), transValues[c]());
  }

  // const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > domain_map = origMatrix_->getDomainMap();
  // const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > range_map  = origMatrix_->getRangeMap();
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > origRowMap = origMatrix_->getRowMap();

  RCP<Teuchos::FancyOStream> out_all = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out_all->setOutputToRootOnly(-1);
  RCP<Teuchos::FancyOStream> out_root = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out_root->setOutputToRootOnly(0);

  *out_root << "transpose before fillComplete()" << std::endl;
  transposeMatrix_->describe(*out_all, Teuchos::VERB_EXTREME);

  // for now, we don't care about the domain/range maps, as we won't be applying this matrix as an operator
  transposeMatrix_->fillComplete(transRowMap, origRowMap, optimizeTranspose_);

  transposeMatrix = transposeMatrix_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWMATRIXTRANSPOSE_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class RowMatrixTransposer< SCALAR, LO , GO , NODE >;


}

