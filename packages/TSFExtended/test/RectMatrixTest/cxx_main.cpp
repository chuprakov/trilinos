//@HEADER
//@HEADER

#include "Teuchos_MPISession.hpp"
#include "TSFVector.hpp"
#include "TSFLinearCombination.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"

using namespace Teuchos;
using namespace TSFExtended;
using namespace TSFExtendedOps;


int main(int argc, void *argv[]) 
{
  try
    {
      int verbosity = 1;

      MPISession::init(&argc, &argv);

      VectorType<double> type = new EpetraVectorType();

      /* create the range (row) space  */
      int nLocalRows = 4;

      int rank = MPISession::getRank();
      int nProc = MPISession::getNProc();

      int rangeDimension = nProc*nLocalRows;
      int lowRow = nLocalRows*rank;

      std::vector<int> localRows(nLocalRows);
      for (int i=0; i<nLocalRows; i++)
        {
          localRows[i] = lowRow + i;
        }

      VectorSpace<double> range = type.createSpace(rangeDimension, 
                                                   nLocalRows, 
                                                   &(localRows[0]));

      /* create the domain (column) space */
      int nLocalCols = nLocalRows;
      if (rank==nProc-1) nLocalCols = nLocalRows + 1;

      int domainDimension = rangeDimension + 1;
      int lowCol = nLocalRows*rank;

      std::vector<int> localCols(nLocalCols);
      for (int i=0; i<nLocalCols; i++)
        {
          localCols[i] = lowCol + i;
        }

      VectorSpace<double> domain = type.createSpace(domainDimension, 
                                                    nLocalCols, 
                                                    &(localCols[0]));
      
      /* create an empty operator, sized to map domain to range */
      LinearOperator<double> A = type.createMatrix(domain, range);
      
      RefCountPtr<LoadableMatrix<double> > mat = A.matrix();

      std::vector<int> colIndices(2);
      std::vector<double> colVals(2);
      colVals[0] = 0.5;
      colVals[1] = 0.5;

      for (int i=0; i<nLocalRows; i++)
        {
          colIndices[0] = localRows[i];
          colIndices[1] = localRows[i] + 1;
          mat->setRowValues(localRows[i], 2, 
                            &(colIndices[0]), &(colVals[0]));
        }

      mat->freezeValues();

      cerr << A << endl;

     //  Vector<double> x = domain.createMember();
//       Vector<double> y = range.createMember();

//       for (int i=0; i<nLocalCols; i++)
//         {
//           x.setElement(localCols[i], localCols[i]);
//         }

//       cerr << "x=" << x << endl;
//       A.apply(x, y);
//       cerr << "y=" << y << endl;
      

    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();
}

