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

      /* create the range space  */
      int nLocalRows = 4;

      int rank = MPISession::getRank();
      int nProc = MPISession::getNProc();

      int spaceDimension = nProc*nLocalRows;
      int lowRow = nLocalRows*rank;

      std::vector<int> localRows(nLocalRows);
      for (int i=0; i<nLocalRows; i++)
        {
          localRows[i] = lowRow + i;
        }

      VectorSpace<double> space = type.createSpace(spaceDimension, 
                                                   nLocalRows, 
                                                   &(localRows[0]));

      /* create an empty operator, sized to map domain to range */
      LinearOperator<double> A = type.createMatrix(space, space);
      
      RefCountPtr<LoadableMatrix<double> > mat = A.matrix();

      for (int i=0; i<nLocalRows; i++)
        {
          Array<int> colIndices;
          Array<double> colVals;
          if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows-1))
            {
              colIndices = tuple(localRows[i]);
              colVals = tuple(1.0);
            }
          else
            {
              colIndices = tuple(localRows[i]-1, localRows[i], localRows[i]+1);
              colVals = tuple(-1.0, 2.0, -1.0);
            }
          mat->setRowValues(localRows[i], colIndices.size(), 
                            &(colIndices[0]), &(colVals[0]));
        }

      mat->freezeValues();

      //      cerr << A << endl;

      Vector<double> x = space.createMember();
      Vector<double> y;

      for (int i=0; i<nLocalRows; i++)
        {
          x.setElement(localRows[i], localRows[i]*localRows[i]);
        }

      cerr << "x=" << x << endl;
      A.apply(x, y);
      cerr << "y=" << y << endl;
      

    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();
}

