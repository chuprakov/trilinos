//@HEADER
//@HEADER

#include "Teuchos_MPISession.hpp"
#include "TSFVector.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFEpetraVectorType.hpp"

using namespace Teuchos;
using namespace TSFExtended;

int main(int argc, void *argv[]) 
{
  try
    {
      MPISession::init(&argc, &argv);

      VectorType<double> type = new EpetraVectorType();

      int n = 10;

      VectorSpace<double> space = type.createSpace(n);

      Vector<double> vec = space.createMember();

      for (int i=0; i<n; i++)
        {
          vec.setElement(i, ::sqrt(i));
        }

      cerr << "vector = " << endl << vec << endl;
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();
}

