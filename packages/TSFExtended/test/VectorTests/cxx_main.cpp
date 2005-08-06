//@HEADER
// ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
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
//@HEADER

#include "Teuchos_MPISession.hpp"
#include "TSFVector.hpp"
#include "TSFLinearCombination.hpp"
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
      int verbosity = 2;

      MPISession::init(&argc, &argv);

      VectorType<double> type = new EpetraVectorType();

      int n = 4;

      int rank = MPISession::getRank();
      int nProc = MPISession::getNProc();

      int dimension = nProc*n;
      int low = n*rank;
      int high = n*(rank+1);
      std::vector<int> localRows(n);
      for (int i=0; i<n; i++)
        {
          localRows[i] = low + i;
        }

      VectorSpace<double> space = type.createSpace(dimension, n, 
                                                   &(localRows[0]));


      cerr << "vector space: " << space.description() << endl;

      Vector<double> u = space.createMember();
      Vector<double> w = space.createMember();
      Vector<double> v = space.createMember();
      Vector<double> x = space.createMember();
      Vector<double> y = space.createMember();
      Vector<double> z = space.createMember();

      for (int i=low; i<high; i++)
        {
          u.setElement(i, 1+i);
          v.setElement(i, i*i);
          w.setElement(i, i*i*i);
          x.setElement(i, 1 + ::sqrt(i));
          y.setElement(i, ::cos(i));
          z.setElement(i, ::sin(i));
        }

      if (verbosity > 1)
        {
          cerr << "u = " << endl << u << endl;
          cerr << "v = " << endl << v << endl;
          cerr << "w = " << endl << w << endl;
          
          cerr << "x = " << endl << x << endl;
          cerr << "y = " << endl << y << endl;
          cerr << "z = " << endl << z << endl;

          cerr << "|u|= " << u.norm2() << endl;
        }

      /* assign into an empty vector */
      Vector<double> a = space.createMember();
      Vector<double> aCheck = space.createMember();
      Vector<double> b = space.createMember();
      Vector<double> bCheck = space.createMember();

      Vector<double> c = space.createMember();
      Vector<double> cCheck = space.createMember();

      Vector<double> d = space.createMember();
      Vector<double> dCheck = space.createMember();

      cerr << "adding... " ;

      a = x + y + z + u + v + w;
      
      cerr << "done" << endl << "doing dotStar... ";

      b = x.dotStar(y);

      cerr << "done" << endl << "doing dotSlash... ";

      c = b.dotSlash(x);

      cerr << "done" << endl << "doing scale..." ;
      d = x * 2.0;
      
      if (verbosity > 1)
        {
          cerr << "a = " << endl << a << endl;
          cerr << "b = " << endl << b << endl;
          cerr << "c = " << endl << c << endl;
          cerr << "d = " << endl << d << endl;
        }


      /* check */

      cerr << "checking..." << endl;
      double aErr = 0.0;
      double bErr = 0.0;

      for (int i=low; i<high; i++)
        {
          double ai = x.getElement(i) + y.getElement(i) + z.getElement(i)
            +  u.getElement(i) + v.getElement(i) + w.getElement(i);
          aCheck.setElement(i, ai);
          bCheck.setElement(i, x.getElement(i) * y.getElement(i) );
          cCheck.setElement(i, b.getElement(i) / x.getElement(i) );
          aErr += ::fabs(ai-a.getElement(i));
          bErr += ::fabs(bCheck.getElement(i)-b.getElement(i));
          if (verbosity > 1)
            {
              cerr << i << " |ai-a[i]| " << ::fabs(ai-a.getElement(i)) << endl;
              cerr << i << " |bi-b[i]| " << ::fabs(bCheck.getElement(i)-b.getElement(i)) << endl;
              cerr << i << " |ci-c[i]| " << ::fabs(cCheck.getElement(i)-c.getElement(i)) << endl;
            }
        }

      
      if (verbosity > 0) cerr << "error norm a  (computed term-by-term) = " 
                              << aErr << endl;

      if (verbosity > 0) cerr << "error norm a (computed with operators) = " 
                              << (a-aCheck).norm2() << endl;

      if (verbosity > 0) cerr << "error norm b  (computed term-by-term) = " 
                              << bErr << endl;

      if (verbosity > 0) cerr << "error norm b (computed with operators) = " 
                              << (b-bCheck).norm2() << endl;

      if (verbosity > 0) cerr << "error norm c (computed with operators) = " 
                              << (c-y).norm2() << endl;

      
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();
}

