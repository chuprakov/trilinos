//@HEADER
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

void timeStats(const vector<double>& opTimings, double& avgTOp,
               double& stdDev)
{
  avgTOp = 0.0;
  for (int i=0; i<opTimings.size(); i++)
    {
      avgTOp += opTimings[i];
    }
  avgTOp = avgTOp/((double) opTimings.size() - 1);

  double varTOp = 0.0;
  for (int i=0; i<opTimings.size(); i++)
    {
      double dt = (opTimings[i] - avgTOp);
      varTOp += dt*dt;
    }
  varTOp = varTOp/((double) opTimings.size() - 1);
  stdDev = ::sqrt(varTOp);
}

int main(int argc, void *argv[]) 
{
  try
    {
      int verbosity = 1;

      MPISession::init(&argc, &argv);

      VectorType<double> type = new EpetraVectorType();

      for (int ln=1; ln<20; ln++)
        {
          int n = 2;
          for (int i=1; i<ln; i++) n *= 2;

          VectorSpace<double> space = type.createSpace(n);

          Vector<double> u = space.createMember();
          Vector<double> w = space.createMember();
          Vector<double> v = space.createMember();
          Vector<double> x = space.createMember();
          Vector<double> y = space.createMember();
          Vector<double> z = space.createMember();

          for (int i=0; i<n; i++)
            {
              u.setElement(i, i);
              v.setElement(i, i*i);
              w.setElement(i, i*i*i);
              x.setElement(i, ::sqrt(i));
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
            }

          /* assign into an empty vector */
          Vector<double> a = space.createMember();

          int nTrials = 10000;

          vector<double> opTimings(nTrials);

          for (int i=0; i<nTrials; i++)
            {
              Time t("op timer");
              t.start();
              a = x + y + z + u + v + w;
              t.stop();
              opTimings[i] = t.totalElapsedTime()/n;
            }
      
          double avgTOp;
          double stdDev;
          timeStats(opTimings, avgTOp, stdDev);
      

          if (verbosity > 1)
            {
              cerr << "a = " << endl << a << endl;
            }

          /* check */
          double err = 0.0;
          for (int i=0; i<n; i++)
            {
              double ai = x.getElement(i) + y.getElement(i) + z.getElement(i)
                +  u.getElement(i) + v.getElement(i) + w.getElement(i);
              err += ::fabs(ai-a.getElement(i));
              if (verbosity > 1)
                {
                  cerr << i << " |ai-a[i]| " << ::fabs(ai-a.getElement(i)) << endl;
                }
            }

          Time t2("TSFCore timer");
          t2.start();

          for (int i=0; i<nTrials; i++)
            {
              Time t("TSFCore timer");
              t.start();
              TSFCore::assign(a.ptr().get(), *(x.ptr()));
              TSFCore::Vp_StV(a.ptr().get(), 1.0, *(y.ptr()));
              TSFCore::Vp_StV(a.ptr().get(), 1.0, *(z.ptr()));
              TSFCore::Vp_StV(a.ptr().get(), 1.0, *(u.ptr()));
              TSFCore::Vp_StV(a.ptr().get(), 1.0, *(v.ptr()));
              TSFCore::Vp_StV(a.ptr().get(), 1.0, *(w.ptr()));
              t.stop();
              opTimings[i] = t.totalElapsedTime()/n;
            }

          double avgTOp2;
          double stdDev2;
          timeStats(opTimings, avgTOp2, stdDev2);
          
          printf("%8d %8d %12.5g %12.5g %12.5g %12.5g\n",
                 n, nTrials, avgTOp, stdDev, avgTOp2, stdDev2);
      
        }
      
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();
}

