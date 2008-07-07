#include <iostream>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_TestForException.hpp"

// Inefficient vector/tensor objects
#include "Special_AlgebraicTypes.hpp"

// TVMET - efficient expression template vector/tensor objects
#include "tvmet/Vector.h"
#include "tvmet/Matrix.h"

int main(int argc, char* argv[])
{

  using namespace std;
  using namespace Teuchos;

  const int num_samples = 1;
  const int num_loops = 500000;
  const int size = 10;
  const int num_vectors = 3;

  TEST_FOR_EXCEPTION(num_loops * size != 5000000, std::logic_error,
		     "Work amount is not constant!");

  // Make all vectors in a contiguous block

  MyVector<double>* vector_array = new MyVector<double>[num_vectors * size];
  ArrayRCP< MyVector<double> > a = 
    arcp< MyVector<double> >(vector_array, 0, size, false);
  ArrayRCP< MyVector<double> > b = 
    arcp< MyVector<double> >(&vector_array[size], 0, size, false);
  ArrayRCP< MyVector<double> > c = 
    arcp< MyVector<double> >(&vector_array[2*size], 0, size, false);

  tvmet::Vector<double, 3>* tvmet_array = 
    new tvmet::Vector<double, 3>[num_vectors * size];
  ArrayRCP< tvmet::Vector<double, 3> > d = 
    arcp< tvmet::Vector<double, 3> >(tvmet_array, 0, size, false);
  ArrayRCP< tvmet::Vector<double, 3> > e = 
    arcp< tvmet::Vector<double, 3> >(&tvmet_array[size], 0, size, false);
  ArrayRCP< tvmet::Vector<double, 3> > f = 
    arcp< tvmet::Vector<double, 3> >(&tvmet_array[2*size], 0, size, false);

  double* raw_array = new double[num_vectors * size * 3];

  double* raw_a = raw_array;
  double* raw_b = &raw_array[size];
  double* raw_c = &raw_array[2*size];

  for (int i=0; i < a.size(); ++i)
    a[i] = 1.0;
  for (int i=0; i < b.size(); ++i)
    b[i] = 2.0;
  for (int i=0; i < c.size(); ++i)
    c[i] = 3.0;
  for (int i=0; i < d.size(); ++i)
    d[i] = 1.0;
  for (int i=0; i < e.size(); ++i)
    e[i] = 2.0;
  for (int i=0; i < f.size(); ++i)
    f[i] = 3.0;
  for (int i=0; i < size; ++i) {
    int offset = i * 3;
    for (int j=0; j < 3; ++j) {
      raw_a[offset + j] = 1.0;
      raw_b[offset + j] = 2.0;
      raw_c[offset + j] = 3.0;
    }
  }

  RCP<Time> vector_time = TimeMonitor::getNewTimer("Vector Time");
  RCP<Time> update_time = TimeMonitor::getNewTimer("Update Time");
  RCP<Time> tvmet_time = TimeMonitor::getNewTimer("TVMET Time");
  RCP<Time> raw_time = TimeMonitor::getNewTimer("Raw Time");
  RCP<Time> raw2_time = TimeMonitor::getNewTimer("Raw2 Time");

  for (int sample = 0; sample < num_samples; ++sample) {
    
    cout << "Vector" << endl;
    {
      TimeMonitor t(*vector_time);
      for (int i=0; i < num_loops; ++i)
	for (int j=0; j < c.size(); ++j)
	  c[j] = a[j] * b[j];
    } 
    
    cout << "Update" << endl;
    {
      TimeMonitor t(*update_time);
      for (int i=0; i < num_loops; ++i)
	for (int j=0; j < c.size(); ++j)
	  c[j].update_multiply(a[j], b[j]);
    }
    
    cout << "TVMET" << endl;
    {
      TimeMonitor t(*tvmet_time);
      for (int i=0; i < num_loops; ++i)
	for (int j=0; j < d.size(); ++j)
	  f[j] = d[j] * e[j];
    }
    
    cout << "Raw" << endl;
    {
      TimeMonitor t(*raw_time);
      for (int i=0; i < num_loops; ++i) {
	for (int j=0; j < size; ++j) {
	  int offset = j * 3;
	  for (int k=0; k < 3; ++k)
	    raw_c[offset + k] = raw_a[offset + k] * raw_b[offset + k];
	}
      }
    }
    
    cout << "Raw2" << endl;
    {
      TimeMonitor t(*raw2_time);
      const int raw_size = 3 * size;  // 3 vector components
      for (int i=0; i < num_loops; ++i) {
	for (int j=0; j < raw_size; ++j) {
	  raw_c[j] = raw_a[j] * raw_b[j];
	}
      }
    }
    
  } // end loop over samples

  TimeMonitor::summarize();
  
  double f_vector = vector_time->totalElapsedTime() / raw_time->totalElapsedTime();
  double f_update = update_time->totalElapsedTime() / raw_time->totalElapsedTime();
  double f_tvmet = tvmet_time->totalElapsedTime() / raw_time->totalElapsedTime();
  double f_raw = raw_time->totalElapsedTime() / raw_time->totalElapsedTime();

  double f_raw2 = raw2_time->totalElapsedTime() / raw_time->totalElapsedTime();

  std::cout << "vector = " << f_vector << std::endl;
  std::cout << "update = " << f_update << std::endl;
  std::cout << "tvmet  = " << f_tvmet << std::endl;
  std::cout << "raw    = " << f_raw << std::endl;
  std::cout << "raw2   = " << f_raw2 << std::endl;

  delete [] vector_array;
  delete [] tvmet_array;
  delete [] raw_array;

  std::cout << "\nTest passed!\n" << std::endl; 
    
  return 0;
}

