#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"

// From test/Utilities directory
#include "Traits.hpp"

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of Field Tag Testing
    // *********************************************************************
    {

      // Dummy data layouts
      RCP<DataLayout> node4 = 
	rcp(new Generic<MyTraits::MY_SCALAR>("Q1_Nodes", 4));
      RCP<DataLayout> quad4 = 
	rcp(new Generic<MyTraits::MY_SCALAR>("Q1_QuadPoints", 4));
      RCP<DataLayout> gradQuad4 = 
	rcp(new Generic<MyTraits::MY_VECTOR>("Q1_QuadPoints", 4));
      
      // Tags with same name but different topology
      FieldTag nodal_density("density", node4);
      FieldTag qp_density("density", quad4);
      FieldTag grad_qp_density("density", gradQuad4);
      

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Ctors
      cout << "Testing ctor with FieldTag...";
      Field<double> a(nodal_density);
      Field< MyVector<double> > b(grad_qp_density);
      cout << "passed!" << endl;
      
      cout << "Testing ctor with individual data...";
      Field<MyTraits::FadType> c("density", node4);
      Field< MyVector<MyTraits::FadType> > d("density", gradQuad4);
      cout << "passed!" << endl;
      
      cout << "Testing empty ctor...";
      Field<double> e;
      Field< MyVector<MyTraits::FadType> > f;
      cout << "passed!" << endl;
      
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // FieldTag accessor
      cout << "Testing fieldTag() accessor...";

      const FieldTag& test_a = a.fieldTag();
      TEST_FOR_EXCEPTION( !(test_a == nodal_density),
			 std::logic_error,
			 "fieldTag() accessor failed!");
      
      const FieldTag& test_b = b.fieldTag();
      TEST_FOR_EXCEPTION( !(test_b == grad_qp_density),
			 std::logic_error,
			 "fieldTag() accessor failed!");
      
      const FieldTag& test_c = c.fieldTag();
      TEST_FOR_EXCEPTION( !(test_c == nodal_density),
			 std::logic_error,
			 "fieldTag() accessor failed!");

      const FieldTag& test_d = d.fieldTag();
      TEST_FOR_EXCEPTION( !(test_d == grad_qp_density),
			 std::logic_error,
			 "fieldTag() accessor failed!");

      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // setFieldTag()
      cout << "Testing setFieldTag() accessor...";
      e.setFieldTag(nodal_density);
      f.setFieldTag(grad_qp_density);
      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // setFieldData()
      cout << "Testing getFieldData() accessor...";
      const int size = 100;
      ArrayRCP<double> a_scalar_scalar = 
	arcp<double>(size);
      ArrayRCP< MyVector<double> > b_vector_scalar = 
	arcp< MyVector<double> >(size);
      ArrayRCP<MyTraits::FadType> c_scalar_fad = 
	arcp<MyTraits::FadType>(size);
      ArrayRCP< MyVector<MyTraits::FadType> > d_vector_fad = 
	arcp< MyVector<MyTraits::FadType> >(size);
      ArrayRCP<double> e_scalar_scalar = 
	arcp<double>(size);
      ArrayRCP< MyVector<MyTraits::FadType> > f_vector_fad = 
	arcp< MyVector<MyTraits::FadType> >(size);

      a.setFieldData(a_scalar_scalar);
      b.setFieldData(b_vector_scalar);
      c.setFieldData(c_scalar_fad);
      d.setFieldData(d_vector_fad);
      e.setFieldData(e_scalar_scalar);
      f.setFieldData(f_vector_fad);
      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // size()
      cout << "Testing size() method...";
      TEST_FOR_EXCEPTION( a.size() != size , std::logic_error, "Size of array a is not equal to requested size.");
      TEST_FOR_EXCEPTION( b.size() != size , std::logic_error, "Size of array b is not equal to requested size.");
      TEST_FOR_EXCEPTION( c.size() != size , std::logic_error, "Size of array c is not equal to requested size.");
      TEST_FOR_EXCEPTION( d.size() != size , std::logic_error, "Size of array d is not equal to requested size.");
      TEST_FOR_EXCEPTION( e.size() != size , std::logic_error, "Size of array e is not equal to requested size.");
      TEST_FOR_EXCEPTION( f.size() != size , std::logic_error, "Size of array f is not equal to requested size.");

      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // operator[]
      cout << "Testing operator[]() accessor...";
      
      for (int i = 0; i < a.size(); ++i) {
	a[i] = 5.0;
	b[i] = 5.0;
	c[i] = 5.0;
	d[i] = MyVector<MyTraits::FadType>(5.0);
	e[i] = 5.0;
	f[i] = MyVector<MyTraits::FadType>(5.0);
      }

      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // ostream
      cout << "Testing operator<<()...";
      ostringstream output;
      output << a << endl;
      cout << "passed!" << endl;
      //cout << output.str() << endl; 
    }

    // *********************************************************************
    // *********************************************************************
    std::cout << "\nTest passed!\n" << std::endl; 
    // *********************************************************************
    // *********************************************************************

  }
  catch (const std::exception& e) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Exception Caught!" << endl;
    std::cout << "Error message is below\n " << e.what() << endl;
    std::cout << "************************************************" << endl;
  }
  catch (...) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Unknown Exception Caught!" << endl;
    std::cout << "************************************************" << endl;
  }

  TimeMonitor::summarize();
    
  return 0;
}
