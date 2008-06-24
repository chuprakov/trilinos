#include "Phalanx_ConfigDefs.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"

// From test/Utilities directory
#include "Traits.hpp"

#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_Generic.hpp"

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of Data Layout Testing
    // *********************************************************************
    {

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // ctor
      cout << "\nTesting constructor...";
      
      RCP<DataLayout> node4 = 
	rcp(new Generic<MyTraits::MY_SCALAR>("Q1_Nodes", 4));
      
      RCP<DataLayout> quad4 = 
	rcp(new Generic<MyTraits::MY_SCALAR>("Q1_QuadPoints", 4));
      
      RCP<DataLayout> gradQuad4 = 
	rcp(new Generic<MyTraits::MY_VECTOR>("Q1_QuadPoints", 4));
      
      cout << "passed!" << endl;

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // name()
      {
	cout << "Testing name() accessor...";
	RCP< Generic<MyTraits::MY_SCALAR> > g_node4 = 
	  rcp_dynamic_cast< Generic<MyTraits::MY_SCALAR> >(node4);
	TEST_FOR_EXCEPTION(is_null(g_node4), 
			   std::logic_error,
			   "dynamic cast from DataLayout to Generic failed!");
	
	TEST_FOR_EXCEPTION(g_node4->name() != std::string("Q1_Nodes"), 
			   std::logic_error,
			   "name() accessor failed!");
	cout << "passed!" << endl;
      }
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // size()
      {
	cout << "Testing size() accessor...";
	TEST_FOR_EXCEPTION(node4->size() != 4, 
			   std::logic_error,
			   "size() accessor failed!");
	cout << "passed!" << endl;
      }
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // operator==()
      {
	cout << "Testing operator==()...";
	
	RCP<DataLayout> unique_node4_copy = 
	  rcp(new Generic<MyTraits::MY_SCALAR>("Q1_Nodes", 4));
	
	
	// same data layout, different object
	TEST_FOR_EXCEPTION( !(*node4 == *unique_node4_copy), 
			    std::logic_error,
			    "operator==() failed!");
	
	// different data layouts - name differentiation
	TEST_FOR_EXCEPTION( *node4 == *quad4, 
			    std::logic_error,
			    "operator==() failed!");
	
	// same name, different algebraic type 
	TEST_FOR_EXCEPTION( *quad4 == *gradQuad4, 
			    std::logic_error,
			    "operator==() failed!");
	
	cout << "passed!" << endl;
      }

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // getAlgebraicTypeInfo
      {
	cout << "Testing getAlgebraicTypeInfo()...";

	RCP< Generic<MyTraits::MY_SCALAR> > g_node4 = 
	  rcp_dynamic_cast< Generic<MyTraits::MY_SCALAR> >(node4);
	
	TEST_FOR_EXCEPTION( node4->getAlgebraicTypeInfo().name() != 
			    quad4->getAlgebraicTypeInfo().name(), 
			    std::logic_error,
			    "getAlgebraicTypeInfo() comparison failed!");

	TEST_FOR_EXCEPTION( node4->getAlgebraicTypeInfo().name() == 
			    gradQuad4->getAlgebraicTypeInfo().name(), 
			    std::logic_error,
			    "getAlgebraicTypeInfo() comparison failed!");

	cout << "passed!" << endl;
      }
      
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // ostream
      cout << "Testing ostream...";
      ostringstream output;
      output << *node4;
      cout << "...passed: " << output.str() << endl;

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
