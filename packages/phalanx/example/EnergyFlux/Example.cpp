#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "CellData.hpp"
#include "Traits.hpp"
#include "FactoryTraits.hpp"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template<typename ScalarT>
void compareScalarFields(PHX::Field<ScalarT>& field1, 
			 PHX::Field<ScalarT>& field2,
			 double tol)
{
  std::cout << "Comparing scalar fields\n" << field1.fieldTag() << "\n" 
	    << field2.fieldTag() << std::endl; 

  TEST_FOR_EXCEPTION(field1.size() != field2.size(), std::logic_error,
		     "Fields for comparison do not have the same size!");

  double error = 0.0;
  for (int i=0; i < field1.size(); ++i)
    error += fabs(field1[i]-field2[i]);
  
  TEST_FOR_EXCEPTION(error > tol, std::runtime_error,
		     "Fields are not equal in comparison!");

  std::cout << "Passed: " << error << " < " << tol << std::endl; 
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template<typename ScalarT>
void compareVectorFields(PHX::Field< MyVector<ScalarT> >& field1, 
			 PHX::Field< MyVector<ScalarT> >& field2,
			 double tol)
{
  std::cout << "Comparing vector fields\n" << field1.fieldTag() << "\n" 
	    << field2.fieldTag() << std::endl; 

  TEST_FOR_EXCEPTION(field1.size() != field2.size(), std::logic_error,
		     "Fields for comparison do not have the same size!");

  double error = 0.0;
  for (int i=0; i < field1.size(); ++i)
    for (int j=0; j < 3; ++j)
      error += fabs(field1[i][j]-field2[i][j]);
  
  TEST_FOR_EXCEPTION(error > tol, std::runtime_error,
		     "Fields are not equal in comparison!");

  std::cout << "Passed: " << error << " < " << tol << std::endl; 
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of FieldManager testing
    // *********************************************************************
    {
      cout << "\nStarting FieldManager Testing !" << endl;

      RCP<DataLayout> scalar_qp = 
	rcp(new Generic<MyTraits::MY_SCALAR>("Q1_QP", 4));
      RCP<DataLayout> vector_qp = 
	rcp(new Generic<MyTraits::MY_VECTOR>("Q1_QP", 4));
      RCP<DataLayout> scalar_node = 
	rcp(new Generic<MyTraits::MY_SCALAR>("Q1_NODE", 4));

      // Parser will build parameter list that determines the field
      // evaluators to build
      map<string, RCP<ParameterList> > evaluators_to_build;
      
      { // Temperature
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Temperature");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", scalar_node);
	evaluators_to_build["DOF_Temperature"] = p;
      }
      { // Density
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_density;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", scalar_qp);
	evaluators_to_build["Density"] = p;
      }

      { // Constant Diffusion Coefficient
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Diffusion Coefficient");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", scalar_qp);
	evaluators_to_build["Diffusion Coefficient"] = p;
      }
      
      { // Nonlinear Source
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_nonlinearsource;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", scalar_qp);
	evaluators_to_build["Nonlinear Source"] = p;
      }

      { // Fourier Energy Flux
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_fourier;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Scalar Data Layout", scalar_qp);
	p->set< RCP<DataLayout> >("Vector Data Layout", vector_qp);
	evaluators_to_build["Energy Flux"] = p;
      }

      { // FE Interpolation
	RCP<ParameterList> p = rcp(new ParameterList);

	int type = MyFactoryTraits<MyTraits>::id_feinterpolation;
	p->set<int>("Type", type);

	p->set<string>("Node Variable Name", "Temperature");
	p->set<string>("QP Variable Name", "Temperature");
	p->set<string>("Gradient QP Variable Name", "Temperature Gradient");

	p->set< RCP<DataLayout> >("Node Data Layout", scalar_node);
	p->set< RCP<DataLayout> >("QP Data Layout", scalar_qp);
	p->set< RCP<DataLayout> >("Gradient QP Data Layout", vector_qp);

	evaluators_to_build["FE Interpolation"] = p;
      }

      // Build Field Evaluators
      EvaluatorFactory<MyTraits,MyFactoryTraits<MyTraits> > factory;
      RCP< vector< RCP<Evaluator_TemplateManager<MyTraits> > > > 
	evaluators;
      evaluators = factory.buildEvaluators(evaluators_to_build);
 
          
      // Request quantities to assemble PDE operators
      FieldManager<MyTraits> vm;
      FieldTag energy_flux("Energy_Flux", vector_qp);
      vm.requireFieldForAllTypes(energy_flux);
      FieldTag source("Nonlinear Source", scalar_qp);
      vm.requireFieldForAllTypes(source);
      
      // Register all Evaluators 
      registerEvaluators(evaluators, vm);

      const std::size_t num_cells = 10;
      const std::size_t num_eval_loops = 1;

      RCP<Time> registration_time = 
	TimeMonitor::getNewTimer("Post Registration Setup Time");
      {
	TimeMonitor t(*registration_time);
	vm.postRegistrationSetup(num_cells);
      }

      cout << vm << endl;
      
      std::vector<CellData> cells(num_cells);

      RCP<Time> eval_time = TimeMonitor::getNewTimer("Evaluation Time");

      vm.preEvaluate<double>(NULL);
      {
	TimeMonitor t(*eval_time);
	for (std::size_t i=0; i < num_eval_loops; ++i)
	  vm.evaluateFields<double>(cells);
      }
      vm.postEvaluate<double>(NULL);

      // Test data retrieval
      cout << "Testing data members" << endl;
      FieldTag d_var("Density", scalar_qp);
      Field<double> den(d_var); 
      vm.getFieldData(den);
      cout << "size of density = " << den.size() << ", should be " 
	   << num_cells * d_var.dataLayout()->size() << "." << endl;
      TEST_FOR_EXCEPTION(den.size() != static_cast<Teuchos::ArrayRCP<double>::Ordinal>(num_cells * d_var.dataLayout()->size()),
			 std::runtime_error, 
			 "Returned arrays are not sized correctly!");
      
      
      cout << endl;

      // Compare temperature fields, should be 2.0
      Field<double> temp("Temperature", scalar_node);
      vm.getFieldData(temp);
      
      Field<double> temp_base("Temperature Baseline", scalar_node);
      ArrayRCP<double> temp_base_data = 
	arcp<double>(num_cells * scalar_node->size());
      temp_base.setFieldData(temp_base_data);
      for (int i=0; i<temp_base.size(); ++i)
	temp_base[i] = 2.0;
      
      compareScalarFields(temp, temp_base, 1.0e-12);

      cout << endl;

      // Compare temperature gradient fields, should be 2.0
      Field< MyVector<double> > tg("Temperature Gradient", vector_qp);
      vm.getFieldData(tg);

      Field< MyVector<double> > 
	tg_base("Temperature Gradient Baseline", vector_qp);
      ArrayRCP< MyVector<double> > tg_base_data = 
	arcp< MyVector<double> >(num_cells * vector_qp->size());
      tg_base.setFieldData(tg_base_data);
      for (int i=0; i<tg_base.size(); ++i)
	tg_base[i] = 2.0;
      
      compareVectorFields(tg, tg_base, 1.0e-12);

      cout << endl;

      // Compare energy flux fields, should be -16.0
      Field< MyVector<double> > ef("Energy_Flux", vector_qp);
      vm.getFieldData(ef);

      Field< MyVector<double> > 
	ef_base("Energy_Flux Baseline", vector_qp);
      ArrayRCP< MyVector<double> > ef_base_data = 
	arcp< MyVector<double> >(num_cells * vector_qp->size());
      ef_base.setFieldData(ef_base_data);
      for (int i=0; i<ef_base.size(); ++i)
	ef_base[i] = -16.0;
      
      compareVectorFields(ef, ef_base, 1.0e-12);

      cout << endl;

    }
    
    // *********************************************************************
    // Finished all testing
    // *********************************************************************
    std::cout << "\nRun has completed successfully!\n" << std::endl; 
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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
