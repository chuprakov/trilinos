#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"

//**********************************************************************
template< typename ScalarT, typename Traits>
FEInterpolation<ScalarT, Traits>::
FEInterpolation(const Teuchos::ParameterList& p) :
  val_node(p.get<std::string>("Node Variable Name"), 
	   p.get< Teuchos::RCP<PHX::DataLayout> >("Node Data Layout") ),
  val_qp(p.get<std::string>("QP Variable Name"), 
	 p.get< Teuchos::RCP<PHX::DataLayout> >("QP Data Layout") ),
  val_grad_qp(p.get<std::string>("Gradient QP Variable Name"), 
	      p.get< Teuchos::RCP<PHX::DataLayout> >("Gradient QP Data Layout") )
{ 
  this->addDependentField(val_node);
  this->addEvaluatedField(val_qp);
  this->addEvaluatedField(val_grad_qp);
  
  this->setName("FEInterpolation");
}

//**********************************************************************
template<typename ScalarT, typename Traits>
FEInterpolation<ScalarT, Traits>::~FEInterpolation()
{ }

//**********************************************************************
template<typename ScalarT, typename Traits> 
void FEInterpolation<ScalarT, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& vm)
{
  vm.setFieldData(val_node);
  vm.setFieldData(val_qp);
  vm.setFieldData(val_grad_qp);
}

//**********************************************************************
template<typename ScalarT, typename Traits>
void FEInterpolation<ScalarT, Traits>::
evaluateFields(typename Traits::EvalData cell_data)
{ 
  
  const int nodes_per_cell = val_node.fieldTag().dataLayout()->size();
  const int qp_per_cell = val_qp.fieldTag().dataLayout()->size();

  // Loop over number of cells
  for (std::size_t cell = 0; cell < cell_data.size(); ++cell) {
    
    std::vector<double>& phi = cell_data[cell].getBasisFunctions();
    std::vector< MyVector<double> >& grad_phi = 
      cell_data[cell].getBasisFunctionGradients();
    int node_offset = cell * nodes_per_cell;
    int qp_offset = cell * qp_per_cell;
    
    // Loop over quad points of cell
    for (int qp = 0; qp < qp_per_cell; ++qp) {
      
      val_qp[qp_offset + qp] = 0.0;
      val_grad_qp[qp_offset + qp] = MyVector<ScalarT>(0.0, 0.0, 0.0);
      
      // Sum nodal contributions to qp
      for (int node = 0; node < nodes_per_cell; ++node) {
	
	val_qp[qp_offset + qp] += phi[node] * val_node[node_offset + node];

	val_grad_qp[qp_offset + qp] += 
	  grad_phi[node] * val_node[node_offset + node];
      }      
    }
    
  }
    
}

//**********************************************************************
