//**********************************************************************
template<typename ScalarT, typename Traits> 
Fourier<ScalarT, Traits>::
Fourier(const Teuchos::ParameterList& p) :
  flux("Energy_Flux", 
       p.get< Teuchos::RCP<PHX::DataLayout> >("Vector Data Layout") ),
  density("Density", 
	  p.get< Teuchos::RCP<PHX::DataLayout> >("Scalar Data Layout") ),
  dc("Diffusion Coefficient", 
     p.get< Teuchos::RCP<PHX::DataLayout> >("Scalar Data Layout") ),
  grad_temp("Temperature Gradient", 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Vector Data Layout") )
{ 
  this->addEvaluatedField(flux);
  this->addDependentField(density);
  this->addDependentField(dc);
  this->addDependentField(grad_temp);

  this->setName("Fourier");
}

//**********************************************************************
template<typename ScalarT, typename Traits> 
Fourier<ScalarT, Traits>::~Fourier()
{ }

//**********************************************************************
template<typename ScalarT, typename Traits> 
void Fourier<ScalarT, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& vm)
{
  vm.getFieldData(flux);
  vm.getFieldData(density);
  vm.getFieldData(dc);
  vm.getFieldData(grad_temp);

  data_layout_size = flux.fieldTag().dataLayout()->size();
}

//**********************************************************************
template<typename ScalarT, typename Traits>
void Fourier<ScalarT, Traits>::evaluateFields(typename Traits::EvalData d)
{ 
  std::size_t size = d.size() * data_layout_size;

  for (std::size_t i = 0; i < size; ++i)
    flux[i] = - density[i] * dc[i] * grad_temp[i];
}

//**********************************************************************
