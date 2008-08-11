// -*- c++ -*-

%module PySundance

%feature("autodoc");

%exception 
{
  try
    {
      $action
    }
  catch (std::exception& e)
    {
      PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(e.what()));
      return NULL;
    }
}


%{
#include "Sundance.hpp"
#include "SundancePathUtils.hpp"
  %}


%inline %{
  bool passFailTest(double err, double tol)
  {
    return SundanceStdFwk::Sundance::passFailTest(err, tol);
  }


  void skipTimingOutput() {Sundance::skipTimingOutput()=true;}

  %}


%include "std_string.i"
namespace SundanceUtils
{
  std::string searchForFile(const std::string& name);
}


%include Mesh.i

%include Utils.i

%include ParameterList.i

%include CellFilter.i

%include TSF.i

%include Quadrature.i

%include Basis.i

%include Spectral.i

%include Symbolics.i

%include Integral.i

%include LinearProblem.i

%include NonlinearProblem.i

%include Viz.i

%include Discrete.i

%include AToC.i  
