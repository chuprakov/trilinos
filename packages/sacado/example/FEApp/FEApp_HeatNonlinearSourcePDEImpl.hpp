// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

template <typename EvalT>
FEApp::HeatNonlinearSourcePDE<EvalT>::
HeatNonlinearSourcePDE(
  const Teuchos::RCP< const FEApp::AbstractFunction<EvalT> >& mat_func,
  const Teuchos::RCP< const FEApp::AbstractSourceFunction<EvalT> >& src_func) :
  mat(mat_func),
  source(src_func),
  num_qp(0),
  num_nodes(0),
  phi(),
  dphi(),
  jac(),
  u(),
  du(),
  udot(),
  a(),
  f()
{
}

template <typename EvalT>
FEApp::HeatNonlinearSourcePDE<EvalT>::
~HeatNonlinearSourcePDE()
{
}

template <typename EvalT>
unsigned int 
FEApp::HeatNonlinearSourcePDE<EvalT>::
numEquations() const
{
  return 1;
}

template <typename EvalT>
void
FEApp::HeatNonlinearSourcePDE<EvalT>::
init(unsigned int numQuadPoints, unsigned int numNodes)
{
  num_qp = numQuadPoints;
  num_nodes = numNodes;

  phi.resize(num_qp);
  dphi.resize(num_qp);
  jac.resize(num_qp);
  u.resize(num_qp);
  du.resize(num_qp);
  udot.resize(num_qp);
  a.resize(num_qp);
  f.resize(num_qp);

  for (unsigned int i=0; i<num_qp; i++) {
    phi[i].resize(num_nodes);
    dphi[i].resize(num_nodes);
  }
}

template <typename EvalT>
void
FEApp::HeatNonlinearSourcePDE<EvalT>::
evaluateElementResidual(const FEApp::AbstractQuadrature& quadRule,
                        const FEApp::AbstractElement& element,
                        const std::vector<ScalarT>* dot,
                        const std::vector<ScalarT>& solution,
                        std::vector<ScalarT>& residual)
 {
  
  // Quadrature points
  const std::vector<double>& xi = quadRule.quadPoints();

  // Weights
  const std::vector<double>& w = quadRule.weights();

  // Evaluate shape functions
  element.evaluateShapes(xi, phi);

  // Evaluate shape function derivatives
  element.evaluateShapeDerivs(xi, dphi);

  // Evaluate Jacobian of transformation to standard element
  element.evaluateJacobian(xi, jac);

  // Evaluate material values
  mat->evaluate(xi, a);

  // Compute u
  for (unsigned int qp=0; qp<num_qp; qp++) {
    u[qp] = 0.0;
    du[qp] = 0.0;
    udot[qp] = 0.0;
    for (unsigned int node=0; node<num_nodes; node++) {
      u[qp] += solution[node] * phi[qp][node];
      du[qp] += solution[node] * dphi[qp][node];
      if (dot != NULL)
        udot[qp] += (*dot)[node] * phi[qp][node];
    }
  }

  // Evaluate source function
  source->evaluate(u, f);

  // Evaluate residual
  for (unsigned int node=0; node<num_nodes; node++) {
    residual[node] = 0.0;
    for (unsigned int qp=0; qp<num_qp; qp++) {
      residual[node] += 
        w[qp]*jac[qp]*(-(1.0/(jac[qp]*jac[qp]))*a[qp]*du[qp]*dphi[qp][node] + 
                       phi[qp][node]*(f[qp] - udot[qp]));
    }
  }

}
