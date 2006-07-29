#! /bin/env python
try:
  import setpath
  import Epetra, Amesos
except:
  from PyTrilinos import Epetra, Amesos

Comm = Epetra.PyComm()

NumGlobalElements = 10
Map = Epetra.Map(NumGlobalElements, 0, Comm)
MyGlobalElements = Map.MyGlobalElements()

Matrix = Epetra.CrsMatrix(Epetra.Copy, Map, 0)

for i in MyGlobalElements:
  Matrix[i, i] = 1.0

Matrix.FillComplete()

ExactSolution = Epetra.MultiVector(Map, 1)
ExactSolutionView = ExactSolution.ExtractView()
n = ExactSolution.MyLength()
for i in xrange(0, n):
  ExactSolutionView[0][i] = sin(i * 3.1415 / n) * sin(i * 3.1415 / n)

LHS = Epetra.MultiVector(Map, 1)
RHS = Epetra.MultiVector(Map, 1)
Matrix.Multiply(False, ExactSolution, RHS)

Problem = Epetra.LinearProblem(Matrix, LHS, RHS);

Factory = Amesos.Factory();
Solver = Factory.Create("Amesos_Klu", Problem);

AmesosList = {
  "PrintTiming":         True,
  "PrintStatus":         True,
}
Solver.SetParameters(AmesosList);

Solver.SymbolicFactorization()
Solver.NumericFactorization()
Solver.Solve();
      
del Solver
