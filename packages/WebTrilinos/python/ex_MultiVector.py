#! /bin/env python
from PyTrilinos import Epetra, AztecOO, Amesos, IFPACK, ML
from Numeric import *

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
