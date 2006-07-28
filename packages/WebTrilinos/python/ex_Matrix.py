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

print Matrix
