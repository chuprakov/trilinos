#! /bin/env python
try:
  import setpath
  import Epetra, IFPACK, AztecOO
except:
  from PyTrilinos import Epetra, AztecOO, IFPACK

Comm = Epetra.PyComm()

NumGlobalElements = 10
Map = Epetra.Map(NumGlobalElements, 0, Comm)
MyGlobalElements = Map.MyGlobalElements()

Matrix = Epetra.CrsMatrix(Epetra.Copy, Map, 0)

for i in MyGlobalElements:
  Matrix[i, i] = 1.0

Matrix.FillComplete()


Factory = IFPACK.Factory()
Prec = Factory.Create("ILU", Matrix)
IFPACKList = {
  "fact: level-of-fill": 1
}

Prec.SetParameters(IFPACKList)
Prec.Initialize()
Prec.Compute()

