#! /bin/env python
try:
  import setpath
  import Epetra
except:
  from PyTrilinos import Epetra

Comm = Epetra.PyComm()

NumGlobalElements = 10
Map = Epetra.Map(NumGlobalElements, 0, Comm)
MyGlobalElements = Map.MyGlobalElements()

print MyGlobalElements
