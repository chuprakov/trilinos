#! /bin/env python
try:
  import setpath
  import Epetra
except:
  from PyTrilinos import Epetra

Comm = Epetra.PyComm()

print Comm
