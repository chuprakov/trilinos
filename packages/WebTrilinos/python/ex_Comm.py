#! /bin/env python
from PyTrilinos import Epetra, AztecOO, Amesos, IFPACK, ML
from Numeric import *

Comm = Epetra.PyComm()

print Comm
