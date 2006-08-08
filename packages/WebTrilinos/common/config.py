# All the settings below refer to MatrixPortal automatically
# generated problems. They can be used to limit the problem size 
# and prevent users to demand to build huge matrices, which 
# may crash the server.
MAX_MATRIX_ROWS = 100000
MAX_MATRIX_NONZEROS = MAX_MATRIX_ROWS * 100
MAX_ITERATIONS = 1550
MAX_KSPACE = 60
# This is the location of HB linear systems
HB_REPOSITORY = "/var/www/html/WebTrilinos/data/HB/"
# This is the location of MM linear systems
MM_REPOSITORY = "/var/www/html/WebTrilinos/data/MM/"
# This is the location of XML linear systems
XML_REPOSITORY = "/var/www/html/WebTrilinos/data/XML/"
