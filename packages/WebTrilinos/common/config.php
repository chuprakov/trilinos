<? 

# ---------------------------------------------- #
# This file contains the location of WebTrilinos #
# scripts and files.                             #
#
# NOTE: You must manually modify this file!      #
# ---------------------------------------------- #

# Specifies the location of the WebTrilinos files. This directory
# is contained in your httpdocs directory; these files must be 
# accessible from your web server. 

$WebTrilinosDirectory = "/var/www/html/MatrixPortal/NewMatrixPortal/WebTrilinos/";

# Specifies the location of the image files, as accessible from
# a web browser. This must be a valid HTTP location.

$HTMLImageDirectory = "http://kythira.ethz.ch/MatrixPortal/NewMatrixPortal/WebTrilinos/tmp";

# Specifies the location of the tmp directory that will be used
# to generate image files. This is typically /tmp/

$TempDirectory = "/tmp/";

# Specifies the location of the Trilinos headers. 

$INCLUDES = "-I/home/msala/Trilinos/LINUX_SERIAL/include";

# Specifies the LDFLAGS to be used to link Trilinos examples.

$LDFLAGS = "-L/home/msala/lib -L/home/msala/Trilinos/LINUX_SERIAL/lib";

# Specifies the LIBS to be used to link Trilinos examples.

$LIBS = "-lml -lifpack -laztecoo -lamesos -lgaleri -lepetraext -lepetra -lteuchos -llapack -lblas -lf2c -lm -lMAIN__";

# Specifies the locatin of PyTrilinos, NumPy, and other Python packages.

$PYTHONPATH = "/home/msala/Trilinos/LINUX_SERIAL/lib/python2.4/site-packages/";

# Specifies the C++ compiler to be used.

$CXX = "g++";

# Specifies the host of this version of WebTrilinos

$WebHost = "COLAB/D-INFK of the ETH Zurich";

# ------------------------------------------------ #
# The following variables are already set for most #
# users. Modify the variable below only if you     #
# really need it!                                  #
# ------------------------------------------------ #

# Specifies the location of the HB matrix data. This is typically
# a subdirectory of WebTrilinosDirectory.

$HBMatrixDirectory = $WebTrilinosDirectory . "/data/HB/";

# Specifies the location of the XML matrix data. This is typically
# a subdirectory of WebTrilinosDirectory.

$XMLMatrixDirectory = $WebTrilinosDirectory . "/data/XML/";

# Specifies the location of the tmp directory that will be used
# to generate image files. This is typically /tmp/

$ImageDirectory = $WebTrilinosDirectory . "/tmp/";

?>
