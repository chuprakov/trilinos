#!/bin/bash

#----------------------------------------------------------------------
# Paths

export PHDMESH_PATH="${HOME}/Trilinos/packages/phdmesh"

#----------------------------------------------------------------------
# The configuration is dependent upon this file

export PHDMESH_CONFIG_DEPS=$0

#----------------------------------------------------------------------
# Compiler and linker configuration:

GNU_PATH="/usr/local/gcc/64Bit/4.0.1/bin"

export PATH="${GNU_PATH}:${PATH}"

source ${PHDMESH_PATH}/config/gnu opt 

#----------------------------------------------------------------------
# Compiler and linker configuration:
#
# GNU_PATH="/usr/local/gcc/32Bit/3.4.3/bin"
# MPI_PATH="/usr/local/mpi/sierra/32Bit/1.2.7/gcc-3.4.3"
#
# export PATH="${MPI_PATH}/bin:${GNU_PATH}:${PATH}"
#
# source ${PHDMESH_PATH}/config/gnu 32bit static opt mpich
#
#----------------------------------------------------------------------
# Compiler and linker configuration:
#
# GNU_PATH="/usr/local/gcc/32Bit/3.4.3/bin"
# MPI_PATH="/usr/local/mpi/sierra/32Bit/1.2.7/gcc-3.4.3"
#
# export PATH="${MPI_PATH}/bin:${GNU_PATH}:${PATH}"
#
# source ${PHDMESH_PATH}/config/gnu 32bit static debug \
#				mpich_purify ${MPI_PATH}/lib
#
#----------------------------------------------------------------------
# Compiler and linker configuration:
#
# GNU_PATH="/usr/local/gcc/32Bit/3.4.3/bin"
# MPI_PATH="/usr/local/mpi/sierra/32Bit/1.2.6..14b-gm-2.0.24/gcc-3.4.3"
#
# export PATH="${MPI_PATH}/bin:${GNU_PATH}:${PATH}"
#
# source ${PHDMESH_PATH}/config/gnu 32bit static opt mpigm /opt/gm/lib
#
#----------------------------------------------------------------------
# SNL ACCESS ExodusII and NemesisI configuration, only for 32Bit :
#
# source ${PHDMESH_PATH}/config/exodusII  /usr/local/eng_sci/struct/i686/current-gcc
#
#----------------------------------------------------------------------

echo `which ${CXX}` ${CXXFLAGS}
echo `which ${CC}` ${CFLAGS}

make -f ${PHDMESH_PATH}/Make.in $*

#----------------------------------------------------------------------

