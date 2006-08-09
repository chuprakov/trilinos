# @HEADER
# ************************************************************************
#
#                WebTrilinos: A Web Interface to Trilinos
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
#
# ************************************************************************
# @HEADER

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
# This is the location of MatrixMarket linear systems
MM_REPOSITORY = "/var/www/html/WebTrilinos/data/MM/"
# This is the location of XML linear systems
XML_REPOSITORY = "/var/www/html/WebTrilinos/data/XML/"
