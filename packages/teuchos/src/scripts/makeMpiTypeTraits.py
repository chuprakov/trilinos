# @HEADER
# ***********************************************************************
#
#                    Teuchos: Common Tools Package
#                 Copyright (2004) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ***********************************************************************
# @HEADER

'''Generate C++ code for converting from a C++ type to its MPI_Datatype.

Author: Mark Hoemmen <mhoemme@sandia.gov>
Date: Aug, Sep 2012

Introduction
============

This module generates C++ code for MpiTypeTraits, a traits class
which, for a given C++ type Packet, returns the corresponding
MPI_Datatype.  The MPI_Datatype tells MPI (Message Passing Interface)
how to send and receive objects of type Packet.  The traits class also
tells users whether or not the MPI_Datatype must be freed (via
MPI_Type_free()) after use.

Usage
=====

If you run this code as a script with no arguments, it will write two
files to the current working directory:

Teuchos_MpiTypeTraits.hpp
Teuchos_MpiTypeTraits.cpp

The first file contains the generic declaration and documentation of
MpiTypeTraits, declarations of full specializations of MpiTypeTraits,
and declarations and definitions of some partial specializations.
The second file contains definitions of the full specializations of
MpiTypeTraits.

After creating these C++ files, one should copy them manually into the
parent directory (teuchos/src).  They should also be listed in
teuchos/src/CMakeLists.txt, conditionally if MPI is enabled.  Look for
"IF (${PACKAGE_NAME}_ENABLE_MPI)", and append Teuchos_MpiTypeTraits.hpp
to the HEADERS set and Teuchos_MpiTypeTraits.cpp to the SOURCES set.'''

# We actually reuse templates, which is why we use string.Template
# instead of "new"-style format statements.  This is hopefully more
# efficient.
from string import Template
# Redirect stdout to a file.
import sys

def emitCopyrightNotice ():
    '''Print the Teuchos package's copyright notice.
    This belongs at the top of every Teuchos header or source file.'''
    
    print '''// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
'''

def makeTypeMap ():
    '''Return the map from C++ Packet type to the corresponding MPI_Datatype.

Each key is a C++ type Packet for which we are generating a specialization
of MpiTypeTraits<Packet>.  Its corresponding value is a (string, integer) pair:
first, the MPI_Datatype which tells MPI how to send and receive data of that
type, and second, the count of consecutive MPI_Datatype to send or receive.

The map only includes keys for C++ types for which the corresponding
MPI_Datatype is a "basic" (not derived) type.  A derived MPI_Datatype
must be freed after use via MPI_Type_free(); a "basic" MPI_Datatype
need not and must not be freed after use, and comes built-in.'''

    # Start with an empty dictionary; fill it below.
    d = {}

    # Many built-in C++ types have a simple conversion to their
    # standard MPI_Datatype values: Capitalize, replace each space
    # with an underscore, and prepend MPI_.
    simpleNames = ['char', 'signed char', 'unsigned char',
		   'short', 'unsigned short',
		   'int', 'unsigned int',
		   'long', 'unsigned long',
		   'long long', 'unsigned long long', \
		   'size_t', 'ptrdiff_t',
		   'float', 'double', 'long double']
    for name in simpleNames:
        d.update ({name, ('MPI_' + name.upper ().replace (' ', '_'), 1)})

    # Other C++ types don't have a straightforward conversion.
    # MPI_CXX_COMPLEX is in the MPI 3.0 standard (thanks in small part
    # to Jeff Hammond and Jed Brown thinking of me!) but
    # MPI_C_COMPLEX, MPI_C_FLOAT_COMPLEX, MPI_C_DOUBLE_COMPLEX, and
    # MPI_C_LONG_DOUBLE_COMPLEX work perfectly well, in case your MPI
    # implementation does not support the 3.0 standard.  We deal with
    # bool elsewhere as a separate case, since there is no basic MPI
    # datatype which standardly supports bool.
    d.update ({'std::complex<float>': ('MPI_C_FLOAT_COMPLEX', 1), \
	       'std::complex<double>': ('MPI_C_DOUBLE_COMPLEX', 1), \
	       'std::complex<long double>': ('MPI_C_LONG_DOUBLE_COMPLEX', 1)})
    for name in simpleNames:
        d.update ({'std::pair<' + name + ', ' + name + '>', d[name]})
    return d

def emit (s, indent=0):
    '''Print code to stdout, with the proper indent level.

    s (string): the string to print.  It may contain one or more lines.
    indent (integer): the number of spaces to indent each line.'''
    for line in s.split('\n'):
	print ' '*indent + line

def emitGenericDecl (indent=0):
    '''Emit the (generic) declaration of MpiTypeTraits, with documentation.

    indent (integer): the number of spaces to indent.'''
    
    with open('MpiTypeTraits.hpp', 'r') as infile:
	for line in infile:
	    emit (line, indent)
	
def emitSpecializationDecls (d, indent=0):
    '''Emit declarations of the specializations of MpiTypeTraits.

    d (dictionary): return value of makeTypeMap().
    indent (integer): the number of spaces to indent.'''
    
    # It's a template template.  Whee!!!
    # Fill it in for each C++ type key of the dictionary.
    tmpl = Template('''// Specialization of MpiTypeTraits<T> for T=${key}.
template<>
class MpiTypeTraits< ${key} > {
 public:
  typedef ${key} packet_type;
  static const bool needFree = false;
  static const bool globallyConsistent = true;
  static const bool globallyConsistentType = true;

  static std::pair<MPI_Datatype, size_t> makeType (const ${key}& example);
};
''')
    for k,v in d.iteritems():
        if k != 'long long' and k != 'unsigned long long' and k != 'size_t' and k != 'ptrdiff_t': # we handle these specially
	    emit (tmpl.substitute (key=k), indent)
    emit ('''#if defined(HAVE_TEUCHOS_LONG_LONG_INT)

// Partial specializations for long long and unsigned long long.
// On platforms with sizeof(size_t) <= sizeof(unsigned long long), 
// this should take care of the size_t and ptrdiff_t specializations
// as well, since we've covered all built-in unsigned integer types
// above with size <= sizeof(unsigned long long).
''', indent)
    emit (tmpl.substitute (key='long long', value=d['long long'][0]), indent)
    emit (tmpl.substitute (key='unsigned long long', \
			   value=d['unsigned long long'][0]), indent)
    emit ('''// The C preprocessor does not allow "sizeof(T)" expressions in #if
// statements, even if T is a built-in type.  Otherwise, we could test
// for 'sizeof(size_t) > sizeof(unsigned long int)'.  The constants
// below are defined in the <cstdint> header file.
#elif SIZE_MAX > ULONG_MAX

// We already have an unsigned long int specialization above.  If
// Teuchos support for "long long" is enabled, then we've taken care
// of all possible lengths of size_t: unsigned (char, short, int,
// long, long long).  If "long long" is _not_ enabled, we need to
// check if sizeof(size_t) > sizeof(unsigned long).  If so, then we
// need a specialization for size_t.  Ditto for ptrdiff_t (which is a
// signed type of the same length as size_t).
''', indent)
    emit (tmpl.substitute (key='ptrdiff_t', value=d['ptrdiff_t'][0]), indent)
    emit (tmpl.substitute (key='size_t', value=d['size_t'][0]), indent)
    emit ('#endif // HAVE_TEUCHOS_LONG_LONG_INT\n', indent)

def emitSpecializationDefs (d, indent=0):
    '''Emit definitions of the specializations of MpiTypeTraits.

    d (dictionary): return value of makeTypeMap().
    indent (integer): the number of spaces to indent.'''
    
    # It's a template template.  Whee!!!
    # Fill it in for each (C++ type, MPI_Datatype) pair in the dictionary.
    tmpl = Template('''// Specialization of MpiTypeTraits<T> for T=${key}.
template<>
MPI_Datatype MpiTypeTraits< ${key} >::makeType () {
  return std::make_pair (${type}, static_cast<size_t> (${count}));
}
''')
    for k,v in d.iteritems():
        if k != 'long long' and k != 'unsigned long long' and k != 'size_t' and k != 'ptrdiff_t': # we handle these specially
	    emit (tmpl.substitute (key=k, type=v[0], count=v[1]), indent)
    emit ('''#if defined(HAVE_TEUCHOS_LONG_LONG_INT)

// Partial specializations for long long and unsigned long long.
// On platforms with sizeof(size_t) <= sizeof(unsigned long long), 
// this should take care of the size_t and ptrdiff_t specializations
// as well, since we've covered all built-in unsigned integer types
// above with size <= sizeof(unsigned long long).
''', indent)
    emit (tmpl.substitute (key='long long', \
                           type=d['long long'][0], \
                           count=d['long long'][1]), indent)
    emit (tmpl.substitute (key='unsigned long long', \
                           type=d['unsigned long long'][0]),
                           count=d['unsigned long long'][1]), indent)
    emit ('''// The C preprocessor does not allow "sizeof(T)" expressions in #if
// statements, even if T is a built-in type.  Otherwise, we could test
// for 'sizeof(size_t) > sizeof(unsigned long int)'.  The constants
// below are defined in the <cstdint> header file.
#elif SIZE_MAX > ULONG_MAX

// We already have an unsigned long int specialization above.  If
// Teuchos support for "long long" is enabled, then we've taken care
// of all possible lengths of size_t: unsigned (char, short, int,
// long, long long).  If "long long" is _not_ enabled, we need to
// check if sizeof(size_t) > sizeof(unsigned long).  If so, then we
// need a specialization for size_t.  Ditto for ptrdiff_t (which is a
// signed type of the same length as size_t).
''', indent)
    emit (tmpl.substitute (key='ptrdiff_t', \
			   type=d['ptrdiff_t'][0], \
	                   count=d['ptrdiff_t'][1]), indent)
    emit (tmpl.substitute (key='size_t', \
			   type=d['size_t'][0], \
			   count=d['size_t'][1]), indent)
    emit ('#endif // HAVE_TEUCHOS_LONG_LONG_INT\n', indent)

def headerizeFilename (filename):
    '''Convert filename to a suitable string for a header file's #ifndef ... #endif.

    filename (string): Name of the header file.'''
    return filename.upper ().replace ('.', '_')

def emitStdPairDeclAndDef (indent=0):
    '''Emit decl. and def. of MpiTypeTraits<std::pair<X,Y> > partial specialization.

    indent (integer): the number of spaces to indent.'''

    with open('MpiTypeTraitsPair.hpp', 'r') as infile:
	for line in infile:
	    emit (line, indent)

def emitBoolDecl (indent=0):
    '''Emit declaration of the specialization of MpiTypeTraits for bool.

    indent (integer): the number of spaces to indent.'''

    with open('MpiTypeTraitsBool.hpp', 'r') as infile:
	for line in infile:
	    emit (line, indent)

def emitDeclBody (d, indent=0):
    emitGenericDecl (indent)
    emitSpecializationDecls (d, indent)
    emitBoolDecl (indent)
    emitStdPairDeclAndDef (indent)

def emitBoolDef (indent=0):
    '''Emit definition of the specialization of MpiTypeTraits for bool.

    indent (integer): the number of spaces to indent.'''

    with open('MpiTypeTraitsBool.cpp', 'r') as infile:
	for line in infile:
	    emit (line, indent)

def emitDefBody (d, indent=0):
    emitSpecializationDefs (d, indent)
    emitBoolDef (indent)

def emitDeclHeader (d, filename, indent=0):
    '''Emit the entire file of MpiTypeTraits declarations (a .hpp file).

    d (dictionary): return value of makeTypeMap().
    filename (string): name of the header file.
    indent (integer): the number of spaces to indent.'''

    hf = headerizeFilename (filename)
    
    emitCopyrightNotice ()
    topMatter = '''#ifndef {0}
#define {0}

#include <mpi.h>
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <utility> // std::pair

/// \\file {1}
/// \\brief Declaration of MpiTypeTraits and its specializations
/// \\warning This file is automatically generated from a script.
///   Any desired changes should be made to the script, not to this
///   file.  Changes to this file will be lost the next time the
///   script is run.
'''.format (hf, filename) + '''namespace Teuchos {
namespace details {
'''
    emit (topMatter, indent)
    emitDeclBody (d, indent)
    emit ('} // namespace details\n} // namespace Teuchos\n', indent)
    emit ('#endif // #ifndef {0}'.format (hf), indent)

def emitDefSource (d, filename, indent=0):
    '''Emit the entire file of MpiTypeTraits definitions (a .cpp file).

    d (dictionary): return value of makeTypeMap().
    filename (string): name of the header file.
    indent (integer): the number of spaces to indent.'''

    emitCopyrightNotice ()
    topMatter = '''/// \\file {0}
/// \\brief Definition of specializations of MpiTypeTraits
/// \\warning This file is automatically generated from a script.
///   Any desired changes should be made to the script, not to this
///   file.  Changes to this file will be lost the next time the
///   script is run.
'''.format (filename) + '''namespace Teuchos {
namespace details {
'''
    emit (topMatter, indent)
    emitDefBody (d, indent)
    emit ('} // namespace details\n} // namespace Teuchos\n', indent)

def run ():
    '''Generate the two header files mentioned in the module's documentation.

    This writes the header file of MpiTypeTraits declarations
    'Teuchos_MpiTypeTraits.hpp', and the source file of definitions
    'Teuchos_MpiTypeTraits.cpp', for all specializations of
    MpiTypeTraits that this module knows how to generate.  Both files
    are written to the current working directory.'''

    d = makeTypeMap ()
    declsFilename = 'Teuchos_MpiTypeTraits.hpp'
    defsFilename = 'Teuchos_MpiTypeTraits.cpp'

    with open(declsFilename, 'w') as outfile:
        try:
	    sys.stdout = outfile # Rebind stdout to the output file.
	    emitDeclHeader (d, declsFilename)
	finally:
	    # If something went wrong, be sure to restore the original
	    # binding of stdout.
	    sys.stdout = sys.__stdout__

    with open(defsFilename, 'w') as outfile:
	try:
	    sys.stdout = outfile
	    emitDefSource (d, defsFilename)
	finally:
	    sys.stdout = sys.__stdout__


# If running the module as an executable script, just call run().
if __name__ == "__main__":
    if len (sys.argv) > 1:
        raise ValueError ('This script does not currently take any command-line arguments.')
    else:
        run ()
