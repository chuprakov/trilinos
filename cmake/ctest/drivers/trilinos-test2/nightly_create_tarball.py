#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER

from optparse import OptionParser
import os
import sys
import shutil
import re

#Make reporting errors easier to make consistent between all system calls
def reportError(errorCode, errorMessage, errorFileName="", verboseError=False):
  if errorCode != 0:
    print errorMessage
    if verboseError and errorFileName != "":
      errorFile = open(errorFileName)
      for line in errorFile:
        print line.strip()
      errorFile.close()
    sys.exit(1)
#end reportError

workingDir = os.getcwd()

parser = OptionParser()

repositoryDefault="software.sandia.gov:/space/git/nightly/Trilinos"

parser.add_option("-r", "--repository", dest="repository", action="store", default=repositoryDefault, 
  help="Sets the repository to pull Trilinos from. Default=%default")
parser.add_option("--cxx-compiler", dest="cxx", action="store", default="g++", 
  help="Sets the C++ compiler to use. Default=%default")
parser.add_option("--c-compiler", dest="cc", action="store", default="gcc", 
  help="Sets the C compiler to use. Default=%default")
parser.add_option("--fortran-compiler", dest="fortran", action="store", default="gfortran", 
  help="Sets the Fortran compiler to use. Default=%default")
parser.add_option("--enable-mpi", dest="enableMpi", action="store_true", default=False, 
  help="Enable an MPI build of Trilinos. Default=%default")
parser.add_option("--mpi-base-path", dest="mpiBaseDir", action="store", default="", 
  help="Sets the base path where we can find the mpi compiler wrappers, run environmnet and headers libraries. Default=%default")
parser.add_option("--boost-dir", dest="boostDir", action="store", default="/usr/lib", 
  help="Sets the location of the boost library. Default=%default")
parser.add_option("--netcdf-dir", dest="netcdfDir", action="store", default="/usr/lib", 
  help="Sets the location of the netcdf library. Default=%default")
parser.add_option("--no-dashboard", dest="doDashboardBuild", action="store_false", default=True, 
  help="Disables the dashboard submission build. The default is to do the dashboard build")
parser.add_option("--install-dir", dest="installDir", action="store", default=workingDir+"/install", 
  help="Sets the location of where Trilinos will be installed to. Default=%default")
parser.add_option("--working-dir", dest="workingDir", action="store", default=workingDir, 
  help="Sets the directory that this script will run from. Default=%default")
parser.add_option("--branch", dest="branch", action="store", default="master", 
  help="Sets the branch that should be used to build the tarball and installation from. Default=%default")
parser.add_option("--shared", dest="shared", action="store_true", default=False, 
  help="Build shared libraries instead of static. Default=%default")
parser.add_option("--enable-secondary-stable", dest="secondaryStable", action="store_true", default=False, 
  help="Enable secondary stable code. Default=master")
parser.add_option("--verbose-errors", dest="verboseErrors", action="store_true", default=False, 
  help="Enable verbose error reporting. This will cause the output of the command that failed to be sent to stdout. Default=master")

(options, args) = parser.parse_args()

options.configurePath = "configure"
options.tarballBuildPath = "tarball_build"
#strip any relative bits out of the path
options.workingDir = os.path.abspath(options.workingDir)

#create the working dir if it doesn't exist
if not os.path.exists(options.workingDir):
  os.mkdir(options.workingDir)

#switch to working dir. This is essentially a noop if the working dir = current dir.
os.chdir(options.workingDir)

print "Running nightly_create_tarball.py"
print "Installing to " + options.installDir

#choosing which compilers to use, mpi vs standard
compilerPathString = ""
if options.enableMpi:
  compilerPathString = """-D TPL_ENABLE_MPI=ON """
  if options.mpiBaseDir != "":
    baseDirString = "-D MPI_BASE_DIR=" + options.mpiBaseDir + " "
    compilerPathString = compilerPathString + baseDirString
else:
  compilerPathString = """\
-D CMAKE_CXX_COMPILER:FILEPATH="""     + options.cxx + """ \
-D CMAKE_C_COMPILER:FILEPATH="""       + options.cc + """ \
-D CMAKE_Fortran_COMPILER:FILEPATH=""" + options.fortran

if options.shared:
  sharedBuildOption = "-DBUILD_SHARED_LIBS:BOOL=ON"
else:
  sharedBuildOption = ""

if options.secondaryStable:
  #have to manually disable the teuchos float and complex options since they are enabled by default when not in development mode.
  secondaryStableOption = "-D Trilinos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON"
else:
  secondaryStableOption = ""

#Base command for configuring
#  Trikota, STK, and Optika are disabled due to explicit dependencies which can't always be met.
# mesquite is disabled because it doesn't generate a MesquiteConfig.cmake file.
baseConfigureCmd = """rm -rf CMakeCache.txt CMakefiles;cmake """ \
+ compilerPathString + """ \
-D HAVE_GCC_ABI_DEMANGLE:BOOL=ON \
-D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE """ \
+ secondaryStableOption + """ \
-DTrilinos_ENABLE_Amesos=ON \
-DTrilinos_ENABLE_Amesos2=OFF \
-DTrilinos_ENABLE_Anasazi=ON \
-DTrilinos_ENABLE_AztecOO=ON \
-DTrilinos_ENABLE_Belos=ON \
-DTrilinos_ENABLE_CTrilinos=OFF \
-DTrilinos_ENABLE_Didasko=OFF \
-DTrilinos_ENABLE_Epetra=ON \
-DTrilinos_ENABLE_EpetraExt=ON \
-DTrilinos_ENABLE_FEI=ON \
-DTrilinos_ENABLE_ForTrilinos=OFF \
-DTrilinos_ENABLE_Galeri=ON \
-DTrilinos_ENABLE_GlobiPack=ON \
-DTrilinos_ENABLE_Ifpack=ON \
-DTrilinos_ENABLE_Ifpack2=ON \
-DTrilinos_ENABLE_Intrepid=ON \
-DTrilinos_ENABLE_Isorropia=ON \
-DTrilinos_ENABLE_Kokkos=ON \
-DTrilinos_ENABLE_Komplex=ON \
-DTrilinos_ENABLE_Mesquite=OFF \
-DTrilinos_ENABLE_ML=ON \
-DTrilinos_ENABLE_Moertel=ON \
-DTrilinos_ENABLE_MOOCHO=ON \
-DTrilinos_ENABLE_NOX=ON \
-DTrilinos_ENABLE_Optika=OFF \
-DTrilinos_ENABLE_OptiPack=ON \
-DTrilinos_ENABLE_Pamgen=ON \
-DTrilinos_ENABLE_Phalanx=ON \
-DTrilinos_ENABLE_Piro=OFF \
-DTrilinos_ENABLE_Pliris=ON \
-DTrilinos_ENABLE_PyTrilinos=OFF \
-DTrilinos_ENABLE_RBGen=OFF \
-DTrilinos_ENABLE_RTOp=ON \
-DTrilinos_ENABLE_Rythmos=ON \
-DTrilinos_ENABLE_Sacado=ON \
-DTrilinos_ENABLE_SEACAS=OFF \
-DTrilinos_ENABLE_Shards=ON \
-DTrilinos_ENABLE_STK=OFF \
-DTrilinos_ENABLE_Stokhos=OFF \
-DTrilinos_ENABLE_Stratimikos=ON \
-DTrilinos_ENABLE_Sundance=OFF \
-DTrilinos_ENABLE_Teko=OFF \
-DTrilinos_ENABLE_Teuchos=ON \
-DTrilinos_ENABLE_ThreadPool=ON \
-DTrilinos_ENABLE_Thyra=ON \
-DTrilinos_ENABLE_Tpetra=ON \
-DTrilinos_ENABLE_TriKota=OFF \
-DTrilinos_ENABLE_TrilinosCouplings=OFF \
-DTrilinos_ENABLE_Triutils=ON \
-DTrilinos_ENABLE_Zoltan=ON """ \
+ sharedBuildOption + """ \
-D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-D TPL_ENABLE_Boost=ON \
-D Boost_INCLUDE_DIRS=""" + options.boostDir + """ \
-D TPL_ENABLE_Netcdf=ON \
-D Netcdf_LIBRARY_DIRS=""" + options.netcdfDir + """/lib \
-D Netcdf_INCLUDE_DIRS=""" + options.netcdfDir + """/include \
-D CMAKE_INSTALL_PREFIX:PATH=""" + options.installDir + " "


#removing install directory first so that subsequent installs aren't polluted by old installs
print "attempting to remove the old install dir"
try:
  shutil.rmtree(options.installDir, True)
except:
  print "execption while removing " + options.installDir
  print sys.exc_info()[1]

print "done trying to remove old install dir"

gitBranchError = 0
#clone repo
if os.path.exists("Trilinos"):
  print "Repository already exists so updating instead of cloning."
  os.chdir("Trilinos")
  if options.branch != "master":
    gitBranchCmd = "git checkout " + options.branch
    gitBranchError = os.system(gitBranchCmd)
  gitPullCmd = "git pull"
  gitError = os.system(gitPullCmd)
  os.chdir(options.workingDir)
else:
  print "Repository doesn't exist so cloning."
  gitCloneCmd = "git clone " + options.repository + " Trilinos"
  gitError = os.system(gitCloneCmd)
  os.chdir("Trilinos")
  if options.branch != "master":
    gitBranchCmd = "git checkout --track origin/" + options.branch
    gitBranchError = os.system(gitBranchCmd)
  os.chdir(options.workingDir)
  

reportError(gitError, "Could not retrieve repository " + options.repository)
reportError(gitBranchError, "Could not checkout branch " + options.branch)

#configure and create tarball
print "Configuring to create a tarball."
print "config dir exists? " + str(os.path.exists(options.configurePath))
print "cwd is " + os.getcwd()
if not os.path.exists(options.configurePath):
  print "Creating configure path in dir " + options.workingDir
  os.mkdir(options.configurePath)

os.chdir(options.configurePath)

configureCmd = baseConfigureCmd + options.workingDir + "/Trilinos &> configure.out"
configureError = os.system(configureCmd)

configureErrorMessage = "Failure while configuring for tarball creation see configure.out for details." + os.linesep + configureCmd
reportError(configureError, configureErrorMessage, "configure.out", options.verboseErrors)
#if configureError != 0:
#  print "Failure while configuring for tarball creation see configure.out for details."
#  print configureCmd
#  sys.exit(1)

print "Creating a tarball."
makeTarballCmd = "make package_source &> make_tarball.out"
makeTarballError = os.system(makeTarballCmd)

reportError(makeTarballError, "Failure while creating tarball, see make_tarball.out for details.", "make_tarball.out", options.verboseErrors)

#expand tarball
os.chdir(options.workingDir)

if os.path.exists(options.tarballBuildPath):
  print "Deleting directory " + options.workingDir + "/" + options.tarballBuildPath
  shutil.rmtree(options.tarballBuildPath)

os.mkdir(options.tarballBuildPath)
os.chdir(options.tarballBuildPath)
fileList = os.listdir(options.workingDir + "/" + options.configurePath)

#find the tarball we need to copy
for entry in fileList:
  if re.match("trilinos-.*-Source\.tar\.gz", entry):
    tarballFileName = entry
    print "found tarball: '" + tarballFileName + "'"
    break

print "Expanding the tarball."
shutil.copy(options.workingDir + "/" + options.configurePath + "/" + tarballFileName, ".")

untarCmd = "tar -xzvf " + tarballFileName + " &> tar.out"
untarError = os.system(untarCmd)

reportError(untarError, "Could not expand tarball, see tar.out for details", "tar.out", options.verboseErrors)

#Make install
print "Building and installing from the tarball."
trilinosDir = tarballFileName[:-7]
os.mkdir("build")
os.chdir("build")

tarballConfigureCmd = baseConfigureCmd + "-D Trilinos_ENABLE_DEVELOPMENT_MODE:BOOL=OFF " + \
  options.workingDir + "/" + options.tarballBuildPath + "/" + trilinosDir + " &> tarball_configure.out"
tarballConfigureError = os.system(tarballConfigureCmd)

reportError(tarballConfigureError, "Configuring from tarball failed, see tarball_configure.out for details.",
            "tarball_configure.out", options.verboseErrors)

tarballMakeInstallCmd = "make install &> tarball_make_install.out"
tarballMakeInstallError = os.system(tarballMakeInstallCmd)

reportError(tarballMakeInstallError, "Make install from the tarball failed, see tarball_make_install.out for details.",
            "tarball_make_install.out", options.verboseErrors)

if options.doDashboardBuild:
  makeExperimentalBuildCmd = "make dashboard &> tarball_make_experimental.out"
  makeExperimentalBuildError = os.system(makeExperimentalBuildCmd)

  reportError(makeExperimentalBuildError, "Make experimental from the tarball failed, see tarball_make_experimental.out for details.",
              "tarball_make_experimental.out", options.verboseErrors)

print "Installation from a tarball completed successfully."
