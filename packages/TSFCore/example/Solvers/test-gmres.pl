# ************************************************************************
# 
#               TSFCore: Trilinos Solver Framework Core
#                 Copyright (2004) Sandia Corporation
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************

#!/usr/bin/perl -w
use strict;
use strict 'refs';
use File::Basename;
my $success = 1;
my $result;
#
# Parse options
#
#my $though_tests=0;
#my $verbose=1;
#
#my $num_args = scalar(@ARGV);
#my ( $script_name, $script_path, $script_ext ) = fileparse($0);
#for (my $i = 0; $i < $num_args; ++$i) {
#  $_ = $ARGV[$i];
#  if (/^-h$/) {
#    print STDERR
#      "Use : ${script_name} [--quiet, --verbose] [\n",
#        $g_use_msg_opts;
#    exit(0);
#  } elsif (/^-t/) {
#    @xml_module_files = ( @xml_module_files, $ARGV[++$i] );
#  } else {
#    @argv_trimmed = (@argv_trimmed, $_);
#  }
#}
# Run tests
$result = system ('./TSFCoreSolversGmresTest.exe ../../../../../packages/belos/test/BlockGmres/orsirr1.hb 500 1e-10 -v');
$success = 0 if ($result!=0);
# RAB: 2003/11/21: Commented out because they are very expensive (todo: ad flag for including them back)
#$result = system ('./TSFCoreSolversGmresTest.exe ../../../../../packages/belos/test/BlockGmres/fidap036.hb 1500 1e-10 -v');
#$success = 0 if ($result!=0);
#$result = system ('./TSFCoreSolversGmresTest.exe ../../../../../packages/belos/test/BlockCG/bcsstk14.hb 1200 1e-10 -v');
#$success = 0 if ($result!=0);
#$result = system ('./TSFCoreSolversGmresTest.exe ../../../../../packages/belos/test/BlockCG/bcsstk15.hb 2100 1e-10 -v');
#$success = 0 if ($result!=0);
#
if($success) {
  printf "\nFinal: Congratulations! all of the linear systems where solved to the given tolerances!\n";
}
else {
  printf "\nFinal: Oh no! at least one of the linear systems was not solved to the given tolerances!\n";
}
exit ($success ? 0 : -1 );
