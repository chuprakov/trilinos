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
