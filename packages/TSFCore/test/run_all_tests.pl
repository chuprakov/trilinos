#!/usr/bin/perl -w
#
# This perl scipt should be run from the directory:
#
#  Trilinos/$BUILD_DIR/package/TSFCore/test/.
#
#

use strict;

my $ierr = 0;
my $iresult;

my $cmnd;

$cmnd = "./epetra-adapters/epetra_adapters.exe --max-flop-rate=1";
print "\n${cmnd}\n";
$iresult = system($cmnd);
$ierr += $iresult;
print "\n***\n*** Above testing program passed!\n***\n" if($iresult==0);
print "\n***\n*** Above testing program failed!\n***\n" if($iresult!=0);

$cmnd = "./product-space/product_space.exe";
print "\n${cmnd}\n";
$iresult = system($cmnd);
$ierr += $iresult;
print "\n***\n*** Above testing program passed!\n***\n" if($iresult==0);
print "\n***\n*** Above testing program failed!\n***\n" if($iresult!=0);

$cmnd = "./scalar-product/scalar_product.exe";
print "\n${cmnd}\n";
$iresult = system($cmnd);
$ierr += $iresult;
print "\n***\n*** Above testing program passed!\n***\n" if($iresult==0);
print "\n***\n*** Above testing program failed!\n***\n" if($iresult!=0);

$cmnd = "./std-ops/std_ops.exe";
print "\n${cmnd}\n";
$iresult = system($cmnd);
$ierr += $iresult;
print "\n***\n*** Above testing program passed!\n***\n" if($iresult==0);
print "\n***\n*** Above testing program failed!\n***\n" if($iresult!=0);

$cmnd = "../example/Core/sillyCgSolve_epetra.exe";
print "\n${cmnd}\n";
$iresult = system($cmnd);
$ierr += $iresult;
print "\n***\n*** Above testing program passed!\n***\n" if($iresult==0);
print "\n***\n*** Above testing program failed!\n***\n" if($iresult!=0);

$cmnd = "../example/Core/sillyCgSolve_serial.exe";
print "\n${cmnd}\n";
$iresult = system($cmnd);
$ierr += $iresult;
print "\n***\n*** Above testing program passed!\n***\n" if($iresult==0);
print "\n***\n*** Above testing program failed!\n***\n" if($iresult!=0);

$cmnd = "../example/Core/sillyPowerMethod_epetra.exe";
print "\n${cmnd}\n";
$iresult = system($cmnd);
$ierr += $iresult;
print "\n***\n*** Above testing program passed!\n***\n" if($iresult==0);
print "\n***\n*** Above testing program failed!\n***\n" if($iresult!=0);

$cmnd = "../example/Core/sillyPowerMethod_serial.exe";
print "\n${cmnd}\n";
$iresult = system($cmnd);
$ierr += $iresult;
print "\n***\n*** Above testing program passed!\n***\n" if($iresult==0);
print "\n***\n*** Above testing program failed!\n***\n" if($iresult!=0);

$cmnd = "../example/Nonlin/NP2DSim/TSFCoreNonlinNP2DSimTest.exe";
print "\n${cmnd}\n";
$iresult = system($cmnd);
$ierr += $iresult;
print "\n***\n*** Above testing program passed!\n***\n" if($iresult==0);
print "\n***\n*** Above testing program failed!\n***\n" if($iresult!=0);

$cmnd = "../example/Nonlin/NP4DOpt/TSFCoreNonlinNP4DOptTest.exe";
print "\n${cmnd}\n";
$iresult = system($cmnd);
$ierr += $iresult;
print "\n***\n*** Above testing program passed!\n***\n" if($iresult==0);
print "\n***\n*** Above testing program failed!\n***\n" if($iresult!=0);

if($ierr) {
  print "\n***\n*** Oh no, at least one of the TSFCore testing programs failed!\n***\n";
}
else {
  print "\n***\n*** Congratulations, All tests in TSFCore seemed to have passed!\n***\n";
}

exit($ierr);
