#!/usr/bin/perl -w
#
# This perl scipt should be run from the directory:
#
#  Trilinos/$BUILD_DIR/package/TSFCore/test/.
#
#

use strict;

my $ierr = 0;

my $cmnd;

$cmnd = "./epetra-adapters/epetra_adapters.exe";
$ierr += system($cmnd);

$cmnd = "./product-space/product_space.exe";
$ierr += system($cmnd);

$cmnd = "./scalar-product/scalar_product.exe";
$ierr += system($cmnd);

$cmnd = "./std-ops/std_ops.exe";
$ierr += system($cmnd);

$cmnd = "../example/Core/sillyCgSolve_epetra.exe";
$ierr += system($cmnd);

$cmnd = "../example/Core/sillyCgSolve_serial.exe";
$ierr += system($cmnd);

$cmnd = "../example/Core/sillyPowerMethod_epetra.exe";
$ierr += system($cmnd);

$cmnd = "../example/Core/sillyPowerMethod_serial.exe";
$ierr += system($cmnd);

$cmnd = "../example/Nonlin/NP2DSim/TSFcoreNonlinNP2DSimTest.exe";
$ierr += system($cmnd);

$cmnd = "../example/Nonlin/NP4DOpt/TSFcoreNonlinNP4DOptTest.exe";
$ierr += system($cmnd);

if($ierr) {
	print "\n***\n*** Oh no, at least one of the TSFCore testing programs failed!\n***\n";
}
else {
	print "\n***\n*** Congratulations, All tests in TSFCore seemed to have passed!\n***\n";
}

exit($ierr);
