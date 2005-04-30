#!/usr/bin/perl -w

use strict;

foreach( <STDIN> ) {
	my ($old_file_name, $new_file_name) = split;
	if( $old_file_name ne $new_file_name ) {
		print "\"${old_file_name}\"\n\"${new_file_name}\"\n\n";
	}
}
