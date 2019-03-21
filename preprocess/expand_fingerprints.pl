#!/usr/bin/perl

# Copyright (c) 2018 Los Alamos National Security, LLC.
# All rights reserved.
# 
# Copyright 2018. Los Alamos National Security, LLC. This software was produced under U.S. Government 
# contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos 
# National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, 
# reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC 
# MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If 
# software is modified to produce derivative works, such modified software should be clearly marked, so 
# as not to confuse it with the version available from LANL.
# 
# This work has been supported in part by the Joint Design of Advanced Computing Solutions for Cancer 
# (JDACS4C) program established by the U.S. Department of Energy (DOE) and the National Cancer Institute 
# (NCI) of the National Institutes of Health. 
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
# associated documentation files (the "Software"), to deal in the Software without restriction, including 
# without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to 
# the following conditions: 
# 
# The above copyright notice and this permission notice shall be included in all copies or substantial 
# portions of the Software.
# 
# THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Author: Jason D. Gans (jgans@lanl.gov)

# Expand the OpenBabel hex-encoded fingerprints into tab-delimited binary
# strings

use strict;

my %hex_to_dec = (
	'0' => 0,
	'1' => 1,
	'2' => 2,
	'3' => 3,
	'4' => 4,
	'5' => 5,
	'6' => 6,
	'7' => 7,
	'8' => 8,
	'9' => 9,
	'a' => 10,
	'b' => 11,
	'c' => 12,
	'd' => 13,
	'e' => 14,
	'f' => 15
);

sub main()
{
	my $delimiter = ',';
	
	# Print the header
	print "DRUG" . $delimiter . "FINGERPRINT\n";
	
	LINE: while(my $line = <STDIN>){
		
		chomp($line);
		
		if($line =~ /^#/){
			
			print STDERR "Skipping: $line\n";
			next LINE;
		}
		
		my @data = split /\t/, $line;
		
		if(scalar(@data) != 2){
			die "Did not read the expected two column format!\n";
		}
		
		my $id = $data[1];
		
		# Remove any prefix from the id string
		if($id =~ /[A-Z|a-z]+\.(\d+)/){
			$id = $1;
		}
		
		# Convert the fingerprint to lower case
		$data[0] =~ tr/A-Z/a-z/;
		
		my @fingerprint = split //, $data[0];
		
		print "$id$delimiter";
		
		# Currently, the fingerprint bits are concatinated into a single string
		foreach my $f (@fingerprint){
			
			if( !defined($hex_to_dec{$f}) ){
				die "Unknown symbol in fingerprint ($f)\n";
			}

			# Convert the hex digit to its decimal representation
			my $value = $hex_to_dec{$f};

			# Each hex value expands to four binary values
			print (($value & 1) ? 1 : 0);
			print (($value & 2) ? 1 : 0);
			print (($value & 4) ? 1 : 0);
			print (($value & 8) ? 1 : 0);
		}
		
		print "\n";
	}
}

main();

exit;
