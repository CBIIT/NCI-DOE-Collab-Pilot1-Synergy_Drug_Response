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

# Convert the Merck drug pair data to same CSV file format used for the
# ALMANAC drugs:
# DRUG1, DRUG2, CELL_LINE1, CELL_LINE2, CELL_LINE3, ...
#
# 1) Map Merck drug names to their corresponding drug CID values
# 1.5) Interpolate growth values using the log(concentration).
# 2) Compute synergy using the Bliss independence model (for now, the
#    Loewe synergy definition will have to wait).
# 3) Define synergy as: (observed growth) - (expected growth)
# 4) Return the minimum synergy for each drug combination and cell line

# For now, use the raw "X/X0" values reported in the Merck data as the growth. According
# to the variable definitions, X/X0 == 0 is "no growth", while X/X0 == 1 is "same growth
# as the DMSO control". There are X/X0 values that are greater than 1, and it is not
# clear if these should be clamped to be <= 1 (values are not currently clamped).

# According to the Merck study ("An Unbiased Oncology Compound Screen to Identify Novel Combination Strategies", O'Niel et al.,
# Mol Cancer Ther. 2016 Jun;15(6):1155-62. doi: 10.1158/1535-7163.MCT-15-0843. Epub 2016 Mar 16.)
# "The screen was conducted in 3 batches in terms of cell line growth and drug batches. Batch 1 & 2 are part of the primary
# large screen while batch 3 is a partial repeat of selected drugs and served as a validation screen".
# For only, only include batch 1 and batch 2.

# Update Aug 4, 2018
#	* Updated Merck drug name to CID mapping
#	* Added optional Merck cell line to CCLE mapping
# Update Sept 17, 2018
#	* Added cell lines from the most recent version of the CCLE RNA-seq dataset. These additional
#	  cell lines are: DLD1 (LARGE INTESTINE) and PA1 (OVARY).
#	* Added additional CCLE cell lines that are expected to be closely related (by not exact
#	  matches) to the Merck cell lines:
#		COLO320DM (Merck) --> COLO320 (CCLE)
#		EFM192B (Merck) --> EFM192A (CCLE)
#
# Updated March 1, 2019
#   * Modified file parser to accept both CSV and UTF-8 CSV files exported from Excel
#     without the need to manually edit the CSV files.
#   * Drug and cell names in the output file are now in sorted in ascending order

use Cwd;
use Getopt::Long;
use strict;

# Name to CID mapping obtained by manually searching the PubChem database
my %drug_to_cid = (
	"5-FU" => 3385,
	"ABT-888" => 11960529,
	"AZD1775" => 24856436,
	"BEZ-235" => 11977753,
	"Bortezomib" => 387447,
	"Carboplatin" => 426756,
	"Cyclophosphamide" => 2907,
	"Dasatinib" => 3062316,
	"Dexamethasone" => 5743,
	"Dinaciclib" => 46926350,
	"Doxorubicin" => 32875,
	"Erlotinib" => 176870,
	"Etoposide" => 439525,
	"geldanamycin" => 5288382,
	"Gemcitabine" => 3461,
	"L778123" => 216454,
	"Lapatinib" => 208908,
	"Metformin" => 4091,
	"Methotrexate" => 126941,
	"Mitomycine" => 5746,
	"MK-2206" => 24964624,
	"MK-4541" => 59691338,
	"MK-4827" => 24958200,
	"MK-5108" => 24748204,
	"MK-8669" => 11520894,
	"MK-8776" => 16224745,
	"MRK-003" => 15953832,
	"Oxaliplatin" => 24197464,
	"Paclitaxel" => 36314,
	"PD325901" => 9826528,
	"SN-38" => 104842,
	"Sorafenib" => 216239,
	"Sunitinib" => 5329102,
	"Temozolomide" => 5394,
	"Topotecan" => 60700,
	"Vinblastine" => 13342,
	"Vinorelbine" => 5311497,
	"Zolinza" => 5311
);

# The mapping between Merck cell line names and CCLE cell line names
my %cell_to_CCLE = (
	"A2058" => "CCLE.A2058",
	"A2780" => "CCLE.A2780",
	"A375" => "CCLE.A375",
	"A427" => "CCLE.A427",
	"CAOV3" => "CCLE.CAOV3",
	"COLO320DM" => "CCLE.COLO320", # Related, but not exact match
	"DLD1" => "CCLE.DLD1",
	"EFM192B" => "CCLE.EFM192A", # Related, but not exact match
	"ES2" => "CCLE.ES2",
	"HCT116" => "CCLE.HCT116",
	"HT144" => "CCLE.HT144",
	"HT29" => "CCLE.HT29",
	"KPL1" => "CCLE.KPL1",
	"LNCAP" => "CCLE.LNCAPCLONEFGC",
	"LOVO" => "CCLE.LOVO",
	"MDAMB436" => "CCLE.MDAMB436",
	"MSTO" => "CCLE.MSTO211H",
	"NCIH1650" => "CCLE.NCIH1650",
	"NCIH2122" => "CCLE.NCIH2122",
	"NCIH23" => "CCLE.NCIH23",
	"NCIH460" => "CCLE.NCIH460",
	"NCIH520" => "CCLE.NCIH520",
	"OV90" => "CCLE.OV90",
	"OVCAR3" => "CCLE.NIHOVCAR3",
	"PA1" => "CCLE.PA1",
	"RKO" => "CCLE.RKO",
	"RPMI7951" => "CCLE.RPMI7951",
	"SKMEL30" => "CCLE.SKMEL30",
	"SKMES1" => "CCLE.SKMES1",
	"SKOV3" => "CCLE.SKOV3",
	"SW620" => "CCLE.SW620",
	"SW837" => "CCLE.SW837",
	"T47D" => "CCLE.T47D",
	"UACC62" => "CCLE.UACC62",
	"VCAP" => "CCLE.VCAP",
	"ZR751" => "CCLE.ZR751"
);

# The drugs for which we have all-against-all data
my %complete_drugs = (
	"ABT-888" => 1,
	"AZD1775" => 1,
	"BEZ-235" => 1,
	"Bortezomib" => 1,
	"Dasatinib" => 1,
	"Dinaciclib" => 1,
	"Erlotinib" => 1,
	"geldanamycin" => 1,
	"L778123" => 1,
	"Lapatinib" => 1,
	"MK-2206" => 1,
	"MK-4541" => 1,
	"MK-4827" => 1,
	"MK-5108" => 1,
	"MK-8669" => 1,
	"MK-8776" => 1,
	"MRK-003" => 1,
	"PD325901" => 1,
	"Sorafenib" => 1,
	"Sunitinib" => 1,
	"Temozolomide" => 1,
	"Zolinza" => 1
);

Getopt::Long::Configure ("bundling");

$|                      = '1';

my $opt_single = "";
my $opt_pair = "";
my $opt_help = (scalar(@ARGV) == 0);

my $command_line = $0 . " " . join(' ', @ARGV);

GetOptions(
   'single=s'  		=> \$opt_single,
   'pair=s'  		=> \$opt_pair,
   'h|?|help'   	=> \$opt_help
);

my $usage = "USAGE: $0\n" .
	"\t--single <Merck single agent CSV data file>\n" .
	"\t--pair <Merck combination response CSV data file>\n" .
	"\t[-h|-?|--help] (print usage)\n";

sub main()
{
	if($opt_help){
   		die $usage;
	}
	
	if($opt_single eq ""){
		die "Please specify a file of single agent responses (--single)\n";
	}
	
	if($opt_pair eq ""){
		die "Please specify a file of combination responses (--pair)\n";
	}
	
	my $single_drug_filename = $opt_single;
	my $pair_drug_filename = $opt_pair;
	
	my $use_gi = 0;
	
	# Set use_viability to 1 to enable viability as the independent variable. Otherwise
	# "X/X0" will be used.
	my $use_viability = 1;
	
	# Set use_ccle to convert Merck cell line names to CCLE cell line names
	my $use_ccle = 1;
	
	# Set to 1 to use abs(viability), otherwise some of the pairwise viability values that are
	# slightly below 0.0 will be kept as negative numbers.
	my $clamp_viability = 0;
	
	if($use_viability && $use_gi){
		die "Please select 'use_viability' or 'use_gi', but not both!\n";
	}

	if($use_viability && $clamp_viability){
		print STDERR "Viability values are clamped to be greater than zero\n";
	}
	
	# The duration of the experiment in arbitrary time units
	my $duration = 96.0;
	my $EPSILON = 1.0e-2;
	
	if($use_viability){
		printf STDERR "Using viability\n";
	}
	else{
		printf STDERR "Using X/X0\n";
	}
	
	if($use_ccle){
		printf STDERR "Using CCLE cell lines\n";
	}
	else{
		printf STDERR "Using Merck cell lines\n";
	}
	
	##########################################################################
	# Parse the single drug data
	##########################################################################
	
	print STDERR "Parsing the single drug data\n";
	
	open(my $INPUT, "$single_drug_filename") or
		die "Unable to open $single_drug_filename for reading single agent data\n";
	
	my $line = "";
	
	if( !($line = <$INPUT>) ){
		die "Error reading header from $INPUT\n";
	}
	
	# If the user has exported the excel file as UTF-8, then there will be a single
    	# special character at the start of the file (a "byte mark"). We will attempt
    	# to remove this character (as a precationary step) by removing *all*
    	# non-ASCII characters. Also remove the DOS carriage return
    	$line =~ tr/\x00-\x7F//cd;
    	$line =~ tr/\r\n//d;
	
	my @header = split /,/, $line;
	
	# Extract the column ids
	my $batch_col = get_column_index(\@header, "BatchID");
	my $drug_col = get_column_index(\@header, "drug_name");
	my $cell_line_col = get_column_index(\@header, "cell_line");
    
    	# Note that the wild-card in front of 'M' is to handle the ASCII 'u'
    	# or (UTF-8) greek mu symbol (which depends on the formate of the input file)
	my $conc_col = get_column_index(\@header, 'Drug_concentration \(.*M\)');
	my $viability1_col = get_column_index(\@header, "viability1");
	my $viability2_col = get_column_index(\@header, "viability2");
	my $viability3_col = get_column_index(\@header, "viability3");
	my $viability4_col = get_column_index(\@header, "viability4");
	my $viability5_col = get_column_index(\@header, "viability5");
	my $viability6_col = get_column_index(\@header, "viability6");
	my $mu_muMax_col = get_column_index(\@header, "mu/muMax");
	my $growth_col = get_column_index(\@header, "X/X0");
	
	# single_data: drug -> cell line -> conc -> growth
	# Average over replicates
	my %single_data = ();
	my %norm_data = ();
	
	my $line_number = 1;
	
	while($line = <$INPUT>){
		
		++$line_number;
		
		# Like chomp($line), but remove any DOS \r characters as well.
		$line =~ tr/\n\r//d;
		
		my @data = split /,/, $line;
		
		my $batch = $data[$batch_col];
		my $drug = $data[$drug_col];
		my $cell_line = $data[$cell_line_col];
		my $conc = $data[$conc_col];
		
		my @viability = ($data[$viability1_col],
				 $data[$viability2_col],
				 $data[$viability3_col],
				 $data[$viability4_col],
				 $data[$viability5_col],
				 $data[$viability6_col]);
		
		my $mu_muMax = $data[$mu_muMax_col];
		 
		my $growth = $data[$growth_col];
		
		if( !defined($batch) || ($batch eq "") ){
			die "Unable to read batch from $line\n";
		}
		
		if( !defined($drug) || ($drug eq "") ){
			die "Unable to read drug from $line\n";
		}
		
		if( !defined($cell_line) || ($cell_line eq "") ){
			die "Unable to read cell line from $line\n";
		}
		
		if( !defined($conc) || ($conc eq "") ){
			die "Unable to read conc from $line\n";
		}
		
		foreach my $v (@viability){
			if( !defined($v) || ($v eq "") ){
				die "Unable to read viability from $line\n";
			}
		}
		
		if( !defined($mu_muMax) || ($mu_muMax eq "") ){
			die "Unable to read mu/muMax from $line\n";
		}		
		
		if( !defined($growth) || ($growth eq "") ){
			die "Unable to read growth from $line\n";
		}
		
		# Skip batch 3 for now
		if( ($batch == 1) || ($batch == 2) ){
			
			my $y = 0.0;
			
			if($use_viability){
			
				$y = median(\@viability);
				
				# Should we allow *negative* viability values?
				if( $clamp_viability && ($y < 0.0) ){
					$y *= -1.0;
				}
			}
			else{
				$y = $growth;				
			}
			
			# Clamp the fractional growth at 1.0
			if($y > 1.0){
				$y = 1.0;
			}
			
			if($use_gi){
			
				########## NCI60 %GI ##########
				my $gi = 0.0;
				
				my $muMax = 1.0;
				
				if( abs($mu_muMax - 1.0) > $EPSILON){ # Epsilon is a small number
					$muMax = log( median(\@viability) )/( $duration*($mu_muMax - 1.0) );
				}
				
				my $mu = $mu_muMax*$muMax;

				if($mu >= 0.0){
					$gi = 100.0*(exp($mu*$duration) - 1.0)/(exp($muMax*$duration) - 1.0);
				}
				else{ # mu < 0.0
					$gi = 100.0*(exp($mu*$duration) - 1.0);
				}
				
				# Clamp the value of %GI to <= 100.0
				if($gi > 100.0){
					$gi = 100.0;
				}
				
				$y = $gi;
			}
			
			if( defined($single_data{$drug}{$cell_line}{$conc}) ){

				$single_data{$drug}{$cell_line}{$conc} += $y;
				++$norm_data{$drug}{$cell_line}{$conc};
			}
			else{
				$single_data{$drug}{$cell_line}{$conc} = $y;
				$norm_data{$drug}{$cell_line}{$conc} = 1;
			}
		}
	}
	
	close($INPUT);
	
	# Extract all of the drug names
	my @valid_drugs = keys(%single_data);
	
	print STDERR "Found " . scalar(@valid_drugs) . " drugs\n";
	
	##########################################################################
	# Normalize the single drug data
	##########################################################################
	foreach my $drug (@valid_drugs){
		
		my @valid_cells = keys(%{$single_data{$drug}});
		
		foreach my $cell (@valid_cells){
			
			# Sort the concentrations to make debugging easier
			my @valid_conc = sort {$a <=> $b} keys(%{$single_data{$drug}{$cell}});
			
			foreach my $conc (@valid_conc){
				$single_data{$drug}{$cell}{$conc} /= $norm_data{$drug}{$cell}{$conc};
				
				# DEBUG
				#print STDERR "$conc\t$single_data{$drug}{$cell}{$conc}\n";
			}
					
			#print STDERR "\t\n";
		}
	}
	
	##########################################################################
	# Parse the drug pair data
	##########################################################################
	
	print STDERR "Parsing the paired drug data\n";
	
	open($INPUT, "$pair_drug_filename") or
		die "Unable to open $pair_drug_filename for reading paired agent data\n";
	
	if( !($line = <$INPUT>) ){
		die "Error reading header from $INPUT\n";
	}
	
	# If the user has exported the excel file as UTF-8, then there will be a single
	# special character at the start of the file (a "byte mark"). We will attempt
	# to remove this character (as a precationary step) by removing *all* 
	# non-ASCII characters. Also remove the DOS carriage return
	$line =~ tr/\x00-\x7F//cd;
	$line =~ tr/\r\n//d;

	@header = split /,/, $line;
	
	# Extract the column ids
	$batch_col = get_column_index(\@header, "BatchID");
	my $drug_A_col = get_column_index(\@header, "drugA_name");
	my $drug_B_col = get_column_index(\@header, "drugB_name");
	$cell_line_col = get_column_index(\@header, "cell_line");
    
    	# Note that the wild-card in front of 'M' is to handle the ASCII 'u'
    	# or (UTF-8) greek mu symbol (which depends on the formate of the input file)
	my $conc_A_col = get_column_index(\@header, 'drugA Conc \(.*M\)');
	my $conc_B_col = get_column_index(\@header, 'drugB Conc \(.*M\)');
    
	$viability1_col = get_column_index(\@header, "viability1");
	$viability2_col = get_column_index(\@header, "viability2");
	$viability3_col = get_column_index(\@header, "viability3");
	$viability4_col = get_column_index(\@header, "viability4");
	$mu_muMax_col = get_column_index(\@header, "mu/muMax");
	$growth_col = get_column_index(\@header, "X/X0");
	
	# pair_data: drug A -> drug B -> cell line -> conc A -> conc B -> growth
	# Note that drug A and drug B will be ordered such that (drug A) < (drug B)
	# using Perl's lexographic "lt" operator.
	# Average over replicates
	my %pair_data = ();
	%norm_data = (); # <-- clean the norm_data hash prior to reuse
	
	$line_number = 1;
	
	while($line = <$INPUT>){
		
		++$line_number;
		
		# Like chomp($line), but remove any DOS \r as well ...
		$line =~ tr/\n\r//d;
		
		my @data = split /,/, $line;
		
		my $batch = $data[$batch_col];
		my $drug_A = $data[$drug_A_col];
		my $drug_B = $data[$drug_B_col];
		my $cell_line = $data[$cell_line_col];
		my $conc_A = $data[$conc_A_col];
		my $conc_B = $data[$conc_B_col];
		
		my @viability = ($data[$viability1_col],
				 $data[$viability2_col],
				 $data[$viability3_col],
				 $data[$viability4_col]);
		
		my $mu_muMax = $data[$mu_muMax_col];
			 
		my $growth = $data[$growth_col];
		
		if( !defined($batch) || ($batch eq "") ){
			die "Unable to read batch from $line\n";
		}
		
		if( !defined($drug_A) || ($drug_A eq "") ){
			die "Unable to read drug A from $line\n";
		}
		
		if( !defined($drug_B) || ($drug_B eq "") ){
			die "Unable to read drug B from $line\n";
		}
		
		if( !defined($cell_line) || ($cell_line eq "") ){
			die "Unable to read cell line from $line\n";
		}
		
		if( !defined($conc_A) || ($conc_A eq "") ){
			die "Unable to read conc A from $line\n";
		}
		
		if( !defined($conc_B) || ($conc_B eq "") ){
			die "Unable to read conc B from $line\n";
		}
		
		foreach my $v (@viability){
			if( !defined($v) || ($v eq "") ){
				die "Unable to read viability from $line\n";
			}
		}
		
		if( !defined($mu_muMax) || ($mu_muMax eq "") ){
			die "Unable to read mu/muMax from $line\n";
		}
		
		if( !defined($growth) || ($growth eq "") ){
			die "Unable to read growth from $line\n";
		}
		
		if($drug_A eq $drug_B){
			die "Found record with the same drugs in A and B positions at line $line_number\n";
		}

		# Swap the drugs if needed
		if( ($drug_A gt $drug_B) ){
			
			# Use the perl list context to perform a fancy, "perlish" swap
			($drug_A, $drug_B) = ($drug_B, $drug_A);
			($conc_A, $conc_B) = ($conc_B, $conc_A);
		}
		
		# Skip batch 3 for now
		if( ($batch == 1) || ($batch == 2) ){
			
			my $y = 0.0;
			
			if($use_viability){
			
				$y = median(\@viability);
				
				# Should we allow *negative* viability values?
				if( $clamp_viability && ($y < 0.0) ){
					$y *= -1.0;
				}
			}
			else{
				$y = $growth;
			}
			
			# Clamp the fractional growth at 1.0
			if($y > 1.0){
				$y = 1.0;
			}
			
			if($use_gi){
			
				########## NCI60 %GI ##########
				my $gi = 0.0;
				
				my $muMax = 1.0;
				my $mu = 1.0;
				
				if( abs($mu_muMax - 1.0) > $EPSILON){ # Epsilon is a small number
					
					my $m = median(\@viability);
					
					if($m > $EPSILON){
					
						$muMax = log($m)/( $duration*($mu_muMax - 1.0) );
						$mu = $mu_muMax*$muMax;
					}
					else{
						#print STDERR "median ($m) is <= 0\n";
						
						# For now, set mu to be a large, negative number
						$mu = -10.0*abs($muMax);
					}
				}
				
				#print STDERR "$muMax\t" . ($mu) . "\n";
				
				if($mu >= 0.0){
					$gi = 100.0*(exp($mu*$duration) - 1.0)/(exp($muMax*$duration) - 1.0);
				}
				else{ # mu < 0.0
					$gi = 100.0*(exp($mu*$duration) - 1.0);
				}
				
				# Clamp the value of %GI to <= 100.0
				if($gi > 100.0){
					$gi = 100.0;
				}
				
				#print STDERR "mu_muMax = $mu_muMax\n";
				#print STDERR "mu_muMax - 1 = " . ($mu_muMax - 1) . "\n";
				
				#print STDERR "duration = $duration\n";
				#print STDERR "median(viability) = " . median(\@viability) . "\n";
				#print STDERR "%GI = $gi\n\n";

				$y = $gi;
			}
			
			if( defined($pair_data{$drug_A}{$drug_B}{$cell_line}{$conc_A}{$conc_B}) ){

				$pair_data{$drug_A}{$drug_B}{$cell_line}{$conc_A}{$conc_B} += $y;
				++$norm_data{$drug_A}{$drug_B}{$cell_line}{$conc_A}{$conc_B};
			}
			else{
				$pair_data{$drug_A}{$drug_B}{$cell_line}{$conc_A}{$conc_B} = $y;
				$norm_data{$drug_A}{$drug_B}{$cell_line}{$conc_A}{$conc_B} = 1;
			}
		}
	}
	
	close($INPUT);
	
	# Extract all of the drug A names
	my @valid_drug_A = sort keys(%pair_data);
	
	##########################################################################
	# Normalize the pair drug data
	##########################################################################
	
	# Collect all of the cell line names (that appear in any experiment)
	my %all_cells = ();
	
	print STDERR "Normalizing the paired drug data\n";
	
	foreach my $drug_A (@valid_drug_A){
		
		my @valid_drug_B = keys(%{$pair_data{$drug_A}});
		
		foreach my $drug_B (@valid_drug_B){
		
			my @valid_cells = keys(%{$pair_data{$drug_A}{$drug_B}});

			foreach my $cell (@valid_cells){

				$all_cells{$cell} = 1;
				
				my @valid_conc_A = keys(%{$pair_data{$drug_A}{$drug_B}{$cell}});

				foreach my $conc_A (@valid_conc_A){
					
					my @valid_conc_B = keys(%{$pair_data{$drug_A}{$drug_B}{$cell}{$conc_A}});
					
					foreach my $conc_B (@valid_conc_B){
						$pair_data{$drug_A}{$drug_B}{$cell}{$conc_A}{$conc_B} /= 
							$norm_data{$drug_A}{$drug_B}{$cell}{$conc_A}{$conc_B};
					}
				}
			}
		}
	}
	
	my @cell_line_names = sort keys(%all_cells);
	
	print STDERR "Found " . scalar(@cell_line_names) . " Merck cell line names\n";
	
	if($use_ccle){
		
		print STDERR "Switching to CCLE cell line names\n";
		
		my @tmp = ();
		
		foreach my $c (@cell_line_names){
			
			if( defined($cell_to_CCLE{$c}) ){
				push @tmp, $c;
			}
		}
		
		@cell_line_names = @tmp;
		
		print STDERR "Found " . scalar(@cell_line_names) . 
			" Merck cell lines in the CCLE cell lines\n";
	}
	
	##########################################################################
	# Write the drug pair synergy data
	##########################################################################
	
	# The output format is CSV:
	# DRUG1,DRUG2,CELL_LINE_0,CELL_LINE_1,CELL_LINE_2,...
	
	print STDERR "Writing the paired drug synergy data\n";
	
	# Write the output header
	print STDOUT "DRUG1,DRUG2";
	
	foreach my $c (@cell_line_names){
		
		if($use_ccle){
			
			if( !defined($cell_to_CCLE{$c}) ){
				die "Did not find matching CCLE cell line name for $c\n";
			}
			
			print STDOUT ",$cell_to_CCLE{$c}";
		}
		else{
			print STDOUT ",$c";
		}
	}
	
	print STDOUT "\n";
	
	print STDERR "Writing the drug pair synergy data\n";
	
	DRUG_A: foreach my $drug_A (@valid_drug_A){
		
		my @valid_drug_B = sort keys(%{$pair_data{$drug_A}});
		
		#if( !defined($complete_drugs{$drug_A}) ){
		#	next DRUG_A;
		#}
		
		DRUG_B: foreach my $drug_B (@valid_drug_B){
			
			#if( !defined($complete_drugs{$drug_B}) ){
			#	next DRUG_B;
			#}
		
			# Convert the drug names into CID values
			if( defined($drug_to_cid{$drug_A}) ){
				print STDOUT "$drug_to_cid{$drug_A}";
			}
			else{
				die "Unable to find CID for drug \"$drug_A\" (A)\n";
			}

			if( defined($drug_to_cid{$drug_B}) ){
				print STDOUT ",$drug_to_cid{$drug_B}";
			}
			else{
				die "Unable to find CID for drug \"$drug_B\" (B)\n";
			}

			foreach my $cell (@cell_line_names){
				
				if( !defined($pair_data{$drug_A}{$drug_B}{$cell}) ){
					# The data for this cell line is missing!
					print STDOUT ",";
				}
				else{
					
					# Find the minimum bliss synergy
					my $min_synergy = 1.0e6; # Very large, positive number
					
					my $ave_synergy = 0.0;
					my $norm = 0;
					
					my $best_conc_A = ""; # For debugging
					my $best_conc_B = ""; # For debugging
					
					my @valid_conc_A = keys(%{$pair_data{$drug_A}{$drug_B}{$cell}});
					
					######################################
					# Debugging!!
					#my @debug_conc = sort {$a <=> $b} keys(%{$single_data{$drug_A}{$cell}});

					#foreach my $c (@debug_conc){
					#	print STDERR "$c\t$single_data{$drug_A}{$cell}{$c}\n";
					#}

					#print STDERR "\t\n";

					#for(my $x = 0.001;$x < 10;$x *= 1.01){
					#	print STDERR "$x\t" . interpolate_growth($single_data{$drug_A}{$cell}, $x) . "\n";
					#}

					#print STDERR "\t\n";

					foreach my $conc_A (@valid_conc_A){

						my @valid_conc_B = keys(%{$pair_data{$drug_A}{$drug_B}{$cell}{$conc_A}});
						
						# Interpolate to find the growth at the requested concetration
						my $growth_A = interpolate_growth($single_data{$drug_A}{$cell}, $conc_A);
						
						foreach my $conc_B (@valid_conc_B){
							
							my $growth_B = interpolate_growth($single_data{$drug_B}{$cell}, $conc_B);
						
							my $observed_growth = $pair_data{$drug_A}{$drug_B}{$cell}{$conc_A}{$conc_B};
							
							my $expected_growth = 0.0;
							
							# The bliss model
							if($use_gi){
								
								if(min($growth_A, $growth_B) < 0.0){
									$expected_growth = min($growth_A, $growth_B);
								}
								else{
									$expected_growth = ($growth_A * $growth_B)/100.0;
								}
							}
							else{
								$expected_growth = $growth_A * $growth_B;
							}
							
							# Synergy is defined as: (observed pair growth) - (expected pair growth)
							my $synergy = $observed_growth - $expected_growth;
							
							# DEBUG
							#if( $synergy < -200.0 ){
							#	print STDERR "A = $drug_A\n";
							#	print STDERR "B = $drug_B\n";
							#	print STDERR "cell = $cell\n";
							#	print STDERR "conc A = $conc_A\n";
							#	print STDERR "conc B = $conc_B\n";
							#	print STDERR "growth_A = $growth_A\n";
							#	print STDERR "growth_B = $growth_B\n";
							#	print STDERR "observed_growth = $observed_growth\n";
							#	print STDERR "expected_growth = $expected_growth\n";
							#	print STDERR "synergy = $synergy\n";
							#	print STDERR "$synergy\n";
							#	print STDERR "$observed_growth\t$expected_growth\n\n";
							#}
							# DEBUG
							#if( ($drug_A eq "AZD1775") && ($drug_B eq "MK-8669") ){
							#	print STDERR "$synergy\n";
							#}
							
							if($synergy < $min_synergy){
								
								$min_synergy = $synergy;
								$best_conc_A = $conc_A;
								$best_conc_B = $conc_B;
							}
							
							$ave_synergy += $synergy;
							++$norm;
						}
					}
					
					if($norm > 0){
						$ave_synergy /= $norm;
					}
					
					print STDOUT ",$min_synergy";
					#print STDOUT ",$ave_synergy";
					
					# DEBUG
					#if($drug_A eq "Paclitaxel"){
					#	print STDERR "$drug_A\t$best_conc_A\t$min_synergy\n";
					#}
					
					#if($drug_B eq "Paclitaxel"){
					#	print STDERR "$drug_B\t$best_conc_B\t$min_synergy\n";
					#}
				}
			}
			
			print STDOUT "\n";
		}
	}
}

sub interpolate_growth($$)
{
	my $data = shift;
	my $query = shift;
	
	if($query <= 0.0){
		die "interpolate_growth: Unable to compute log(query)\n";
	}
	
	# Extract all of the concentrations from the data hash
	my @support_conc = keys(%{$data});
	
	# Sort the concetration support points in ascending order
	@support_conc = sort { $a <=> $b } @support_conc;
						
	my $num_support = scalar(@support_conc);
	
	if($num_support == 0){
		die "interpolate_growth: Did not find any concentration values to interpolate!\n";
	}
	
	if($query <= $support_conc[0]){
		return ${$data}{$support_conc[0]};
	}
	
	my $log_query = log($query);
	
	for(my $i = 1;$i < $num_support;++$i){
		
		if($query <= $support_conc[$i]){
			
			if( ($support_conc[$i - 1] <= 0.0) || ($support_conc[$i] <= 0.0) ){
				die "interpolate_growth: Unable to compute the log of the support points\n";
			}
			
			my $x0 = $support_conc[$i - 1];
			my $x1 = $support_conc[$i];
			my $y0 = ${$data}{$x0};
			my $y1 = ${$data}{$x1};
			
			
			# Perform a linear interpolation in log(concentration)
			$x0 = log($x0);
			$x1 = log($x1);
			
			my $y = ($log_query - $x0)*($y1 - $y0)/($x1 - $x0) + $y0;
			
			return $y;
		}
	}
	
	# If we get here, than query is greater than the largest support concentration value
	return ${$data}{$support_conc[$num_support - 1]};
}

sub get_column_index($$)
{
	my $header = shift;
	my $query = shift;
	
	my $num_col = scalar(@{$header});
	
	for(my $i = 0;$i < $num_col;++$i){
		
		if(${$header}[$i] =~ /^$query$/){
			return $i;
		}
	}

	die "Unable to find \"$query\" in header\n";
}

sub median($)
{
	my $input = shift;
	my @data = ();
	
	foreach my $i (@$input){
		if($i ne "NULL"){
			push @data, $i;
		}
	}
	
	my $n = scalar(@data);
	
	if($n == 0){
		die "median: Unable to compute median\n";
	}
	
	@data = sort {$a <=> $b} (@data);
	
	if($n % 2 == 0){		
		return 0.5*($data[ int($n/2) - 1 ] + $data[ int($n/2) ]);
	}
	
	return $data[ int($n/2) ];
}

sub min($$)
{
	my $A = shift;
	my $B = shift;
	
	if($A < $B){
		return $A;
	}
	
	return $B;
}

main();
exit;
