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

# 1) Convert ALMANAC data to use CID instead of NSC drug ids
# 2) Switch from NCI60 to CCLE cell lines

use strict;

# Drug id mapping (this mapping is specific to the NSC drug id's used
# in the ALMANAC study).
my %nsc_to_cid = (
	"740" => "126941",
	"750" => "2478",
	"752" => "2723601",
	"755" => "667490",
	"762" => "4033",
	"1390" => "2094",
	"3053" => "44415057",
	"3088" => "2708",
	"6396" => "5453",
	"8806" => "460612",
	"9706" => "5799",
	"13875" => "2123",
	"14229" => "237",
	"18509" => "137",
	"19893" => "3385",
	"24559" => "163659",
	"25154" => "4842",
	"26271" => "2907",
	"26980" => "5746",
	"27640" => "5790",
	"32065" => "3657",
	"34462" => "6194",
	"38721" => "4211",
	"45388" => "5353562",
	"45923" => "4114",
	"49842" => "13342",
	"63878" => "6253",
	"66847" => "5426",
	"67574" => "5978",
	"71423" => "11683",
	"77213" => "4915",
	"79037" => "3950",
	"82151" => "30323",
	"85998" => "23615975",
	"92859" => "518740",
	"102816" => "1805",
	"105014" => "1546",
	"109724" => "3690",
	"118218" => "3367",
	"119875" => "441203",
	"122758" => "444795",
	"122819" => "54610154",
	"123127" => "32875",
	"125066" => "5360373",
	"125973" => "36314",
	"127716" => "16886",
	"138783" => "65628",
	"141540" => "439525",
	"169780" => "71384",
	"180973" => "2733526",
	"218321" => "40926",
	"226080" => "54600319",
	"241240" => "426756",
	"246131" => "41744",
	"256439" => "3685",
	"256942" => "41867",
	"266046" => "24197464",
	"279836" => "4212",
	"296961" => "2141",
	"362856" => "5394",
	"369100" => "57469",
	"409962" => "2578",
	"606869" => "354624",
	"608210" => "5311497",
	"609699" => "60700",
	"613327" => "3461",
	"628503" => "148124",
	"673596" => "104842",
	"681239" => "387447",
	"686673" => "3011155",
	"698037" => "446556",
	"701852" => "5311",
	"702294" => "259329",
	"707389" => "11354606",
	"712807" => "400633",
	"713563" => "60198",
	"715055" => "123631",
	"718781" => "176870",
	"719276" => "104741",
	"719344" => "2187",
	"719345" => "3902",
	"719627" => "2662",
	"721517" => "68740",
	"732517" => "3062316",
	"733504" => "54608520",
	"737754" => "10113978",
	"743414" => "5291",
	"745750" => "208908",
	"747599" => "644241",
	"747971" => "216239",
	"747972" => "216326",
	"747973" => "6445540",
	"747974" => "5035",
	"749226" => "132971",
	"750690" => "5329102",
	"753082" => "42611257",
	"754143" => "5352062",
	"754230" => "148121",
	"755986" => "24776445",
	"756645" => "54613769",
	"757441" => "6450551",
	"760766" => "3081361",
	# Exclude NC 761431, as this drug appears twice (the other version is
	# NSC 753082) and this version has very little data.
	#"761431" => "42611257",
	"761432" => "9854073",
	"763371" => "25126798"
);

# Cell line name mapping
my %nci60_to_ccle = (
	"786-0" => "CCLE.786O",
	"A498" => "CCLE.A498",
	"A549/ATCC" => "CCLE.A549",
	"ACHN" => "CCLE.ACHN",
	"BT-549" => "CCLE.BT549",
	"CAKI-1" => "CCLE.CAKI1",
	"DU-145" => "CCLE.DU145",
	"EKVX" => "CCLE.EKVX",
	"HCT-116" => "CCLE.HCT116",
	"HCT-15" => "CCLE.HCT15",
	"HL-60(TB)" => "CCLE.HL60",
	"HOP-62" => "CCLE.HOP62",
	"HOP-92" => "CCLE.HOP92",
	"HS 578T" => "CCLE.HS578T",
	"HT29" => "CCLE.HT29",
	"IGROV1" => "CCLE.IGROV1",
	"K-562" => "CCLE.K562",
	"KM12" => "CCLE.KM12",
	"LOX IMVI" => "CCLE.LOXIMVI",
	"MALME-3M" => "CCLE.MALME3M",
	"MCF7" => "CCLE.MCF7",
	"MDA-MB-231/ATCC" => "CCLE.MDAMB231",
	"MDA-MB-435" => "CCLE.MDAMB435S",
	"MDA-MB-468" => "CCLE.MDAMB468",
	"NCI-H226" => "CCLE.NCIH226",
	"NCI-H23" => "CCLE.NCIH23",
	"NCI-H460" => "CCLE.NCIH460",
	"NCI-H522" => "CCLE.NCIH522",
	"OVCAR-3" => "CCLE.NIHOVCAR3",
	"OVCAR-4" => "CCLE.OVCAR4",
	"OVCAR-8" => "CCLE.OVCAR8",
	"PC-3" => "CCLE.PC3",
	"RPMI-8226" => "CCLE.RPMI8226",
	"SF-268" => "CCLE.SF268",
	"SF-295" => "CCLE.SF295",
	"SF-539" => "CCLE.SF539",
	"SK-MEL-28" => "CCLE.SKMEL28",
	"SK-MEL-5" => "CCLE.SKMEL5",
	"SK-OV-3" => "CCLE.SKOV3",
	"SNB-75" => "CCLE.SNB75",
	"SR" => "CCLE.SR786",
	"SW-620" => "CCLE.SW620",
	"T-47D" => "CCLE.T47D",
	"U251" => "CCLE.U251MG",
	"UACC-257" => "CCLE.UACC257",
	"UACC-62" => "CCLE.UACC62",
	"UO-31" => "CCLE.UO31"
);

sub main()
{
	if(scalar(@ARGV) != 1){
		die "Usage: $0 <ALMANAC data file>\n";
	}
	
	my $filename = $ARGV[0];
	
	open(my $INPUT, $filename) or
		die "Unable to open $filename fo reading\n";
	
	my $line = <$INPUT>;
	
	# Remove DOS end-of-line symbols
	$line =~ tr/\r//d;
	
	chomp($line);
	
	my @header = split /,/, $line;
	
	if(scalar(@header) != 29){
		die "Did not read the expected number of columns in the header\n";
	}
	
	my $nsc1_col = -1;
	my $nsc2_col = -1;
	my $cellname_col = -1;
	
	for(my $i = 0;$i < scalar(@header);++$i){
		
		if($header[$i] eq "NSC1"){
			$nsc1_col = $i;
		}
		
		if($header[$i] eq "NSC2"){
			$nsc2_col = $i;
		}
		
		if($header[$i] eq "CELLNAME"){
			$cellname_col = $i;
		}
	}
	
	if($nsc1_col < 0){
		die "Did not read \"NSC1\" in the header\n";
	}
	
	if($nsc2_col < 0){
		die "Did not read \"NSC2\" in the header\n";
	}
	
	if($cellname_col < 0){
		die "Did not read \"CELLNAME\" in the header\n";
	}
	
	# Change the header to indicate that we are using CID, not NSC, drug IDs
	$header[$nsc1_col] = "CID1";
	$header[$nsc2_col] = "CID2";
	
	print STDOUT join(',', @header) . "\n";
	
	my $line_number = 1;
	
	# The ALMANAC data is supposed to be vanilla ASCII, but is contaminated
	my %unmapped_cellname = ();
	
	LINE: while($line = <$INPUT>){
		
		++$line_number;
		
		# Remove DOS end-of-line symbols and the SUB (0x1A) character
		# that appears at least once in the current ALMANAC datafile
		# (likely due to an error when combining multiple database
		# dumps).
		$line =~ tr/\r\x1A//d;
	
		chomp($line);
		
		my @data = split /,/, $line;
		
		if(scalar(@data) != 29){
			die "Did not read the expected number of columns in line $line_number\n";
		}
		
		# Remap the cell name
		if( !defined($nci60_to_ccle{$data[$cellname_col]}) ){
			#die "Unable to lookup CELLNAME ($data[$cellname_col]) to CCLE mapping\n";
			
			++$unmapped_cellname{ $data[$cellname_col] };
			
			# Skip this line
			next LINE;
		}
		
		$data[$cellname_col] = $nci60_to_ccle{$data[$cellname_col]};
		
		# Skip the only sparsely measured version of Vemurafenib (NSC 761431)
		if( ($data[$nsc1_col] eq "761431") || ($data[$nsc2_col] eq "761431") ){
			next LINE;
		}
		
		# Remap drug 1
		if( !defined($nsc_to_cid{$data[$nsc1_col]}) ){
			die "Unable to lookup NSC1 ($data[$nsc1_col]) to CID mapping\n";
		}
		
		$data[$nsc1_col] = $nsc_to_cid{$data[$nsc1_col]};
		
		# Remap drug 2 (if present)
		if($data[$nsc2_col] ne ""){
		
			if( !defined($nsc_to_cid{$data[$nsc2_col]}) ){
				die "Unable to lookup NSC2 ($data[$nsc2_col]) to CID mapping\n";
			}

			$data[$nsc2_col] = $nsc_to_cid{$data[$nsc2_col]};
		}
		
		print STDOUT join(',', @data) . "\n";
	}

	close($INPUT);
	
	my @cellname = keys(%unmapped_cellname);
	
	foreach my $c (@cellname){
		
		print STDERR "Did not find $c : $unmapped_cellname{$c} times\n";
	}
}

main();
exit;

