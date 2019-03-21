#!/bin/sh

#Copyright (c) 2018 Los Alamos National Security, LLC.
#All rights reserved.
#
#Copyright 2018. Los Alamos National Security, LLC. This software was produced under U.S. Government 
#contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos 
#National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, 
#reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC 
#MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If 
#software is modified to produce derivative works, such modified software should be clearly marked, so 
#as not to confuse it with the version available from LANL.
#
#This work has been supported in part by the Joint Design of Advanced Computing Solutions for Cancer 
#(JDACS4C) program established by the U.S. Department of Energy (DOE) and the National Cancer Institute 
#(NCI) of the National Institutes of Health. 
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
#associated documentation files (the "Software"), to deal in the Software without restriction, including 
#without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to 
#the following conditions: 
#
#The above copyright notice and this permission notice shall be included in all copies or substantial 
#portions of the Software.
#
#THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
#OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
#FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR 
#CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
#OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
#EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
#Author: Jason D. Gans (jgans@lanl.gov)

OUTPUT_BLISS_DIR=data_bliss/
OUTPUT_LOEWE_DIR=data_loewe/
MERCK_DIR=.
ALMANAC_DIR=.

merck_cells()
{
	# Extract the per-cell line synergy files for the Merck data
	echo "Extracting per-cell line synergy files for the Merck data (w/ CCLE & CID)"

	for CELL in "CCLE.HT144" "CCLE.A427" "CCLE.CAOV3" "CCLE.COLO320" "CCLE.DLD1" "CCLE.EFM192A" \
		"CCLE.OV90" "CCLE.NCIH520" "CCLE.MDAMB436" \
		"CCLE.LOVO" "CCLE.RPMI7951" "CCLE.A375" "CCLE.HT29" "CCLE.ZR751" "CCLE.T47D" "CCLE.SW620" \
		"CCLE.LNCAPCLONEFGC" "CCLE.MSTO211H" "CCLE.KPL1" "CCLE.SKMEL30" "CCLE.HCT116" "CCLE.RKO" \
		"CCLE.NCIH1650" "CCLE.VCAP" "CCLE.NCIH23" "CCLE.SKOV3" "CCLE.NCIH2122" "CCLE.ES2" \
		"CCLE.NIHOVCAR3" "CCLE.PA1" "CCLE.SW837" "CCLE.SKMES1" "CCLE.UACC62" "CCLE.A2058" \
		"CCLE.NCIH460" "CCLE.A2780"
	do
		# Gave up trying to figure out how to encode "/", so use [:punct:] instead
		CELL_NAME=`echo $CELL | tr "[:blank:][:punct:]" "_"`

		echo "--------Merck: $CELL --------"

		./synergy_search -n 1 \
			-o $OUTPUT_BLISS_DIR/Merck.$CELL_NAME.csv \
			--dependent $MERCK_DIR/merck_pair.csv --include "$CELL"
	done
}

nci60_cells()
{
	# Extract the per-cell line synergy files for the ALMANAC data
	echo "Extracting per-cell line synergy files for the ALMANAC data (w/ CCLE & CID)"

	for CELL in "CCLE.786O" "CCLE.A498" "CCLE.A549" "CCLE.ACHN" "CCLE.BT549" "CCLE.CAKI1" \
		"CCLE.DU145" "CCLE.EKVX" "CCLE.HCT116" "CCLE.HCT15" "CCLE.HL60" "CCLE.HOP62" \
		"CCLE.HOP92" "CCLE.HS578T" "CCLE.HT29" "CCLE.IGROV1" "CCLE.K562" "CCLE.KM12" \
		"CCLE.LOXIMVI" "CCLE.MALME3M" "CCLE.MCF7" "CCLE.MDAMB231" "CCLE.MDAMB435S" \
		"CCLE.MDAMB468" "CCLE.NCIH226" "CCLE.NCIH23" "CCLE.NCIH460" "CCLE.NCIH522" \
		"CCLE.NIHOVCAR3" "CCLE.OVCAR4" "CCLE.OVCAR8" "CCLE.PC3" "CCLE.RPMI8226" "CCLE.SF268" "CCLE.SF295" \
		"CCLE.SF539" "CCLE.SKMEL28" "CCLE.SKMEL5" "CCLE.SKOV3" "CCLE.SNB75" "CCLE.SR786" \
		"CCLE.SW620" "CCLE.T47D" "CCLE.U251MG" "CCLE.UACC257" "CCLE.UACC62" "CCLE.UO31"
	do
		# Gave up trying to figure out how to encode "/", so use [:punct:] instead
		CELL_NAME=`echo $CELL | tr "[:blank:][:punct:]" "_"`

		echo "--------ALMANAC: $CELL --------"

		# Use the Loewe synergy
		./synergy_search -n 1 \
			-o $OUTPUT_LOEWE_DIR/ALMANAC.$CELL_NAME.csv \
			--dependent $ALMANAC_DIR/ALMANAC_CID_CCLE_loewe_ave_synergy.csv --include "$CELL"
		
		# Use the Bliss synergy
		./synergy_search -n 1 \
			-o $OUTPUT_BLISS_DIR/ALMANAC.$CELL_NAME.csv \
			--dependent $ALMANAC_DIR/ALMANAC_CID_CCLE_bliss_ave_synergy.csv --include "$CELL"
		
	done
}

nci60_cells
merck_cells

