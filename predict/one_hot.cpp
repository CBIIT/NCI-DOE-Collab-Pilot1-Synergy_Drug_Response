// Copyright (c) 2018 Los Alamos National Security, LLC.
// All rights reserved.
// 
// Copyright 2018. Los Alamos National Security, LLC. This software was produced under U.S. Government 
// contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos 
// National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, 
// reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC 
// MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If 
// software is modified to produce derivative works, such modified software should be clearly marked, so 
// as not to confuse it with the version available from LANL.
// 
// This work has been supported in part by the Joint Design of Advanced Computing Solutions for Cancer 
// (JDACS4C) program established by the U.S. Department of Energy (DOE) and the National Cancer Institute 
// (NCI) of the National Institutes of Health. 
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
// associated documentation files (the "Software"), to deal in the Software without restriction, including 
// without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
// copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to 
// the following conditions: 
// 
// The above copyright notice and this permission notice shall be included in all copies or substantial 
// portions of the Software.
// 
// THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
// OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR 
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Author: Jason D. Gans (jgans@lanl.gov)

#include "rx.h"
#include "deque_set.h"
#include "keys.hpp"
#include <iostream>

using namespace std;

void remove_prefix(string &m_name, const string &m_prefix);
string cell_to_tissue(const string &m_cell);

// Return the dimensionality of the one hot encoded feature vector
size_t one_hot_encode(MAP<string, vector<float> > &m_features, 
	const vector<string> &m_categories)
{

	// Collect all of the categories from both the existing features and the input
	// categories
	deque<string> cat;
	
	for(MAP<string, vector<float> >::const_iterator i = m_features.begin();i != m_features.end();++i){
		cat.push_back(i->first);
	}
	
	for(vector<string>::const_iterator i = m_categories.begin();i != m_categories.end();++i){
		cat.push_back(*i);
	}

	// Make the categories unique
	make_set(cat);
	
	const size_t num_cat = cat.size();
	
	if(num_cat == 0){
		throw __FILE__ ":one_hot_encode: Did not find any categories to one hot encode!";
	}
	
	MAP<string, vector<float> > ret;
	
	size_t index = 0;
	
	// The ordering of encoded features is arbitrary!
	for(deque<string>::const_iterator i = cat.begin();i != cat.end();++i, ++index){
		
		vector<float> encoded(num_cat);
		
		encoded[index] = 1;
		
		ret[*i] = encoded;
	}
	
	// Replace the old feature set with the new feature set
	m_features = ret;
	
	return num_cat;
}

// Return the list of tissues used to one hot encode
vector<string> one_hot_encode_by_tissue(MAP<string, vector<float> > &m_features, 
	const vector<string> &m_categories)
{

	// Collect all of the cell lines from both the existing features and the input
	// categories
	deque<string> cells;
	
	for(MAP<string, vector<float> >::const_iterator i = m_features.begin();i != m_features.end();++i){
		cells.push_back(i->first);
	}
	
	for(vector<string>::const_iterator i = m_categories.begin();i != m_categories.end();++i){
		cells.push_back(*i);
	}

	make_set(cells);
	
	const size_t num_cells = cells.size();
	
	if(num_cells == 0){
		throw __FILE__ ":one_hot_encode_by_tissue: Did not find any cells to one hot encode!";
	}
	
	MAP<string, size_t> tissue_index;
	
	for(deque<string>::const_iterator i = cells.begin();i != cells.end();++i){
	
		const string t = cell_to_tissue(*i);
		
		if( tissue_index.find(t) == tissue_index.end() ){
			
			// Compute the one hot encoding index before adding it to the
			// tissue_index (adding the index in a single instruction leads to
			// poor readability, as the size of tissue_index gets incremented
			// *before* the insertion).
			const size_t index = tissue_index.size();
			
			tissue_index[t] = index;
		}
	}
	
	const size_t num_tissues = tissue_index.size();
	
	if(num_tissues == 0){
		throw __FILE__ ":one_hot_encode_by_tissue: Did not find any tissues to one hot encode!";
	}
	
	MAP<string, vector<float> > ret;
	
	// The ordering of encoded features is arbitrary!
	for(deque<string>::const_iterator i = cells.begin();i != cells.end();++i){
		
		MAP<string, size_t>::const_iterator iter = tissue_index.find( cell_to_tissue(*i) );
		
		if( iter == tissue_index.end() ){
			throw __FILE__ ":one_hot_encode_by_tissue: Unable to lookup tissue name to encoding index";
		}
		
		if(iter->second >= num_tissues){
			throw __FILE__ ":one_hot_encode_by_tissue: Encoding index overflow";
		}
		
		vector<float> encoded(num_tissues);
		
		encoded[iter->second] = 1;
		
		ret[*i] = encoded;
	}
	
	// Replace the old feature set with the new feature set
	m_features = ret;
	
	return keys(tissue_index);
}

string cell_to_tissue(const string &m_cell)
{
	string ret(m_cell);
	
	// Remove any dataset prefixes
	remove_prefix(ret, "CCLE.");
	remove_prefix(ret, "NCI60.");
	remove_prefix(ret, "CTRP.");
	remove_prefix(ret, "GDSC.");
	remove_prefix(ret, "gCSI.");
	remove_prefix(ret, "GDC.");
	remove_prefix(ret, "NCIPDM.");

	/////////////////////////////////////////////////
	// ALMANAC cell lines
	/////////////////////////////////////////////////
	if( (ret == "SF-539") || (ret == "SF539") ){

		ret = "TISSUE_CNS";
		return ret;
	}

	if( (ret == "MDA-MB-435") || (ret == "MDAMB435S") ){

		ret = "TISSUE_ME";
		return ret;
	}

	if(ret == "NCI-H322M"){

		ret = "TISSUE_LC";
		return ret;
	}

	if( (ret == "SR") || (ret == "SR786") ){

		ret = "TISSUE_LE";
		return ret;
	}

	if(ret == "SN12C"){

		ret = "TISSUE_RE";
		return ret;
	}

	if( (ret == "PC-3") || (ret == "PC3") ){

		ret = "TISSUE_PR";
		return ret;
	}

	if( (ret == "UO-31") || (ret == "UO31") ){

		ret = "TISSUE_RE";
		return ret;
	}

	if( (ret == "HS 578T") || (ret == "HS578T") ){

		ret = "TISSUE_BR";
		return ret;
	}

	if(ret == "TK-10"){

		ret = "TISSUE_RE";
		return ret;
	}

	if(ret == "ACHN"){

		ret = "TISSUE_RE";
		return ret;
	}

	if( (ret == "CAKI-1") || (ret == "CAKI1") ){

		ret = "TISSUE_RE";
		return ret;
	}

	if( (ret == "OVCAR-4") || (ret == "OVCAR4") ){

		ret = "TISSUE_OV";
		return ret;
	}

	if(ret == "RXF 393"){

		ret = "TISSUE_RE";
		return ret;
	}

	if(ret == "CCRF-CEM"){

		ret = "TISSUE_LE";
		return ret;
	}

	if(ret == "SK-MEL-2"){

		ret = "TISSUE_ME";
		return ret;
	}

	if( (ret == "MDA-MB-468") || (ret == "MDAMB468") ){

		ret = "TISSUE_BR";
		return ret;
	}

	if( (ret == "786-0") || (ret == "786O") ){

		ret = "TISSUE_RE";
		return ret;
	}

	if(ret == "M14"){

		ret = "TISSUE_ME";
		return ret;
	}

	if( (ret == "HOP-92") || (ret == "HOP92") ){

		ret = "TISSUE_LC";
		return ret;
	}

	if(ret == "COLO 205"){

		ret = "TISSUE_CO";
		return ret;
	}

	if(ret == "OVCAR-5"){

		ret = "TISSUE_OV";
		return ret;
	}

	if( (ret == "SNB-75") || (ret == "SNB75") ){

		ret = "TISSUE_CNS";
		return ret;
	}

	if(ret == "EKVX"){

		ret = "TISSUE_LC";
		return ret;
	}

	if( (ret == "SK-MEL-28") || (ret == "SKMEL28") ){

		ret = "TISSUE_ME";
		return ret;
	}

	if( (ret == "NCI-H226") || (ret == "NCIH226") ){

		ret = "TISSUE_LC";
		return ret;
	}

	if(ret == "IGROV1"){

		ret = "TISSUE_OV";
		return ret;
	}

	if( (ret == "HOP-62") || (ret == "HOP62") ){

		ret = "TISSUE_LC";
		return ret;
	}

	if( (ret == "SF-295") || (ret == "SF295") ){

		ret = "TISSUE_CNS";
		return ret;
	}

	if( (ret == "K-562") || (ret == "K562") ){

		ret = "TISSUE_LE";
		return ret;
	}

	if( (ret == "SK-MEL-5") || (ret == "SKMEL5") ){

		ret = "TISSUE_ME";
		return ret;
	}

	if(ret == "NCI/ADR-RES"){

		ret = "TISSUE_OV";
		return ret;
	}

	if( (ret == "NCI-H522") || (ret == "NCIH522") ){

		ret = "TISSUE_LC";
		return ret;
	}

	if( (ret == "A549/ATCC") || (ret == "A549") ){

		ret = "TISSUE_LC";
		return ret;
	}

	if(ret == "SW-620"){

		ret = "TISSUE_CO";
		return ret;
	}

	if(ret == "HCT-116"){

		ret = "TISSUE_CO";
		return ret;
	}

	if( (ret == "OVCAR-8") || (ret == "OVCAR8") ){

		ret = "TISSUE_OV";
		return ret;
	}

	if(ret == "OVCAR-3"){

		ret = "TISSUE_OV";
		return ret;
	}

	if(ret == "SK-OV-3"){

		ret = "TISSUE_OV";
		return ret;
	}

	if( (ret == "UACC-257") || (ret == "UACC257") ){

		ret = "TISSUE_ME";
		return ret;
	}

	if(ret == "MCF7"){

		ret = "TISSUE_BR";
		return ret;
	}

	if(ret == "HCC-2998"){

		ret = "TISSUE_CO";
		return ret;
	}

	if( (ret == "RPMI-8226") || (ret == "RPMI8226") ){

		ret = "TISSUE_LE";
		return ret;
	}

	if( (ret == "HL-60(TB)") || (ret == "HL60") ){

		ret = "TISSUE_LE";
		return ret;
	}

	if(ret == "KM12"){

		ret = "TISSUE_CO";
		return ret;
	}

	if( (ret == "MDA-MB-231/ATCC") || (ret == "MDAMB231") ){

		ret = "TISSUE_BR";
		return ret;
	}

	if( (ret == "DU-145") || (ret == "DU145") ){

		ret = "TISSUE_PR";
		return ret;
	}

	if(ret == "NCI-H460"){

		ret = "TISSUE_LC";
		return ret;
	}

	if( (ret == "BT-549") || (ret == "BT549") ){

		ret = "TISSUE_BR";
		return ret;
	}

	if( (ret == "T-47D") || (ret == "T47D") ){

		ret = "TISSUE_BR";
		return ret;
	}

	if(ret == "SNB-19"){

		ret = "TISSUE_CNS";
		return ret;
	}

	if( (ret == "HCT-15") || (ret == "HCT15") ){

		ret = "TISSUE_CO";
		return ret;
	}

	if(ret == "HT29"){

		ret = "TISSUE_CO";
		return ret;
	}

	if( (ret == "MALME-3M") || (ret == "MALME3M") ){

		ret = "TISSUE_ME";
		return ret;
	}

	if(ret == "A498"){

		ret = "TISSUE_RE";
		return ret;
	}

	if( (ret == "U251") || (ret == "U251MG") ){

		ret = "TISSUE_CNS";
		return ret;
	}

	if(ret == "NCI-H23"){

		ret = "TISSUE_LC";
		return ret;
	}

	if( (ret == "SF-268") || (ret == "SF268") ){

		ret = "TISSUE_CNS";
		return ret;
	}

	if( (ret == "UACC-62") || (ret == "UACC62") ){

		ret = "TISSUE_ME";
		return ret;
	}

	if( (ret == "LOX IMVI") || (ret == "LOXIMVI") ){

		ret = "TISSUE_ME";
		return ret;
	}

	if(ret == "MOLT-4"){

		ret = "TISSUE_LE";
		return ret;
	}

	//////////////////////////////////////////////////////////
	// Merck cell lines. Tissue type determined by searching: 
	//	Supplementary online table 1: http://mct.aacrjournals.org/content/15/6/1155.figures-only
	// 	COSMIC database: http://cancer.sanger.ac.uk
	//	Cellosaurus: https://web.expasy.org/cellosaurus
	//	Google	
	//////////////////////////////////////////////////////////

	if(ret == "A2058"){ // COSMIC

		ret = "TISSUE_ME";
		return ret;
	}

	if(ret == "A2780"){ // COSMIC

		ret = "TISSUE_OV";
		return ret;
	}

	if(ret == "A375"){ // COSMIC

		ret = "TISSUE_ME";
		return ret;
	}

	if(ret == "A427"){ // COSMIC

		ret = "TISSUE_LC";
		return ret;
	}

	if(ret == "CAOV3"){ // Google search

		ret = "TISSUE_OV";
		return ret;
	}

	if(ret == "COLO320DM"){ // Cellosaurus

		ret = "TISSUE_CO";
		return ret;
	}

	if(ret == "COLO320"){ // CCLE

		ret = "TISSUE_CO";
		return ret;
	}

	if(ret == "DLD1"){ // Cellosaurus

		ret = "TISSUE_CO";
		return ret;
	}

	if(ret == "EFM192B"){ // Cellosaurus

		ret = "TISSUE_BR";
		return ret;
	}

	if(ret == "ES2"){ // Supplementary table 1

		ret = "TISSUE_OV";
		return ret;
	}

	if(ret == "HCT116"){ // Supplementary table 1

		ret = "TISSUE_CO";
		return ret;
	}

	if(ret == "HT144"){ // Supplementary table 1

		ret = "TISSUE_ME";
		return ret;
	}

	if(ret == "HT29"){ // Supplementary table 1

		ret = "TISSUE_CO";
		return ret;
	}

	if(ret == "KPL1"){ // Supplementary table 1

		ret = "TISSUE_BR";
		return ret;
	}

	if(ret == "LNCAP"){ // Supplementary table 1

		ret = "TISSUE_PR";
		return ret;
	}

	if(ret == "LOVO"){ // Supplementary table 1

		ret = "TISSUE_CO";
		return ret;
	}

	if(ret == "MDAMB436"){ // Supplementary table 1

		ret = "TISSUE_BR";
		return ret;
	}

	if(ret == "MSTO"){ // Supplementary table 1

		ret = "TISSUE_LC";
		return ret;
	}

	if(ret == "NCIH1650"){ // Supplementary table 1

		ret = "TISSUE_LC";
		return ret;
	}

	if(ret == "NCIH2122"){ // Supplementary table 1

		ret = "TISSUE_LC";
		return ret;
	}

	if(ret == "NCIH23"){ // Supplementary table 1

		ret = "TISSUE_LC";
		return ret;
	}

	if(ret == "NCIH460"){ // Cellosaurus

		ret = "TISSUE_LC";
		return ret;
	}

	if(ret == "NCIH520"){ // Supplementary table 1

		ret = "TISSUE_LC";
		return ret;
	}

	if(ret == "OCUBM"){ // Supplementary table 1

		ret = "TISSUE_BR";
		return ret;
	}

	if(ret == "OV90"){ // Supplementary table 1

		ret = "TISSUE_OV";
		return ret;
	}

	if( (ret == "OVCAR3") || (ret == "NIHOVCAR3") ){ // Supplementary table 1

		ret = "TISSUE_OV";
		return ret;
	}

	if(ret == "PA1"){ // Supplementary table 1

		ret = "TISSUE_OV";
		return ret;
	}

	if(ret == "RKO"){ // Supplementary table 1

		ret = "TISSUE_CO";
		return ret;
	}

	if(ret == "RPMI7951"){ // Supplementary table 1

		ret = "TISSUE_ME";
		return ret;
	}

	if(ret == "SKMEL30"){ // Supplementary table 1

		ret = "TISSUE_ME";
		return ret;
	}

	if(ret == "SKMES1"){ // Supplementary table 1

		ret = "TISSUE_LC";
		return ret;
	}

	if(ret == "SKOV3"){ // Supplementary table 1

		ret = "TISSUE_OV";
		return ret;
	}

	if(ret == "SW620"){ // Supplementary table 1

		ret = "TISSUE_CO";
		return ret;
	}

	if(ret == "SW837"){ // Supplementary table 1

		ret = "TISSUE_CO";
		return ret;
	}

	if(ret == "UACC62"){ // Supplementary table 1

		ret = "TISSUE_ME";
		return ret;
	}

	if(ret == "UWB1289"){ // Supplementary table 1

		ret = "TISSUE_OV";
		return ret;
	}

	if(ret == "UWB1289BRCA1"){ // Google

		ret = "TISSUE_OV";
		return ret;
	}

	if(ret == "VCAP"){ // Supplementary table 1

		ret = "TISSUE_PR";
		return ret;
	}

	if(ret == "ZR751"){ // Supplementary table 1

		ret = "TISSUE_BR";
		return ret;
	}

	if(ret == "LNCAPCLONEFGC"){ // Google -> Cellosaurus

		ret = "TISSUE_PR";
		return ret;
	}

	if(ret == "MSTO211H"){ //Cellosaurus

		ret = "TISSUE_LC";
		return ret;
	}

	if(ret == "EFM192A"){ // CCLE

		ret = "TISSUE_BR";
		return ret;
	}

	cerr << "Unknown cell line type: \"" << ret << "\"; can't find corresponding tissue type" << endl;
	throw __FILE__ ":tissue_type: Unable to find cell line name";
	
	return ret;
}

void remove_prefix(string &m_name, const string &m_prefix)
{
	string::size_type loc = m_name.find(m_prefix);
	
	if(loc == 0){
	
		const size_t prefix_len = m_prefix.size();
		
		m_name = m_name.substr(prefix_len, m_name.size() - prefix_len);
	}
}
