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

#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_set>

#include <string.h>

#include "zlib.h"
#include "gemini.h"
#include "parse_csv.h"

using namespace std;

deque<DrugPairData> parse_pair_drug(string &m_filename, vector<string> &m_cell_names,
                                      const float &m_normalize_by /*= 1.0*/)
{
	deque<DrugPairData> ret;
	
	ifstream fin( m_filename.c_str() );
	
	if(!fin){
		throw __FILE__ ":parse_pair_drug: Unable to open drug pair growth file for reading";
	}
	
	string line;
	
	// Read the header
	if( !getline(fin, line) ){
		throw __FILE__ ":parse_pair_drug: Unable to read header";
	}
	
	const vector<string> header = split(line, ',');
	
	const size_t num_col = header.size();
	
	// Make sure that there is at least two drug columns and one cell name column
	if(num_col < 3){
		throw __FILE__ ":parse_pair_drug: Did not read the minimum number of columns";
	}
	
	const size_t drug1_col = get_col(header, "DRUG1");
	const size_t drug2_col = get_col(header, "DRUG2");
	const size_t first_cell_col = drug2_col + 1;
	
	// The header defines the cell line order
	const size_t num_cell = num_col - 2;
	
	m_cell_names.resize(num_cell);
	
	for(size_t i = 0;i < num_cell;++i){
		m_cell_names[i] = header[first_cell_col + i];
	}
	
	while( getline(fin, line) ){
		
		const vector<string> data = split(line, ',');
		
		if(data.size() != num_col){
			throw __FILE__ ":parse_pair_drug: Did not read the expected number of columns";
		}
		
		ret.push_back( DrugPairData() );
		
		DrugPairData &ref = ret.back();
		
		ref.drug_id_1 = data[drug1_col];
		ref.drug_id_2 = data[drug2_col];
		
		ref.growth.resize(num_cell);
		
		for(size_t i = 0;i < num_cell;++i){
			
			// Is this missing data?
			if( data[first_cell_col + i].empty() ){
				ref.growth[i] = MISSING_DATA;
			}
			else{
				ref.growth[i] = atof( data[first_cell_col + i].c_str() )/m_normalize_by;
			}
		}
	}
	
	return ret;
}

void parse_single_drug(MULTIMAP<string, DrugFeatures> &m_single,
	const string &m_filename, 
	const float &m_normalize_by /* = 1.0*/)
{
	ifstream fin( m_filename.c_str() );
	
	if(!fin){
		throw __FILE__ ":parse_single_drug: Unable to open single drug growth file for reading";
	}
	
	string line;
	
	// Read the header
	if( !getline(fin, line) ){
		throw __FILE__ ":parse_single_drug: Unable to read header";
	}
	
	const vector<string> header = split(line, ',');
	
	const size_t num_col = header.size();
	
	const size_t id_col = get_col(header, "DRUG");
    	const size_t fingerprint_col = get_col(header, "FINGERPRINT");
	
	// Is this a fingerprint file?
	if( (id_col != COLUMN_NOT_FOUND) && (fingerprint_col != COLUMN_NOT_FOUND) ){
		
		MAP<string, DrugFeatures> local;
		
		while( getline(fin, line) ){
		
			const vector<string> data = split(line, ',');

			if(data.size() != num_col){
				throw __FILE__ ":parse_single_drug: Did not read the expected number of columns";
			}
			
			// Check for erroneously repeated drugs
			if( local.find(data[id_col]) != local.end() ){
				throw __FILE__ ":parse_single_drug: Found a duplicate drug";
			}
			
			const size_t num_fingerprints = data[fingerprint_col].size();
			
			vector<float> growth(num_fingerprints);
			
			for(size_t i = 0;i < num_fingerprints;++i){
				
				switch(data[fingerprint_col][i]){
					case '0':
						growth[i] = 0.0;
						break;
					case '1':
						growth[i] = 1.0;
						break;
				};
			}
			
			local.insert( make_pair(data[id_col], growth) );
		}
		
		// Since this function needs to return *pairs* of single drugs, combine the single drug data into a
		// a multimap that contains all *pairs* of drugs
		
		const deque<string> drugs = keys(local);
		
		for(deque<string>::const_iterator i = drugs.begin();i != drugs.end();++i){
			
			const MAP<string, DrugFeatures>::const_iterator iter1 = local.find(*i);
				
			if( iter1 == local.end() ){
				throw __FILE__ ":parse_single_drug: Unable to find first drug";
			}
				
			for(deque<string>::const_iterator j = i + 1;j != drugs.end();++j){
								
				// Note that the keys() function sorts the drugs in ascending order,
				// so we will pack the drugs as "i + j" (as opposed to "j + i")
				
				const string name = *i + "+" + *j;
        
				const MAP<string, DrugFeatures>::const_iterator iter2 = local.find(*j);
				
				if( iter2 == local.end() ){
					throw __FILE__ ":parse_single_drug: Unable to find second drug";
				}
				
				m_single.insert( make_pair(name, iter1->second) );
				m_single.insert( make_pair(name, iter2->second) );
			}
		}
		
		return;
	}
	
	const size_t drug1_col = get_col(header, "DRUG1");
    	const size_t drug2_col = get_col(header, "DRUG2");
    	const size_t drug_col = get_col(header, "DRUG");
	const size_t first_cell_col = drug_col + 1;
	
	if(drug_col == COLUMN_NOT_FOUND){
		throw __FILE__ ":parse_single_drug: Unable to find column \"DRUG\"";
	}
	
	if( (drug1_col == COLUMN_NOT_FOUND) || (drug2_col == COLUMN_NOT_FOUND) ){
		
		// Make sure that *both* the DRUG1 and DRUG2 columns are missing
		if(drug1_col != drug2_col){
			throw __FILE__ ":parse_single_drug: Missing either \"DRUG1\" or \"DRUG2\" columns, but not both";
		}
		
		// If we get here, then this is a single drug feature file that only contains features on a per-drug basis
		
		if(num_col <= 1){
			throw __FILE__ ":parse_single_drug: Did not read enough feature columns";
		}

		const size_t num_feature = num_col - 1;
		
		MAP<string, DrugFeatures> local;
		
		while( getline(fin, line) ){
		
			const vector<string> data = split(line, ',');

			if(data.size() != num_col){
				throw __FILE__ ":parse_single_drug: Did not read the expected number of columns";
			}

			// Check for erroneously repeated drugs
			if( local.find(data[drug_col]) != local.end() ){
				throw __FILE__ ":parse_single_drug: Found a duplicate drug";
			}

			// Initialize all growth values as missing data
			vector<float> growth(num_feature, MISSING_DATA);

			for(size_t i = 0;i < num_feature;++i){

				// Since we've initialized all values as missing data, we can just skip
				// missing data
				if( !data[first_cell_col + i].empty() ){
					growth[i] = atof( data[first_cell_col + i].c_str() )/m_normalize_by;
				}
			}

        		local.insert( make_pair(data[drug_col], growth) );
		}
		
		// Since this function needs to return *pairs* of single drugs, combine the single drug data into a
		// a multimap that contains all *pairs* of drugs
		
		const deque<string> drugs = keys(local);
		
		for(deque<string>::const_iterator i = drugs.begin();i != drugs.end();++i){
			
			const MAP<string, DrugFeatures>::const_iterator iter1 = local.find(*i);
				
			if( iter1 == local.end() ){
				throw __FILE__ ":parse_single_drug: Unable to find first drug";
			}
				
			for(deque<string>::const_iterator j = i + 1;j != drugs.end();++j){
								
				// Note that the keys() function sorts the drugs in ascending order,
				// so we will pack the drugs as "i + j" (as opposed to "j + i")
				
				const string name = *i + "+" + *j;
        
				const MAP<string, DrugFeatures>::const_iterator iter2 = local.find(*j);
				
				if( iter2 == local.end() ){
					throw __FILE__ ":parse_single_drug: Unable to find second drug";
				}
				
				m_single.insert( make_pair(name, iter1->second) );
				m_single.insert( make_pair(name, iter2->second) );
			}
		}
		
		return;
	}
	
	if(num_col <= 3){
		throw __FILE__ ":parse_single_drug: Did not read enough feature columns";
	}
	
	const size_t num_feature = num_col - 3;
	
	while( getline(fin, line) ){
		
		const vector<string> data = split(line, ',');
		
		if(data.size() != num_col){
			throw __FILE__ ":parse_single_drug: Did not read the expected number of columns";
		}
		
		// Initialize all growth values as missing data
		vector<float> growth(num_feature, MISSING_DATA);
		
		for(size_t i = 0;i < num_feature;++i){
			
			// Since we've initialized all values as missing data, we can just skip
			// missing data
			if( !data[first_cell_col + i].empty() ){
				growth[i] = atof( data[first_cell_col + i].c_str() )/m_normalize_by;
			}
		}
		
        	const string name = data[drug1_col] + "+" + data[drug2_col];
        
		m_single.insert( make_pair(name, growth) );
	}
}

void parse_fingerprints(MAP<string, DrugFeatures> &m_single, const string &m_filename)
{
	// Allow reading (gzip) compressed or uncompressed files
	gzFile fin = gzopen(m_filename.c_str(), "r");
	
	if(fin == NULL){
		throw __FILE__ ":parse_fingerprints: Unable to open input file for reading";
	}

	const int buffer_len = 4096;
	char buffer[buffer_len];
	
	size_t line_number = 0;
	
	size_t num_col = 0;
	size_t nsc_col = 0;
	size_t fingerprint_col = 0;
	
	deque< pair<size_t, size_t> > cell_line_col_to_index;
	
	size_t num_features = 0;
	
	unordered_set<string> current_drugs;	
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		++line_number;
		
		const size_t len = strlen(buffer);
		
		if( len == (buffer_len - 1) ){
			
			cerr << "Error parsing " << m_filename << " at line #" << line_number << endl;
			throw __FILE__ ":parse_fingerprints: Potential buffer overflow!";
		}
		
		// Note that while gzgets stores the terminal '\n', we do not need to explicitly remove
		// it because the split() function skips any terminal '\n' or '\r'.
		if( num_col == 0 ){
			
			// Read the header
			vector<string> header = split(buffer, ',');
			
			if( header.empty() ){
				throw __FILE__ ":parse_fingerprints: Empty header";
			}
			
			num_col = header.size();
			
			nsc_col = get_col(header, "NSC");
			
			if(nsc_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_fingerprints: Unable to find \"NSC\" column header";
			}
			
			fingerprint_col = get_col(header, "fingerprint");
			
			if(fingerprint_col == COLUMN_NOT_FOUND){
				
				cerr << "buffer = " << buffer << endl;
				throw __FILE__ ":parse_fingerprints: Unable to find \"fingerprint\" column header";
			}
		}
		else{
			
			const vector<string> line = split(buffer, ',');
			
			if(line.size() != num_col){
				
				cerr << "Did not read the correct number of columns @ line " 
					<< line_number << endl;
				cerr << "Found " << line.size() << " columns, but expected " << num_col << endl;
				throw __FILE__ ":parse_fingerprints: I/O error";
			}
			
			// Skip missing drug descriptor data (many drugs do not have any structure-based
			// descriptor data available).
			if(line[fingerprint_col] == "-"){
				continue;
			}
			
			const string drug_id = line[nsc_col];
			
			// Convert the bitstring into vector<bool>
			if( line[fingerprint_col].empty() ){
				throw __FILE__ ":parse_fingerprints: Unable to read fingerprint data";
			}
			
			if(num_features == 0){
			
				num_features = line[fingerprint_col].size();
				
				if(num_features == 0){
					throw __FILE__ ":parse_fingerprints: Found a zero length fingerprint";
				}
			}
			else{
				if( num_features != line[fingerprint_col].size() ){
					throw __FILE__ ":parse_fingerprints: Did not read expected number of descriptor features";
				}
			}
			
			vector<float> data(num_features);
			
			for(size_t i = 0;i < num_features;++i){
				
				switch(line[fingerprint_col][i]){
					case '1':
						data[i] = 1.0;
						break;
					case '0':
						data[i] = 0.0;
						break;
					default:
						throw __FILE__ ":parse_fingerprints: Invalid fingerprint symbol";
				}
			}
			
			// Make sure that we don't have any duplicate entries
			if( current_drugs.find(drug_id) != current_drugs.end() ){
				throw __FILE__ ":parse_fingerprints: Duplicate entry";
			}
			
			current_drugs.insert(drug_id);
			
			MAP<string, DrugFeatures>::iterator iter = m_single.find(drug_id);
			
			if( iter != m_single.end() ){
			
				// Concatinate the features
				iter->second.insert( iter->second.end(), data.begin(), data.end() );
			}
			else{
				m_single[drug_id] = data;
			}
		}
	}
	
	gzclose(fin);
}
