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

#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>

#include "rx.h"
#include "parse_csv.h"

using namespace std;

bool is_missing(const string &m_str);

MAP< string, vector<float> > parse_fingerprint(const string &m_filename,
	vector<string> &m_feature_id,
	const string &m_drug, const string &m_fingerprint,
	const char m_delim /*= ','*/)
{
	MAP< string, vector<float> > ret;
	
	ifstream fin( m_filename.c_str() );
	
	if(!fin){
		throw __FILE__ ":parse_fingerprint: Unable to open input file for reading";
	}
	
	string line;
	
	// Read the header
	if( !getline(fin, line) ){
		throw __FILE__ ":parse_fingerprint: Unable to read header";
	}
	
	const vector<string> header = split(line, m_delim);
	
	const size_t num_col = header.size();
	
	const size_t drug_col = get_col(header, m_drug);
	
	if(drug_col == COLUMN_NOT_FOUND){
		throw __FILE__ ":parse_fingerprint: Did not find the requested drug name column";
	}
	
	const size_t fingerprint_col = get_col(header, m_fingerprint);
	
	if(fingerprint_col == COLUMN_NOT_FOUND){
		throw __FILE__ ":parse_fingerprint: Did not find the requested fingerprint column";
	}
	
	size_t fingerprint_len = 0;
	
	while( getline(fin, line) ){
		
		const vector<string> data = split(line, m_delim);
		
		if(data.size() != num_col){
			throw __FILE__ ":parse_fingerprint: Did not read the expected number of columns";
		}
		
		// Check for duplicate drug names
		const string &drug = data[drug_col];
		const string &fingerprint = data[fingerprint_col];
		
		if( ret.find(drug) != ret.end() ){
			
			cerr << "Drug id \"" << drug << "\" is a duplicate!" << endl;
			throw __FILE__ ":parse_fingerprint: Found a duplicate drug fingerprint entry!";
		}
		
		if( ret.empty() ){
			fingerprint_len = fingerprint.size();
		}
		
		if( fingerprint_len != fingerprint.size() ){
			throw __FILE__ ":parse_fingerprint: Unequal fingerprint lengths";
		}
		
		vector<float> local(fingerprint_len);
		
		for(size_t i = 0;i < fingerprint_len;++i){
				
			switch(fingerprint[i]){
				case '1':
					local[i] = 1.0;
					break;
				case '0':
					local[i] = 0.0;
					break;
				default:
					throw __FILE__ ":parse_fingerprint: Unknown fingerprint symbol";
			};				
		}
		
		ret[drug] = local;
	}

	// Since molecular fingerprint files do not typically provide the definition of
	// each bit, we will create dummy feature ids for each bit position in the fingerprint
	m_feature_id.resize(fingerprint_len);
	
	for(size_t i = 0;i < fingerprint_len;++i){
		
		stringstream ssin;
		
		ssin << "Bit" << i;
		
		m_feature_id[i] = ssin.str();
	}
	
	return ret;
}

bool is_fingerprint(const string &m_filename)
{
	ifstream fin( m_filename.c_str() );
	
	if(!fin){
		throw __FILE__ ":is_fingerprint: Unable to open input file for reading";
	}
	
	string line;
	
	// Read the header
	if( !getline(fin, line) ){
		throw __FILE__ ":is_fingerprint: Unable to read header";
	}
	
	const vector<string> header = split(line, ',');
	
	const size_t num_col = header.size();
	
	// Fingerprint files only have two columns (DRUG,FINGERPRINT)
	if(num_col > 2){
		return false;
	}
	
	if(num_col == 2){
		return true;
	}
	
	throw __FILE__ ":is_fingerprint: Too few columns in input file";
	return false;
		
}

MAP< string, vector<float> > parse_descriptor(const string &m_filename,
	vector<string> &m_feature_id,
	const string &m_drug,
	const char m_delim /*= ','*/)
{
	MAP< string, vector<float> > ret;
	MAP< string, vector<size_t> > norm; // Per-feature normalization
	
	ifstream fin( m_filename.c_str() );
	
	if(!fin){
		throw __FILE__ ":parse_descriptor: Unable to open input file for reading";
	}
	
	string line;
	
	// Read the header
	if( !getline(fin, line) ){
		throw __FILE__ ":parse_descriptor: Unable to read header";
	}
	
	const vector<string> header = split(line, m_delim);
	
	const size_t num_col = header.size();
	
	const size_t drug_col = get_col(header, m_drug);
	
	if(drug_col == COLUMN_NOT_FOUND){
		throw __FILE__ ":parse_descriptor: Did not find the requested drug name column";
	}
	
	// Store the feature names
	m_feature_id.resize(num_col - 1);
	size_t feature_index = 0;
	
	for(size_t i = 0;i < num_col;++i){
		
		if(i != drug_col){
		
			m_feature_id[feature_index] = header[i];
			++feature_index;
		}
	}
	
	while( getline(fin, line) ){
		
		const vector<string> data = split(line, m_delim);
		
		if(data.size() != num_col){
			throw __FILE__ ":parse_descriptor: Did not read the expected number of columns";
		}
		
		// Check for duplicate drug names
		const string &drug = data[drug_col];
		
		vector<float> local;
		
		local.reserve(num_col - 1);
		
		for(size_t i = 0;i < num_col;++i){
				
			if(i != drug_col){
				if( data[i].empty() ){
					local.push_back(MISSING_DATA);
				}
				else{
					local.push_back( string_to_float(data[i]) );
				}
			}				
		}
		
		const size_t num_features = local.size();
		
		// Per-feature normalization for handling missing data
		vector<size_t> local_norm(num_features);

		for(size_t i = 0;i != num_features;++i){
			local_norm[i] = (local[i] == MISSING_DATA) ? 0 : 1;
		}

		// Average duplicate feature vectors (since we are not currently using
		// batch-specific feature vectors
		
		MAP< string, vector<float> >::iterator iter = ret.find(drug);
		
		if( iter == ret.end() ){
		
			ret[drug] = local;
			norm[drug] = local_norm;
		}
		else{
			//cerr << "Duplicate drug name: " << drug << endl;
			//throw __FILE__ ":parse_descriptor: Found a duplicate drug feature vector entry!";
			
			if( num_features != iter->second.size() ){
				throw __FILE__ ":parse_descriptor: Unequal feature vector sizes";
			}
			
			MAP< string, vector<size_t> >::iterator norm_iter = norm.find(drug);
			
			if( norm_iter == norm.end() ){
				throw __FILE__ ":parse_descriptor: Unable to lookup normalization vector";
			}

			for(size_t i = 0;i < num_features;++i){
				
				if(local[i] != MISSING_DATA){
					
					if(iter->second[i] == MISSING_DATA){
						
						// Replace, don't average, missing data
						iter->second[i] = local[i];
						norm_iter->second[i] = local_norm[i];
					}
					else{
						iter->second[i] += local[i];
						norm_iter->second[i] += local_norm[i];
					}
				}
			}
		}
	}

	for(MAP< string, vector<float> >::iterator i = ret.begin();i != ret.end();++i){
		
		MAP< string, vector<size_t> >::iterator norm_iter = norm.find(i->first);

		if( norm_iter == norm.end() ){
			throw __FILE__ ":parse_descriptor: Unable to lookup normalization vector (2)";
		}
		
		const size_t num_features = i->second.size();
		
		if( num_features != norm_iter->second.size() ){
			throw __FILE__ ":parse_descriptor: Normalization vector size mismatch";
		}
		
		for(size_t j = 0;j < num_features;++j){
			
			if(i->second[j] != MISSING_DATA){
				
				if(norm_iter->second[j] == 0){
					throw __FILE__ ":parse_descriptor: Normalization divide by zero error!";
				}
				
				i->second[j] /= norm_iter->second[j];
			}
		}
	}
	
	return ret;
}

MAP<string, string> parse_label(const string &m_filename,
	const string &m_drug, const string &m_label)
{
	MAP<string, string> ret;
	
	ifstream fin( m_filename.c_str() );
	
	if(!fin){
		throw __FILE__ ":parse_label: Unable to open input file for reading";
	}
	
	string line;
	
	// Read the header
	if( !getline(fin, line) ){
		throw __FILE__ ":parse_label: Unable to read header";
	}
	
	const vector<string> header = split(line, ',');
	
	const size_t num_col = header.size();
	
	const size_t drug_col = get_col(header, m_drug);
	
	if(drug_col == COLUMN_NOT_FOUND){
		throw __FILE__ ":parse_label: Did not find the requested drug name column";
	}
	
	const size_t label_col = get_col(header, m_label);
	
	if(label_col == COLUMN_NOT_FOUND){
		throw __FILE__ ":parse_label: Did not find the requested label column";
	}
	
	while( getline(fin, line) ){
		
		const vector<string> data = split(line, ',');
		
		if(data.size() != num_col){
			throw __FILE__ ":parse_label: Did not read the expected number of columns";
		}
		
		// Check for duplicate drug names
		const string &drug = data[drug_col];
		const string &label = data[label_col];
		
		if( ret.find(drug) != ret.end() ){
			throw __FILE__ ":parse_label: Found a duplicate drug label entry!";
		}
		
		ret[drug] = label;
	}
	
	return ret;
}

MAP<string, float> parse_response(const string &m_filename,
	const string &m_drug, const string &m_response)
{
	MAP<string, float> ret;
	
	ifstream fin( m_filename.c_str() );
	
	if(!fin){
		throw __FILE__ ":parse_response: Unable to open input file for reading";
	}
	
	string line;
	
	// Read the header
	if( !getline(fin, line) ){
		throw __FILE__ ":parse_response: Unable to read header";
	}
	
	const vector<string> header = split(line, ',');
	
	const size_t num_col = header.size();
	
	const size_t drug_col = get_col(header, m_drug);
	
	if(drug_col == COLUMN_NOT_FOUND){
		throw __FILE__ ":parse_response: Did not find the requested drug name column";
	}
	
	const size_t response_col = get_col(header, m_response);
	
	if(response_col == COLUMN_NOT_FOUND){
		throw __FILE__ ":parse_response: Did not find the requested response column";
	}
	
	while( getline(fin, line) ){
		
		const vector<string> data = split(line, ',');
		
		if(data.size() != num_col){
			throw __FILE__ ":parse_response: Did not read the expected number of columns";
		}
		
		// Check for duplicate drug names
		const string &drug = data[drug_col];
		const string &response = data[response_col];
		
		if( ret.find(drug) != ret.end() ){
			throw __FILE__ ":parse_response: Found a duplicate drug label entry!";
		}
		
		ret[drug] = string_to_float(response);
	}
	
	return ret;
}

MAP< string, MAP<string, float> > parse_matrix(const string &m_filename)
{
	MAP< string, MAP<string, float> > ret;
	
	ifstream fin( m_filename.c_str() );
	
	if(!fin){
	
		cerr << "Error reading matrix file: " << m_filename << endl;
		throw __FILE__ ":parse_matrix: Unable to open input file for reading";
	}
	
	string line;
	
	// Read the header
	if( !getline(fin, line) ){
		
		cerr << "Error reading matrix file: " << m_filename << endl;
		throw __FILE__ ":parse_matrix: Unable to read header";
	}
	
	const vector<string> header = split(line, ',');
	
	const size_t num_col = header.size();
	
	while( getline(fin, line) ){
		
		const vector<string> data = split(line, ',');
		
		if(data.size() != num_col){
			throw __FILE__ ":parse_matrix: Did not read the expected number of columns";
		}
		
		// The row ID
		const string row_id = data[0];
		
		for(size_t i = 1;i < num_col;++i){
			
			// Check for duplicate entries
			if( ret.find(row_id) != ret.end() ){
				if( ret[row_id].find(header[i]) != ret[row_id].end() ){
					throw __FILE__ ":parse_matrix: Duplicate entry found!";
				}
			}
			
			if( is_missing(data[i]) ){
				ret[row_id][ header[i] ] = MISSING_DATA;
			}
			else{
				const float value = atof( data[i].c_str() );
				
				// A little sanity checking ...
				if(fabs(value) > 1000.0){
				
					cerr << "Suspicious value: " << row_id << '\t' << header[i] << '\t' << value << endl;
					throw __FILE__ ":parse_matrix: Found a questionable synergy value!";
				}
				
				ret[row_id][ header[i] ] = value;
			}
		}
	}
	
	return ret;
}

bool is_missing(const string &m_str)
{
	for(string::const_iterator i = m_str.begin();i != m_str.end();++i){
		if( !isspace(*i) ){
			return false;
		}
	}
	
	return true;
}

MAP<string, string> parse_drug_target(const string &m_filename)
{
	MAP<string, string> ret;
	
	if( m_filename.empty() ){
		return ret;
	}
	
	ifstream fin( m_filename.c_str() );
	
	if(!fin){
		throw __FILE__ ":parse_drug_target: Unable to open input file for reading";
	}
	
	string line;
	
	// Read the header
	if( !getline(fin, line) ){
		throw __FILE__ ":parse_drug_target: Unable to read header";
	}
	
	const vector<string> header = split(line, '\t');
	
	const size_t num_col = header.size();
	
	const size_t drug_col = get_col(header, "DRUG");
	
	if(drug_col == COLUMN_NOT_FOUND){
		throw __FILE__ ":parse_drug_target: Did not find the requested drug name column";
	}
	
	const size_t target_col = get_col(header, "TARGET");
	
	if(target_col == COLUMN_NOT_FOUND){
		throw __FILE__ ":parse_drug_target: Did not find the requested target column";
	}
	
	while( getline(fin, line) ){
		
		const vector<string> data = split(line, '\t');
		
		if(data.size() != num_col){
			throw __FILE__ ":parse_drug_target: Did not read the expected number of columns";
		}
		
		const string& drug = data[drug_col];
		const string& target = data[target_col];
		
		if( ret.find(drug) != ret.end() ){
			throw __FILE__ ":parse_drug_target: Duplicate drug target";
		}
		
		ret[drug] = target;
	}
	return ret;
}
