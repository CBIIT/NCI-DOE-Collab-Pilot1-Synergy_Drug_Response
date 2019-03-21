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
#include <unordered_set>
#include <string.h>
#include <math.h>

#include "gemini.h"
#include "zlib.h"
#include "parse_csv.h"

using namespace std;

struct sort_by_plate_and_num_drug
{
    inline bool operator()(const GrowthRecord &m_a, const GrowthRecord &m_b) const{
        
        if(m_a.plate == m_b.plate){
            
            return m_a.num_drug() < m_b.num_drug();
        }
        
        return m_a.plate < m_b.plate;
    };
};

char guess_delimeter(const string &m_filename);

bool parse_combo(const string &m_filename, deque<GrowthRecord> &m_data,
                 const int &m_lab_mask)
{
	// Since the ALMANAC datafile went from tab-delimited to comma delimited, we need to guess the
	// delimiter before parsing.
	const char delim = guess_delimeter(m_filename);
	
	// Allow reading (gzip) compressed or uncompressed files
	gzFile fin = gzopen(m_filename.c_str(), "r");

	if(fin == NULL){
	
		cerr <<":parse_combo: Unable to open input file for reading" << endl;
		return false;
	}

	const int buffer_len = 2048;
	char buffer[buffer_len];
	size_t line_number = 0;
	
	size_t num_col = 0;
	size_t id_col = 0;
	size_t screener_col = 0;
	size_t study_col = 0;
	size_t testdate_col = 0;
	size_t plate_col = 0;
	size_t nsc1_col = 0;
	size_t conc_index1_col = 0;
	size_t conc1_col = 0;
	size_t concunit1_col = 0;
	size_t nsc2_col = 0;
	size_t conc_index2_col = 0;
	size_t conc2_col = 0;
	size_t concunit2_col = 0;
	size_t percentgrowth_col = 0;
	//size_t percentgrowth_notz_col = 0;
	size_t testvalue_col = 0;
	size_t controlvalue_col = 0;
	size_t tzvalue_col = 0;
	size_t expected_growth_col = 0;
	size_t valid_col = 0;
	size_t cellname_col = 0;
	
	cerr << "Parsing growth data from: " << m_filename << " ... ";
	
	size_t num_invalid = 0;
	
	MAP<string, float> bliss_error;
	MAP<string, float> bliss_error2;
	MAP<string, size_t> bliss_count;
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		++line_number;
		
		const size_t actual_len = strlen(buffer);
		
		if( actual_len == (buffer_len - 1) ){
		
			cerr << __FILE__ ":parse_combo: Potential buffer overflow!" << endl;
			return false;
		}
		
		if(actual_len == 0){
			continue;
		}

		// Replace the terminal '\n' with a '\0'
		if(buffer[actual_len - 1] == '\n'){
			buffer[actual_len - 1] = '\0';
		}
		
		if( (actual_len == 1) && (buffer[0] == 26) ){
			
			cerr << "Skipping terminal SUB (0x1A) character" << endl;
			continue;
		}
		
		if( num_col == 0 ){

			// Read the header
			const vector<string> header = split(buffer, delim);

			num_col = header.size();
			
			id_col = get_col(header, "COMBODRUGSEQ");
			
			if(id_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"COMBODRUGSEQ\" in the header";
			}
			
			screener_col = get_col(header, "SCREENER");
			
			if(screener_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"SCREENER\" in the header";
			}
			
			study_col = get_col(header, "STUDY");
			
			if(study_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"STUDY\" in the header";
			}
			
			testdate_col = get_col(header, "TESTDATE");
			
			if(testdate_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"TESTDATE\" in the header";
			}
			
			plate_col = get_col(header, "PLATE");
			
			if(plate_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"PLATE\" in the header";
			}
			
			nsc1_col = get_col(header, "NSC1");
			
			if(nsc1_col == COLUMN_NOT_FOUND){
				
				// If we are unable to find a NSC1, look for CID1 before giving up
				nsc1_col = get_col(header, "CID1");
				
				if(nsc1_col == COLUMN_NOT_FOUND){
					throw __FILE__ ":parse_combo: Unable to find \"NSC1\" or \"CID1\" in the header";
				}
			}
			
			conc_index1_col = get_col(header, "CONCINDEX1");
			
			if(conc_index1_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"CONCINDEX1\" in the header";
			}
			
			conc1_col = get_col(header, "CONC1");
			
			if(conc1_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"CONC1\" in the header";
			}
			
			concunit1_col = get_col(header, "CONCUNIT1");
			
			if(concunit1_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"CONCUNIT1\" in the header";
			}
			
			nsc2_col = get_col(header, "NSC2");
			
			if(nsc2_col == COLUMN_NOT_FOUND){
				
				// If we are unable to find a NSC2, look for CID2 before giving up
				nsc2_col = get_col(header, "CID2");
				
				if(nsc2_col == COLUMN_NOT_FOUND){
					throw __FILE__ ":parse_combo: Unable to find \"NSC2\" or \"CID2\" in the header";
				}
			}
			
			conc_index2_col = get_col(header, "CONCINDEX2");
			
			if(conc_index2_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"CONCINDEX2\" in the header";
			}
			
			conc2_col = get_col(header, "CONC2");
			
			if(conc2_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"CONC2\" in the header";
			}
			
			concunit2_col = get_col(header, "CONCUNIT2");
			
			if(concunit2_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"CONCUNIT2\" in the header";
			}
			
			percentgrowth_col = get_col(header, "PERCENTGROWTH");
			
			if(percentgrowth_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"PERCENTGROWTH\" in the header";
			}
			
			//percentgrowth_notz_col = get_col(header, "PERCENTGROWTHNOTZ");
			testvalue_col = get_col(header, "TESTVALUE");
			
			if(testvalue_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"TESTVALUE\" in the header";
			}
			
			controlvalue_col = get_col(header, "CONTROLVALUE");
			
			if(controlvalue_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"CONTROLVALUE\" in the header";
			}
			
			tzvalue_col = get_col(header, "TZVALUE");
			
			if(tzvalue_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"TZVALUE\" in the header";
			}
			
			expected_growth_col = get_col(header, "EXPECTEDGROWTH");
			
			if(expected_growth_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"EXPECTEDGROWTH\" in the header";
			}
			
			valid_col = get_col(header, "VALID");
			
			if(valid_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"VALID\" in the header";
			}
			
			cellname_col = get_col(header, "CELLNAME");
			
			if(cellname_col == COLUMN_NOT_FOUND){
				throw __FILE__ ":parse_combo: Unable to find \"CELLNAME\" in the header";
			}
			
		}
		else{

			const vector<string> line = split(buffer, delim);

			if(line.size() != num_col){

				cerr << "Did not read the correct number of columns @ line " 
					<< line_number << endl;
				cerr << "actual_len = " << actual_len << endl;
				cerr << "buffer[0] = " << int(buffer[0]) << endl;
				cerr << "r = " << int('\r') << endl;
				cerr << "n = " << int('\n') << endl;
				cerr << "b = " << int('\b') << endl;
				cerr << buffer << endl;
				cerr << __FILE__ ":parse_combo: IO error" << endl;
				return false;
			}
			
			// There can be multiple headers in the combo data file
			if(line[nsc1_col] == "NSC1"){
			
				// This is another header line. Skip it.
				continue;
			}
			
			if(line[valid_col] != "Y"){
			
				//cerr << "Found a combo record marked as invalid on line " 
				//	<< line_number << endl;
				++num_invalid;
				continue;
			}

			GrowthRecord rec;

			rec.record_id = line[id_col];
			rec.parse_screener(line[screener_col]);
			rec.study = line[study_col];
			rec.testdate = Date(line[testdate_col]);
			rec.plate = line[plate_col];

			rec.nsc_1 = line[nsc1_col];
			rec.conc_index_1 = atoi( line[conc_index1_col].c_str() );
			rec.conc_1 = atof( line[conc1_col].c_str() );

			if(line[concunit1_col].size() != 1){
				throw __FILE__ ":parse_combo: Unable to parse CONCUNIT1";
			}

			rec.conc_unit_1 = line[concunit1_col][0];

			rec.nsc_2 = line[nsc2_col];
			rec.conc_index_2 = atoi( line[conc_index2_col].c_str() );
			rec.conc_2 = atof( line[conc2_col].c_str() );

			if( !rec.nsc_2.empty() && (line[concunit2_col].size() != 1) ){
				throw __FILE__ ":parse_combo: Unable to parse CONCUNIT2";
			}

			rec.conc_unit_2 = line[concunit2_col][0];
			
			rec.__percent_growth = atof( line[percentgrowth_col].c_str() );
			rec.test_value = atof( line[testvalue_col].c_str() );
			rec.control_value = atof( line[controlvalue_col].c_str() );
			rec.tz_value = atof( line[tzvalue_col].c_str() );
			rec.cellname = line[cellname_col];
            
			// Exclude records that do not belong to the requested laboratories
			if(rec.screener & m_lab_mask){
				m_data.push_back(rec);
			}
			
			// The expected growth (from the Bliss model of independent 
			// drug activity) is effectively a prediction of drug growth.
			// Compute the screening site-dependent mean squared error for
			// this prediction.
			if( false == line[expected_growth_col].empty() ){
				
				float delta_growth = rec.__percent_growth - 
					atof( line[expected_growth_col].c_str() );
				
				// Compute the mean *squared* error
				delta_growth *= delta_growth;
				
				bliss_error[ line[screener_col] ] += delta_growth;
				bliss_error2[ line[screener_col] ] += delta_growth*delta_growth;
				++bliss_count[ line[screener_col] ];
			}
		}
	}
	
	gzclose(fin);
	
	cerr << "done" << endl;

	cerr << "Found " << m_data.size() << " valid growth records" << endl;
    	cerr << "(Found " << num_invalid << " invalid growth records)" << endl;
	
	const deque<string> sites = keys(bliss_error);
	
	float total_bliss_error = 0.0;
	float total_bliss_error2 = 0.0;
	size_t total_bliss_count = 0;
	
	cerr << "Bliss model MSE for *all* concentrations (per screening site):" << endl;
	
	for(deque<string>::const_iterator i = sites.begin();i != sites.end();++i){
		
		MAP<string, float>::const_iterator error_iter = bliss_error.find(*i);
		
		if( error_iter == bliss_error.end() ){
			
			// Skip screening sites with no data
			continue;
		}
		
		MAP<string, float>::const_iterator error2_iter = bliss_error2.find(*i);
		
		if( error2_iter == bliss_error2.end() ){
			
			throw __FILE__ ":parse_combo: Unable to look up bliss error 2nd moment";
		}
		
		MAP<string, size_t>::const_iterator count_iter = bliss_count.find(*i);
		
		if( count_iter == bliss_count.end() ){
			throw __FILE__ ":parse_combo: Unable to look up bliss count";
		}
		
		if(count_iter->second > 0){
			
			total_bliss_error += error_iter->second;
			total_bliss_error2 += error2_iter->second;
			total_bliss_count += count_iter->second;
			
			const float mse = error_iter->second/count_iter->second;
			const float stdev = 
				sqrt(error2_iter->second/count_iter->second - mse*mse);
			
			cerr << '\t' << *i << '\t' << mse << " +/- " << stdev << endl;
		}
	}

	if(total_bliss_count > 0){
		
		total_bliss_error /= total_bliss_count;
		total_bliss_error2 = 
			sqrt(total_bliss_error2/total_bliss_count - 
				total_bliss_error*total_bliss_error);
		
		cerr << "Total Bliss model MSE for *all* concentrations = " << total_bliss_error 
			<< " +/- " << total_bliss_error2 << endl;
	}

#ifdef NCI_ALMANAC_QC_RULES
	
	// Sort the records by "PLATE" and then by single vs paired drug
	// data (so that the single drug data is first).
	//sort( m_data.begin(), m_data.end(), sort_by_plate_and_num_drug() );

	// Apply the following filter rule from the Almanac paper:
	//      "All data for a drug pair/cell line were removed from the
	//      dataset if the growth percent increased by 50% or greater
	//      between any 2 adjacent doses of either drug."

	cerr << "Applying Almanac QC filter rules:" << endl;
	
	MULTIMAP<string /*drug pair + cell line*/, deque<GrowthRecord>::iterator> filter;
	const float NCI_Alamanac_fraction_increase = 0.5f;
	
	for(deque<GrowthRecord>::iterator i = m_data.begin();i != m_data.end();++i){
        
		if( i->is_pair() ){
			filter.insert( make_pair(i->nsc_1 + ":" + i->nsc_2 + ":" + i->cellname, i) );
		}
    	}
    
	const deque<string> filter_keys = keys(filter);

	// The number of records that passed and failed the filter
	size_t num_filter_pass = 0;
	size_t num_filter_fail = 0;
	
	for(deque<string>::const_iterator i = filter_keys.begin();i != filter_keys.end();++i){
		
		// To enable easy testing, build a 2D matrix of records, ordered by the concentration
		// for each drug
		
		int max_index_1 = 0;
		int max_index_2 = 0;
		int min_index = 0;
		
		typedef MULTIMAP<string /*drug pair + cell line*/, deque<GrowthRecord>::iterator>::iterator I;

		pair<I, I> range = filter.equal_range(*i);

		// Find the matrix size
		for(I j = range.first;j != range.second;++j){

			max_index_1 = max(max_index_1, j->second->conc_index_1);
			max_index_2 = max(max_index_2, j->second->conc_index_2);
			
			min_index = min(min_index, min(j->second->conc_index_1, j->second->conc_index_2) );
		}
		
		if(min_index < 0){
			throw __FILE__ ":parse_combo: Found an index out of bounds!";
		}
		
		vector< vector<deque<GrowthRecord>::iterator> > m(max_index_1 + 1, 
        		vector<deque<GrowthRecord>::iterator>(max_index_2 + 1, m_data.end() ) );
		
		// Pack the matrix
		for(I j = range.first;j != range.second;++j){

			if( m[j->second->conc_index_1][j->second->conc_index_2] != m_data.end() ){
				
				cerr << "Filter key: " << *i << endl;
				throw __FILE__ ":parse_combo: Found an unexpected replicate";
			}
			
			m[j->second->conc_index_1][j->second->conc_index_2] = j->second;
		}
		
		bool qc_valid = true;
		
		// Test the *rows* of the matrix
		for(int j = 0;j < max_index_1;++j){
			
			for(int k = 1;k < max_index_2;++k){
				
				if( ( m[j][k - 1] == m_data.end() ) || ( m[j][k] == m_data.end() ) ){
					continue;
				}
				
				const float delta = m[j][k]->percent_growth() - m[j][k - 1]->percent_growth();
				
				if( delta > NCI_Alamanac_fraction_increase*fabs( m[j][k - 1]->percent_growth() ) ){
					qc_valid = false;
				}
			}
		}
		
		// Test the *columns* of the matrix
		for(int k = 0;k < max_index_2;++k){
			
			for(int j = 1;j < max_index_1;++j){
				
				if( ( m[j - 1][k] == m_data.end() ) || ( m[j][k] == m_data.end() ) ){
					continue;
				}
				
				const float delta = m[j][k]->percent_growth() - m[j - 1][k]->percent_growth();
				
				if( delta > NCI_Alamanac_fraction_increase*fabs( m[j - 1][k]->percent_growth() ) ){
					qc_valid = false;
				}
			}
		}
		
		// If there were any failures to pass the fileter, qc_valid will be false. When qc_valid is
		// false we must remove *all* of the dose measurements for this drug pair/cell line combination.
		if(qc_valid == false){
			
			for(I j = range.first;j != range.second;++j){
				
				j->second->passed_almanac_qc = false;
				++num_filter_fail;
			}
		}
		else{
			for(I j = range.first;j != range.second;++j){
				++num_filter_pass;
			}
		}
	}
    
    	cerr << "\t" << num_filter_pass << " drug pair records passed the QC process ("
		<< (100.0*num_filter_pass)/(num_filter_pass + num_filter_fail) << "%)" << endl;
	cerr << "\t" << num_filter_fail << " drug pair records failed the QC process ("
		<< (100.0*num_filter_fail)/(num_filter_pass + num_filter_fail) << "%)" << endl;

#endif // NCI_ALMANAC_QC_RULES
    
    	deque<GrowthRecord> clean_data;
    
    	for(deque<GrowthRecord>::const_iterator i = m_data.begin();i != m_data.end();++i){
        
        	// Only accept "clean", non-control measurements
        	if( i->passed_almanac_qc && !i->is_control() ){
            		clean_data.push_back(*i);
        	}
    	}

    	// Return the cleaned data
    	m_data.swap(clean_data);
    
    	cerr << "Imposing a constant drug 1 & drug 2 order ... ";
		
	size_t num_reordered = 0;

	for(deque<GrowthRecord>::iterator i = m_data.begin();i != m_data.end();++i){

		if( !i->is_pair() ){
			continue;
		}

		if(i->nsc_1 > i->nsc_2){

			i->swap_drugs();
			++num_reordered;
		}
	}

	cerr << "done (reordered " << num_reordered << " records)." << endl;

	return true;
}

char guess_delimeter(const string &m_filename)
{
	// Sample a small number of lines
	const size_t sample_size = 20;
	
	// Allow reading (gzip) compressed or uncompressed files
	gzFile fin = gzopen(m_filename.c_str(), "r");

	if(fin == NULL){
	
		cerr <<":guess_delimeter: Unable to open input file for reading" << endl;
		return false;
	}

	const int buffer_len = 2048;
	char buffer[buffer_len];
	size_t line_number = 0;
	
	MAP<char, size_t> delim_count;
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		++line_number;
		
		const size_t actual_len = strlen(buffer);
		
		if( actual_len == (buffer_len - 1) ){
		
			cerr << __FILE__ ":guess_delimeter: Potential buffer overflow!" << endl;
			return false;
		}
		
		if(actual_len == 0){
			continue;
		}

		for(size_t i = 0;i < actual_len;++i){
			
			switch(buffer[i]){
				case '\t':
					++delim_count['\t'];
					break;
				case ',':
					++delim_count[','];
					break;
			};
		}
		
		if(line_number > sample_size){
			break;
		}
	}
	
	gzclose(fin);
	
	char best_delim = '?';
	size_t best_count = 0;
	
	for(MAP<char, size_t>::const_iterator i = delim_count.begin();i != delim_count.end();++i){
		
		if(best_count < i->second){
			
			best_delim = i->first;
			best_count = i->second;
		}
	}
	
	if(best_count == 0){
		throw __FILE__ ":guess_delimeter: Unable to guess the record delimiter";
	}
	
	return best_delim;
}
