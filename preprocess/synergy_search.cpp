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

// Rank different synergy definitions and tissue types by likelihood of occuring by chance
// J. D. Gans
// Bioscience Division, B-11
// Los Alamos National Laboratory
// Wed Sep 13 08:58:19 2017

#include <iostream>
#include <fstream>
#include <unordered_set>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#include "gemini.h"
#include "parse_csv.h"
#include "shuffle.h"
#include "multistream.h"
#include "correlation.hpp"

using namespace std;

enum {
	TISSUE_NONE = 0,
	TISSUE_CNS = (1 << 0),
	TISSUE_ME =  (1 << 1),
	TISSUE_LC =  (1 << 2),
	TISSUE_LE =  (1 << 3),
	TISSUE_RE =  (1 << 4),
	TISSUE_PR =  (1 << 5),
	TISSUE_BR =  (1 << 6),
	TISSUE_OV =  (1 << 7),
	TISSUE_CO =  (1 << 8),
	TISSUE_ALL =  (1 << 9) - 1
};

// Allow up to 64 different tissue types
typedef size_t TissueType;

#define	USE_GRID_SEARCH				0xFFFFFFFF
#define	REQUIRED_FRACTION_OF_VALID_DATA		0.5

TissueType tissue_type(string m_cell /*copy*/);
size_t upper_diagonal_index(const size_t &m_i, const size_t &m_j, const size_t &m_N);
vector< pair<float, string> > compute_fractional_synergy(const vector<int> &m_synergy, 
	const deque<string> &m_drugs);
vector<float> compute_expected_fractional_synergy(const vector<int> &m_synergy, const deque<string> &m_drugs,
	gsl_rng *m_rand_gen);
void remove_prefix(string &m_name, const string &m_prefix);
	
int main(int argc, char *argv[])
{
	try{
		
		// Command line options:
		// -o <output file>
		// --dependent <csv file of drug pair values to predict>
		// [-n <number of cell lines needed for synergy>]
		// [--seed <random number seed>]
		// [-x <drug id to exclude, can be repeated>]
		// [-t <synergy threshold> (default is a grid search of potential thresholds)]
		// Tissue subsets (default is all)
		// [--CNS | --ME | --LC | --LE | --RE | --PR | --BR | --OV | --CO]
		// Cell line subsets (can be repeated)
		// [--include <celll line name>]
		// [--include.shared (include the ALMANAC/Merck shared cell lines)]
		
		const char* options = "o:n:x:t:?h";
		int config_opt = 0;
		int long_index = 0;

		struct option long_opts[] = {
			{"dependent", true, &config_opt, 1},
			{"seed", true, &config_opt, 2},
			{"CNS", false, &config_opt, 3},
			{"ME", false, &config_opt, 4},
			{"LC", false, &config_opt, 5},
			{"LE", false, &config_opt, 6},
			{"RE", false, &config_opt, 7},
			{"PR", false, &config_opt, 8},
			{"BR", false, &config_opt, 9},
			{"OV", false, &config_opt, 10},
			{"CO", false, &config_opt, 11},
			{"include", true, &config_opt, 12},
			{"include.shared", false, &config_opt, 13},
			{0,0,0,0} // Terminate options list
		};
		
		int opt_code;
		opterr = 0;
		
		bool print_usage = (argc == 1);
		string filename_pair_data;
		string filename_output;
		unsigned int seed = 0; // Default is time-based
		unsigned int num_cell_line_synergy_threshold = 1;
		TissueType tt = TISSUE_NONE;
		deque<string> exclude_drugs;
		float synergy_threshold = USE_GRID_SEARCH;
		unordered_set<string> cell_line_white_list;
		
		while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){
	
			switch( opt_code ){
				case 0:
					
					if(config_opt == 1){ // dependent
					
						filename_pair_data = optarg;
						break;
					}
					
					if(config_opt == 2){ // seed
					
						seed = string_to_uint(optarg);
						break;
					}
					
					if(config_opt == 3){ // CNS
					
						tt |= TISSUE_CNS;
						break;
					}
					
					if(config_opt == 4){ // ME
					
						tt |= TISSUE_ME;
						break;
					}
					
					if(config_opt == 5){ // LC
					
						tt |= TISSUE_LC;
						break;
					}
					
					if(config_opt == 6){ // LE
					
						tt |= TISSUE_LE;
						break;
					}
					
					if(config_opt == 7){ // RE
					
						tt |= TISSUE_RE;
						break;
					}
					
					if(config_opt == 8){ // PR
					
						tt |= TISSUE_PR;
						break;
					}
					
					if(config_opt == 9){ // BR
					
						tt |= TISSUE_BR;
						break;
					}
					
					if(config_opt == 10){ // OV
					
						tt |= TISSUE_OV;
						break;
					}
					
					if(config_opt == 11){ // CO
					
						tt |= TISSUE_CO;
						break;
					}
					
					if(config_opt == 12){ // include
					
						cell_line_white_list.insert(optarg);
						break;
					}
					
					if(config_opt == 13){ // include.shared
					
						cell_line_white_list.insert("HCT-116"); cell_line_white_list.insert("HCT116");
						cell_line_white_list.insert("NCI-H23"); cell_line_white_list.insert("NCIH23");
						cell_line_white_list.insert("NCI-H460"); cell_line_white_list.insert("NCIH460");
						cell_line_white_list.insert("OVCAR-3"); cell_line_white_list.insert("OVCAR3");
						cell_line_white_list.insert("SK-OV-3"); cell_line_white_list.insert("SKOV3");
						cell_line_white_list.insert("SW-620"); cell_line_white_list.insert("SW620");
						cell_line_white_list.insert("T-47D"); cell_line_white_list.insert("T47D");
						cell_line_white_list.insert("UACC-62"); cell_line_white_list.insert("UACC62");
						break;
					}
					
					cerr << "Unknown flag!" << endl;
					break;
				case 'o':
					filename_output = optarg;
					break;
				case 'n':
					num_cell_line_synergy_threshold = string_to_uint(optarg);
					break;
				case 'x':
					exclude_drugs.push_back(optarg);
					break;
				case  't':
					synergy_threshold = atof(optarg);
					break;
				case 'h':
				case '?':
					print_usage = true;
					break;
				default:
					cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
					break;
			};
		}
		
		if(tt == TISSUE_NONE){
			tt = TISSUE_ALL;
		}
		
		if(print_usage){
			
			cerr << "Usage:" << endl;
			cerr << "\t-o <output file>" << endl;
			cerr << "\t--dependent <csv file of drug pair values to predict>" << endl;
			cerr << "\t[-n <number of cell lines needed for synergy>] (default is 1)" << endl;
			cerr << "\t[--seed <random number seed>" << endl;
			cerr << "\t[-x <drug id to exclude, can be repeated>]" << endl;
			cerr << "\t[-t <synergy threshold> (default is a grid search of potential thresholds)]" << endl;
			cerr << "\tTissue subsets to include (default is all)" << endl;
			cerr << "\t[--CNS | --ME | --LC | --LE | --RE | --PR | --BR | --OV | --CO]" << endl;
			cerr << "\tCell line subsets (can be repeated)" << endl;
			cerr << "\t[--include <celll line name>]" << endl;
			cerr << "\t[--include.shared (include the ALMANAC/Merck shared cell lines)]" << endl;
			return EXIT_SUCCESS;
		}
		
		if( filename_output.empty() ){
		
			cerr << "Please specify an output file (-o)" << endl;
			return EXIT_FAILURE;
		}
		
		if( filename_pair_data.empty() ){
		
			cerr << "Please specify the dependent variable file (--dependent)" << endl;
			return EXIT_FAILURE;
		}
		
		if(num_cell_line_synergy_threshold == 0){
			
			cerr << "Please specify a cell line synergy threshold greater than zero (-n)" << endl;
			return EXIT_FAILURE;
		}
		
		if(seed == 0){
			
			// Use a time-based seed
			seed = time(NULL);
		}
		
		ofstream fout( filename_output.c_str() );
		
		if(!fout){
		
			cerr << "Unable to open output file \"" << filename_output << "\" for writing" << endl;
			return EXIT_FAILURE;
		}
		
		// Write the command line arguments
		//fout << "# Command line:";
		//
		//for(int i = 0;i < argc;++i){
		//	fout << ' ' << argv[i];
		//}
		//
		//fout << endl;
		
		mostream mout;
		
		// A single object to handle writing to multiple streams
		mout.push_back(&fout);
		mout.push_back(&cerr);

		gsl_rng *rand_gen = gsl_rng_alloc(gsl_rng_ranlux389);
		
		if(rand_gen == NULL){
			throw __FILE__ ":main: Unable to initialize random number generator";
		}

		gsl_rng_set(rand_gen, seed);
		
		if(synergy_threshold == USE_GRID_SEARCH){
			cerr << "Using a grid search to identify the best synergy threshold value" << endl;
		}
		else{
			cerr << "Using synergy threshold value = " << synergy_threshold << endl;
		}
		
		cerr << "Requiring synergy in " << num_cell_line_synergy_threshold 
			<< " cell lines to define an interaction as synergistic" << endl;
		
		cerr << "Including tissues:";
		
		if(tt & TISSUE_CNS){
			cerr << " CNS";
		}
		
		if(tt & TISSUE_ME){
			cerr << " ME";
		}
		
		if(tt & TISSUE_LC){
			cerr << " LC";
		}
		
		if(tt & TISSUE_LE){
			cerr << " LE";
		}
		
		if(tt & TISSUE_RE){
			cerr << " RE";
		}
		
		if(tt & TISSUE_PR){
			cerr << " PR";
		}
		
		if(tt & TISSUE_BR){
			cerr << " BR";
		}
		
		if(tt & TISSUE_OV){
			cerr << " OV";
		}
		
		if(tt & TISSUE_CO){
			cerr << " CO";
		}
		
		cerr << endl;
		
		
		if( !cell_line_white_list.empty() ){
			
			cerr << "Including the following cell lines:" << endl;
			
			for(unordered_set<string>::const_iterator i = cell_line_white_list.begin();i != cell_line_white_list.end();++i){
				
				cerr << '\t' << *i << endl;
			}
		}
		
		vector<string> cell_names;
		
		// Non-const variable so we can randomize for pair-based cross validation
		deque<DrugPairData> pair_drug = parse_pair_drug(filename_pair_data, cell_names,
                                                                1.0 /*normalization_factor*/);
		
		// There are two ALMANAC drugs that were not tested (mostly).
		// These are the corresponding NSC values
		exclude_drugs.push_back("707389"); // Eribulin
		exclude_drugs.push_back("761431"); //Vemurafenib
		
		// These are the corresponding CID values
		exclude_drugs.push_back("54611489"); // Eribulin sulfate
		exclude_drugs.push_back("11354606"); // Eribulin
		
		// Don't exclude the CID for Vemurafenib, sicne this drug appears
		// *twice* in the ALMANAC dataset. The other version, NSC 75308,
		// has complete data (and will be represented by CID 42611257).
		//exclude_drugs.push_back("42611257"); //Vemurafenib
		
		make_set(exclude_drugs);
		
		deque<DrugPairData> tmp;
		
		// Remove the excluded drugs
		for(deque<DrugPairData>::const_iterator i = pair_drug.begin();i != pair_drug.end();++i){
			
			if( set_contains(exclude_drugs, i->drug_id_1) || set_contains(exclude_drugs, i->drug_id_2) ){
				continue;
			}
			
			tmp.push_back(*i);
		}
		
		pair_drug = tmp;
		
		tmp.clear();
		
		deque<string> drugs;
		
		for(deque<DrugPairData>::const_iterator i = pair_drug.begin();i != pair_drug.end();++i){
			
			drugs.push_back(i->drug_id_1);
			drugs.push_back(i->drug_id_2);
		}
		
		make_set(drugs);
		
		const size_t N = drugs.size();
		
		cerr << "Found a total of " << N << " valid drugs" << endl;
		
		float min_synergy = 1.0e25;
		float max_synergy = -1.0e25;
		float delta_synergy = 0.0;
		unsigned int num_step = 100;
		
		unordered_set<string> matching_cells;

		// Compute the range of scores
		for(deque<DrugPairData>::const_iterator i = pair_drug.begin();i != pair_drug.end();++i){

		    for(vector<float>::const_iterator j = i->growth.begin();j != i->growth.end();++j){

			// Select by cell type
			if( (tissue_type(cell_names[j - i->growth.begin()]) & tt) == 0 ){
				continue;
			}

			// Skip all of the cells that are *not* white-listed
			if( !cell_line_white_list.empty() ){
			
				if( cell_line_white_list.find(cell_names[j - i->growth.begin()]) == 
				    cell_line_white_list.end() ){
					continue;
				}
			}
			
			if(*j != MISSING_DATA){

				min_synergy = min(min_synergy, *j);
				max_synergy = max(max_synergy, *j);

				matching_cells.insert(cell_names[j - i->growth.begin()]);
			}
		    }
		}

		cerr << "Synergy score range: " << min_synergy << " <= synergy <= " << max_synergy << endl;

		cerr << "Found " << matching_cells.size() << " matching cell lines" << endl;

		for(unordered_set<string>::const_iterator i = matching_cells.begin();i != matching_cells.end();++i){
			cerr << '\t' << *i << endl;
		}

		if( !cell_line_white_list.empty() && ( matching_cells.size() != cell_line_white_list.size() ) ){
		
			cerr << "** Warning ** Did not find all of the requested cell line names" << endl;
			cerr << "Missing cell line names:" << endl;
			
			for(unordered_set<string>::const_iterator i = cell_line_white_list.begin();i != cell_line_white_list.end();++i){
				
				if( matching_cells.find(*i) == matching_cells.end() ){
					cerr << '\t' << *i << endl;
				}
			}
		}

		if( matching_cells.empty() ){
		
			cerr << "** Unable to find any matching cell lines **" << endl;
			return EXIT_FAILURE;
		}
			
		delta_synergy = (max_synergy - min_synergy)/num_step;
		
		if(synergy_threshold != USE_GRID_SEARCH){
		
			min_synergy = max_synergy = synergy_threshold;
			delta_synergy = 1.0; // Any value greater than zero will do
			num_step = 1;
		}
		
		vector< pair<float, string> > best_fractional_synergy;
		vector<float> best_synergy;
		float best_score = 0.0;
		float best_threshold = 0.0;
		
		// synergy threshold grid search
		for(float synergy_threshold = min_synergy;
			//synergy_threshold <= min(0.0f, max_synergy);
			synergy_threshold <= max_synergy;
			synergy_threshold += delta_synergy){
			
			vector<float> synergy( (N*(N - 1))/2, 0 /* 0 means "not tested" */);
			
			for(size_t i = 0;i < N;++i){

				for(size_t j = 0;j < N;++j){

					if(i == j){
						continue;
					}

					size_t number_synergetic = 0;
					bool tested = false;
					bool valid = false;
					float min_synergy = 1.0e6;
					
					for(deque<DrugPairData>::const_iterator k = pair_drug.begin();k != pair_drug.end();++k){

						if( ( (k->drug_id_1 == drugs[i]) && (k->drug_id_2 == drugs[j]) ) ||
						    ( (k->drug_id_2 == drugs[i]) && (k->drug_id_1 == drugs[j]) ) ){

						    tested = true;

						    for(vector<float>::const_iterator iter = k->growth.begin();iter != k->growth.end();++iter){

							// Select by cell type
							if( (tissue_type(cell_names[iter - k->growth.begin()]) & tt) == 0 ){
								continue;
							}
							
							// Skip all of the cells that are *not* white listed
							if( !cell_line_white_list.empty() ){
							
								if( cell_line_white_list.find(cell_names[iter - k->growth.begin()]) == 
								    cell_line_white_list.end() ){
									continue;
								}
							}
							
							if( (*iter != MISSING_DATA) && (*iter <= synergy_threshold) ){
								++number_synergetic;
							}
							
							if( (*iter != MISSING_DATA) && (*iter < min_synergy) ){
								
								valid = true;
								min_synergy = *iter;
							}

						    }
						}
					}
					
					#ifdef SYNERGY_BY_COUNTING
					const bool has_synergy = (number_synergetic >= num_cell_line_synergy_threshold);
					
					synergy[ upper_diagonal_index(i, j, N) ] = (tested ? (has_synergy ? 1 : -1) : 0);
					#else
					synergy[ upper_diagonal_index(i, j, N) ] = (tested && valid) ? min_synergy : MISSING_DATA;
					#endif // SYNERGY_BY_COUNTING
				}
			}
			
			#ifndef SYNERGY_BY_COUNTING
			best_synergy = synergy;
			
			cerr << "Breaking out of search loop to report minimum synergy values" << endl;
			break;
			#endif // SYNERGY_BY_COUNTING
			
			//const vector< pair<float, string> > fractional_synergy = 
			//	compute_fractional_synergy(synergy, drugs);
			const vector< pair<float, string> > fractional_synergy;
			
			//const vector<float> expected_synergy = 
			//	compute_expected_fractional_synergy(synergy, drugs, rand_gen);
			const vector<float> expected_synergy;
			
			float score = 0.0;
			
			float norm = 0.0;
			
			for(size_t i = 0;i < N;++i){
				
				// Don't include drugs with too much missing data
				if( (fractional_synergy[i].first >= 0.0) && (expected_synergy[i] >= 0.0) ){
					
					const float delta = fractional_synergy[i].first - expected_synergy[i];
				
					score += delta*delta;
					++norm;
				}
			}
			
			if(norm > 0){
				score = sqrt(score/norm);
			}
			
			if(score > best_score){
			
				best_score = score;
				best_fractional_synergy = fractional_synergy;
				best_threshold = synergy_threshold;
				best_synergy = synergy;
			}
			
			// DEBUG
			//mout << synergy_threshold << '\t' << score << endl;
			cerr << synergy_threshold << '\t' << score << endl;
		}
				
		cerr << "Best score is " << best_score << endl;
		cerr << "Best threshold is " << best_threshold << endl;
		
		#ifdef FRACTION_SYNERGISTIC
		fout << "DRUG,Fraction synergistic" << endl;
		
		// Print in descending (reverse) order
		for(vector< pair<float, string> >::const_reverse_iterator i = best_fractional_synergy.rbegin();
			i != best_fractional_synergy.rend();++i){

			fout << i->second << ',' << i->first << endl;
		}
		#endif // FRACTION_SYNERGISTIC
		
		fout << "DRUG";
		
		for(size_t i = 0;i < N;++i){
			fout << ',' << drugs[i];
		}
		
		fout << endl;
		
		// Print in descending (reverse) order
		for(size_t i = 0;i < N;++i){
			
			fout << drugs[i];
			
			for(size_t j = 0;j < N;++j){
				
				fout << ',';
				
				if(i == j){
					//fout << '0'; // Not tested
					//Missing data
				}
				else{
					if(best_synergy[ upper_diagonal_index(i, j, N) ] != MISSING_DATA){
						fout << best_synergy[ upper_diagonal_index(i, j, N) ];
					}
				}
			}
			
			fout << endl;
		}
		
		// Deallocate the random number generator
		gsl_rng_free(rand_gen);
	}
	catch(const char *error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(const string error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){
		cerr << "Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}

vector<float> compute_expected_fractional_synergy(const vector<int> &m_synergy, const deque<string> &m_drugs,
	gsl_rng *m_rand_gen)
{
	const size_t N = m_drugs.size();
	
	const size_t num_iter = 100;
	
	vector<float> ret(N);
	vector<size_t> norm(N);

	vector<int> synergy(m_synergy); // Make a copy
		
	for(size_t iter = 0;iter < num_iter;++iter){

		randomize(synergy.begin(), synergy.end(), m_rand_gen);

		const vector< pair<float, string> > f = compute_fractional_synergy(synergy, m_drugs);

		for(size_t i = 0;i < N;++i){
			
			if(f[i].first >= 0.0){
			
				ret[i] += f[i].first;
				++norm[i];
			}
		}
	}

	for(size_t i = 0;i < N;++i){
		
		if(norm[i] > 0){
			ret[i] /= norm[i];
		}
	}
		
	return ret;
}

vector< pair<float, string> > compute_fractional_synergy(const vector<int> &m_synergy, const deque<string> &m_drugs)
{
	const size_t N = m_drugs.size();
	
	vector< pair<float, string> > ret(N);
	
	for(size_t i = 0;i < N;++i){

		size_t num_synergy = 0;
		size_t num_no_synergy = 0;
		size_t num_no_test = 0;

		for(size_t j = 0;j < N;++j){

			if(i != j){
				switch(m_synergy[ upper_diagonal_index(i, j, N) ]){
					case -1:
						++num_no_synergy;
						break;
					case 0:
						++num_no_test;
						break;
					case 1:
						++num_synergy;
						break;
					default:
						throw __FILE__ ":compute_fractional_synergy: Unknown state!";
				};
			}
		}

		float f = -1.0;
		
		// Only count the fraction of synergistic interactions from the drugs that have
		// more than REQUIRED_FRACTION_OF_VALID_DATA fraction of valid measurements.
		// Some drug combination data sets have lots of missing data (such that there may
		// only be a small number of mostly complete drugs and a large number of 
		// mostly inconplete drugs).
		if( float(num_synergy + num_no_synergy)/N > REQUIRED_FRACTION_OF_VALID_DATA ){
			f = float(num_synergy)/(num_synergy + num_no_synergy);
		}
		
		ret[i] = make_pair(f, m_drugs[i]);
	}
	
	sort( ret.begin(), ret.end() );
	
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
	
TissueType tissue_type(string m_cell /*copy*/)
{
	// Remove dataset prefixes
	remove_prefix(m_cell, "CCLE.");
	remove_prefix(m_cell, "NCI60.");
	remove_prefix(m_cell, "CTRP.");
	remove_prefix(m_cell, "GDSC.");
	remove_prefix(m_cell, "gCSI.");
	remove_prefix(m_cell, "GDC.");
	remove_prefix(m_cell, "NCIPDM.");
	
	/////////////////////////////////////////////////
	// ALMANAC cell lines
	/////////////////////////////////////////////////
	if( (m_cell == "SF-539") || (m_cell == "SF539") ){
		return TISSUE_CNS;
	}
	
	if( (m_cell == "MDA-MB-435") || (m_cell == "MDAMB435S") ){
		return TISSUE_ME;
	}
	
	if(m_cell == "NCI-H322M"){
		return TISSUE_LC;
	}
	
	if( (m_cell == "SR") || (m_cell == "SR786") ){
		return TISSUE_LE;
	}
	
	if(m_cell == "SN12C"){
		return TISSUE_RE;
	}
	
	if( (m_cell == "PC-3") || (m_cell == "PC3") ){
		return TISSUE_PR;
	}
	
	if( (m_cell == "UO-31") || (m_cell == "UO31") ){
		return TISSUE_RE;
	}
	
	if( (m_cell == "HS 578T") || (m_cell == "HS578T") ){
		return TISSUE_BR;
	}

	if(m_cell == "TK-10"){
		return TISSUE_RE;
	}
	
	if(m_cell == "ACHN"){
		return TISSUE_RE;
	}
	
	if( (m_cell == "CAKI-1") || (m_cell == "CAKI1") ){
		return TISSUE_RE;
	}
	
	if( (m_cell == "OVCAR-4") || (m_cell == "OVCAR4") ){
		return TISSUE_OV;
	}
	
	if(m_cell == "RXF 393"){
		return TISSUE_RE;
	}
	
	if(m_cell == "CCRF-CEM"){
		return TISSUE_LE;
	}
	
	if(m_cell == "SK-MEL-2"){
		return TISSUE_ME;
	}
	
	if( (m_cell == "MDA-MB-468") || (m_cell == "MDAMB468") ){
		return TISSUE_BR;
	}
	
	if( (m_cell == "786-0") || (m_cell == "786O") ){
		return TISSUE_RE;
	}
	
	if(m_cell == "M14"){
		return TISSUE_ME;
	}
	
	if( (m_cell == "HOP-92") || (m_cell == "HOP92") ){
		return TISSUE_LC;
	}
	
	if(m_cell == "COLO 205"){
		return TISSUE_CO;
	}
	
	if(m_cell == "OVCAR-5"){
		return TISSUE_OV;
	}
	
	if( (m_cell == "SNB-75") || (m_cell == "SNB75") ){
		return TISSUE_CNS;
	}
	
	if(m_cell == "EKVX"){
		return TISSUE_LC;
	}
		
	if( (m_cell == "SK-MEL-28") || (m_cell == "SKMEL28") ){
		return TISSUE_ME;
	}
	
	if( (m_cell == "NCI-H226") || (m_cell == "NCIH226") ){
		return TISSUE_LC;
	}
	
	if(m_cell == "IGROV1"){
		return TISSUE_OV;
	}
	
	if( (m_cell == "HOP-62") || (m_cell == "HOP62") ){
		return TISSUE_LC;
	}
	
	if( (m_cell == "SF-295") || (m_cell == "SF295") ){
		return TISSUE_CNS;
	}
		
	if( (m_cell == "K-562") || (m_cell == "K562") ){
		return TISSUE_LE;
	}
	
	if( (m_cell == "SK-MEL-5") || (m_cell == "SKMEL5") ){
		return TISSUE_ME;
	}
	
	if(m_cell == "NCI/ADR-RES"){
		return TISSUE_OV;
	}
	
	if( (m_cell == "NCI-H522") || (m_cell == "NCIH522") ){
		return TISSUE_LC;
	}
	
	if( (m_cell == "A549/ATCC") || (m_cell == "A549") ){
		return TISSUE_LC;
	}
	
	if(m_cell == "SW-620"){
		return TISSUE_CO;
	}
	
	if(m_cell == "HCT-116"){
		return TISSUE_CO;
	}
	
	if( (m_cell == "OVCAR-8") || (m_cell == "OVCAR8") ){
		return TISSUE_OV;
	}
	
	if(m_cell == "OVCAR-3"){
		return TISSUE_OV;
	}
	
	if(m_cell == "SK-OV-3"){
		return TISSUE_OV;
	}
	
	if( (m_cell == "UACC-257") || (m_cell == "UACC257") ){
		return TISSUE_ME;
	}
	
	if(m_cell == "MCF7"){
		return TISSUE_BR;
	}
	
	if(m_cell == "HCC-2998"){
		return TISSUE_CO;
	}
	
	if( (m_cell == "RPMI-8226") || (m_cell == "RPMI8226") ){
		return TISSUE_LE;
	}
	
	if( (m_cell == "HL-60(TB)") || (m_cell == "HL60") ){
		return TISSUE_LE;
	}
	
	if(m_cell == "KM12"){
		return TISSUE_CO;
	}
	
	if( (m_cell == "MDA-MB-231/ATCC") || (m_cell == "MDAMB231") ){
		return TISSUE_BR;
	}
	
	if( (m_cell == "DU-145") || (m_cell == "DU145") ){
		return TISSUE_PR;
	}
	
	if(m_cell == "NCI-H460"){
		return TISSUE_LC;
	}
	
	if( (m_cell == "BT-549") || (m_cell == "BT549") ){
		return TISSUE_BR;
	}
	
	if( (m_cell == "T-47D") || (m_cell == "T47D") ){
		return TISSUE_BR;
	}
	
	if(m_cell == "SNB-19"){
		return TISSUE_CNS;
	}
	
	if( (m_cell == "HCT-15") || (m_cell == "HCT15") ){
		return TISSUE_CO;
	}
	
	if(m_cell == "HT29"){
		return TISSUE_CO;
	}
	
	if( (m_cell == "MALME-3M") || (m_cell == "MALME3M") ){
		return TISSUE_ME;
	}
	
	if(m_cell == "A498"){
		return TISSUE_RE;
	}
	
	if( (m_cell == "U251") || (m_cell == "U251MG") ){
		return TISSUE_CNS;
	}
	
	if(m_cell == "NCI-H23"){
		return TISSUE_LC;
	}
	
	if( (m_cell == "SF-268") || (m_cell == "SF268") ){
		return TISSUE_CNS;
	}
	
	if( (m_cell == "UACC-62") || (m_cell == "UACC62") ){
		return TISSUE_ME;
	}
	
	if( (m_cell == "LOX IMVI") || (m_cell == "LOXIMVI") ){
		return TISSUE_ME;
	}
	
	if(m_cell == "MOLT-4"){
		return TISSUE_LE;
	}
	
	//////////////////////////////////////////////////////////
	// Merck cell lines. Tissue type determined by searching: 
	//	Supplementary online table 1: http://mct.aacrjournals.org/content/15/6/1155.figures-only
	// 	COSMIC database: http://cancer.sanger.ac.uk
	//	Cellosaurus: https://web.expasy.org/cellosaurus
	//	Google	
	//////////////////////////////////////////////////////////
	
	if(m_cell == "A2058"){ // COSMIC
		return TISSUE_ME;
	}
	
	if(m_cell == "A2780"){ // COSMIC
		return TISSUE_OV;
	}
	
	if(m_cell == "A375"){ // COSMIC
		return TISSUE_ME;
	}
	
	if(m_cell == "A427"){ // COSMIC
		return TISSUE_LC;
	}
	
	if(m_cell == "CAOV3"){ // Google search
		return TISSUE_OV;
	}

	if(m_cell == "COLO320DM"){ // Cellosaurus
		return TISSUE_CO;
	}
	
	if(m_cell == "COLO320"){ // CCLE
		return TISSUE_CO;
	}
	
	if(m_cell == "DLD1"){ // Cellosaurus
		return TISSUE_CO;
	}

	if(m_cell == "EFM192B"){ // Cellosaurus
		return TISSUE_BR;
	}

	if(m_cell == "ES2"){ // Supplementary table 1
		return TISSUE_OV;
	}

	if(m_cell == "HCT116"){ // Supplementary table 1
		return TISSUE_CO;
	}

	if(m_cell == "HT144"){ // Supplementary table 1
		return TISSUE_ME;
	}

	if(m_cell == "HT29"){ // Supplementary table 1
		return TISSUE_CO;
	}

	if(m_cell == "KPL1"){ // Supplementary table 1
		return TISSUE_BR;
	}

	if(m_cell == "LNCAP"){ // Supplementary table 1
		return TISSUE_PR;
	}

	if(m_cell == "LOVO"){ // Supplementary table 1
		return TISSUE_CO;
	}

	if(m_cell == "MDAMB436"){ // Supplementary table 1
		return TISSUE_BR;
	}

	if(m_cell == "MSTO"){ // Supplementary table 1
		return TISSUE_LC;
	}

	if(m_cell == "NCIH1650"){ // Supplementary table 1
		return TISSUE_LC;
	}

	if(m_cell == "NCIH2122"){ // Supplementary table 1
		return TISSUE_LC;
	}

	if(m_cell == "NCIH23"){ // Supplementary table 1
		return TISSUE_LC;
	}

	if(m_cell == "NCIH460"){ // Cellosaurus
		return TISSUE_LC;
	}

	if(m_cell == "NCIH520"){ // Supplementary table 1
		return TISSUE_LC;
	}

	if(m_cell == "OCUBM"){ // Supplementary table 1
		return TISSUE_BR;
	}

	if(m_cell == "OV90"){ // Supplementary table 1
		return TISSUE_OV;
	}

	if( (m_cell == "OVCAR3") || (m_cell == "NIHOVCAR3") ){ // Supplementary table 1
		return TISSUE_OV;
	}

	if(m_cell == "PA1"){ // Supplementary table 1
		return TISSUE_OV;
	}

	if(m_cell == "RKO"){ // Supplementary table 1
		return TISSUE_CO;
	}

	if(m_cell == "RPMI7951"){ // Supplementary table 1
		return TISSUE_ME;
	}

	if(m_cell == "SKMEL30"){ // Supplementary table 1
		return TISSUE_ME;
	}

	if(m_cell == "SKMES1"){ // Supplementary table 1
		return TISSUE_LC;
	}

	if(m_cell == "SKOV3"){ // Supplementary table 1
		return TISSUE_OV;
	}

	if(m_cell == "SW620"){ // Supplementary table 1
		return TISSUE_CO;
	}

	if(m_cell == "SW837"){ // Supplementary table 1
		return TISSUE_CO;
	}

	if(m_cell == "UACC62"){ // Supplementary table 1
		return TISSUE_ME;
	}

	if(m_cell == "UWB1289"){ // Supplementary table 1
		return TISSUE_OV;
	}
	
	if(m_cell == "UWB1289BRCA1"){ // Google
		return TISSUE_OV;
	}

	if(m_cell == "VCAP"){ // Supplementary table 1
		return TISSUE_PR;
	}

	if(m_cell == "ZR751"){ // Supplementary table 1
		return TISSUE_BR;
	}

	if(m_cell == "LNCAPCLONEFGC"){ // Google -> Cellosaurus
		return TISSUE_PR;
	}
	
	if(m_cell == "MSTO211H"){ //Cellosaurus
		return TISSUE_LC;
	}
	
	if(m_cell == "EFM192A"){ // CCLE
		return TISSUE_BR;
	}
	
	cerr << "Unknown cell type: \"" << m_cell << "\"" << endl;
	
	throw __FILE__ ":tissue_type: Unable to find cell line name";
	return TISSUE_NONE;
}

size_t upper_diagonal_index(const size_t &m_i, const size_t &m_j, const size_t &m_N)
{
	if(m_i == m_j){
		throw __FILE__ ":upper_diagonal_index: Illegal diagonal element access";
	}

	if( (m_i >= m_N) || (m_j >= m_N) ){
		throw __FILE__ ":upper_diagonal_index: Index out of bounds!";
	}

	const size_t x = min(m_i, m_j);
	const size_t y = max(m_i, m_j);

	return ( x*(2*m_N - 1 - x) )/2 + y - x - 1;
}

