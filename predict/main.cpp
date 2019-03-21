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

// Regress drugs
// J. D. Gans
// Bioscience Division, B-11
// Los Alamos National Laboratory
// Fri Dec  1 10:18:27 2017

#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
#include <getopt.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <deque>
#include <vector>
#include <algorithm>
#include <unordered_set>

#include "rx.h"
#include "mpi_util.h"
#include "random_forest.h"
#include "shuffle.h"
#include "deque_set.h"

#include "correlation.hpp"
#include "keys.hpp"

using namespace std;

// Global variables for MPI
int mpi_numtasks;
int mpi_rank;
double start_time;

void terminate_program(int m_sig);
string report_run_time();
float compute_distance(const vector<float> &m_a, const vector<float> &m_b);
pair<float /*value*/, float /*weight*/> synergy_score(const vector<string> &m_drugs, 
	const MAP<string, float> &m_synergy, const ScoreFunction &m_func,
	const float &m_synergy_threshold);

MAP< string, MAP<string, float> > NSC_to_CID(const MAP< string, MAP<string, float> > &m_synergy_matrix);
string NSC_to_CID(const string &m_nsc);

void randomize_independent(MAP<string, vector<float> > &m_data, unsigned int *m_ptr_seed);
void randomize_dependent(MAP<string, MAP< string, MAP<string, float> > > &m_data, unsigned int *m_ptr_seed);
vector<bool> find_invariant_features(const MAP<string, vector<float> > &m_drug_features, 
	const MAP<string, vector<float> > &m_blind_test_drug_features);
size_t count(const vector<bool> &m_x);
template<class T>
void remove_invariant_features(vector<T> &m_x, const vector<bool> &m_invariant_mask);
void remove_drug(MAP< string, MAP<string, float> > &m_matrix, const string &m_drug);

pair<float, float> average_and_stdev(const deque<float> &m_data);
double compute_enrichment(const deque< pair<float, bool> > &m_synergy_data);
double compute_AUROC(const deque< pair<float, bool> > &m_synergy_data);
double compute_gini(const double &m_auroc);
void normalize(MAP<string, vector<float> > &m_features);
double compute_p_value(const float &m_x, const deque<float> &m_data);
float fraction_synergy(const MAP< string, MAP< string, MAP<string, float> > > &m_synergy_matrix,
	const float &m_synergy_threshold);
string probability_heat_map(const MAP< string, MAP< string, MAP<string, float> > > &m_synergy_matrix,
	const float &m_synergy_threshold);
	
vector<string> set_difference(const vector<string> &m_a, const vector<string> &m_b);

int main(int argc, char *argv[])
{
	try{
		
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_numtasks);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		
		start_time = MPI_Wtime();
		
		signal( SIGINT, terminate_program );
		signal( SIGTERM, terminate_program );
		signal( SIGSEGV, terminate_program );
		
		Options opt;
		
		if(mpi_rank == 0){
			opt.parse(argc, argv);
		}
		
		broadcast(opt, mpi_rank, 0);
		
		if(opt.print_usage){
			
			MPI_Finalize();
			return EXIT_SUCCESS;
		}

		if(mpi_rank != 0){
		
			// Redirect the mpi_rank != 0 output to /dev/null
			opt.file_name_output = "/dev/null";
			
			// Each rank gets a unique seed
			opt.seed += mpi_rank;
		}
		
		ofstream fout;
		
		if( !opt.file_name_output.empty() ){
			
			fout.open( opt.file_name_output.c_str() );
			
			if(!fout){
				throw __FILE__ ":main: Unable to open output file for writing";
			}
		}
		
		ostream &out = ( fout.is_open() ) ? fout : cout;
		
		out << "Drug Synergy Prediction Framework version " << RX_VERSION << endl;
		out << "Running with " << mpi_numtasks << " MPI ranks" << endl;
		
		out << "Input parameters:" << endl;
		
		out << "\tDrug feature filename = " << opt.file_name_drug << endl;
		out << "\tCell line feature filename = " << opt.file_name_cell << endl;
		
		if(opt.randomize_drug){
			out << "\t** Randomizing drug features **" << endl;
		}
		
		if(opt.randomize_cell){
			out << "\t** Randomizing cell line features **" << endl;
		}
		
		out << "\tFound " << opt.train_cell_to_synergy_file.size() << " training synergy files to read" << endl;
		out << "\tFound " << opt.test_cell_to_synergy_file.size() << " testing synergy files to read" << endl;
		
		if( !opt.test_cell_to_synergy_file.empty() ){
			out << "\tSetting number of cv folds to 1 to enable use of the explicit test set" << endl;
		}
		
		if( !opt.file_name_blind_drug.empty() && !opt.file_name_blind_cell.empty() ){

			out << "\tBlind drug features will be from " << opt.file_name_blind_drug << endl;
			out << "\tBlind cell features will be from " << opt.file_name_blind_cell << endl;
		}
		
		if(opt.use_pair_prediction){
			out << "\tPerforming both drug pair and single drug-based synergy prediction" << endl;
		}
		else{
			out << "\tPerforming single drug-based synergy prediction" << endl;
		}
		
		if(opt.randomize_test_synergy){
			out << "\t** Randomizing testing synergy variables **" << endl;
		}
		
		if(opt.randomize_train_synergy){
			out << "\t** Randomizing training synergy variables **" << endl;
		}
		
		if(opt.NSC_to_CID_y){
			out << "\t** Converting synergy variable drug IDs from NSC to CID **" << endl;
		}
		
		if(opt.one_hot_cell){
			out << "\tUsing one-hot-encoded cell lines" << endl;
		}
		
		if(opt.one_hot_tissue){
			out << "\tUsing one-hot-encoded cell line tissue types" << endl;
		}
		
		if(opt.one_hot_drug){
			out << "\tUsing one-hot-encoded drugs" << endl;
		}
		
		out << "\tRandom Forest Parameters:" << endl;
		out << "\t\tnum_tree = " << opt.num_tree << endl;
		out << "\t\tleaf_size = " << opt.leaf_size << endl;
		out << "\t\tbag_fraction = " << opt.bag_fraction << endl;
		
		out << "\tRandom number seed = " << opt.seed << endl;
		out << "\tNumber of cross validation folds = " << opt.num_fold << endl;
		
		switch(opt.score_func){
			case AVERAGE_SCORE:

				out << "\tPredicting average synergy scores" << endl;
				break;
			case MEDIAN_SCORE:

				out << "\tPredicting median synergy scores" << endl;
				break;
			case FILTERED_MEAN_SCORE:

				out << "\tPredicting filtered mean synergy scores" << endl;
				break;
			case COUNTING_SCORE:

				out << "\tPredicting probability of being synergistic (w/ synergy threshold = " 
					<< opt.synergy_threshold << ")" << endl;
				break;
			default:
				throw __FILE__ ":synergy_score: Unknown synergy scoring function";
		};
		
		// The name of each drug feature (for summarizing informative features at the end)
		vector<string> cell_feature_id;
		MAP<string, vector<float> > cell_features;
		MAP<string, vector<float> > blind_test_cell_features;
		MULTIMAP<string, string> blind_test_cell_prefix_map; // For grouping NCIPDM models by prefix
		vector<string> blind_test_prefix;
		
		// Real valued descriptor file and make sure that we have cell line features for 
		// all of the synergy matricies
		bool valid_cell_features = true;
		
		if(mpi_rank == 0){
			
			// Cell line features are optional, but only attempt to parse them if the user
			// has specified an input file
			if( !opt.file_name_cell.empty() ){
			
				cell_features = parse_descriptor(opt.file_name_cell, cell_feature_id, "Sample", '\t');
				
				if(opt.normalize_cell){
					normalize(cell_features);
				}
				
				for(MAP<string, string>::const_iterator i = opt.train_cell_to_synergy_file.begin();
					i != opt.train_cell_to_synergy_file.end();++i){

					if( cell_features.find(i->first) == cell_features.end() ){

						cerr << "Unable to find cell features for training cell line: " << i->first << endl;
						valid_cell_features = false;
					}
				}

				for(MAP<string, string>::const_iterator i = opt.test_cell_to_synergy_file.begin();
					i != opt.test_cell_to_synergy_file.end();++i){

					if( cell_features.find(i->first) == cell_features.end() ){

						cerr << "Unable to find cell features for testing cell line: " << i->first << endl;
						valid_cell_features = false;
					}
				}
				
				out << "Found features for " << cell_features.size() << " cell lines" << endl;
			}
			
			if(opt.one_hot_cell){
				
				one_hot_encode( cell_features, keys(opt.train_cell_to_synergy_file) );
				
				const size_t one_hot_dim = 
					one_hot_encode( cell_features, keys(opt.test_cell_to_synergy_file) );
				
				out << "One hot encoded cell line features in " << one_hot_dim 
					<< " dimensions" << endl;
			}
			
			if(opt.one_hot_tissue){
				
				one_hot_encode_by_tissue( cell_features, keys(opt.train_cell_to_synergy_file) );
				
				const vector<string> one_hot_tissues = 
					one_hot_encode_by_tissue( cell_features, keys(opt.test_cell_to_synergy_file) );				
					
				out << "One hot encoded cell line features by tissue in " << one_hot_tissues.size() 
					<< " dimensions:" << endl;
				
				// List the tissues
				for(vector<string>::const_iterator i = one_hot_tissues.begin();
					i != one_hot_tissues.end();++i){
					
					out << '\t' << *i << endl;
				}
			}
		}
		
		broadcast(valid_cell_features, mpi_rank, 0);
		
		if(valid_cell_features == false){
			
			MPI_Finalize();
			return EXIT_FAILURE;
		}
		
		broadcast(cell_features, mpi_rank, 0);
						
		if( (mpi_rank == 0) && !opt.file_name_blind_cell.empty() ){
			
			vector<string> dummy_id;
			
			blind_test_cell_features = parse_descriptor(opt.file_name_blind_cell, dummy_id, "Sample", '\t');
			
			if(opt.normalize_cell){
				normalize(blind_test_cell_features);
			}

			for(MAP<string, vector<float> >::const_iterator i = blind_test_cell_features.begin();
				i != blind_test_cell_features.end();++i){
				
				// Extract the prefix (the PDM models have a name the is delimeted by '~')
				const string::size_type loc = i->first.find('~');
				
				if(loc == string::npos){
					
					cerr << "Warning: Did not find a blind test cell name delimeter" << endl;
					
					blind_test_cell_prefix_map.insert( make_pair(i->first, i->first) );
				}
				else{
					const string prefix = i->first.substr(0, loc);
					
					blind_test_cell_prefix_map.insert( make_pair(prefix, i->first) );
				}
			}
			
			blind_test_prefix = keys(blind_test_cell_prefix_map);
			
			out << "Found " << blind_test_cell_features.size() << " blind test cells" << endl;
			out << "Found " << blind_test_prefix.size() << " blind test cell prefixes" << endl;
		}
		
		broadcast(blind_test_cell_features, mpi_rank, 0);
		broadcast(blind_test_prefix, mpi_rank, 0);
		
		MAP< string /*cell*/, 
			MAP< string /*drug*/, 
				MAP<string /*drug*/, float> > > train_synergy_matrix;
		
		MAP< string /*cell*/, 
			MAP< string /*drug*/, 
				MAP<string /*drug*/, float> > > test_synergy_matrix;
		
		if(mpi_rank == 0){

			// Load all of the training synergy matricies
			for(MAP<string, string>::const_iterator i = opt.train_cell_to_synergy_file.begin();
				i != opt.train_cell_to_synergy_file.end();++i){

				train_synergy_matrix[i->first] = parse_matrix(i->second);
			}

			// Load all of the testing synergy matricies
			for(MAP<string, string>::const_iterator i = opt.test_cell_to_synergy_file.begin();
				i != opt.test_cell_to_synergy_file.end();++i){

				test_synergy_matrix[i->first] = parse_matrix(i->second);
			}
			
			/////////////////////////////////////////////////////////////////////////////////
			// Convert the drug names from NSC to CID if requested
			/////////////////////////////////////////////////////////////////////////////////
			if(opt.NSC_to_CID_y){

				out << "Converting training synergy matrix drug names from NSC to CID" << endl;

				for(MAP<string, MAP<string, MAP<string, float> > >::iterator i = train_synergy_matrix.begin();
					i != train_synergy_matrix.end();++i){

					i->second = NSC_to_CID(i->second);
				}

				if(!test_synergy_matrix.empty() ){

					out << "Converting testing synergy matrix drug names from NSC to CID" << endl;

					for(MAP<string, MAP<string, MAP<string, float> > >::iterator i = test_synergy_matrix.begin();
						i != test_synergy_matrix.end();++i){

						i->second = NSC_to_CID(i->second);
					}
				}
			}
		}
		
		broadcast(train_synergy_matrix, mpi_rank, 0);
		broadcast(test_synergy_matrix, mpi_rank, 0);
		
		out << 100.0*fraction_synergy(train_synergy_matrix, opt.synergy_threshold)
			<< "% of training measurement are synergistic" << endl;
			
		if( !test_synergy_matrix.empty() ){
		
			out << 100.0*fraction_synergy(test_synergy_matrix, opt.synergy_threshold)
				<< "% of testing measurement are synergistic" << endl;
		}
		
		vector<string> training_cells = keys(train_synergy_matrix);
		vector<string> testing_cells = keys(test_synergy_matrix);
		
		out << "Found " << training_cells.size() << " training cell lines" << endl;
		out << "Found " << testing_cells.size() << " testing cell lines" << endl;
		
		// The dependent variable drugs from the *training* synergy matricies
		vector<string> training_drugs;
		
		for(MAP<string, MAP<string, MAP<string, float> > >::iterator i = train_synergy_matrix.begin();
			i != train_synergy_matrix.end();++i){
			
			if( training_drugs.empty() ){
				
				// Extract the drugs from the first synergy matrix
				training_drugs = keys(i->second);
			}
			else{
				
				// Make sure that all synergy matricies have the same set of drugs
				if( training_drugs != keys(i->second) ){
					throw __FILE__ ":main: Inconsistent drugs in different training synergy matricies";
				}
			}
		}
		
		//#define EXTRACT_CELL_LINE_SYNERGY_DISTRIBUTION
		#ifdef EXTRACT_CELL_LINE_SYNERGY_DISTRIBUTION
		// Extract the distribution of cell lines that show synergy for a given drug pair
		
		vector<size_t> histogram(train_synergy_matrix.size() + 1);
		
		for(vector<string>::const_iterator i = training_drugs.begin();i != training_drugs.end();++i){
			for(vector<string>::const_iterator j = i + 1;j != training_drugs.end();++j){
				
				size_t cell_count = 0; // The number of valid synergy measurements
				size_t synergy_count = 0; // The number of cell lines that display synergy
				                          // for this drug pair
				
				for(MAP<string, MAP<string, MAP<string, float> > >::iterator k = 
					train_synergy_matrix.begin();
					k != train_synergy_matrix.end();++k){
					
					MAP<string, MAP<string, float> >::const_iterator row_iter = 
						k->second.find(*i);
					
					if( row_iter == k->second.end() ){
						continue;
					}
					
					MAP<string, float>::const_iterator col_iter = row_iter->second.find(*j);
					
					if( col_iter == row_iter->second.end() ){
						continue;
					}
					
					if(col_iter->second == MISSING_DATA){
						continue;
					}
					
					++cell_count;
					
					synergy_count += (col_iter->second <= opt.synergy_threshold);
				}
				
				//if(mpi_rank == 0){
					//cout << synergy_count << '\t' << cell_count << endl;
				//}
					
				if(synergy_count >= 20){

					// Print the identities of "pan-cell line" synergistic drugs
					//if(mpi_rank == 0){
					//	cout << synergy_count << '\t' << *i << '\t' << *j << endl;
					//}
					// **Remove** these "pan-cell line" synergistic drug pairs!
					for(MAP<string, MAP<string, MAP<string, float> > >::iterator k = 
						train_synergy_matrix.begin();
						k != train_synergy_matrix.end();++k){

						k->second[*i][*j] = k->second[*j][*i] = 1.0;
					}
				}
				
				++histogram[synergy_count];
			}
		}
		
		//if(mpi_rank == 0){
		//	for(size_t i = 0;i < histogram.size();++i){
		//		cout << i << '\t' << histogram[i] << endl;
		//	}
		//}
		
		#endif // EXTRACT_CELL_LINE_SYNERGY_DISTRIBUTION
		
		// Uncomment to produce a heat map (i.e. matrix) that plots the singe drug synergy
		// probability for each drug/cell line pair
		//out << probability_heat_map(train_synergy_matrix, opt.synergy_threshold) << endl;
		
		// The dependent variable drugs from the *testing* synergy matricies
		vector<string> testing_drugs;
		
		for(MAP<string, MAP<string, MAP<string, float> > >::iterator i = test_synergy_matrix.begin();
			i != test_synergy_matrix.end();++i){
			
			if( testing_drugs.empty() ){
				
				// Extract the drugs from the first test synergy matrix
				testing_drugs = keys(i->second);
			}
			else{
				
				// Make sure that all testing synergy matricies have the same set of drugs
				if( testing_drugs != keys(i->second) ){
					throw __FILE__ ":main: Inconsistent drugs in different testing synergy matricies";
				}
			}
		}
		
		// The name of each drug feature (for summarizing informative features at the end)
		vector<string> drug_feature_id;
		MAP<string, vector<float> > drug_features;
		MAP<string, vector<float> > blind_test_drug_features;
		
		if(mpi_rank == 0){
		
			// Drug features are optional and only attempt to parse them if the user
			// has specified an input file
			if( !opt.file_name_drug.empty() ){
			
				if( is_fingerprint(opt.file_name_drug) ){

					out << "Reading drug features as fingerprints" << endl;

					// Binary fingerprint file
					drug_features = parse_fingerprint(opt.file_name_drug, drug_feature_id,
						"DRUG", "FINGERPRINT");

					if( !opt.file_name_blind_drug.empty() ){

						if(is_fingerprint(opt.file_name_blind_drug) == false){
							throw __FILE__ ":main: Detected a fingerprint-based training set and a descriptor-based blind test set!";
						}

						vector<string> dummy_id;

						blind_test_drug_features = parse_fingerprint(opt.file_name_blind_drug, dummy_id,
							"DRUG", "FINGERPRINT");
					}
				}
				else{

					out << "Reading drug features as descriptors" << endl;

					// Real valued descriptor file
					drug_features = parse_descriptor(opt.file_name_drug, drug_feature_id, "DRUG");

					if( !opt.file_name_blind_drug.empty() ){

						if(is_fingerprint(opt.file_name_blind_drug) == true){
							throw __FILE__ ":main: Detected a descriptor-based training set and a fingerprint-based blind test set!";
						}

						vector<string> dummy_id;

						blind_test_drug_features = parse_descriptor(opt.file_name_blind_drug, dummy_id,
							"DRUG");
					}
				}
				
				#ifdef REMOVE_INVARIANT_DRUG_FEATURES
				// For many types of molelcular fingerprints and descriptors, there are elements of the feature vector
				// that are invariant across all input drugs. Remove these invariant drug features.
				const vector<bool> feature_mask = 
					find_invariant_features(drug_features, blind_test_drug_features);
				
				remove_invariant_features(drug_feature_id, feature_mask);
				
				for(MAP<string, vector<float> >::iterator i = drug_features.begin();
					i != drug_features.end();++i){
					
					remove_invariant_features(i->second, feature_mask);
				}
				
				for(MAP<string, vector<float> >::iterator i = blind_test_drug_features.begin();
					i != blind_test_drug_features.end();++i){
					
					remove_invariant_features(i->second, feature_mask);
				}
				
				out << "Found " << feature_mask.size() << " initial drug features" << endl;
				out << '\t' << count(feature_mask) << " of these features are invariant" << endl;
				out << '\t' << feature_mask.size() - count(feature_mask) << " variable features remain" << endl;
				#endif // REMOVE_INVARIANT_DRUG_FEATURES
			}
			
			if(opt.one_hot_drug){
				
				one_hot_encode(drug_features, training_drugs);
				const size_t one_hot_dim = one_hot_encode(drug_features, testing_drugs);
				
				out << "One hot encoded drug features in " << one_hot_dim 
					<< " dimensions" << endl;
					
				drug_feature_id = keys(drug_features);
			}

			// Has the user specified any drug features? If not, create a set of dummy drug features
			if( drug_features.empty() ){
			
				for(vector<string>::const_iterator i = training_drugs.begin();i != training_drugs.end();++i){
					drug_features[*i] = vector<float>(1, 1.0);
				}

				for(vector<string>::const_iterator i = testing_drugs.begin();i != testing_drugs.end();++i){   
                                        drug_features[*i] = vector<float>(1, 1.0);
                                }

				drug_feature_id = vector<string>(1, "dummy");
			}
		}
		
		broadcast(drug_feature_id, mpi_rank, 0);
		broadcast(drug_features, mpi_rank, 0);
		
		// Parse the drug target information (if available). This target information is only used for
		// the blind test drugs
		MAP<string, string> drug_target;
		
		if(mpi_rank == 0){
		
			parse_drug_target(opt.file_name_drug_target);
		
			if( !blind_test_drug_features.empty() ){

				MAP<string, vector<float> > non_redundant_blind_test_drug_features;

				out << "Found " << blind_test_drug_features.size() << " blind test drugs" << endl;

				// Merge drugs with identical feature vectors
				while(true){

					bool found_match = false;

					for(MAP<string, vector<float> >::const_iterator i = blind_test_drug_features.begin();
						i != blind_test_drug_features.end();++i){

						deque<string> matches;

						for(MAP<string, vector<float> >::const_iterator j = blind_test_drug_features.begin();
							j != blind_test_drug_features.end();++j){

							if(i == j){
								continue;
							}

							if(i->second == j->second){
								matches.push_back(j->first);
							}
						}

						if( matches.empty() ){

							non_redundant_blind_test_drug_features[i->first] = i->second;
							continue;
						}

						found_match = true;

						string composite_name = i->first;

						for(deque<string>::const_iterator j = matches.begin();j != matches.end();++j){
							composite_name +=  "|" + *j;
						}

						non_redundant_blind_test_drug_features[composite_name] = i->second;

						// Create composite records in the drug_target map
						MAP<string, string>::const_iterator target_iter = drug_target.find(i->first);
						string consensus_target;

						if( target_iter != drug_target.end() ){

							consensus_target = target_iter->second;
							drug_target.erase(target_iter);
						}

						for(deque<string>::const_iterator j = matches.begin();j != matches.end();++j){

							target_iter = drug_target.find(*j);

							if( target_iter != drug_target.end() ){

								// We if don't have a consensus target, set it now 
								if( consensus_target.empty() ){
									consensus_target = target_iter->second;
								}

								if(consensus_target != target_iter->second){

									out << "**Warning** Inconsistent blind drug targets: "
										<< consensus_target << " != " << target_iter->second << endl;
								}

								drug_target.erase(target_iter);
							}
						}

						if( !consensus_target.empty() ){
							drug_target.insert( make_pair(composite_name, consensus_target) );
						}

						// Remove the cluster of identical drugs from the blind target list
						blind_test_drug_features.erase(i->first);

						for(deque<string>::const_iterator j = matches.begin();j != matches.end();++j){
							blind_test_drug_features.erase(*j);
						}

						break;
					}

					if(found_match == false){
						break;
					}
				}

				blind_test_drug_features = non_redundant_blind_test_drug_features;

				out << "Found " << blind_test_drug_features.size() << " *non-redundant* blind test drugs" << endl;
			}		
		}
		
		broadcast(drug_target, mpi_rank, 0);
		broadcast(blind_test_drug_features, mpi_rank, 0);
		
		// The independent variable drugs from the fingerprint/descriptor data
		vector<string> independent_drugs = keys(drug_features);
		
		out << "Found training synergy data for " << training_drugs.size() << " drugs" << endl;
		
		if( !testing_drugs.empty() ){
			out << "Found testing synergy data for " << testing_drugs.size() << " drugs" << endl;
		}
		
		out << "Found fingerprint/descriptor data for " << independent_drugs.size() << " drugs" << endl;
		
		// Since we need both dependent and independent data, we will restrict the analysis 
		// to the *intersection* of these two drug lists
		
		vector<string> drugs = intersection(training_drugs, independent_drugs);
		 
		out << "Found " << drugs.size() << " drugs in the intersection of training drugs and drugs with descriptors/fingerprints" << endl;
		 
		// Remove drugs from the dependent *training* variables (i.e. training synergy matrix) that do *not* have fingerprint data
		for(vector<string>::const_iterator i = training_drugs.begin();i != training_drugs.end();++i){
			
			if( find(drugs.begin(), drugs.end(), *i) == drugs.end() ){
				
				// This drug is *NOT* in the set of final drugs
				
				// For each synergy matrix ...
				for(MAP<string, MAP< string, MAP<string, float> > >::iterator j = train_synergy_matrix.begin();
					j != train_synergy_matrix.end();++j){
					
					// Remove the row corresponding to this drug
					j->second.erase(*i);
					
					// Remove the columns corresponding to this drug
					for(MAP< string, MAP<string, float> >::iterator k = j->second.begin();k != j->second.end();++k){
						k->second.erase(*i);
					}
				}
			}
			else{
				// This drug *is* present in the training synergy matrix. Make sure we have enough valid measurements
				
				// For each synergy matrix ...
				for(MAP<string, MAP< string, MAP<string, float> > >::iterator j = train_synergy_matrix.begin();
					j != train_synergy_matrix.end();++j){
					
					MAP< string, MAP<string, float> >::iterator row_iter = j->second.find(*i);

					if( row_iter == j->second.end() ){
						throw __FILE__ ":main: Unable to look up training synergy row by drug name";
					}

					size_t num_valid = 0;
					size_t num_invalid = 0;

					for(MAP<string, float>::iterator col_iter = row_iter->second.begin();
						col_iter != row_iter->second.end();++col_iter){

						if(col_iter->second == MISSING_DATA){
							++num_invalid;
						}
						else{
							++num_valid;
						}
					}

					// For now, only remove rows/column that are *entirely* made up of missing data!
					if(num_valid == 0){
					
						// Remove the row corresponding to this drug					
						j->second.erase(*i);

						// Remove the columns corresponding to this drug
						for(MAP< string, MAP<string, float> >::iterator k = j->second.begin();k != j->second.end();++k){
                                                	k->second.erase(*i);
                                        	}	
					}
				}				
			}
		}
		
		// Recompute the drugs for which we have both dependent training (i.e. training synergy) 
		// and independent (fingerprint) data		
		unordered_set<string> drug_union;
		
		for(MAP<string, MAP<string, MAP<string, float> > >::iterator i = train_synergy_matrix.begin();
			i != train_synergy_matrix.end();++i){
			
			const vector<string> local = keys(i->second);
			
			// The union of all drugs from all synergy matricies
			for(vector<string>::const_iterator j = local.begin();j != local.end();++j){
				drug_union.insert(*j);
			}
		}
		
		// The final set of dependent drugs is the *union* of the valid drugs from
		// each training synergy matrix
		training_drugs.assign( drug_union.begin(), drug_union.end() );
		
		drug_union.clear();
		
		out << "Found " << training_drugs.size() << " drugs with training synergy AND fingerprint data" << endl;
		
		if( !testing_drugs.empty() ){
		
			drugs = intersection(testing_drugs, independent_drugs);
		 
			out << "Found " << drugs.size() << " drugs in the intersection of testing drugs and drugs with descriptors/fingerprints" << endl;

			// Remove drugs from the dependent *testing* variables (i.e. testing synergy matrix) that do *not* have fingerprint data
			for(vector<string>::const_iterator i = testing_drugs.begin();i != testing_drugs.end();++i){

				if( find(drugs.begin(), drugs.end(), *i) == drugs.end() ){

					// This drug is *NOT* in the set of final drugs

					// For each synergy matrix ...
					for(MAP<string, MAP< string, MAP<string, float> > >::iterator j = test_synergy_matrix.begin();
						j != test_synergy_matrix.end();++j){

						// Remove the row corresponding to this drug
						j->second.erase(*i);

						// Remove the columns corresponding to this drug
						for(MAP< string, MAP<string, float> >::iterator k = j->second.begin();k != j->second.end();++k){
							k->second.erase(*i);
						}
					}
				}
				else{
					// This drug *is* present in the testing synergy matrix. Make sure we have enough valid measurements

					// For each synergy matrix ...
					for(MAP<string, MAP< string, MAP<string, float> > >::iterator j = test_synergy_matrix.begin();
						j != test_synergy_matrix.end();++j){

						MAP< string, MAP<string, float> >::iterator row_iter = j->second.find(*i);

						if( row_iter == j->second.end() ){
							throw __FILE__ ":main: Unable to look up testing synergy row by drug name";
						}

						size_t num_valid = 0;
						size_t num_invalid = 0;

						for(MAP<string, float>::iterator col_iter = row_iter->second.begin();
							col_iter != row_iter->second.end();++col_iter){

							if(col_iter->second == MISSING_DATA){
								++num_invalid;
							}
							else{
								++num_valid;
							}
						}

						// For now, only remove rows/column that are *entirely* made up of missing data!
						if(num_valid == 0){
							
							// Remove the row corresponding to this drug
							j->second.erase(*i);

							// Remove the columns corresponding to this drug	
							for(MAP< string, MAP<string, float> >::iterator k = j->second.begin();k != j->second.end();++k){
                                                        	k->second.erase(*i);
                                                	}	
						}
					}
				}
			}

			// Recompute the drugs for which we have both dependent testing (i.e. testing synergy) 
			// and independent (fingerprint) data		
			for(MAP<string, MAP<string, MAP<string, float> > >::iterator i = test_synergy_matrix.begin();
				i != test_synergy_matrix.end();++i){

				const vector<string> local = keys(i->second);

				// The union of all drugs from all synergy matricies
				for(vector<string>::const_iterator j = local.begin();j != local.end();++j){
					drug_union.insert(*j);
				}
			}

			// The final set of dependent drugs is the *union* of the valid drugs from
			// each testing synergy matrix
			testing_drugs.assign( drug_union.begin(), drug_union.end() );

			drug_union.clear();
		
			out << "Found " << testing_drugs.size() << " drugs with testing synergy AND fingerprint data" << endl;
		}
		
		size_t num_blind_with_target = 0;
		size_t num_blind_no_target = 0;
		
		for(MAP<string, vector<float> >::const_iterator i = blind_test_drug_features.begin();
			i != blind_test_drug_features.end();++i){
			
			if( drug_target.find(i->first) == drug_target.end() ){
				++num_blind_no_target;
			}
			else{
				++num_blind_with_target;
			}
		}
		
		if( !blind_test_drug_features.empty() ){
		
			out << num_blind_with_target << " blind drugs have target information" << endl;
			out << num_blind_no_target << " blind drugs do not have target information" << endl;
		}
		
		if( opt.overlap_to_train && !testing_cells.empty() ){
			
			testing_cells = set_difference(testing_cells, training_cells);
			testing_drugs = set_difference(testing_drugs, training_drugs);
			
			out << "Assigned overlapping cell lines and drugs to the *training* set" << endl;
			out << "\t|test cell lines| = " << testing_cells.size() << endl;
			out << "\t|test drugs| = " << testing_drugs.size() << endl;
			
			// Update the synergy matricies to ensure that we don't reintroduce the overlap
			// during synergy matrix permutation
			for(vector<string>::const_iterator i = training_cells.begin();i != training_cells.end();++i){
				
				// Remove the training cell lines from the test synergy values
				if( test_synergy_matrix.find(*i) != test_synergy_matrix.end() ){
				
					test_synergy_matrix.erase(*i);
					
					out << "\t\tRemoved " << *i << " from the test set" << endl;
				}
			}
			
			for(MAP<string, MAP< string, MAP<string, float> > >::iterator i = test_synergy_matrix.begin();
				i != test_synergy_matrix.end();++i){
				
				// Remove the training drugs from each of the individual test synergy matricies
				for(vector<string>::const_iterator j = training_drugs.begin();j != training_drugs.end();++j){
					remove_drug(i->second, *j);
				}
			}
		}
		
		if( opt.overlap_to_test && !testing_cells.empty() ){
			
			training_cells = set_difference(training_cells, testing_cells);
			training_drugs = set_difference(training_drugs, testing_drugs);
			
			out << "Assigned overlapping cell lines and drugs to the *test* set" << endl;
			out << "\t|training drugs| = " << training_drugs.size() << endl;
			out << "\t|training cell lines| = " << training_cells.size() << endl;
			
			// Update the synergy matricies to ensure that we don't reintroduce the overlap
			// during synergy matrix permutation
			for(vector<string>::const_iterator i = testing_cells.begin();i != testing_cells.end();++i){
				
				// Remove the test cell lines from the training synergy values
				if( train_synergy_matrix.find(*i) != train_synergy_matrix.end() ){
				
					train_synergy_matrix.erase(*i);
					
					out << "\t\tRemoved " << *i << " from the training set" << endl;
				}
			}
			
			for(MAP<string, MAP< string, MAP<string, float> > >::iterator i = train_synergy_matrix.begin();
				i != train_synergy_matrix.end();++i){
				
				// Remove the testing drugs from each of the individual training synergy matricies
				for(vector<string>::const_iterator j = testing_drugs.begin();j != testing_drugs.end();++j){
					remove_drug(i->second, *j);
				}
			}
		}
		
		// Compute the range of values in each *training* synergy matrix
		for(MAP< string, MAP< string, MAP<string, float> > >::const_iterator i = train_synergy_matrix.begin();i != train_synergy_matrix.end();++i){
		
			float min_synergy = 0.0;
			float max_synergy = 0.0;
			float ave_synergy = 0.0;
			size_t synergy_norm = 0;
			size_t num_synergistic = 0;
			size_t num_non_synergistic = 0;

			for(MAP< string, MAP<string, float> >::const_iterator j = i->second.begin();j != i->second.end();++j){

				for(MAP<string, float>::const_iterator k = j->second.begin();k != j->second.end();++k){

					if(k->second == MISSING_DATA){
						continue;
					}

					if(synergy_norm == 0){
						min_synergy = max_synergy = k->second;
					}
					else{
						min_synergy = min(min_synergy, k->second);
						max_synergy = max(max_synergy, k->second);
					}

					if(k->second <= opt.synergy_threshold){
						++num_synergistic;
					}
					else{
						++num_non_synergistic;
					}

					ave_synergy += k->second;
					++synergy_norm;
				}
			}

			if(synergy_norm > 0){
				ave_synergy /= synergy_norm;
			}

			out << "For training cell: " << i->first << endl;
			out << '\t' << synergy_norm << " synergy values range from " << min_synergy << " to " << max_synergy 
				<< " with an average of " << ave_synergy << endl;

			if(opt.score_func == COUNTING_SCORE){

				out << '\t' << (100.0*num_synergistic)/(num_synergistic + num_non_synergistic) 
					<< "% of drug pairs are defined as synergistic" << endl;
				out << '\t' << (100.0*num_non_synergistic)/(num_synergistic + num_non_synergistic) 
					<< "% of drug pairs are defined as non-synergistic" << endl;

				if(opt.synergy_threshold <= min_synergy){
					out << '\t' << "!!Warning!! All synergy values are greater than the specificy synergy threshold (all interactions are non-synergistic)" << endl;
				}

				if(opt.synergy_threshold >= max_synergy){
					out << '\t' << "!!Warning!! All synergy values are less than the specificy synergy threshold (all interactions are synergistic)" << endl;
				}
			}
		}
		
		// Compute the range of values in each *testing* synergy matrix
		for(MAP< string, MAP< string, MAP<string, float> > >::const_iterator i = test_synergy_matrix.begin();i != test_synergy_matrix.end();++i){
		
			float min_synergy = 0.0;
			float max_synergy = 0.0;
			float ave_synergy = 0.0;
			size_t synergy_norm = 0;
			size_t num_synergistic = 0;
			size_t num_non_synergistic = 0;

			for(MAP< string, MAP<string, float> >::const_iterator j = i->second.begin();j != i->second.end();++j){

				for(MAP<string, float>::const_iterator k = j->second.begin();k != j->second.end();++k){

					if(k->second == MISSING_DATA){
						continue;
					}

					if(synergy_norm == 0){
						min_synergy = max_synergy = k->second;
					}
					else{
						min_synergy = min(min_synergy, k->second);
						max_synergy = max(max_synergy, k->second);
					}

					if(k->second <= opt.synergy_threshold){
						++num_synergistic;
					}
					else{
						++num_non_synergistic;
					}

					ave_synergy += k->second;
					++synergy_norm;
				}
			}

			if(synergy_norm > 0){
				ave_synergy /= synergy_norm;
			}

			out << "For testing cell: " << i->first << endl;
			out << '\t' << synergy_norm << " synergy values range from " << min_synergy << " to " << max_synergy 
				<< " with an average of " << ave_synergy << endl;

			if(opt.score_func == COUNTING_SCORE){

				out << '\t' << (100.0*num_synergistic)/(num_synergistic + num_non_synergistic) 
					<< "% of drug pairs are defined as synergistic" << endl;
				out << '\t' << (100.0*num_non_synergistic)/(num_synergistic + num_non_synergistic) 
					<< "% of drug pairs are defined as non-synergistic" << endl;

				if(opt.synergy_threshold <= min_synergy){
					out << '\t' << "!!Warning!! All synergy values are less than the specificy synergy threshold (all interactions are synergistic)" << endl;
				}

				if(opt.synergy_threshold >= max_synergy){
					out << '\t' << "!!Warning!! All synergy values are greater than the specificy synergy threshold (all interactions are non-synergistic)" << endl;
				}
			}
		}
		
		// Randomize the independent variables 
		if(opt.randomize_drug){
			
			out << "** Randomizing independent (feature) variables **" << endl;
			
			randomize_independent(drug_features, &opt.seed);
			
			// rank 0 defines the ordering
			broadcast(drug_features, mpi_rank, 0);
		}
		
		if(opt.randomize_train_synergy){
			
			out << "** Randomizing dependent (training synergy matrix) variables **" << endl;
			out << "** This will also randomize the synergy matricies used for testing in cross validation **" << endl;
			
			randomize_dependent(train_synergy_matrix, &opt.seed);
			
			// rank 0 defines the ordering
			broadcast(train_synergy_matrix, mpi_rank, 0);
		}
		
		// Randomly shuffle the drugs prior to cross validation
		randomize(training_drugs.begin(), training_drugs.end(), &opt.seed);
		
		// rank 0 defines the ordering of training drugs
		broadcast(training_drugs, mpi_rank, 0);
				
		// Randomly shuffle the training cells prior to cross validation
		randomize(training_cells.begin(), training_cells.end(), &opt.seed);
		
		// rank 0 defines the ordering of training cells
		broadcast(training_cells, mpi_rank, 0);
		
		// The mean squared error for single drug synergy prediction
		double mse = 0.0;
		double average_mse = 0.0;
		double mse_norm = 0.0;
		
		// The mean squared error for single drug synergy prediction
		double pair_mse = 0.0;
		double pair_average_mse = 0.0;
		double pair_mse_norm = 0.0;
		
		// For per-fold classification/counting-based assessments of enrichment
		deque<float> fractional_enrichment;
		deque<float> single_auroc;
		deque<float> single_gini;
		
		deque<float> fractional_pair_enrichment;
		deque<float> pair_auroc;
		deque<float> pair_gini;
		
		deque< pair<float, float> > actual_vs_predicted;
		
		// Compute the correlation between predicted and observed drug *pair* 
		// synergy values for the test set
		deque< pair<float, float> > pair_synergy_comparison;
			
		time_t profile = time(NULL);
		
		const bool explicit_test_provided = 
			(opt.num_fold == 1) && 
			!testing_drugs.empty() && 
			!testing_cells.empty();
		
		for(size_t fold = 0;fold < opt.num_fold;++fold){
			
			out << "Fold " << fold + 1 << endl;
			
			// Split the drugs in to a training set and testing set
			vector<string> cv_training_drugs;
			vector<string> cv_testing_drugs;
			
			vector<string> cv_training_cells;
			vector<string> cv_testing_cells;
			
			if(explicit_test_provided){
			
				// The user has explicitly provided a testing set
				cv_testing_drugs = testing_drugs;
				cv_training_drugs = training_drugs;
				
				cv_testing_cells = testing_cells;
				cv_training_cells = training_cells;
			}
			else{
				// Standard cross validation that splits the training set
				// into test and training folds
				const size_t num_train_drugs = training_drugs.size();
				
				for(size_t i = 0;i < num_train_drugs;++i){

					if(i%opt.num_fold == fold){

						cv_testing_drugs.push_back(training_drugs[i]);
					}
					else{
						cv_training_drugs.push_back(training_drugs[i]);
					}
				}

				const size_t num_training_cells = training_cells.size();
				
				for(size_t i = 0;i < num_training_cells;++i){

					if(i%opt.num_fold == fold){

						cv_testing_cells.push_back(training_cells[i]);
					}
					else{
						cv_training_cells.push_back(training_cells[i]);
					}
				}
				
				// Fewer cells than cross validation folds are provided, use *all*
				// cells for both training and testing
				if(num_training_cells < opt.num_fold){
					
					out << "\t*Using all cells for test and training*" << endl;
					
					cv_testing_cells = training_cells;
					cv_training_cells = training_cells;
				}
				
			}

			out << "\t|training drugs| = " << cv_training_drugs.size() << endl;
			out << "\t|training cells| = " << cv_training_cells.size() << endl;
			
			out << "\t|testing drugs| = " << cv_testing_drugs.size() << endl;
			out << "\t|testing cells| = " << cv_testing_cells.size() << endl;
			
			// Pack the training data for regression. 
			vector<LabeledData> train;
			vector<LabeledData> train_pair;
			size_t num_features = 0;
			
			train.reserve( cv_training_drugs.size() * cv_training_cells.size() );
			
			if(opt.use_pair_prediction){
				
				// NumDrugs*(NumDrugs - 1)*NumCells/2
				train_pair.reserve( cv_training_drugs.size() * ( cv_training_drugs.size() - 1) * 
					cv_training_cells.size()/2 );
			}
			
			////////////////////
			// Training
			////////////////////
			for(vector<string>::const_iterator c = cv_training_cells.begin();c != cv_training_cells.end();++c){

				MAP<string, MAP< string, MAP<string, float> > >::const_iterator synergy_iter = train_synergy_matrix.find(*c);
					
				if( synergy_iter == train_synergy_matrix.end() ){
					throw __FILE__ ":main: Unable to look up training synergy matrix";
				}
					
				for(vector<string>::const_iterator d_1 = cv_training_drugs.begin();d_1 != cv_training_drugs.end();++d_1){

					MAP<string, vector<float> >::const_iterator drug_features_iter = drug_features.find(*d_1);
				
					if( drug_features_iter == drug_features.end() ){
						throw __FILE__ ":main: Unable to look up fingerprint data for selected training drug (1)";
					}
					
					MAP< string, MAP<string, float> >::const_iterator row_iter = synergy_iter->second.find(*d_1);

					if( row_iter == synergy_iter->second.end() ){
						
						// Synergy matricies may be missing drugs
						continue;
					}

					if( d_1 == cv_training_drugs.begin() ){
						num_features = drug_features_iter->second.size();
					}

					if( num_features != drug_features_iter->second.size() ){
						throw __FILE__ ":main: Variable number of features detected";
					}

					// Only use the TRAINING drugs to compute the synergy score.
					// The synergy score can be:
					// (a) fraction of  synergistic interactions
					// (b) median synergy value (over all drug pairs that contain this drug)
					const pair<float, float> score = synergy_score(
						cv_training_drugs, 
						row_iter->second, 
						opt.score_func, 
						opt.synergy_threshold);

					if(score.second <= 0.0){
					
						// A score weight <= 0.0 indicates too much missing
						// data to compute the score
						continue;
					}

					train.push_back( LabeledData() );

					LabeledData &ref = train.back();

					ref.response = score.first;
					ref.drug_id = *d_1;
					ref.cell_id = *c;
					
					if(opt.use_pair_prediction){
						
						// Only train on unique pairs of drugs (i.e. once we have trained on {A,B}, don't
						// train on {B, A}
						for(vector<string>::const_iterator d_2 = d_1 + 1;d_2 != cv_training_drugs.end();++d_2){

							MAP<string, vector<float> >::const_iterator drug_features_iter = drug_features.find(*d_2);

							if( drug_features_iter == drug_features.end() ){
								throw __FILE__ ":main: Unable to look up fingerprint data for selected training drug (2)";
							}

							if( num_features != drug_features_iter->second.size() ){
								throw __FILE__ ":main: Variable number of features detected";
							}
							
							MAP<string, float>::const_iterator col_iter = row_iter->second.find(*d_2);

							if( ( col_iter == row_iter->second.end() ) || (col_iter->second == MISSING_DATA) ){

								// Synergy matricies may be missing drug combinations
								continue;
							}
							
							train_pair.push_back( LabeledData() );

							LabeledData &ref = train_pair.back();
							
							// We will be predicting the actual synergy value for this pair
							ref.response = col_iter->second;
							
							ref.drug_id = *d_1;
							ref.drug_id_2 = *d_2;
							ref.cell_id = *c;
						}
					}
				}
			}
			
			// Train the machine learning models
			RandomForest rf(opt.num_tree, 
				opt.leaf_size, 
				opt.bag_fraction,
				1.0, // Don't bag the data
				&opt.seed);
			
			rf.build(train, cell_features, drug_features);
			
			RandomForest pair_rf(opt.num_tree, 
				opt.leaf_size, 
				opt.bag_fraction,
				opt.bag_fraction, // *Do* bag the data
				&opt.seed);
			
			if(opt.use_pair_prediction){
				pair_rf.build(train_pair, cell_features, drug_features);
			}
			
			double local_mse = 0.0; 	// For predicting single drug synergy
			double local_norm = 0.0;
			double local_average_mse = 0.0;
			
			double local_pair_mse = 0.0;	// For predicting drug pair synergy
			double local_pair_norm = 0.0;
			double local_pair_average_mse = 0.0;

			deque< pair<float, float> > local_actual_vs_predicted;
			
			// Can we predict pair synergy from single drug synergies?
			MAP<string /*cell*/, MAP<string /*drug*/, float> > test_drug_synergy;
			
			// We will only perform permutation testing when the user has requested randomization of
			// the test synergy and has provided an explicit test set
			const size_t num_permutation = 1 + 
				( (opt.randomize_test_synergy && explicit_test_provided) ? 2000 : 0);
			
			////////////////////
			// Testing
			////////////////////
			
			// For permutation testing of the enrichment results
			deque<float> permuted_fractional_enrichment;
			deque<float> permuted_single_auroc;
			deque<float> permuted_single_gini;
			
			deque<float> permuted_fractional_pair_enrichment;
			deque<float> permuted_pair_auroc;
			deque<float> permuted_pair_gini;
		
			for(size_t permutation = 0;permutation < num_permutation;++permutation){

				// Pack the testing data for regression. 
				vector<LabeledData> test;
				vector<LabeledData> test_pair;

				test.reserve( cv_testing_drugs.size() * cv_testing_cells.size() );

				if(opt.use_pair_prediction){

					// NumDrugs*(NumDrugs - 1)*NumCells/2
					test_pair.reserve( cv_testing_drugs.size() * (cv_testing_drugs.size() - 1)* 
						cv_testing_cells.size()/2 );
				}

				for(vector<string>::const_iterator c = cv_testing_cells.begin();c != cv_testing_cells.end();++c){

					MAP<string, MAP< string, MAP<string, float> > >::const_iterator synergy_iter;

					if(explicit_test_provided){

						synergy_iter = test_synergy_matrix.find(*c);

						if( synergy_iter == test_synergy_matrix.end() ){
							throw __FILE__ ":main: Unable to look up explicit test synergy matrix";
						}
					}
					else{
						synergy_iter = train_synergy_matrix.find(*c);

						if( synergy_iter == train_synergy_matrix.end() ){
							throw __FILE__ ":main: Unable to look up test value from training synergy matrix";
						}
					}

					for(vector<string>::const_iterator d_1 = cv_testing_drugs.begin();d_1 != cv_testing_drugs.end();++d_1){

						MAP<string, vector<float> >::const_iterator drug_features_iter = drug_features.find(*d_1);

						if( drug_features_iter == drug_features.end() ){
							throw __FILE__ ":main: Unable to look up fingerprint data for selected testing drug (1)";
						}

						MAP< string, MAP<string, float> >::const_iterator row_iter = synergy_iter->second.find(*d_1);

						if( row_iter == synergy_iter->second.end() ){

							// Synergy matricies may have missing drugs
							continue;
						}

						if( num_features != drug_features_iter->second.size() ){
							throw __FILE__ ":main: Variable number of testing features detected";
						}


						pair<float, float> score(0.0f, 0.0f);

						// Fully disjoint cross validation: Only include pairs of *test* drugs 
						// in the test set (*both* are absent from the training set)
						score = synergy_score(
								cv_testing_drugs, 
								row_iter->second,
								opt.score_func, 
								opt.synergy_threshold);

						if(score.second <= 0.0){

							// A score weight less than zero indicates too much missing data
							continue;
						}

						test.push_back( LabeledData() );

						LabeledData &ref = test.back();

						ref.response = score.first;
						ref.drug_id = *d_1;
						ref.cell_id = *c;

						if(opt.use_pair_prediction){

							// Only test on unique pairs of drugs (i.e. once we have tested {A,B}, don't
							// test {B, A}
							for(vector<string>::const_iterator d_2 = d_1 + 1;d_2 != cv_testing_drugs.end();++d_2){

								MAP<string, vector<float> >::const_iterator drug_features_iter = drug_features.find(*d_2);

								if( drug_features_iter == drug_features.end() ){
									throw __FILE__ ":main: Unable to look up fingerprint data for selected testing drug (2)";
								}

								if( num_features != drug_features_iter->second.size() ){
									throw __FILE__ ":main: Variable number of testing features detected";
								}

								MAP<string, float>::const_iterator col_iter = row_iter->second.find(*d_2);

								if( ( col_iter == row_iter->second.end() ) || (col_iter->second == MISSING_DATA) ){

									// Synergy matricies may have missing drug pairs
									continue;
								}

								test_pair.push_back( LabeledData() );

								LabeledData &ref = test_pair.back();

								// We will be predicting the actual synergy value for this pair
								ref.response = col_iter->second;

								ref.drug_id = *d_1;
								ref.drug_id_2 = *d_2;
								ref.cell_id = *c;
							}
						}
					}
				}

				if( num_features != drug_feature_id.size() ){

					cerr << "num_features = " << num_features << endl;
					cerr << "|drug_feature_id| = " << drug_feature_id.size() << endl;
					throw __FILE__ ":main: num_features != |drug_feature_id|";
				}

				if(permutation == 0){
					
					out << "\t|test| = " << test.size() << endl;
					out << "\t|train| = " << train.size() << endl;
				}
				
				if(opt.use_pair_prediction){

					if(permutation == 0){
					
						out << "\t|pair test| = " << test_pair.size() << endl;
						out << "\t|pair train| = " << train_pair.size() << endl;
					}
				}

				// Compute the *weighted* average response for the base line MSE
				double ave_response = 0.0;
				double norm = 0.0;

				for(vector<LabeledData>::const_iterator i = test.begin();i != test.end();++i){

					ave_response += i->response;
					norm += 1.0;
				}

				if(norm > 0.0){
					ave_response /= norm;
				}

				if(permutation == 0){
					out << "\tAverage test response = " << ave_response << endl;
				}

				double pair_ave_synergy = 0.0;
				norm = 0.0;

				for(vector<string>::const_iterator i = cv_testing_drugs.begin();i != cv_testing_drugs.end();++i){

					for(vector<string>::const_iterator j = cv_testing_cells.begin();j != cv_testing_cells.end();++j){

						MAP<string, MAP< string, MAP<string, float> > >::const_iterator synergy_iter;

						if(explicit_test_provided){

							synergy_iter = test_synergy_matrix.find(*j);

							if( synergy_iter == test_synergy_matrix.end() ){
								throw __FILE__ ":main: Unable to look up explicit test synergy matrix";
							}
						}
						else{
							synergy_iter = train_synergy_matrix.find(*j);

							if( synergy_iter == train_synergy_matrix.end() ){
								throw __FILE__ ":main: Unable to look up test value from training synergy matrix";
							}
						}

						MAP< string, MAP<string, float> >::const_iterator row_iter = synergy_iter->second.find(*i);

						if( row_iter == synergy_iter->second.end() ){

							// Synergy matricies may have missing drugs
							continue;
						}

						for(vector<string>::const_iterator j = i + 1;j != cv_testing_drugs.end();++j){

							MAP<string, float>::const_iterator col_iter = row_iter->second.find(*j);

							if( ( col_iter == row_iter->second.end() ) || (col_iter->second == MISSING_DATA) ){

								// Synergy matricies may have missing drugs
								continue;
							}

							pair_ave_synergy += col_iter->second;
							norm += 1.0;
						}
					}
				}

				if(norm > 0.0){
					pair_ave_synergy /= norm;
				}

				if(permutation == 0){
					out << "\tpair average test synergy = " << pair_ave_synergy << endl;
				}

				for(vector<LabeledData>::const_iterator i = test.begin();i != test.end();++i){

					MAP<string, vector<float> >::const_iterator cell_iter = 
						cell_features.find(i->cell_id);

					if( !cell_features.empty() && ( cell_iter == cell_features.end() ) ){
						throw __FILE__ ":main: Unable to look up cell iterator";
					}

					MAP<string, vector<float> >::const_iterator drug_iter = 
						drug_features.find(i->drug_id);

					if( !drug_features.empty() && ( drug_iter == drug_features.end() ) ){
						throw __FILE__ ":main: Unable to look up drug iterator";
					}

					float y = 0.0;
					
					if( cell_iter == cell_features.end() ){

						// No cell features
						if( drug_iter == drug_features.end() ){

							// No drug features
							y = rf.predict( vector<float>(), vector<float>() );
						}
						else{

							// Drug features provided
							y = rf.predict(vector<float>(), drug_iter->second);
						}
					}
					else{
						// Cell features are provided
						if( drug_iter == drug_features.end() ){

							// No drug features
							y = rf.predict( cell_iter->second, vector<float>() );
						}
						else{

							// Drug features provided
							y = rf.predict(cell_iter->second, drug_iter->second);
						}
					}

					actual_vs_predicted.push_back( make_pair(i->response, y) );
					local_actual_vs_predicted.push_back( make_pair(i->response, y) );

					test_drug_synergy[i->cell_id][i->drug_id] = y;

					float delta = (y - i->response);

					local_mse += delta*delta;

					delta = (ave_response - i->response);

					local_average_mse += delta*delta;

					local_norm += 1.0;
				}
				
				// Compute the correlation between predicted and observed drug *pair* 
				// synergy values for the test set
				deque< pair<float, float> > local_pair_synergy_comparison;

				// For measuring rank-based enrichment across all testing cell lines
				deque< 
					pair<float /*synergy score used for ranking*/, 
						bool /*actual synergy label*/> > local_synergy_rank;

				deque< 
					pair<float /*synergy score used for ranking*/, 
						bool /*actual synergy label*/> > local_pair_synergy_rank;

				for(vector<string>::const_iterator c = cv_testing_cells.begin();c != cv_testing_cells.end();++c){

					// For measuring rank-based enrichment for the current test cell line
					deque< 
						pair<float /*synergy score used for ranking*/, 
							bool /*actual synergy label*/> > local_cell_synergy_rank;

					deque< 
						pair<float /*synergy score used for ranking*/, 
							bool /*actual synergy label*/> > local_cell_pair_synergy_rank;

					MAP<string, MAP< string, MAP<string, float> > >::const_iterator synergy_iter;

					if(explicit_test_provided){

						synergy_iter = test_synergy_matrix.find(*c);

						if( synergy_iter == test_synergy_matrix.end() ){
							throw __FILE__ ":main: Unable to look up explicit test synergy matrix";
						}
					}
					else{
						synergy_iter = train_synergy_matrix.find(*c);

						if( synergy_iter == train_synergy_matrix.end() ){
							throw __FILE__ ":main: Unable to look up test value from training synergy matrix";
						}
					}

					MAP<string /*cell*/, MAP<string /*drug*/, float> >::const_iterator test_synergy_iter 
						= test_drug_synergy.find(*c);

					if( test_synergy_iter == test_drug_synergy.end() ){
						throw __FILE__ ":main: Unable to look up test drug synergy";
					}

					// Only needed for pair-based prediction
					MAP<string, vector<float> >::const_iterator cell_iter = cell_features.end();

					if(opt.use_pair_prediction){

						cell_iter = cell_features.find(*c);

						if( !cell_features.empty() && ( cell_iter == cell_features.end() ) ){
							throw __FILE__ ":main: Unable to look up cell iterator";
						}
					}

					for(vector<string>::const_iterator i = cv_testing_drugs.begin();i != cv_testing_drugs.end();++i){

						MAP< string, MAP<string, float> >::const_iterator row_iter = synergy_iter->second.find(*i);

						if( row_iter == synergy_iter->second.end() ){

							// Synergy matricies may have missing drugs
							continue;
						}

						MAP<string /*drug*/, float>::const_iterator synergy_1_iter = test_synergy_iter->second.find(*i);

						if( synergy_1_iter == test_synergy_iter->second.end() ){

							// Synergy matricies may have missing drugs
							continue;
							//throw __FILE__ ":main: Unable to find single drug synergy 1";
						}

						// Only needed for pair-based prediction
						MAP<string, vector<float> >::const_iterator drug_iter_1 = drug_features.end();

						if(opt.use_pair_prediction){

							drug_iter_1 = drug_features.find(*i);

							if( !drug_features.empty() && (drug_iter_1 == drug_features.end() ) ){
								throw __FILE__ ":main: Unable to look up drug iterator 1";
							}
						}

						// Don't compare a drug against itself and don't double count drug pairs
						for(vector<string>::const_iterator j = i + 1;j != cv_testing_drugs.end();++j){

							MAP<string, float>::const_iterator col_iter = row_iter->second.find(*j);

							if( ( col_iter == row_iter->second.end() ) || (col_iter->second == MISSING_DATA) ){

								// Synergy matricies may have missing drugs
								continue;
							}

							MAP<string /*drug*/, float>::const_iterator synergy_2_iter = test_synergy_iter->second.find(*j);

							if( synergy_2_iter == test_synergy_iter->second.end() ){

								// Synergy matricies may have missing drugs
								continue;
							}

							// Only needed for pair-based prediction
							MAP<string, vector<float> >::const_iterator drug_iter_2 = drug_features.end();

							if(opt.use_pair_prediction){

								drug_iter_2 = drug_features.find(*j);

								if( !drug_features.empty() && ( drug_iter_2 == drug_features.end() ) ){
									throw __FILE__ ":main: Unable to look up drug iterator 2";
								}
							}

							const bool is_synergistic = (col_iter->second <= opt.synergy_threshold);

							if(opt.score_func == COUNTING_SCORE){

								// **************************************************************************
								// The predicted synergy is defined here as **minus** the product of the single 
								// drug predicted probabilities. Why minus you ask?
								// The minus sign allows us to rank single drug and drug pair prediction with
								// the same code when computing the fractional enrichment (more negative means
								// a stronger prediction of synergy).
								// **************************************************************************
								const float predicted_synergy = -synergy_1_iter->second*synergy_2_iter->second;

								local_synergy_rank.push_back( make_pair(predicted_synergy, is_synergistic) );
								local_cell_synergy_rank.push_back( make_pair(predicted_synergy, is_synergistic) );
							}
							else{
								// A quick test (using cell line "786-0") suggests that the {independent, average, strongest} 
								// methods for combining single drug synergy values yeild similar correlation values:
								// [method] 	[pearson] 	[spearman]	[mse]
								// independent	0.257963	0.224935	-97.231
								// ave		0.247426	0.212759	5.55936
								// strongest	0.215334	0.175944	3.96274
								// weakest	0.232046	0.204585	4.27527
								//
								// However, the independent method does not do a good job predicting absolute
								// synergy values.

								// Assume that the individual drugs act independently
								//const float predicted_synergy = synergy_1_iter->second + synergy_2_iter->second + 
								//	synergy_1_iter->second*synergy_2_iter->second;

								// ... strongest effect dominates
								//const float predicted_synergy = min(synergy_1_iter->second, synergy_2_iter->second);

								// ... weakest effect dominates
								//const float predicted_synergy = max(synergy_1_iter->second, synergy_2_iter->second);

								// ... average effect
								//const float predicted_synergy = 0.5*(synergy_1_iter->second + synergy_2_iter->second);

								// Geometric mean (thanks to Nick Hengartner for pointing this out!)
								const float predicted_synergy = 
									-sqrt( min(0.0f, synergy_1_iter->second)*min(0.0f, synergy_2_iter->second) );

								local_pair_synergy_comparison.push_back( make_pair(col_iter->second, predicted_synergy) );

								pair_synergy_comparison.push_back( make_pair(col_iter->second, predicted_synergy) );

								local_synergy_rank.push_back( make_pair(predicted_synergy, is_synergistic) );
								local_cell_synergy_rank.push_back( make_pair(predicted_synergy, is_synergistic) );

								float delta = (col_iter->second - predicted_synergy);

								local_pair_mse += delta*delta;

								delta = (col_iter->second - pair_ave_synergy);

								local_pair_average_mse += delta*delta;

								local_pair_norm += 1.0;
							}

							if(opt.use_pair_prediction){

								float pair_synergy = 0.0;

								if( cell_iter == cell_features.end() ){

									// No cell features
									if( drug_iter_1 == drug_features.end() ){

										// No drug features
										pair_synergy = pair_rf.predict( vector<float>(), 
											vector<float>(), vector<float>() );
									}
									else{

										// Drug features provided
										pair_synergy = pair_rf.predict(vector<float>(), 
											drug_iter_1->second, drug_iter_2->second);
									}
								}
								else{
									// Cell features are provided
									if( drug_iter_1 == drug_features.end() ){

										// No drug features
										pair_synergy = pair_rf.predict( cell_iter->second, 
											vector<float>(), vector<float>() );
									}
									else{

										// Drug features provided
										pair_synergy = pair_rf.predict(cell_iter->second, 
											drug_iter_1->second, drug_iter_2->second);
									}
								}

								local_pair_synergy_rank.push_back( make_pair(pair_synergy, is_synergistic) );
								local_cell_pair_synergy_rank.push_back( make_pair(pair_synergy, is_synergistic) );
							}
						}
					}

					if(permutation == 0){

						const float local_enrichment = compute_enrichment(local_cell_synergy_rank);
						const float local_auroc = compute_AUROC(local_cell_synergy_rank);
						const float local_gini = compute_gini(local_auroc);
						
						// Compute the enrichment for this test cell line
						out << "\tPer test cell line single drug-based enrichment|AUROC|Gini: " 
							<< *c << '\t' 
							<< local_enrichment << '\t'
							<< local_auroc << '\t'
							<< local_gini << endl;
													
						if(opt.use_pair_prediction){

							const float local_enrichment = compute_enrichment(local_cell_pair_synergy_rank);
							const float local_auroc = compute_AUROC(local_cell_pair_synergy_rank);
							const float local_gini = compute_gini(local_auroc);
							
							out << "\tPer test cell line drug pair-based enrichment|AUROC|Gini: " 
								<< *c << '\t' 
								<< local_enrichment << '\t'
								<< local_auroc << '\t'
								<< local_gini << endl;
						}
					}
				}

				if(permutation == 0){

					mse += local_mse;
					average_mse += local_average_mse;
					mse_norm += local_norm;

					pair_mse += local_pair_mse;
					pair_average_mse += local_pair_average_mse;
					pair_mse_norm += local_pair_norm;

					if(local_norm > 0.0){

						local_mse /= local_norm;
						local_average_mse /= local_norm;
					}

					if(local_pair_norm > 0.0){

						local_pair_mse /= local_pair_norm;
						local_pair_average_mse /= local_pair_norm;
					}

					out << "\tMSE = " << local_mse << endl;
					out << "\tAverage MSE = " << local_average_mse << endl;

					out << "\tpearson r = " << pearson_correlation(local_actual_vs_predicted) << endl;
					out << "\tspearman r = " << spearman_correlation(local_actual_vs_predicted) << endl;

					if(opt.use_pair_prediction){

						const float local_enrichment = compute_enrichment(local_pair_synergy_rank);
						const float local_auroc = compute_AUROC(local_pair_synergy_rank);
						const float local_gini = compute_gini(local_auroc);
						
						fractional_pair_enrichment.push_back(local_enrichment);
						pair_auroc.push_back(local_auroc);
						pair_gini.push_back(local_gini);
						
						out << "\tPair synergy enrichment|AUROC|Gini = " 
							<< local_enrichment << '\t'
							<< local_auroc << '\t'
							<< local_gini << endl;
						
						if(mpi_rank == 0){	
						
							cerr << "Fold " << (fold + 1) 
								<< ": drug pair enrichment|AUROC|Gini = " 
								<< local_enrichment << '\t'
								<< local_auroc << '\t'
								<< local_gini << endl;
						}
					}

					if(opt.score_func == COUNTING_SCORE){

						const float local_enrichment = compute_enrichment(local_synergy_rank);
						const float local_auroc = compute_AUROC(local_synergy_rank);
						const float local_gini = compute_gini(local_auroc);
						
						fractional_enrichment.push_back(local_enrichment);
						single_auroc.push_back(local_auroc);
						single_gini.push_back(local_gini);
						
						out << "\tsingle drug synergy enrichment|AUROC|Gini = " 
							<< local_enrichment << '\t'
							<< local_auroc << '\t'
							<< local_gini << endl;
						
						if(mpi_rank == 0){
						
							cerr << "Fold " << (fold + 1) 
								<< ": single drug enrichment|AUROC|Gini = " 
								<< local_enrichment << '\t'
								<< local_auroc << '\t'
								<< local_gini << endl;
						}
						
						#ifdef WRITE_ENRICHMENT_CURVE
						// Print the detailed enrichment curve 
						if( fout.is_open() ){

							// Rank our predictions in order of *ascending* predicted synergy values
							// (more negative means more synergistic).
							sort( local_synergy_rank.begin(), local_synergy_rank.end() );

							// For the random (i.e. null model), the expected number of synergistic pairs
							// that would be randomly selected in a subset of M pairs is M*S/N,
							// where N is the total number of pairs and S is the total number of synergistic pairs.

							const size_t total_num_pair = local_synergy_rank.size();
							size_t total_num_synergistic = 0;

							for(size_t i = 0;i < total_num_pair;++i){
								total_num_synergistic += local_synergy_rank[i].second ? 1 : 0;
							}

							fout << "\tSelection size\tNumber Synergistic\tExpected Synergistic" << endl;

							// Start at one (since we always need to have at least one prediction)
							for(size_t i = 1;i <= total_num_pair;++i){

								// Count the number of correctly predicted synergistic pairs
								size_t num_synergistic = 0;

								for(size_t j = 0;j < i;++j){
									num_synergistic += (local_synergy_rank[j].second == true) ? 1 : 0;				
								}

								const float expected_num_synergistic = float(i*total_num_synergistic)/total_num_pair;

								fout << "ENRICHMENT[" << fold << "]\t" 
									<< float(i)/total_num_pair 						// Fractional selection size
									<< '\t' << float(num_synergistic)/total_num_synergistic 		// Fraction predicted synergistic
									<< '\t' << float(expected_num_synergistic)/total_num_synergistic	// Expected Synergistic
									<< endl;
							}
						}
						#endif // WRITE_ENRICHMENT_CURVE
					}
					else{

						out << "\tpair MSE = " << local_pair_mse << endl;
						out << "\tpair average MSE = " << local_pair_average_mse << endl;

						out << "\tpair synergy correlation (pearson) = " << pearson_correlation(local_pair_synergy_comparison) << endl;
						out << "\tpair synergy correlation (spearman) = " << spearman_correlation(local_pair_synergy_comparison) << endl;
					}

					out << endl;
				}
				else{ // permutation > 0
					
					if(opt.use_pair_prediction){

						const float local_enrichment = compute_enrichment(local_pair_synergy_rank);
						const float local_auroc = compute_AUROC(local_pair_synergy_rank);
						const float local_gini = compute_gini(local_auroc);
						
						permuted_fractional_pair_enrichment.push_back(local_enrichment);
						permuted_pair_auroc.push_back(local_auroc);
						permuted_pair_gini.push_back(local_gini);
						
						if(mpi_rank == 0){
						
							cerr << "\tPermutation " << permutation 
								<< ": drug pair enrichment|AUROC|Gini = " 
								<< local_enrichment << '\t'
								<< local_auroc << '\t'
								<< local_gini << endl;
						}
					}

					if(opt.score_func == COUNTING_SCORE){

						const float local_enrichment = compute_enrichment(local_synergy_rank);
						const float local_auroc = compute_AUROC(local_synergy_rank);
						const float local_gini = compute_gini(local_auroc);
						
						permuted_fractional_enrichment.push_back(local_enrichment);
						permuted_single_auroc.push_back(local_auroc);
						permuted_single_gini.push_back(local_gini);
						
						if(mpi_rank == 0){
						
							cerr << "\tPermutation " << permutation 
								<< ": single drug enrichment|AUROC|Gini = " 
								<< local_enrichment << '\t'
								<< local_auroc << '\t'
								<< local_gini << endl;
						}
					}
				}
				
				// Randomize the test synergy matricies for permutation testing
				if(num_permutation > 1){

					if(mpi_rank == 0){
						cerr << "Permuting test matricies (permutation round " << permutation << ")" << endl;
					}
					
					randomize_dependent(test_synergy_matrix, &opt.seed);
					
					// rank 0 defines the ordering
					broadcast(test_synergy_matrix, mpi_rank, 0);
				}
			}
			
			// If we were performing permutation testing, then we need to output the enrichment p-values
			if(num_permutation > 1){
				
				out << "Performed " << num_permutation 
					<< " rounds of permutation testing for the fractional synergy enrichment" << endl;
				
				if( !fractional_enrichment.empty() ){
					
					out << "Fractional single drug-based synergy enrichment p-value (two-tailed) = " 
						<< compute_p_value(fractional_enrichment.back(), permuted_fractional_enrichment) 
						<< endl;
				}
				
				if( !single_auroc.empty() ){
					
					out << "Single drug-based synergy AUROC p-value (two-tailed) = " 
						<< compute_p_value(single_auroc.back(), permuted_single_auroc) 
						<< endl;
				}
				
				if( !single_gini.empty() ){
					
					out << "Single drug-based synergy Gini p-value (two-tailed) = " 
						<< compute_p_value(single_gini.back(), permuted_single_gini) 
						<< endl;
				}
				
				if( !fractional_pair_enrichment.empty() ){
					
					out << "Fractional drug pair-based synergy enrichment p-value (two-tailed) = " 
						<< compute_p_value(fractional_pair_enrichment.back(), permuted_fractional_pair_enrichment) 
						<< endl;
				}
				
				if( !pair_auroc.empty() ){
					
					out << "Fractional drug pair-based synergy AUROC p-value (two-tailed) = " 
						<< compute_p_value(pair_auroc.back(), permuted_pair_auroc) 
						<< endl;
				}
				
				if( !pair_gini.empty() ){
					
					out << "Fractional drug pair-based synergy Gini p-value (two-tailed) = " 
						<< compute_p_value(pair_gini.back(), permuted_pair_gini) 
						<< endl;
				}
				
				out << endl;
			}
		}
		
		if(opt.num_fold >= 1){
		
			if(mse_norm > 0.0){
			
				mse /= mse_norm;
				average_mse /= mse_norm;
			}
			
			if(pair_mse_norm > 0.0){
			
				pair_mse /= pair_mse_norm;
				pair_average_mse /= pair_mse_norm;
			}
			
			out << "Final MSE = " << mse << endl;
			out << "Final Average MSE = " << average_mse << endl;

			out << "Final improvement over average MSE = " << 100.0*(average_mse - mse)/average_mse << " %" << endl;
			out << "Final pearson r = " << pearson_correlation(actual_vs_predicted) << endl;
			out << "Final spearman r = " << spearman_correlation(actual_vs_predicted) << endl;
			out << endl;
			
			if(opt.use_pair_prediction){
				
				pair<float, float> info = average_and_stdev(fractional_pair_enrichment);
				
				out << "Final fractional drug pair-based synergy enrichment = " 
					<< info.first << " +/- " << info.second << endl;
					
				info = average_and_stdev(pair_auroc);
				
				out << "Final drug pair-based synergy AUROC = " 
					<< info.first << " +/- " << info.second << endl;
					
				info = average_and_stdev(pair_gini);
				
				out << "Final drug pair-based synergy Gini = " 
					<< info.first << " +/- " << info.second << endl;
			}
			
			if(opt.score_func == COUNTING_SCORE){
				
				pair<float, float> info = average_and_stdev(fractional_enrichment);
				
				out << "Final fractional single drug-based synergy enrichment = " 
					<< info.first << " +/- " << info.second << endl;
					
				info = average_and_stdev(single_auroc);
				
				out << "Final single drug-based synergy AUROC = " 
					<< info.first << " +/- " << info.second << endl;
					
				info = average_and_stdev(single_gini);
				
				out << "Final single drug-based synergy Gini = " 
					<< info.first << " +/- " << info.second << endl;
			}
			else{
				
				out << "Final pair MSE = " << pair_mse << endl;
				out << "Final pair average MSE = " << pair_average_mse << endl;
				out << "Final improvement over pair average MSE = " << 100.0*(pair_average_mse - pair_mse)/pair_average_mse << " %" << endl;

				out << "Final pair synergy correlation (pearson) = " << pearson_correlation(pair_synergy_comparison) << endl;
				out << "Final pair synergy correlation (spearman) = " << spearman_correlation(pair_synergy_comparison) << endl;
			}
		}
		
		if( !blind_test_drug_features.empty() ){
			
			// Pack the training data for regression. 
			vector<LabeledData> train;
			size_t num_features = 0;
			
			train.reserve( training_cells.size() * drugs.size() );
			
			for(vector<string>::const_iterator c = training_cells.begin();c != training_cells.end();++c){
			
				MAP<string, MAP< string, MAP<string, float> > >::const_iterator synergy_iter = train_synergy_matrix.find(*c);

				if( synergy_iter == train_synergy_matrix.end() ){
					throw __FILE__ ":main: Unable to look up cell in synergy matrix";
				}

				for(vector<string>::const_iterator i = drugs.begin();i != drugs.end();++i){

					MAP<string, vector<float> >::const_iterator drug_features_iter = drug_features.find(*i);

					if( drug_features_iter == drug_features.end() ){
						throw __FILE__ ":main: Unable to look up fingerprint data for selected training drug";
					}

					MAP< string, MAP<string, float> >::const_iterator row_iter = synergy_iter->second.find(*i);

					if( row_iter == synergy_iter->second.end() ){
					
						// Synergy matricies may have missing drugs
						continue;
					}

					if( i == drugs.begin() ){
						num_features = drug_features_iter->second.size();
					}

					if( num_features != drug_features_iter->second.size() ){
						throw __FILE__ ":main: Variable number of features detected";
					}

					// Only use the TRAINING drugs to compute the synergy score.
					// The synergy score can be:
					// (a) fraction of  synergistic interactions
					// (b) median synergy value (over all drug pairs that contain this drug)
					const pair<float, float> score = synergy_score(
						drugs, 
						row_iter->second, 
						opt.score_func, 
						opt.synergy_threshold);

					if(score.second <= 0.0){

						// A score weight less than zero indicates too much missing data
						continue;
					}

					train.push_back( LabeledData() );

					LabeledData &ref = train.back();

					ref.response = score.first;
					ref.drug_id = *i;
					ref.cell_id = *c;
				}
			}
			
			out << "\t|blind test| = " << blind_test_drug_features.size() << endl;
			out << "\t|train| = " << train.size() << endl;
			out << endl;
			
			RandomForest rf(opt.num_tree, 
				opt.leaf_size, 
				opt.bag_fraction, 
				1.0, // Don't bag the data
				&opt.seed);
			
			rf.build(train, cell_features, drug_features);
					
			for(vector<string>::const_iterator i = blind_test_prefix.begin();
				i != blind_test_prefix.end();++i){
			
				typedef MULTIMAP<string, string>::const_iterator I;
				
				const pair<I, I> range = blind_test_cell_prefix_map.equal_range(*i);
				
				size_t norm = 0;
					
				for(I k = range.first;k != range.second;++k){
					++norm;
				}
					
				deque<SynergyPrediction> blind_test_drug_synergy;
				
				for(MAP<string, vector<float> >::const_iterator j = blind_test_drug_features.begin();
					j != blind_test_drug_features.end();++j){
					
					// Average over all of the blind test cells that share the same
					// prefix
					SynergyPrediction pred;
					
					for(I k = range.first;k != range.second;++k){
						
						MAP<string, vector<float> >::const_iterator iter = 
							blind_test_cell_features.find(k->second);
						
						if( iter == blind_test_cell_features.end() ){
							throw __FILE__ ":main: Unable to lookup blind test cell features";
						}
						
						const float y = rf.predict(iter->second, j->second);
						
						pred.ave += y;
						pred.stdev += y*y;
					}
					
					if(norm > 1){
						
						pred.ave /= norm;
						pred.stdev = sqrt(pred.stdev/(norm - 1.0) - pred.ave*pred.ave);
					}
					
					pred.drug_id = j->first;
					
					// Store the per-drug synergy values for the current cell
					blind_test_drug_synergy.push_back(pred);
				}
				
				// Sort by single drug synergy propensity
				sort( blind_test_drug_synergy.begin(), blind_test_drug_synergy.end() );

				// Sort in descending order
				reverse( blind_test_drug_synergy.begin(), blind_test_drug_synergy.end() );
				
				fout << "Predicted single drug synergy propensities for " << norm 
					<< " blind test cells belonging to: " << *i << endl;
				
				for(deque<SynergyPrediction>::const_iterator j = blind_test_drug_synergy.begin();
					j != blind_test_drug_synergy.end();++j){

					fout << j->drug_id << '\t' << j->ave << '\t' << j->stdev << endl;
				}

				fout << endl;
				
				if(num_blind_with_target > 0){
					
					const size_t num_best = 50;
					
					fout << "Best " << num_best << " blind predictions for drug pairs with target information" << endl;
					
					deque< pair< float, pair<string, string> > > drug_pairs;
					
					for(deque<SynergyPrediction>::const_iterator j = blind_test_drug_synergy.begin();
						j != blind_test_drug_synergy.end();++j){
						
						MAP<string, string>::const_iterator target_j_iter = drug_target.find(j->drug_id);
						
						if( target_j_iter == drug_target.end() ){
							continue;
						}
						
						// Don't combine a drug with itself
						for(deque<SynergyPrediction>::const_iterator k = j + 1;
							k != blind_test_drug_synergy.end();++k){
							
							MAP<string, string>::const_iterator target_k_iter = drug_target.find(k->drug_id);

							if( target_k_iter == drug_target.end() ){
								continue;
							}
							
							// Don't combine two drugs that share the same target
							if(target_j_iter->second == target_k_iter->second){
								continue;
							}
							
							const float score = sqrt(j->ave*k->ave);
							
							drug_pairs.push_back( make_pair( score, make_pair(j->drug_id, k->drug_id) ) );
						}
					}
					
					// Sort by drug pair synergy propensity
					sort( drug_pairs.begin(), drug_pairs.end() );
					
					// Sort in descending order
					reverse( drug_pairs.begin(), drug_pairs.end() );
					
					fout << "Score\tDrug 1 ID\tDrug 2 ID\tDrug 1 target\tDrug 2 target" << endl;
					
					size_t index = 0;
					
					for(deque< pair< float, pair<string, string> > >::const_iterator j = drug_pairs.begin();
						( j != drug_pairs.end() ) && (num_best > index) ;++j,++index){
						
						MAP<string, string>::const_iterator target_1_iter = drug_target.find(j->second.first);
						
						if( target_1_iter == drug_target.end() ){
							throw __FILE__ ":main: Unable to lookup target 1";
						}
						
						MAP<string, string>::const_iterator target_2_iter = drug_target.find(j->second.second);
						
						if( target_2_iter == drug_target.end() ){
							throw __FILE__ ":main: Unable to lookup target 2";
						}
						
						fout << j->first << '\t' 
							<< j->second.first << '\t' 
							<< j->second.second << '\t'
							<< target_1_iter->second << '\t'
							<< target_2_iter->second << endl;
					}
					
					fout << endl;
				}
				
				#ifdef OUTPUT_SYNERGY_MATRIX
				fout << "Predicted synergy matrix for blind test set: " << i->first << endl;

				// Column headers
				for(deque< pair<float, string> >::const_iterator j = blind_test_drug_synergy.begin();
					j != blind_test_drug_synergy.end();++j){
					fout << '\t' << j->second;
				}

				fout << endl;
				
				for(deque< pair<float, string> >::const_iterator j = blind_test_drug_synergy.begin();
					j != blind_test_drug_synergy.end();++j){

					// Row header
					fout << j->second;

					for(deque< pair<float, string> >::const_iterator k = blind_test_drug_synergy.begin();
						k != blind_test_drug_synergy.end();++k){

						fout << '\t';

						if(j != k){

							const float predicted_synergy = sqrt(j->first*k->first);

							fout << predicted_synergy;
						}
					}

					fout << endl;
				}
				#endif // OUTPUT_SYNERGY_MATRIX
			}
		}
		
		profile = time(NULL) - profile;
		
		if(opt.num_fold > 0){
			out << "Cross validation complete in " << profile << " sec" << endl;
		}
		
		MPI_Finalize();
	}
	catch(const char *error){
		cerr << "[" << mpi_rank << "] Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(const string error){
		cerr << "[" << mpi_rank << "] Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){
		cerr << "[" << mpi_rank << "] Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}
	
	
	return EXIT_SUCCESS;
}

float compute_distance(const vector<float> &m_a, const vector<float> &m_b)
{
	const size_t len = m_a.size();
	
	if( len != m_b.size() ){
		throw __FILE__ ":compute_distance: Size mismatch!";
	}
	
	float ret = 0.0;
	size_t overlap = 0;
	
	for(size_t i = 0;i < len;++i){
		
		if( (m_a[i] == MISSING_DATA) || (m_b[i] == MISSING_DATA) ){
			continue;
		}
		
		++overlap;
		
		const float delta = m_a[i] - m_b[i];
		
		ret += delta*delta;
	}
	
	if(overlap == 0){
		throw __FILE__ ":compute_distance: Unable to compute distance between non-overlapping vectors";
	}
	
	return sqrt(ret);
}

pair<float /*value*/, float /*weight*/> synergy_score(const vector<string> &m_drugs, const MAP<string, float> &m_synergy,
	const ScoreFunction &m_func, const float &m_synergy_threshold)
{
	if(m_func == COUNTING_SCORE){

		size_t num_synergistic = 0;
		size_t num_non_synergistic = 0;

		for(vector<string>::const_iterator i = m_drugs.begin();i != m_drugs.end();++i){

			MAP<string, float>::const_iterator iter = m_synergy.find(*i);

			if( iter == m_synergy.end() ){
				continue; // skip missing data
				//throw __FILE__ ":synergy_score: Unable to find query drug";
			}
			
			// Skip missing data
			if(iter->second == MISSING_DATA){
				continue;
			}
			
			if(iter->second <= m_synergy_threshold){
				++num_synergistic;
			}
			else{
				++num_non_synergistic;
			}
		}

		if( (num_synergistic + num_non_synergistic) == 0 ){

			// Use a score weight less than zero to indicate too much missing data
			return make_pair(0.0f, -1.0f);
			//throw __FILE__ ":synergy_score: Unable to find a single valid drug pair to count!";
		}

		return make_pair( float(num_synergistic)/(num_synergistic + num_non_synergistic),
			1.0); // Equal weights for now
	}
	
	deque<float> data;
	
	for(vector<string>::const_iterator i = m_drugs.begin();i != m_drugs.end();++i){

		MAP<string, float>::const_iterator iter = m_synergy.find(*i);

		if( iter == m_synergy.end() ){
			throw __FILE__ ":synergy_score: Unable to find query drug";
		}

		if(iter->second != MISSING_DATA){
			data.push_back(iter->second);
		}
	}

	const size_t N = data.size();
	
	if(N == 0){
		throw __FILE__ ":synergy_score: Unable to compute a score (no valid data points!)";
	}
	
	float ret = 0.0;
	float w = 0.0;
	
	switch(m_func){
		case AVERAGE_SCORE:
			
			for(deque<float>::const_iterator i = data.begin();i != data.end();++i){
				ret += *i;
			}

			ret /= N;
			w = sqrt(N);
			break;
		case MEDIAN_SCORE:

			sort( data.begin(), data.end() );

			ret = (N%2 == 1) ? data[N/2] // Odd number of points
				: 0.5*(data[N/2 - 1] + data[N/2]); // Even number of points

			w = 1.0; // Equal weights for now
			break;
		case FILTERED_MEAN_SCORE:
			
			if(N < 3){
				throw __FILE__ ":synergy_score: Unable to compute a filtered mean score (not enough valid data points!)";
			}

			// Sort to compute the filtererd mean
			sort( data.begin(), data.end() );

			for(size_t i = 1;i < (N - 1);++i){
				ret += data[i];
			}

			ret /= (N - 2);
			w = sqrt(N - 2); // <-- How should we actually weight the filtered mean?
			break;
		default:
			throw __FILE__ ":synergy_score: Unknown synergy scoring function";
	};
	
	return make_pair( ret, w );
}

// Convert the drug names from NSC to CID
MAP< string, MAP<string, float> > NSC_to_CID(const MAP< string, MAP<string, float> > &m_synergy_matrix)
{
	MAP< string, MAP<string, float> > ret;
	
	for(MAP< string, MAP<string, float> >::const_iterator i = m_synergy_matrix.begin();i != m_synergy_matrix.end();++i){
		
		const string name_i = NSC_to_CID(i->first);
		
		for(MAP<string, float>::const_iterator j = i->second.begin();j != i->second.end();++j){
		
			string name_j = NSC_to_CID(j->first);
			
			ret[name_i][name_j] = j->second;
		}
	}
	
	return ret;
}

string NSC_to_CID(const string &m_nsc)
{
	if(m_nsc == "740") return "126941";
	if(m_nsc == "750") return "2478";
	if(m_nsc == "752") return "2723601";
	if(m_nsc == "755") return "667490";
	if(m_nsc == "762") return "5935";
	if(m_nsc == "1390") return "2094";
	if(m_nsc == "3053") return "44415057";
	if(m_nsc == "3088") return "2708";
	if(m_nsc == "6396") return "5453";
	if(m_nsc == "8806") return "9927978";
	if(m_nsc == "9706") return "5799";
	if(m_nsc == "13875") return "2123";
	if(m_nsc == "14229") return "6239";
	if(m_nsc == "18509") return "123608";
	if(m_nsc == "19893") return "3385";
	if(m_nsc == "24559") return "163659";
	if(m_nsc == "25154") return "4842";
	if(m_nsc == "26271") return "2907";
	if(m_nsc == "26980") return "5746";
	if(m_nsc == "27640") return "5790";
	if(m_nsc == "32065") return "3657";
	if(m_nsc == "34462") return "6194";
	if(m_nsc == "38721") return "4211";
	if(m_nsc == "45388") return "5353562";
	if(m_nsc == "45923") return "4114";
	if(m_nsc == "49842") return "241902";
	if(m_nsc == "63878") return "6252";
	if(m_nsc == "66847") return "5426";
	if(m_nsc == "67574") return "249332";
	if(m_nsc == "71423") return "11683";
	if(m_nsc == "77213") return "4915";
	if(m_nsc == "79037") return "3950";
	if(m_nsc == "82151") return "30323";
	if(m_nsc == "85998") return "23615975";
	if(m_nsc == "92859") return "518740";
	if(m_nsc == "102816") return "1805";
	if(m_nsc == "105014") return "1546";
	if(m_nsc == "109724") return "3690";
	if(m_nsc == "118218") return "3367";
	if(m_nsc == "119875") return "5702198";
	if(m_nsc == "122758") return "444795";
	if(m_nsc == "122819") return "54610154";
	if(m_nsc == "123127") return "32874";
	if(m_nsc == "125066") return "54608728";
	if(m_nsc == "125973") return "36314";
	if(m_nsc == "127716") return "16886";
	if(m_nsc == "138783") return "77082";
	if(m_nsc == "141540") return "439525";
	if(m_nsc == "169780") return "71384";
	if(m_nsc == "180973") return "2733525";
	if(m_nsc == "218321") return "40926";
	if(m_nsc == "226080") return "54600319";
	if(m_nsc == "241240") return "426756";
	if(m_nsc == "246131") return "41744";
	if(m_nsc == "256439") return "3685";
	if(m_nsc == "256942") return "32874";
	if(m_nsc == "266046") return "24197464";
	if(m_nsc == "279836") return "4212";
	if(m_nsc == "296961") return "2141";
	if(m_nsc == "362856") return "5394";
	if(m_nsc == "369100") return "57469";
	if(m_nsc == "409962") return "2578";
	if(m_nsc == "606869") return "354624";
	if(m_nsc == "608210") return "25136944";
	if(m_nsc == "609699") return "60699";
	if(m_nsc == "613327") return "356653";
	if(m_nsc == "628503") return "148124";
	if(m_nsc == "673596") return "104842";
	if(m_nsc == "681239") return "387447";
	if(m_nsc == "686673") return "3011155";
	if(m_nsc == "698037") return "446556";
	if(m_nsc == "701852") return "5311";
	if(m_nsc == "702294") return "54611422";
	if(m_nsc == "707389") return "54611489";
	if(m_nsc == "712807") return "400633";
	if(m_nsc == "713563") return "60198";
	if(m_nsc == "715055") return "123631";
	if(m_nsc == "718781") return "176871";
	if(m_nsc == "719276") return "104741";
	if(m_nsc == "719344") return "2187";
	if(m_nsc == "719345") return "3902";
	if(m_nsc == "719627") return "2662";
	if(m_nsc == "721517") return "68740";
	if(m_nsc == "732517") return "3062316";
	if(m_nsc == "733504") return "54608520";
	if(m_nsc == "737754") return "11525740";
	if(m_nsc == "743414") return "5291";
	if(m_nsc == "745750") return "208908";
	if(m_nsc == "747599") return "644241";
	if(m_nsc == "747971") return "216239";
	if(m_nsc == "747972") return "216326";
	if(m_nsc == "747973") return "6445540";
	if(m_nsc == "747974") return "5035";
	if(m_nsc == "749226") return "132971";
	if(m_nsc == "750690") return "5329102";
	if(m_nsc == "753082") return "42611257";
	if(m_nsc == "754143") return "5352062";
	if(m_nsc == "754230") return "148121";
	if(m_nsc == "755986") return "24776445";
	if(m_nsc == "756645") return "54613769";
	if(m_nsc == "757441") return "6450551";
	if(m_nsc == "760766") return "3081361";
	if(m_nsc == "761431") return "42611257";
	if(m_nsc == "761432") return "9854073";
	if(m_nsc == "763371") return "25126798";

	cerr << "Did not find a match for putative NSC value: " << m_nsc << endl;
	throw __FILE__ ":NSC_to_CID: Error mapping NSC to CID string";
	return string();
}

void randomize_independent(MAP<string, vector<float> > &m_data, unsigned int *m_ptr_seed)
{
	// The independent variables (i.e. features) are stored in an associative array:
	// drug name -> feature vector. Shuffle this array so that the drug
	// name to feature vector mapping is randomized.
	
	const size_t num_drugs = m_data.size();
	
	vector<string> drug_names;
	
	drug_names.reserve(num_drugs);
	
	for(MAP<string, vector<float> >::const_iterator i = m_data.begin();i != m_data.end();++i){
		drug_names.push_back(i->first);
	}
	
	if(drug_names.size() != num_drugs){
		throw __FILE__ ":randomize_independent: Did not find the expected number of drug names";
	}
	
	vector<string> shuffled_drug_names = drug_names;
	
	randomize(shuffled_drug_names.begin(), shuffled_drug_names.end(), m_ptr_seed);
	
	// The new drug name to feature vector mapping
	MAP<string, vector<float> > ret;
	
	for(size_t i = 0;i < num_drugs;++i){
		
		MAP<string, vector<float> >::const_iterator iter = m_data.find(drug_names[i]);
		
		if( iter == m_data.end() ){
			throw __FILE__ ":randomize_independent: Unable to find drug name!";
		}
		
		ret[ shuffled_drug_names[i] ] = iter->second;
	}
	
	m_data = ret;
}

void randomize_dependent(MAP<string, MAP< string, MAP<string, float> > > &m_data, unsigned int *m_ptr_seed)
{
	// The dependent variables (i.e. synergy matrix elements) are stored in a symmetric, 
	// 3D associative array: cell name -> drug name row -> drug name col -> synergy value. 
	// Shuffle column and row labels of this 2D array so that the matrix is still
	// symmetric, but each drug pair is now randomly mapped to the synergy value of a 
	// different drug pair. In addition, a drug against the same drug (i.e. the diagonal
	// elements of the synergy matrix) are unchanged.
	
	const size_t num_cells = m_data.size();
	
	// The new drug name to feature vector mapping
	MAP<string, MAP< string, MAP<string, float> > > ret;
	
	vector<string> cell_names;
	
	cell_names.reserve(num_cells);
	
	for(MAP<string, MAP< string, MAP<string, float> > >::const_iterator i = m_data.begin();i != m_data.end();++i){
		cell_names.push_back(i->first);
	}
	
	vector<string> shuffled_cell_names = cell_names;
	
	randomize(shuffled_cell_names.begin(), shuffled_cell_names.end(), m_ptr_seed);
	
	for(size_t k = 0;k < num_cells;++k){

		MAP<string, MAP< string, MAP<string, float> > >::const_iterator src_iter = m_data.find(cell_names[k]);
		
		if( src_iter == m_data.end() ){
			throw __FILE__ ":randomize_dependent: Unable to look up synergy data by cell line name";
		}

		const size_t num_drugs = src_iter->second.size();

		vector<string> drug_names;

		drug_names.reserve(num_drugs);

		for(MAP< string, MAP<string, float> >::const_iterator i = src_iter->second.begin();i != src_iter->second.end();++i){
			drug_names.push_back(i->first);
		}

		if(drug_names.size() != num_drugs){
			throw __FILE__ ":randomize_dependent: Did not find the expected number of drug names";
		}

		vector<string> shuffled_drug_names = drug_names;

		randomize(shuffled_drug_names.begin(), shuffled_drug_names.end(), m_ptr_seed);

		for(size_t i = 0;i < num_drugs;++i){

			MAP< string, MAP<string, float> >::const_iterator row_iter = src_iter->second.find(drug_names[i]);

			if( row_iter == src_iter->second.end() ){
				throw __FILE__ ":randomize_dependent: Unable to find drug row name!";
			}

			for(size_t j = 0;j < num_drugs;++j){

				MAP<string, float>::const_iterator col_iter = row_iter->second.find(drug_names[j]);

				if( col_iter == row_iter->second.end() ){
					throw __FILE__ ":randomize_dependent: Unable to find drug column name!";
				}

				ret[ shuffled_cell_names[k] ][ shuffled_drug_names[i] ][ shuffled_drug_names[j] ] = col_iter->second;
			}
		}
	}
	
	m_data = ret;
}

pair<float, float> average_and_stdev(const deque<float> &m_data)
{
	pair<float, float> ret(0.0f, 0.0f);

	if( m_data.empty() ){
		return ret;
	}
	
	for(deque<float>::const_iterator i = m_data.begin();i != m_data.end();++i){
		ret.first += *i;
	}
	
	ret.first /= m_data.size();
	
	for(deque<float>::const_iterator i = m_data.begin();i != m_data.end();++i){
		
		const float delta = ret.first - *i;
		
		ret.second += delta*delta;
	}
	
	if(m_data.size() > 1){
		ret.second = sqrt( ret.second/(m_data.size() - 1) );
	}
	
	return ret;
}

double compute_enrichment(const deque< pair<float, bool> > &m_synergy_data)
{
	// Make a copy of the input data that we can sort
	deque< pair<float, bool> > local_synergy_data(m_synergy_data);
	
	double local_fractional_enrichment = 0.0; // Compared to the null model of random
	
	// Rank our predictions in order of *ascending* predicted synergy.
	sort( local_synergy_data.begin(), local_synergy_data.end() );

	// For the random (i.e. null model), the expected number of synergistic pairs
	// that would be randomly selected in a subset of M pairs is M*S/N,
	// where N is the total number of pairs and S is the total number of synergistic pairs.

	const size_t total_num_pair = local_synergy_data.size();
	size_t total_num_synergistic = 0;

	for(size_t i = 0;i < total_num_pair;++i){
		total_num_synergistic += local_synergy_data[i].second ? 1 : 0;
	}

	// Start at one (since we always need to have at least one prediction)
	for(size_t i = 1;i <= total_num_pair;++i){

		// Count the number of correctly predicted synergistic pairs
		size_t num_synergistic = 0;

		for(size_t j = 0;j < i;++j){
			num_synergistic += (local_synergy_data[j].second == true) ? 1 : 0;						
		}

		const double expected_num_synergistic = float(i*total_num_synergistic)/total_num_pair;

		local_fractional_enrichment += float(num_synergistic) - float(expected_num_synergistic);
	}

	// Normalize the local_fractional_enrichment
	if(total_num_synergistic > 0){
		local_fractional_enrichment /= 0.5*total_num_synergistic*total_num_pair;					
	}

	return local_fractional_enrichment;
}

// Compute the area under the reciever operator characteristic curve
double compute_AUROC(const deque< pair<float, bool> > &m_synergy_data)
{
	if( m_synergy_data.empty() ){
		throw __FILE__ ":compute_AUROC: no data!";
	}
	
	// Make a copy of the input data that we can sort
	deque< pair<float, bool> > local_synergy_data(m_synergy_data);
	
	// Rank our predictions in order of *ascending* predicted synergy
	// (i.e. from most to least predicted synergy).
	sort( local_synergy_data.begin(), local_synergy_data.end() );

	double auroc = 0.0;
	
	// The ROC curve is the true positive rate as a function of the false positive rate
	// 	true positive rate = sensitivity = recall = TP/(TP + FN)
	//	false positive rate = FP/(FP + TN)
	
	size_t TP = 0;
	size_t TN = 0;
	size_t FP = 0;
	size_t FN = 0;
	
	// Initialize the counts by first predicting all points to be 
	// "not synergistic" (i.e. negatives)
	for(deque< pair<float, bool> >::const_iterator i = local_synergy_data.begin();
		i != local_synergy_data.end();++i){
		
		if(i->second){
			// Since this point is actually synergistic, this is a false negative
			++FN;
		}
		else{
			// Since this point is actually *not* synergistic, this is a true negative
			++TN; 
		}
	}
	
	// Start integrating the curve at zero
	double last_tpr = 0.0;
	double last_fpr = 0.0;
	
	double last_rank = local_synergy_data.front().first;
	
	for(deque< pair<float, bool> >::const_iterator i = local_synergy_data.begin();
		i != local_synergy_data.end();++i){

		if(i->second){
		
			// When we encounter a synergistic point, it gets "moved" from a
			// false negative to a true positive.
			++TP;
			--FN;
		}
		else{
			
			// When we encounter a non synergistic point, it gets "moved" from a
			// true negative to a false positive.
			++FP;
			--TN;
		}
		
		// Only update the integral of the ROC when the ranking value changes
		// or we reach the last point
		if( (last_rank != i->first) || ( (i + 1) == local_synergy_data.end() ) ){
			
			// Make sure we don't try to compute 0/0 when TP + FN == 0 or FP + TN == 0
			const double tpr = (TP == 0) ? 0.0 : double(TP)/(TP + FN);
			const double fpr = (FP == 0) ? 0.0 : double(FP)/(FP + TN);
			
			const double delta_fpr = fpr - last_fpr;
			
			// Use the trapezoid rule to compute the integral
			auroc += 0.5*delta_fpr*(tpr + last_tpr);
			
			last_rank = i->first;
			last_tpr = tpr;
			last_fpr = fpr;
		}
	}
	
	return auroc;
}

double compute_gini(const double &m_auroc)
{
	return 2.0*m_auroc - 1.0;
}

// Report our rank, the signal we caught and the time spent running
void terminate_program(int m_sig)
{
	cerr << "[" << mpi_rank << "] caught signal " << m_sig << endl;
	cerr << report_run_time() << endl;
	
	MPI_Abort(MPI_COMM_WORLD, 0);
}

// Run time computes the total run time. The results are formatted as a string.
string report_run_time()
{
	double elapsed_time = MPI_Wtime() - start_time; // In sec
	
	const double elapsed_sec = fmod(elapsed_time, 60.0);
	
	elapsed_time = (elapsed_time - elapsed_sec)/60.0; // In min
	
	const double elapsed_min = fmod(elapsed_time, 60.0);
	elapsed_time = (elapsed_time - elapsed_min)/60.0; // In hour
	
	const double elapsed_hour = fmod(elapsed_time, 24.0);
	elapsed_time = (elapsed_time - elapsed_hour)/24.0; // In day
	
	stringstream sout;
	
	sout << "Run time is " 
		<< elapsed_time 
		<< " days, "
		<< elapsed_hour
		<< " hours, "
		<< elapsed_min
		<< " min and "
		<< elapsed_sec
		<< " sec";
	
	return sout.str();
}

void normalize(MAP<string, vector<float> > &m_features)
{
	for(MAP<string, vector<float> >::iterator i = m_features.begin();i != m_features.end();++i){
					
		const size_t num_features = i->second.size();
		double ave = 0.0;
		double stdev = 0.0;

		for(vector<float>::const_iterator j = i->second.begin();j != i->second.end();++j){
			ave += *j;
		}

		if(num_features > 0){
			ave /= num_features;
		}

		for(vector<float>::const_iterator j = i->second.begin();j != i->second.end();++j){

			const float local = ave - *j;

			stdev += local*local;
		}

		if(num_features > 1){
			stdev /= num_features - 1;
		}

		stdev = sqrt(stdev - ave*ave);

		if(stdev > 0.0){

			for(vector<float>::iterator j = i->second.begin();j != i->second.end();++j){
				*j = (*j - ave)/stdev;
			}
		}
	}
}


// Compute a two-tailed p-value
double compute_p_value(const float &m_x, const deque<float> &m_data)
{
	size_t num_greater = 0;
	
	for(deque<float>::const_iterator i = m_data.begin();i != m_data.end();++i){
		
		// Count the number of permuted enrichment scores that are *greater* than our target score
		num_greater += ( fabs(*i) > fabs(m_x) );
	}
	
	if( m_data.empty() ){
		throw __FILE__ ":compute_p_value: |data| == 0, can't compute p-value!";
	}
	
	// Return the fraction of permuted values that have an enirchment score greater than the target value
	return double(num_greater)/m_data.size();
}

vector<bool> find_invariant_features(const MAP<string, vector<float> > &m_drug_features, 
	const MAP<string, vector<float> > &m_blind_test_drug_features)
{
	if( m_drug_features.empty() ){
		return vector<bool>();
	}
	
	const size_t num_features = m_drug_features.begin()->second.size();
	
	// Assume that all features are constant 
	vector<bool> invariant_feature_mask(num_features, true);
	
	for(size_t i = 0;i < num_features;++i){
		
		float feature_value = 0.0;
		
		for(MAP<string, vector<float> >::const_iterator j = m_drug_features.begin();
			j != m_drug_features.end();++j){
			
			if(j->second.size() != num_features){
				throw __FILE__ ":find_invariant_features: Unexpected number of features (1)";
			}
			
			if( j == m_drug_features.begin() ){
				feature_value = j->second[i];
			}
			else{
				if(feature_value != j->second[i]){
				
					// This feature varies!
					invariant_feature_mask[i] = false;
					
					break; // Stop testing once we find the first example of variability
				}
			}
		}
		
		// If the feature appears to be invariant, we will need to test the blind features (if present)
		if(invariant_feature_mask[i] && !m_blind_test_drug_features.empty() ){
			
			for(MAP<string, vector<float> >::const_iterator j = m_blind_test_drug_features.begin();
				j != m_blind_test_drug_features.end();++j){
			
				if(j->second.size() != num_features){
					throw __FILE__ ":find_invariant_features: Unexpected number of features (2)";
				}


				if(feature_value != j->second[i]){
				
					// This feature varies!
					invariant_feature_mask[i] = false;

					break; // Stop testing once we find the first example of variability
				}
			}
		}
	}
	
	return invariant_feature_mask;
}


size_t count(const vector<bool> &m_x)
{
	size_t ret = 0;
	
	for(vector<bool>::const_iterator i = m_x.begin();i != m_x.end();++i){
		
		ret += (*i ? 1 : 0);
	}
	
	return ret;
}

template<class T>
void remove_invariant_features(vector<T> &m_x, const vector<bool> &m_invariant_mask)
{
	deque<T> ret;
	
	const size_t len = m_invariant_mask.size();
	
	if( len != m_x.size() ){
		throw __FILE__ ":remove_invariant_features: size mismatch!";
	}

	for(size_t i = 0;i < len;++i){
		
		// if m_invariant_mask is false, then we need to *keep* this variable!
		if(m_invariant_mask[i] == false){
			ret.push_back(m_x[i]);
		}
	}
	
	// Copy the list of variable features back to the calling function
	m_x.assign( ret.begin(),ret.end() );
}

// Subtract set m_b from set m_a and return the resulting set
vector<string> set_difference(const vector<string> &m_a, const vector<string> &m_b)
{
	vector<string> ret;
	
	for(vector<string>::const_iterator i = m_a.begin();i != m_a.end();++i){
		
		// Use a simple, brute-force search since we don't expect a large number of set elements
		// (i.e. drugs or cell lines).
		if( find(m_b.begin(), m_b.end(), *i) == m_b.end() ){
			
			// Only save the elements of m_a that are not found in m_b
			ret.push_back(*i);
		}
	}
	
	return ret;
}

void remove_drug(MAP< string, MAP<string, float> > &m_matrix, const string &m_drug)
{
	m_matrix.erase(m_drug);
	
	for(MAP< string, MAP<string, float> >::iterator i = m_matrix.begin();i != m_matrix.end();++i){
		i->second.erase(m_drug);
	}
}

float fraction_synergy(const MAP< string, MAP< string, MAP<string, float> > > &m_synergy_matrix,
	const float &m_synergy_threshold)
{

	size_t num_synergy = 0;
	size_t total = 0;
	
	for(MAP<string, MAP<string, MAP<string, float> > >::const_iterator i = m_synergy_matrix.begin();
		i != m_synergy_matrix.end();++i){

		for(MAP<string, MAP<string, float> >::const_iterator j = i->second.begin();
			j != i->second.end();++j){
				
			for(MAP<string, float>::const_iterator k = j->second.begin();
				k != j->second.end();++k){

				if(k->second == MISSING_DATA){
					continue;
				}

				num_synergy += (k->second <= m_synergy_threshold);
				++total;
			}
		}
	}
	
	return float(num_synergy)/total;
}

string probability_heat_map(const MAP< string, MAP< string, MAP<string, float> > > &m_synergy_matrix,
	const float &m_synergy_threshold)
{

	const vector<string> cells = keys(m_synergy_matrix);
	deque<string> drugs;
	
	for(MAP< string, MAP< string, MAP<string, float> > >::const_iterator i = m_synergy_matrix.begin();
		i != m_synergy_matrix.end();++i){
		
		for(MAP< string, MAP<string, float> >::const_iterator j = i->second.begin();
			j != i->second.end();++j){
			
			drugs.push_back(j->first);
		}
	}
	
	make_set(drugs);
	
	stringstream sout;
	
	// Print the header of drug names
	for(deque<string>::const_iterator i = drugs.begin();i != drugs.end();++i){
		sout << '\t' << *i;
	}
	
	sout << '\n';
	
	for(vector<string>::const_iterator c = cells.begin();c != cells.end();++c){
		
		// The row names are cells
		sout << *c;
		
		MAP< string, MAP< string, MAP<string, float> > >::const_iterator cell_iter = 
			m_synergy_matrix.find(*c);
			
		if( cell_iter == m_synergy_matrix.end() ){
			throw __FILE__ ":probability_heat_map:Unable to lookup cell name";
		}
		
		for(deque<string>::const_iterator d = drugs.begin();d != drugs.end();++d){
			
			MAP< string, MAP<string, float> >::const_iterator drug_iter = 
				cell_iter->second.find(*d);
				
			if( drug_iter == cell_iter->second.end() ){
				throw __FILE__ ":probability_heat_map: Unable to lookup drug name";
			}
			
			size_t num_valid = 0;
			size_t num_synergistic = 0;
			
			for(MAP<string, float>::const_iterator i = drug_iter->second.begin();
				i != drug_iter->second.end();++i){
				
				if(i->first != *d){ // Don't compare a drug against itself!
					
					if(i->second != MISSING_DATA){
						
						++num_valid;
						
						num_synergistic += (i->second <= m_synergy_threshold);
					}
				}
			}
			
			sout << '\t';
			
			if(num_valid > 0){
				sout << float(num_synergistic)/num_valid;
			}
		}
		
		sout << '\n';
	}
	
	return sout.str();
}
