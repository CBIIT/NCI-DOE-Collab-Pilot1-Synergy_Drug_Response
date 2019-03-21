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

#include <iostream>

#include <getopt.h>

using namespace std;

pair<string, string> parse_pair(const string &m_str);
string trim_flanking_space(const string &m_str);

void Options::parse(int argc, char *argv[])
{

	// Command line options:
	// [-o <output file>] (default is stdout)
	// [--o. <output file>] (default is stdout)
	// --drug <file of drug features>
	// [--drug.random (randomize the drug features)]
	// [--drug.blind (blind test set of drug features)]
	// [--drug.target (per-drug target/mechanism of action)]
	// --cell <file of cell features>
	// [--cell.random (randomize the cell features)]
	// [--cell.blind (blind test set of cell features)]
	// [--cell.normalize (normalize the cell features to zero mean and unit variance)]
	// Optional one hot encoding
	// [--hot.cell (encode cell lines as one-hot features)]
	// [--hot.tissue (encode cell line tissue type as one-hot features)]
	// [--hot.drug (encode drugs as one-hot features)]
	// Synergy variable to regress
	// --synergy | --synergy.train <"cell name->file of training synergy measurements"> (can be repeated)
	// [--synergy.test <"cell name->file of testing synergy measurements"> (can be repeated)]
	// [--synergy.NSC_to_CID (convert drug IDs from NSC to CID)]
	// [--synergy.test.random (randomize the *dependent* testing variables to compute a p-value)]
	// [--synergy.train.random (randomize the *dependent* training variables)]
	// [-s <random number seed> (default is time based)]
	// [--enable-pair-prediction (perform both single drug and drug *pair* prediction)]
	// Cross validation
	// [--cv.folds <number of cross validation folds>]
	// Generalization (when test and train are different datasets)
	// [--overlap.to_train]
	// [--overlap.to_test]
	// Synergy variable to regress
	// --synergy.average | --synergy.median | --synergy.fmean | --synergy.threshold <threshold value>
	// Random forest regression algorithm parameters
	// [--forest.size <number of trees>]
	// [--forest.bag <fraction of variables to bag>]
	// [--forest.leaf <minimum leaf size>]

	const char* options = "o:s:?h";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"drug", true, &config_opt, 1},
		{"drug.random", false, &config_opt, 2},
		{"synergy", true, &config_opt, 3}, // synonym for --synergy.train
		{"synergy.train", true, &config_opt, 3},
		{"synergy.NSC_to_CID", false, &config_opt, 4},
		{"synergy.test.random", false, &config_opt, 5},
		{"forest.size", true, &config_opt, 6},
		{"forest.bag", true, &config_opt, 7},
		{"forest.leaf", true, &config_opt, 8},
		{"shared-to-none", false, &config_opt, 21},
		{"synergy.average", false, &config_opt, 22},
		{"synergy.median", false, &config_opt, 23},
		{"synergy.fmean", false, &config_opt, 24},
		{"synergy.threshold", true, &config_opt, 25},
		{"drug.blind", true, &config_opt, 27},
		{"cell", true, &config_opt, 28},
		{"cell.random", false, &config_opt, 29},
		{"cell.blind", true, &config_opt, 30},
		{"cv.folds", true, &config_opt, 31},
		{"drug.target", true, &config_opt, 32},
		{"synergy.test", true, &config_opt, 33},
		{"enable-pair-prediction", false, &config_opt, 34},
		{"hot.cell", false, &config_opt, 35},
		{"hot.tissue", false, &config_opt, 36},
		{"hot.drug", false, &config_opt, 37},
		{"cell.normalize", false, &config_opt, 38},
		{"synergy.train.random", false, &config_opt, 39},
		{"overlap.to_train", false, &config_opt, 40},
		{"overlap.to_test", false, &config_opt, 41},
		{0,0,0,0} // Terminate options list
	};

	int opt_code;
	opterr = 0;

	// Set the defaults
	print_usage = (argc == 1);
	num_tree = DEFAULT_NUM_TREE;
	leaf_size = DEFAULT_LEAF_SIZE;
	bag_fraction = DEFAULT_BAG_FRACTION;

	num_fold = DEFAULT_NUM_FOLD;

	randomize_drug = false;
	randomize_cell = false;
	randomize_test_synergy = false;
	randomize_train_synergy = false;
	NSC_to_CID_y = false;
	score_func = UNDEFINED_SCORE;
	synergy_threshold = 0.0;

	seed = 0; // 0 -> use a time based seed

	// By default, we only perform single drug prediction. Turning on this
	// flag enables *both* single drug and drug pair-based prediction.
	use_pair_prediction = false;
	
	// By default, don't use one hot encoding
	one_hot_cell = false;
	one_hot_drug = false;
	one_hot_tissue = false;
	
	// By default, overlap is allowed
	overlap_to_train = false;
	overlap_to_test = false;
	
	while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

		switch( opt_code ){
			case 0:

				if(config_opt == 1){ // drug

					file_name_drug = optarg;
					break;
				}

				if(config_opt == 2){ // drug.random

					randomize_drug = true;
					break;
				}

				if(config_opt == 3){ // synergy or synergy.train

					const pair<string, string> local = parse_pair(optarg);

					if(train_cell_to_synergy_file.insert(local).second == false){

						cerr << "Duplicate training cell line key: " << optarg << endl;
						throw __FILE__ ":main: Unable to insert training key/value";
					}

					break;
				}

				if(config_opt == 4){ // synergy.NSC_to_CID

					NSC_to_CID_y = true;
					break;
				}

				if(config_opt == 5){ // synergy.test.random

					randomize_test_synergy = true;
					break;
				}

				if(config_opt == 6){ // forest.size

					const int temp = atoi(optarg);

					if(temp <= 0){

						cerr << "Please specify forest.size > 0" << endl;
						print_usage = true;
						return;
					}

					num_tree = size_t(temp);

					break;
				}

				if(config_opt == 7){ // forest.bag

					bag_fraction = atof(optarg);

					if( (bag_fraction <= 0.0) || (bag_fraction >= 1.0) ){

						cerr << "Please specify 0.0 < forest.bag < 1.0" << endl;
						
						print_usage = true;
						return;
					}

					break;
				}

				if(config_opt == 8){ // forest.leaf

					const int temp = atoi(optarg);

					if(temp < 1){

						cerr << "Please specify forest.leaf >= 1" << endl;
						
						print_usage = true;
						return;
					}

					leaf_size = size_t(temp);

					break;
				}					

				if(config_opt == 22){ // synergy.average

					score_func = AVERAGE_SCORE;
					break;
				}

				if(config_opt == 23){ // synergy.median

					score_func = MEDIAN_SCORE;
					break;
				}

				if(config_opt == 24){ // synergy.fmean

					score_func = FILTERED_MEAN_SCORE;
					break;
				}

				if(config_opt == 25){ // synergy.threshold

					score_func = COUNTING_SCORE;
					synergy_threshold = atof(optarg);
					break;
				}

				if(config_opt == 27){ // drug.blind

					file_name_blind_drug = optarg;
					break;
				}

				if(config_opt == 28){ // cell

					file_name_cell = optarg;
					break;
				}

				if(config_opt == 29){ // cell.random

					randomize_cell = true;
					break;
				}

				if(config_opt == 30){ // cell.blind

					file_name_blind_cell = optarg;
					break;
				}

				if(config_opt == 31){ // cv.folds

					const int temp = atoi(optarg);

					if(temp < 0){

						cerr << "Please specify number of cross-validation folds (--cv.folds) >= 0" << endl;
						
						print_usage = true;
						
						return;
					}

					num_fold = size_t(temp);

					break;
				}

				if(config_opt == 32){ // drug.target

					file_name_drug_target = optarg;
					break;
				}

				if(config_opt == 33){ // synergy.test

					const pair<string, string> local = parse_pair(optarg);

					if(test_cell_to_synergy_file.insert(local).second == false){

						cerr << "Duplicate test cell line key: " << optarg << endl;
						
						print_usage = true;
						
						return;
					}

					break;
				}

				if(config_opt == 34){ // enable-pair-prediction

					use_pair_prediction = true;
					break;
				}
				
				if(config_opt == 35){ // hot.cell

					one_hot_cell = true;
					break;
				}
				
				if(config_opt == 36){ // hot.tissue

					one_hot_tissue = true;
					break;
				}
				
				if(config_opt == 37){ // hot.drug

					one_hot_drug = true;
					break;
				}
				
				if(config_opt == 38){ // cell.normalize

					normalize_cell = true;
					break;
				}
				
				if(config_opt == 39){ // synergy.train.random

					randomize_train_synergy = true;
					break;
				}
				
				if(config_opt == 40){ // overlap.to_train

					overlap_to_train = true;
					break;
				}
				
				if(config_opt == 41){ // overlap.to_test

					overlap_to_test = true;
					break;
				}
				
				cerr << "Unknown flag!" << endl;
				break;
			case 'o':
				file_name_output = optarg;
				break;
			case 's':
				seed = abs( atoi(optarg) );
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

	if(print_usage){

		cerr << "Usage: predict_synergy version " << RX_VERSION << ":" << endl;
		cerr << "\t\t[-o <output file>] (default is stdout)" << endl;
		cerr << "\t\t[--drug <file of drug features>]" << endl;
		cerr << "\t\t[--drug.random (randomize the drug features)]" << endl;
		cerr << "\t\t[--drug.blind (specific test set of drug features)]" << endl;
		cerr << "\t\t[--drug.target (per-drug target/mechanism of action)]" << endl;
		cerr << "\t\t[--cell <file of cell features>]" << endl;
		cerr << "\t\t[--cell.random (randomize the cell features)]" << endl;
		cerr << "\t\t[--cell.blind (specific test set of cell features)]" << endl;
		cerr << "\t\t[--cell.normalize (normalize the cell features to zero mean and unit variance)]" << endl;
		cerr << "\t\t[--enable-pair-prediction (perform both single drug and drug *pair* prediction)]" << endl;
		cerr << "\tOptional one hot encoding" << endl;
		cerr << "\t\t[--hot.cell (encode cell lines as one-hot features)]" << endl;
		cerr << "\t\t[--hot.tissue (encode cell line tissue type as one-hot features)]" << endl;
		cerr << "\t\t[--hot.drug (encode drugs as one-hot features)]" << endl;
		cerr << "\tSynergy variable to regress" << endl;
		cerr << "\t\t--synergy.train | --synergy <\"cell name->file of training synergy measurements\"> (can be repeated)" << endl;
		cerr << "\t\t[--synergy.test <\"cell name->file of testing synergy measurements\"> (can be repeated)]" << endl;
		cerr << "\t\t--synergy.average | --synergy.median | --synergy.fmean | --synergy.threshold <threshold value>" << endl;
		cerr << "\t\t[--synergy.NSC_to_CID (convert drug IDs from NSC to CID)]" << endl;
		cerr << "\t\t[--synergy.test.random (randomize the testing synergy variables to compute permutation p-value)]" << endl;
		cerr << "\t\t[--synergy.train.random (randomize the training synergy variables)]" << endl;
		cerr << "\t\t[-s <random number seed> (default is time based)]" << endl;
		cerr << "\tCross-validation" << endl;
		cerr << "\t\t[--cv.folds <number of cross validation folds>] (default is " << DEFAULT_NUM_FOLD << ")" << endl;
		cerr << "\tGeneralization (when test and training are from different data sets)" << endl;
		cerr << "\t\t[--overlap.to_train] (shared cell lines/drugs assigend to train)" << endl;
		cerr << "\t\t[--overlap.to_test] (shared cell lines/drugs assigend to test)" << endl;
		cerr << "\tRegression parameters" << endl;
		cerr << "\t\t[--forest.size <number of trees>] (default is " << DEFAULT_NUM_TREE << ")" << endl;
		cerr << "\t\t[--forest.bag <fraction of variables to bag>] (default is " << DEFAULT_BAG_FRACTION << ")" << endl;
		cerr << "\t\t[--forest.leaf <minimum leaf size>] (default is " << DEFAULT_LEAF_SIZE << ")" << endl;
		
		return;
	}

	if(score_func == UNDEFINED_SCORE){

		cerr << "Please define the synergy score to regress on" << endl;
		
		print_usage = true;
		return;
	}

	if(seed == 0){
		seed = time(NULL);
	}
	
	// The user is not allowed to specifcy both cell-based and tissue-based one hot encoding
	if(one_hot_cell && one_hot_tissue){
		
		cerr << "One hot encoding of both cell (--hot.cell) and tissue (--hot.tissue) is not allowed!" << endl;
		print_usage = true;
		
		return;
	}
	
	// The user is not allowed to specify both one-hot and explicit cell line features
	if( (one_hot_cell || one_hot_tissue) && !file_name_cell.empty() ){
		
		cerr << "Please select either one-hot or explcit cell line features, but not both" << endl;
		print_usage = true;
		
		return;
	}

	// The user is not allowed to specify both one-hot and explicit drug features
	if( one_hot_drug && !file_name_drug.empty() ){
		
		cerr << "Please select either one-hot or explcit drug features, but not both" << endl;
		print_usage = true;
		
		return;
	}
	
	if( !file_name_blind_drug.empty() || !file_name_blind_cell.empty() ){
		
		if( !file_name_drug.empty() && file_name_blind_drug.empty() ){
			
			cerr << "Please provide --drug.blind features" << endl;
			
			print_usage = true;
			return;
		}
		
		if( !file_name_cell.empty() && file_name_blind_cell.empty() ){
			
			cerr << "Please provide --cell.blind features" << endl;
			
			print_usage = true;
			return;
		}
		
		if( file_name_output.empty() ){
			cerr << "\t**Warning** Blind test set output requires an output file (not stdout)!" << endl;
		}

		num_fold = 0;
	}
	
	if( !test_cell_to_synergy_file.empty() ){
	
		// When the user has specified an explicit test set, we must set
		// the number of folds to 1.
		num_fold = 1;
	}
	
	if(overlap_to_train && overlap_to_test){
		
		cerr << "Both --overlap.to_train and --overlap.to_test cannot be specified together" << endl;
			
		print_usage = true;
		return;
	}
}

pair<string, string> parse_pair(const string &m_str)
{
	const char delim = '=';
	
	const string::size_type loc = m_str.find(delim);
	
	if(loc == string::npos){
		throw __FILE__ ":parse_pair: Unable to find key/value delimeter";
	}
	
	return make_pair( 
		trim_flanking_space( m_str.substr(0, loc) ), 
		trim_flanking_space( m_str.substr( loc + 1, m_str.size() ) ) );
}


string trim_flanking_space(const string &m_str)
{
	if( m_str.empty() ){
		return m_str;
	}
	
	int start = 0;
	int stop = m_str.size() - 1;
	
	while(start <= stop){
		
		if( isspace(m_str[start]) ){
			++start;
		}
		else{
			break;
		}
	}
	
	while(start <= stop){
		
		if( isspace(m_str[stop]) ){
			--stop;
		}
		else{
			break;
		}
	}
	
	if(start > stop){
		return string();
	}
	
	return m_str.substr(start, stop - start + 1);
}
