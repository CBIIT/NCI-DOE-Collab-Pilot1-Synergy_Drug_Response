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

// Extract dose repsonse data for both drug pairs and individual drugs.
// J. D. Gans
// Bioscience Division, B-11
// Los Alamos National Laboratory
// Wed Sep  6 12:06:23 2017

#include <iostream>
#include <fstream>
#include <tuple>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#include <parallel/algorithm> // For __gnu_parallel::sort

#include "gemini.h"
#include "dose_response_function.h"
#include "stats.hpp"

using namespace std;

struct Options{

	bool print_usage;
	int lab_selection;
	
	string filename_input;
	string filename_bliss_ave_synergy;
	string filename_bliss_stdev_synergy;
	string filename_loewe_ave_synergy;
	string filename_loewe_stdev_synergy;
	string filename_pair_ave_growth;
	string filename_pair_stdev_growth;
	ComboScoreType dependent_variable_operator;
	
	Options(int argc, char* argv[]);
};

float loewe_pair_growth(const Concentration &m_d1, const float &m_alpha1, const float &m_beta1,
	const Concentration &m_d2, const float &m_alpha2, const float &m_beta2);
double loewe_func(const double &m_f,
		const Concentration &m_d1, const float &m_alpha1, const float &m_beta1,
		const Concentration &m_d2, const float &m_alpha2, const float &m_beta2);
void write_header(ofstream &m_fout, const deque<string> &m_cell_names);
deque<GrowthRecord>::const_iterator find_by_conc(const deque<GrowthRecord>::const_iterator &m_begin,
	const deque<GrowthRecord>::const_iterator &m_end, const Concentration &m_conc);
deque< pair<Concentration, float> >log_conc_growth(const deque<GrowthRecord> &m_data);

deque< tuple<Concentration, Concentration, float> > compute_bliss_synergy(const deque<GrowthRecord> &m_pair_data, 
	const deque<GrowthRecord> &m_drug_1_data, const deque<GrowthRecord> &m_drug_2_data);
deque< tuple<Concentration, Concentration, float> > compute_loewe_synergy(const deque<GrowthRecord> &m_pair_data, 
	const deque<GrowthRecord> &m_drug_1_data, const deque<GrowthRecord> &m_drug_2_data);
deque< tuple<Concentration, Concentration, float> > compute_pair_growth(const deque<GrowthRecord> &m_pair_data);
void accumulate_data(MAP<Concentration, MAP<Concentration, deque<float> > > &m_data, 
	const deque< tuple<Concentration, Concentration, float> > &m_sample);
pair<float, float> compute_score(const MAP<Concentration, MAP<Concentration, deque<float> > > &m_data, 
        const ComboScoreType &m_op);
void write_results(ofstream &m_fave, ofstream &m_fstdev, 
	const MAP< string, MAP<string, MAP<Concentration, MAP<Concentration, deque<float> > > > > &m_data, 
	const deque<string> &m_cell_names, const ComboScoreType &m_op);

int main(int argc, char *argv[])
{
	try{
		cerr << "gemini_prep version " << VERSION << endl;
		
		// Command line options:
		// -i <input file>
		// [-p <output file>]
		// [--pair.min (minimum dependent variable over all concentrations; default)]
		// [--pair.median (median dependent variable over all concentrations)]
		// [--pair.average (average dependent variable over all concentrations)]
		// [--exclude.1A] (exclude data from site 1A)
		// [--exclude.FF] (exclude data from site FF)
		// [--exclude.FG] (exclude data from site FG)
		Options opt(argc, argv);
		
		if(opt.print_usage){
			
			cerr << "Usage:" << endl;
			cerr << "\t-i <input file>" << endl;
			cerr << "\t[-p <output file prefix>] (default is no prefix)" << endl;
			cerr << "\t[--pair.min] (default; minimum dependent variable over all concentrations)" << endl;
			cerr << "\t[--pair.median] (median dependent variable over all concentrations)" << endl;
			cerr << "\t[--pair.average] (average dependent variable over all concentrations)" << endl;
			cerr << "\t[--exclude.1A] (exclude data from site 1A)" << endl;
			cerr << "\t[--exclude.FF] (exclude data from site FF)" << endl;
			cerr << "\t[--exclude.FG] (exclude data from site FG)" << endl;
			
			return EXIT_SUCCESS;
		}
		
		if( opt.filename_input.empty() ){
		
			cerr << "Please specify an input filename (-i)" << endl;
			return EXIT_FAILURE;
		}
		
		if( !( opt.lab_selection & (S_FF | S_FG | S_1A) ) ){
		
			cerr << "All site data has been excluded!" << endl;
			return EXIT_FAILURE;
		}

		cerr << "Dependent variable (synergy) pair data operation: ";
		
		switch(opt.dependent_variable_operator){
			case COMBO_AVERAGE:
				cerr << "Average over all concentrations tested" << endl;
				break;
			case COMBO_MEDIAN:
				cerr << "Median over all concentrations tested" << endl;
				break;
			case COMBO_MIN:
				cerr << "Minimum over all concentrations tested" << endl;
				break;
			default:
				throw __FILE__ ":main: Unknown dependent_variable_operator";
		};
		
		time_t profile = time(NULL);
		
		cerr << "Writing Bliss drug average synergy data to " << opt.filename_bliss_ave_synergy << endl;
		cerr << "Writing Bliss drug standard deviation synergy data to " << opt.filename_bliss_stdev_synergy << endl;
		
		cerr << "Writing Loewe drug average synergy data to " << opt.filename_loewe_ave_synergy << endl;
		cerr << "Writing Loewe drug standard deviation synergy data to " << opt.filename_loewe_stdev_synergy << endl;
		
		cerr << "Writing drug pair average growth data to " << opt.filename_pair_ave_growth << endl;
		cerr << "Writing drug pair standard deviation growth data to " << opt.filename_pair_stdev_growth << endl;
		
        	deque<GrowthRecord> data;
        
		if( !parse_combo(opt.filename_input, data, opt.lab_selection) ){
			
			cerr << "Unable to parse combo data file" << endl;
			return EXIT_FAILURE;
		}
		
		cerr << "Found a total of " << data.size() << " raw growth records" << endl;
		
		////////////////////////////////////////////////////////////////////
		// Extract the cell name
		deque<string> cell_names;

        	// Extract all of the cell names
        	for(deque<GrowthRecord>::const_iterator i = data.begin();i != data.end();++i){

        	    if( i->is_pair() ){
                	cell_names.push_back(i->cellname);
        	    }
        	}

        	make_set(cell_names);

        	const size_t num_cell = cell_names.size();

        	cerr << "Found " << num_cell << " cell line names" << endl;
		
		////////////////////////////////////////////////////////////////////
		// Prepare the headers for output files
		ofstream fbliss_ave_synergy( opt.filename_bliss_ave_synergy.c_str() );

		if(!fbliss_ave_synergy){
			throw __FILE__ ":main: Unable to open Bliss average synergy file for writing";
		}

		write_header(fbliss_ave_synergy, cell_names);
        
		ofstream fbliss_stdev_synergy( opt.filename_bliss_stdev_synergy.c_str() );

		if(!fbliss_stdev_synergy){
			throw __FILE__ ":main: Unable to open Bliss standard deviation synergy file for writing";
		}

		write_header(fbliss_stdev_synergy, cell_names);
		
		//*****************************
		
		ofstream floewe_ave_synergy( opt.filename_loewe_ave_synergy.c_str() );

		if(!floewe_ave_synergy){
			throw __FILE__ ":main: Unable to open Loewe average synergy file for writing";
		}

		write_header(floewe_ave_synergy, cell_names);
		
		ofstream floewe_stdev_synergy( opt.filename_loewe_stdev_synergy.c_str() );

		if(!floewe_stdev_synergy){
			throw __FILE__ ":main: Unable to open Loewe standard deviation synergy file for writing";
		}

		write_header(floewe_stdev_synergy, cell_names);
		
		//*****************************
		
		// Write the parsed data to disk
		ofstream fpair_ave_growth( opt.filename_pair_ave_growth.c_str() );

		if(!fpair_ave_growth){
			throw __FILE__ ":main: Unable to open paired average drug growth file for writing";
		}

		write_header(fpair_ave_growth, cell_names);
		
		// Write the parsed data to disk
		ofstream fpair_stdev_growth( opt.filename_pair_stdev_growth.c_str() );

		if(!fpair_stdev_growth){
			throw __FILE__ ":main: Unable to open paired standard deviation drug growth file for writing";
		}

		write_header(fpair_stdev_growth, cell_names);
		
		////////////////////////////////////////////////////////////////////
		// Group the growth records into "experiments", where each experiment is a collection of
		// drug pair growth measurements and associated single drug growth records
		
		// Step 1: Split the raw growth records into drug pair and single drug data
		deque<GrowthRecord> drug_pair_data;
		deque<GrowthRecord> single_drug_data;
		
		while( !data.empty() ){
			
			const GrowthRecord rec = data.back();
			data.pop_back();
			
			if( rec.is_pair() ){
				drug_pair_data.push_back(rec);
			}
			else{
				single_drug_data.push_back(rec);
			}
			
		}
		
		const size_t num_drug_pair = drug_pair_data.size();
		const size_t num_single_drug = single_drug_data.size();
		
		cerr << "\tFound a total of " << num_drug_pair << " drug pair growth records" << endl;
		cerr << "\tFound a total of " << num_single_drug << " single drug growth records" << endl;
		
		// Step 1 is to stort all of the Growth records by:
		// screener -> study -> cellname -> drug 1 -> drug 2
		__gnu_parallel::sort( drug_pair_data.begin(), drug_pair_data.end() );
		__gnu_parallel::sort( single_drug_data.begin(), single_drug_data.end() );
		
		size_t num_valid = 0;
		size_t num_invalid = 0;
		
		// Store the final data (including replicates) with the following multidimensional map:
		// 	[drug1,drug2][cellname][drug 1 concentration][drug 2 concentration] -> (sum, replicate count)
		MAP< string, 
			MAP<string, 
				MAP<Concentration, 
					MAP<Concentration, 
						deque<float>
					> 
				>
			>
		> bliss_synergy;
		
		MAP< string, 
			MAP<string, 
				MAP<Concentration, 
					MAP<Concentration, 
						deque<float>
					>
				>
			>
		> loewe_synergy;
		
		MAP< string, 
			MAP<string, 
				MAP<Concentration, 
					MAP<Concentration, 
						deque<float>
					>
				>
			>
		> pair_growth;
		
		while( !drug_pair_data.empty() ){
			
			// Extract the data for a single pair experiment
			deque<GrowthRecord> pair_data;
			deque<GrowthRecord> drug_1_data;
			deque<GrowthRecord> drug_2_data;
			
			// 1) Make sure the pair growth values are greater than zero
			// 2) Make sure that we have one, and only one, single drug record for each drug and
			// concentration found in the pair data
			bool is_valid = true;
			
			while( !drug_pair_data.empty() ){
				
				if( pair_data.empty() ){
				
					pair_data.push_back( drug_pair_data.back() );
					drug_pair_data.pop_back();
					continue;
				}
				
				const GrowthRecord &curr = pair_data.back();
				const GrowthRecord &next = drug_pair_data.back();
				
				if( (curr.screener != next.screener) ||
				    (curr.cellname != next.cellname) ||
				    (curr.nsc_1 != next.nsc_1) ||
				    (curr.nsc_2 != next.nsc_2) ||
				    (curr.study != next.study) ){
				    
				    	break;
				}
				
				// Impose an additional plate-based criteria on measurements 
				// from screening sites *other* than S_1A
				if( (curr.screener != S_1A) && 
				    (curr.plate != next.plate) ){
				
					break;
				}

				pair_data.push_back( drug_pair_data.back() );
				drug_pair_data.pop_back();
			}
			
			if( pair_data.empty() ){
				throw __FILE__ ":main: Failed to find any pair data ...";
			}
			
			// There are a handfull of measurements from screening site 1A that have negative 
			// growth values due to negative test_value's.
			for(deque<GrowthRecord>::const_iterator i = pair_data.begin();i != pair_data.end();++i){
				
				const float growth = i->test_value/i->control_value;
				
				if(growth < 0.0){
					
					is_valid = false;
					
					cerr << '\t' << i->screener << '\t' << i->study << '\t' 
						<< i->plate << '\t' << i->cellname << '\t' << i->nsc_1
						<< '\t' << i->nsc_2 << " <- negative pair growth value" << endl;
				}
			}
			
			if(is_valid == false){
        		    	continue;
			}

			// Please note that the following test code (for |pair_data| != 9 and
			// |pair_data| != 15) finds 45 experiments from screening site S_1A
			// that have 30 pairs each. These pairs have the same cellname and drug names
			// but different plate names. However, splitting by plate name does not help,
			// as these appear to actually be two different experiments (involving slightly
			// different concentration ranges) split over three plates!
			// All of these pair measurments involve:
			//	study = 1307AC42
			//	drug 1 = 3685
			//	drug 2 = 41867
			//
			// For now (since drug pairs are treated mostly independently), we will let
			// this anomally slide and include all of these pair measurements
			//if( (pair_data.size() != 9) && (pair_data.size() != 15) ){
			//	
			//	cout << "Found " << pair_data.size() << " pairs" << endl;
			//	
			//	for(deque<GrowthRecord>::const_iterator i = pair_data.begin();
			//		i != pair_data.end();++i){
			//		
			//		cout << '\t' << i->screener << '\t' << i->study << '\t' 
			//			<< i->plate << '\t' << i->cellname << '\t' << i->nsc_1
			//			<< '\t' << i->nsc_2 << endl;
			//	}
			//}
			
			const GrowthRecord &first_pair = pair_data.front();
			
			// Extract all of the single drug records for this experiment
			deque<GrowthRecord>::const_iterator iter = 
				lower_bound( single_drug_data.begin(), single_drug_data.end(), first_pair);
				
			if( iter == single_drug_data.end() ){
				
				cerr << '\t' << first_pair.screener << '\t' << first_pair.study << '\t' 
						<< first_pair.plate << '\t' << first_pair.cellname << '\t' << first_pair.nsc_1
						<< '\t' << first_pair.nsc_2 << endl;
				throw __FILE__ ":main: Unable to find pair experiment in single drug data!";
			}
			
			while( iter != single_drug_data.end() ){
				
				if( (first_pair.screener != iter->screener) ||
				    (first_pair.cellname != iter->cellname) ||
				    (first_pair.study != iter->study) ){
				    	break;
				}
				
				if(first_pair.nsc_1 == iter->nsc_1){
					drug_1_data.push_back(*iter);
				}
				else{
					if(first_pair.nsc_2 == iter->nsc_1){
						drug_2_data.push_back(*iter);
					}
				}
				
				++iter;
			}
		
			// Make sure that all of the single drug values refer to pair measurement
			const size_t num_drug_1 = drug_1_data.size();
			const size_t num_drug_2 = drug_2_data.size();
			
			vector<bool> drug_1_mask(num_drug_1, false);
			vector<bool> drug_2_mask(num_drug_2, false);
			
			for(deque<GrowthRecord>::const_iterator i = pair_data.begin();
				i != pair_data.end();++i){
				
				const Concentration conc_1 = i->concentration_1();
				const Concentration conc_2 = i->concentration_2();
				
				size_t num_match_1 = 0;
				size_t num_match_2 = 0;
				
				for(size_t j = 0;j < num_drug_1;++j){
					
					if(drug_1_data[j].concentration_1() == conc_1){
					
						++num_match_1;
						drug_1_mask[j] = true;
					}
				}
				
				if(num_match_1 != 1){
					
					cerr << '\t' << i->screener << '\t' << i->study << '\t' 
						<< i->plate << '\t' << i->cellname << '\t' << i->nsc_1
						<< '\t' << i->nsc_2 << " <- missing drug 1 @ " 
						<< conc_1.value << conc_1.units << endl;

					is_valid = false;
				}
				
				for(size_t j = 0;j < num_drug_2;++j){
					
					if(drug_2_data[j].concentration_1() == conc_2){
					
						++num_match_2;
						drug_2_mask[j] = true;
					}
				}
				
				if(num_match_2 != 1){
					
					cerr << '\t' << i->screener << '\t' << i->study << '\t' 
						<< i->plate << '\t' << i->cellname << '\t' << i->nsc_1
						<< '\t' << i->nsc_2 << " <- missing drug 2 @ "
						<< conc_2.value << conc_2.units << endl;

					is_valid = false;
				}
			}
			
			// Check that all of the single drugs were refered to at least once
			for(size_t i = 0;i < num_drug_1;++i){

				if(drug_1_mask[i] == false){

					cerr << '\t' << drug_1_data[i].screener << '\t' << drug_1_data[i].study << '\t' 
						<< drug_1_data[i].plate << '\t' << drug_1_data[i].cellname 
						<< '\t' << drug_1_data[i].nsc_1 << " <- unreferenced drug 1" << endl;
					is_valid = false;
				}
			}
			
			for(size_t i = 0;i < num_drug_2;++i){

				if(drug_2_mask[i] == false){

					cerr << '\t' << drug_2_data[i].screener << '\t' << drug_2_data[i].study << '\t' 
						<< drug_2_data[i].plate << '\t' << drug_2_data[i].cellname 
						<< '\t' << drug_2_data[i].nsc_1 << " <- unreferenced drug 2" << endl;
					is_valid = false;
				}
			}
			
			if(!is_valid){
				
				++num_invalid;
				continue;
			}
			
			++num_valid;
			
			// Whew! We are finally ready to compute growth and synergy values!
			const deque< tuple<Concentration, Concentration, float> > bliss = 
				compute_bliss_synergy(pair_data, drug_1_data, drug_2_data);
				
			const deque< tuple<Concentration, Concentration, float> > loewe = 
				compute_loewe_synergy(pair_data, drug_1_data, drug_2_data);
				
			const deque< tuple<Concentration, Concentration, float> > growth = 
				compute_pair_growth(pair_data);
			
			// Make sure that we have the same number of synergy and growth results
			const size_t num_pair = pair_data.size();
			
			if(bliss.size() != num_pair){
				throw __FILE__ ":main: |bliss| != |pairs|";
			}
			
			if(loewe.size() != num_pair){
				throw __FILE__ ":main: |loewe| != |pairs|";
			}
			
			if(growth.size() != num_pair){
				throw __FILE__ ":main: |growth| != |pairs|";
			}
			
			// Plot the bliss vs the loewe synergy values
			//for(size_t i = 0;i < num_pair;++i){
			//	cout << bliss[i] << '\t' << loewe[i] << endl;
			//}
			
			// DEBUG
			#ifdef DISCORDANT_SYNERGY
			bool debug = false;
			
			for(size_t i = 0;i < num_pair;++i){
				
				if( (bliss[i] == 0.0) && ( fabs(loewe[i] - bliss[i]) > 0.1) ){
					debug = true;
				}
			}
			
			if(debug){
				
				for(size_t i = 0;i < num_pair;++i){
					cout << bliss[i] << '\t' << loewe[i] << endl;
				}
			
				cout << '\t' << endl;
				
				for(deque<GrowthRecord>::const_iterator i = pair_data.begin();
					i != pair_data.end();++i){
					
					cout << i->screener << '\t' << i->study << '\t' 
						<< i->plate << '\t' << i->cellname << '\t' 
						<< i->nsc_1 << '\t' << i->nsc_2 << '\t'
						<< i->conc_1 << i->conc_unit_1 << '\t' 
						<< i->conc_2 << i->conc_unit_2 << '\t'
						<< min(1.0f, i->test_value/i->control_value) << endl;
				}
				
				cout << '\t' << endl;
				
				for(size_t i = 0;i < num_drug_1;++i){

					cout << drug_1_data[i].screener << '\t' << drug_1_data[i].study << '\t' 
						<< drug_1_data[i].plate << '\t' << drug_1_data[i].cellname 
						<< '\t' << drug_1_data[i].nsc_1 << '\t' 
						<< drug_1_data[i].conc_1 << drug_1_data[i].conc_unit_1 << '\t'
						<< min(1.0f, drug_1_data[i].test_value/drug_1_data[i].control_value) << endl;
				}
				
				cout << '\t' << endl;
				
				DoseResponseFunction sigmoid1(0.0 /*no regularization*/);
				sigmoid1.fit( log_conc_growth(drug_1_data) );
				
				cout << "Drug 1 alpha = " << sigmoid1.get_alpha() << endl;
				cout << "Drug 1 beta = " << sigmoid1.get_beta() << endl;
				
				cout << '\t' << endl;
				
				for(size_t i = 0;i < num_drug_2;++i){

					cout << drug_2_data[i].screener << '\t' << drug_2_data[i].study << '\t' 
						<< drug_2_data[i].plate << '\t' << drug_2_data[i].cellname 
						<< '\t' << drug_2_data[i].nsc_1 << '\t'
						<< drug_2_data[i].conc_1 << drug_2_data[i].conc_unit_1 << '\t'
						<< min(1.0f, drug_2_data[i].test_value/drug_2_data[i].control_value) << endl;
				}
				
				cout << '\t' << endl;
				
				DoseResponseFunction sigmoid2(0.0 /*no regularization*/);
				sigmoid2.fit( log_conc_growth(drug_2_data) );
				
				cout << "Drug 2 alpha = " << sigmoid2.get_alpha() << endl;
				cout << "Drug 2 beta = " << sigmoid2.get_beta() << endl;
				
				cout << '\t' << endl;
				
				for(deque<GrowthRecord>::const_iterator i = pair_data.begin();
					i != pair_data.end();++i){
					
					const double eps = 1.0e-3;
					
					cout << i->screener << '\t' << i->study << '\t' 
						<< i->plate << '\t' << i->cellname << '\t' 
						<< i->nsc_1 << '\t' << i->nsc_2 << '\t'
						<< i->conc_1 << i->conc_unit_1 << '\t'
						<< i->conc_2 << i->conc_unit_2 << '\t'
						<< loewe_func(eps,
							i->concentration_1(), sigmoid1.get_alpha(), sigmoid1.get_beta(), 
							i->concentration_2(), sigmoid2.get_alpha(), sigmoid2.get_beta() ) << '\t'
						<< loewe_func(1.0 - eps,
							i->concentration_1(), sigmoid1.get_alpha(), sigmoid1.get_beta(), 
							i->concentration_2(), sigmoid2.get_alpha(), sigmoid2.get_beta() ) << endl;
				}
				
				cout << "------------------------------------------------------------" << endl;
			}
			
			#endif // DISCORDANT_SYNERGY
			
			const string drug_pair_key = first_pair.nsc_1 + "," + first_pair.nsc_2;
			
			accumulate_data(bliss_synergy[drug_pair_key][first_pair.cellname], bliss);
			accumulate_data(loewe_synergy[drug_pair_key][first_pair.cellname], loewe);
			accumulate_data(pair_growth[drug_pair_key][first_pair.cellname], growth);
		}
		
		// Write the processed data to disk
		write_results(fbliss_ave_synergy, fbliss_stdev_synergy, 
			bliss_synergy, cell_names, opt.dependent_variable_operator);
		write_results(floewe_ave_synergy, floewe_stdev_synergy, 
			loewe_synergy, cell_names, opt.dependent_variable_operator);
		write_results(fpair_ave_growth, fpair_stdev_growth, 
			pair_growth, cell_names, opt.dependent_variable_operator);
	
		cerr << "Found " << num_valid << " valid pairs (" 
			<< (100.0*num_valid)/(num_valid + num_invalid) << "%) and " 
			<< num_invalid << " invalid pairs ("
			<< (100.0*num_invalid)/(num_valid + num_invalid) << "%)" << endl;
		
		profile = time(NULL) - profile;
		
		cerr << "Processed data in " << profile << " sec" << endl;
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

deque< pair<Concentration, float> >log_conc_growth(const deque<GrowthRecord> &m_data)
{
	deque< pair<Concentration, float> > ret;
	
	// Convert the concentration to the log concentration
	for(deque<GrowthRecord>::const_iterator i = m_data.begin();i != m_data.end();++i){
	
		Concentration conc = i->concentration_1();
		
		if(conc.value <= 0.0){
			throw __FILE__ ":log_conc_growth: Unable to compute log of negative number";
		}
		
		conc.value = log(conc.value);
		
		const float growth = min(1.0f, i->test_value/i->control_value);
		
		ret.push_back( make_pair(conc, growth) );
	}
	
	sort( ret.begin(), ret.end() );
	
	return ret;
}

double loewe_func(const double &m_f,
		const Concentration &m_d1, const float &m_alpha1, const float &m_beta1,
		const Concentration &m_d2, const float &m_alpha2, const float &m_beta2)
{

	const double max_value = 1.0e6;
	
	const double local = log(1.0/m_f - 1.0);
	
	//return m_d1.value/exp(local/m_alpha1 - m_beta1) + m_d2.value/exp(local/m_alpha2 - m_beta2) - 1.0;
	
	// Rewrite the Loewe function to cancel terms in the exponential *before* calculating the exp() function.
	const double arg1 = log(m_d1.value) - local/m_alpha1 + m_beta1;
	const double arg2 = log(m_d2.value) - local/m_alpha2 + m_beta2;
	
	// DEBUG
	//cerr << "m_f = " << m_f << endl;
	//cerr << "m_alpha1 = " << m_alpha1 << endl;
	//cerr << "m_beta1 = " << m_beta1 << endl;
	//cerr << "m_alpha2 = " << m_alpha2 << endl;
	//cerr << "m_beta2 = " << m_beta2 << endl;
	
	//cerr << "arg1 = " << arg1 << endl;
	//cerr << "arg2 = " << arg2 << endl;
	
	//cerr << "exp(arg1) = " << exp(arg1) << endl;
	//cerr << "exp(arg2) = " << exp(arg2) << endl;
	//cerr << "exp(arg1) + exp(arg2) - 1.0 = " << exp(arg1) + exp(arg2) - 1.0 << endl;
	
	// Guard against overflow
	const double max_arg = 1.0e3;
	
	if( (arg1 >= max_arg) || (arg2 >= max_arg) ){
		return max_value;
	}
	
	return min(max_value, exp(arg1) + exp(arg2) - 1.0);
}
	
float loewe_pair_growth(const Concentration &m_d1, const float &m_alpha1, const float &m_beta1,
	const Concentration &m_d2, const float &m_alpha2, const float &m_beta2)
{
	// Use root finding to solve for the pair growth value that is expected under the Loewe synergy model
	const double eps = 1.0e-3;
	
	double lower_bound = eps;
	double lower_bound_value = loewe_func(lower_bound,
		m_d1, m_alpha1, m_beta1, m_d2, m_alpha2, m_beta2);
		
	double upper_bound = 1.0 - eps;
	double upper_bound_value = loewe_func(upper_bound,
		m_d1, m_alpha1, m_beta1, m_d2, m_alpha2, m_beta2);
	
	// DEBUG
	//cerr << "lower_bound_value = " << lower_bound_value << endl;
	//cerr << "upper_bound_value = " << upper_bound_value << endl;
	
	if( (   (lower_bound_value < 0.0) && (upper_bound_value < 0.0) ) ||
		( (lower_bound_value > 0.0) && (upper_bound_value > 0.0) ) ){
		
		//cerr << "lower_bound_value = " << lower_bound_value << endl;
		//cerr << "upper_bound_value = " << upper_bound_value << endl;
		//throw __FILE__ ":loewe_pair_growth: Both values are on the same side as zero";
		
		// From Kashif et al. (2017) Oncotarget. 2017 Nov 28; 8(61): 103952-103967
		// Return the argmin over f of the loewe function
		if( fabs(lower_bound_value) > fabs(upper_bound_value) ){
		
			// The upper_bound_value is closer to zero (the root), so return the
			// maximum growth. This is equivalent to giving up on root finding, and instead
			// returning argmin over f for loewe_func(f)
			return 1.0;
		}
		else{
			// The lower_bound_value is closer to zero (the root), so return the
			// maximum growth. This is equivalent to giving up on root finding, and instead
			// returning argmin over f for loewe_func(f)
			return 0.0;
		}
	}
	
	// DEBUG
	//cerr << "m_d1.value = " << m_d1.value << endl;
	//cerr << "m_alpha1 = " << m_alpha1 << endl;
	//cerr << "m_beta1 = " << m_beta1 << endl;
	
	//cerr << "m_d2.value = " << m_d2.value << endl;
	//cerr << "m_alpha2 = " << m_alpha2 << endl;
	//cerr << "m_beta2 = " << m_beta2 << endl;
	
	//size_t iteration = 0;
	
	const float convergence_threshold = 1.0e-4;
	
	while( fabs(upper_bound - lower_bound) > convergence_threshold){
		
		const double mid = 0.5*(upper_bound + lower_bound);
		
		const double mid_value = loewe_func(mid,
			m_d1, m_alpha1, m_beta1, m_d2, m_alpha2, m_beta2);
		
		// DEBUG
		//cerr << "\titer = " << ++iteration << endl;
		
		//cerr << "\tupper_bound = " << upper_bound << endl;
		//cerr << "\tupper_bound_value = " << upper_bound_value << endl;
		
		//cerr << "\tlower_bound = " << lower_bound << endl;
		//cerr << "\tlower_bound_value = " << lower_bound_value << endl;
		
		//cerr << "\tmid = " << mid << endl;
		//cerr << "\tmid_value = " << mid_value << endl;
		
		//cerr << endl;
		
		if(mid_value < 0.0){
			
			if(lower_bound_value <= 0.0){
				
				if( fabs(mid_value) <= fabs(lower_bound_value) ){
					
					lower_bound = mid;
					lower_bound_value = mid_value;
				}
				else{
					throw __FILE__ ":loewe_pair_growth: Lower bound (-) step did not converge";
				}
			}
			else{ // upper_bound_value < 0.0
				
				if( fabs(mid_value) <= fabs(upper_bound_value) ){
					
					upper_bound = mid;
					upper_bound_value = mid_value;
				}
				else{
					throw __FILE__ ":loewe_pair_growth: Upper bound (-) step did not converge";
				}
			}
		}
		else{ // mid_value >= 0.0
		
			if(lower_bound_value >= 0.0){
				
				if( fabs(mid_value) <= fabs(lower_bound_value) ){
					
					lower_bound = mid;
					lower_bound_value = mid_value;
				}
				else{
					throw __FILE__ ":loewe_pair_growth: Lower bound (+) step did not converge";
				}
			}
			else{ // upper_bound_value >= 0.0
				
				if( fabs(mid_value) <= fabs(upper_bound_value) ){
					
					upper_bound = mid;
					upper_bound_value = mid_value;
				}
				else{
					const double local = log(1.0/upper_bound - 1.0);
					
					cerr << "local = " << local << endl;
					
					cerr << "arg func_A = " << log(m_d1.value) - local/m_alpha1 + m_beta1 << endl;
					cerr << "arg func_B = " << log(m_d2.value) - local/m_alpha2 + m_beta2 << endl;
					
					cerr << "func_A = " << exp(log(m_d1.value) - local/m_alpha1 + m_beta1) << endl;
					cerr << "func_B = " << exp(log(m_d2.value) - local/m_alpha2 + m_beta2) << endl;
					
					cerr << "m_d1.value = " << m_d1.value << endl;
					cerr << "m_alpha1 = " << m_alpha1 << endl;
					cerr << "m_beta1 = " << m_beta1 << endl;

					cerr << "m_d2.value = " << m_d2.value << endl;
					cerr << "m_alpha2 = " << m_alpha2 << endl;
					cerr << "m_beta2 = " << m_beta2 << endl;

					cerr << "upper_bound = " << upper_bound << endl;
					cerr << "upper_bound_value = " << upper_bound_value << endl;
					
					cerr << "lower_bound = " << lower_bound << endl;
					cerr << "lower_bound_value = " << lower_bound_value << endl;
					
					cerr << "mid = " << mid << endl;
					cerr << "mid_value = " << mid_value << endl;
					
					throw __FILE__ ":loewe_pair_growth: Upper bound (+) step did not converge";
				}
			}
		}
		
	}
	
	return 0.5*(upper_bound + lower_bound);
}

Options::Options(int argc, char* argv[])
{

	const char* options = "i:p:?h";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"exclude.1A", false, &config_opt, 5},
		{"exclude.FF", false, &config_opt, 6},
		{"exclude.FG", false, &config_opt, 7},
		{"pair.min", false, &config_opt, 11},
		{"pair.median", false, &config_opt, 12},
		{"pair.average", false, &config_opt, 13},
		{0,0,0,0} // Terminate options list
	};

	int opt_code;
	opterr = 0;

	string output_prefix = "";
	
	print_usage = (argc == 1);
	filename_input = "";
	filename_bliss_ave_synergy = "bliss_ave_synergy.csv";
	filename_bliss_stdev_synergy = "bliss_stdev_synergy.csv";
	filename_loewe_ave_synergy = "loewe_ave_synergy.csv";
	filename_loewe_stdev_synergy = "loewe_stdev_synergy.csv";
	filename_pair_ave_growth = "pair_ave_growth.csv";
	filename_pair_stdev_growth = "pair_stdev_growth.csv";
	lab_selection = S_FF | S_1A | S_FG; // All laboratories by default
        dependent_variable_operator = COMBO_MIN;

	while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

		switch( opt_code ){
			case 0:
				if(config_opt == 5){ // exclude.1A

					// Turn off the bits associated with S_1A
					lab_selection &= ~S_1A;
					break;
				}

				if(config_opt == 6){ // exclude.FF

					// Turn off the bits associated with S_FF
					lab_selection &= ~S_FF;
					break;
				}

				if(config_opt == 7){ // exclude.FG

					// Turn off the bits associated with S_FG
					lab_selection &= ~S_FG;
					break;
				}

				if(config_opt == 11){ // pair.min

					dependent_variable_operator = COMBO_MIN;
					break;
				}

				if(config_opt == 12){ // pair.median

					dependent_variable_operator = COMBO_MEDIAN;
					break;
				}

				if(config_opt == 13){ // pair.average

					dependent_variable_operator = COMBO_AVERAGE;
					break;
				}

				cerr << "Unknown flag!" << endl;
				break;
			case 'i':
				filename_input = optarg;
				break;
			case 'p':
				output_prefix = optarg;
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
	
	// Prepend the filename prefix:
	filename_bliss_ave_synergy = output_prefix + filename_bliss_ave_synergy;
	filename_bliss_stdev_synergy = output_prefix + filename_bliss_stdev_synergy;
	
	filename_loewe_ave_synergy = output_prefix + filename_loewe_ave_synergy;
	filename_loewe_stdev_synergy = output_prefix + filename_loewe_stdev_synergy;
	
	filename_pair_ave_growth = output_prefix + filename_pair_ave_growth;
	filename_pair_stdev_growth = output_prefix + filename_pair_stdev_growth;
	
	cerr << "Using samples collected from:";
		
	if(lab_selection & S_FF){
		cerr << " FF";
	}

	if(lab_selection & S_1A){
		cerr << " 1A";
	}

	if(lab_selection & S_FG){
		cerr << " FG";
	}

	cerr << endl;
}

void write_header(ofstream &m_fout, const deque<string> &m_cell_names)
{
	// For the pair drug data, use the following CSV format:
	// DRUG1,DRUG2,CELL1,CELL2,CELL3,...,CELL60
	m_fout << "DRUG1,DRUG2";

	for(deque<string>::const_iterator i = m_cell_names.begin();i != m_cell_names.end();++i){
		m_fout << ',' << *i;
	}

	m_fout << endl;
}

deque< tuple<Concentration, Concentration, float> > compute_bliss_synergy(const deque<GrowthRecord> &m_pair_data, 
	const deque<GrowthRecord> &m_drug_1_data, const deque<GrowthRecord> &m_drug_2_data)
{
	deque< tuple<Concentration, Concentration, float> > ret;
	
	for(deque<GrowthRecord>::const_iterator pair_iter = m_pair_data.begin();pair_iter != m_pair_data.end();++pair_iter){
		
		const Concentration conc_1 = pair_iter->concentration_1();
        	const Concentration conc_2 = pair_iter->concentration_2();
		
		const deque<GrowthRecord>::const_iterator drug_1_iter = find_by_conc(m_drug_1_data.begin(), m_drug_1_data.end(), conc_1);
		
		if( drug_1_iter == m_drug_1_data.end() ){
			throw __FILE__ ":compute_bliss_synergy: Unable to find drug 1 by concentration";
		}
		
		const deque<GrowthRecord>::const_iterator drug_2_iter = find_by_conc(m_drug_2_data.begin(), m_drug_2_data.end(), conc_2);
		
		if( drug_2_iter == m_drug_2_data.end() ){
			throw __FILE__ ":compute_bliss_synergy: Unable to find drug 2 by concentration";
		}
		
		// All growth values are clamped
		const float y1 = min(1.0f, drug_1_iter->test_value/drug_1_iter->control_value);
        	const float y2 = min(1.0f, drug_2_iter->test_value/drug_2_iter->control_value);

        	const float bliss_z = y1*y2;
		
		const float growth = min(1.0f, pair_iter->test_value/pair_iter->control_value);
		
		// This is the "combo score" as defined by the ALMANAC paper
		const float delta_growth = growth - bliss_z;
		
		ret.push_back( make_tuple(conc_1, conc_2, delta_growth) );
	}
	
	return ret;
}

deque< tuple<Concentration, Concentration, float> > compute_loewe_synergy(const deque<GrowthRecord> &m_pair_data, 
	const deque<GrowthRecord> &m_drug_1_data, const deque<GrowthRecord> &m_drug_2_data)
{
	deque< tuple<Concentration, Concentration, float> > ret;
	
	// Fit each set of single drug growth measurements
	DoseResponseFunction sigmoid1(0.0 /*no regularization*/);
	DoseResponseFunction sigmoid2(0.0 /*no regularization*/);

	sigmoid1.fit( log_conc_growth(m_drug_1_data) );
	sigmoid2.fit( log_conc_growth(m_drug_2_data) );
	
	for(deque<GrowthRecord>::const_iterator pair_iter = m_pair_data.begin();pair_iter != m_pair_data.end();++pair_iter){
		
		const Concentration conc_1 = pair_iter->concentration_1();
        	const Concentration conc_2 = pair_iter->concentration_2();
		
		// All growth values are clamped
		const float loewe_z = loewe_pair_growth(conc_1, sigmoid1.get_alpha(), sigmoid1.get_beta(),
			conc_2, sigmoid2.get_alpha(), sigmoid2.get_beta() );

		// Clamp the pair growth
		const float growth = min(1.0f, pair_iter->test_value/pair_iter->control_value);
		
		// This is the "combo score" as defined by the ALMANAC paper (but with the Loewe predicted growth
		// instead of the Bliss predicted growth).
		const float delta_growth = growth - loewe_z;
		
		ret.push_back( make_tuple(conc_1, conc_2, delta_growth) );
	}
	
	return ret;
}

deque< tuple<Concentration, Concentration, float> > compute_pair_growth(const deque<GrowthRecord> &m_pair_data)
{
	deque< tuple<Concentration, Concentration, float> > ret;
	
	for(deque<GrowthRecord>::const_iterator pair_iter = m_pair_data.begin();pair_iter != m_pair_data.end();++pair_iter){
		
		// Clamp the pair growth
		const float growth = min(1.0f, pair_iter->test_value/pair_iter->control_value);
		
		ret.push_back( make_tuple(pair_iter->concentration_1(), pair_iter->concentration_2(), growth) );
	}
	
	return ret;
}

deque<GrowthRecord>::const_iterator find_by_conc(const deque<GrowthRecord>::const_iterator &m_begin,
	const deque<GrowthRecord>::const_iterator &m_end, const Concentration &m_conc)
{
	for(deque<GrowthRecord>::const_iterator i = m_begin;i != m_end;++i){
		
		if(i->concentration_1() == m_conc){
			return i;
		}
	}
	
	return m_end;
}

void write_results(ofstream &m_fave, ofstream &m_fstdev, 
	const MAP< string, MAP<string, MAP<Concentration, MAP<Concentration, deque<float> > > > >  &m_data, 
	const deque<string> &m_cell_names, const ComboScoreType &m_op)
{
	// Sort the output by drug pair (to make sure that all output files have the same row ordering).
	deque<string> keys;
	
	for(MAP< string, MAP<string, MAP<Concentration, MAP<Concentration, deque<float> > > > >::const_iterator 
		i = m_data.begin();i != m_data.end();++i){
		keys.push_back(i->first);
	}
	
	make_set(keys);
	
	for(deque<string>::const_iterator i = keys.begin();i != keys.end();++i){
		
		MAP< string, MAP<string, MAP<Concentration, MAP<Concentration, deque<float> > > > >::const_iterator 
			iter = m_data.find(*i);
		
		if( iter == m_data.end() ){
			throw __FILE__ ":write_results: Unable to look up drug name";
		}
		
		m_fave << iter->first; // The drug pair names
		m_fstdev << iter->first; // The drug pair names

		for(deque<string>::const_iterator j = m_cell_names.begin();j != m_cell_names.end();++j){

			MAP<string, MAP<Concentration, MAP<Concentration, deque<float> > > >::const_iterator 
				k = iter->second.find(*j);

			if( k == iter->second.end() ){

				// There is no data for this drug pair/cell line combination
				m_fave << ',';
				m_fstdev << ',';
			}
			else{
				const pair<float, float> score = compute_score(k->second, m_op);

				// DEBUG
				//if( (fabs(score.first) > 1.0e-3) && 
				//   (fabs( fabs(score.second/score.first) - sqrt(2.0) ) < 1.0e-3) ){
				//	cout << score.first << '\t' << score.second << "\t#"
				//		<< *i << '\t' << *j << endl;
				//}
				
				m_fave << ',' << score.first;
				m_fstdev << ',' << score.second;
			}
		}

		m_fave << endl;
		m_fstdev << endl;
	}
}

pair<float, float> compute_score(const MAP<Concentration, MAP<Concentration, deque<float> > > &m_data, 
	const ComboScoreType &m_op)
{
	float ave = 0.0;
	float stdev = 0.0;

	deque< pair<float, float> > x;

	for(MAP<Concentration, MAP<Concentration, deque<float> > >::const_iterator i = m_data.begin();
		i != m_data.end();++i){
		
		for(MAP<Concentration, deque<float> >::const_iterator j = i->second.begin();
			j != i->second.end();++j){
		
			x.push_back( make_pair( 
				average<float>( j->second.begin(), j->second.end() ), 
				standard_deviation<float>( j->second.begin(), j->second.end() ) ) );
		}
	}

	const size_t len = x.size();
	
	if(len == 0){
		throw __FILE__ ":compute_score: len == 0";
	}
		
	if(m_op == COMBO_MEDIAN){
		
		sort( x.begin(), x.end() );

		if(len%2 == 1){ // Odd
		
                	ave = x[len/2].first;
			stdev = x[len/2].second;
		}
		else{
			ave = 0.5*(x[len/2 - 1].first + x[len/2].first);
			stdev = 0.5*(x[len/2 - 1].second + x[len/2].second);
		}
	}
	else{
		
		if(m_op == COMBO_MIN){
			
			sort( x.begin(), x.end() );
			
			tie(ave, stdev) = x.front();
		}
		else{
			
			if(m_op == COMBO_MAX){
			
				sort( x.begin(), x.end() );

				tie(ave, stdev) = x.back();
			}
			else{
				if(m_op == COMBO_AVERAGE){
					
					for(deque< pair<float, float> > ::const_iterator i = x.begin();
						i != x.end();++i){
						
						ave += i->first;
						stdev += i->second*i->second;
					}
					
					ave /= len;
					stdev = sqrt(stdev/len);
				}
				else{
					throw __FILE__ ":compute_score: Unknown combo score!";
				}
			}
		}
	}
	
	return make_pair(ave, stdev);
}

void accumulate_data(MAP<Concentration, MAP<Concentration, deque<float> > > &m_data, 
	const deque< tuple<Concentration, Concentration, float> > &m_sample)
{
	for(deque< tuple<Concentration, Concentration, float> >::const_iterator i = m_sample.begin();i != m_sample.end();++i){
				
		deque<float> &ref = m_data[get<0>(*i)][get<1>(*i)];

		ref.push_back( get<2>(*i) );
	}
}
