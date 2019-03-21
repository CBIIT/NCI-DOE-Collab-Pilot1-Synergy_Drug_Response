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

#include "random_forest.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <list>
#include <math.h>
#include <mpi.h>

#include "shuffle.h"
#include "mpi_util.h"

// A way to potentially speed up the sorting and partitioning routines (which can consume a non-trivial
// amount of CPU time) is to use the g++-specific __gnu_parallel:: funcions instead
// of the std:: functions. Since __gnu_parallel:: is a drop in replacement for
// std::, include:
#include <parallel/algorithm>
// which defines _GLIBCXX_PARALLEL_PARALLEL_H and enables parallel sorting
// Please note that approach will not work on non-g++ compilers (i.e. clang or intel).

#if defined(_GLIBCXX_PARALLEL_PARALLEL_H) && defined(_OPENMP)

	#define	SORT		__gnu_parallel::sort
	#define	PARTITION	__gnu_parallel::partition
#else

	#define	SORT		std::sort
	#define	PARTITION	std::partition
#endif // _GLIBCXX_PARALLEL_PARALLEL_H

// If USE_MIN_MAX_PAIR_FEATURES is defined, we represent drug pair features with
// the per-element min and max of each input feature vaule. If USE_MIN_MAX_PAIR_FEATURES is
// *NOT* defined, then the the drug pair features are represented with the per-element, average 
// feature value.
#define	USE_MIN_MAX_PAIR_FEATURES

using namespace std;

// Global variables for MPI
extern int mpi_numtasks;
extern int mpi_rank;

float average_response(
	const vector<LabeledData>::const_iterator m_begin,
	const vector<LabeledData>::const_iterator m_end);	

float best_split(pair<unsigned int, float> &m_boundary, vector<unsigned int> &m_left, 
	RandomForest::FeatureType &m_best_feature,
	const size_t &m_leaf, 
	const vector<LabeledData>::const_iterator &m_begin, 
	const vector<LabeledData>::const_iterator &m_end, 
	const vector< vector<float> > &m_cell_features,
	const vector< vector<float> > &m_drug_features,
	const vector< pair<unsigned int, RandomForest::FeatureType> > &m_feature_indicies);
			
void RandomForest::build(vector<LabeledData> &m_data, 
	const unordered_map<string, vector<float> > &m_cell_features,
	const unordered_map<string, vector<float> > &m_drug_features)
{
	
	if(forest_size == 0){
		
		// when number of trees is less than the number if workers, an
		// individual work may have nothing to do.
		return;
	}
	
	if( (forest_feature_bag <= 0.0) || (forest_feature_bag > 1.0) ){
		throw __FILE__ ":RandomForest::build: Please specify 0 < forest_feature_bag <= 1.0";
	}
	
	if( (forest_data_bag <= 0.0) || (forest_data_bag > 1.0) ){
		throw __FILE__ ":RandomForest::build: Please specify 0 < forest_data_bag <= 1.0";
	}
	
	// Allow integer-based indexing of the feature vectors (hash table lookups are too slow!).
	unordered_map<string, unsigned int> cell_index;
	unordered_map<string, unsigned int> drug_index;
	
	const size_t num_cells = m_cell_features.size();
	const size_t num_drugs = m_drug_features.size();
	const size_t num_data = m_data.size();
	
	// Make a local copy of the drug and cell features. The extra overhead of this copy is very small
	// compared to the benefit of integer-based indexing (and not having to perform string-based
	// hash table lookups).
	vector< vector<float> > cell_features(num_cells);
	vector< vector<float> > drug_features(num_drugs);
	
	unsigned int index = 0;
	
	for(unordered_map<string, vector<float> >::const_iterator i = m_cell_features.begin();
		i != m_cell_features.end();++i){
		
		cell_index[i->first] = index;
		cell_features[index] = i->second;
		++index;
	}
	
	index = 0;
	
	for(unordered_map<string, vector<float> >::const_iterator i = m_drug_features.begin();
		i != m_drug_features.end();++i){
		
		drug_index[i->first] = index;
		drug_features[index] = i->second;
		++index;
	}
	
	// Make sure that all of the data has the same number of features (for both cell and drug)
	// and the same number of response values
	unsigned int num_cell_features = 0;
	unsigned int num_drug_features = 0;
	
	size_t num_drug_pair = 0;
	
	for(vector<LabeledData>::iterator i = m_data.begin();i != m_data.end();++i){
		
		// Convert from the string-based indexing (which requires expensive hash-table lookups)
		if( !cell_index.empty() ){
		
			unordered_map<string, unsigned int>::const_iterator index_iter = 
				cell_index.find(i->cell_id);

			if( index_iter == cell_index.end() ){
				throw __FILE__ ":RandomForest::build: Unable to lookup cell index";
			}

			i->cell_index = index_iter->second;
			
			if( i == m_data.begin() ){
				num_cell_features = cell_features[i->cell_index].size();
			}
			
			if( num_cell_features != cell_features[i->cell_index].size() ){
				throw __FILE__ ":RandomForest::build: Variable number of cell features not allowed";
			}
		}
		
		if( !drug_index.empty() ){
		
			unordered_map<string, unsigned int>::const_iterator index_iter = 
				drug_index.find(i->drug_id);
		
			if( index_iter == drug_index.end() ){
				throw __FILE__ ":RandomForest::build: Unable to lookup drug index";
			}
		
			i->drug_index = index_iter->second;
		
			if( i == m_data.begin() ){
				num_drug_features = drug_features[i->drug_index].size();
			}
		
			if( num_drug_features != drug_features[i->drug_index].size() ){
				throw __FILE__ ":RandomForest::build: Variable number of drug features not allowed";
			}
			
			// Does this data point refer to a pair of points?
			if( !i->drug_id_2.empty() ){

				unordered_map<string, unsigned int>::const_iterator index_iter = 
					drug_index.find(i->drug_id_2);

				if( index_iter == drug_index.end() ){
					throw __FILE__ ":RandomForest::build: Unable to lookup drug index (2)";
				}

				i->drug_index_2 = index_iter->second;

				++num_drug_pair;
			}
		}
	}
	
	if( (num_drug_pair != 0) && ( num_drug_pair != m_data.size() ) ){
		throw __FILE__ ":RandomForest::build: Mixed pair and single drug data not allowed!";
	}

	const size_t num_single_drug_features = num_drug_features;
	
	#ifdef USE_MIN_MAX_PAIR_FEATURES
	if(num_drug_pair > 0){
		
		// Drug pairs require twice as many features
		num_drug_features *= 2;
	}
	#endif // USE_MIN_MAX_PAIR_FEATURES
	
        unsigned int num_bagged_cell_features = 
		max( (unsigned int)(1), (unsigned int)(forest_feature_bag*num_cell_features) );
	
	// Since the number of bagged cell features was clamped to be >= 1, we need to 
	// explicitly test for the case that num_cell_features is zero!
	if(num_cell_features == 0){
		num_bagged_cell_features = 0;
	}
	
	unsigned int num_bagged_drug_features = 
		max( (unsigned int)(1), (unsigned int)(forest_feature_bag*num_drug_features) );
	
	// Since the number of bagged drug features was clamped to be >= 1, we need to 
	// explicitly test for the case that num_drug_features is zero!
	if(num_drug_features == 0){
		num_bagged_drug_features = 0;
	}
	
	const size_t num_bagged_data = 
		max( (size_t)(1), (size_t)(forest_data_bag*num_data) );
	
	vector<unsigned char> cell_mask(num_cell_features, false);
	vector<unsigned char> drug_mask(num_drug_features, false);

	for(unsigned int i = 0;i < num_bagged_cell_features;++i){
		cell_mask[i] = true;
	}	

	for(unsigned int i = 0;i < num_bagged_drug_features;++i){
		drug_mask[i] = true;
	}
	
	vector< pair<unsigned int, FeatureType> > feature_indicies;
	
	feature_indicies.reserve(num_bagged_cell_features
		+ num_bagged_drug_features);
	
	// The number of trees that this rank is responsible for building
	const size_t local_forest_size = forest_size/mpi_numtasks + 
		( int(forest_size%mpi_numtasks) > mpi_rank ? 1 : 0) ;

	forest.resize(local_forest_size);

	// Keep the user informed
	string info_buffer;
	
	if(mpi_rank == 0){
		cerr << "\tForest growth: ";
	}
	
	for(size_t i = 0;i < local_forest_size;++i){
		
		randomize(cell_mask.begin(), cell_mask.end(), ptr_seed);
		randomize(drug_mask.begin(), drug_mask.end(), ptr_seed);
		
		// Randomize the *input* data for every tree if we are
		// bagging data
		if(num_bagged_data < num_data){
			randomize(m_data.begin(), m_data.end(), ptr_seed);
		}
		
		// While we could just pass the feature masks, passing
		// the actual target indicies was intended to improve
		// OpenMPI-based multi-threaded load balancing
		feature_indicies.clear();
		
		for(unsigned int j = 0;j < num_cell_features;++j){
			
			if(cell_mask[j]){
				feature_indicies.push_back( make_pair(j, CELL_FEATURE) );
			}
		}
		
		for(unsigned int j = 0;j < num_drug_features;++j){
			
			if(drug_mask[j]){
			
				if(num_drug_pair > 0){
					
					// Drug pair
					#ifdef USE_MIN_MAX_PAIR_FEATURES
					if(j < num_single_drug_features){
						feature_indicies.push_back( 
							make_pair(j, 
								DRUG_PAIR_MIN_FEATURE) );
					}
					else{
						feature_indicies.push_back( 
							make_pair(j - num_single_drug_features, 
								DRUG_PAIR_MAX_FEATURE) );
					}
					#else
					feature_indicies.push_back( 
							make_pair(j, 
								DRUG_PAIR_AVE_FEATURE) );
					#endif // USE_MIN_MAX_PAIR_FEATURES
					
				}
				else{
					// Single drug
					feature_indicies.push_back( make_pair(j, DRUG_FEATURE) );
				}
			}
		}
		
		RandomForest::build_tree(forest[i], 
			m_data.begin(), m_data.begin() + num_bagged_data, 
			cell_features, drug_features,
			feature_indicies);

		if(mpi_rank == 0){
		
			// Update the forest building process
			for(string::const_iterator j = info_buffer.begin();j != info_buffer.end();++j){
				cerr << '\b';
			}

			for(string::const_iterator j = info_buffer.begin();j != info_buffer.end();++j){
                        	cerr << ' ';
                	}

			for(string::const_iterator j = info_buffer.begin();j != info_buffer.end();++j){
                        	cerr << '\b';
                	}

			stringstream ssin;

			ssin << (100.0*i)/local_forest_size << '%';

			info_buffer = ssin.str();

			cerr << info_buffer;
		}
	}
	
	if(mpi_rank == 0){
		cerr << endl;
	}
}

struct IsLeft
{
	inline bool operator()(const LabeledData &m_x) const
	{
		return m_x.partition_left;
	};
};

void RandomForest::build_tree(Tree &m_tree, 
	vector<LabeledData>::iterator m_begin,
	vector<LabeledData>::iterator m_end,
	const vector< vector<float> > &m_cell_features,
	const vector< vector<float> > &m_drug_features,
	const vector< pair<unsigned int, FeatureType> > &m_feature_indicies)
{	
	if(m_end <= m_begin){
		throw __FILE__ ":RandomForest::build_tree: No data!";
	}
	
	const unsigned int num_data = m_end - m_begin;
	
	TreeNode local;
	
	m_tree.push_back(local);
	
	// Do we have enough data to make a split?
	if(num_data <= forest_leaf){
		
		// Return the average over all leaf members
		m_tree.back().prediction = average_response(m_begin, m_end);
		return;
	}
	
	// Search for the partition that obtains the smallest *weighted* mean square error
	pair<unsigned int, float> best_boundary;
	vector<unsigned int> best_left;
	FeatureType best_feature = RandomForest::UNKNOWN_FEATURE;
	
	best_split(best_boundary, best_left, best_feature, forest_leaf, 
		m_begin, m_end, 
		m_cell_features, m_drug_features, 
		m_feature_indicies);
	
	// We could not find a valid split
	if(best_feature == RandomForest::UNKNOWN_FEATURE){

		// Return the average over all leaf members
		m_tree.back().prediction = average_response(m_begin, m_end);
		return;
	}
	
	// Partition the data into left and right branches.
	for(vector<unsigned int>::iterator i = best_left.begin();i != best_left.end();++i){
		(m_begin + *i)->partition_left = true;
	}
	
	PARTITION( m_begin, m_end, IsLeft() );
	
	vector<LabeledData>::iterator boundary_iter = m_begin + best_left.size();
	
	// Unset the partition flag bit so we can use the cell_id variable normally
	for(vector<LabeledData>::iterator i = m_begin;i != boundary_iter;++i){
		i->partition_left = false;
	}
	
	const unsigned int node = m_tree.size() - 1;
	
	m_tree[node].boundary = best_boundary;
	m_tree[node].boundary_feature = best_feature;
	
	m_tree[node].left = m_tree.size();

	build_tree(m_tree, 
		m_begin, boundary_iter,
		m_cell_features, m_drug_features,
		m_feature_indicies);

	m_tree[node].right = m_tree.size();

	build_tree(m_tree, 
		boundary_iter, m_end,
		m_cell_features, m_drug_features,
		m_feature_indicies);
}

float RandomForest::predict(
	const vector<float> &m_cell_features, 
	const vector<float> &m_drug_features) const
{
	// forest_size refers to the *global* (not local) forest. This size
	// should never be zero (evern when the number of trees assigned to
	// this MPI rank is zero).
	if(forest_size == 0){
		throw __FILE__ ":RandomForest::predict: forest_size == 0";
	}

	// Input variables are taken from rank 0
	vector<float> cell_features(m_cell_features);
	vector<float> drug_features(m_drug_features);
	
	broadcast(cell_features, mpi_rank, 0);
	broadcast(drug_features, mpi_rank, 0);
	
	double sum = 0.0;
	
	for(vector<Tree>::const_iterator i = forest.begin();i != forest.end();++i){
		sum += predict_tree(*i, cell_features, drug_features);
	}

	MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
					
	return sum/forest_size;
}

float RandomForest::predict(
	const vector<float> &m_cell_features, 
	const vector<float> &m_drug_features,
	const vector<float> &m_drug_features_2) const
{
	
	// forest_size refers to the *global* (not local) forest. This size
	// should never be zero (evern when the number of trees assigned to
	// this MPI rank is zero).
	if(forest_size == 0){
		throw __FILE__ ":RandomForest::predict: forest_size == 0";
	}
	
	// Input variables are taken from rank 0
	vector<float> cell_features(m_cell_features);
	vector<float> drug_features(m_drug_features);
	vector<float> drug_features_2(m_drug_features_2);
	
	broadcast(cell_features, mpi_rank, 0);
	broadcast(drug_features, mpi_rank, 0);
	broadcast(drug_features_2, mpi_rank, 0);
	
	double sum = 0.0;
	
	for(vector<Tree>::const_iterator i = forest.begin();i != forest.end();++i){
		sum += predict_tree(*i, cell_features, drug_features, drug_features_2);
	}
	
	MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return sum/forest_size;
}

float RandomForest::predict_tree(const Tree &m_tree,
	const vector<float> &m_cell_features, 
	const vector<float> &m_drug_features,
	const vector<float> &m_drug_features_2) const
{
	if( m_tree.empty() ){
		throw __FILE__ ":RandomForest::predict_tree: Empty tree!";
	}
	
	const unsigned int num_cell_features = m_cell_features.size();
	const unsigned int num_drug_features = m_drug_features.size();
	const unsigned int num_drug_features_2 = m_drug_features_2.size();
	
	if(num_drug_features != num_drug_features_2){
		throw __FILE__ ":RandomForest::predict_tree: |drug 1 features| != |drug 2 features|";
	}
	
	unsigned int index = 0;
	
	while(true){

		const TreeNode &node = m_tree[index];

		if( node.is_leaf() ){
			return node.prediction;
		}
		
		switch(node.boundary_feature){
			case CELL_FEATURE:
			
				if(num_cell_features <= node.boundary.first){
					throw __FILE__ ":RandomForest::predict_tree: Cell feature index out of bounds!";
				}

				index = (m_cell_features[node.boundary.first] < node.boundary.second) ?
					node.left :
					node.right;					
				break;
			case DRUG_PAIR_MIN_FEATURE:
			
				if( (num_drug_features <= node.boundary.first) || 
				    (num_drug_features_2 <= node.boundary.first) ){
					throw __FILE__ ":RandomForest::predict_tree: Min drug pair feature index out of bounds!";
				}

				index = ( min(m_drug_features[node.boundary.first], m_drug_features_2[node.boundary.first]) 
						< node.boundary.second) ?
					node.left :
					node.right;
				break;
			case DRUG_PAIR_MAX_FEATURE:
			
				if( (num_drug_features <= node.boundary.first) || 
				    (num_drug_features_2 <= node.boundary.first) ){
					throw __FILE__ ":RandomForest::predict_tree: Max drug pair feature index out of bounds!";
				}

				index = ( max(m_drug_features[node.boundary.first], m_drug_features_2[node.boundary.first]) 
						< node.boundary.second) ?
					node.left :
					node.right;				
				break;
			case DRUG_PAIR_AVE_FEATURE:
			
				if(num_drug_features <= node.boundary.first){
					throw __FILE__ ":RandomForest::predict_tree: Drug feature index out of bounds!";
				}

				index = ( 0.5*(m_drug_features[node.boundary.first] + m_drug_features_2[node.boundary.first])
					< node.boundary.second) ?
					node.left :
					node.right;			
				break;
			default:
				throw __FILE__ ":RandomForest::predict_tree: Unknown boundary feature type!";
		};
	}
	
	throw __FILE__ ":RandomRegressionForest::predict_tree: Should never get here!";
	
	return 0.0f;
}

float RandomForest::predict_tree(const Tree &m_tree,
	const vector<float> &m_cell_features, 
	const vector<float> &m_drug_features) const
{
	if( m_tree.empty() ){
		throw __FILE__ ":RandomForest::predict_tree: Empty tree!";
	}
	
	const unsigned int num_cell_features = m_cell_features.size();
	const unsigned int num_drug_features = m_drug_features.size();
	
	unsigned int index = 0;

	while(true){

		const TreeNode &node = m_tree[index];

		if( node.is_leaf() ){
			return node.prediction;
		}
		
		switch(node.boundary_feature){
			case CELL_FEATURE:
			
				if(num_cell_features <= node.boundary.first){
					throw __FILE__ ":RandomForest::predict_tree: Cell feature index out of bounds!";
				}

				index = (m_cell_features[node.boundary.first] < node.boundary.second) ?
					node.left :
					node.right;
				break;
			case DRUG_FEATURE:
			
				if(num_drug_features <= node.boundary.first){
					throw __FILE__ ":RandomForest::predict_tree: Drug feature index out of bounds!";
				}

				index = (m_drug_features[node.boundary.first] < node.boundary.second) ?
					node.left :
					node.right;

				break;
			default:
				throw __FILE__ ":RandomForest::predict_tree: Unknown boundary feature type!";
		};
	}
	
	throw __FILE__ ":RandomRegressionForest::predict_tree: Should never get here!";
	
	return 0.0f;
}

float average_response(
	const vector<LabeledData>::const_iterator m_begin,
	const vector<LabeledData>::const_iterator m_end)
{	
	if(m_end <= m_begin){
		throw __FILE__ ":average_response: No data!";
	}
	
	const unsigned int N = m_end - m_begin;
	
	// Accumulate as double to avoid loss of precision for
	// large datasets
	double ret = 0.0;
	
	#pragma omp parallel for \
		reduction(+:ret)
	for(unsigned int i = 0;i < N;++i){
		ret += (m_begin + i)->response;
	}
	
	ret /= N;
	
	return ret;
}

float best_split(pair<unsigned int, float> &m_boundary, vector<unsigned int> &m_left, 
	RandomForest::FeatureType &m_best_feature,
	const size_t &m_leaf, 
	const vector<LabeledData>::const_iterator &m_begin, 
	const vector<LabeledData>::const_iterator &m_end, 
	const vector< vector<float> > &m_cell_features,
	const vector< vector<float> > &m_drug_features,
	const vector< pair<unsigned int, RandomForest::FeatureType> > &m_feature_indicies)
{
	#define		REALLY_LARGE_VALUE	1.0e20
	
	float best_score = REALLY_LARGE_VALUE; // <-- A really large value!
	
	if(m_end <= m_begin){
		throw __FILE__ ":best_split: No data!";
	}
	
	const unsigned int num_data = m_end - m_begin;
	
	if(m_leaf < 1){
		throw __FILE__ ":best_split: m_leaf < 1";
	}

	const unsigned int num_features_to_test = m_feature_indicies.size();
	
	#pragma omp parallel
	{
		// Store the values for the i^th independent feature. We can allocate
		// this memory outside the for loop, since the size does not change
		vector< pair<float, unsigned int> > feature_slice(num_data);

		float local_best_score = REALLY_LARGE_VALUE; // <-- A really large value!
		RandomForest::FeatureType local_best_feature = RandomForest::UNKNOWN_FEATURE;
		pair<unsigned int, float> local_boundary;
		vector<unsigned int> local_left;
		
		// Test each feature and every possible boundary value within a feature 
		#pragma omp for
		for(unsigned int f = 0;f < num_features_to_test;++f){

			const pair<unsigned int, RandomForest::FeatureType> &curr_feature = m_feature_indicies[f];
			
			switch(curr_feature.second){
				case RandomForest::CELL_FEATURE:
					
					for(unsigned int i = 0;i < num_data;++i){
						feature_slice[i] = 
							make_pair(m_cell_features[(m_begin + i)->cell_index][curr_feature.first], i);
					}
					break;
				case RandomForest::DRUG_FEATURE:
				
					for(unsigned int i = 0;i < num_data;++i){
						
						feature_slice[i] = 
							make_pair(m_drug_features[ (m_begin + i)->drug_index ][curr_feature.first], i);
					}
					break;
				case RandomForest::DRUG_PAIR_MIN_FEATURE:
				
					for(unsigned int i = 0;i < num_data;++i){
						
						const float min_value = 
							min(
								m_drug_features[ (m_begin + i)->drug_index ][curr_feature.first],
								m_drug_features[ (m_begin + i)->drug_index_2 ][curr_feature.first]);
						
						feature_slice[i] = 
							make_pair(min_value, i);
					}
					break;
				case RandomForest::DRUG_PAIR_MAX_FEATURE:
				
					for(unsigned int i = 0;i < num_data;++i){
						
						const float max_value = 
							max(
								m_drug_features[ (m_begin + i)->drug_index ][curr_feature.first],
								m_drug_features[ (m_begin + i)->drug_index_2 ][curr_feature.first]);
						
						feature_slice[i] = 
							make_pair(max_value, i);
					}
					break;
				case RandomForest::DRUG_PAIR_AVE_FEATURE:
				
					for(unsigned int i = 0;i < num_data;++i){
						
						const float ave_value = 
							0.5*(
								m_drug_features[ (m_begin + i)->drug_index ][curr_feature.first] +
								m_drug_features[ (m_begin + i)->drug_index_2 ][curr_feature.first]);
						
						feature_slice[i] = 
							make_pair(ave_value, i);
					}
					break;
				default:
					throw __FILE__ ":best_split: Unknown feature type!";
			};

			// Sort the feature values in ascending order. Since we are already in
			// a parallel section, use a serial sort
			sort( feature_slice.begin(), feature_slice.end() );

			// To make the calculation efficient, track the running sum of y values in the left and 
			// right branches. Use double to accumulate the floating point moments
			double sum_left = 0.0f;
			double sum_right = 0.0f;

			double sum_left2 = 0.0;
			double sum_right2 = 0.0;

			unsigned int num_left = 0;
			unsigned int num_right = 0;

			// By default, all of the data starts in the *right* branch
			
			num_right = num_data;
			
			for(unsigned int i = 0;i < num_data;++i){

				const LabeledData &ref = *(m_begin + i);
				
				sum_right += ref.response;
				sum_right2 += ref.response*ref.response;
			}
			
			for(unsigned int i = 0;i < num_data;++i){

				// Move data point i from the right branch to the left branch
				const float y = (m_begin + feature_slice[i].second)->response;

				// Since the sums are unnormalized, we can simply remove a point from the right and
				// add it to the left
				sum_right -= y;
				sum_left += y;

				const float y2 = y*y;
				
				sum_right2 -= y2;
				sum_left2 += y2;

				--num_right;
				++num_left;

				if( (i < m_leaf) || (i >= (num_data - m_leaf) ) ){
					continue;
				}

				// Don't split on equal values! We can access element i + 1 of the
				// feature slice vector since m_leaf must be greater than 0 (and the
				// above test on i will trigger a "continue" at the end of the vector range)
				if(feature_slice[i].first == feature_slice[i + 1].first){
					continue;
				}
				
				if( (num_left <= 0) || (num_right <= 0) ){
					throw __FILE__ ":best_split: Unable to normalize split variance";
				}

				const float left_mean_square_error = (sum_left2 - sum_left*sum_left/num_left);
				const float right_mean_square_error = (sum_right2 - sum_right*sum_right/num_right);

				// Here is the original, more readable code:
				//left_mean_square_error /= num_left;
				//right_mean_square_error /= num_right;
				//
				//const float trial_mean_square_error = (num_left*left_mean_square_error + 
				//	num_right*right_mean_square_error)/num_data;

				// And here is the code when we cancel the factors of L and R
				const float trial_mean_square_error = 
					(left_mean_square_error + right_mean_square_error)/num_data;

				if(trial_mean_square_error < local_best_score){

					local_best_score = trial_mean_square_error;
					local_best_feature = curr_feature.second;
					
					// Place the boundary at the midpoint between the current and the
					// next feature value.
					local_boundary = make_pair(curr_feature.first, 
						0.5*(feature_slice[i].first + feature_slice[i + 1].first) );

					local_left.resize(i + 1);

					for(unsigned int j = 0;j <= i;++j){
						local_left[j] = feature_slice[j].second;
					}
					
					// If a data point is not in the left branch, it must be in the
					// right branch (so we don't need to explicitly store the data points
					// that belong to the right hand branch).
				}
			}
		}
		
		#pragma omp critical
		if(local_best_score < best_score){
			
			best_score = local_best_score;
			m_best_feature = local_best_feature;
			m_boundary = local_boundary;
			m_left = local_left;
		}
	}
	
	return best_score;
}
