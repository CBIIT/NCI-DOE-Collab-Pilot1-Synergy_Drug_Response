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

#ifndef __RANDOM_FOREST
#define __RANDOM_FOREST

#include <vector>
#include <deque>
#include <random>
#include <unordered_map>

struct LabeledData
{
	LabeledData()
	{
		partition_left = false;
	};
	
	float response;
	bool partition_left;
	
	// An integer index for fast lookup of the drug and
	// cell line features
	unsigned int drug_index;
	unsigned int drug_index_2;
	unsigned int cell_index;
	
	// The string-based name of each drug and cell
	std::string drug_id;
	std::string drug_id_2;
	std::string cell_id;
};

class RandomForest
{	
	public:
		typedef enum {
			CELL_FEATURE,
			DRUG_FEATURE,
			DRUG_PAIR_MIN_FEATURE,
			DRUG_PAIR_MAX_FEATURE,
			DRUG_PAIR_AVE_FEATURE,
			UNKNOWN_FEATURE
		} FeatureType;

		struct TreeNode
		{

			TreeNode()
			{
				left = right = 0;
			};


			std::pair<unsigned int /*index*/, 
				float /*threshold*/> boundary;
			FeatureType boundary_feature;

			unsigned int left;
			unsigned int right;

			float prediction;

			inline bool is_leaf() const
			{
				return left == right;
			};
		};


	private:
		
		typedef std::deque<TreeNode> Tree;
		
		std::vector<Tree> forest;
		
		size_t forest_size;
		size_t forest_leaf;
		float forest_feature_bag;
		float forest_data_bag;
		unsigned int *ptr_seed;
		
		void build_tree(Tree &m_tree, 
			std::vector<LabeledData>::iterator m_begin,
			std::vector<LabeledData>::iterator m_end,
			const std::vector< std::vector<float> > &m_cell_features,
			const std::vector< std::vector<float> > &m_drug_features,
			const std::vector< std::pair<unsigned int, 
				FeatureType> > &m_feature_indicies);
			
		float predict_tree(const Tree &m_tree,
			const std::vector<float> &m_cell_features, 
			const std::vector<float> &m_drug_features) const;
			
		float predict_tree(const Tree &m_tree,
			const std::vector<float> &m_cell_features, 
			const std::vector<float> &m_drug_features,
			const std::vector<float> &m_drug_features_2) const;
		
	public:
	
		RandomForest(const size_t &m_forest_size, 
			const size_t &m_forest_leaf, 
			const float &m_forest_feature_bag,
			const float &m_forest_data_bag,
			unsigned int *m_seed_ptr) :
			forest_size(m_forest_size),
			forest_leaf(m_forest_leaf),
			forest_feature_bag(m_forest_feature_bag),
			forest_data_bag(m_forest_data_bag),
			ptr_seed(m_seed_ptr)
		{
		};
		
		// m_data needs to be mutable to allow the tree building process
		// to reorder the data points "in place" as individual trees
		// are built.
		void build(std::vector<LabeledData> &m_data,
			const std::unordered_map<std::string, std::vector<float> > &m_cell_features,
			const std::unordered_map<std::string, std::vector<float> > &m_drug_features);
		
		float predict(
			const std::vector<float> &m_cell_features, 
			const std::vector<float> &m_drug_features) const;
		
		float predict(
			const std::vector<float> &m_cell_features, 
			const std::vector<float> &m_drug_features,
			const std::vector<float> &m_drug_features_2) const;
};

#endif // __RANDOM_FOREST
