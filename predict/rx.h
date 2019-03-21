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

#ifndef __RX
#define __RX

#include <unordered_map>
#include <string>
#include <vector>
#include <deque>

#define	MAP 		std::unordered_map
#define	MULTIMAP 	std::unordered_multimap

// Allowed values for the synery matrix (in counting mode)

#define	SYNERGY_NOT_TESTED	0.0
#define	SYNERGY_YES		1.0
#define	SYNERGY_NO		-1.0

#define	MISSING_DATA		0xFFFFFFFF

// Keep track of changes to this program
#define	RX_VERSION		"0.6"

#define	DEFAULT_NUM_TREE	100000
#define	DEFAULT_LEAF_SIZE	3
#define	DEFAULT_BAG_FRACTION	0.3
#define	DEFAULT_NUM_FOLD	5


typedef enum {
	CLUSTER_NONE,
	CLUSTER_SINGLE_LINKAGE,
	CLUSTER_TOTAL_LINKAGE
} ClusteringAlgorithm;

typedef enum {
	AVERAGE_SCORE,
	MEDIAN_SCORE,
	FILTERED_MEAN_SCORE,
	COUNTING_SCORE,
	UNDEFINED_SCORE
} ScoreFunction;

struct SynergyPrediction
{
	float ave;
	float stdev;
	std::string drug_id;
	
	SynergyPrediction()
	{
		ave = stdev = 0.0;
	};
	
	inline bool operator<(const SynergyPrediction &m_rhs) const
	{
		return (ave < m_rhs.ave);
	};
};

struct Options
{
	bool print_usage;
	size_t num_tree;
	size_t leaf_size;
	float bag_fraction;

	size_t num_fold;

	std::string file_name_drug;
	std::string file_name_blind_drug;
	std::string file_name_drug_target;
	bool randomize_drug;

	std::string file_name_cell;
	std::string file_name_blind_cell;
	bool randomize_cell;

	MAP<std::string, std::string> train_cell_to_synergy_file;
	MAP<std::string, std::string> test_cell_to_synergy_file;
	bool randomize_test_synergy;
	bool randomize_train_synergy;
	bool NSC_to_CID_y;
	ScoreFunction score_func;
	float synergy_threshold;

	std::string file_name_output;

	unsigned int seed; // 0 -> use a time based seed
	
	bool use_pair_prediction;
	
	bool one_hot_cell; // One hot encode cell lines
	bool one_hot_drug; // One hot encode drugs
	bool one_hot_tissue; // One hot encode cell line tissue type
	
	bool normalize_cell;
	
	// When the test and training sets are different (i.e. ALMANAC
	// and Merck), what should we do with the overlapping cell lines
	// and drugs?
	bool overlap_to_train;
	bool overlap_to_test;
	
	Options(int argc, char *argv[])
	{
		parse(argc, argv);
	};
	
	Options()
	{
	
	};
	
	void parse(int argc, char *argv[]);
	
	template<class T> friend 
		size_t mpi_size(const T &m_obj);
	template<class T> friend 
		unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);
	template<class T> friend 
		unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
};

template<> size_t mpi_size(const Options &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const Options &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, Options &m_obj);

// In io.cpp
MAP< std::string, std::vector<float> > parse_fingerprint(const std::string &m_filename,
	std::vector<std::string> &m_feature_id,
	const std::string &m_drug, const std::string &m_fingerprint,
	const char m_delim = ',');

MAP< std::string, std::vector<float> > parse_descriptor(const std::string &m_filename,
	std::vector<std::string> &m_feature_id, const std::string &m_drug,
	const char m_delim = ',');
	
MAP<std::string, std::string> parse_label(const std::string &m_filename,
	const std::string &m_drug, const std::string &m_label);

MAP<std::string, float> parse_response(const std::string &m_filename,
	const std::string &m_drug, const std::string &m_response);

MAP< std::string, MAP<std::string, float> > parse_matrix(const std::string &m_filename);

bool is_fingerprint(const std::string &m_filename);

MAP<std::string, std::string> parse_drug_target(const std::string &m_filename);

// In cluster.cpp
std::deque< std::deque<std::string> > cluster(const std::vector<std::string> &m_drugs, 
	const MAP< std::string, std::vector<float> > &m_data,
	const float &m_threshold, const ClusteringAlgorithm &m_alg);

// In one_hot.cpp
size_t one_hot_encode(MAP<std::string, std::vector<float> > &m_features, 
	const std::vector<std::string> &m_categories);
	
std::vector<std::string> one_hot_encode_by_tissue(MAP<std::string, std::vector<float> > &m_features, 
	const std::vector<std::string> &m_categories);

#endif // __RX
