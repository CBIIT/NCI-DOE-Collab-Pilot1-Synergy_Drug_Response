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

#include <math.h>

#include <iostream>

using namespace std;

deque< deque<string> > linkage_cluster(const vector<string> &m_drugs, 
	const MAP<string, vector<float> > &m_data,
	const float &m_threshold,
	const ClusteringAlgorithm &m_alg);

deque< deque<string> > separate_clusters(const vector<string> &m_drugs);

float cluster_distance(const deque<string> &m_A, const deque<string> &m_B, 
	const MAP<string, vector<float> > &m_data, const ClusteringAlgorithm &m_alg);

float drug_distance(const vector<float> &m_A, const vector<float> &m_B);

deque< deque<string> > cluster(const vector<string> &m_drugs, 
	const MAP<string, vector<float> > &m_data,
	const float &m_threshold,
	const ClusteringAlgorithm &m_alg)
{
	switch(m_alg){
		case CLUSTER_NONE:
			
			// Every drug is in a separate cluster
			return separate_clusters(m_drugs);

		case CLUSTER_SINGLE_LINKAGE:
		case CLUSTER_TOTAL_LINKAGE:
			return linkage_cluster(m_drugs, m_data, m_threshold, m_alg);
		default:
			throw __FILE__ ":cluster: Unknown clustering algorithm requested";
	};
	
	throw __FILE__ ":cluster: Shouldn't get here!";
	return deque< deque<string> >();
}

deque< deque<string> > separate_clusters(const vector<string> &m_drugs)
{
	deque< deque<string> > ret;
	
	// Each drug starts in its own cluster
	for(vector<string>::const_iterator i = m_drugs.begin();i != m_drugs.end();++i){
		
		ret.push_back( deque<string>() );
		ret.back().push_back(*i);
	}
	
	return ret;
}

deque< deque<string> > linkage_cluster(const vector<string> &m_drugs, 
	const MAP<string, vector<float> > &m_data, const float &m_threshold,
	const ClusteringAlgorithm &m_alg)
{
	deque< deque<string> > ret;
	
	// Each drug starts in its own cluster
	for(vector<string>::const_iterator i = m_drugs.begin();i != m_drugs.end();++i){
		
		ret.push_back( deque<string>() );
		ret.back().push_back(*i);
	}
	
	// Iteratively merge the closest clusters
	while(true){
		
		const size_t num_cluster = ret.size();
		pair<size_t, size_t> best_pair(0, 0);
		float best_distance = 0.0;
		
		for(size_t i = 0;i < num_cluster;++i){
			
			for(size_t j = i + 1;j < num_cluster;++j){
				
				const float d = cluster_distance(ret[i], ret[j], m_data, m_alg);
				
				if( (d < best_distance) || (best_pair.first == best_pair.second) ){
					
					best_distance = d;
					best_pair = make_pair(i, j);
				}
			}
		}
		
		if( (best_distance > m_threshold) || (best_pair.first == best_pair.second) ){
			
			// No clusters were merged
			break;
		}
		
		// Merge best_pair.second into best_pair.first
		for(deque<string>::const_iterator i = ret[best_pair.second].begin();i != ret[best_pair.second].end();++i){
			ret[best_pair.first].push_back(*i);
		}
		
		// Erase cluster best_pair.second
		ret.erase(ret.begin() + best_pair.second);
	}
	
	return ret;
}

float cluster_distance(const deque<string> &m_A, const deque<string> &m_B, 
	const MAP<string, vector<float> > &m_data, const ClusteringAlgorithm &m_alg)
{
	float min_distance = 0.0;
	float max_distance = 0.0;
	
	for(deque<string>::const_iterator i = m_A.begin();i != m_A.end();++i){
		
		MAP<string, vector<float> >::const_iterator iter_A = m_data.find(*i);
		
		if( iter_A == m_data.end() ){
			throw __FILE__ ":cluster_distance: Unable to find drug A";
		}
		
		for(deque<string>::const_iterator j = m_B.begin();j != m_B.end();++j){
			
			MAP<string, vector<float> >::const_iterator iter_B = m_data.find(*j);

			if( iter_B == m_data.end() ){
				throw __FILE__ ":cluster_distance: Unable to find drug B";
			}

			const float d = drug_distance(iter_A->second, iter_B->second);
			
			if( ( i == m_A.begin() ) && ( j == m_B.begin() ) ){
				min_distance = max_distance = d;
			}
			else{
				min_distance = min(d, min_distance);
				max_distance = max(d, max_distance);
			}
		}
	}
	
	switch(m_alg){
		case CLUSTER_SINGLE_LINKAGE:
			return min_distance;
		case CLUSTER_TOTAL_LINKAGE:
			return max_distance;
		default:
			throw __FILE__ ":cluster_distance: Invalid/unknown clustering algorithm";
	};
	
	throw __FILE__ ":cluster_distance: Unknown clustering algorithm";
	return 0.0;
}

float drug_distance(const vector<float> &m_A, const vector<float> &m_B)
{
	float ret = 0.0;
	
	const size_t N = m_A.size();
	
	if(m_B.size() != N){
		throw __FILE__ ":drug_distance: |A| != |B|";
	}
	
	for(size_t i = 0;i < N;++i){
		
		const float diff = m_A[i] - m_B[i];
		
		ret += diff*diff;
	}
	
	return sqrt(ret);
}
