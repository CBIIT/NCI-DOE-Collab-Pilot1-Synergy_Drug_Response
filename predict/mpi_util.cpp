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

#include "mpi_util.h"
#include "rx.h"
#include <string.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for std::string
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<string>(const string &m_str)
{
	return sizeof(size_t) + m_str.size();
}

template<>
unsigned char* mpi_pack<string>(unsigned char* m_ptr, const string &m_str)
{
	size_t len = m_str.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	memcpy(m_ptr, m_str.c_str(), len);
	m_ptr += len;
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<string>(unsigned char* m_ptr, string &m_str)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_str.assign( (char*)m_ptr, len );
	m_ptr += len;
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for __uint128_t
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<__uint128_t>(const __uint128_t &m_obj)
{
	return sizeof(__uint128_t);
}

template<>
unsigned char* mpi_pack<__uint128_t>(unsigned char* m_ptr, const __uint128_t &m_obj)
{
	memcpy( m_ptr, &m_obj, sizeof(__uint128_t) );
	m_ptr += sizeof(__uint128_t);
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<__uint128_t>(unsigned char* m_ptr, __uint128_t &m_obj)
{
	memcpy( &m_obj, m_ptr, sizeof(__uint128_t) );
	m_ptr += sizeof(__uint128_t);
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for Options
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const Options &m_obj)
{
	return 
	mpi_size(m_obj.print_usage) +
	mpi_size(m_obj.num_tree) +
	mpi_size(m_obj.leaf_size) +
	mpi_size(m_obj.bag_fraction) +
	mpi_size(m_obj.num_fold) +
	mpi_size(m_obj.file_name_drug) +
	mpi_size(m_obj.file_name_blind_drug) +
	mpi_size(m_obj.file_name_drug_target) +
	mpi_size(m_obj.randomize_drug) +
	mpi_size(m_obj.file_name_cell) +
	mpi_size(m_obj.file_name_blind_cell) +
	mpi_size(m_obj.randomize_cell) +
	mpi_size(m_obj.train_cell_to_synergy_file) +
	mpi_size(m_obj.test_cell_to_synergy_file) +
	mpi_size(m_obj.randomize_test_synergy) +
	mpi_size(m_obj.randomize_train_synergy) +
	mpi_size(m_obj.NSC_to_CID_y) +
	mpi_size(m_obj.score_func) +
	mpi_size(m_obj.synergy_threshold) +
	mpi_size(m_obj.file_name_output) +
	mpi_size(m_obj.seed) + 
	mpi_size(m_obj.use_pair_prediction) + 
	mpi_size(m_obj.one_hot_cell) + 
	mpi_size(m_obj.one_hot_drug) + 
	mpi_size(m_obj.one_hot_tissue) + 
	mpi_size(m_obj.normalize_cell) +
	mpi_size(m_obj.overlap_to_train) + 
	mpi_size(m_obj.overlap_to_test);
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, Options &m_obj)
{
	m_ptr = mpi_unpack(m_ptr, m_obj.print_usage);
	m_ptr = mpi_unpack(m_ptr, m_obj.num_tree);
	m_ptr = mpi_unpack(m_ptr, m_obj.leaf_size);
	m_ptr = mpi_unpack(m_ptr, m_obj.bag_fraction);
	m_ptr = mpi_unpack(m_ptr, m_obj.num_fold);
	m_ptr = mpi_unpack(m_ptr, m_obj.file_name_drug);
	m_ptr = mpi_unpack(m_ptr, m_obj.file_name_blind_drug);
	m_ptr = mpi_unpack(m_ptr, m_obj.file_name_drug_target);
	m_ptr = mpi_unpack(m_ptr, m_obj.randomize_drug);
	m_ptr = mpi_unpack(m_ptr, m_obj.file_name_cell);
	m_ptr = mpi_unpack(m_ptr, m_obj.file_name_blind_cell);
	m_ptr = mpi_unpack(m_ptr, m_obj.randomize_cell);
	m_ptr = mpi_unpack(m_ptr, m_obj.train_cell_to_synergy_file);
	m_ptr = mpi_unpack(m_ptr, m_obj.test_cell_to_synergy_file);
	m_ptr = mpi_unpack(m_ptr, m_obj.randomize_test_synergy);
	m_ptr = mpi_unpack(m_ptr, m_obj.randomize_train_synergy);
	m_ptr = mpi_unpack(m_ptr, m_obj.NSC_to_CID_y);
	m_ptr = mpi_unpack(m_ptr, m_obj.score_func);
	m_ptr = mpi_unpack(m_ptr, m_obj.synergy_threshold);
	m_ptr = mpi_unpack(m_ptr, m_obj.file_name_output);
	m_ptr = mpi_unpack(m_ptr, m_obj.seed);
	m_ptr = mpi_unpack(m_ptr, m_obj.use_pair_prediction);	
	m_ptr = mpi_unpack(m_ptr, m_obj.one_hot_cell);
	m_ptr = mpi_unpack(m_ptr, m_obj.one_hot_drug);
	m_ptr = mpi_unpack(m_ptr, m_obj.one_hot_tissue);
	m_ptr = mpi_unpack(m_ptr, m_obj.normalize_cell);
	m_ptr = mpi_unpack(m_ptr, m_obj.overlap_to_train);
	m_ptr = mpi_unpack(m_ptr, m_obj.overlap_to_test);
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const Options &m_obj)
{

	m_ptr = mpi_pack(m_ptr, m_obj.print_usage);
	m_ptr = mpi_pack(m_ptr, m_obj.num_tree);
	m_ptr = mpi_pack(m_ptr, m_obj.leaf_size);
	m_ptr = mpi_pack(m_ptr, m_obj.bag_fraction);
	m_ptr = mpi_pack(m_ptr, m_obj.num_fold);
	m_ptr = mpi_pack(m_ptr, m_obj.file_name_drug);
	m_ptr = mpi_pack(m_ptr, m_obj.file_name_blind_drug);
	m_ptr = mpi_pack(m_ptr, m_obj.file_name_drug_target);
	m_ptr = mpi_pack(m_ptr, m_obj.randomize_drug);
	m_ptr = mpi_pack(m_ptr, m_obj.file_name_cell);
	m_ptr = mpi_pack(m_ptr, m_obj.file_name_blind_cell);
	m_ptr = mpi_pack(m_ptr, m_obj.randomize_cell);
	m_ptr = mpi_pack(m_ptr, m_obj.train_cell_to_synergy_file);
	m_ptr = mpi_pack(m_ptr, m_obj.test_cell_to_synergy_file);
	m_ptr = mpi_pack(m_ptr, m_obj.randomize_test_synergy);
	m_ptr = mpi_pack(m_ptr, m_obj.randomize_train_synergy);
	m_ptr = mpi_pack(m_ptr, m_obj.NSC_to_CID_y);
	m_ptr = mpi_pack(m_ptr, m_obj.score_func);
	m_ptr = mpi_pack(m_ptr, m_obj.synergy_threshold);
	m_ptr = mpi_pack(m_ptr, m_obj.file_name_output);
	m_ptr = mpi_pack(m_ptr, m_obj.seed);
	m_ptr = mpi_pack(m_ptr, m_obj.use_pair_prediction);
	m_ptr = mpi_pack(m_ptr, m_obj.one_hot_cell);
	m_ptr = mpi_pack(m_ptr, m_obj.one_hot_drug);
	m_ptr = mpi_pack(m_ptr, m_obj.one_hot_tissue);
	m_ptr = mpi_pack(m_ptr, m_obj.normalize_cell);
	m_ptr = mpi_pack(m_ptr, m_obj.overlap_to_train);
	m_ptr = mpi_pack(m_ptr, m_obj.overlap_to_test);
	
	return m_ptr;
}
