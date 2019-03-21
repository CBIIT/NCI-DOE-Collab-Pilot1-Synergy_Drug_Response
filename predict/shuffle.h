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

#ifndef __RANDOM_SHUFFLE
#define __RANDOM_SHUFFLE

#include <stdlib.h>

// A random_shuffle-like function that uses the re-entrant random number generator
template <class T>
void randomize(const T &m_begin, const T &m_end, unsigned int *m_ptr_seed)
{
	const size_t len = m_end - m_begin;
	
	const double norm = double(len)/RAND_MAX;
	
	for(size_t i = 0;i < len;++i){
	
		// Generate a random number between [0, len)
		size_t index = size_t( rand_r(m_ptr_seed)*norm );
		
		while(index == len){
			index = size_t( rand_r(m_ptr_seed)*norm );
		}
		
		std::swap( *(m_begin + i), *(m_begin + index) );
	}
}

template<class PAIR_CONTAINER>
void randomize_pairs(const PAIR_CONTAINER &m_begin, const PAIR_CONTAINER &m_end, 
	unsigned int *m_ptr_seed)
{
	const size_t len = m_end - m_begin;
	
	const double norm = double(len)/RAND_MAX;
	
	// Randomly shuffle the first element of the pair
	for(size_t i = 0;i < len;++i){
	
		// Generate a random number between [0, len)
		size_t index = size_t( rand_r(m_ptr_seed)*norm );
		
		while(index == len){
			index = size_t( rand_r(m_ptr_seed)*norm );
		}
		
		std::swap( (m_begin + i)->first, (m_begin + index)->first );
	}
	
	// Randomly shuffle the second element of the pair
	for(size_t i = 0;i < len;++i){
	
		// Generate a random number between [0, len)
		size_t index = size_t( rand_r(m_ptr_seed)*norm );
		
		while(index == len){
			index = size_t( rand_r(m_ptr_seed)*norm );
		}
		
		std::swap( (m_begin + i)->second, (m_begin + index)->second );
	}
}

#endif // __RANDOM_SHUFFLE
