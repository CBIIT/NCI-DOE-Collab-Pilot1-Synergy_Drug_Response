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

#ifndef __DEQUE_SET
#define __DEQUE_SET

#include <deque>
#include <algorithm>

// A way to potentially speed up the sorting routines (which can consume a non-trivial
// amount of CPU time) is to use the g++-specific __gnu_parallel::sort funcion instead
// of the std::sort function. Since __gnu_parallel::sort is a drop in replacement for
// std::sort, we could include:
//#include <parallel/algorithm>
// and define
//#define	SORT	__gnu_parallel::sort
// or
#define	SORT	std::sort
// to switch between the two versions. This approach still needs benchmarking on the target machine
// (for speed and correctness). This approach will not work on non-g++ compilers (i.e. clang or intel).


template <class T>
void make_set(std::deque<T> &m_set)
{
	SORT( m_set.begin(), m_set.end() );
	m_set.erase( unique( m_set.begin(), m_set.end() ), m_set.end() );
}

template <class T>
std::deque<T> intersection(const std::deque<T> &m_A, const std::deque<T> &m_B)
{
	std::deque<T> ret;
	
	// All inputs must be sorted in ascending order and contain unique elements
	typename std::deque<T>::const_iterator a = m_A.begin();
	typename std::deque<T>::const_iterator b = m_B.begin();
	
	while( ( a != m_A.end() ) && ( b != m_B.end() ) ){
		
		if(*a == *b){
		
			ret.push_back(*a);
			++a;
			++b;
		}
		else{
			if(*a < *b){
				++a;
			}
			else{ // *a > *b
				++b;
			}
		}
	}
	
	return ret;
}

template<class T>
bool set_contains(const std::deque<T> &m_set, const T &m_query)
{
	typename std::deque<T>::const_iterator i = std::lower_bound(m_set.begin(), m_set.end(), m_query);
	
	return ( i != m_set.end() ) && (*i == m_query);
}

template<class T>
size_t get_index(const std::deque<T> &m_set, const T &m_query)
{
	typename std::deque<T>::const_iterator iter = 
		std::lower_bound(m_set.begin(), m_set.end(), m_query);
	
	if( ( iter == m_set.end() ) || (*iter != m_query) ){
		throw __FILE__ ":get_index: Unable to find index!";
	}
	
	return ( iter - m_set.begin() );
};

template<class T>
bool operator==(const std::deque<T> &m_a, const std::deque<T> &m_b)
{
	const size_t len = m_a.size();
	
	if( len != m_b.size() ){
		return false;
	}
	
	for(size_t i = 0;i < len;++i){
		if(m_a[i] != m_b[i]){
			return false;
		}
	}
	
	return true;
}

template<class T>
bool operator!=(const std::deque<T> &m_a, const std::deque<T> &m_b)
{
	const size_t len = m_a.size();
	
	if( len != m_b.size() ){
		return true;
	}
	
	for(size_t i = 0;i < len;++i){
		if(m_a[i] != m_b[i]){
			return true;
		}
	}
	
	return false;
}

template <class T> 
bool find_in_set(const std::deque<T> &m_set, const T &m_query)
{
	typename std::deque<T>::const_iterator iter = 
		lower_bound(m_set.begin(), m_set.end(), m_query);
	
	if( ( iter == m_set.end() ) || (*iter != m_query) ){
		return false;
	}
	
	return true;
}

#endif // __DEQUE_SET
