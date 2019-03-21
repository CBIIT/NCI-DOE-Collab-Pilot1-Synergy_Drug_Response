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

template<class K, class V>
std::vector<K> keys(const MAP<K, V> &m_map)
{
	std::vector<K> ret;
	
	ret.reserve( m_map.size() );
	
	for(typename MAP<K, V>::const_iterator i = m_map.begin();i != m_map.end();++i){
		ret.push_back(i->first);
	}
	
	std::sort( ret.begin(), ret.end() );
	
	return ret;
};

template<class K, class V>
std::vector<K> keys(const MULTIMAP<K, V> &m_map)
{
	std::vector<K> ret;
	
	ret.reserve( m_map.size() );
	
	for(typename MULTIMAP<K, V>::const_iterator i = m_map.begin();i != m_map.end();++i){
		ret.push_back(i->first);
	}
	
	std::sort( ret.begin(), ret.end() );
	ret.erase( unique( ret.begin(), ret.end() ), ret.end() );
	
	return ret;
};

template<class T>
std::vector<T> intersection(std::vector<T> m_A, std::vector<T> m_B) // Make copies of the input vectors
{
	std::vector<T> ret;
	
	// Sort the copies of the inputs to make the set intersection easier to compute
	std::sort( m_A.begin(), m_A.end() );
	std::sort( m_B.begin(), m_B.end() );
	
	typename std::vector<T>::const_iterator a = m_A.begin();
	typename std::vector<T>::const_iterator b = m_B.begin();
	
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
