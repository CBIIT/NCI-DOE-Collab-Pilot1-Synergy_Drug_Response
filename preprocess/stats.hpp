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

#ifndef __STATS
#define __STATS

#include <deque>
#include <algorithm>

template<class T>
T average(typename std::deque<T>::const_iterator m_begin, 
	typename std::deque<T>::const_iterator m_end)
{
	T ret = 0.0;
	size_t count = 0;

	for(typename std::deque<T>::const_iterator i = m_begin;i != m_end;++i){
		
		ret += *i;
		++count;
	}

	if(count == 0){
		throw __FILE__ ":average: count == 0";
	}

	ret /= count;

	return ret;
}

template<class T>
T standard_deviation(typename std::deque<T>::const_iterator m_begin, 
	typename std::deque<T>::const_iterator m_end)
{       
        T ave = 0.0;
	T stdev = 0.0;
        size_t count = 0;

        for(typename std::deque<T>::const_iterator i = m_begin;i != m_end;++i){
                
                ave += *i;
                ++count;
        }

        if(count == 0){
                throw __FILE__ ":average: count == 0";
        }

        ave /= count;

	for(typename std::deque<T>::const_iterator i = m_begin;i != m_end;++i){

		const T local = *i - ave;

                stdev += local*local;
        }

	if(count > 1){
		stdev = sqrt(stdev);
	}
	else{
		stdev = 0.0;
	}

        return stdev;
}

template<class T>
T median(typename std::deque<T>::const_iterator m_begin, 
	typename std::deque<T>::const_iterator m_end)
{
	std::deque<T> local(m_begin, m_end);
	
	std::sort( local.begin(), local.end() );
	
	const size_t len = local.size();
	
	if(len == 0){
		throw __FILE__ ":median: Median not defined (len == 0)!";
	}
	
	if(len%2 == 1){ // Odd
                return local[len/2];
	}
	
	// Even
        return 0.5*(local[len/2 - 1] + local[len/2]);
};

template<class T>
T minimum(typename std::deque<T>::const_iterator m_begin, 
	typename std::deque<T>::const_iterator m_end)
{
	if(m_begin == m_end){
		throw __FILE__ ":minimum: Minimum not defined (len == 0)!";
	}
	
	T ret = *m_begin;
	
	for(typename std::deque<T>::const_iterator i = m_begin;i != m_end;++i){
		ret = min(ret, *i);
	}
	
        return ret;
};

template<class T>
T maximum(typename std::deque<T>::const_iterator m_begin, 
	typename std::deque<T>::const_iterator m_end)
{
	if(m_begin == m_end){
		throw __FILE__ ":maximum: Minimum not defined (len == 0)!";
	}
	
	T ret = *m_begin;
	
	for(typename std::deque<T>::const_iterator i = m_begin;i != m_end;++i){
		ret = max(ret, *i);
	}
	
        return ret;
};

#endif // __STATS
