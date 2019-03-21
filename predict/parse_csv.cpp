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

#include <algorithm>
#include <iostream>
#include <limits.h>
#include <stdlib.h> // atof
#include "parse_csv.h"

using namespace std;

vector<string> split(const string &m_line, const char &m_delim)
{
	// Delimiters are litteral if they are contained within matching protect characters
	const char protect = '"';
	
	vector<string> ret;
	
	size_t len = 1;
	
	size_t protect_count = 0;
	
	for(string::const_iterator i = m_line.begin();i != m_line.end();++i){
		
		// A '\n' or DOS/Windows '\r' symbol forces the end of the line
		if( (*i == '\r') || (*i == '\n') ){
			break;
		}
		
		protect_count += (*i == protect);
		
		len += (*i == m_delim) && (protect_count%2 == 0);
	}
	
	if(protect_count%2 == 1){
		throw __FILE__ ":split: Unmatched protection symbol";
	}
	
	ret.resize(len);
	
	for(size_t i = 0;i < len;++i){
		ret[i].reserve(32);
	}
	
	size_t index = 0;
	
	protect_count = 0;
	
	for(string::const_iterator i = m_line.begin();i != m_line.end();++i){	
		
		// A '\n' or DOS/Windows '\r' symbol forces the end of the line
		if( (*i == '\r') || (*i == '\n') ){
			break;
		}
		
		protect_count += (*i == protect);
		
		if( (*i == m_delim) && (protect_count%2 == 0) ){
			++index;
		}
		else{
			ret[index].push_back(*i);
		}
	}
	
	return ret;
}

size_t get_col(const vector<string> &m_header, const string &m_key)
{
	vector<string>::const_iterator iter = find(m_header.begin(), m_header.end(), m_key);

	if( ( iter == m_header.end() ) || (*iter != m_key) ){
		return COLUMN_NOT_FOUND;
	}
	
	return ( iter - m_header.begin() );
}

size_t string_to_size_t(const string &m_str)
{
	size_t ret = 0;
	
	for(string::const_iterator i = m_str.begin();i != m_str.end();++i){
		
		switch(*i){
			case '0':
				ret = ret*10;
				break;
			case '1':
				ret = ret*10 + 1;
				break;
			case '2':
				ret = ret*10 + 2;
				break;
			case '3':
				ret = ret*10 + 3;
				break;
			case '4':
				ret = ret*10 + 4;
				break;
			case '5':
				ret = ret*10 + 5;
				break;
			case '6':
				ret = ret*10 + 6;
				break;
			case '7':
				ret = ret*10 + 7;
				break;
			case '8':
				ret = ret*10 + 8;
				break;
			case '9':
				ret = ret*10 + 9;
				break;
			default:
				cerr << "Unable to parse: \"" << m_str << "\"" << endl;
				throw __FILE__ ":string_to_size_t: Illegal symbol!";
				break;
		};
	}
	
	return ret;
}

unsigned int string_to_uint(const string &m_str)
{
	const size_t ret = string_to_size_t(m_str);
	
	if(ret > UINT_MAX){
		throw __FILE__ ":string_to_uint: Overflow!";
	}
	
	return (unsigned int)ret;
}

float string_to_float(const string &m_str)
{
	return atof( m_str.c_str() );
}

bool has_digit(const string &m_str)
{
	for(string::const_iterator i = m_str.begin();i != m_str.end();++i){
		if( isdigit(*i) ){
			return true;
		}
	}
	
	return false;
}
