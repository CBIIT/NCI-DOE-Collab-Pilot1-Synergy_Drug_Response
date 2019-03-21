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

#ifndef __MPI_UTIL
#define	__MPI_UTIL

#include <mpi.h>
#include <limits.h>
#include <deque>
#include <vector>
#include <string>
#include <string.h> // memcpy
#include <unordered_map>

// Use a template for simple objects. Specialize as needed for more complex types
template <class T>
size_t mpi_size(const T &m_obj)
{
	// Force a *compile* time test of whether this is a native or derived type
	static_assert(std::is_fundamental<T>::value || std::is_enum<T>::value,
		":mpi_size: Non-fundamental or non-enum type passed as template");
		
	return sizeof(m_obj);
}

template <class T>
unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj)
{
	// Force a *compile* time test of whether this is a native or derived type
	static_assert(std::is_fundamental<T>::value || std::is_enum<T>::value,
		":mpi_pack: Non-fundamental or non-enum type passed as template");
		
	memcpy( m_ptr, &m_obj, sizeof(m_obj) );
	m_ptr += sizeof(m_obj);
	
	return m_ptr;
}

template <class T>
unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj)
{
	// Force a *compile* time test of whether this is a native or derived type
	static_assert(std::is_fundamental<T>::value || std::is_enum<T>::value,
		":mpi_unpack: Non-fundamental or non-enum type passed as template");
		
	memcpy( &m_obj, m_ptr, sizeof(m_obj) );
	m_ptr += sizeof(m_obj);
	
	return m_ptr;
}

// Specialization for string
template<>
	size_t mpi_size(const std::string &m_str);
template<>
	unsigned char* mpi_pack(unsigned char* m_ptr, const std::string &m_str);
template<>
	unsigned char* mpi_unpack(unsigned char* m_ptr, std::string &m_str);

// Specialization for __uint128_t (which is not recognized by g++ as either a fundamental
// integral or scalar type!). __uint128_t is used to store both accessions and taxonomy ids ...
template<>
	size_t mpi_size(const __uint128_t &m_obj);
template<>
	unsigned char* mpi_pack(unsigned char* m_ptr, const __uint128_t &m_obj);
template<>
	unsigned char* mpi_unpack(unsigned char* m_ptr, __uint128_t &m_obj);
	
/////////////////////////////////////////////////////////////////////////////////////////
// Overload for std::deque
/////////////////////////////////////////////////////////////////////////////////////////
template<class T>
size_t mpi_size(const std::deque<T> &m_obj)
{
	size_t len = sizeof(size_t);
	
	for(typename std::deque<T>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		len += mpi_size(*i);
	}
	
	return len;
}

template<class T>
unsigned char* mpi_pack(unsigned char* m_ptr, const std::deque<T> &m_obj)
{
	size_t len = m_obj.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	for(typename std::deque<T>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		m_ptr = mpi_pack(m_ptr, *i);
	}
	
	return m_ptr;
}

template<class T>
unsigned char* mpi_unpack(unsigned char* m_ptr, std::deque<T> &m_obj)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_obj.resize(len);
	
	for(size_t i = 0;i < len;++i){
		m_ptr = mpi_unpack(m_ptr, m_obj[i]);
	}
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Overload for std::vector
/////////////////////////////////////////////////////////////////////////////////////////
template<class T>
size_t mpi_size(const std::vector<T> &m_obj)
{
	size_t len = sizeof(size_t);
	
	for(typename std::vector<T>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		len += mpi_size(*i);
	}
	
	return len;
}

template<class T>
unsigned char* mpi_pack(unsigned char* m_ptr, const std::vector<T> &m_obj)
{
	size_t len = m_obj.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	for(typename std::vector<T>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		m_ptr = mpi_pack(m_ptr, *i);
	}
	
	return m_ptr;
}

template<class T>
unsigned char* mpi_unpack(unsigned char* m_ptr, std::vector<T> &m_obj)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_obj.resize(len);
	
	for(size_t i = 0;i < len;++i){
		m_ptr = mpi_unpack(m_ptr, m_obj[i]);
	}
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Overload for std::unordered_map
/////////////////////////////////////////////////////////////////////////////////////////
template<class A, class B>
size_t mpi_size(const std::unordered_map<A,B> &m_obj)
{
	size_t len = sizeof(size_t);
	
	for(typename std::unordered_map<A,B>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		
		len += mpi_size(i->first);
		len += mpi_size(i->second);
	}
	
	return len;
}

template<class A, class B>
unsigned char* mpi_pack(unsigned char* m_ptr, const std::unordered_map<A,B> &m_obj)
{
	size_t len = m_obj.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	for(typename std::unordered_map<A,B>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		
		m_ptr = mpi_pack(m_ptr, i->first);
		m_ptr = mpi_pack(m_ptr, i->second);
	}
	
	return m_ptr;
}

template<class A, class B>
unsigned char* mpi_unpack(unsigned char* m_ptr, std::unordered_map<A,B> &m_obj)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_obj.clear();
	
	for(size_t i = 0;i < len;++i){
		
		std::pair<A,B> local;
		
		m_ptr = mpi_unpack(m_ptr, local.first);
		m_ptr = mpi_unpack(m_ptr, local.second);
		
		m_obj.insert(local);
	}
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Overload for std::unordered_multimap
/////////////////////////////////////////////////////////////////////////////////////////
template<class A, class B>
size_t mpi_size(const std::unordered_multimap<A, B> &m_obj)
{
	size_t len = sizeof(size_t);
	
	for(typename std::unordered_multimap<A, B>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		
		len += mpi_size(i->first);
		len += mpi_size(i->second);
	}
	
	return len;
}

template<class A, class B>
unsigned char* mpi_pack(unsigned char* m_ptr, const std::unordered_multimap<A, B> &m_obj)
{
	size_t len = m_obj.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	for(typename std::unordered_multimap<A, B>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		
		m_ptr = mpi_pack(m_ptr, i->first);
		m_ptr = mpi_pack(m_ptr, i->second);
	}
	
	return m_ptr;
}

template<class A, class B>
unsigned char* mpi_unpack(unsigned char* m_ptr, std::unordered_multimap<A, B> &m_obj)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_obj.clear();
	
	for(size_t i = 0;i < len;++i){
		
		std::pair<A, B> local;
		
		m_ptr = mpi_unpack(m_ptr, local.first);
		m_ptr = mpi_unpack(m_ptr, local.second);
		
		m_obj.insert(local);
	}
	
	return m_ptr;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Generic broadcast (from rank zero to all other ranks)
//////////////////////////////////////////////////////////////////////////////////////////////
template <class T>
void broadcast(T &m_obj, const int &m_my_rank, const int &m_src_rank)
{
	size_t len = (m_my_rank == m_src_rank) ? mpi_size(m_obj) : 0;

	if(len >= INT_MAX){
		throw __FILE__ ":broadcast: Max size exceeded";
	}
	
	MPI_Bcast(&len, sizeof(len), MPI_BYTE, m_src_rank, MPI_COMM_WORLD);
	
	unsigned char *buffer = new unsigned char[len];
	
	if(buffer == NULL){
		throw __FILE__ ":broadcast: Unable to allocate buffer";
	}
	
	if(m_my_rank == m_src_rank){
		mpi_pack(buffer, m_obj);		
	}

	MPI_Bcast(buffer, len, MPI_BYTE, m_src_rank, MPI_COMM_WORLD);
		
	if(m_my_rank != m_src_rank){
		mpi_unpack(buffer, m_obj);
	}
		
	delete [] buffer;
}
#endif // __MPI_UTIL
