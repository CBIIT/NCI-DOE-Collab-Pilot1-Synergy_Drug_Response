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

#ifndef __DOSE_RESPONSE_FUNCTION
#define __DOSE_RESPONSE_FUNCTION

#include <deque>
#include <string.h> // Needed for memcpy
#include "concentration.h"

//#define	MAX_RESPONSE	100.0
//#define	MIN_RESPONSE	-100.0

#define	MAX_RESPONSE	1.0
#define	MIN_RESPONSE	0.0

#define	DEFAULT_DOSE_RESPONSE_REGULARIZATION	10.0

// Currently a two parameter sigmoid function: 
// 	100*(2/( 1 + exp( alpha*(x + beta) ) ) - 1)
namespace DoseResponse{
	double sigmoid(const double &m_x, const double &m_alpha, const double &m_beta);
}

class DoseResponseFunction
{
	private:
	
		double alpha;
		double beta;
		double regularization;
		
		// Was the curve fitting sucessfull?
		bool valid;		
	public:

		DoseResponseFunction(const double &m_regularization)
		{
			alpha = beta = 0.0;
			valid = false;
			regularization = m_regularization;
		};
		
		inline operator bool() const
		{
			return valid;
		};
		
		inline double get_alpha() const
		{
			return alpha;
		};
		
		inline void set_alpha(const double &m_alpha)
		{
			alpha = m_alpha;
		};
		
		inline double get_beta() const
		{
			return beta;
		};
		
		inline void set_beta(const double &m_beta)
		{
			beta = m_beta;
		};
		
		inline double operator()(const Concentration &m_x) const
		{
			return DoseResponse::sigmoid(m_x.value, alpha, beta);
		};
		
		void fit(const std::deque< std::pair<Concentration, float> > &m_xy);
		
		inline size_t mpi_size() const
		{
			return sizeof(alpha) + sizeof(beta) + 
				sizeof(valid);
		};
		
		inline unsigned char* mpi_pack(unsigned char* m_ptr) const
		{
			memcpy( m_ptr, &alpha, sizeof(alpha) );
			m_ptr += sizeof(alpha);
			
			memcpy( m_ptr, &beta, sizeof(beta) );
			m_ptr += sizeof(beta);
			
			memcpy( m_ptr, &valid, sizeof(valid) );
			m_ptr += sizeof(valid);
			
			return m_ptr;
		};
		
		inline unsigned char* mpi_unpack(unsigned char* m_ptr)
		{
			memcpy( &alpha, m_ptr, sizeof(alpha) );
			m_ptr += sizeof(alpha);
			
			memcpy( &beta, m_ptr, sizeof(beta) );
			m_ptr += sizeof(beta);
			
			memcpy( &valid, m_ptr, sizeof(valid) );
			m_ptr += sizeof(valid);
			
			return m_ptr;
		};
};

#endif // __DOSE_RESPONSE_FUNCTION
