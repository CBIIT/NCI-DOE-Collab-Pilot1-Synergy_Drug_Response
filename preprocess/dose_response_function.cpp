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

#include "dose_response_function.h"

#include <iostream> // DEBUG
#include <algorithm> // DEBUG

extern "C"
{
	#include <gsl/gsl_multimin.h>
	#include <gsl/gsl_blas.h>
}

using namespace std;

// Parameter indicies for curving fitting with the GSL
#define	LN_ALPHA	0
#define	LN_BETA		1

namespace DoseResponse
{	
	double f_min(const gsl_vector *x, void *params);
	void df_min(const gsl_vector *x, void *params, gsl_vector *df);
	void fdf_min(const gsl_vector *x, void *params, double *f, gsl_vector *df);
}

struct MinimizationParam
{
	double regularization;
	const deque< pair<Concentration, float> > *ptr_xy;
};

// Curve fitting of a sigmoid function to the dose repsonse data is handled by
// the GSL multidimensional minimization functions (as opposed to the non-linear
// curve fitting routines). This is needed to allow for a regularization 
// (i.e. constraint) term in alpha (the slope of the curve at the inflection point).
// If we don't constrain alpha to be "small", noise in the data can lead to large
// values of alpha (that fit the data well, but drop sharply after the highest
// concentration data point) that make it hard to compare functions.
void DoseResponseFunction::fit(const deque< pair<Concentration, float> > &m_xy)
{
	if( m_xy.empty() ){
		throw __FILE__ ":DoseResponseFunction::fit: No data provided!";
	}
	
	const size_t num_param = 2;
	
	gsl_vector *x = gsl_vector_alloc(num_param);
	
	if(x == NULL){
		throw __FILE__ ":DoseResponseFunction::fi: Unable to allocate x vector";
	}
	
	gsl_multimin_fdfminimizer *minimizer = 
		gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, num_param);
	
	if(minimizer == NULL){
		throw __FILE__ ":DoseResponseFunction::fi: Unable to allocate minimizer";
	}
	
	MinimizationParam p;
	
	p.regularization = regularization;
	p.ptr_xy = &m_xy;
	
	gsl_multimin_function_fdf my_func;
	
	my_func.n = num_param;
	my_func.f = DoseResponse::f_min;
	my_func.df = DoseResponse::df_min;
	my_func.fdf = DoseResponse::fdf_min;
	my_func.params = &p;
	
	// Find the range of the *independent* variable
	pair<Concentration, Concentration> range(m_xy.front().first, m_xy.front().first);
	
	for(deque< pair<Concentration, float> >::const_iterator i = m_xy.begin();i != m_xy.end();++i){
		
		range.first = min(range.first, i->first);
		range.second = max(range.second, i->first);
	}
	
	gsl_vector_set( x, LN_ALPHA, 0.0 );
	gsl_vector_set( x, LN_BETA, 
		log( max(1.0f, 
			fabsf( -0.5*(range.second.value + range.first.value) ) 
			) ) );
	
	gsl_multimin_fdfminimizer_set(minimizer, &my_func, x, 0.1, 1e-4);
	
	size_t iteration = 0;
	const size_t max_iteration = 100;
	
	int status = 0;
	
	do
	{
		++iteration;
		status = gsl_multimin_fdfminimizer_iterate(minimizer);

		if(status){
			break;
		}
		
		// A convergence threshold of 1.0e-3 is to large when we need to invert
		// the fitted function to compute the Loewe synergy.
		//status = gsl_multimin_test_gradient(minimizer->gradient, 1e-3);
		status = gsl_multimin_test_gradient(minimizer->gradient, 1e-5);

	}
	while( (status == GSL_CONTINUE) && (iteration < max_iteration) );
	
  	alpha = exp( gsl_vector_get(minimizer->x, LN_ALPHA) );
	beta = exp( gsl_vector_get(minimizer->x, LN_BETA) );
	valid = true; // For now, all drug fits are valid (regardless of convergence)
	
	//////////////////////////////////
	// DEBUG
	
	//bool debug = true;
	
	//for(deque< pair<Concentration, float> >::const_iterator i = m_xy.begin();
	//	i != m_xy.end();++i){
	//	
	//	if(i->second != 1.0){
	//		debug = false;
	//	}
	//}
	//
	//if(debug){
	//
	//	cerr << "alpha = " << alpha << endl;
	//	cerr << "beta = " << beta << endl;
	//	cerr << "iteration = " << iteration << endl;
	//
	//	deque< pair<Concentration, float> > xy = m_xy;
	//
	//	sort( xy.begin(), xy.end() );
	//
	//	for(deque< pair<Concentration, float> >::const_iterator i = xy.begin();
	//		i != xy.end();++i){
	//
	//		cerr << i->first.value << '\t' << i->second << endl;
	//	}
	//
	//	cerr << '\t' << endl;
	//
	//	const float min_c = xy.front().first.value;
	//	const float max_c = xy.back().first.value;
	//	const float delta_c = (max_c - min_c)/20.0;
	//
	//	for(float c = min_c;c <= max_c;c += delta_c){
	//		cerr << c << '\t' << DoseResponse::sigmoid(c, alpha, beta) << endl;
	//	}
	//
	//	cerr << '\t' << endl;
	//}
	
	////////////////////////////////////
	
	gsl_multimin_fdfminimizer_free(minimizer);
	gsl_vector_free(x);
}

double DoseResponse::sigmoid(const double &m_x, const double &m_alpha, const double &m_beta)
{
	//return 100.0*(2.0/( 1.0 + exp( m_alpha*(m_x + m_beta) ) ) - 1.0);
	return (MAX_RESPONSE - MIN_RESPONSE)/( 1.0 + exp( m_alpha*(m_x + m_beta) ) ) + MIN_RESPONSE;
}

double DoseResponse::f_min(const gsl_vector *x, void *params)
{
	double ret = 0.0;
	MinimizationParam* p = (MinimizationParam*)params;
	
	const double alpha = exp( gsl_vector_get(x, LN_ALPHA) );
	const double beta = exp( gsl_vector_get(x, LN_BETA) );
	
	const double norm = p->ptr_xy->empty() ? 1.0 : 1.0/p->ptr_xy->size();
	
	for(deque< pair<Concentration, float> >::const_iterator i = p->ptr_xy->begin();i != p->ptr_xy->end();++i){
		
		const double delta = sigmoid(i->first.value, alpha, beta) - i->second;
		
		ret += delta*delta;
	}
	
	ret *= norm;
	
	ret += p->regularization*alpha*alpha;
	
	return ret;
}

void DoseResponse::df_min(const gsl_vector *x, void *params, gsl_vector *df)
{
	MinimizationParam* p = (MinimizationParam*)params;
	
	const double alpha = exp( gsl_vector_get(x, LN_ALPHA) );
	const double beta = exp( gsl_vector_get(x, LN_BETA) );

	const double norm = p->ptr_xy->empty() ? 1.0 : 1.0/p->ptr_xy->size();
	
	double d_ln_alpha = 0.0;
	double d_ln_beta = 0.0;
	
	for(deque< pair<Concentration, float> >::const_iterator i = p->ptr_xy->begin();i != p->ptr_xy->end();++i){
		
		const double local = exp( alpha*(i->first.value + beta) );
		const double delta = i->second - sigmoid(i->first.value, alpha, beta);
		const double denom = (1.0 + local)*(1.0 + local);
		
		d_ln_alpha += 2.0*delta*(MAX_RESPONSE - MIN_RESPONSE)*local*(i->first.value + beta)/denom;
		
		d_ln_beta += 2.0*delta*(MAX_RESPONSE - MIN_RESPONSE)*local*alpha/denom;
	}
	
	d_ln_alpha *= norm;
	d_ln_beta *= norm;
	
	d_ln_alpha += 2.0*p->regularization*alpha;
	
	// The log correction
	d_ln_alpha *= alpha;
	d_ln_beta *= beta;
	
	gsl_vector_set(df, LN_ALPHA, d_ln_alpha);
	gsl_vector_set(df, LN_BETA, d_ln_beta);
}

void DoseResponse::fdf_min(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{

	*f = DoseResponse::f_min(x, params);
	DoseResponse::df_min(x, params, df);
}
