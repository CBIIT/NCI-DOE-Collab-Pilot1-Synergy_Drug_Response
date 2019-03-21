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

#ifndef __GEMINI
#define __GEMINI

#include <string>
#include <sstream>
#include <vector>
#include <deque>
#include <unordered_map>

#include <string.h> // memcpy

#include "deque_set.h"
#include "parse_csv.h"
#include "concentration.h"
#include "missing_data.h"
#include "date.h"

#define MAP		std::unordered_map
#define MULTIMAP	std::unordered_multimap

#define	DEFAULT_FOLD			5
#define	DEFAULT_FOREST_SIZE		100
#define	DEFAULT_FOREST_BAG_SAMPLE	1.0
#define	DEFAULT_FOREST_BAG_FEATURE	1.0
#define	DEFAULT_FOREST_LEAF		25

#define	VERSION			"0.06: loewe & bliss"

struct DrugPairData
{
	std::string drug_id_1;
	std::string drug_id_2;
	
	std::vector<float> growth;
};

typedef std::vector<float> DrugFeatures;

typedef enum {
	COMBO_AVERAGE, 
	COMBO_MEDIAN, 
	COMBO_MIN, 
	COMBO_MAX
} ComboScoreType;

typedef enum {
    S_NOT_SET = 0,
    S_FF = (1 << 0), // SRI International
    S_1A = (1 << 1), // National Cancer Institue
    S_FG = (1 << 2), // University of Pittsburg
} Laboratory;

struct GrowthRecord
{
    // COMBODRUGSEQ -- column 1
    std::string record_id; // For debugging
    
    // SCREENER -- column 2
    Laboratory screener;
    
    // STUDY -- column 3
    std::string study;
    
    // TESTDATE -- column 4
    Date testdate;
    
    // PLATE -- column 5
    std::string plate;
    
    // PANELNBR -- column 6
    // Not stored
    
    // CELLNBR -- column 7
    // Not stored
    
    // PREFIX1 -- column 8
    // Not stored
    
    // NSC1 -- column 9
    std::string nsc_1;
    
    // SAMPLE1 -- column 10
    // Not stored
    
    // CONCINDEX1 -- column 11
    int conc_index_1;
    
    // CONC1 -- column 12
    float conc_1;
    
    // CONCUNIT1 -- column 13
    char conc_unit_1;
    
    // PREFIX2 -- column 14
    // Not stored
    
    // NSC2 -- column 15
    std::string nsc_2;
    
    // SAMPLE2 -- column 16
    // Not stored
    
    // CONCINDEX2 -- column 17
    int conc_index_2;
    
    // CONC2 -- column 18
    float conc_2;
    
    // CONCUNIT2 -- column 19
    char conc_unit_2;
    
    // PERCENTGROWTH -- column 20
    float __percent_growth; // Stored to compare with computed value
    
    // PERCENTGROWTHNOTZ -- column 21
    // Not stored
    
    // TESTVALUE -- column 22
    float test_value;
    
    // CONTROLVALUE -- column 23
    float control_value;
    
    // TZVALUE -- column 24
    float tz_value;
    
    // EXPECTEDGROWTH -- column 25
    // Not stored
    
    // SCORE -- column 26
    // Not stored
    
    // VALID -- column 27
    // Not stored
    
    // PANEL -- column 28
    // Not stored
    
    // CELLNAME -- column 29
    std::string cellname;
    
    // Did this record pass the QC process defined by the Almanac paper?
    bool passed_almanac_qc;
    
    GrowthRecord()
    {
        screener = S_NOT_SET;
        conc_index_1 = -1;
        conc_1 = 0.0;
        conc_unit_1 = '?';
        conc_index_2 = -1;
        conc_2 = 0.0;
        conc_unit_2 = '?';
	__percent_growth = 0.0;
        test_value = 0.0;
        control_value = 0.0;
        tz_value = 0.0;
        passed_almanac_qc = true;
    };
    
    void parse_screener(const std::string &m_str)
    {
        if(m_str == "FF"){
            screener = S_FF;
        }
        else{
            if(m_str == "1A"){
                screener = S_1A;
            }
            else{
                if(m_str == "FG"){
                    screener = S_FG;
                }
                else{
                    throw __FILE__ ":parse_screener: Unknown screener!";
                }
            }
        }
    };
    
    inline bool is_control() const
    {
        if(conc_index_1 == -1){
            return true;
        }
        
        return false;
    };
    
    inline bool is_single() const
    {
        return nsc_2.empty();
    };
    
    inline bool is_pair() const
    {
        return (nsc_2.empty() == false);
    };
    
    inline size_t num_drug() const
    {
        return !nsc_1.empty() + !nsc_2.empty();
    };
    
    inline float percent_growth() const
    {
        if(test_value < tz_value){
            return 100.0f*(test_value - tz_value)/tz_value;
        }
        
        return 100.0f*(test_value - tz_value)/(control_value - tz_value);
    };
    
    inline Concentration concentration_1() const
    {
        return Concentration(conc_1, conc_unit_1);
    };
    
    inline Concentration concentration_2() const
    {
        return Concentration(conc_2, conc_unit_2);
    };
    
    void swap_drugs()
    {
    	std::swap(nsc_1, nsc_2);
	std::swap(conc_index_1, conc_index_2);
	std::swap(conc_1, conc_2);
	std::swap(conc_unit_1, conc_unit_2);
    };
    
    // This is used for merging replicate records
    inline GrowthRecord& operator+=(const GrowthRecord &m_rhs)
    {
    	__percent_growth += m_rhs.__percent_growth;
	test_value += m_rhs.test_value;
	control_value += m_rhs.control_value;
	tz_value += m_rhs.tz_value;
	
	return *this;
    };
    
    // For sorting
    inline bool operator<(const GrowthRecord &m_rhs) const
    {
    	// Sort growth recording using the following:
	if(screener != m_rhs.screener){
		return screener < m_rhs.screener;
	}
	
	if(cellname != m_rhs.cellname){
		return cellname < m_rhs.cellname;
	}
	
	if(testdate != m_rhs.testdate){
		return testdate < m_rhs.testdate;
	}
	
	// There can be multiple (i.e. replicate) studies with the same cellname,
	// drug 1 and drug 2 names on the *same* date!
	if(study != m_rhs.study){
		return study < m_rhs.study;
	}
	
	// Drug pair records are *less* than single drug records
	if( is_pair() ){
		
		if( m_rhs.is_pair() ){
			
			// There can be multiple (i.e. replicate) studies with the same cellname,
			// drug 1 and drug 2 names on the *same* date and *same* study!
			// However, this rule does not apply to screening site S_1A, which has
			// many experiments split across multiple plates ...
			if( (screener != S_1A) && (plate != m_rhs.plate) ){
				return plate < m_rhs.plate;
			}

			if(nsc_1 != m_rhs.nsc_1){
				return nsc_1 < m_rhs.nsc_1;
			}
	
			return nsc_2 < m_rhs.nsc_2;
		}
		
		// m_rhs is a single drug record
		return true;
	}
	
	// This is a single drug record
	if( m_rhs.is_single() ){
		return nsc_1 < m_rhs.nsc_1;
	}
	
	return false;
    };
};

struct ExperimentRecord
{
	// An "experiment" is a collection of GrowthRecords that correspond
	// to a particular drug pair/cell line/study	
	std::deque<GrowthRecord> pair_data;
	std::deque<GrowthRecord> drug_1_data;
	std::deque<GrowthRecord> drug_2_data;
};

template <class K, class V>
std::deque<K> keys(const MAP<K, V> &m_data)
{
	std::deque<K> ret;
	
	for(typename MAP<K, V>::const_iterator i = m_data.begin();i != m_data.end();++i){
		ret.push_back(i->first);
	}
	
	make_set(ret);
	
	return ret;
}

template <class K, class V>
std::deque<K> keys(const MULTIMAP<K, V> &m_data)
{
	std::deque<K> ret;
	
	for(typename MULTIMAP<K, V>::const_iterator i = m_data.begin();i != m_data.end();++i){
		ret.push_back(i->first);
	}
	
	make_set(ret);
	
	return ret;
}

// parse_combo.cpp
bool parse_combo(const std::string &m_filename, 
    std::deque<GrowthRecord> &m_data,
    const int &m_lab_mask);

// io.cpp
std::deque<DrugPairData> parse_pair_drug(std::string &m_filename,
                                           std::vector<std::string> &m_cell_names,
                                           const float &m_normalize_by = 1.0);

void parse_single_drug(MULTIMAP<std::string, DrugFeatures> &m_single,
	const std::string &m_filename, 
    	const float &m_normalize_by = 1.0);

void parse_fingerprints(MAP<std::string, DrugFeatures> &m_single, 
	const std::string &m_filename);

// nearest_neighbor.cpp
std::vector<float> nearest_neighbor_growth_prediction(
	const std::deque< std::pair< std::vector<float>, std::vector<float> > > &m_training,
	const std::vector<float> &m_x,
	const size_t &m_num_nearest_neighbors);
							
#endif // __GEMINI
