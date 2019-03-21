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

#ifndef __DATE
#define __DATE

#include <string>

struct Date 
{    

    // Start the months at 1
    enum {JAN = 1, FEB, MAR, APR, MAY, JUN, JUL,
        AUG, SEP, OCT, NOV, DEC, NO_MONTH} month;
    unsigned int day;
    unsigned int year;
    
    Date()
    {
        
    };
    
    Date(const std::string &m_str)
    {
    	// Date strings can be in two different formats!
	// The original format looks like: Mon Jul 18 00:00:00 2011
	// The new format looks like: 06/23/2011
        
	// Use the following heuristic: If the input string only contains numbers and '/',
	// then it belongs to the new format. Otherwise, treat it as the old format.
	bool use_new_format = true;
	
	for(std::string::const_iterator i = m_str.begin();i != m_str.end();++i){
		
		if( isalpha(*i) ){
			
			// This belongs to the old format!
			use_new_format = false;
			break;
		}
	}
	
	if(use_new_format){
		parse_month_day_year(m_str);
	}
	else{
		parse_day_month_day_time_year(m_str);
	}
	
    };
    
    inline bool operator==(const Date &m_rhs) const
    {
        return (year == m_rhs.year) &&
            (month == m_rhs.month) &&
            (day == m_rhs.day);
    };
    
    inline bool operator!=(const Date &m_rhs) const
    {
    	return (year != m_rhs.year) ||
            (month != m_rhs.month) ||
            (day != m_rhs.day);
    };
    
    inline bool operator<(const Date &m_rhs) const
    {
    	if(year != m_rhs.year){
		return year < m_rhs.year;
	}
	
	if(month != m_rhs.month){
		return month < m_rhs.month;
	}
	
	return day < m_rhs.day;
    };
    
    // Compute the *approximate* difference between two dates (in days)
    inline int operator-(const Date &m_rhs) const
    {
    	int ret = 0;
	
	ret += 365*( int(year) - int(m_rhs.year) );
	ret += 30*( int(month) - int(m_rhs.month) );
	ret += int(day) - int(m_rhs.day);
	
	return ret;
    };
    
    inline std::string str() const
    {
    	std::stringstream ssout;
	
	switch(month){
		case JAN:
			ssout << "01";
			break;
		case FEB:
			ssout << "02";
			break;
		case MAR:
			ssout << "03";
			break;
		case APR:
			ssout << "04";
			break;
		case MAY:
			ssout << "05";
			break;
		case JUN:
			ssout << "06";
			break;
		case JUL:
			ssout << "07";
			break;
		case AUG:
			ssout << "08";
			break;
		case SEP:
			ssout << "09";
			break;
		case OCT:
			ssout << "10";
			break;
		case NOV:
			ssout << "11";
			break;
		case DEC:
			ssout << "12";
			break;
		default:
			throw __FILE__ ":Date::str: Unknown month!";
	};
	
	ssout << '/' << day << '/' << year;
	
	return ssout.str();
    };
    
private:
    
    void parse_month_day_year(const std::string &m_str)
    {
    	const size_t len = m_str.size();
        size_t start = 0;
	
	// Skip any leading white space
        while( (start < len) && isspace(m_str[start]) ){
            ++start;
        }
	
	size_t stop = start;
	
	// Read the month
	while( (stop < len) && isdigit(m_str[stop]) ){
		++stop;
	}
	
	unsigned int tmp = string_to_uint( m_str.substr(start, stop - start) );
	
	switch(tmp){
		case 1:
			month = JAN;
			break;
		case 2:
			month = FEB;
			break;
		case 3:
			month = MAR;
			break;
		case 4:
			month = APR;
			break;
		case 5:
			month = MAY;
			break;
		case 6:
			month = JUN;
			break;
		case 7:
			month = JUL;
			break;
		case 8:
			month = AUG;
			break;
		case 9:
			month = SEP;
			break;
		case 10:
			month = OCT;
			break;
		case 11:
			month = NOV;
			break;
		case 12:
			month = DEC;
			break;
		default:
			throw __FILE__ ":parse_month_day_year: Invalid month!";
			break;
	};
	
	// Skip the delimiter
	if( (stop == len) || (m_str[stop] != '/') ){
		throw __FILE__ ":parse_month_day_year: Failed to read month-day delimiter";
	}
	
	++stop;
	start = stop;
	
	// Read the day
	while( (stop < len) && isdigit(m_str[stop]) ){
		++stop;
	}
	
	parse_day( m_str.substr(start, stop - start) );
		
	// Skip the delimiter
	if( (stop == len) || (m_str[stop] != '/') ){
		throw __FILE__ ":parse_month_day_year: Failed to read day-year delimiter";
	}
	
	++stop;
	start = stop;
	
	// Read the year
	while( (stop < len) && isdigit(m_str[stop]) ){
		++stop;
	}
	
	parse_year( m_str.substr(start, stop - start) );
    };
    
    void parse_day_month_day_time_year(const std::string &m_str)
    {
        // We expect a space delimited date string
        const size_t len = m_str.size();
        size_t start = 0;
        
        // Skip leading white space
        while( (start < len) && isspace(m_str[start]) ){
            ++start;
        }
        
        // Skip the day of the week
        while( (start < len) && !isspace(m_str[start]) ){
            ++start;
        }
        
        // Skip the white space between day of the week and month
        while( (start < len) && isspace(m_str[start]) ){
            ++start;
        }
        
        size_t stop = start;
        
        // Extract the month
        while( (stop < len) && !isspace(m_str[stop]) ){
            ++stop;
        }
        
        parse_month( m_str.substr(start, stop - start) );
        
        start = stop;
        
        // Skip the white space between month and day
        while( (start < len) && isspace(m_str[start]) ){
            ++start;
        }
        
        stop = start;
        
        // Extract the day
        while( (stop < len) && !isspace(m_str[stop]) ){
            ++stop;
        }
        
        parse_day( m_str.substr(start, stop - start) );
        
        start = stop;
        
        // Skip the white space between day and time
        while( (start < len) && isspace(m_str[start]) ){
            ++start;
        }
        
        // Skip the time
        while( (start < len) && !isspace(m_str[start]) ){
            ++start;
        }
        
        // Skip the white space between time and year
        while( (start < len) && isspace(m_str[start]) ){
            ++start;
        }
        
        stop = start;
        
        // Extract the year
        while( (stop < len) && !isspace(m_str[stop]) ){
            ++stop;
        }
        
        parse_year( m_str.substr(start, stop - start) );
    };
    
    void parse_day(const std::string &m_str)
    {
        day = string_to_uint(m_str);
        
        if( (day < 1) || (day > 31) ){
            throw __FILE__ ":parse_day: Day out of bounds";
        }
    };
    
    void parse_year(const std::string &m_str)
    {
        year = string_to_uint(m_str);
        
        // There are some growth measurements from 2009 (!)
        if( (year < 2009) || (year > 2017) ){
            throw __FILE__ ":parse_year: Year out of bounds";
        }
    };
    
    void parse_month(std::string m_str)
    {
        // Make the month string all upper case
        for(std::string::iterator i = m_str.begin();i != m_str.end();++i){
            *i = toupper(*i);
        }
        
        month = NO_MONTH;
        
        if(m_str == "JAN"){
            month = JAN;
        }
        else{
            if(m_str == "FEB"){
                month = FEB;
            }
            else{
                if(m_str == "MAR"){
                    month = MAR;
                }
                else{
                    if(m_str == "APR"){
                        month = APR;
                    }
                    else{
                        if(m_str == "MAY"){
                            month = MAY;
                        }
                        else{
                            if(m_str == "JUN"){
                                month = JUN;
                            }
                            else{
                                if(m_str == "JUL"){
                                    month = JUL;
                                }
                                else{
                                    if(m_str == "AUG"){
                                        month = AUG;
                                    }
                                    else{
                                        if(m_str == "SEP"){
                                            month = SEP;
                                        }
                                        else{
                                            if(m_str == "OCT"){
                                                month = OCT;
                                            }
                                            else{
                                                if(m_str == "NOV"){
                                                    month = NOV;
                                                }
                                                else{
                                                    if(m_str == "DEC"){
                                                        month = DEC;
                                                    }
                                                    else{
                                                        throw __FILE__ ":parse_month: Unknown month";
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        if(month == NO_MONTH){
            throw __FILE__ ":parse_month: Did not find a valid month!";
        }
    };
};

namespace std
{   
    template<> struct hash<Date>
    {
        size_t operator()(const Date& s) const
        {
            size_t ret = 0;
            
	    // Max allowed year is 2048 (11 bits)
	    ret = (ret << 11) | s.year;
	    
	    // Max allowed month is 16 (4 bits)
	    ret = (ret << 4) | (unsigned int)(s.month);
	    
            // Max allowed day is 32 (5 bits)
	    ret = (ret << 5) | s.day;
	    
	    return ret;
        };
    };
}


#endif // __DATE
