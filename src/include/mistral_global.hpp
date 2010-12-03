
/*
  Mistral is a constraint satisfaction and optimisation library
  Copyright (C) 2003-2005  Emmanuel Hebrard
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  The author can be contacted electronically at emmanuel.hebrard@gmail.com.
*/


/** \file mistral_global.hpp
    \brief Global Definitions.
*/


#ifndef __GLOBAL_HPP
#define __GLOBAL_HPP


#include <string>
#include <iostream>


namespace Mistral {



  typedef int Event;
  typedef int Outcome;
  

  
  const int getlast[256] = {-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7};
  
  const int MIN_CAPACITY   = 16;  
  const int NOVAL          = (int)((~(unsigned int)0)/2);
  const int MAXINT         =  NOVAL;
  const int MININT         = -NOVAL;

//   //#ifdef _STATIC_CAST

// #define BITSET_VAR  0xe0000000;
// #define LIST_VAR    0xc0000000;
// #define RANGE_VAR   0xa0000000;
// #define VIRTUAL_VAR 0x60000000;
// #define CONST_VAR   0x40000000;

//   //#else

//   const int NOTYPE         = 0;
//   const int BITSET         = 1;
//   const int RANGE          = 2;
//   const int LIST           = 4;
//   const int VIRTUAL_REV    = 8;
//   const int VIRTUAL_DOM    = 16;
//   const int CONSTANT       = 32;
  
//   const int BITSET_VAR     = (BITSET  | RANGE | VIRTUAL_REV);
//   const int BOOL_VAR       = (BITSET  | RANGE);
//   const int RANGE_VAR      = (RANGE   | VIRTUAL_REV);
//   const int LIST_VAR       = (BITSET  | RANGE | VIRTUAL_REV | VIRTUAL_DOM | LIST);
//   const int VIRTUAL_VAR    = (VIRTUAL_REV | VIRTUAL_DOM | RANGE);
//   const int CONST_VAR      = (BITSET  | RANGE | CONSTANT);
//   const int DYN_VAR        = 1111111

  const int CONST_VAR      = 1;
  const int BOOL_VAR       = 2;
  const int RANGE_VAR      = 4;
  const int BITSET_VAR     = 8;
  const int LIST_VAR       = 16;
  const int VIRTUAL_VAR    = 0;
  const int DYN_VAR        = (CONST_VAR | BOOL_VAR | RANGE_VAR | BITSET_VAR | LIST_VAR);
  const int EXPRESSION     = 3;


    //#endif

  //       fnvrlu
  //fail   000000
  //no     010000
  //domain 100000
  //range  100100
  //lb     100110
  //ub     100101
  // value 101100

  const Event NO_EVENT     = 0;
  const Event DOMAIN_EVENT = 1;
  const Event RANGE_EVENT  = 1+2;
  const Event UB_EVENT     = 1+2+4;
  const Event LB_EVENT     = 1+2+8;
  const Event VALUE_EVENT  = 1+2+16;
  const Event FAIL_EVENT   = 32;

  const Outcome SAT = 1;
  const Outcome OPT = 3;
  const Outcome UNSAT = 0;
  const Outcome UNKNOWN = 2;
  const Outcome LIMITOUT = 4;

//   const Event NO_EVENT     = 32;
//   const Event DOMAIN_EVENT = 1;
//   const Event RANGE_EVENT  = 3;
//   const Event UB_EVENT     = 7;
//   const Event LB_EVENT     = 11;
//   const Event VALUE_EVENT  = 19;
//   const Event FAIL_EVENT   = 0;

  const int _value_  = 0;
  const int _range_  = 1;
  const int _domain_ = 2;

  inline bool is_domain(Event e) {
    return ((e & DOMAIN_EVENT) == DOMAIN_EVENT);
  }

  inline bool is_range(Event e) {
    return ((e & RANGE_EVENT) == RANGE_EVENT);
  }

  inline bool is_upper_bound(Event e) {
    return ((e & UB_EVENT) == UB_EVENT);
  }

  inline bool is_lower_bound(Event e) {
    return ((e & LB_EVENT) == LB_EVENT);
  }

  inline bool is_value(Event e) {
    return ((e & VALUE_EVENT) == VALUE_EVENT);
  }

  /**********************************************
   * Timing Memory and Command line utilities 
   *********************************************/

  double getRunTime();
  unsigned long int getMemory();
  void getCommandLine(const char**,int*,int,const char**,const char**,int,char**,int);


//   std::string toString(const int x);


//   std::string toString(const IntStack& x);


//   std::string toString(const Queue& x);


//   std::string toString(const MultiSet& x);


//   std::string toString(const ConstraintTrigger& x);


//   std::string toString(const Constraint* x);


//   std::string toString(const SolverStatistics& x);

// //   class VariableInt;
// //   std::string toString(const VariableInt* x);

//   //class Goal;
//   //std::string toString(const Goal* x);


  class IntStack;
  std::ostream& operator<< (std::ostream& os, const IntStack& x);
  
  class Queue;
  std::ostream& operator<< (std::ostream& os, const Queue& x);

  class MultiSet;
  std::ostream& operator<< (std::ostream& os, const MultiSet& x);

  class ConstraintTrigger;
  std::ostream& operator<< (std::ostream& os, const ConstraintTrigger& x);

  class Constraint;
  std::ostream& operator<< (std::ostream& os, const Constraint& x);

  class Variable;
  std::ostream& operator<< (std::ostream& os, const Variable& x);

  class BitsetDomain;
  std::ostream& operator<< (std::ostream& os, const BitsetDomain& x);

  class SolverStatistics;
  std::ostream& operator<< (std::ostream& os, const SolverStatistics& x);


  std::ostream& operator<< (std::ostream& os, const IntStack* x);

  std::ostream& operator<< (std::ostream& os, const Queue* x);

  std::ostream& operator<< (std::ostream& os, const MultiSet* x);

  std::ostream& operator<< (std::ostream& os, const ConstraintTrigger* x);

  std::ostream& operator<< (std::ostream& os, const Constraint* x);

  std::ostream& operator<< (std::ostream& os, const Variable* x);

  std::ostream& operator<< (std::ostream& os, const BitsetDomain* x);

  std::ostream& operator<< (std::ostream& os, const SolverStatistics* x);

  //std::ostream& operator<< (std::ostream& os, const VariableInt* x);

  //std::ostream& operator<< (std::ostream& os, const Goal* x);


  template <class WORD_TYPE>
  void showUint(WORD_TYPE n, std::ostream& os) {
    WORD_TYPE mask=1;
    while(mask){
      if(mask & n) os<<1;
      else os<<0;
      mask = mask << 1;
    }
  }
  template <class WORD_TYPE>
  void printBitset(WORD_TYPE n, const int idx, std::ostream& os) {

//     os << std::endl;
//     showUint(n, os);
//     os << std::endl;

    int offset = 8*sizeof(WORD_TYPE)*idx;
    WORD_TYPE mask=1;
    int last, cur, serie=0, k=0;
    os << "{";
    while(mask){
      if(mask & n) {
	last = (offset + k);
	os << last;
	break;
      }
      mask = mask << 1;
      ++k;
    }
    
    ++k;
    mask <<= 1;
    while(mask){
      if(mask & n) {
	cur = (k + offset);
	if(cur == last+1) ++serie;
	else {
	  if(serie > 1) os << ".." << last ;
	  else if(serie > 0) os << "," << last;
	  os << "," << cur;
	  serie = 0;
	}
	last = cur;
      }
      ++k;
      mask = mask << 1;
    }
    if(serie > 1) os << ".." << last ;
    else if(serie > 0) os << "," << last;
    
    os << "}";
  }

  /**********************************************
   * Knuth's Random number generator (code from sp-1.4)
   **********************************************/

  void usrand (unsigned seed);
  unsigned urand0 (void);
  unsigned urand (void);
  int randint(int upto);
  double randreal();


  /// Other utils:

  int log2( const unsigned int v );

}

#define MAX_URAND 0xFFFFFFFFL

   
#endif // __GLOBAL_HPP


/**
   Boolean: 1 bit for the domain, 31 bits for the assignment level
   Range: 32 bits for the min, 32 bits for the max
   Word: 5 bits for the min, 5 bits for the max, 5 bits for the size, 32 bits for the domain, 17 bits for the offset
   Integer: 32 bits for the min, 32 bits for the max, 32 bits for a ptr to the rest of the domain 

00
01
10
11
*/

/**
class Intvar {
public:
  char type;
  void *_ptr;





  int min(, Variable *X) {
    if(type == BOOLEAN) {
      int x = (int)(*X);
      return (x&1 && x>>1 <= level); 
    } else {
      VariableFDomain *x = (VariableFDomain*)X;
      return X->min;
    }
}
  
}

-- Variables are read only structs.
-- WriteManager class takes care of the reversibility.
*/




