
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


#ifndef _MISTRAL_GLOBAL_HPP
#define _MISTRAL_GLOBAL_HPP


#include <string>
#include <iostream>

#include <mistral_structure.hpp>


namespace Mistral {



  typedef int Event;
  typedef int Outcome;
  typedef int PropagationOutcome;


  typedef unsigned int Lit;
  typedef unsigned int Atom;
  typedef unsigned int Value;

  typedef Array<Lit> Clause;
  
  
#define NORESTART 0
#define GEOMETRIC 1
#define LUBY 2

#define LARGE_VALUE NOVAL/16384
#define INFTY  NOVAL/2
#define MAXINT NOVAL
#define MININT -NOVAL
#define MIN_CAPACITY 16
  
#define BOUND_CONSISTENCY 1
#define FORWARD_CHECKING 0

#define CONST_VAR   1
#define BOOL_VAR    2
#define RANGE_VAR   4
#define BITSET_VAR  8
#define LIST_VAR    16
#define VIRTUAL_VAR 0
#define DYN_VAR     27
  //#define DYN_VAR     31
#define REMOVED_VAR 512
#define EXPRESSION  3
  
#define DOMAIN_C     1
#define RANGE_C      2
#define UB_C         4
#define LB_C         8
#define VALUE_C     16

#define NO_EVENT     0
#define DOMAIN_EVENT 1
  // DOMAIN_C+RANGE_C
#define RANGE_EVENT  3
  // DOMAIN_C+RANGE_C+UB_C
#define UB_EVENT     7
  // DOMAIN_C+RANGE_C+LB_C
#define LB_EVENT     11
  // DOMAIN_C+RANGE_C+VALUE_C
#define VALUE_EVENT  31
#define FAIL_EVENT   32

#define DOMAIN_CHANGED(e) (bool)((e)&DOMAIN_C)
#define RANGE_CHANGED(e)  (bool)((e)&RANGE_C)
#define LB_CHANGED(e)     (bool)((e)&LB_C)
#define UB_CHANGED(e)     (bool)((e)&UB_C)
#define ASSIGNED(e)       (bool)((e)&VALUE_C)

#define IS_FAIL(e) ((e)&FAIL_EVENT)

#define SAT      1
#define OPT      3
#define UNSAT    0
#define UNKNOWN  2
#define LIMITOUT 4

#define IS_OK(o) (o<0)
#define FAILURE(x) x
#define CONSISTENT -1
  
#define _VALUE_ 0
#define _RANGE_ 1
#define _DOMAIN_ 2

#define EVENT_TYPE(e) (2-(RANGE_CHANGED(e))-(ASSIGNED(e)))


  /**********************************************
   * Timing Memory and Command line utilities 
   *********************************************/

  double getRunTime();
  unsigned long int getMemory();
  void getCommandLine(const char**,int*,int,const char**,const char**,int,char**,int);

  template <class WORD_TYPE>
  void printBitset(WORD_TYPE n, const int idx, std::ostream& os) {
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

  // class Random {

  // public:
  //   unsigned int mistral_rand_x[56];
  //   unsigned int mistral_rand_y[256];
  //   unsigned int mistral_rand_z;

  //   int mistral_rand_j;
  //   int mistral_rand_k;

  //   Random(unsigned int seed=12345) { usrand(seed); }
  //   virtual ~Random() {}

    void usrand (unsigned seed);
    unsigned urand0 (void);
    unsigned urand (void);
    int randint(int upto);
    double randreal();
  //  };

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




