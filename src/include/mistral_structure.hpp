
/*
  Mistral 2.0 is a constraint satisfaction and optimisation library
  Copyright (C) 2009  Emmanuel Hebrard
  
	This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.

  The author can be contacted electronically at emmanuel.hebrard@gmail.com.
*/


/*! \file mistral_structure.hpp
  \brief Header file for the basic data-structures.
*/


#include <iostream>
#include <iomanip>
#include <stdlib.h> 
#include <sstream> 
#include <string.h>
#include <limits.h>
#include <assert.h>
//#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

//#include <mistral_global.hpp>

//Temporary
#include <vector>
#include <algorithm>

#ifndef __STRUCTURE_HPP
#define __STRUCTURE_HPP


namespace Mistral {


template <class WORD_TYPE>
void showUint(WORD_TYPE n, std::ostream& os) {
  WORD_TYPE mask=1;
  while(mask){
    if(mask & n) os<<1;
    else os<<0;
    mask = mask << 1;
  }
}

const int getlast[256] = {-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7};

const int NOVAL = (int)((~(unsigned int)0)/2);
#define INFTY  NOVAL/2

#define IN_COVER(x)  (2*x)
#define IN_STABLE(x) (2*x+1)
#define NODE(e) (e/2)
#define STATUS(e) (e&1)
#define COV 0
#define STA 1


#define _VALGRIND_


//#define _STORE_NOGOOD
#define _FORGET_NOGOOD
//#define _TRACE_ true
//#define _DEBUG_BACKJUMPS true //true //(nb_nodes >= 48430)
//#define _DEBUG_REDUCTION true
#define _USE_CLIQUES true
#define _COVER_SET true
//#define _NOGOOD_STAT true
//#define _VERIF_TRAIL true
//#define _VERIF_REASON true
//#define _VERIF_WATCHERS true
//#define _DEBUG_UPDATE
//#define _DEBUG_CLIQUEDOM true //(node.size<=4)
//#define _DEBUG_ADJCC true //(V.size<=4)
//#define _COUNT_OP			
//#define _DEBUG_CLIQUECOV
//#define _DEBUG_NEIGHCLIQUE true
//#define _DEBUG_KERNCLIQUE true
//#define _DEBUG_KERNCLIQUEINCR true
//#define _DEBUG_BUSS true
//#define _DEBUG_UNWATCH true
//#define _VERIFIED_RCG true
//#define _DEBUG_UPDATE true
//#define _UNOGOOD_


#define _DEGREE_
// #define _DEG_SAT_
// #define _SAT_DEG_
// #define _DEG_PLUS_SAT_

#define _NOT_SO_COMPACT



// WARNING: Maintaining the dual would make conflict extraction and dominance detection wrong! (we need the orginal degree/edges)
//#define	_MAINTAIN_DUAL		


#ifdef _STORE_NOGOOD				
#define NO_REASON NULL
#else
#define NO_REASON -1
#define DEDUCTION -2
#endif	






//   /**********************************************
//    * Vector
//    **********************************************/
//   /*! \class Vector
//     \brief Simple vector class     
//   */
  
//   template <class DATA_TYPE>
//   class Vector {
//   public:

//     /*!@name Parameters*/
//     //@{
//     DATA_TYPE* stack_;
//     unsigned int capacity;
//     unsigned int size;
//     //@}

//     /*!@name Constructor*/
//     //@{
//     Vector()
//     {
//       capacity = 0;
//       size = 0;
//       stack_ = NULL;
//     }

//     Vector(const Vector<DATA_TYPE>& s)
//     {
//       capacity = s.size;
//       stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
//       //stack_ = new DATA_TYPE[capacity];
//       for(size = 0; size<capacity; ++size)
// 	stack_[size] = s[size];
//     }

//     Vector(const int n)
//     {
//       capacity = n;
//       size = n;
//       if( capacity ) {
// 	stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
// 	//stack_ = new DATA_TYPE[capacity];
// 	int f = sizeof(DATA_TYPE)/sizeof(int);
// 	std::fill((int*)stack_, (int*)stack_+(capacity*f), 0);
//       }
//       else stack_ = NULL;
//     }
//     //@}

//     /*!@name Destructor*/
//     //@{  
//     virtual ~Vector()
//     {
// #ifdef _DEBUG_MEMORY
//       std::cout << "c delete vector: " << size << " " << capacity << " " // ;
//       // display(std::cout);
//       // std::cout 
// 	<< std::endl;
// #endif
//       free( stack_ );
//       //delete [] stack_;
//     }
//     //@}

//     /*!@name Initialisation*/
//     //@{
//     void initialise(const unsigned int c)
//     {
//       size = 0;
//       capacity = c;

//       //std::cout << "init(c): " << capacity << std::endl;

//       stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      
//       // stack_ = new DATA_TYPE[capacity];
      
//       int f = sizeof(DATA_TYPE)/sizeof(int);
//       std::fill((int*)stack_, (int*)stack_+(capacity*f), 0);

//       //DATA_TYPE x(0);
//       //std::fill(stack_, stack_+capacity, x);
//     }

//     void initialise(const unsigned int s, const unsigned int c)
//     {
//       size = s;
//       capacity = c;
//       stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));

//       //std::cout << "init(s,c): " << capacity << std::endl;

//       // stack_ = new DATA_TYPE[capacity];
      
//       int f = sizeof(DATA_TYPE)/sizeof(int);
//       std::fill((int*)stack_, (int*)stack_+(capacity*f), 0);
      
//       //DATA_TYPE x(0);
//       //std::fill(stack_, stack_+capacity, x);
//     }

//     void extendStack( const unsigned int l=0 )
//     {

//       //std::cout << "extend stack!! " << this << std::endl;

//       unsigned int increment = (l ? l : (capacity+1) << 1);
//       capacity += increment;

//       // DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, capacity*sizeof(DATA_TYPE));
//       // std::cout << (int*)new_stack << " " << (int*)stack_ << std::endl;
//       // stack_ = new_stack;


//       DATA_TYPE* new_stack = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
//       // for(int i=0; i<capacity-increment; ++i)
//       // 	new_stack[i] = stack_[i];

//       memcpy(new_stack, stack_, (capacity-increment)*sizeof(DATA_TYPE));

//       //memcpy((int*)new_stack, (int*)stack_, (capacity-increment)*f);
//       free(stack_); 
//       stack_ = new_stack;

//       // std::cout << "extend: " << capacity << std::endl;
     
//       // DATA_TYPE* new_stack_ = new DATA_TYPE[capacity];
//       // memcpy(new_stack_, stack_, capacity-increment);
//       // delete [] stack_; 
//       // stack_ = new_stack_;
      
//       int f = sizeof(DATA_TYPE)/sizeof(int);
//       std::fill((int*)stack_+(capacity-increment)*f, (int*)stack_+(capacity*f), 0);


//       //std::cout << " ==> " << this << std::endl;

//       //DATA_TYPE x(0);
//       //std::fill(stack_+capacity-increment, stack_+capacity, x);
//     }

//     void resize( const unsigned int l )
//     {
//       if( capacity < l ) {
// 	capacity = l;
// 	DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, capacity*sizeof(DATA_TYPE));
// 	stack_ = new_stack;

// 	// std::cout << "resize: " << l << std::endl;
	
// 	// DATA_TYPE* new_stack_ = new DATA_TYPE[l];
// 	// memcpy(new_stack_, stack_, capacity);
// 	// delete [] stack_;
// 	// stack_ = new_stack_;

// 	// capacity = l;
//       }
//       size = l;
//     }
//     //@}

//     /*!@name Accessors*/
//     //@{
//     inline int empty() const
//     {
//       return !size;
//     }
  
//     inline void add(DATA_TYPE x)
//     {
//       if( capacity == size ) 
// 	extendStack();
//       stack_[size++] = x;
//     }

//     inline void secure(const DATA_TYPE x) 
//     {
//       if( capacity == size ) 
// 	extendStack();
//     }

//     inline void fast_add(DATA_TYPE x)
//     {
//       stack_[size++] = x;
//     }

//     inline void push_back(DATA_TYPE x)
//     {
//       if( capacity == size ) 
// 	extendStack();
//       stack_[size++] = x;
//     }

//     inline DATA_TYPE pop_until(const unsigned int level)
//     {
//       size = level;
//       return stack_[size];
//     }

//     inline DATA_TYPE pop()
//     {
//       return stack_[--size];
//     }

//     inline void pop(DATA_TYPE& x)
//     {
//       x = stack_[--size];
//     }

//     inline void clear()
//     {
//       size = 0;
//     }

//     inline void remove(const unsigned int i)
//     {  
//       stack_[i] = stack_[--size];
//     }

//     inline void remove_elt(DATA_TYPE& elt)
//     {
//       unsigned int j=size;
//       while(j && stack_[--j] != elt);
//       stack_[j] = stack_[--size];
//     }

//     inline void setBack(const DATA_TYPE& x, const int k=1)
//     {
//       stack_[size-k] = x;
//     }

//     inline DATA_TYPE& front(const int k=0)
//     {
//       return stack_[k];
//     }

//     inline const DATA_TYPE front(const int k=0) const
//     {
//       return stack_[k];
//     }

//     inline DATA_TYPE& back(const int k=1)
//     {
//       return stack_[size-k];
//     }

//     inline const DATA_TYPE back(const int k=1) const
//     {
//       return stack_[size-k];
//     }

//     inline DATA_TYPE& operator[](const int i)
//     {
//       return stack_[i];
//     }

//     inline const DATA_TYPE operator[](const int i) const
//     {
//       return stack_[i];
//     }

//     inline Vector< DATA_TYPE >& operator=(const Vector< DATA_TYPE >& x)
//     {
//       initialise(0, x.capacity);
//       for(unsigned int i=0; i<x.size; ++i)
// 	add(x[i]);
//       return *this;
//     }
//     //@}

//     /*!@name Printing*/
//     //@{
//     std::ostream& display(std::ostream& os) const {
//       os << "[";
//       if(size) os << stack_[0] ;
//       for(unsigned int i=1; i<size; ++i)
// 	os << " " << stack_[i];
//       os << "]";
//       return os;
//     }
//     //@}

//   };



  /**********************************************
   * Vector
   **********************************************/
  /*! \class Vector
    \brief Simple vector class     
  */
  //int global = 0;
  
  template <class DATA_TYPE>
  int increasing_order(const void *x, const void *y) {
    DATA_TYPE& x_ = *((DATA_TYPE*)x);
    DATA_TYPE& y_ = *((DATA_TYPE*)y);
    return(x_ < y_ ? -1 : (x_ > y_ ? 1 : 0));
  }
  
  template <class DATA_TYPE>
  class Vector {
  public:

    typedef DATA_TYPE* iterator;
    
    /*!@name Parameters*/
    //@{
    DATA_TYPE* stack_;
    unsigned int capacity;
    unsigned int size;
    //@}

    /*!@name Constructor*/
    //@{
    Vector()
    {
      capacity = 0;
      size = 0;
      stack_ = NULL;
    }

    Vector(const Vector<DATA_TYPE>& s)
    {
      capacity = s.size;
      //stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      stack_ = new DATA_TYPE[capacity];
      for(size = 0; size<capacity; ++size)
	stack_[size] = s[size];
    }

    Vector(const int n)
    {
      capacity = n;
      size = n;
      if( capacity ) {
	//stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
	stack_ = new DATA_TYPE[capacity];
	//int f = sizeof(DATA_TYPE)/sizeof(int);
	//std::fill((int*)stack_, (int*)stack_+(capacity*f), 0);
	std::fill(stack_, stack_+capacity, (DATA_TYPE)0);
      }
      else stack_ = NULL;
    }
    //@}

    /*!@name Destructor*/
    //@{  
    virtual ~Vector()
    {
#ifdef _DEBUG_MEMORY
      std::cout << "c delete vector: " << size << " " << capacity << " " // ;
      // display(std::cout);
      // std::cout 
	<< std::endl;
#endif
      //free( stack_ );
      delete [] stack_;
    }
    //@}

    /*!@name Initialisation*/
    //@{
    void initialise(const unsigned int c)
    {
      size = 0;
      capacity = c;
      stack_ = new DATA_TYPE[capacity];
      std::fill(stack_, stack_+capacity, DATA_TYPE());
    }

    void initialise(const unsigned int s, const unsigned int c)
    {
      size = s;
      capacity = c;
      stack_ = new DATA_TYPE[capacity];
      std::fill(stack_, stack_+capacity, DATA_TYPE());
    }

    void initialise(const unsigned int s, const unsigned int c, DATA_TYPE x)
    {
      size = s;
      capacity = c;
      stack_ = new DATA_TYPE[capacity];
      std::fill(stack_, stack_+capacity, x);
    }

    void copy(Vector<DATA_TYPE>& vec) {
      size = vec.size;
      capacity = vec.capacity;
      stack_ = vec.stack_;
    }

    void neutralise() {
      stack_ = NULL;
    }

    void extendStack( const unsigned int l=0 )
    {

      // if(size > 1400) {
      // std::cerr << "extend " << size << std::endl;
      // //display(std::cout);
      // std::cout << std::endl;

      unsigned int increment = (l ? l : (capacity+1) << 1);
      capacity += increment;

      DATA_TYPE* new_stack = new DATA_TYPE[capacity];
      for(unsigned int i=0; i<capacity-increment; ++i)
       	new_stack[i] = stack_[i];

      delete [] stack_;
      stack_ = new_stack;
      
      std::fill(stack_+capacity-increment, stack_+capacity, DATA_TYPE());


      // std::cout << "extended " ;
      // //display(std::cout);
      // std::cout << std::endl << std::endl;


    }

    void resize( const unsigned int l )
    {
      if( capacity < l ) {
	DATA_TYPE* new_stack_ = new DATA_TYPE[l];
	memcpy(new_stack_, stack_, capacity);
	delete [] stack_;
	stack_ = new_stack_;
	std::fill(stack_+capacity, stack_+l, (DATA_TYPE)0);

	capacity = l;
      }
      size = l;
    }
    //@}

    /*!@name Accessors*/
    //@{

    void sort() {
      qsort(stack_, size, sizeof(DATA_TYPE), increasing_order<DATA_TYPE>);
    }

    inline iterator begin() const {
      return stack_;
    }
    inline iterator end() const {
      return stack_+size;
    }

    inline int empty() const
    {
      return !size;
    }
  
    inline void add(DATA_TYPE x)
    {
      if( capacity == size ) 
	extendStack();
      stack_[size++] = x;
    }

    inline void secure(const DATA_TYPE x) 
    {
      if( capacity == size ) 
	extendStack();
    }

    inline void fast_add(DATA_TYPE x)
    {
      stack_[size++] = x;
    }

    inline void push_back(DATA_TYPE x)
    {
      if( capacity == size ) 
	extendStack();
      stack_[size++] = x;
    }

    inline DATA_TYPE pop_until(const unsigned int level)
    {
      size = level;
      return stack_[size];
    }

    inline DATA_TYPE pop()
    {
      return stack_[--size];
    }

    inline void pop(DATA_TYPE& x)
    {
      x = stack_[--size];
    }

    inline void clear()
    {
      size = 0;
    }

    inline void remove(const unsigned int i)
    {  
      stack_[i] = stack_[--size];
    }

    inline void remove_elt(DATA_TYPE& elt)
    {
      unsigned int j=size;
      while(j && stack_[--j] != elt);
      stack_[j] = stack_[--size];
    }

    inline void set_back(const DATA_TYPE& x, const int k=1)
    {
      stack_[size-k] = x;
    }

    inline DATA_TYPE& front(const int k=0)
    {
      return stack_[k];
    }

    inline const DATA_TYPE front(const int k=0) const
    {
      return stack_[k];
    }

    inline DATA_TYPE& back(const int k=1)
    {
      return stack_[size-k];
    }

    inline const DATA_TYPE back(const int k=1) const
    {
      return stack_[size-k];
    }

    inline DATA_TYPE& operator[](const int i)
    {
      return stack_[i];
    }

    inline const DATA_TYPE operator[](const int i) const
    {
      return stack_[i];
    }

    inline Vector< DATA_TYPE >& operator=(const Vector< DATA_TYPE >& x)
    {
      if(x.capacity && !stack_) {
	initialise(0, x.capacity);
      } else if(capacity<x.size) {
	extendStack(x.capacity-capacity);
      }

      clear();

      for(unsigned int i=0; i<x.size; ++i)
	add(x[i]);
      
      return *this;
    }
    //@}

    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const {
      os << "[";
      if(size) os << stack_[0] ;
      for(unsigned int i=1; i<size; ++i)
	os << " " << stack_[i];
      os << "]";
      return os;
    }
    //@}

  };
	
	
	template< class T >
	class TwoWatchedThing {
	
	public:

		Vector<T> * watched_by;
		Vector<int> * watched_other;
		Vector<int> * watched_index;
	
	
	
		TwoWatchedThing() {
		}
		TwoWatchedThing(const int n) {
			initialise(n);
		}
	
		virtual ~TwoWatchedThing() {
			delete [] watched_by;
			delete [] watched_other;
			delete [] watched_index;
		}
	
		void initialise(const int n) {
			watched_by = new Vector<T>[n];
			watched_other = new Vector<int>[n];
			watched_index = new Vector<int>[n];
		}

		void watch(const T x, const int w1, const int w2) {
		
			watched_index[w1].add(watched_by[w2].size);
			watched_index[w2].add(watched_by[w1].size);
			watched_other[w1].add(w2);
			watched_other[w2].add(w1);
			watched_by[w1].add(x);
			watched_by[w2].add(x);
		
		}
	
	
		void unwatch(const int ix1, const int w1) {
							
			// second watcher of x, [OK]
			int  w2 = watched_other[w1][ix1];
			int x = watched_by[w1][ix1];	


			
#ifdef _DEBUG_UNWATCH
			
			
			if(_DEBUG_UNWATCH) {
				std::cout << " stop watching " << w1 << " and " << w2 << " with " << watched_by[w1][ix1] << std::endl;

				std::cout << w1 << ":\n";
				std::cout << watched_by[w1] << std::endl;
				std::cout << watched_index[w1] << std::endl;
				std::cout << watched_other[w1] << std::endl;
				std::cout << std::endl;

				std::cout << w2 << ":\n";
				std::cout << watched_by[w2] << std::endl;
				std::cout << watched_index[w2] << std::endl;
				std::cout << watched_other[w2] << std::endl;
				std::cout << std::endl;
			}
#endif
			
			
			// elt to swap with x in w1's list, [OK] 
			T l1 = watched_by[w1].pop();
			// other watcher of l1 (besides w1), [OK]
			int o2 = watched_other[w1].pop();		
			// elt to swap with x in w2's list, [OK] 
			T l2 = watched_by[w2].pop();
			// other watcher of l2 (besides w2), [OK]
			int o1 = watched_other[w2].pop();
			
			
			// index of x in w2's list [MUT!!]
			int ix2 = watched_index[w1][ix1];
			
			
			if(l1!=x) {
#ifdef _DEBUG_UNWATCH
				if(_DEBUG_UNWATCH) {		
					std::cout << " move " << l1 << " to "<< watched_by[w1][ix1] << "'s slot in " << w1 << "'s list" << std::endl << "=>\n";
				}
#endif	

				// move l1 to x's slot in w1
				watched_by[w1][ix1] = l1;
				watched_other[w1][ix1] = o2;
				int il1o2 = watched_index[w1].pop();
				watched_index[w1][ix1] = il1o2;
			
#ifdef _DEBUG_UNWATCH
				if(_DEBUG_UNWATCH) {		
					std::cout << w1 << ":\n";
					std::cout << watched_by[w1] << std::endl;
					std::cout << watched_index[w1] << std::endl;
					std::cout << watched_other[w1] << std::endl;
				}
#endif			
			
				// change its pointer in o2
				watched_index[o2][il1o2] = ix1;
			
#ifdef _DEBUG_UNWATCH
				if(_DEBUG_UNWATCH) {	
					std::cout << " and repoint its counterpart in " << o2 << "'s list" << std::endl;
					std::cout << o2 << ":\n";
					std::cout << watched_by[o2] << std::endl;
					std::cout << watched_index[o2] << std::endl;
					std::cout << watched_other[o2] << std::endl;
					std::cout << std::endl;
				}
#endif			
			} else watched_index[w1].pop();
	
			if(l2!=x) {		
#ifdef _DEBUG_UNWATCH
				if(_DEBUG_UNWATCH) {		
					std::cout << " move " << l2 << " to "<< watched_by[w2][ix2] << "'s slot in " << w2 << "'s list" << std::endl << "=>\n";
				}
#endif	
			
				// move l2 to x's slot in w2
				watched_by[w2][ix2] = l2;
				watched_other[w2][ix2] = o1;
				int il2o1 = watched_index[w2].pop();
				watched_index[w2][ix2] = il2o1;
			
#ifdef _DEBUG_UNWATCH
				if(_DEBUG_UNWATCH) {		
					std::cout << w2 << ":\n";
					std::cout << watched_by[w2] << std::endl;
					std::cout << watched_index[w2] << std::endl;
					std::cout << watched_other[w2] << std::endl;
				}
#endif			
			
				// change its pointer in o1
				watched_index[o1][il2o1] = ix2;
			
#ifdef _DEBUG_UNWATCH
				if(_DEBUG_UNWATCH) {	
					std::cout << " and repoint its counterpart in " << o1 << "'s list" << std::endl;
					std::cout << o1 << ":\n";
					std::cout << watched_by[o1] << std::endl;
					std::cout << watched_index[o1] << std::endl;
					std::cout << watched_other[o1] << std::endl;
					std::cout << std::endl;
				}
#endif		
			} else watched_index[w2].pop();
			
			
	//
//
//
//
// 			int  w2 = watched_other[w1][ix1];
//
// #ifdef _DEBUG_UNWATCH
// 			if(_DEBUG_UNWATCH) {
// 			std::cout << " stop watching " << w1 << " and " << w2 << " with " << watched_by[w1][ix1] << std::endl;
//
// 			std::cout << w1 << ":\n";
// 			std::cout << watched_by[w1] << std::endl;
// 			std::cout << watched_index[w1] << std::endl;
// 			std::cout << watched_other[w1] << std::endl;
// 			std::cout << std::endl;
//
// 			std::cout << w2 << ":\n";
// 			std::cout << watched_by[w2] << std::endl;
// 			std::cout << watched_index[w2] << std::endl;
// 			std::cout << watched_other[w2] << std::endl;
// 			std::cout << std::endl;
// 		}
// #endif
//
//
// 			int ix2 = watched_index[w1][ix1];
//
// 			T l1 = watched_by[w1].pop();
// 			T l2 = watched_by[w2].pop();
//
// 			int o1 = watched_other[w2].pop();
// 			int o2 = watched_other[w1].pop();
//
// 			int il1o2 = watched_index[w1].pop();
// 			int il2o1 = watched_index[w2].pop();
//
//
//
// #ifdef _DEBUG_UNWATCH
// 			if(_DEBUG_UNWATCH) {
// 			std::cout << "the other watched of " << l1 << " and " << l2 << " are, respectively, " << o2 << " and " << o1 << std::endl;
//
//
// 			std::cout << o2 << ":\n";
// 			std::cout << watched_by[o2] << std::endl;
// 			std::cout << watched_index[o2] << std::endl;
// 			std::cout << watched_other[o2] << std::endl;
// 			std::cout << std::endl;
//
// 			std::cout << o1 << ":\n";
// 			std::cout << watched_by[o1] << std::endl;
// 			std::cout << watched_index[o1] << std::endl;
// 			std::cout << watched_other[o1] << std::endl;
// 			std::cout << std::endl;
//
// 			std::cout << "the index of " << l1 << " in " << w1 <<"'s list has changed to " << ix1 << std::endl;
// 			std::cout << "the index of " << l2 << " in " << w2 <<"'s list has changed to " << ix2 << std::endl;
// 		}
// #endif
//
// 			//
// 			// change the index of l1 in watched[o2]
// 			watched_index[o2][il1o2] = ix1;
//
// 			// change the index of l2 in watched[o1]
// 			watched_index[o1][il2o1] = ix2;
//
//
// #ifdef _DEBUG_UNWATCH
// 			if(_DEBUG_UNWATCH) {
// 			std::cout << o2 << ":\n";
// 			std::cout << watched_by[o2] << std::endl;
// 			std::cout << watched_index[o2] << std::endl;
// 			std::cout << watched_other[o2] << std::endl;
// 			std::cout << std::endl;
//
// 			std::cout << o1 << ":\n";
// 			std::cout << watched_by[o1] << std::endl;
// 			std::cout << watched_index[o1] << std::endl;
// 			std::cout << watched_other[o1] << std::endl;
// 			std::cout << std::endl;
// 		}
// #endif
//
//
//
// 			// move l1 in watched[w1]
// 			watched_by[w1][ix1] = l1;
// 			watched_other[w1][ix1] = o2;
// 			watched_index[w1][ix1] = il1o2;
//
//
// 			// move l2 in watched[w2]
// 			watched_by[w2][ix2] = l2;
// 			watched_other[w2][ix2] = o1;
// 			watched_index[w2][ix2] = il2o1;
//
//
// #ifdef _DEBUG_UNWATCH
// 			if(_DEBUG_UNWATCH) {
// 			std::cout << " swap " << watched_by[w1][ix1] << " with " << l1 << " and " << l2 << ", respectively" << std::endl << "=>\n";
// 			std::cout << w1 << ":\n";
// 			std::cout << watched_by[w1] << std::endl;
// 			std::cout << watched_index[w1] << std::endl;
// 			std::cout << watched_other[w1] << std::endl;
// 			std::cout << std::endl;
//
// 			std::cout << w2 << ":\n";
// 			std::cout << watched_by[w2] << std::endl;
// 			std::cout << watched_index[w2] << std::endl;
// 			std::cout << watched_other[w2] << std::endl;
// 			std::cout << std::endl;
// 		}
// #endif


			
			
		}


	};


template < int N, class T >
  
  class Tuple {
  
 public:
  
  T data[N];
  
  Tuple() {}
  Tuple(T a) {data[0] = a;}
  Tuple(T a, T b) {data[0] = a; data[1] = b;}
  Tuple(T a, T b, T c) {data[0] = a; data[1] = b; data[2] = c;}
  Tuple(T a, T b, T c, T d) {data[0] = a; data[1] = b; data[2] = c; data[3] = d;}
  Tuple(T a, T b, T c, T d, T e) {data[0] = a; data[1] = b; data[2] = c; data[3] = d; data[4] = e;}

  T operator[](const int i) const { return data[i]; }
  T& operator[](const int i) { return data[i]; }
  
  /*!@name Printing*/
  //@{
  std::ostream& display(std::ostream& os) const {
    os << "<";
    if(N) os << data[0] ;
    for(unsigned int i=1; i<N; ++i)
      os << " " << data[i];
    os << ">";
    return os;
  }
  //@}

};


  typedef Tuple< 2, int > Pair;
  
  template <class T1, class T2, class T3>
  class Triplet {
  public:
    T1 first;
    T2 second;
    T3 third;

    Triplet() {}
    Triplet(T1 t1, T2 t2, T3 t3) : first(t1), second(t2), third(t3) {}
    ~Triplet() {}

    operator T1() { return first; }

   std::ostream& display(std::ostream& os) const {
      os << "<" << first << ", " << second << ", " << third << ">";
      return os;
    }    

  };


  template <class T>
  class TwoWayStack {
  public:

    /*!@name Parameters*/
    //@{
    /** <stack_> is an array of object T such that:
	- The first element of the stack is the one at index <start>
	- There are <size> elements;
	- The current size of the stack is <capacity>
	Moreover, we assume that elements T can be casted into int, and that
	this integer is between 0 and <capacity>, then we have:
	- <index_[i]> is the index_ of element i in <stack_>
     */

    T   *stack_;
    int *index_;

    // index
    unsigned int start;
    unsigned int capacity;
    unsigned int size;
    //@}

    /*!@name Constructor*/
    //@{
    TwoWayStack()
    {
      start = 0;
      capacity = 0;
      size = 0;
      stack_ = NULL;
      index_ = NULL;
    }
    //@}

    /*!@name Destructor*/
    //@{  
    virtual ~TwoWayStack()
    {
      free( stack_ );
      free( index_ );
    }
    //@}

    /*!@name Initialisation*/
    //@{
    void initialise(const unsigned int c)
    {
      start = 0;
      size = 0;
      capacity = c;

      stack_ = (T  *) malloc(capacity*sizeof(T  ));
      index_ = (int*) malloc(capacity*sizeof(int));
      
      int f = sizeof(T)/sizeof(int);

      //int out = (start+size+1)%capacity;

      std::fill((int*)stack_, (int*)stack_+(capacity*f), 0);
      std::fill(index_, index_+capacity, -1);



      // std::cout << "CREATE " << std::endl;
      // for(int i=0; i<capacity; ++i)
      //  	std::cout << index_[i] << " ";
      // std::cout << std::endl;
    }


    // void check_intergrity() {
    //   // check that elements 
    // }


    void extend_stack( const unsigned int l=0 )
    {
      unsigned int increment = (l ? l : (capacity+1) << 1);
      int f = sizeof(T)/sizeof(int);
      //capacity += increment;

      int* new_index = (int*) malloc((capacity+increment)*sizeof(int));
      for(unsigned int i=0; i<capacity; ++i) {
	new_index[i] = (index_[i]>=0 ? (index_[i]+capacity-start)%capacity : -1);
      }
      // for(int i=start; i<capacity; ++i) {
      // 	new_index[i] = index_[i]-start;
      // }
      //int out = (start+size+1)%(capacity+increment);
      std::fill(new_index+capacity, new_index+capacity+increment, -1);

      T* new_stack = (T*) malloc((capacity+increment)*sizeof(T));
      memcpy(new_stack, stack_+start, (capacity-start)*sizeof(T));
      memcpy(new_stack+capacity-start, stack_, (start)*sizeof(T));
      std::fill((int*)new_stack+(capacity*f), (int*)new_stack+(capacity+increment)*f, 0);

      free(stack_); 
      stack_ = new_stack;

      free(index_); 
      index_ = new_index;

      capacity += increment;
      start = 0;

      // std::cout << "EXTEND" << std::endl;
      // for(int i=0; i<capacity; ++i)
      //  	std::cout << index_[i] << " ";
      // std::cout << std::endl;

    }

    void declare( const int x ) {
      if(x >= (int)capacity) extend_stack();
    }
    //@}

    /*!@name Accessors*/
    //@{
    inline int empty() const
    {
      return !size;
    }

    inline void probe() { if(size == capacity) extend_stack(); }
  
    inline void push_back(T x)
    {

      //if((int)x == 355) std::cout << "PUSH B x355" << std::endl; 

      int idx = (start+(size++))%capacity;
      stack_[idx] = x;
      index_[(int)x] = idx;
    }

    inline void push_front(T x)
    {

      //if((int)x == 355) std::cout << "PUSH F x355" << std::endl; 

      ++size;
      start = ((start+capacity-1)%capacity);
      stack_[start] = x;

      index_[(int)x] = start;
    }

    inline bool contain(const int x) {

      // if(x == 355) {

      // 	std::cout << "contain 355? " << index_[x] << std::endl;

      // }


      return index_[x] >= 0;

      //return (int)(stack_[index_[x]]) == x;

      // return (( capacity+
      // 		index_[x]-start)%capacity) < size;
    }

    inline T& operator[](const int x) {
      return stack_[index_[x]];
    }

    inline T pop_back()
    {
      int idx = (start+--size)%capacity;
      index_[(int)stack_[idx]] = -1;
      return stack_[idx];
    }

    inline T pop_front()
    {
      --size;
      T rval = stack_[start];
      index_[(int)stack_[start]] = -1;
      start = (start+1)%capacity;
      return rval;
    }

    inline void pop_back(T& x)
    {
      int idx = (start+--size)%capacity;
      index_[(int)stack_[idx]] = -1;
      //index_[idx] = -1;
      x = stack_[idx];
    }

    inline void pop_front(T& x)
    {
      --size;
      x = stack_[start];
      index_[(int)stack_[start]] = -1;
      //index_[start] = -1;
      start = (start+1)%capacity;
    }

    inline void clear()
    {
      for(unsigned int i=0; i<size; ++i) {
       	index_[stack_[(start+i)%capacity]] = -1;
      }
      size = 0;
    }

    //@}

    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const {
      os << "[";
      if(size) os << stack_[start]  ;
      for(unsigned int i=1; i<size; ++i)
	os << " " << stack_[(start+i)%capacity] ;
      os << "]";
      return os;
    }
    //@}

  };




  template <class MAIN_TYPE, class AUX_TYPE>
  class BiStack {
  public:

    /*!@name Parameters*/
    //@{
    MAIN_TYPE* main_stack_;
    AUX_TYPE* aux_stack_;
    unsigned int start;
    unsigned int capacity;
    unsigned int size;
    //@}

    /*!@name Constructor*/
    //@{
    BiStack()
    {
      start = 0;
      capacity = 0;
      size = 0;
      main_stack_ = NULL;
      aux_stack_ = NULL;
    }
    //@}

    /*!@name Destructor*/
    //@{  
    virtual ~BiStack()
    {
      free( main_stack_ );
      free( aux_stack_ );
    }
    //@}

    /*!@name Initialisation*/
    //@{
    void initialise(const unsigned int c)
    {
      start = 0;
      size = 0;
      capacity = c;

      main_stack_ = (MAIN_TYPE*) malloc(capacity*sizeof(MAIN_TYPE));
      aux_stack_ = (AUX_TYPE*) malloc(capacity*sizeof(AUX_TYPE));
      
      int f = sizeof(MAIN_TYPE)/sizeof(int);
      std::fill((int*)main_stack_, (int*)main_stack_+(capacity*f), 0);
      f = sizeof(AUX_TYPE)/sizeof(int);
      std::fill((int*)aux_stack_, (int*)aux_stack_+(capacity*f), 0);
    }


    void extend_stack( const unsigned int l=0 )
    {
      unsigned int increment = (l ? l : (capacity+1) << 1);
      capacity += increment;

      MAIN_TYPE* new_main_stack = (MAIN_TYPE*) malloc(capacity*sizeof(MAIN_TYPE));
      memcpy(new_main_stack, main_stack_, (capacity-increment)*sizeof(MAIN_TYPE));

      AUX_TYPE* new_aux_stack = (AUX_TYPE*) malloc(capacity*sizeof(AUX_TYPE));
      memcpy(new_aux_stack, aux_stack_, (capacity-increment)*sizeof(AUX_TYPE));

      free(main_stack_); 
      main_stack_ = new_main_stack;

      free(aux_stack_); 
      aux_stack_ = new_aux_stack;

      int f = sizeof(MAIN_TYPE)/sizeof(int);
      std::fill((int*)main_stack_, (int*)main_stack_+(capacity*f), 0);
      f = sizeof(AUX_TYPE)/sizeof(int);
      std::fill((int*)aux_stack_, (int*)aux_stack_+(capacity*f), 0);
    }
    //@}

    /*!@name Accessors*/
    //@{
    inline int empty() const
    {
      return !size;
    }

    inline void probe() { if(size == capacity) extend_stack(); }
  
    inline void push_back(MAIN_TYPE x, AUX_TYPE y)
    {
      int idx = (start+(size++))%capacity;
      main_stack_[idx] = x;
      aux_stack_[idx] = y;
    }

    inline void push_front(MAIN_TYPE x, AUX_TYPE y)
    {
      ++size;
      start = (start+capacity-1)%capacity;
      main_stack_[start] = x;
      aux_stack_[start] = y;
    }


    inline void pop_back(MAIN_TYPE& x, AUX_TYPE& y)
    {
      int idx = (start+--size)%capacity;
      x = main_stack_[idx];
      y = aux_stack_[idx];
    }

    inline void pop_front(MAIN_TYPE& x, AUX_TYPE& y)
    {
      --size;
      x = main_stack_[start];
      y = aux_stack_[start];
      start = (start+1)%capacity;
    }

    inline void clear()
    {
      size = 0;
    }

    //@}

    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const {
      for(unsigned int i=0; i<capacity; ++i)
	os << main_stack_[i] << " ";
      os << std::endl;

      os << "[";
      if(size) os << main_stack_[start] << "/" << aux_stack_[start] ;
      for(unsigned int i=1; i<size; ++i)
	os << " " << main_stack_[(start+i)%capacity] << "/" << aux_stack_[(start+i)%capacity];
      os << "]";
      return os;
    }
    //@}

  };



	class IntList {
	public: 

		int capacity;
		int *next;
		int *prev;
		int size;
	
		//Vector<int> removed;
	
		IntList() {}
		IntList(const int n, const bool filled=false) {
			initialise(n,filled);
		}
	
		virtual ~IntList() {
			delete [] next;
			delete [] prev;
		}
	
		void initialise(const int n, const bool filled=false) {
			capacity = n;
			next = new int[capacity+1];
			prev = new int[capacity+1];
		
			if(filled) {
				fill();
			} else {
				clear();
			}
			
#ifdef _DEBUG_BACKJUMPS			
			verify("after init");
#endif
		}
	
		bool contain(const int elt) const {
			return next[elt] != elt;
		}
	
		int get_first() const { return next[capacity]; }
		int get_last () const { return prev[capacity]; }
		
		void add_after(const int elt, const int ref) {
#ifdef _DEBUG_BACKJUMPS				
			verify("before add_after");
#endif			
			int n = next[ref];			
			prev[n] = elt;
			next[elt] = n;
			next[ref] = elt;
			prev[elt] = ref;
			++size;
#ifdef _DEBUG_BACKJUMPS				
			verify("after add_after");
#endif
		}
		
		void add_before(const int elt, const int ref) {
#ifdef _DEBUG_BACKJUMPS				
			verify("before add_before");
#endif			
			int p = prev[ref];			
			next[p] = elt;
			prev[elt] = p;
			next[elt] = ref;
			prev[ref] = elt;
			++size;
#ifdef _DEBUG_BACKJUMPS				
			verify("after add_before");
#endif
		}
		
		void add_end(const int elt) {
#ifdef _DEBUG_BACKJUMPS				
			verify("before add_end");
#endif			
			int p = prev[capacity];			
			next[p] = elt;
			prev[elt] = p;
			next[elt] = capacity;
			prev[capacity] = elt;
			++size;		
#ifdef _DEBUG_BACKJUMPS				
			verify("after add_end");
#endif
		}
		
		void add_beg(const int elt) {
#ifdef _DEBUG_BACKJUMPS				
			verify("before add_beg");
#endif			
			int n = next[capacity];			
			prev[n] = elt;
			next[elt] = n;
			next[capacity] = elt;
			prev[elt] = capacity;
			++size;			
#ifdef _DEBUG_BACKJUMPS				
			verify("after add_beg");
#endif
		}
	
		void remove(const int elt) {
#ifdef _DEBUG_BACKJUMPS				
			verify("before remove");
#endif			
			int p = prev[elt], n=next[elt];
			next[p] = n;
			prev[n] = p;
			next[elt] = elt;
			--size;
#ifdef _DEBUG_BACKJUMPS				
			verify("after remove");
#endif			
			//removed.add(elt);
		}
	
		void fill() {
			size = capacity;
			for(int i=0; i<=capacity; ++i) {
				next[i] = (i+1)%(capacity+1);
				prev[i] = (capacity+i)%(capacity+1);
			}
			// int i, elt, x;
			// for(i=0; i<removed.size; ++i) {
			// 	elt = removed[i];
			//
			// 	x = (elt+1)%(capacity+1);
			// 	next[elt] = x;
			// 	prev[x] = elt;
			//
			// 	x = (capacity+elt)%(capacity+1);
			// 	prev[elt] = x;
			// 	next[x] = elt;
			// }
			// removed.clear();
		}
		
		
		void clear() {
			size = 0;
			for(int i=0; i<=capacity; ++i) {
				next[i] = i;
				prev[i] = i;
			}
		}
		
		void quick_clear() {
#ifdef _DEBUG_BACKJUMPS				
			verify("before quick clear");
#endif			
			size = 0;
			next[capacity] = capacity;
			prev[capacity] = capacity;
#ifdef _DEBUG_BACKJUMPS				
			verify("after quick clear");
#endif
		}
		
		
		void verify(const char* msg) {
			int rsize = 0;
			int elt = capacity;
			do {
				
				elt = next[elt];
				if(elt == capacity) break;
				
				if(rsize>size) {
					std::cout << "INIFINITE LIST!! " << msg << std::endl ;
					display(std::cout);
					std::cout << std::endl;
					exit(1);
				}
				++rsize;
				
				
			} while(true);
			
			
			if(rsize<size) {
				std::cout << "LIST SIZE INCORRECT!! (" << size << ") " << msg << std::endl ;
					display(std::cout);
					std::cout  << std::endl;
				exit(1);
			}
			
		}
		
		

    std::ostream& display(std::ostream& os) const {
			
			int rsize = 0;
			int elt = capacity;
      os << "[";
			do {
				elt = next[elt];
				if(elt == capacity) break;
				if(prev[elt] != capacity) os << " ";
				os << elt;
				
				if(rsize++>size) break;
			} while(true);
			os << "]";
						
      return os;
    }
		

	};
	
	std::ostream& operator<< (std::ostream& os, const IntList& x);
	std::ostream& operator<< (std::ostream& os, const IntList* x);



  /**********************************************
   * Vector2
   **********************************************/
  /*! \class Vector2
    \brief Simple vector class     
  */
  
  template <class DATA_TYPE>
  class Vector2 {
  public:

    /*!@name Parameters*/
    //@{
    DATA_TYPE* stack_;
    unsigned int capacity;
    unsigned int size;
    //@}

    /*!@name Constructor*/
    //@{
    Vector2()
    {
      capacity = 0;
      size = 0;
      stack_ = NULL;
    }

    Vector2(const Vector2<DATA_TYPE>& s)
    {
      capacity = s.size;
      stack_ = new DATA_TYPE[capacity];
      //(DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      for(size = 0; size<capacity; ++size)
  	stack_[size] = s[size];
    }

    Vector2(const int n)
    {
      capacity = n;
      size = n;
      if( capacity )
  	//stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
  	stack_ = new DATA_TYPE[capacity]; //(DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      else stack_ = NULL;
    }
    //@}

    /*!@name Destructor*/
    //@{  
    virtual ~Vector2()
    {

#ifdef _DEBUG_MEMORY
      std::cout << "c delete vector: " << size << " " << capacity << " " // ;
      // display(std::cout);
      // std::cout 
	<< std::endl;
#endif

      delete [] stack_;
      //free( stack_ );
    }
    //@}

    /*!@name Initialisation*/
    //@{
    void initialise(const unsigned int c)
    {
      size = 0;
      capacity = c;
      //stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      stack_ = new DATA_TYPE[capacity];

      DATA_TYPE x(0);
      std::fill(stack_, stack_+capacity, x);
    }

    void initialise(const unsigned int s, const unsigned int c)
    {
      size = s;
      capacity = c;
      //stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      stack_ = new DATA_TYPE[capacity];

      DATA_TYPE x(0);
      std::fill(stack_, stack_+capacity, x);
    }

    //     void extendStack()
    //     {
    //       unsigned int new_capacity = ((capacity+1) << 1);
    //       DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, new_capacity*sizeof(DATA_TYPE));
    //       stack_ = new_stack;

    //       DATA_TYPE x;
    //       std::fill(stack_+capacity, stack_+new_capacity, x);
      
    //       capacity = new_capacity;
    //     }

    void extendStack( const unsigned int l=0 )
    {
      unsigned int increment = (l ? l : (capacity+1) << 1);
      capacity += increment;
      //DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, capacity*sizeof(DATA_TYPE));
      DATA_TYPE* new_stack_ = new DATA_TYPE[capacity];
      memcpy(new_stack_, stack_, capacity-increment);
      delete [] stack_;
      stack_ = new_stack_;

      DATA_TYPE x(0);
      std::fill(stack_+capacity-increment, stack_+capacity, x);
    }

    void resize( const unsigned int l )
    {
      if( capacity < l ) {
	int old_capacity = capacity;
  	capacity = l;
  	// DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, capacity*sizeof(DATA_TYPE));
  	// stack_ = new_stack;

	DATA_TYPE* new_stack_ = new DATA_TYPE[capacity];
	memcpy(new_stack_, stack_, old_capacity);
	delete [] stack_;
	stack_ = new_stack_;
	
	DATA_TYPE x(0);
	std::fill(stack_+old_capacity, stack_+capacity, x);
      }
      size = l;
    }
    //@}

    /*!@name Accessors*/
    //@{
    inline int empty() const
    {
      return !size;
    }
  
    inline void add(DATA_TYPE x)
    {
      if( capacity == size ) 
  	extendStack();
      stack_[size++] = x;
    }

    inline void secure(const DATA_TYPE x) 
    {
      if( capacity == size ) 
  	extendStack();
    }

    inline void fast_add(DATA_TYPE x)
    {
      stack_[size++] = x;
    }

    inline void push_back(DATA_TYPE x)
    {
      if( capacity == size ) 
  	extendStack();
      stack_[size++] = x;
    }

    inline DATA_TYPE pop_until(const unsigned int level)
    {
      size = level;
      return stack_[size];
    }

    inline DATA_TYPE pop()
    {
      return stack_[--size];
    }

    inline void pop(DATA_TYPE& x)
    {
      x = stack_[--size];
    }

    inline void clear()
    {
      size = 0;
    }

    inline void remove(const unsigned int i)
    {  
      stack_[i] = stack_[--size];
    }

    inline void remove_elt(DATA_TYPE& elt)
    {
      unsigned int j=size;
      while(j && stack_[--j] != elt);
      stack_[j] = stack_[--size];
    }

    inline void setBack(const DATA_TYPE& x, const int k=1)
    {
      stack_[size-k] = x;
    }

    inline DATA_TYPE& front(const int k=0)
    {
      return stack_[k];
    }

    inline const DATA_TYPE front(const int k=0) const
    {
      return stack_[k];
    }

    inline DATA_TYPE& back(const int k=1)
    {
      return stack_[size-k];
    }

    inline const DATA_TYPE back(const int k=1) const
    {
      return stack_[size-k];
    }

    inline DATA_TYPE& operator[](const int i)
    {
      return stack_[i];
    }

    inline const DATA_TYPE operator[](const int i) const
    {
      return stack_[i];
    }

    inline Vector2< DATA_TYPE >& operator=(const Vector2< DATA_TYPE >& x)
    {
      initialise(0, x.capacity);
      for(unsigned int i=0; i<x.size; ++i)
  	add(x[i]);
      return *this;
    }
    //@}

    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const {
      os << "[";
      if(size) os << stack_[0] ;
      for(unsigned int i=1; i<size; ++i)
  	os << " " << stack_[i];
      os << "]";
      return os;
    }
    //@}

  };


  /**********************************************
   * Array
   **********************************************/    
  /*! \class Array
    \brief Simple array class 
  */

  typedef unsigned int Atom;
  typedef unsigned int Literal;

  //class Decision;
  class Explanation {

  public:
    
    typedef Literal* iterator;

 
    //typedef Decision* iterator;

    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end) = 0;

    virtual std::ostream& display(std::ostream& os) const = 0;

    virtual bool is_clause() {return true;}
    
  };


  
  template <class DATA_TYPE>
  class Array : public Explanation {

  public: 

    /*!@name Parameters*/
    //@{
    unsigned int size;
    DATA_TYPE data[0];
    //@}

    Array(const Vector<DATA_TYPE>& ps) 
    {
      size = ps.size;
      for (unsigned int i=0; i<ps.size; ++i) 
			data[i] = ps[i];
    }
		
    Array(const IntList& ps) 
    {
      size = ps.size;
			int i=0;
			for(int x=ps.get_first(); x!=ps.get_last(); x=ps.next[x]) {
				data[i] = x;
				++i;
			}
			data[i] = ps.get_last();
    }

    virtual ~Array() {}

    // virtual Explanation::iterator begin(Atom a) { return &(data[0]); }
    // virtual Explanation::iterator end  (Atom a) { return &(data[size]); }
    
    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end) { end = (unsigned*)(&(data[size])); return (unsigned*)&(data[0]); }

    static Array<DATA_TYPE>* Array_new(const Vector<DATA_TYPE>& ps)
    {
      void* mem = malloc(sizeof(Array<DATA_TYPE>) + sizeof(DATA_TYPE)*(ps.size));
      return new (mem) Array<DATA_TYPE>(ps); 
    }
		
    static Array<int>* Array_new(const IntList& ps)
    {
      void* mem = malloc(sizeof(Array<int>) + sizeof(int)*(ps.size));
      return new (mem) Array<int>(ps); 
    }

    inline DATA_TYPE& operator [] (const int i) { return data[i]; }
    inline DATA_TYPE operator [] (const int i) const { return data[i]; }

    std::ostream& display(std::ostream& os) const {
      os << "[";
      if(size) os << data[0];
      for(unsigned int i=1; i<size; ++i)
	os << " " << data[i];
      os << "]";
      return os;
    }

  };


  /**********************************************
   * MultiSet
   **********************************************/
  /// Sparse multi-set representation

  class MultiSet 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    int *values_; // contiguous list of values (repeated)
    unsigned int **index_; // index[v][o] is the index of occurrence 'o' of value 'v' in 'values_'
    unsigned int *occurrence_; // self explanatory
    
    /// current max capacity
    unsigned int capacity;
    /// current size (total number of occurrences)
    unsigned int size;

    /// minimum and maximum value 
    int lb_;
    int ub_;
    //@}

    /*!@name Constructors*/
    //@{
    MultiSet()
    {
      size = 0;
      capacity = 0;
      lb_ = 0;
      ub_ = -1;

      values_ = NULL;
      index_ = NULL;
      occurrence_ = NULL;
    }

    MultiSet(const int lb, const int ub, const int* occ=NULL, const int* mocc=NULL)
    {
      if(!occ) {
	int tmp_occ[ub-lb+1];
	int tmp_mocc[ub-lb+1];
	for(int i=0; i<ub-lb+1; ++i) {
	  tmp_occ[i] = 0;
	  tmp_mocc[i] = 4;
	}
	initialise(lb, ub, tmp_occ, tmp_mocc);
      } else if(!mocc) {
	int tmp_mocc[ub-lb+1];
	for(int i=0; i<ub-lb+1; ++i) 
	  tmp_mocc[i] = occ[i];
	initialise(lb, ub, occ, tmp_mocc);
      } else {
	initialise(lb, ub, occ, mocc);
      }
    }

    virtual ~MultiSet()
    {
      delete [] values_;
      occurrence_ += lb_;
      delete [] occurrence_;
      for(int v=lb_; v<=ub_; ++v)
	delete [] index_[v];
      index_ += lb_;
      delete [] index_;
    }

    void initialise(const int lb, const int ub, const int* occ, const int* mocc)
    {
      int span = ub-lb+1, i=0, j=0, v, o;

      lb_ = lb;
      ub_ = ub;

      index_ = new unsigned int*[span];
      occurrence_ = new unsigned int[span];
      capacity = 0;
      size = 0;

      for(v=0; v<span; ++v) {
	capacity += mocc[v];
	size += occ[v];
	occurrence_[v] = occ[v];
	index_[v] = new unsigned int[mocc[v]];
      }
      
      occurrence_ -= lb_;
      index_ -= lb_;
      values_ = new int[capacity];
      
      for(v=lb_; v<=ub_; ++v) {
	for(o=0; o<occ[v-lb_]; ++o) {
	  values_[i] = v;
	  index_[v][o] = i++;
	}
	for(o=occ[v-lb_]; o<mocc[v-lb_]; ++o) {
	  values_[size+j] = v;
	  index_[v][o] = size+j++;
	}
      }
    }
    //@}    


    inline void remove(const int elt)
    {
      --size;
      --occurrence_[elt];

      int v = values_[size];
      int i = index_[elt][occurrence_[elt]];

      index_[v][occurrence_[v]-1] = i;
      index_[elt][occurrence_[elt]] = size;

      values_[i] = values_[size];
      values_[size] = elt;
    }

    inline int pop()
    {
      int v = values_[--size];
      --occurrence_[v];
      return v;
    }

    inline int pop_head()
    {
      const int elt = values_[0];

      --size;
      --occurrence_[elt];

      index_[values_[size]][occurrence_[values_[size]]-1] = 0;
      index_[elt][occurrence_[elt]] = size;

      values_[0] = values_[size];
      values_[size] = elt;

      return elt;
    }

    inline int head()
    {
      return *values_;
    }
    
    inline int back()
    {
      return values_[size-1];
    }

    inline void add(const int elt)
    {
      int i = index_[elt][occurrence_[elt]];
      int v = values_[size]; 

      index_[v][occurrence_[v]] = i;
      index_[elt][occurrence_[elt]] = size;

      values_[i] = v;
      values_[size] = elt;

      ++occurrence_[elt];
      ++size;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    //     std::string getString() const {
    //       std::string return_str = "{{";
    //       bool first = true;
    //       for(int v=lb_; v<=ub_; ++v) {
    // 	if(occurrence_[v] > 0) {
    // 	  if(first) first = false;
    // 	  else return_str += ",";
    // 	  return_str += ("("+toString((int)(occurrence_[v]))+"*"+toString(v)+")");
    // 	}
    //       }
    //       return_str += "}}";
    //       return return_str;
    //     }

    std::ostream& display(std::ostream& os) const {
      os << "{{";
      bool first = true;
      for(int v=lb_; v<=ub_; ++v) {
	if(occurrence_[v] > 0) {
	  if(first) first = false;
	  else os << ",";
	  os << "(" << occurrence_[v] << "*" << v << ")" ;
	}
      }
      os << "}}";
      return os;
    }
    //@}
  };




  /**********************************************
   * IntStack
   **********************************************/
  /// Sparse set representation

  class IntStack 
  {
  public:

    //typedef int* iterator;

    /*!@name Parameters*/
    //@{
    /// list of values
    int *list_;
    /// current max capacity
    unsigned int index_capacity;
    unsigned int list_capacity;
    /// current size
    unsigned int size;
    /// values' indices
    unsigned int *index_;
    unsigned int *start_;
    //@}

    /*!@name Constructors*/
    //@{
    IntStack()
    {
      size = 0;
      index_capacity = 0;
      list_capacity = 0;

      list_ = NULL;
      index_ = NULL;
      start_ = NULL;
    }

    IntStack(const int lb, const int ub, const bool full=true)
    {
      initialise(lb, ub, ub-lb+1, full);
    }

    IntStack(const int lb, const int ub, const int sz, const bool full)
    {
      initialise(lb, ub, sz, full);
    }

    IntStack(IntStack& shared, const int sz)
    {
      initialise(shared, sz);
    }

    void  neutralise() {
      start_ = NULL;
    }
    
    virtual ~IntStack()
    {
      delete [] list_;
      delete [] start_;
    }

    virtual void initialise(IntStack& shared, const int sz);
    virtual void initialise(const int lb, const int ub, const int sz, const bool full);


    void extend_list();
    void extend(const int new_elt);


    /*!@name Accessors*/
    //@{ 
     int get_min() const;

     int get_max() const;

     bool safe_contain(const int elt) const ;

     bool contain(const int elt) const ;
		 
		 bool contain(const int elt, const int min_idx, const int max_idx) const ;
  
     bool empty()const ;

     int next(const int elt) const ;

     int prev(const int elt) const ;
    
     int operator[](const unsigned int idx) const ;

     int& operator[](const unsigned int idx) ;
    //@}

    /*!@name List Manipulation*/
    //@{
     int* begin() ;

     int* end() ;

     int* end_mem() ;

     void fill() ;

     void clear() ;
  
     void set_to(const int elt) ;

     void remove(const int elt) ;
		 void remove(const int elt, int& rk) ; // put elt at rank rk and increase rk

     int next() ;

     int pop() ;

     int pop_head() ;

     int head() const ;
    
     int back() const ;
		 
		 int back(const int offset) const ;

     void init_add(const int elt) ;

     void add(const int elt) ;
		 void add(const int elt, int& rk) ; // put elt at rank rk-1 and decrease rk

     void safe_add(const int elt) ;

    // create a new element that can potentially be outside the bounds
     void create(const int elt) ;

     void ordered_add(const int elt) ;

     void revert_to(const int level) ;

     void index() ;
		 
		 int get_index(const int elt) const ;
		 
		 void move(const int elt, const int idx) ;
		 
		 void set(const int elt, const int idx) ;
    //@}


    /*!@name Miscellaneous*/
    //@{
    std::ostream& display(std::ostream& os) const;
		std::ostream& display(std::ostream& os, const int rank) const;
		std::string to_str() const;
    // std::ostream& display(std::ostream& os) const {
    //   int min =  INFTY;
    //   int max = -INFTY;

    //   for(unsigned int i=0; i<size; ++i) {
    // 	if(list_[i] < min) min = list_[i];
    // 	if(list_[i] > max) max = list_[i];
    //   }

    //   Bitset< unsigned long int > elts(min, max, BitSet::empt);
    //   for(unsigned int i=0; i<size; ++i) {
    // 	elts.add(list_[i]);
    //   }

    //   os << elts;
    //   // os << "(";
    //   // bool not_empty = (size>0);

    //   // if(not_empty) os << list_[0];
    //   // for(unsigned int i=1; i<size; ++i)
    //   // 	os << " " << list_[i];
    //   // os << ")";
    //   return os;
    // }
    //@}
  };
  
	 std::ostream& operator<< (std::ostream& os, const IntStack& x);
  
  




//   //#define TO_ELEMENT(x) ((x) << 3);
//   //#define TO_INDEX(x) ((x) << 11);
//   //#define SIZE 0xfffffff8

// #define SIZE 0xfffffffc
// #define CURRENT 0x3fffffff

  
//   const int ELEMENT = {
//     /*    
// 	  0xffffffe7,
// 	  0xffffff9f,
// 	  0xfffffe7f,
// 	  0xfffff9ff
//     */

//     0xfffffff3,
//     0xffffffcf,
//     0xffffff3f,
//     0xfffffcff

//   };

//   const int INDEX = {
//     /*
//     0xffffe7ff,
//     0xffff9fff,
//     0xfffe7fff,
//     0xfff9ffff
//     */

//     0xfffff3ff,
//     0xffffcfff,
//     0xffff3fff,
//     0xfffcffff

//   };


//   const int HISTORY = {

//     0xfff3ffff,
//     0xffcfffff,
//     0xff3fffff,
//     0xfcffffff

//   };





//>>>>>>> 72038e9f7f75621d70d570affe4773a0b30a46c0
//   /*
// #define ELT0 0xffffffe7
// #define ELT1 0xffffff9f
// #define ELT2 0xfffffe7f
// #define ELT3 0xfffff9ff
// #define IDX0 0xffffe7ff
// #define IDX1 0xffff9fff
// #define IDX2 0xfffe7fff
// #define IDX3 0xfff9ffff
// */

//   /**********************************************
//    * HexStack
//    **********************************************/
//   /// Sparse set representation

//   class HexStack 
//   {
//   public:

//     /*!@name Parameters*/
//     //@{
//     int _data_;

//     // the first 2 bits stand for the current size (1-4)
//     // the following 8 bits store the values (0-3)
//     // the following 8 bits store the index (0-3)
//     // the following 8 bits store the suucessive sizes (1-4)
//     // the following 4 bits don't do anything
//     // the last 2 bits store the number of changes
//     //@}

//     /*!@name Constructors*/
//     //@{
//     HexStack()
//     {
//       _data_ = 0;
//     }

//     HexStack(const int n=1)
//     {
//       initialise(n);
//     }

//     virtual ~HexStack()
//     {
//     }

//     inline int size() { return (_data_ & 3)+1; }
//     inline int element(const int i) { return (_data_ >> (2+i)) & 3; }
//     inline int index(const int e) { return (_data_ >> (10+e)) & 3; }

//     virtual void initialise(const int n)
//     {
//       _data_ = n-1;

//       for(int i=0, int e=0; e<4; ++i, ++e) 
//   	{
//   	  // store the value
//   	  _data_ |= (e << (i+2));
//   	  // store the index
//   	  _data_ |= (i << (e+10));
//   	}
//     }
//     //@}    

//     /*!@name Accessors*/
//     //@{ 
//     inline bool contain(const int elt) const 
//     {
//       return ((_data_ >> (elt+10)) & 3)<=(_data_ & 3);
//     } 
  
//     inline bool singleton() const 
//     {
//       return !(_data_ & 3);
//     } 

//     inline int next(const int elt) const
//     {
//       unsigned int idx = ((_data_ >> (elt+10)) & 3);
//       return ((_data_ >> (2+idx+(idx < (_data_ & 3)))) & 3);
//     }
    
//     inline int operator[](const unsigned int idx) const
//     {
//       return (_data_ >> (2+idx)) & 3;
//     }
//     //@}

//     /*!@name List Manipulation*/
//     //@{
//     inline void fill()
//     {
//       _data_ &= SIZE;
//       _data_ |= 3;
//     }

//     inline void clear()
//     {
//       _data_ &= SIZE;
//     }
  
//     inline void remove(const int elt)
//     {
//       // get the index of the last element
//       unsigned int last_idx = (_data_ & 3);

//       // decrease the size
//       --_data_;

//       // idx of elt: 
//       unsigned int elt_idx = ((_data_ >> (10+elt)) & 3);

//       // idx of last: 
//       unsigned int last_elt = ((_data_ >> (2+last_idx)) & 3);
      
//       // set the index of the last element to the index of elt
//       _data_ &= INDEX[last_elt];
//       _data_ |= (elt_idx << (10+last_elt));


//       index_[list_[size]] = index_[elt];
//       list_[index_[elt]] = list_[size];
//       list_[size] = elt;
//       index_[elt] = size;
//     }

//     inline int next()
//     {
//       return list_[size];
//     }

//     inline int pop()
//     {
//       return list_[--size];
//     }

//     inline int pop_head()
//     {
//       --size;
//       index_[list_[size]] = 0;
//       const int elt = *list_;
//       *list_ = list_[size];
//       list_[size] = elt;
//       index_[elt] = size;
//       return elt;
//     }

//     inline int head()
//     {
//       return *list_;
//     }
    
//     inline int back()
//     {
//       return list_[size-1];
//     }

//     inline void add(const int elt)
//     {
//       //       std::cout << elt << ", " << size << " <= " << capacity << std::endl; 
//       //       std::cout << index_[elt] << " <= " << capacity << std::endl; 

//       index_[list_[size]] = index_[elt];
//       list_[index_[elt]] = list_[size];
//       list_[size] = elt;
//       index_[elt] = size;
//       ++size;
//     }

//     // create a new element that can potentially be outside the bounds
//     inline void create(const int elt)
//     {
//       extend(elt);
//       add(elt);
//     }

//     inline void ordered_add(const int elt)
//     {
//       // the first non-element goes where elt was
//       index_[list_[size]] = index_[elt];
//       list_[index_[elt]] = list_[size];

//       int idx = size;
//       while(idx && list_[idx-1] > elt) { // push every values greater than elt above elt
//   	list_[idx] = list_[idx-1];
//   	index_[list_[idx-1]] = idx;
//   	--idx;
//       }

//       list_[idx] = elt;
//       index_[elt] = idx;
//       ++size;
//     }


//     inline void revert_to(const int level)
//     {
//       size = level;
//     }

//     inline void index()
//     {
//       for(unsigned int i=0; i<capacity; ++i)
//   	index_[list_[i]] = i;
//     }
//     //@}

//     /*!@name Miscellaneous*/
//     //@{
//     //     std::string getString() const {
//     //       std::string return_str = "(";
//     //       if(size) return_str += toString(list_[0]);
//     //       for(unsigned int i=1; i<size; ++i)
//     // 	return_str += (" "+toString(list_[i]));
//     //       return_str += ")";
      
//     //       return return_str;
//     //     }

//     std::ostream& display(std::ostream& os) const {
//       os << "(";
//       if(size) os << list_[0];
//       for(unsigned int i=1; i<size; ++i)
//   	os << " " << list_[i];
//       os << ")";
//       return os;
//     }
//     //@}
//   };



  // template < int ARITY >
  // class ActiveStack {

  //   int _data_;
  //   //int _level_;

  //   ActiveStack(const int n) {
  //     _data_ = (1 << n)-1;
  //   }

  //   ~ActiveStack() {}

  //   void remove(const int x, const int lvl) 
  //   { 
  //     if(lvl == _level_) {
  // 	_data_ ^= (1 << x);
  //     } else {
  // 	_data_ = (_data_ << ARITY) | (_data_ ^ (1 << x)); }
  //   }


  // };


  // /**********************************************
  //  * Stack
  //  **********************************************/
  // /// Sparse set representation

  // template< class PTR_TYPE >
  // class Stack 
  // {
  // public:

  //   /*!@name Parameters*/
  //   //@{
  //   /// list of values
  //   PTR_TYPE *list_;
  //   /// current max capacity
  //   unsigned int capacity;
  //   /// current size
  //   unsigned int size;
  //   /// values' indices
  //   unsigned int *index_;
  //   int offset;
  //   //unsigned int *start_;
  //   //@}

  //   /*!@name Constructors*/
  //   //@{
  //   Stack()
  //   {
  //     size = 0;
  //     capacity = 0;
  //     list_ = NULL;
  //     offset = 0;
  //     index_ = NULL;
  //   }

  //   Stack(Vector< PTR_TYPE >& obj, bool full=true)
  //   {
  //     initialise(obj, full);
  //   }

  //   virtual ~Stack()
  //   {
  //     delete [] list_;
  //     index_  += offset;
  //     delete [] index_;
  //   }

  //   void initialise(Vector< PTR_TYPE >& obj, const bool full=true)
  //   {
  //     assert((obj.size == 0) || ((unsigned int)(obj.back()->id - obj[0]->id + 1) == obj.size));

  //     capacity = (obj.size ? obj.size : obj.capacity);
  //     list_ = new PTR_TYPE[capacity];
  //     offset = (obj.size ? obj[0]->id : 0);
  //     index_ = new unsigned int[capacity];
  //     index_ -= offset;
  //     for(unsigned int i=0; i<capacity; ++i) 
  // 	{
  // 	  index_[i+offset] = i;
  // 	  list_[i] = obj[i];
  // 	}
      
  //     size = (full ? capacity : 0);
  //   }
  //   //@}    

  //   /*!@name Accessors*/
  //   //@{  
  //   inline bool contain(const PTR_TYPE elt) const 
  //   {
  //     return index_[elt->id]<size;
  //   } 
  //   inline bool contain(const int elt) const 
  //   {
  //     return index_[elt]<size;
  //   } 
  
  //   inline bool empty()const 
  //   {
  //     return !size;
  //   } 

  //   inline PTR_TYPE next(const PTR_TYPE elt) const
  //   {
  //     unsigned int idx = index_[elt->id]+1;
  //     return (idx < size ? list_[idx] : elt);
  //   }
  //   inline PTR_TYPE next(const int elt) const
  //   {
  //     unsigned int idx = index_[elt]+1;
  //     return (idx < size ? list_[idx] : elt);
  //   }
    
  //   inline PTR_TYPE operator[](const unsigned int idx) const
  //   {
  //     return list_[idx];
  //   }

  //   inline PTR_TYPE& operator[](const unsigned int idx)
  //   {
  //     return list_[idx];
  //   }
  //   //@}

  //   /*!@name List Manipulation*/
  //   //@{
  //   inline void fill()
  //   {
  //     size = capacity;
  //   }

  //   inline void clear()
  //   {
  //     size = 0;
  //   }
  
  //   inline void set_to(const PTR_TYPE elt)
  //   {
  //     int idx = elt->id;
  //     size=1;
  //     index_[(*list_)->id] = index_[idx];
  //     list_[index_[idx]] = *list_;
  //     *list_ = elt;
  //     index_[idx] = 0;
  //   }

  //   inline void remove(const PTR_TYPE elt)
  //   {
  //     int idx = elt->id;
  //     --size;
  //     index_[list_[size]->id] = index_[idx];
  //     list_[index_[idx]] = list_[size];
  //     list_[size] = elt;
  //     index_[idx] = size;
  //   }

  //   inline void remove(const int idx)
  //   {
  //     PTR_TYPE elt = list_[index_[idx]];
  //     --size;
  //     index_[list_[size]->id] = index_[idx];
  //     list_[index_[idx]] = list_[size];
  //     list_[size] = elt;
  //     index_[idx] = size;
  //   }

  //   inline PTR_TYPE next()
  //   {
  //     return list_[size];
  //   }

  //   inline PTR_TYPE pop()
  //   {
  //     return list_[--size];
  //   }

  //   inline PTR_TYPE pop_head()
  //   {
  //     --size;
  //     index_[list_[size]->id] = 0;
  //     const PTR_TYPE elt = *list_;
  //     *list_ = list_[size];
  //     list_[size] = elt;
  //     index_[elt->id] = size;
  //     return elt;
  //   }

  //   inline PTR_TYPE head()
  //   {
  //     return *list_;
  //   }
    
  //   inline PTR_TYPE back()
  //   {
  //     return list_[size-1];
  //   }

  //   inline void add(const PTR_TYPE elt)
  //   {
  //     int idx = elt->id;
  //     index_[list_[size]->id] = index_[idx];
  //     list_[index_[idx]] = list_[size];
  //     list_[size] = elt;
  //     index_[idx] = size;
  //     ++size;
  //   }

  //   inline void ordered_add(const PTR_TYPE elt)
  //   {
  //     int idx = elt->id;

  //     // the first non-element goes where elt was
  //     index_[list_[size]->id] = index_[idx];
  //     list_[index_[idx]] = list_[size];

  //     int rank = size;
  //     while(idx && list_[rank-1]->id > elt->id) { // push every values greater than elt above elt
  // 	list_[rank] = list_[rank-1];
  // 	index_[list_[rank-1]->id] = rank;
  // 	--rank;
  //     }

  //     list_[rank] = elt;
  //     index_[idx] = rank;
  //     ++size;
  //   }

  //   inline void revert_to(const int level)
  //   {
  //     size = level;
  //   }

  //   inline void index()
  //   {
  //     for(unsigned int i=0; i<capacity; ++i)
  // 	index_[list_[i]->id] = list_[i]->id;
  //   }
  //   //@}

  //   /*!@name Miscellaneous*/
  //   //@{
  //   //     std::string getString() const {
  //   //       std::string return_str = "(";
  //   //       if(size) return_str += toString(list_[0]);
  //   //       for(unsigned int i=1; i<size; ++i)
  //   // 	return_str += (" "+toString(list_[i]));
  //   //       return_str += ")";
      
  //   //       return return_str;
  //   //     }

  //   std::ostream& display(std::ostream& os) const {
  //     os << "(";
  //     if(size) os << list_[0];
  //     for(unsigned int i=1; i<size; ++i)
  // 	os << " " << list_[i];
  //     os << ")";
  //     return os;
  //   }
  //   //@}
  // };


  /**********************************************
   * VarStack
   **********************************************/
  /// Sparse set representation

  class Solver;
  template< class VAR_TYPE, class SIZE_TYPE >
  class VarStack 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    VAR_TYPE *list_;
    /// current max capacity
    unsigned int capacity;
    /// current size
    SIZE_TYPE size;
    /// values' indices
    unsigned int *index_;
    int offset;
    //unsigned int *start_;
    //@}

    /*!@name Constructors*/
    //@{
    VarStack() : size(0)
    {
      //size = 0;
      capacity = 0;
      list_ = NULL;
      offset = 0;
      index_ = NULL;
    }

    VarStack(Vector< VAR_TYPE >& obj, bool full=true)
    {
      initialise(obj, full);
    }

    virtual ~VarStack()
    {
      delete [] list_;
      index_  += offset;
      delete [] index_;
    }

    void point_to(VarStack< VAR_TYPE, SIZE_TYPE >& vs) {
      capacity = vs.capacity;
      size = vs.size;
      list_ = vs.list_;
      offset = vs.offset;
      index_ = vs.index_;
    }

    void initialise(Solver *s) {
      size.initialise(s);
    }

    void initialise(const int n)
    {
      capacity = n;
      list_ = new VAR_TYPE[capacity];
      VAR_TYPE x;
      std::fill(list_, list_+capacity, x);
      
      offset = 0;
      index_ = new unsigned int[capacity];
      for(unsigned int i=0; i<capacity; ++i) 
	{
	  index_[i] = i;
	}
      index_ -= offset;
      size = 0;
    }

    void extend(const int idx) {
      //int idx = elt.id();
      int new_lb = offset;
      int new_ub = capacity-offset-1;
      
      if(idx < new_lb || idx > new_ub) {

	if(idx < new_lb) new_lb = idx;
	else if(idx > new_ub) new_ub = idx;
	
	unsigned int new_capacity = new_ub-new_lb+1;
	if(new_capacity < 2*capacity) new_capacity = 2*capacity;
	
	unsigned int *aux_index = index_+offset;
	index_ = new unsigned int[new_capacity];
	memcpy(index_, aux_index, capacity*sizeof(unsigned int));
	for(unsigned int i=capacity; i<new_capacity; ++i) 
	  {
	    index_[i] = i;
	  }
	delete [] aux_index;
	
	VAR_TYPE *aux_list = list_;
	list_ = new VAR_TYPE[new_capacity];
	memcpy(list_, aux_list, capacity*sizeof(VAR_TYPE));

	VAR_TYPE x;
	std::fill(list_+capacity, list_+new_capacity, x);

	delete [] aux_list;
	
	index_ -= new_lb;
	capacity = new_capacity;
	offset = new_lb;
      }
    }

    void declare(VAR_TYPE elt) {
      int idx = elt.id();
      extend(idx);
      
      if(idx < offset || idx >= offset+(int)size) {
	list_[index_[idx+offset]] = list_[size];
	list_[size] = elt;
	index_[idx+offset] = size;
	++size;
      } else {
	add(elt);
      }
    }

    void initialise(Vector< VAR_TYPE >& obj, const bool full=true)
    {
      assert((obj.size == 0) || ((unsigned int)(obj.back().id() - obj[0].id() + 1) == obj.size));

      capacity = (obj.size ? obj.size : obj.capacity);
      list_ = new VAR_TYPE[capacity];
      offset = (obj.size ? obj[0].id() : 0);
      index_ = new unsigned int[capacity];
      index_ -= offset;
      for(unsigned int i=0; i<capacity; ++i) 
	{
	  index_[i+offset] = i;
	  list_[i] = obj[i];
	}
      
      size = (full ? capacity : 0);
    }
    //@}    

    /*!@name Accessors*/
    //@{  
    inline bool contain(const VAR_TYPE elt) const 
    {
      return (int)(index_[elt.id()])<(int)size;
    } 
    inline bool safe_contain(const VAR_TYPE elt) const 
    {
      return (elt.id()-offset <= (int)capacity && (int)(index_[elt.id()])<(int)size);
    } 
    inline bool contain(const int elt) const 
    {
      return index_[elt]<(unsigned int)size;
    } 
    inline bool safe_contain(const int elt) const 
    {
      return (elt-offset <= capacity && index_[elt]<(unsigned int)size);
    } 
  
    inline bool empty()const 
    {
      return !(int)size;
    } 

    inline VAR_TYPE next(const VAR_TYPE elt) const
    {
      unsigned int idx = index_[elt.id()]+1;
      return (idx < size ? list_[idx] : elt);
    }
    inline VAR_TYPE next(const int elt) const
    {
      unsigned int idx = index_[elt]+1;
      return (idx < size ? list_[idx] : elt);
    }
    
    inline VAR_TYPE operator[](const unsigned int idx) const
    {
      return list_[idx];
    }

    inline VAR_TYPE& operator[](const unsigned int idx)
    {
      return list_[idx];
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline void fill()
    {
      size = capacity;
    }

    inline void clear()
    {
      size = 0;
    }
  
    inline void set_to(const VAR_TYPE elt)
    {
      int idx = elt.id();
      size=1;
      index_[(*list_).id()] = index_[idx];
      list_[index_[idx]] = *list_;
      *list_ = elt;
      index_[idx] = 0;
    }

    inline void remove(const VAR_TYPE elt)
    {
      int idx = elt.id();
      --size;
      index_[list_[size].id()] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
    }

    inline void remove(const int idx)
    {
      VAR_TYPE elt = list_[index_[idx]];
      --size;
      index_[list_[size].id()] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
    }

    inline VAR_TYPE next()
    {
      return list_[size];
    }

    inline VAR_TYPE pop()
    {
      return list_[--size];
    }

    inline VAR_TYPE pop_head()
    {
      --size;
      index_[list_[size].id()] = 0;
      const VAR_TYPE elt = *list_;
      *list_ = list_[size];
      list_[size] = elt;
      index_[elt.id()] = size;
      return elt;
    }

    inline VAR_TYPE head()
    {
      return *list_;
    }
    
    inline VAR_TYPE back()
    {
      return list_[size-1];
    }

    inline void add(const VAR_TYPE elt)
    {
      
      // std::cout << "add " << elt <<  " to " ;
      // display(std::cout);
      // std::cout << " " << (int*)list_ << " " << index_ << std::endl;

      int idx = elt.id();
      index_[list_[size].id()] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
      ++size;
    }

    inline void ordered_add(const VAR_TYPE elt)
    {
      int idx = elt.id();

      // the first non-element goes where elt was
      index_[list_[size].id()] = index_[idx];
      list_[index_[idx]] = list_[size];

      int rank = size;
      while(idx && list_[rank-1].id() > elt.id()) { // push every values greater than elt above elt
	list_[rank] = list_[rank-1];
	index_[list_[rank-1].id()] = rank;
	--rank;
      }

      list_[rank] = elt;
      index_[idx] = rank;
      ++size;
    }

    inline void revert_to(const int level)
    {
      size = level;
    }

    inline int index(const int elt_idx) {
      if(elt_idx>((int)capacity-offset)) return -1;
      return index_[elt_idx];
    }

    inline void index()
    {
      for(unsigned int i=0; i<capacity; ++i)
	index_[list_[i].id()] = list_[i].id();
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    //     std::string getString() const {
    //       std::string return_str = "(";
    //       if(size) return_str += toString(list_[0]);
    //       for(unsigned int i=1; i<size; ++i)
    // 	return_str += (" "+toString(list_[i]));
    //       return_str += ")";
      
    //       return return_str;
    //     }

    std::ostream& display(std::ostream& os) const {
      os << "(";
      if(size) os << list_[0];
      for(int i=1; i<size; ++i)
	os << " " << list_[i];
      os << ")";
      return os;
    }
    //@}
  };



  /**********************************************
   * ConStack
   **********************************************/
  /// Sparse set representation

  template< class CON_TYPE >
  class ConStack 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    CON_TYPE *list_;
    /// current max capacity
    unsigned int capacity;
    /// current size
    unsigned int size;
    /// values' indices
    unsigned int *index_;
    int offset;
    //unsigned int *start_;
    //@}

    /*!@name Constructors*/
    //@{
    ConStack()
    {
      size = 0;
      capacity = 0;
      list_ = NULL;
      offset = 0;
      index_ = NULL;
    }

    ConStack(Vector< CON_TYPE >& obj, bool full=true)
    {
      initialise(obj, full);
    }

    virtual ~ConStack()
    {
      delete [] list_;
      index_  += offset;
      delete [] index_;
    }

    void initialise(const int n)
    {
      capacity = n;
      list_ = new CON_TYPE[capacity];
      offset = 0;
      index_ = new unsigned int[capacity];
      for(unsigned int i=0; i<capacity; ++i) 
	{
	  index_[i] = i;
	}
      index_ -= offset;
      size = 0;
    }

    void extend(const int idx) {
      //int idx = elt.id();
      int new_lb = offset;
      int new_ub = capacity-offset-1;
      
      if(idx < new_lb || idx > new_ub) {

	if(idx < new_lb) new_lb = idx;
	else if(idx > new_ub) new_ub = idx;
	
	unsigned int new_capacity = new_ub=new_lb+1;
	if(new_capacity < 2*capacity) new_capacity = capacity;
	
	unsigned int *aux_index = index_+offset;
	index_ = new unsigned int[new_capacity];
	memcpy(index_, aux_index, capacity*sizeof(unsigned int));
	delete [] aux_index;
	
	CON_TYPE *aux_list = list_;
	list_ = new CON_TYPE[new_capacity];
	memcpy(list_, aux_list, capacity*sizeof(CON_TYPE));
	delete [] aux_list;

	index_ -= new_lb;
	capacity = new_capacity;
      }
    }

    void declare(CON_TYPE elt) {
      int idx = elt->id;
      extend(idx);
      
      if(idx < offset || idx >= offset+(int)size) {
	list_[index_[idx+offset]] = list_[size];
	list_[size] = elt;
	index_[idx+offset] = size;
	++size;
      } else {
	add(elt);
      }
    }

    void initialise(Vector< CON_TYPE >& obj, const bool full=true)
    {
      assert((obj.size == 0) || ((unsigned int)(obj.back()->id - obj[0]->id + 1) == obj.size));

      capacity = (obj.size ? obj.size : obj.capacity);
      list_ = new CON_TYPE[capacity];
      offset = (obj.size ? obj[0]->id : 0);
      index_ = new unsigned int[capacity];
      index_ -= offset;
      for(unsigned int i=0; i<capacity; ++i) 
	{
	  index_[i+offset] = i;
	  list_[i] = obj[i];
	}
      
      size = (full ? capacity : 0);
    }
    //@}    

    /*!@name Accessors*/
    //@{  
    inline bool contain(const CON_TYPE elt) const 
    {
      return index_[elt->id]<size;
    } 
    inline bool contain(const int elt) const 
    {
      return index_[elt]<size;
    } 
  
    inline bool empty()const 
    {
      return !size;
    } 

    inline CON_TYPE next(const CON_TYPE elt) const
    {
      unsigned int idx = index_[elt->id]+1;
      return (idx < size ? list_[idx] : elt);
    }
    inline CON_TYPE next(const int elt) const
    {
      unsigned int idx = index_[elt]+1;
      return (idx < size ? list_[idx] : elt);
    }
    
    inline CON_TYPE operator[](const unsigned int idx) const
    {
      return list_[idx];
    }

    inline CON_TYPE& operator[](const unsigned int idx)
    {
      return list_[idx];
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline void fill()
    {
      size = capacity;
    }

    inline void clear()
    {
      size = 0;
    }
  
    inline void set_to(const CON_TYPE elt)
    {
      int idx = elt->id;
      size=1;
      index_[(*list_)->id] = index_[idx];
      list_[index_[idx]] = *list_;
      *list_ = elt;
      index_[idx] = 0;
    }

    inline void remove(const CON_TYPE elt)
    {
      int idx = elt->id;
      --size;
      index_[list_[size]->id] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
    }

    inline void remove(const int idx)
    {
      CON_TYPE elt = list_[index_[idx]];
      --size;
      index_[list_[size]->id] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
    }

    inline CON_TYPE next()
    {
      return list_[size];
    }

    inline CON_TYPE pop()
    {
      return list_[--size];
    }

    inline CON_TYPE pop_head()
    {
      --size;
      index_[list_[size]->id] = 0;
      const CON_TYPE elt = *list_;
      *list_ = list_[size];
      list_[size] = elt;
      index_[elt->id] = size;
      return elt;
    }

    inline CON_TYPE head()
    {
      return *list_;
    }
    
    inline CON_TYPE back()
    {
      return list_[size-1];
    }

    inline void add(const CON_TYPE elt)
    {
      int idx = elt->id;
      index_[list_[size]->id] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
      ++size;
    }

    inline void ordered_add(const CON_TYPE elt)
    {
      int idx = elt->id;

      // the first non-element goes where elt was
      index_[list_[size]->id] = index_[idx];
      list_[index_[idx]] = list_[size];

      int rank = size;
      while(idx && list_[rank-1]->id > elt->id) { // push every values greater than elt above elt
	list_[rank] = list_[rank-1];
	index_[list_[rank-1]->id] = rank;
	--rank;
      }

      list_[rank] = elt;
      index_[idx] = rank;
      ++size;
    }

    inline void revert_to(const int level)
    {
      size = level;
    }

    inline void index()
    {
      for(unsigned int i=0; i<capacity; ++i)
	index_[list_[i]->id] = list_[i]->id;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    //     std::string getString() const {
    //       std::string return_str = "(";
    //       if(size) return_str += toString(list_[0]);
    //       for(unsigned int i=1; i<size; ++i)
    // 	return_str += (" "+toString(list_[i]));
    //       return_str += ")";
      
    //       return return_str;
    //     }

    std::ostream& display(std::ostream& os) const {
      os << "(";
      if(size) os << list_[0];
      for(unsigned int i=1; i<size; ++i)
	os << " " << list_[i];
      os << ")";
      return os;
    }
    //@}
  };

 
  // /**********************************************
  //  * IdStack
  //  **********************************************/
  // /// Sparse set representation
  // template< class DATA_TYPE >
  // class IdStack 
  // {
  // public:

  //   /*!@name Parameters*/
  //   //@{
  //   /// list of values
  //   DATA_TYPE *list_;
  //   /// current max capacity
  //   unsigned int capacity;
  //   /// current size
  //   unsigned int size;
  //   //@}

  //   /*!@name Constructors*/
  //   //@{
  //   IdStack()
  //   {
  //     size = 0;
  //     capacity = 0;
  //     list_ = NULL;
  //   }

  //   IdStack(Vector< DATA_TYPE >& obj, bool full=true)
  //   {
  //     initialise(obj, full);
  //   }

  //   virtual ~IdStack()
  //   {
  //     delete [] list_;
  //   }

  //   void initialise(Vector< DATA_TYPE >& obj, const bool full=true)
  //   {
  //     capacity = obj.size;
  //     size = (full ? capacity : 0);
  //     list_ = new DATA_TYPE[capacity];
  //     for(unsigned int i=0; i<capacity; ++i) 
  // 	{
  // 	  list_[i] = obj[i];
  // 	  list_[i].set_stack_id(i);
  // 	}
  //   }
  //   //@}    

  //   /*!@name Accessors*/
  //   //@{  
  //   inline bool contain(const DATA_TYPE elt) const 
  //   {
  //     return elt.get_stack_id()<size;
  //   } 
  
  //   inline bool empty()const 
  //   {
  //     return !size;
  //   } 

  //   inline DATA_TYPE next(const DATA_TYPE elt) const
  //   {
  //     unsigned int idx = elt.get_stack_id()+1;
  //     return (idx < size ? list_[idx] : elt);
  //   }
    
  //   inline DATA_TYPE operator[](const unsigned int idx) const
  //   {
  //     return list_[idx];
  //   }

  //   inline DATA_TYPE& operator[](const unsigned int idx)
  //   {
  //     return list_[idx];
  //   }
  //   //@}

  //   /*!@name List Manipulation*/
  //   //@{
  //   inline void fill()
  //   {
  //     size = capacity;
  //   }

  //   inline void clear()
  //   {
  //     size = 0;
  //   }
  
  //   inline void set_to(DATA_TYPE elt)
  //   {
  //     int idx = elt.get_stack_id();
  //     size=1;
      
  //     (*list_).set_stack_id(idx);
  //     elt.set_stack_id(0);

  //     list_[idx] = *list_;
  //     *list_ = elt;
  //   }

  //   inline void remove(DATA_TYPE elt)
  //   {
  //     int idx = elt.get_stack_id();

  //     list_[--size].set_stack_id(idx);
  //     elt.set_stack_id(size);
      
  //     list_[idx] = list_[size];
  //     list_[size] = elt;
  //   }

  //   inline DATA_TYPE pop()
  //   {
  //     return list_[--size];
  //   }

  //   inline DATA_TYPE pop_head()
  //   {
  //     --size;
  //     const DATA_TYPE elt = *list_;

  //     list_[size].set_stack_id(0);
  //     elt.set_stack_id(size);

  //     *list_ = list_[size];
  //     list_[size] = elt;
  //   }

  //   inline DATA_TYPE head()
  //   {
  //     return *list_;
  //   }
    
  //   inline DATA_TYPE back(const int offset=0)
  //   {
  //     return list_[size-1+offset];
  //   }

  //   inline void add(DATA_TYPE elt)
  //   {
  //     int idx = elt.get_stack_id();
      
  //     elt.set_stack_id(size);
  //     list_[size].set_stack_id(idx);

  //     list_[idx] = list_[size];
  //     list_[size] = elt;
      
  //     ++size;
  //   }

  //   inline void revert_to(const int level)
  //   {
  //     size = level;
  //   }

  //   std::ostream& display(std::ostream& os) const {
  //     os << "(";
  //     if(size) os << list_[0];
  //     for(unsigned int i=1; i<size; ++i)
  // 	os << " " << list_[i];
  //     os << ")";
  //     return os;
  //   }
  //   //@}
  // };


  template< class DATA_TYPE >
  class Indexed 
  {
    
  public:
    
    DATA_TYPE element;
    int index;

    Indexed() {};
    Indexed(DATA_TYPE x) { element = x; index = INFTY; };
    virtual ~Indexed() {};
  };

  /**********************************************
   * Stack
   **********************************************/
  /// Sparse set representation
  template< class DATA_TYPE >
  class Stack 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    Indexed< DATA_TYPE > *list_;
    /// current max capacity
    unsigned int capacity;
    /// current size if full
    unsigned int fsize;
    /// current size
    unsigned int size;
    //@}

    /*!@name Constructors*/
    //@{
    Stack()
    {
      size = 0;
      fsize = 0;
      capacity = 0;
      list_ = NULL;
    }

    Stack(Vector< DATA_TYPE >& obj, bool full=true)
    {
      initialise(obj, full);
    }

    virtual ~Stack()
    {
      delete [] list_;
    }

    void extend()
    {
      capacity = 2*fsize;
      Indexed<DATA_TYPE> *new_list_ = new Indexed<DATA_TYPE>[capacity];
      memcpy(new_list_, list_, sizeof(Indexed<DATA_TYPE>)*fsize);
      delete [] list_;
      list_ = new_list_;
    }

    void initialise(Vector< DATA_TYPE >& obj, const bool full=true)
    {
      fsize = obj.size;
      size = (full ? fsize : 0);
      capacity = 2*fsize;
      list_ = new Indexed<DATA_TYPE>[capacity];
      for(unsigned int i=0; i<fsize; ++i) 
	{
	  list_[i] = Indexed<DATA_TYPE>(obj[i]);
	  list_[i].index = i;
	}
    }
    //@}    

    /*!@name Accessors*/
    //@{  
    inline bool contain(const Indexed<DATA_TYPE> elt) const 
    {
      return elt.index<size;
    } 
  
    inline bool empty()const 
    {
      return !size;
    } 

    inline DATA_TYPE next(const Indexed<DATA_TYPE> elt) const
    {
      unsigned int idx = elt.index+1;
      return (idx < size ? list_[idx] : elt);
    }
    
    inline DATA_TYPE operator[](const unsigned int idx) const
    {
      return list_[idx].element;
    }

    inline DATA_TYPE& operator[](const unsigned int idx)
    {
      return list_[idx].element;
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline void fill()
    {
      size = capacity;
    }

    inline void clear()
    {
      size = 0;
    }
  
    inline void set_to(Indexed<DATA_TYPE> elt)
    {
      int idx = elt.index;
      size=1;
      
      list_->index = idx;
      elt.index = 0;

      list_[idx] = *list_;
      *list_ = elt;
    }

    inline void remove(Indexed<DATA_TYPE> elt)
    {
      int idx = elt.index;

      list_[--size].index = idx;
      elt.index = size;
      
      list_[idx] = list_[size];
      list_[size] = elt;
    }

    inline DATA_TYPE pop()
    {
      return list_[--size].element;
    }

    inline DATA_TYPE pop_head()
    {
      --size;
      const DATA_TYPE elt = list_->element;

      list_[size].index = 0;
      elt.index = size;

      *list_ = list_[size];
      list_[size] = elt;
    }

    inline DATA_TYPE head()
    {
      return list_->element;
    }
    
    inline DATA_TYPE back(const int offset=0)
    {
      return list_[size-1+offset].element;
    }

    inline void add(Indexed<DATA_TYPE> elt)
    {
      int idx = elt.index;
      
      elt.index = size;
      list_[size].index = idx;

      list_[idx] = list_[size];
      list_[size] = elt;
      
      ++size;
    }

    inline void declare(DATA_TYPE elt)
    {
      Indexed<DATA_TYPE> x(elt, fsize);
      if(fsize == capacity) extend();
      list_[fsize++] = x;
      
      add(x);
    }

    inline void revert_to(const int level)
    {
      size = level;
    }

    std::ostream& display(std::ostream& os) const {
      os << "(";
      if(size) os << list_[0].element;
      for(unsigned int i=1; i<size; ++i)
	os << " " << list_[i].element;
      os << ")";
      return os;
    }
    //@}
  };


  /**********************************************
   * Node
   **********************************************/
  /// Node of Multilist

  template < class DATA_TYPE >
  class Node {
  public:

    /*!@name Parameters*/
    //@{
    int prev;
    int next;
    DATA_TYPE elt;
    //@}

    Node(DATA_TYPE e=(DATA_TYPE(0)), int p=-1, int n=-1) {elt = e; prev=p; next=n;}
    inline operator const DATA_TYPE() const { return elt; }

    //inline Node& operator=(const Node& nd) const { prev = nd.prev; return elt; }

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    //     std::string getString() const {
    //       std::string return_str = "("+toString(prev)+"|"+toString(elt)+"|"+toString(next)+")";
    //       return return_str;
    //     }
    std::ostream& display(std::ostream& os) const {
      //os << "(" << prev << "|" << elt << "|" << next << ")";
      os << elt;
      return os;
    }
    //@}
  };


  /**********************************************
   * Multilist
   **********************************************/
  /// List with multiple entry points

  template < class DATA_TYPE, int NUM_HEAD >
  class MultiList 
  {

  public:
    /*!@name Parameters*/
    //@{
    Vector< Node< DATA_TYPE > > data;
    int head[NUM_HEAD+1]; // index of the first element of the list
    unsigned int degree; // number of element in the entire list
    //@}

    /*!@name Constructors*/
    //@{
    MultiList() {
      data.initialise(0, std::max(2*NUM_HEAD,16));
      for(int k=0; k<=NUM_HEAD; ++k) {
	head[k] = k; 
	Node< DATA_TYPE > x((DATA_TYPE)0, k-1, k+1); 
	data.add(x);
      }
      degree=0;
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline int create(DATA_TYPE elt, const int k=0) {
      Node<DATA_TYPE> x(elt);
      data.add(x);
      add(data.size-1, k);
      return data.size-1;
    }

    inline void add(const int idx, const int k=0) {
      int succ = data.stack_[head[k]].next;
      data.stack_[succ].prev = idx;
      data.stack_[head[k]].next = idx;
      data.stack_[idx].next = succ;
      data.stack_[idx].prev = head[k];
      ++degree;
    }

    inline void remove(const int idx, const int k=0) {
      int succ = data.stack_[idx].next;
      int prec = data.stack_[idx].prev;   

      data.stack_[succ].prev = prec;
      data.stack_[prec].next = succ;

      --degree;
    }

    inline Node<DATA_TYPE>& first(const int h) {
      return data.stack_[head[h]];
    }

    inline Node<DATA_TYPE>& last(const int h) {
      return data.stack_[data.stack_[head[h+1]].prev];
    }

    inline DATA_TYPE pop(const int h) {
      int idx = data.stack_[head[h+1]].prev;
      remove(idx, h);
      return data.stack_[idx].elt;
    }

    inline bool empty(const int h) {
      return (data.stack_[head[h]].next > NUM_HEAD);
    }

    inline bool next(Node< DATA_TYPE >& node) const {
      unsigned int idx_next;
      do { 
	idx_next = node.next;
	if(idx_next == NUM_HEAD) return false;
	node = data.stack_[idx_next];
      } while(!node.elt);
      return true;
    }

    inline bool prev(Node< DATA_TYPE >& node) const {
      unsigned int idx_prev;
      do { 
	idx_prev = node.prev;
	if(!idx_prev) return false;
	node = data.stack_[idx_prev];
      } while(!node.elt);
      return true;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    //     std::string getString() const {
    //       std::string return_str = "";
    //       for(unsigned int i=0; i<data.size; ++i)
    // 	return_str += toString(data[i]);
    //       return return_str;
    //     }
    std::ostream& display(std::ostream& os) const;//  {

    void debug_print(std::ostream& os) const 
    {
      for(unsigned int i=head[0]; i<data.size; ++i)
	os << "(" << data.stack_[i].prev << "|" 
	   << data.stack_[i].elt << "|" 
	   << data.stack_[i].next << ")";
      //o << data.stack_[i].print(o);
    }
    void flat_print(std::ostream& o) const 
    {
      o << "[ ";
      for(unsigned int i=NUM_HEAD+1; i<data.size; ++i)
	o << data.stack_[i].elt << " ";
      o << "] ";
    }
    void print(std::ostream& o) const 
    {
      for(int k=0; k<NUM_HEAD; ++k) {
	o << "[ ";

	Node< DATA_TYPE > nd = data.stack_[head[k]];
	while(next(nd)) {
	  o << (DATA_TYPE)nd << " ";
	}
	o << "] ";
      }
    }
    //@}

  };


  /**********************************************
   * BitSet
   **********************************************/

  /*! \class BitSet
    \brief A representation of sets using a vector of bits.

    The sets have a static capacity. 
    Template are used so that ReversibleWords can used instead of unsigned int
  */
  template< class WORD_TYPE, class FLOAT_TYPE >
  class Bitset
  { 

  public:
    
    /*!@name Class attributes*/
    //@{
    static const WORD_TYPE empt = 0;
    static const WORD_TYPE full = ~0;
    static const unsigned int EXP = (sizeof(empt) == 4 ? 5 /*32 bits*/ : 6 /*64 bits*/);
    static const unsigned int size_word_bit = (1 << EXP);
    static const unsigned int size_word_byte = (size_word_bit >> 3);
    static const unsigned int CACHE = (size_word_bit - 1);
    static const unsigned int LASTCHAR = (size_word_bit - 8);
    static const unsigned int mantissa = (sizeof(empt) == 4 ? 23 /*32 bits*/ : 52 /*64 bits*/);
    static const unsigned int float_offset = (sizeof(empt) == 4 ? (0x7f) /*32 bits*/ : (0x3ff) /*64 bits*/);
    static const WORD_TYPE mask_first_char = 0xff;
    static const WORD_TYPE mask_last_char = (mask_first_char << ((size_word_bit) - 8));
    //@}

    /*!@name Parameters*/
    //@{
    /// index of the first word used to represent the set
    int pos_words;
    /// 1 + index of the last word used to represent the set
    int neg_words;
    /// A vector of bits 
    WORD_TYPE* table;
    //@}

    Bitset()
    {
      initialise();
    }

    Bitset(int sz)
    {
      if(sz>0) {
	initialise(sz,0);
      } else {
      	initialise();
      }
    }

    void initialise()
    {
      pos_words = 0;
      neg_words = 0;
      table = NULL;
    }

    Bitset(const int sz, const int* elt) 
    {
      int lb =  NOVAL;
      int ub = -NOVAL;
      for(int i=0; i<sz; ++i) {
	if(elt[i] > ub) ub = elt[i];
	if(elt[i] < lb) lb = elt[i];
      }

      initialise(lb,ub,empt);

      for(int i=0; i<sz; ++i) 
	add( elt[i] );
    }

    void initialise(const Vector< int >& elt) 
    {
      int min = elt.front();
      // int max = elt.front();
      // for(unsigned int i=1; i<elt.size; ++i) {
      // 	if(elt[i] < min) min = elt[i];
      // 	if(elt[i] > max) max = elt[i];
      // }
      int max = elt.back();

      initialise(min,max,empt);
      
      for(unsigned int i=0; i<elt.size; ++i) 
	add( elt[i] );
    }

    void initialise(const int lb, const int ub, const Vector< int >& elt, int&sz) 
    {

      // std::cout << "initialise bitset with " << lb << ".." << ub << ": " ;
      // elt.display(std::cout);
      // std::cout << std::endl; 

      initialise(lb,ub,empt);
			sz = 0;

      // display(std::cout);
      // std::cout << std::endl;

      for(unsigned int i=0; i<elt.size; ++i) {
				if(!contain(elt[i])) {
					add( elt[i] );
					++sz;
				}
	// display(std::cout);
	// std::cout << std::endl;
      }
      
    }

    Bitset(const int lb, const int ub, const WORD_TYPE p)
    {
      initialise(lb,ub,p,NULL);
    }

    inline int word_index(const int elt) const
    {
      return (elt >> EXP);
    }

		

    bool operator==(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) {
      return equal(s);
    }

    bool operator!=(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) {
      return !equal(s);
    }

    Bitset<WORD_TYPE,FLOAT_TYPE>& operator=(const Bitset<WORD_TYPE,FLOAT_TYPE>& q) 
    {
      if(!table)
	clone(q);
      else
	copy(q);
      return *this;
    }

    void reinitialise(const int lb, const int ub, const WORD_TYPE p) 
    {
      table += neg_words;
      delete [] table;
      initialise(lb, ub, p, NULL);
    }

    void initialise(const int sz, const WORD_TYPE p) 
    {
      pos_words = sz;
      neg_words = 0;

      if( sz>=0 ) {
	table = new WORD_TYPE[pos_words];
	for(int i=0; i<pos_words; ++i) 
	  table[i]=p;
      } else table = NULL;
    }

    void initialise(const int lb, const int ub, const WORD_TYPE p, WORD_TYPE *pool=NULL) 
    {
      neg_words = (lb >> EXP);
      pos_words = (ub >> EXP)+1;
      if(pool==NULL) table = new WORD_TYPE[pos_words-neg_words];
      else table = pool;
      for(int i=0; i<pos_words-neg_words; ++i) 
	table[i]=p;
      table[pos_words-neg_words-1] &= 
	(p >> (size_word_bit-1-(ub & CACHE)));
      table[0] &= (p << (lb & CACHE));
      table -= neg_words;
    }

    void initialise(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      pos_words = s.pos_words;
      neg_words = s.neg_words;
  
      table = new WORD_TYPE[pos_words-neg_words];
      table -= neg_words;
      for(int i=neg_words; i<pos_words; ++i) 
	//table[i].initialise(s.table+i, s.size(i));
	table[i] = s.table[i];
    }

    inline void declare(const int elt)
    {
      int i = (elt >> EXP);
      if( (i < neg_words) ||
	  (i >=  pos_words) ) 
	{
	  extend(elt);
	}
      fast_add(elt);
    }

    void extend(const int elt) 
    {
      int nval = (elt >> EXP);
      if( (nval < neg_words) ||
	  (nval >=  pos_words) ) 
	{
	  int new_neg_words = neg_words;
	  //nval;
	  int new_pos_words = pos_words;
	  //nval+1;
	  bool need_to_extend = false;
	  if(nval < new_neg_words) {
	    new_neg_words = nval;
	    need_to_extend = true;
	  }
	  if(nval >= new_pos_words) {
	    new_pos_words = nval+1;
	    need_to_extend = true;
	  }

	  if(need_to_extend) {
	    WORD_TYPE *aux = table;
	    table = new WORD_TYPE[new_pos_words-new_neg_words];
	    table -= new_neg_words;
	    
	    memcpy(table+neg_words, aux+neg_words, 
		   (pos_words-neg_words)*sizeof(WORD_TYPE));

	    if(new_neg_words < neg_words)
	      std::fill(table+new_neg_words, table+neg_words, 0);

	    if(new_pos_words > pos_words)
	      std::fill(table+pos_words, table+new_pos_words, 0);
	    
	    aux += neg_words;
	    delete [] aux;
	    
	    pos_words = new_pos_words; 
	    neg_words = new_neg_words; 
	  }
	}
    }

    Bitset(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {      
      initialise();
      clone( s );
    }

    void clone(const Bitset<WORD_TYPE,FLOAT_TYPE>& s)
    {
      if(table) {
	table += neg_words;
	delete [] table;
      }
      neg_words = s.neg_words;
      pos_words = s.pos_words;
      table = new WORD_TYPE[pos_words-neg_words];
      memcpy(table, s.table+neg_words,
	     size_word_byte*(pos_words-neg_words));
      table -= neg_words;
    }

    void point_to(Bitset<WORD_TYPE,FLOAT_TYPE>& s)
    {
      neg_words = s.neg_words;
      pos_words = s.pos_words;
      table = s.table;
    }

    void point_to(WORD_TYPE *t)
    {
      neg_words = 0;
      pos_words = 1;
      table = t;
    }

		void copy(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
		{
			int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
			int k, j = (neg_words < s.neg_words ? s.neg_words : neg_words);
			for( k=neg_words; k<j; ++k)
				table[k] = empt;
			for( k=i; k<pos_words; ++k)
				table[k] = empt;
			if( i>j )
				memcpy(table+j,s.table+j,size_word_byte*(i-j));
		}

    virtual ~Bitset() 
    {
      table += neg_words;
			// if(table)
			delete [] table; 
    }

    void destroy() 
    {		
      table += neg_words;
      neg_words = 0;
      delete [] table; 
      table = NULL;
    }

    bool is_built()
    {
      return (table != NULL);
    }

    inline void swap(Bitset<WORD_TYPE,FLOAT_TYPE>& s)
    {	
      WORD_TYPE *aux = s.table;
      s.table = table;
      table = aux;
    }


    void iterate_into_b(const int size, int *buffer) {
      int elt;
      int idx;
      int nval;

      union {FLOAT_TYPE f; WORD_TYPE i; } t;
      WORD_TYPE b;
      WORD_TYPE v; 

      nval = 1;
      elt = buffer[0];
      idx = ((elt+1) >> EXP);
      v = (table[idx] & (full << ((elt+1) & CACHE)));  
      
      while(nval < size) {
	// find the next word that is not null
	while(!v) v = table[++idx];

	// find the first element in the set:
	// remove all other element
	b = v & -v;

	// cast into float, which will be coded as 1*2^exp, and 'exp' is precisely the index of the first element 
	t.f = (FLOAT_TYPE)b; // cast the least significant bit in v to a float
	// keep only the exponant part
	elt = t.i >> mantissa;

	elt += idx * size_word_bit - float_offset;

	do {
	  // put the element in the buffer
	  buffer[nval] = elt;
	  ++nval;
	
	  // remove it from v
	  v ^= b;
	  
	  do {
	    
	    // try the next element
	    b <<= 1;
	    ++elt;

	  } while( b && !(v & b) );

	} while(v);
	
      }
      
    }


    void iterate_into(const int size, int *buffer) {
      int elt;
      int idx;
      int nval;

      union {FLOAT_TYPE f; WORD_TYPE i; } t;
      WORD_TYPE b;
      WORD_TYPE v; 

      nval = 1;
      elt = buffer[0];
      idx = ((elt+1) >> EXP);
      v = (table[idx] & (full << ((elt+1) & CACHE)));  
      
      elt = (idx * size_word_bit - float_offset);

      while(nval < size) {
	// find the next word that is not null
	while(!v) {
	  v = table[++idx];
	  elt += size_word_bit;
	}

	// find the first element in the set:
	do {
	  // remove all other element
	  b = v & -v;
	  
	  // cast into float, which will be coded as 1*2^exp, and 'exp' is precisely the index of the first element 
	  t.f = (FLOAT_TYPE)b; // cast the least significant bit in v to a float
	  	  	  
	  // put the element in the buffer
	  buffer[nval++] = (t.i >> mantissa)+elt;
	  
	  // remove it from v
	  v ^= b;
	  
	} while(v);      
      }
    }


    inline int lsb_mantissa(const WORD_TYPE v) const {
      union {FLOAT_TYPE f; WORD_TYPE i; } t;
      WORD_TYPE b = v & -v;

      t.f = (FLOAT_TYPE)b; // cast the least significant bit in v to a float
      b = t.i >> mantissa;
      return b - float_offset;
    }

#ifdef _BIT64
    inline int lsb_gcc(const WORD_TYPE v) const { return __builtin_ctzl(v); }
    inline int msb_gcc(const WORD_TYPE v) const { return __builtin_clzl(v); }
#else
    inline int lsb_gcc(const WORD_TYPE v) const { return __builtin_ctz(v); }
    inline int msb_gcc(const WORD_TYPE v) const { return __builtin_clz(v); }
#endif

      
		inline int minimum_element(int idx, WORD_TYPE v, const int def=NOVAL) const
		{
           
			while(v == 0) {
				if( ++idx >= pos_words )
					return def;
				v = table[idx];
			}

      
#ifdef _VALGRIND_
      // CODE THAT PASSES VALGRIND
			union {FLOAT_TYPE f; WORD_TYPE i; } t;
			WORD_TYPE b = v & -v;

			t.f = (FLOAT_TYPE)b; // cast the least significant bit in v to a float
			b = t.i >> mantissa;
			
			return b + idx * size_word_bit - float_offset;

#else

			// CODE THAT DOES NOT
			return __builtin_ctz(v) //__builtin_ffs(v) - 1
				+ (idx * size_word_bit);

#endif
      
		}

    /*!
      Minimum element in the set [O(N/32)]. 
    */
    inline int min() const
    { 
      int idx = neg_words;
      WORD_TYPE v = table[idx];
      return minimum_element(idx,v);
    }

    /*!
      Maximum element in the set [O(N/8)]. 
    */
		inline int max() const
		{ 
			WORD_TYPE tab;
			int i=pos_words, j, k;
    
			while( i-- > neg_words )
			if( (tab = table[i]) ) {
				j = size_word_byte;
				while( j-- ) {
					if( (k = getlast[(tab & mask_last_char) >> LASTCHAR]) >= 0 ) 
						return ( (i<<EXP)+(j<<3)+k );	
					tab = (tab << 8);
				}
			}
			return NOVAL;
		}

		inline void  remove(const int elt) 
		{
			int i = (elt >> EXP);
			if( (i >= neg_words) && 
				(i <  pos_words) )
					table[i] &= (full ^ ((WORD_TYPE)1 << (elt & CACHE)));
		}

    inline void  fast_remove(const int elt) 
    {
      table[(elt >> EXP)] &= (full ^ ((WORD_TYPE)1 << (elt & CACHE)));
    }

    inline void  word_remove(const int elt) 
    {
      table[neg_words] &= (full ^ ((WORD_TYPE)1 << (elt & CACHE)));
    }


		inline int next(const int elt) const {
			int idx = ((elt+1) >> EXP);
			if(idx >= pos_words) return elt;
			if(idx < neg_words) return min();
			WORD_TYPE v = (table[idx] & (full << ((elt+1) & CACHE)));
			while(v == 0) {
				if(++idx >= pos_words) return elt;
				v = table[idx];
			}
      
#ifdef _VALGRIND_
			return lsb_mantissa(v) + (idx * size_word_bit);
#else
			return lsb_gcc(v) + (idx * size_word_bit);
#endif
		}


    /*
    inline int next(const int elt) const {
      int idx = ((elt+1) >> EXP);
      if(idx >= pos_words) return elt;
      WORD_TYPE v = (table[idx] & (full << ((elt+1) & CACHE)));
      return minimum_element(idx,v,elt);
    }
    */

    inline int prev(const int elt) const {

			WORD_TYPE tab;
			int i = ((elt-1) >> EXP);
			if(i>=pos_words) return max();
			
			int SHFT = size_word_byte;

			if( i >= neg_words ) {
				int e = ((elt-1) & CACHE), k;
				int j = 1+(e >> 3);

				if( (tab = ((table[i] & (full >> (CACHE - e))) << ((SHFT-j) << 3))) ) 
				while( j-- ) {
					if( (k = getlast[(tab & mask_last_char) >> LASTCHAR]) >= 0 )
						return ( (i<<EXP)+(j<<3)+k );
					tab = (tab << 8);
				}
				while( i-- > neg_words ) 
				if( (tab = table[i]) ) {
					j = size_word_byte;
					while( j-- ) {
						if( (k = getlast[(tab & mask_last_char) >> LASTCHAR]) >= 0 )
							return ( (i<<EXP)+(j<<3)+k );
						tab = (tab << 8);
					}
				}
			}

			return elt;
		}

		inline void xor_to(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
		{
			int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
			int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
			while( i-- > j )
				s.table[i] ^= table[i];
		}

		inline void fast_xor_to(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
		{
			int i = pos_words;
			while( i-- > neg_words )
				s.table[i] ^= table[i];
		}

		inline void xor_with(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
		{
			int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
			int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
			while( i-- > j )
				table[i] ^= s.table[i];
		}

		inline void fast_xor_with(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
		{
			int i = pos_words;
			while( i-- > neg_words )
				table[i] ^= s.table[i];
		}
		
		inline void xor_with(const int lb, const int ub) 
		{
			int wlb = (lb >> EXP), i;
			int wub = (ub >> EXP), j;
			
			if(wlb<neg_words) wlb=neg_words;
			if(wub>=pos_words) wub=pos_words-1;
			
			if(wlb==wub) {
				table[wlb] ^= ((full << (lb & CACHE)) & (full >> (size_word_bit-1-(ub & CACHE))));
			} else {
				++wub;
				
				if(wub>neg_words && wlb<pos_words) {
					if(wub<=pos_words) {
						i = wub-1;
	


						// std::cout << std::endl;
						// showUint(table[wub-1], std::cout);
						// std::cout << std::endl << " xor " << std::endl ;
						// showUint((full >> (size_word_bit-1-(ub & CACHE))), std::cout);
						// std::cout << std::endl << " = " << std::endl;
						
						table[i] ^= (full >> (size_word_bit-1-(ub & CACHE)));
						
						// showUint(table[wub-1], std::cout);
						// std::cout << std::endl ;

					
						
						if(wlb>=neg_words) {
							table[wlb] ^= (full << (lb & CACHE));
							j = wlb+1;
						} else {
							j = wlb;
						}
					} else {
						i = wub;
						if(wlb>=neg_words) {
							table[wlb] ^= (full << (lb & CACHE));
							j = wlb+1;
						} else {
							j = wlb;
						}
					}
					while( i-- > j )
						table[i] ^= full;
				}
			}
		}

    inline void union_with(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	table[i] |= s.table[i];
    }

    inline void union_with(const int s) 
    {
      if(pos_words>0 && neg_words<=0) table[0] |= s;
    }
		
		inline void union_with(const int lb, const int ub) 
		{
			int wlb = (lb >> EXP), i;
			int wub = (ub >> EXP), j;
			
			if(wlb==wub) {
				table[wlb] |= ((full << (lb & CACHE)) & (full >> (size_word_bit-1-(ub & CACHE))));
			} else {
				++wub;
				
				if(wub<=pos_words) {
					i = wub-1;
					table[i] |= (full >> (size_word_bit-1-(ub & CACHE)));
					if(wlb>=neg_words) {
						table[wlb] |= (full << (lb & CACHE));
						j = wlb+1;
					} else {
						j = wlb;
					}
				} else {
					i = wub;
					if(wlb>=neg_words) {
						table[wlb] |= (full << (lb & CACHE));
						j = wlb+1;
					} else {
						j = wlb;
					}
				}
				while( i-- > j )
					table[i] = full;
			}
		}

    inline void union_to(Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
    {
      s.union_with( *this );
    }

		inline void intersect_with(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
		{
			int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
			int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
			int k = pos_words;
			while( k > i ) {
				--k;
				table[k] = empt;
			}
			while( k > j ) {
				--k;
				table[k] &= s.table[k];
			}
			while( k-- > neg_words )
				table[k] = empt;
		}
		
		inline bool intersect_with_check(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
		{
			bool not_empty = false;
			int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
			int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
			int k = pos_words;
			while( k > i ) {
				--k;
				table[k] = empt;
			}
			while( k > j ) {
				--k;
				table[k] &= s.table[k];
				not_empty |= table[k];
			}
			while( k-- > neg_words )
				table[k] = empt;
			
			return not_empty;
		}

    inline void intersect_with(const int s) 
    {
      int i=pos_words;
      while(i-- > 1) table[i] = empt;
      i = 0;
      while(i-- > neg_words) table[i] = empt;
      if(pos_words>0 && neg_words<=0) table[0] &= s;
    }

    inline void intersect_to(Bitset<WORD_TYPE,FLOAT_TYPE>& s) const
    {
      s.intersect_with( *this );
    }

    inline void setminus_with (const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	table[i] &= (~(s.table[i]));
    }

    inline void setminus_to (Bitset<WORD_TYPE,FLOAT_TYPE>& s) const
    {
      s.setminus_with( *this );
    }

    inline void xor_with(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) 
    {
      xor_with(*s);
    }

    inline void union_with  (const Bitset<WORD_TYPE,FLOAT_TYPE>* s) 
    {
      union_with(*s);
    }

    inline void intersect_with(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) 
    {
      intersect_with(*s);
    }

    inline void setminus_with (const Bitset<WORD_TYPE,FLOAT_TYPE>* s) 
    {
      setminus_with(*s);
    }
    inline void union_to  (Bitset<WORD_TYPE,FLOAT_TYPE>* s) const
    {
      s->union_with(*this);
    }

    inline void intersect_to(Bitset<WORD_TYPE,FLOAT_TYPE>* s) const
    {
      s->intersect_with(*this);
    }

    inline void setminus_to (Bitset<WORD_TYPE,FLOAT_TYPE>* s) const
    {
      s->setminus_with(*this);
    }

		inline bool equal(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
		{
			int i=pos_words;
			int j=s.pos_words;
			int k;
  
			while( j > i )
				if( s.table[--j] ) return false; 
			while( i > j )
				if( table[--i] ) return false;

			j=neg_words;
			k=s.neg_words;

			while( j > k )
				if( s.table[k++] ) return false;
			while( k > j )
				if( table[j++] ) return false;

			while( i-- > j )
				if( table[i] != s.table[i] ) return false;
  
			return true;
		}
    
		inline bool includes(const WORD_TYPE s) const 
		{
			return( pos_words && neg_words<1 && (table[0] & s) == s );
		}

		inline bool included(const WORD_TYPE s) const 
		{
			bool inc = true;
			int k = pos_words;
			if(neg_words>0 || pos_words<1) {
				while(k>neg_words && inc) inc = !(table[--k]);
			} else {
				while(k>1 && inc) inc = !(table[--k]);
				inc = ((table[--k] & s) == table[0]);
				while(k>neg_words && inc) inc = !(table[--k]);
			}
			return inc;
		}

		inline bool included(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
		{
			int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
			int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
			int k = pos_words;
			while( k > i ) {
				--k;
				if( table[k] ) return false;
			}
			while( k > j ) {
				--k;
				if( table[k] != (table[k] & s.table[k]) ) return false;
			}
			while( k-- > neg_words )
				if( table[k] ) return false;
			return true;
		}

		inline bool includes(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
		{
			int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
			int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
			int k = s.pos_words;
			while( k > i ) {
				--k;
				if( s.table[k] ) return false;
			}
			while( k > j ) {
				--k;
				if( s.table[k] != (table[k] & s.table[k]) ) {
					return false;
				}
			}
			while( k-- > s.neg_words ) {
				if( s.table[k] ) return false;
			}

			return true;
		}

    inline bool included(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) const 
    {
      return included( *s );
    }

    inline bool includes(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) const 
    {
      return includes( *s );
    }

    inline bool intersect(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) const 
    {
      return intersect( *s );
    }

		inline bool included(const int lb, const int ub) const 
		{
			int neg_int = lb >> EXP;
			int pos_int = ub >> EXP;
			int k = pos_words;
			while( k > pos_int ) 
				if( table[--k] ) return false;
			k = neg_words;
			while( k < neg_int )
				if( table[k++] ) return false;
			if(neg_int == pos_int) {
				k = ((full << (lb & CACHE)) & (full >> (CACHE - (ub & CACHE))));
				return (k & table[neg_int]) == table[neg_int];
			} else {
				return ((((full << (lb & CACHE)) & table[neg_int]) == table[neg_int]) &&
					(((full >> (CACHE - (ub & CACHE))) & table[pos_int]) == table[pos_int]));
			}
		}

		inline bool includes(const int lb, const int ub) const
		{
      int neg_int = lb >> EXP;
      int pos_int = ub >> EXP;

			if(neg_int < neg_words || pos_int >= pos_words) 
				return false;
			
      int k = pos_int - 1;
      unsigned int u, l;
      while (k > neg_int) {
          if (table[k] != full)
              return false;
          --k;
      }
      if (neg_int == pos_int) {
          u = ((full << (lb & CACHE)) & (full >> (CACHE - (ub & CACHE))));
          return (u & table[neg_int]) == u;
      } else {
          u = (full >> (CACHE - (ub & CACHE)));
          l = (full << (lb & CACHE));
          return (((l & table[neg_int]) == l) && ((u & table[pos_int]) == u));
      }
		}



    /*!
     * Returns the number of bits set in v.
     * For a derivation of this algorithm, see
     * "Algorithms and data structures with applications to 
     *  graphics and geometry", by Jurg Nievergelt and Klaus Hinrichs,
     *  Prentice Hall, 1993.
     */
    inline unsigned int word_size(WORD_TYPE v) const 
    {
#ifdef _BIT32
        return __builtin_popcount(v);
#else
				return __builtin_popcountl(v);
#endif
      /*
      v = v - ((v >> 1) & (WORD_TYPE)~(WORD_TYPE)0/3);                           // temp
      v = (v & (WORD_TYPE)~(WORD_TYPE)0/15*3) + ((v >> 2) & (WORD_TYPE)~(WORD_TYPE)0/15*3);      // temp
      v = (v + (v >> 4)) & (WORD_TYPE)~(WORD_TYPE)0/255*15;                      // temp
      return (WORD_TYPE)(v * ((WORD_TYPE)~(WORD_TYPE)0/255)) >> (sizeof(v) - 1) * CHAR_BIT; // count
      */
    }
	
		inline unsigned int size() const 
		{  
			int i=pos_words;
			unsigned int c=0;
			WORD_TYPE v;
			while( i-- > neg_words ) 
				if( (v = table[i]) ) 
					c += word_size(v);
			return c;  
		}

		inline unsigned int word_size() const 
		{  
			unsigned int v, c=0;
			if( (v = table[neg_words]) ) 
				c = word_size(v);
			return c;  
		}

		inline unsigned int size( const int i ) const
		{  
			WORD_TYPE v;
			unsigned int c=0;
			if( (v = table[i]) ) 
				c = word_size(v);
			return c;  
		}



    /*!
      Check if element elt belong to the set [O(1)]
    */
    inline  bool contain(const int elt)const 
    {
      int i = (elt >> EXP);
      return ( (i >= neg_words) && 
	       (i <  pos_words) && 
	       (table[i] & ((WORD_TYPE)1 << (elt & CACHE))) );
    }

    inline  bool fast_contain(const int elt)const 
    {
      return ( (table[(elt >> EXP)] & ((WORD_TYPE)1 << (elt & CACHE))) );
    }

    inline  bool word_contain(const int elt)const 
    {
      return ( (table[neg_words] & ((WORD_TYPE)1 << (elt & CACHE))) );
    }
    /*!
      Add element elt into the set [O(1)]
    */

    inline  void add(const int elt)
    {
      int i = (elt >> EXP);
      if( (i >= neg_words) && 
	  (i <  pos_words) ) 
	table[i] |= ((WORD_TYPE)1 << (elt & CACHE));
    }

    inline  void fast_add(const int elt)
    {
      table[(elt >> EXP)] |= ((WORD_TYPE)1 << (elt & CACHE));
    }

    inline   void word_add(const int elt)
    {
      table[neg_words] |= ((WORD_TYPE)1 << (elt & CACHE));
    }

    /*!
      Add element elt into the set or remove it if it is already contain [O(1)]
    */
    inline  void invert(const int elt)
    {
      int i = (elt >> EXP);
      if( (i >= neg_words) && 
	  (i <  pos_words) )
	table[i] ^= ((WORD_TYPE)1 << (elt & CACHE));
    }

    inline  void fast_invert(const int elt)
    {
      table[(elt >> EXP)] ^= ((WORD_TYPE)1 << (elt & CACHE));
    }

    inline  void word_invert(const int elt)
    {
      table[neg_words] ^= ((WORD_TYPE)1 << (elt & CACHE));
    }

    /*!
      Return true iff the set is empty [O(N/32)]
    */
    inline  bool empty()const
    { 
      int i = pos_words;
      while( i-- > neg_words ) 
	if(table[i]) return false;
      return true;  
    }

    /*!
      Return true iff the calling object intersect s [O(N/32)]
    */
    inline bool intersect(const Bitset<WORD_TYPE,FLOAT_TYPE>& s)const
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	if(table[i] & s.table[i]) return true;
      return false;
    }

    /*!
      Return true iff the calling object intersect s (s is assumed to be a bitset in {0,..,31}) [O(N/32)]
    */
    inline bool intersect(const int s) const
    {
      return(pos_words && neg_words<1 && (table[0]&s));
    }

    inline bool word_intersect(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
    {
      return ( table[neg_words] & s.table[neg_words] ) ;
    }

    inline  bool fast_intersect(const Bitset<WORD_TYPE,FLOAT_TYPE>& s, int& idx)const
    {
      if( table[idx] & s.table[idx] ) return true; 
      if( pos_words > neg_words ) {
	idx = pos_words;
	while( idx > neg_words ) {
	  --idx;
	  if(table[idx] & s.table[idx]) return true;
	}
      }
      return false;
    }

    /*!
      Return true iff the calling object intersect [lo..up] [O(N/32)]
    */
    inline   bool intersect(const int lb, const int ub)const
    {
      int i = (ub >> EXP);
      int j = (lb >> EXP);

      if( i < neg_words || j >= pos_words ) 
	return false;

      WORD_TYPE masked_lb = (full << (lb & CACHE));
      WORD_TYPE masked_ub = (full >> (CACHE - (ub & CACHE)));

      if( i == j ) {
	if( table[i] & (masked_lb & masked_ub) ) return true;
	else return false;
      }
      if( i >= pos_words )
	i = pos_words-1;
      else if( table[i--] & masked_ub ) return true;
  
      if( j < neg_words ) 
	j = neg_words;
      else if( table[j++] & masked_lb ) return true;

      while( i >= j )
	if(table[i--]) return true;
      return false;
    }

    /*!
      Increment by x all elements in the set.
      Any element greater or equal than the capacity 
      of the set minus x is removed [O(N/32)]
    */
    inline void increment(const int v)
    {
      int step = (v >> EXP); 
      int i = pos_words;
      int e = (v & CACHE);
      int f = size_word_bit-e;
      int j = neg_words+step;
      WORD_TYPE mask = ((WORD_TYPE)~0 << f);
      while( --i > j ) 
	table[i] = ((table[i-step] << e) | ((table[i-step-1] & mask) >> f));
      if( i >= neg_words+step ) table[i] = (table[i-step] << e);
      while( i > neg_words ) table[--i] = 0;
    }

    /*!
      Decrement by x all elements in the set.
      Any element lower than x is removed [O(N/32)]
    */
    inline  void decrement(const int v)
    {
      int step = (v >> EXP); 
      int i = neg_words-1;
      int e = (v & CACHE);
      int f = size_word_bit-e;
      int j = pos_words-step-1;
      WORD_TYPE mask = ((WORD_TYPE)~0 >> e);
      while( ++i < j )
	table[i] = ((table[i+step] >> e) | ((table[i+step+1] & mask) << f));
      if( i < pos_words-step ) table[i] = (table[i+step] >> e);
      while( ++i < pos_words ) table[i] = 0;
    }

    /*!
      Changes every value to its arythmetic negation
    */
    inline void negate( Bitset<WORD_TYPE,FLOAT_TYPE>& s ) const
    {
      int i = (pos_words > -s.neg_words ? -s.neg_words : pos_words);
      int j = (neg_words < -s.pos_words ? -s.pos_words : neg_words);

      unsigned int a; 
      WORD_TYPE mask, v, aux, rest = ( i < pos_words && (table[i] & 1) );

      while( i-- > j ) {
	aux = (table[i] & 1);
	v = (table[i] >> 1);
	mask = ~0;         
	a = sizeof(v) * CHAR_BIT; // bit size; must be power of 2 
	while ((a >>= 1) > 0) 
	  {
	    mask ^= (mask << a);
	    v = ((v >> a) & mask) | ((v << a) & ~mask);
	  }	
	s.table[-i-1] = (v | rest);
	rest = aux;
      }
      if(rest)
	s.table[i+1] |= rest;
    }

    /*!
      Add all elements between 0 to capacity [O(N/32)]
    */
    inline  void fill()
    {
      int i = pos_words;
      while( i > neg_words )
	table[--i] = full;
    }

    inline  void fill(const int lb, const int ub)
    {
      int i = (ub >> EXP);
      int j = (lb >> EXP);
      
      if( i >= neg_words || j < pos_words ) {
		
	WORD_TYPE masked_lb = (full << (lb & CACHE));
	WORD_TYPE masked_ub = (full >> (CACHE - (ub & CACHE)));
	
	if( i == j ) {

	  table[i] |= (masked_lb & masked_ub);
	  
	} else {
	  
	  if( i >= pos_words ) {
	    i = pos_words-1;
	  } else {
	    table[i--] |= masked_ub;
	  }
	  
	  if( j < neg_words ) {
	    j = neg_words;
	  } else {
	    table[j++] |= masked_lb;
	  }
	  
	  while( i >= j )
	    table[i--] |= full;
	}
      }
    }

    /*!
      Remove all elements [O(N/32)]
    */
    inline  void clear()
    {
      int i = pos_words;
      while( i > neg_words )
	table[--i] = empt;
    }

    /*!
      Remove all elements but v [O(N/32)]
    */
    inline  void set_to( const int v )
    {
      int i, j = (v >> EXP);
      for(i=neg_words; i<j; ++i)
	table[i] = empt;
      table[j] = ((WORD_TYPE)1 << v);
      for(i=j+1; i<pos_words; ++i) 
	table[i] = empt;
    }

    /*!
      flip all elements [O(N/32)]
    */
    inline  void flip()
    {
      int i = pos_words;
      while( i > neg_words )
	table[--i] ^= full;
    }

    /*!
      Remove all elements strictly lower than l [O(N/32)]
    */
    inline void set_min(const int bound)
    {
      int ith_word=(bound >> EXP);
      if( ith_word >= neg_words ) {
	if( ith_word <  pos_words ) {      
	  int i=ith_word;
	  while( i-- > neg_words ) table[i]=0;
	  table[ith_word] &= (full << (bound & CACHE));
	} else clear();
      }
    }

    /*!
      Remove all elements strictly greater than u [O(N/32)]
    */
    inline void set_max(const int bound)
    {
      int ith_word=(bound >> EXP);
      if( ith_word <  pos_words ) {
	if( ith_word >= neg_words ) {
	  int i=pos_words;
	  while( --i > ith_word ) table[i]=0;
	  table[ith_word] &= (full >> (CACHE - (bound & CACHE)));
	} else clear();
      }
    }

    /*!
      Remove all elements in the interval [l..u] [O(N/32)]
    */
    inline  void remove_interval(const int lb, const int ub)
    {
      if( lb <= ub ) {
	int lb_word = lb >> EXP;
	int ub_word = ub >> EXP;

	lb_word = ( lb_word < neg_words ? neg_words : lb_word );
	ub_word = ( ub_word >= pos_words ? pos_words-1 : ub_word );

	WORD_TYPE masked_lb = 0;
	WORD_TYPE masked_ub = 0;

	if( lb_word >= neg_words ) 
	  // add a '0' on the 32nd bit, because >> 32 does nothing
	  masked_lb = ((full/2) >> (CACHE - (lb & CACHE)));
	if( ub_word < pos_words ) 
	  masked_ub = ((full-1) << (ub & CACHE));

	if( lb_word == ub_word ) {
	  table[lb_word] &= (masked_lb | masked_ub);
	} else {
	  table[lb_word] &= masked_lb;
	  table[ub_word] &= masked_ub;
	  while( --ub_word > lb_word )
	    table[ub_word] = 0;
	}
      }
    }

    /*!
      Add all elements in the interval [l..u] [O(N/32)]
    */
    inline void add_interval(int lb, int ub)
    {

      if( lb <= ub ) {
	int lb_word = lb >> EXP;
	int ub_word = ub >> EXP;
    
	lb_word = ( lb_word < neg_words ? neg_words : lb_word );
	ub_word = ( ub_word >= pos_words ? pos_words-1 : ub_word );

	WORD_TYPE masked_lb = full;
	WORD_TYPE masked_ub = full;
	if( lb_word >= neg_words ) 
	  //masked_lb ^= (full >> (CACHE - (lb & CACHE) + 1));
	  masked_lb ^= ((full/2) >> (CACHE - (lb & CACHE)));
	if( ub_word < pos_words ) 
	  //masked_ub ^= ((full-1) << (ub & CACHE) );
	  masked_ub ^= ((full-1) << (ub & CACHE));

	if( lb_word == ub_word ) 
	  table[lb_word] |= (masked_lb & masked_ub);
	else {
	  table[lb_word] |= masked_lb;
	  table[ub_word] |= masked_ub;
	  while( --ub_word > lb_word )
	    table[ub_word] = full;
	}
      }
    }

    inline bool operator[](const int i)
    {
      return fast_contain(i);
    }
		
    std::string to_str() const {
			    std::ostringstream oss;
				oss << "{";
				//std::string rstr = std::string("{");
      if( !empty() ) {
	int last = NOVAL, cur=min(), aft;

	bool flag=false;
	do{
	  aft = next(cur);

	  if(aft != cur+1 || cur != last+1) {
	    if( flag )
	      //rstr += std::string(",");
			oss << ",";

  	    oss << (int)cur;

	    flag = true;
	  } else if(flag) {
	    //rstr += std::string("..");
		  oss << "..";
	    flag = false;
	  }
	  last = cur;
	  cur = aft;
	} while( cur != NOVAL && cur != last );
      }
      //rstr += std::string("}");
	  oss << "}";
      return oss.str();
    }

    std::ostream& display(std::ostream& os) const {
      os << "{";
      if( !empty() ) {
	int last = NOVAL, cur=min(), aft;

	bool flag=false;
	do{
	  aft = next(cur);

	  if(aft != cur+1 || cur != last+1) {
	    if( flag )
	      os << ",";
	    os << (int)cur;
	    flag = true;
	  } else if(flag) {
	    os << "..";
	    flag = false;
	  }
	  last = cur;
	  cur = aft;
	} while( cur != NOVAL && cur != last );
      }
      os << "}";
      return os;
    }

    void  print_bits(std::ostream & os) const 
    {
      os << "[";
      for(int i=neg_words; i<pos_words; ++i) {
				if(i) os << " ";
				showUint( table[i], os );
			}
      os << "]";
    }

  };


  /**********************************************
   * InlinedBitSet
   **********************************************/
  /// The data array is stored directly in the object

  template< class WORD_TYPE, class FLOAT_TYPE >
  class InlinedBitset : public Bitset<WORD_TYPE,FLOAT_TYPE> 
  {

  public:

    /**@name Parameters*/
    //@{    
    WORD_TYPE data[0];
    //@}
    
    static InlinedBitset<WORD_TYPE,FLOAT_TYPE>* make_new(const int lb, const int ub, const WORD_TYPE p)
    {
      int n_w = (lb >> Bitset<WORD_TYPE,FLOAT_TYPE>::EXP);
      int p_w = (ub >> Bitset<WORD_TYPE,FLOAT_TYPE>::EXP)+1;
      int size = p_w - n_w;
      
      void* mem = malloc(sizeof(InlinedBitset<WORD_TYPE,FLOAT_TYPE>) + sizeof(WORD_TYPE)*(size));
      return new (mem) InlinedBitset<WORD_TYPE,FLOAT_TYPE>(lb, ub, p);
    }

    InlinedBitset(const int lb, const int ub, const WORD_TYPE p)
    {
      initialise(lb,ub,p,&(data[0]));
    }
    
  };

  typedef Bitset< unsigned long long int, double > Bitset64;
  typedef Bitset< unsigned int, float > Bitset32;
  typedef InlinedBitset< unsigned long long int, double > iBitset64;
  typedef InlinedBitset< unsigned int, float > iBitset32;

#ifdef _BIT64
  
  typedef Bitset64 BitSet;

#else

  typedef Bitset32 BitSet;

#endif


  /**********************************************
   * Queue
   **********************************************/
  /// List of constraints for computing the GAC closure (fifo).
  
  class Queue {
  public:
    
    /**@name Parameters*/
    //@{
    int *next; // next element in the list (one for each elt + one for the head)
    int _head; // last_index+1
    int _tail; // last element
    int offset; // first index
    //@}

    /**@name Constructors*/
    //@{
    Queue() {
      next = NULL;
      _head = NOVAL;
      _tail = NOVAL;
      offset = 0;
    }
    ~Queue()
    {
      next += offset;
      delete [] next;
    }
    void cancel()
    {
      next = NULL;
      offset = 0;
    }
    bool is_initialised() 
    {
      return (next != NULL);
    }
    void initialise(const int n)
    {
      next = new int[n+1];
      std::fill( next, next+n, NOVAL );
      _head = _tail = n;
      next[n] = n;
      offset = 0;
    }
    void initialise(const int lb, const int ub)
    {
      next = new int[ub-lb+2];
      next-=lb;
      std::fill( next+lb, next+ub+1, NOVAL );
      _head = _tail = ub+1;
      next[ub+1] = ub+1;
      offset = lb;
    }
    void extend(const int elt)
    {
      int new_lb = (elt < offset ? elt : offset);
      int new_ub = (elt >= _head ? elt : _head-1);

      if(new_lb < offset || new_ub >= _head) {

	// 	std::cout << "EXTEND QUEUE" << std::endl;
	
	int *old_next = next;
	int old_head = _head;
	int old_tail = _tail;
	int old_offset = offset;

	// 	std::cout << "before: "
	// 		  << _head << " " 
	// 		  << _tail << " "
	// 		  << next[_head] 
	// 		  << std::endl;

	initialise(new_lb, new_ub);

	// 	std::cout << "afteri: "
	// 		  << _head << " " 
	// 		  << _tail << " "
	// 		  << next[_head] 
	// 		  << std::endl;


	// 	std::cout << " saved: "
	// 		  << old_head << " " 
	// 		  << old_tail << " "
	// 		  << old_next[old_head] 
	// 		  << std::endl;
	
	if(old_tail != old_head) {
	  int elt = old_next[old_head];
	  while(elt != old_head) {
	    //std::cout << "ADD " << elt << std::endl;
	    add(elt);
	    elt = old_next[elt];
	  }
	}

	// 	std::cout << " after: "
	// 		  << _head << " " 
	// 		  << _tail << " "
	// 		  << next[_head] 
	// 		  << std::endl;


	// 	int *aux = next;
	// 	next = new int[new_ub-new_lb+2];
	// 	std::fill(next, next+new_ub-new_lb+2, NOVAL);
	// 	next-=new_lb;
	// 	//memcpy(next+offset, aux+offset, (_head-offset+1)*sizeof(int));

	
	

	// // 	for(int i=new_lb; i<=new_ub+1; ++i)
	// // 	  if(next[i] != NOVAL)
	// // 	    std::cout << std::setw(3) << next[i] ;
	// // 	  else 
	// // 	    std::cout << " . ";
	// // 	std::cout << std::endl;
	// // 	for(int i=new_lb; i<=new_ub+1; ++i)	  
	// // 	  if(i>=offset && i<=_head && aux[i] != NOVAL)
	// // 	    std::cout << std::setw(3) << next[i] ;
	// // 	  else 
	// // 	    std::cout << " . ";
	// // 	std::cout << std::endl;


	// 	for(int i=new_lb; i<offset; ++i)
	// 	  next[i] = NOVAL;

	// 	for(int i=_head+1; i<new_ub+1; ++i)
	// 	  next[i] = NOVAL;

	// 	if(_tail == _head) {
	// 	  _tail = new_ub+1;
	// 	}
	
	// 	next[_tail] = new_ub+1;
	// 	next[new_ub+1] = next[_head];
	

	// 	for(int i=new_lb; i<=new_ub+1; ++i)
	// 	  if(next[i] != NOVAL)
	// 	    std::cout << std::setw(3) << next[i] ;
	// 	  else 
	// 	    std::cout << " . ";
	// 	std::cout << std::endl;
	// 	for(int i=new_lb; i<=new_ub+1; ++i)	  
	// 	  if(i>=offset && i<=_head && aux[i] != NOVAL)
	// 	    std::cout << std::setw(3) << next[i] ;
	// 	  else 
	// 	    std::cout << " . ";
	// 	std::cout << std::endl;

	// 	//exit(1);

	
	// // 	if(_head <= new_ub) {
	// // 	  next[new_ub+1] = aux[_head];
	// // 	  next[_head] = NOVAL;
	// // 	}

	// 	for(int i=new_lb; i<=new_ub+1; ++i)
	// 	  if(next[i] != NOVAL)
	// 	    std::cout << std::setw(3) << next[i] ;
	// 	  else 
	// 	    std::cout << " . ";
	// 	std::cout << std::endl;
	// 	for(int i=new_lb; i<=new_ub+1; ++i)	  
	// 	  if(i>=offset && i<=_head && aux[i] != NOVAL)
	// 	    std::cout << std::setw(3) << aux[i] ;
	// 	  else 
	// 	    std::cout << " . ";
	// 	std::cout << std::endl;

	// // 	//exit(1);

	// 	aux += offset;
	// 	delete [] aux;

	// 	offset = new_lb;
	// 	_head = new_ub+1;

	old_next += old_offset;
	delete [] old_next;

	//	display(std::cout);
      }
    }
    void declare(const int elt) {
      extend(elt);
      add(elt);
    }
    //@}
  
    /*!@name Accessors*/
    //@{
    /// Add an element to the queue
    inline int first() {
      return next[_head];
    }

    inline void add(const int elt)
    {
      next[_tail] = elt;
      _tail = elt;
      next[elt] = _head;
    }

    /// Pop the first event of the queue
    inline int pop()
    {
      int elt=next[_head];
      next[_head] = next[elt];
      if( next[_head] == _head ) _tail = _head;

      return elt;
    }

    /// Is the queue empty?
    inline bool empty() const
    {
      return _head == _tail;
    }

    /// Clear all elements of the queue
    inline void clear()
    {
      _tail = _head;
    }
    //@}
    
    Queue& operator=(const Queue& q) {
      _head = q._head;
      _tail = q._tail;
      next = q.next;
      offset = q.offset;
      return *this;
    }

    /*!@name Miscellanous*/
    //@{ 
    std::ostream& display(std::ostream& os) const {

      os << _head << " " << _tail << " " << offset << ": ";
      if(_head != NOVAL){
	for(int i=offset; i<=_head; ++i) {
	  if(next[i] != NOVAL)
	    os << next[i] << " ";
	  else 
	    os << ". ";
	}
      }
      os << "\n" ;

      if(_head != NOVAL){
	int i;
	os << _head << "(";
	os.flush();
	i = next[_head];
	if(i != _head) {
	  os << i;
	  i = next[i];
	  while( i != _head ) {
	    os << " " << i;
	    
	    if(i == NOVAL) break;
	    
	    i = next[i];
	  }
	}
	os << ")";
      }


      return os;
    }
    //@}

  };



	template <class T>
	class BinaryMinHeap {

	private:
    
		inline void sift_up(unsigned int index, const T x) {
			++num_operations;
      
			unsigned int ascendant = (index-1)/2;    
			while(index && x < data[ascendant]) {
				++num_operations;
	
				data[index] = data[ascendant];
				data[ascendant] = x;
				index = ascendant;
				ascendant /= 2;
			}
		}
    
		inline void sift_down(unsigned int index, const T x) {
			unsigned int descendant = index*2+1;
			unsigned int smallest = index;
			++num_operations;
      
			while(descendant < data.size) {
				++num_operations;
	
				if(data[descendant] < data[smallest]) {
					smallest = descendant;
				}
				++descendant;
				if(descendant < data.size && data[descendant] < data[smallest]) {
					smallest = descendant;
				}
				if(smallest == index) break;
	
				data[index] = data[smallest];
				data[smallest] = x;
	
				index = smallest;
				descendant = smallest*2+1;
			}
		}
    
	public: 
    
		int num_operations;
		Vector< T > data;
    
		BinaryMinHeap() { num_operations=0; }
		virtual ~BinaryMinHeap() { }
    
		void add(const T x) {
			data.add(x);
			sift_up(data.size-1, x);
		}
    
		T pop_min() {
			T the_min = data[0];
			data[0] = data[--data.size];
			data[data.size] = the_min;
			sift_down(0, data[0]);
			return the_min;
		}
    
		/*!@name Printing*/
		//@{
		std::ostream& display(std::ostream& os) const {
      

			unsigned int n_level = 0;
			while((unsigned int)(1<<n_level)<=data.size) ++n_level;
      
			unsigned int offset = 0;
			unsigned int step = 1;

			unsigned int i, j, l, first, last;

			for(i=0; i<30; ++i) {
				os << " " << data[i];
			}
			os << std::endl;
	
			for(l=0; l<n_level; ++l) {
				first = (1<<(n_level-l-1))-1;
				last = (1<<(n_level-l))-1;
				if(last>data.size) last = data.size;
				for(j=0; j<offset; ++j) os << "   ";
				for(i=first; i<last; ++i) {
					os << std::setw(3) << data[i];
					if(i<last-1) for(j=0; j<step; ++j) os << "   ";
				}
				os << std::endl;
				offset = offset+(1+step)/2;
				step = 2*step+1;
			}
			return os;
		}
		//@}
    
	};





	template <class T>
	class IndexedBinaryMaxHeap {

	private:
    
		// starting at position 'rank' in the tree, move up to the root swapping elements until this branch is sorted
		inline void sift_up(unsigned int rank) {
			unsigned int ascendant = (rank-1)/2;
			int x = heap[rank];
			while(rank && value[x] > value[heap[ascendant]]) {
				index[x] = ascendant;
				index[heap[ascendant]] = rank;
	
				heap[rank] = heap[ascendant];
				heap[ascendant] = x;

				rank = ascendant--;
				ascendant /= 2;
			}
		}
    
		// starting at position 'rank' in the tree, move down to a leaf swapping elements until this branch is sorted
		inline void sift_down(unsigned int rank) {
			unsigned int descendant = rank*2+1;
			unsigned int largest = rank;
			int x = heap[rank];
			while(descendant < heap.size) {
				if(value[heap[descendant]] > value[heap[largest]]) {
					largest = descendant;
				}
				++descendant;
				if(descendant < heap.size && value[heap[descendant]] > value[heap[largest]]) {
					largest = descendant;
				}
				if(largest == rank) break;
	
				index[heap[largest]] = rank;
				index[x] = largest;

				heap[rank] = heap[largest];
				heap[largest] = x;
	
				rank = largest;
				descendant = largest*2+1;
			}
		}
    
	public: 

		// the value of each element (indexed by element names)
		Vector< T >   value;

		// the data structure of the heap (populated with element names)
		Vector< int >  heap;

		// the current position of each element in the heap (indexed by element names)
		Vector< int > index;

    
		IndexedBinaryMaxHeap() {}
		virtual ~IndexedBinaryMaxHeap() { }
    
		int add(const T x) {
			int i = value.size;
			value.add(x);
			index.add(heap.size);
			heap.add(i);
			sift_up(heap.size-1);
      
			return i;
		}

		void change(const int i, const T new_value) {
			if(new_value > value[i]) {
				value[i] = new_value;
				sift_up(index[i]);
			} else {
				value[i] = new_value;
				sift_down(index[i]);
			}
		}

		T pop_max() {
			int the_root = heap[0];

			heap[0] = heap[--heap.size];
			heap[heap.size] = the_root;

			index[the_root] = heap.size;
			index[heap[0]] = 0;

			sift_down(0);
			return value[the_root];
		}



    
		/*!@name Printing*/
		//@{
		std::ostream& display(std::ostream& os) const {
      

			unsigned int n_level = 0;
			while((unsigned int)(1<<n_level)<=heap.size) ++n_level;
      
			unsigned int offset = 0;
			unsigned int step = 1;

			unsigned int i, j, l, first, last;

   
	
			for(l=0; l<n_level; ++l) {
				first = (1<<(n_level-l-1))-1;
				last = (1<<(n_level-l))-1;
				if(last>heap.size) last = heap.size;
				for(j=0; j<offset; ++j) os << "   ";
				for(i=first; i<last; ++i) {
					if(index[heap[i]] != i) {
						os << "!!!!" << std::endl;
						exit(1);
					}
					os << std::setw(3) << value[heap[i]] ;
					if(i<last-1) for(j=0; j<step; ++j) os << "   ";
				}
				os << std::endl;
				offset = offset+(1+step)/2;
				step = 2*step+1;
			}
			return os;
		}
		//@}
    
	};


  int __modulo_fct__(const int x, const int m) ;

  //class Variable;
  class BiInterval;
  class Interval {

  public:
    
    int min;
    int max;
    
    Interval();
    Interval(const BiInterval b);
    Interval(const int _min, const int _max);
    virtual ~Interval();
    
    Interval get_union(Interval arg);
    
    bool contain(const int x) const;
    bool empty() const;
    
    Interval operator*(const Interval);
    Interval operator/(const Interval);
    Interval anti_mul(const Interval);
    Interval operator-() const;
    Interval operator%(const int);
    Interval positive_modulo(const int);
    Interval operator_modulo(const int);

    Interval target_positive_modulo(const int, const Interval);
    Interval target_modulo(const int, const Interval);
    Interval target_c_modulo(const int, const Interval);

    void operator+=(const int x);
    void operator-=(const int x);

    std::ostream& display(std::ostream& os) const;
    
  };




  class Variable;
  class PositiveHalfDomain;
  class NegativeHalfDomain : public Interval {

public:
  
  NegativeHalfDomain();
  NegativeHalfDomain(const int _min, const int _max);
  virtual ~NegativeHalfDomain();


  Interval operator*(const PositiveHalfDomain arg);
  Interval operator*(const NegativeHalfDomain arg);
  Interval operator*(const int arg);

  Interval anti_mul(const PositiveHalfDomain arg);
  Interval anti_mul(const NegativeHalfDomain arg);
  Interval anti_mul(const int arg);

  Interval operator/(const PositiveHalfDomain arg);
  Interval operator/(const NegativeHalfDomain arg);
  Interval operator/(const int arg);
  Interval divided_by(const PositiveHalfDomain arg, const Variable target);
  Interval divided_by(const NegativeHalfDomain arg, const Variable target);

  Interval anti_div_X(const PositiveHalfDomain arg);
  Interval anti_div_X(const NegativeHalfDomain arg);
  Interval anti_div_X_pos(const int arg);
  Interval anti_div_X_neg(const int arg);

  Interval anti_div_Y(const PositiveHalfDomain arg);
  Interval anti_div_Y(const NegativeHalfDomain arg);
  Interval anti_div_Y_pos(const int arg);
  Interval anti_div_Y_neg(const int arg);
	

  // Interval operator%(const int mod);
  // Interval operator%(const int mod);

  // // the return value is the interval I such that I%mod = this
  // Interval anti_modulo(const int mod);
  // Interval anti_modulo(const int mod);

  //void operator=(const Interval I) {min = I.min; max = I.max;}
    //Interval operator-() { return Interval(-max, -min); }
  
};


class PositiveHalfDomain : public Interval {

public:
  
  PositiveHalfDomain();
  PositiveHalfDomain(const int _min, const int _max);
  virtual ~PositiveHalfDomain();

  Interval operator*(const PositiveHalfDomain arg);
  Interval operator*(const NegativeHalfDomain arg);
  Interval operator*(const int arg);

  Interval anti_mul(const PositiveHalfDomain arg);
  Interval anti_mul(const NegativeHalfDomain arg);
  Interval anti_mul(const int arg);

  Interval operator/(const PositiveHalfDomain arg);
  Interval operator/(const NegativeHalfDomain arg);
  Interval divided_by(const PositiveHalfDomain arg, const Variable target, const bool pos=true);
  Interval divided_by(const NegativeHalfDomain arg, const Variable target);
  Interval operator/(const int arg);

  Interval anti_div_X(const PositiveHalfDomain arg);
  Interval anti_div_X(const NegativeHalfDomain arg);
  Interval anti_div_X_pos(const int arg);
  Interval anti_div_X_neg(const int arg);

  Interval anti_div_Y(const PositiveHalfDomain arg);
  Interval anti_div_Y(const NegativeHalfDomain arg);
  Interval anti_div_Y_pos(const int arg);
  Interval anti_div_Y_neg(const int arg);
  // Interval operator%(const int mod);
  // Interval operator%(const int mod);

  // // the return value is the interval I such that I%mod = this
  // Interval anti_modulo(const int mod);
  // Interval anti_modulo(const int mod);


  //void operator=(const Interval I) {min = I.min; max = I.max;}
  //Interval operator-() { return Interval(-max, -min); }
};


  class Variable;
  class BiInterval {

  public:
    PositiveHalfDomain positive;
    NegativeHalfDomain negative;
    bool zero;

    BiInterval();
    BiInterval(const Variable x);
    BiInterval(const int _min, const int _max);
    BiInterval(const Interval I);
    BiInterval(const Interval neg, const Interval pos, const bool z);
    BiInterval(const Interval neg, const Interval pos, const Interval z);
    BiInterval(const int n_min, const int n_max, const int p_min, const int p_max, const bool z);
    void initialise(const int _min, const int _max);
    virtual ~BiInterval();
    

    int get_min() const;
    int get_max() const;

    int get_min_abs() const;
    int get_max_abs() const;

    int is_hollow() const;
    //Interval get_hole() const;

    BiInterval operator*(const BiInterval arg);
    BiInterval anti_mul(const BiInterval arg);
    BiInterval operator/(const BiInterval arg);
    BiInterval divided_by(const BiInterval arg, const Variable target);
    BiInterval anti_div_X(const BiInterval arg);
    BiInterval anti_div_Y(const BiInterval arg);

    BiInterval operator*(const int arg);
    BiInterval operator/(const int arg);
    BiInterval anti_div_X(const int arg);
    BiInterval anti_div_Y(const int arg);
    // BiInterval operator%(const int mod);
    // BiInterval anti_modulo(const int mod);

    bool operator==(const int x) const;
    void operator=(const int x);

    std::ostream& display(std::ostream& os) const;
    
  };


  class IntervalList : public Vector<Interval> {
    
  public:
    IntervalList();
    virtual ~IntervalList();

    void union_with(const IntervalList& with, IntervalList& into) const;
    void intersect_with(const IntervalList& with, IntervalList& into) const;

    void operator=(const IntervalList& l);

    void push(const int lb, const int ub);
    void push(const Interval& I);

  };

 


  // /**********************************************
  //  * VariableQueue
  //  **********************************************/
  // /// List of integer, used because we can do fifo push/pop
  // class VariableQueue {
  // public:

  //   /**@name Constructors*/
  //   //@{
  //   // index of the following element 
  //   int *next;
  //   // type of event
  //   int *trigger;
  //   // index of the first element
  //   int head;
  //   // last element
  //   int tail;
  //   // size of the structure
  //   int capacity;
  //   //@}

  //   /**@name Constructors*/
  //   //@{
  //   VariableQueue() {
  //     trigger = NULL;
  //     next = NULL;
  //     capacity = 0;
  //   }
  //   ~VariableQueue()
  //   {
  //     delete [] trigger;
  //     delete [] next;
  //   }
  //   void initialise(int n)
  //   {
  //     //      X = x;
  //     trigger = new int[n];
  //     std::fill( trigger, trigger+n, 0 );
  //     next = new int[n+2];
  //     std::fill( next, next+n, NOVAL );
  //     head = tail = n;
  //     next[n] = n;
  //     capacity = n;
  //   }

  //   void declare(int x) {
  //     if(x>=capacity) {
  // 	int new_capacity = 2*capacity;

  // 	int *new_trigger = new int[new_capacity];
  // 	memcpy(new_trigger, trigger, capacity*sizeof(int));
  // 	std::fill( new_trigger+capacity, new_trigger+new_capacity, 0 );

  // 	int *new_next = new int[new_capacity+2];
  // 	memcpy(new_next, next, capacity*sizeof(int));
  // 	std::fill( new_next+capacity, new_next+new_capacity, 0 );
  // 	std::fill( next+capacity, next+new_capacity, NOVAL );

  // 	if(head == tail) {
  // 	  head = tail = new_capacity;
  // 	  next[new_capacity] = new_capacity;
  // 	}

  // 	delete [] trigger;
  // 	delete [] next;
	
  // 	trigger = new_trigger;
  // 	next = new_trigger;
  // 	capacity = new_capacity;
  //     }
  //   }
  //   //@}
  
  //   /*!@name Accessors*/
  //   //@{
  //   /// Add an element to the queue
  //   inline void add(const int i, const int event)
  //   {
  //     //int i=v->id;
  //     if( !trigger[i] ) {
  // 	next[tail] = i;
  // 	tail = i;
  // 	next[i] = head;
  //     }
  //     trigger[i] |= event;
  //   }

  //   /// Pop the first event of the queue
  //   inline int pop( int& event )
  //   {
  //     int i=next[head];
  //     next[head] = next[i];
  //     event = trigger[i];
  //     trigger[i] = 0;

  //     if( next[head] == head )
  // 	tail = head;
  //     return i;
  //   }

  //   /// Is the queue empty?
  //   inline bool empty() const
  //   {
  //     return (next[head] == head);
  //   }

  //   /// Clear all elements of the queue
  //   inline void clear()
  //   {
  //     while( next[head] != head ) {
  // 	trigger[next[head]] = 0;
  // 	next[head] = next[next[head]];
  //     }
  //     tail = head;
  //   }
  //   //@}

  //   /*!@name Miscellanous*/
  //   //@{ 
  //   std::ostream& display(std::ostream& os) const
  //   {
  //     int i;
  //     os << "{";
  //     i = next[head];
  //     while( i != head ) {
  // 	os << i << " ";
  // 	i = next[i];
  //     }
  //     os << "}";
  //     return os;
  //   }
  //   //@}

  // };



  //   template < class DATA_TYPE > 
  //   std::string toString(const Vector< DATA_TYPE >& x) {
  //     std::ostringstream os;
  //     os << x;
  //     return os.str();
  //   }
  
  //   template < class DATA_TYPE > 
  //   std::string toString(const Stack< DATA_TYPE >& x) {
  //     std::ostringstream os;
  //     os << x;
  //     return os.str();
  //   }
  
  //   template < class DATA_TYPE > 
  //   std::string toString(const Array< DATA_TYPE >& x) {
  //     return x.getString();
  //   }

  //   template < class DATA_TYPE > 
  //   std::string toString(const Node< DATA_TYPE >& x) {
  //     return x.getString();
  //   }

  //   template < class DATA_TYPE, int NUM_HEAD > 
  //   std::string toString(const MultiList< DATA_TYPE, NUM_HEAD >& x) {
  //     return x.getString();
  //   }

  //   template< class WORD_TYPE, class FLOAT_TYPE >
  //   std::string toString(const Bitset< WORD_TYPE, FLOAT_TYPE >& x) {
  //     return x.getString();
  //   }



  std::ostream& operator<< (std::ostream& os, const Explanation& x);  

  std::ostream& operator<< (std::ostream& os, const Interval& x);  

  std::ostream& operator<< (std::ostream& os, const BiInterval& x);

  std::ostream& operator<< (std::ostream& os, const MultiSet& x);

  std::ostream& operator<< (std::ostream& os, const Queue& x);

  template < int N, class T > 
  std::ostream& operator<< (std::ostream& os, const Tuple< N, T >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const TwoWayStack< DATA_TYPE >& x) {
    return x.display(os);
  }

  template <class T1, class T2, class T3>
  std::ostream& operator<< (std::ostream& os, const Triplet<T1, T2, T3>& x) {
    return x.display(os);
  }


  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const BinaryMinHeap< DATA_TYPE >& x) {
    return  x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const IndexedBinaryMaxHeap< DATA_TYPE >& x) {
    return  x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Vector< DATA_TYPE >& x) {
    return x.display(os);
  }

  // template < class DATA_TYPE > 
  // std::ostream& operator<< (std::ostream& os, const Vector< DATA_TYPE >& x) {
  //   return x.display(os);
  // }

  template < class MAIN_TYPE, class AUX_TYPE > 
  std::ostream& operator<< (std::ostream& os, const BiStack< MAIN_TYPE, AUX_TYPE >& x) {
    return x.display(os);
  }

  // template < class DATA_TYPE > 
  // std::ostream& operator<< (std::ostream& os, const Stack< DATA_TYPE >& x) {
  //   return x.display(os);
  // }

  template < class DATA_TYPE, class SIZE_TYPE > 
  std::ostream& operator<< (std::ostream& os, const VarStack< DATA_TYPE, SIZE_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const ConStack< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Stack< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Array< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Node< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE, int NUM_HEAD > 
  std::ostream& operator<< (std::ostream& os, const MultiList< DATA_TYPE, NUM_HEAD >& x) {
    return x.display(os);
  }

  template< class WORD_TYPE, class FLOAT_TYPE >
  std::ostream& operator<< (std::ostream& os, const Bitset< WORD_TYPE, FLOAT_TYPE >& x) {
    return x.display(os);
  }


  std::ostream& operator<< (std::ostream& os, const Explanation* x);  

  std::ostream& operator<< (std::ostream& os, const IntStack* x);

  std::ostream& operator<< (std::ostream& os, const Queue* x);

  std::ostream& operator<< (std::ostream& os, const MultiSet* x);

  template < int N, class T > 
  std::ostream& operator<< (std::ostream& os, const Tuple< N, T >* x) {
    return x->display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const TwoWayStack< DATA_TYPE >* x) {
    return x->display(os);
  }

  template <class T1, class T2, class T3>
  std::ostream& operator<< (std::ostream& os, const Triplet<T1, T2, T3>* x) {
    return x->display(os);
  }


  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const BinaryMinHeap< DATA_TYPE >* x) {
    return (x ? x->display(os) : os);
  }

 template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const IndexedBinaryMaxHeap< DATA_TYPE >* x) {
    return (x ? x->display(os) : os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Vector< DATA_TYPE >* x) {
    return (x ? x->display(os) : (os << "nill"));
  }

  // template < class DATA_TYPE > 
  // std::ostream& operator<< (std::ostream& os, const Vector< DATA_TYPE >* x) {
  //   return (x ? x->display(os) : (os << "nill"));
  // }

  template < class MAIN_TYPE, class AUX_TYPE > 
  std::ostream& operator<< (std::ostream& os, const BiStack< MAIN_TYPE, AUX_TYPE >* x) {
    return x->display(os);
  }

  // template < class DATA_TYPE > 
  // std::ostream& operator<< (std::ostream& os, const Stack< DATA_TYPE >* x) {
  //   return (x ? x->display(os) : os);
  // }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Array< DATA_TYPE >* x) {
    return (x ? x->display(os) : os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Node< DATA_TYPE >* x) {
    return (x ? x->display(os) : os);
  }

  //   template < class DATA_TYPE, int NUM_HEAD > 
  //   std::ostream& operator<< (std::ostream& os, const MultiList< DATA_TYPE, NUM_HEAD >* x) {
  //     return x->display(os);
  //   }

  template< class WORD_TYPE, class FLOAT_TYPE >
  std::ostream& operator<< (std::ostream& os, const Bitset< WORD_TYPE, FLOAT_TYPE >* x) {
    return (x ? x->display(os) : os);
  }

  //   template < class DATA_TYPE, int NUM_HEAD >
  //   std::ostream& MultiList< DATA_TYPE, NUM_HEAD >::display(std::ostream& os) const {
  //     //for(unsigned int i=0; i<data.size; ++i)
  //     //os << data[i];
    
  //     for(int k=0; k<NUM_HEAD; ++k) {
      
  //       int min_elt =  NOVAL;
  //       int max_elt = -NOVAL;
      
  //       Node< DATA_TYPE > nd = data.stack_[head[k]];
  //       while(next(nd)) {
  // 	if(min_elt > (int)nd) min_elt = (int)nd;
  // 	if(max_elt < (int)nd) max_elt = (int)nd;
  //       }
      
  //       if(min_elt == NOVAL) {
  // 	min_elt = 0;
  // 	max_elt = 1;
  //       }
      
  //       BitSet tmp(min_elt, max_elt, BitSet::empt);
      
  //       nd = data.stack_[head[k]];
  //       while(next(nd)) {
  // 	tmp.add((int)nd);
  //       }
      
  //       os << tmp;
  //     }
  //     return os;
  //   }



	class Graph {
	
	public:
		int       capacity;
		int      num_edges;
		IntStack      node;
		IntStack *neighbor;
		
		BitSet     *matrix;
		BitSet    node_set;
		
		std::vector<int> matching;
		int nmatch;
		bool matching_is_max;
		
		std::vector<int> indexlist;
	
		Graph() {}
	
		Graph(const int n, const int full_nodes=true, const int full_edges=false) {
			initialise(n, full_nodes, full_edges);
		}
	
		Graph(const Graph& g) {
				
			num_edges = g.num_edges;
			//initialise(g.capacity,g.size() == g.capacity);
			initialise(g.capacity, true, false);
			// if(g.size() != g.capacity)
			// 	for(int i=0; i<g.node.size; ++i) { node.add(g.node[i]); }
			//if(g.size() != g.capacity)
			for(int i=g.size(); i<g.capacity; ++i) { node.remove(g.node[i]); }
			for(int x=0; x<capacity; ++x) {
				matrix[x].copy(g.matrix[x]);
				for(unsigned int i=0; i<g.neighbor[x].size; ++i) {
					neighbor[x].add(g.neighbor[x][i]);
				}
			}
			node_set.copy(g.node_set);
			
			assert(num_edges = g.num_edges);
	
		}
	
		virtual void initialise(const int n, const bool full_nodes=true, const bool full_edges=false) {
			
			//std::cout << "initialise " << n << " " << full_nodes << " " << full_edges << std::endl;
			
			capacity = n;
			node.initialise(0,capacity-1,capacity,full_nodes);
			//node_set.initialise(0,capacity-1, BitSet::full);
			node_set.initialise(0,capacity-1,(full_nodes ? ~0 : 0));
				//(full ? BitSet::full : BitSet::empt));
			neighbor = new IntStack[capacity];
			matrix = new BitSet[capacity];
			for(int i=0; i<capacity; ++i) {
				//matrix[i].initialise(0, capacity-1, (_full ? BitSet::full : BitSet::empt));				
				matrix[i].initialise(0, capacity-1, (full_edges ? ~0 : 0));
				neighbor[i].initialise(0, capacity-1, capacity, true);
				neighbor[i].remove(i);
				if(!full_edges)
					neighbor[i].clear();
			}
			
			num_edges = (full_edges ? (n*(n-1)/2) : 0);
			
			indexlist.reserve(capacity);
			matching.assign(capacity,-1);
			nmatch=0;
			matching_is_max = false;
		}
	
		virtual ~Graph() {
			delete [] neighbor;
			delete [] matrix;
		}
	
		int size() const { return node.size; }
	
		void clear() { 
			while(!node.empty()) { 
				neighbor[node.pop()].clear(); 
			} 
		}
	
		void add_node(const int x) {
			node.add(x);
			node_set.add(x);
		}
	
		// void add_directed(const int x, const int y) {
		// 	neighbor[x].add(y);
		// 	matrix[x].fast_add(y);
		// }
	
		void add_undirected(const int x, const int y) {
			neighbor[x].add(y);
			neighbor[y].add(x);

			matrix[x].fast_add(y);
			matrix[y].fast_add(x);
			
			++num_edges;
		}
		
		void rem_undirected(const int x, const int y) {
			neighbor[x].remove(y);
			neighbor[y].remove(x);

			matrix[x].fast_remove(y);
			matrix[y].fast_remove(x);
			
			--num_edges;
		}
	
		bool exist_arc(const int x, const int y) const {
			return neighbor[x].contain(y);
		}
	
		int degree(const int x) const {
			return neighbor[x].size;
		}
	
		int get_neighbor(const int x, const int i) const {
			return neighbor[x][i];
		}
	
		int get_next_neighbor(const int x, const int y) const {
			return neighbor[x].next(y);
		}
	
		int get_prev_neighbor(const int x, const int y) const {
			return neighbor[x].prev(y);
		}
		
		void find_maximal_matching() {
		
			indexlist.clear();
			for (int i=0; i<size(); i++) indexlist.push_back(i);
			std::random_shuffle (indexlist.begin(),indexlist.end());
		
			int u,v;
			matching.assign(capacity,-1);
			nmatch=0;
			for (int i=0; i<size(); i++) {
				u = node[indexlist[i]];
				if (!is_matched(u)) {
					for (unsigned int j=0; j<neighbor[u].size; j++) {
						v = neighbor[u][j];
						if (!is_matched(v)) {
							matching[u] = v;
							matching[v] = u;
							nmatch++;
							break;
						}
					}
				}
			}
		}
		
		bool is_matched(const int u) {
			return (matching[u] != -1);
		}
	
		std::ostream& display(std::ostream& os) const {
			for(unsigned int i=0; i<node.size; ++i) {
				os << node[i] << ": ";
				neighbor[node[i]].display(os);
				os << std::endl;
			}
			return os;
		}
	
	};


	std::ostream& operator<< (std::ostream& os, const Graph& x);

	std::ostream& operator<< (std::ostream& os, const Graph* x);
	
	class CompactGraph {
	
	public:
		int          capacity;
		int         num_edges;
		IntStack         node;
		Vector<int> *neighbor;
		Vector<int> *nb_index;
		
		BitSet       node_set;
		
#ifdef _NOT_SO_COMPACT
		BitSet        *matrix;
#endif
				
		std::vector<int> matching;
		int nmatch;
		
		std::vector<int> indexlist;
	
		CompactGraph() {}
	
		CompactGraph(const int n) {
			initialise(n);
		}
	
		CompactGraph(const CompactGraph& g) {
			initialise(g.capacity);
			num_edges = g.num_edges;
			node_set.copy(g.node_set);
			for(int i=g.size(); i<g.capacity; ++i) { node.remove(g.node[i]); }
			for(int x=0; x<capacity; ++x) {
				for(unsigned int i=0; i<g.neighbor[x].size; ++i) {
					neighbor[x].add(g.neighbor[x][i]);
				}
				for(unsigned int i=0; i<g.nb_index[x].size; ++i) {
					nb_index[x].add(g.nb_index[x][i]);
				}
				
#ifdef _NOT_SO_COMPACT
				matrix[x].copy(g.matrix[x]);
#endif
				
			}
		}
	
		virtual void initialise(const int n) {
			capacity = n;
			
			node.initialise(0,capacity-1,capacity,true);
			
			node_set.initialise(0, capacity-1, BitSet::empt);
						
			neighbor = new Vector<int>[capacity];
			nb_index = new Vector<int>[capacity];
			
			num_edges = 0;
			
			indexlist.reserve(capacity);
			matching.assign(capacity,-1);
			nmatch=0;
			
#ifdef _NOT_SO_COMPACT
			matrix = new BitSet[capacity];
			for(int i=0; i<capacity; ++i) {
				matrix[i].initialise(0, capacity-1, BitSet::empt);
			}
#endif
			
			
		}
	
		virtual ~CompactGraph() {
			delete [] neighbor;
			delete [] nb_index;
			
#ifdef _NOT_SO_COMPACT
			delete [] matrix;
#endif
			
		}
	
		int size() const { return node.size; }
		
		
		void flip() {
			
			//std::cout << "flip graph" << std::endl << " 1/ clear nodes and edges" << std::endl;
			
			num_edges = 0;
			//node.clear();
			//node_set.clear();
			
			for(int x=0; x<capacity; ++x) {
				
				//std::cout << " " << (x+2) << "/ clear neighbors of " << x << std::endl;
				
				neighbor[x].clear();
				nb_index[x].clear(); 
			}
			for(int x=0; x<capacity; ++x) {
				
				//std::cout << "  flip matrix[" << x << "] = " << matrix[x] << std::endl; 
				
				matrix[x].flip();
				matrix[x].set_max(capacity-1);
				matrix[x].remove(x);
				
				//std::cout << "                => " << matrix[x] << std::endl; 
				
				int y=x, nxt = matrix[x].next(x);
				while(y!=nxt) {
					y = nxt;
					
					//std::cout << "   -> add " << x << "," << y << std::endl;
					
					//add_undirected(x,y);
					add_undirected_nm(x,y);
					
					nxt = matrix[x].next(y);
				}
			}
		}
		
	
		void clear() { 
			int x;
			while(!node.empty()) { 
				x = node.pop();
				neighbor[x].clear(); 
				nb_index[x].clear(); 
			} 
		}
	
		void add_node(const int x) {
			node.add(x);
		}
	
		void add_undirected(const int x, const int y) {
			nb_index[y].add(neighbor[x].size);
			nb_index[x].add(neighbor[y].size);
			neighbor[x].add(y);
			neighbor[y].add(x);
			
#ifdef _NOT_SO_COMPACT			
			matrix[x].add(y);
			matrix[y].add(x);
#endif
			
			++num_edges;
		}
		
		void add_undirected_nm(const int x, const int y) {
			nb_index[y].add(neighbor[x].size);
			nb_index[x].add(neighbor[y].size);
			neighbor[x].add(y);
			neighbor[y].add(x);
	
			++num_edges;
		}
	
		int degree(const int x) const {
			return neighbor[x].size;
		}
		
		void find_maximal_matching() {
		
			indexlist.clear();
			for (int i=0; i<size(); i++) indexlist.push_back(i);
			std::random_shuffle (indexlist.begin(),indexlist.end());
		
			int u,v;
			matching.assign(capacity,-1);
			nmatch=0;
			for (int i=0; i<size(); i++) {
				u = node[indexlist[i]];
				if (!is_matched(u)) {
					for (unsigned int j=0; j<neighbor[u].size; j++) {
						v = neighbor[u][j];
						if (!is_matched(v)) {
							matching[u] = v;
							matching[v] = u;
							nmatch++;
							break;
						}
					}
				}
			}
		}
		
		bool is_matched(const int u) {
			return (matching[u] != -1);
		}
	
		std::ostream& display(std::ostream& os) const {
			for(unsigned int i=0; i<node.size; ++i) {
				int x = node[i];
				
				
				os << x << ": ";
				neighbor[x].display(os);
				// os << " ";
				// matrix[x].display(os);
				
				
				// for(int j=0; j<nb_index[x].size; ++j) {
				// 	int y = neighbor[x][j];
				// 	os << " " << neighbor[y][nb_index[x][j]];
				// }
				// os << " " ;
				// nb_index[node[i]].display(os);
				
				os << std::endl;
			}
			return os;
		}
	
	};


	std::ostream& operator<< (std::ostream& os, const CompactGraph& x);

	std::ostream& operator<< (std::ostream& os, const CompactGraph* x);
	
	
	/* Temporary bipartite graph class. Supports maximum matching algorithms.
	   TempBipartiteGraph does not inherit from Graph since it uses very different data structures.
	 
	class TempBipartiteGraph {
	
	public:
	
		// col_ids contains all adjacency vectors. The adjacency vector of a left-side vertex v starts at index col_ptr[v] and ends at index col_ptr[v+1]-1.
		int *col_ids;
		int *col_ptrs;
		
		// If v is a left-side vertex, then match[v] is the right-side vertex matched with v, and -1 otherwise.
		int *match;
		// Same thing as match, but for right-side vertices.
		int *row_match;
		
		// (n,m) : number of (left-side,right-side) vertices
		int n;
		int m;
		
		// Structures for the matching algorithm
		int* visited;
		int* stack;
		int* colptrs;
		int* lookahead;
		int* unmatched;
			
		int leftsize;
		int rightsize;
		
		int nmatched;
		vector<int> matching;
		
		bool is_matched(int u) {
			return (matching[u] >= 0);
		}
		
		void greedy_matching()
		{
			int v;
  			for(size_t u = 0; u != leftsize; u++) {
  				for(int i = 0; i != neighbor[u].size; i++) {
				v = neighbor[u][i];
  				if (!is_matched(v)) {
  					matching[v] = u;
  					matching[u] = v;
  					++nmatched;
  					break;
  				}
  			}
		}
		
		TempBipartiteGraph(const int* col_ptrs_in, const int* col_ids_in, const int n_in, const int m_in) {
			col_ids = col_ids_in;
			col_ptrs_in = col_ptrs_in;
			n = n_in;
			m = m_in;
			match = (int*)malloc(sizeof(int)*n);
			row_match = (int*)malloc(sizeof(int)*m);
			visited = (int*)malloc(sizeof(int)*m);
			stack = (int*)malloc(sizeof(int)*n);
			colptrs = (int*)malloc(sizeof(int)*n);
			lookahead = (int*)malloc(sizeof(int)*n);
			unmatched = (int*)malloc(sizeof(int)*n);
			
			memset(match, -1, sizeof(int)*n);
		}
		
		void find_matching() {
			old_cheap();
			match_pf_fair();
		}
		
		void old_cheap() {
			int ptr;
			for(i=0; i < n; i++) {
				int s_ptr = col_ptrs[i];
				int e_ptr = col_ptrs[i + 1];
				for(ptr = s_ptr; ptr < e_ptr; ptr++) {
					int r_id = col_ids[ptr];
					if(row_match[r_id] == -1) {
						match[i] = r_id;
						row_match[r_id] = i;
						break;
					}
				}
			}
		}

		
		void match_pf_fair() {

			int	i, j, row, col, stack_col, temp, ptr, eptr, stack_last,
				stop = 0, pcount = 1, stack_end_ptr, nunmatched = 0, nextunmatched = 0,
				current_col, inc = 1;

			memset(visited, 0, sizeof(int) * m);
			memcpy(lookahead, col_ptrs, sizeof(int) * n);

			for(i = 0; i < n; i++) {
				if(match[i] == -1 && col_ptrs[i] != col_ptrs[i+1]) {
					unmatched[nunmatched++] = i;
				}
			}

			while(!stop) {
				stop = 1; stack_end_ptr = n;
				if(inc) {
					for(i = 0; i < nunmatched; i++) {
						current_col = unmatched[i];
						stack[0] = current_col; stack_last = 0; colptrs[current_col] = col_ptrs[current_col];

						while(stack_last > -1) {
							stack_col = stack[stack_last];

							eptr = col_ptrs[stack_col + 1];
							for(ptr = lookahead[stack_col]; ptr < eptr && row_match[col_ids[ptr]] != -1; ptr++){}
							lookahead[stack_col] = ptr + 1;

							if(ptr >= eptr) {
								for(ptr = colptrs[stack_col]; ptr < eptr; ptr++) {
									temp = visited[col_ids[ptr]];
									if(temp != pcount && temp != -1) {
										break;
									}
								}
								colptrs[stack_col] = ptr + 1;

								if(ptr == eptr) {
									if(stop) {stack[--stack_end_ptr] = stack_col;}
									--stack_last;
									continue;
								}

								row = col_ids[ptr]; visited[row] = pcount;
								col = row_match[row]; stack[++stack_last] = col; colptrs[col] = col_ptrs[col];
							
							} else {
							
								row = col_ids[ptr]; visited[row] = pcount;
								while(row != -1){
									col = stack[stack_last--];
									temp = match[col];
									match[col] = row; row_match[row] = col;
									row = temp;
								}
								stop = 0;
								break;
							}
						}

						if(match[current_col] == -1) {
							if(stop) {
								for(j = stack_end_ptr + 1; j < n; j++) {
									visited[match[stack[j]]] = -1;
								}
								stack_end_ptr = n;
							} else {
								unmatched[nextunmatched++] = current_col;
							}
						}
					}
				} else {
					for(i = 0; i < nunmatched; i++) {
						current_col = unmatched[i];
						stack[0] = current_col; stack_last = 0; colptrs[current_col] = col_ptrs[current_col + 1] - 1;

						while(stack_last > -1) {
							stack_col = stack[stack_last];

							eptr = col_ptrs[stack_col + 1];
							for(ptr = lookahead[stack_col]; ptr < eptr && row_match[col_ids[ptr]] != -1; ptr++){}
							lookahead[stack_col] = ptr + 1;

							if(ptr >= eptr) {
								eptr = col_ptrs[stack_col] - 1;
								for(ptr = colptrs[stack_col]; ptr > eptr; ptr--) {
									temp = visited[col_ids[ptr]];
									if(temp != pcount && temp != -1) {
										break;
									}
								}
								colptrs[stack_col] = ptr - 1;

								if(ptr == eptr) {
									if(stop) {stack[--stack_end_ptr] = stack_col;}
									--stack_last;
									continue;
								}

								row = col_ids[ptr]; visited[row] = pcount;
								col = row_match[row]; stack[++stack_last] = col;
								colptrs[col] = col_ptrs[col + 1] - 1;

							} else {
								row = col_ids[ptr]; visited[row] = pcount;
								while(row != -1){
									col = stack[stack_last--];
									temp = match[col];
									match[col] = row; row_match[row] = col;
									row = temp;
								}
								stop = 0;
								break;
							}
						}

						if(match[current_col] == -1) {
							if(stop) {
								for(j = stack_end_ptr + 1; j < n; j++) {
									visited[match[stack[j]]] = -1;
								}
								stack_end_ptr = n;
							} else {
								unmatched[nextunmatched++] = current_col;
							}
						}
					}
				}
				pcount++; nunmatched = nextunmatched; nextunmatched = 0; inc = !inc;
			}
		}
	
	}
	*/
	
	class BipartiteGraph : public Graph {
	
	public:
	
		IntStack left;
		IntStack right;
		
		IntStack vcover;
		
		// Left-side vertices that belong to ALL maximum matchings
		IntStack lpmatched;
		
		// structures for matching
		std::vector<int> leftfrontier;
		std::vector<int> rightfrontier;
		std::vector<int> rightvisited;
		std::vector<int> leftvisited;
		std::vector<int> rightvisited_toclear;
		std::vector<int> leftvisited_toclear;
		std::vector<int> leftbackp;
		std::vector<int> rightbackp;
		
		// structures for vertex cover
		std::vector<int> leftaccessed;
		std::vector<int> rightaccessed;
		
		Vector<int> trail;
		Vector<int> lefttrail;
		Vector<int> righttrail;
		
		BipartiteGraph() {}
		
		BipartiteGraph(const int n)
		{
			initialise(n);
		}
		
		void initialise(const int n, const bool full_nodes=true, const bool full_edges=false) 
		{
			capacity = n;
			node.initialise(0,capacity-1,capacity,full_nodes);
			node_set.initialise(0,capacity-1,(full_nodes ? ~0 : 0));
			neighbor = new IntStack[capacity];
			matrix = new BitSet[capacity];
			for(int i=0; i<capacity; ++i) {		
				matrix[i].initialise(0, capacity-1, (full_edges ? ~0 : 0));
				neighbor[i].initialise(0, capacity-1, capacity, true);
				neighbor[i].remove(i);
				if(!full_edges)
					neighbor[i].clear();
			}
			
			num_edges = (full_edges ? (n*(n-1)/2) : 0);
			
			indexlist.reserve(capacity);
			matching.assign(capacity,-1);
			nmatch=0;
			
			left.initialise(0,capacity-1,capacity,false);
			right.initialise(0,capacity-1,capacity,false);
			
			leftvisited.resize(capacity,false);
			rightvisited.resize(capacity,false);
			leftaccessed.resize(capacity,false);
			rightaccessed.resize(capacity,false);
			leftbackp.assign(capacity,-1);
			rightbackp.assign(capacity,-1);
			
			vcover.initialise(0,capacity-1,capacity,false);
			lpmatched.initialise(0,capacity-1,capacity,false);
		}
		
		void set_left(const int v) {
			left.add(v);
		}
		
		void set_right(const int v) {
			right.add(v);
		}
		
		void remove_and_save(const int v) {
		
			if (is_matched(v)) {
				matching[matching[v]] = -1;
				matching[v] = -1;
				nmatch--;
				matching_is_max = false;
				if (vcover.contain(v)) vcover.remove(v);
			}
			
			node.remove(v);
			//EH CHANGE?? for (int i=0; i<neighbor[v].size; i++) neighbor[v][i].remove(v);
			for (unsigned int i=0; i<neighbor[v].size; i++) neighbor[neighbor[v][i]].remove(v);
			
			trail.add(v);
			
			if (left.contain(v)) {
				left.remove(v);
				lefttrail.add(true);
			} else {
				right.remove(v);
				lefttrail.add(false);
			}
			
		}

		void remove_extended_neighbourhood_and_save(const int v) {
		
			while (neighbor[v].size > 0) {
				remove_and_save(neighbor[v][0]);
			}
			remove_and_save(v);
		
		}
		
		
		void restore(const int v, const bool isleft) {
			
			node.add(v);
			isleft ? left.add(v) : right.add(v);
			//EH CHANGE?? for (int i=0; i<neighbor[v].size; i++) neighbor[v][i].add(v);
			for (unsigned int i=0; i<neighbor[v].size; i++) neighbor[neighbor[v][i]].add(v);
			matching_is_max = false; // Maybe we should keep track of the old matching and restore it ?
			
		}
		
		void check() {
			int u,v;
			for (unsigned int i=0; i<left.size; i++) {
				u = left[i];
				for (unsigned int j=0; j<neighbor[u].size; j++) {
					v = neighbor[u][j];
					assert(right.contain(v));
				}
			}
			for (unsigned int i=0; i<right.size; i++) {
				u = right[i];
				for (unsigned int j=0; j<neighbor[u].size; j++) {
					v = neighbor[u][j];
					assert(left.contain(v));
				}
			}
			for (unsigned int i=0; i<node.size; i++) {
				assert(left.contain(node[i]) || right.contain(node[i]));
			}
		}
				
		
		void clear_visited()
		{
  			for(size_t i = 0; i != rightvisited_toclear.size(); ++i) rightvisited[rightvisited_toclear[i]] = false;
  			for(size_t i = 0; i != leftvisited_toclear.size(); ++i) leftvisited[leftvisited_toclear[i]] = false;
  			rightvisited_toclear.clear();
  			leftvisited_toclear.clear();
		}
		
		void maximum_matching()
		{
			greedy_matching();
			find_matching();
		}
		
		void persistently_matched_left_vertices()
		{
			maximum_matching();
			
			leftaccessed.clear();
			leftfrontier.clear();
    		rightfrontier.clear();
    		leftvisited.assign(capacity,false);
			rightvisited.assign(capacity,false);
			int u,v;
			
			// Identify all vertices reacheable from an unmatched vertex via an alternating path of even length (using breadth-first search)
    			for (unsigned int i=0; i<left.size; i++) {
    				if (is_matched(left[i])) continue;
    				leftfrontier.push_back(left[i]);
    				leftvisited[left[i]] = true;
    				leftaccessed[left[i]] = true;
    				do {
    					while (!leftfrontier.empty()) {
    						u = leftfrontier.back();
    						leftfrontier.pop_back();
    						for (unsigned int q=0; q<neighbor[u].size; q++) {
    							v = neighbor[u][q];
    							if (matching[u] == v) continue;
    							if (rightvisited[v]) continue;
    							assert(is_matched(v)); // otherwise we would have an augmenting path !
    							rightfrontier.push_back(v);
    							rightvisited[v] = true;
    						}
    					}
    					while (!rightfrontier.empty()) {
						v = rightfrontier.back();
						rightfrontier.pop_back();
						int next = matching[v];
						if ((next == -1) || (leftvisited[next])) continue;
						leftfrontier.push_back(next);
						leftvisited[next] = true;
						leftaccessed[next] = true;
					}
				} while (!leftfrontier.empty());
				clear_visited();
			}
			
			for (unsigned int i=0; i<left.size; i++) {
				if (!leftaccessed[left[i]]) lpmatched.add(left[i]);
			}
		}
		
		void vertex_cover()
		{
			greedy_matching();
			find_matching();
			matching_to_vertex_cover();
			assert((int)(vcover.size) == nmatch);
		}
		
		void check_vertex_cover()
		{
			int u;
			for (unsigned int i=0; i<left.size; i++) {
				u = left[i];
				if (vcover.contain(u)) continue;
				for (unsigned int j=0; j<neighbor[u].size;j++) {
					assert(vcover.contain(neighbor[u][j]));
				}
			}
		}
		
		void check_matching()
		{
			for (int i=0; i<(int)(matching.size()); i++) {
				if (matching[i] != -1) {
					assert(matching[matching[i]] == i);
				}
			}
		}
		
		void output_matching()
		{
			for (unsigned int i=0; i<matching.size(); i++) {
				std::cout << "[" << i << "," << matching[i] << "]" << std::endl;
			}
		}
		
		void greedy_matching()
		{
			matching.assign(capacity,-1);
			nmatch=0;
			int u,v;
  			for(unsigned int i = 0; i < left.size; i++) {
  				u = left[i];
    				for(unsigned int q = 0; q < neighbor[u].size; q++) {
    					v = neighbor[u][q];
      					if( !is_matched(v) ) {
        					matching[u] = v;
        					matching[v] = u;
        					nmatch++;
        					break;
      					}
    				}
  			}
		}
		
		void find_matching()
		{
			int u,v;
			for (unsigned int fleft = 0; (fleft < left.size) && (nmatch < (int)(left.size)) && (nmatch < (int)(right.size)); fleft++) {
				if (is_matched(left[fleft])) continue;
				// do a bfs for an augmenting path
				leftfrontier.push_back(left[fleft]);
				int pathend = -1;
				do {
					while (!leftfrontier.empty()) {
						u = leftfrontier.back();
						assert(u != -1);
						leftfrontier.pop_back();
						for (size_t q = 0; q < neighbor[u].size; q++) {
							v = neighbor[u][q];
							assert(v != -1);
							if (matching[u] == v) continue;
							if (rightvisited[v]) continue;
							rightbackp[v] = u;
							if (!is_matched(v)) {
								pathend = v;
								goto augment_path;
							}
							rightfrontier.push_back(v);
							rightvisited[v] = true;
							rightvisited_toclear.push_back(v);
						}
					}
					while (!rightfrontier.empty()) {
						v = rightfrontier.back();
						assert(v != -1);
						rightfrontier.pop_back();
						int next = matching[v];
						if ((next == -1) || (leftvisited[next])) continue;
						leftfrontier.push_back(next);
						leftbackp[next] = v;
						leftvisited[next] = true;
						leftvisited_toclear.push_back(next);
					}
				} while (!leftfrontier.empty());
				clear_visited();
				continue;
				augment_path:
					int rn = pathend;
					assert(rn != -1);
					int ln = rightbackp[rn];
					assert(ln != -1);
					matching[ln] = rn;
					matching[rn] = ln;
					while (ln != left[fleft]) {
						rn = leftbackp[ln];
      						ln = rightbackp[rn];
      						matching[ln] = rn;
      						matching[rn] = ln;
      					}
      					nmatch++;
      					leftfrontier.clear();
    					rightfrontier.clear();
    					clear_visited();
    			}
    		}
    		
    		void matching_to_vertex_cover()
    		{
    			leftaccessed.clear();
    			rightaccessed.clear();
    			leftfrontier.clear();
    			rightfrontier.clear();
    			leftvisited.assign(capacity,false);
			rightvisited.assign(capacity,false);
			int u,v;
			
			// Identify all vertices reacheable from an unmatched vertex via an alternating path (using breadth-first search)
    			for (size_t i=0; i<left.size; i++) {
    				if (is_matched(left[i])) continue;
    				leftfrontier.push_back(left[i]);
    				leftvisited[left[i]] = true;
    				leftaccessed[left[i]] = true;
    				do {
    					while (!leftfrontier.empty()) {
    						u = leftfrontier.back();
    						leftfrontier.pop_back();
    						for (size_t q=0; q<neighbor[u].size; q++) {
    							v = neighbor[u][q];
    							if (matching[u] == v) continue;
    							if (rightvisited[v]) continue;
    							assert(is_matched(v)); // otherwise we would have an augmenting path !
    							rightfrontier.push_back(v);
    							rightvisited[v] = true;
    							rightaccessed[v] = true;
    						}
    					}
    					while (!rightfrontier.empty()) {
						v = rightfrontier.back();
						rightfrontier.pop_back();
						int next = matching[v];
						if ((next == -1) || (leftvisited[next])) continue;
						leftfrontier.push_back(next);
						leftvisited[next] = true;
						leftaccessed[next] = true;
					}
				} while (!leftfrontier.empty());
				clear_visited();
			}
			
			// Fill the vertex cover according to Koenig's theorem
			for (size_t i=0; i<left.size; i++) {
				if (!leftaccessed[left[i]]) vcover.add(left[i]);
			}
			for (size_t i=0; i<right.size; i++) {
				if (rightaccessed[right[i]]) vcover.add(right[i]);
			}
				
    		}
    		
		std::ostream& display(std::ostream& os) const {
			os << "left: " << left << "\nright: " << right << std::endl; 
			Graph::display(os);
			return os;
		}

	};


	std::ostream& operator<< (std::ostream& os, const BipartiteGraph& x);

	std::ostream& operator<< (std::ostream& os, const BipartiteGraph* x);
	
	
	class ReversibleGraph : public Graph {
	
	public:
		
		//int num_edges;
				
		int min_degree;
		int max_degree;
		
		IntStack *node_of_degree;
		
		int *original_size;
		
		
		ReversibleGraph() {}
	
		ReversibleGraph(const int n) {
			initialise(n, true, false);
			initialise_degree();
		}
	
		ReversibleGraph(const Graph& g) : Graph(g) {	
			initialise_degree();
		}
		
		virtual ~ReversibleGraph() {
			delete [] node_of_degree;
			delete [] original_size;
		}
	
		virtual void initialise_degree() {
			
			//std::cout << "INIT DEGREE " << (capacity-1) << std::endl;
			
			node_of_degree = new IntStack[capacity];
			for(int i=0; i<capacity; ++i) {
				node_of_degree[i].initialise(0, capacity-1, capacity, false);
			}
			min_degree = capacity-1;
			max_degree = 0;
			//num_edges = 0;
			original_size = new int[capacity];
			for(int i=0; i<capacity; ++i) {
				int d = degree(i);
				original_size[i] = d;
								
				assert(d<=capacity-1);
				
				//num_edges+=d;
				if(d < min_degree) min_degree = d;
				if(d > max_degree) max_degree = d;
				node_of_degree[d].add(i);
			}
			
			//num_edges /= 2;
			
			
			
			verify("init");
		}
		
		
		void verify(const char* msg)
		{
			// check that edges are in both directions
			for(size_t x=0; x<node.size; ++x) {
				int y = node[x];
				
#ifdef _MAINTAIN_MATRIX			
				if(neighbor[y].size != matrix[y].size()) {
					std::cout << msg << ": incosistency in representations: " << neighbor[y] << " != " << matrix[y] << std::endl;
					exit(1);
				}
#endif
				
				for(size_t j=0; j<neighbor[y].size; ++j) {
					if(!matrix[y].contain(neighbor[y][j])) {
						std::cout << msg << ": " << "matrix[" << y << "] does not contain " << neighbor[y][j]
							<< " (" << matrix[y] << ")" << std::endl;
						exit(1);
					}
					
					if(!neighbor[neighbor[y][j]].contain(y)) {
						std::cout << msg << ": " << "neighbor[" << neighbor[y][j] << "] does not contain " << y 
							<< " (" << neighbor[neighbor[y][j]] << ")" << std::endl;
						exit(1);
					}
				}
				
				// check that the vertices in [neighbor[y].size, original_size[y][ iif they are not in the graph 
				
				for(int j=neighbor[y].size; j<original_size[y]; ++j) {
					if(node.contain(neighbor[y][j])) {
						std::cout << msg << ": wrong dual! (" << neighbor[y][j] 
							<< " should be at or after position " << original_size[y] << ")" << std::endl;
						exit(1);
					}
				}
				
				for(int j=original_size[y]; j<capacity; ++j) {
					if(!node.contain(neighbor[y][j])) {
						std::cout << msg << ": wrong dual! (" << neighbor[y][j] 
							<< " should be before position " << original_size[y] << ")" << std::endl;
						exit(1);
					}
				}
			}
			
			
			// check min/max degree
			int lb=capacity;
			int ub=0;
			int ne=0;
			for(size_t x=0; x<node.size; ++x) {
				int d = degree(node[x]);
				if(d < lb) lb = d;
				if(d > ub) ub = d;
				
				ne += d;
				
				//assert(node_of_degree[d].contain(x));
				if(!node_of_degree[d].contain(node[x])) {
					std::cout << *this << std::endl << msg << ": node_of_degree[" << d << "] does not contain " 
						<< node[x] << " (" << node_of_degree[d] << ")" << std::endl;
					exit(1);
				}
			}
			
			if(lb!=min_degree) {
				std::cout << *this << std::endl << msg << ": min_degree=" << min_degree << " (should be " << lb << ")" << std::endl;
				exit(1);
			}
			if(ub!=max_degree) {
				std::cout << *this << std::endl << msg << ": max_degree=" << max_degree << " (should be " << ub << ")" << std::endl;
				exit(1);
			}
			assert(lb==min_degree);
			assert(ub==max_degree);
			ne /= 2;
			
			
			if(ne!=num_edges) {
				std::cout << *this << std::endl << msg << ": num_edges=" << num_edges << " (should be " << ne << ")" << std::endl;
				exit(1);
			}
			assert(ne==num_edges);
			
			
			// check degree ordering
			for(int d=min_degree; d<=max_degree; ++d) {
				for(size_t i=0; i<node_of_degree[d].size; ++i) {
					if(degree(node_of_degree[d][i]) != d) {
						std::cout << *this << std::endl << msg << ": degree of " << node_of_degree[d][i] << " is not " << d << std::endl;
						exit(1);
					}
					
					assert(degree(node_of_degree[d][i]) == d); 
				}
			}
					
			
		}
		
		
		

		void add_and_update(const int x) {
			// assert(x>=0);
			// assert(x<capacity);
			// assert(!node.contain(x));

#ifdef _DEBUG_UPDATE
			std::cout << std::endl << *this << "\nadd node " << x << " " ;
			neighbor[x].display(std::cout, original_size[x]);
			std::cout << std::endl;
#endif			
			
			node.add(x);
			node_set.add(x);
			int i = degree(x), y, d;
			node_of_degree[i].add(x);
			num_edges += i;
			
			while(i--) {
				y = neighbor[x][i];
				d = degree(y);
				node_of_degree[d].remove(y);
				node_of_degree[d+1].add(y);
				
				if(max_degree==d)
					max_degree = d+1;
				if(min_degree==d && node_of_degree[d].empty())
					min_degree = d+1;
				
				neighbor[y].add(x);
				
#ifdef _MAINTAIN_MATRIX			
				matrix[y].fast_add(x);
#endif
			}
			
			d = degree(x);
			if(max_degree < d) max_degree = d;
			if(min_degree > d) min_degree = d;
			
#ifdef _MAINTAIN_DUAL			
			for(i=original_size[x]; i<capacity-1; ++i) {
				y = neighbor[x][i]; // neighbor in the dual
				
#ifdef _DEBUG_UPDATE
				std::cout << " - " << y << ": ";
				neighbor[y].display(std::cout,original_size[y]);
				std::cout << " => " ;
#endif
									
				neighbor[y].add(x,original_size[y]);

#ifdef _DEBUG_UPDATE
				neighbor[y].display(std::cout,original_size[y]);
				std::cout << std::endl;
#endif
			
			}
#endif			
			
			
			//verify("add");
		}

		
		void rem_and_update(const int x) {
			// assert(x>=0);
			// assert(x<capacity);
			// assert(node.contain(x));
			
#ifdef _DEBUG_UPDATE			
			std::cout << std::endl << *this << "\nremove node " << x << " " ;
			neighbor[x].display(std::cout, original_size[x]);
			std::cout << std::endl;
#endif			
			
			node.remove(x);
			node_set.remove(x);
			int i = degree(x), y, d;
			node_of_degree[i].remove(x);
			num_edges -= i;
			
			while(i--) {
				y = neighbor[x][i];
				d = degree(y);
				
				
				// std::cout << "remove " << x << " from neighbor[" << y << "] => move " << y << " from "
				// 	<< node_of_degree[d] << " to " << node_of_degree[d-1] << std::endl;
				//			
				node_of_degree[d].remove(y);
				node_of_degree[d-1].add(y);
				
				// std::cout << "  ==> " << node_of_degree[d] << " and " << node_of_degree[d-1] << " ("
				// 	<< node_of_degree[d].empty() << ")" << std::endl;
				//			
				if(min_degree==d)
					min_degree = d-1;
				
				if(max_degree==d && node_of_degree[d].empty())
					max_degree = d-1;
				
				neighbor[y].remove(x);
				
#ifdef _MAINTAIN_MATRIX			
				matrix[y].fast_remove(x);
#endif
				
			}
			
			d = degree(x);
			if(max_degree==d) {
				while(d>=min_degree && node_of_degree[d].empty()) --d;
				max_degree = d;
			}
			if(min_degree==d) {
				while(d<=max_degree && node_of_degree[d].empty()) ++d;
				min_degree = d;
			}
			
#ifdef _MAINTAIN_DUAL
			for(i=original_size[x]; i<capacity-1; ++i) {
				y = neighbor[x][i]; // neighbor in the dual

#ifdef _DEBUG_UPDATE				
				std::cout << " - " << y << ": ";
				neighbor[y].display(std::cout,original_size[y]);
				std::cout << " => " ;
#endif
									
				neighbor[y].remove(x,original_size[y]);

#ifdef _DEBUG_UPDATE				
				neighbor[y].display(std::cout,original_size[y]);
				std::cout << std::endl;
#endif
				
			}
#endif
					
			//verify("rem");
		}
	
		std::ostream& display(std::ostream& os) const {
			Graph::display(os);
			int d = max_degree;
			
			while(d >= min_degree) {
				if(!node_of_degree[d].empty()) {
					os << d << ": [" << node_of_degree[d][0];
					for(size_t i=1; i<node_of_degree[d].size; ++i) {
						os << " " << node_of_degree[d][i];
					}
					os << "] ";
				}
				--d;
			}
			return os;
		}
	
	};
	
	
	std::ostream& operator<< (std::ostream& os, const ReversibleGraph& x);

	std::ostream& operator<< (std::ostream& os, const ReversibleGraph* x);
	
	
	
	template<class GraphType>
	class VCAlgo {
	
	public:
			
		int init_level;
			
		GraphType             _graph;
		int             _lower_bound;
		int             _upper_bound;
		Vector<int>           _trail;
		Vector<int>        _decision;
		IntStack              _cover;
		IntStack         _best_cover;
		int            _last_checked;
		
		// greedy bound 
		int * degree_count;
		
		// explanation stuff
		int*                   _rank;
		int*                  _level;
		IntList         _explanation;
	
		Vector<int>          _nogood;
		int*                 _reason;
		
		
		
		// stats
		double time_start;
		int nb_nodes;
		int nb_restarts;
		int nb_solutions;
		int nb_conflicts;
		long int nogood_size;
		long int backjump_size;
		bool optimal;
		
		// params
		int verbose;
		int fail_limit;
		bool use_backjump;
		bool use_cliques;
		bool use_buss;
		int nbtry;
		
		
		// util
		BitSet util_set;
		bool first_propag;
		bool incremental;
		
	
		VCAlgo(const GraphType& g, const int k)
			: _graph(g), _upper_bound(k)
		{
			
			int n = _graph.capacity;
			
			init_level = 0;
			_last_checked = 0;
			
			nb_nodes = 0;
			nb_restarts = -1;
			nb_solutions = 0;
			nb_conflicts = 0;
			nogood_size = 0;
			backjump_size = 0;
			optimal = false;
			
			degree_count = new int[n];
			_cover.initialise(0, n-1, n, false);
			_best_cover.initialise(0, n-1, n, false);

			_rank = new int[n+1];
			_rank[n] = -1;
			_level = new int[n+1];
			_level[n] = -1;
			_explanation.initialise(n);
			_reason = new int[n];
			
			verbose = true;
			fail_limit = -1;
			use_backjump = false;
			use_cliques = true;
			use_buss = true;
			nbtry = 10;
			
			util_set.initialise(0, n-1, BitSet::empt);
			first_propag = true;
			incremental = true;
		}

		virtual ~VCAlgo()
		{
			delete [] degree_count;
			delete [] _rank;
			delete [] _level;
			delete [] _reason;
		}
		
		double cpu_time() {
			struct rusage ru;
			getrusage(RUSAGE_SELF, &ru);
			return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000; 
		}
		
		void verify_trail(const char* msg) {
		
			int x;
			int *reason;
			int rsize;
			int dec_idx = init_level;
			for(int i=0; i<_trail.size; ++i) {
				x = _trail[i];
				
				assert(_rank[x]==i);
				if(i) assert(_level[x]>=_level[_trail[i-1]]);
				
				if(dec_idx<_decision.size && x==_decision[dec_idx]) {
					assert(_reason[x]==NO_REASON);
					if(i) assert(_level[x]==_level[_trail[i-1]]+1);
					else assert(_level[x]<=init_level+1);
					++dec_idx; 
				} else {
					
					if(_reason[x]==NO_REASON) {
						std::cout << "ERROR: reason/trail of " << x << " (" << msg << ")"<< std::endl;
						print_trail(true);
					}
						
					assert(_reason[x]!=NO_REASON);
					if(_reason[x]!=DEDUCTION) {
						
						//std::cout << " (" << msg << ") reason[" << x << "] = " << _reason[x] << std::endl;
						
						if(_reason[x]<0) {
							std::cout << " (" << msg << ") reason[" << x << "] = " << _reason[x] << std::endl;
						}
						
						assert(_reason[x]>=0);
						
						if(_reason[x]>=_nogood.size) {
							std::cout << " (" << msg << ") reason[" << x << "] = " << _reason[x] << " nogoods = " << _nogood << std::endl;
						}
						
						assert(_reason[x]<_nogood.size);
						
						if(_reason[x]+_nogood[_reason[x]]>=_nogood.size) {
							std::cout << " (" << msg << ") reason[" << x << "] = " << _reason[x] << ", |reason[" 
								<< x << "]| = " << _nogood[_reason[x]] << " nogoods = " << _nogood << std::endl;
						}
						
						
						assert(_reason[x]+_nogood[_reason[x]]<_nogood.size);
					}
				}
				
				if(_reason[x]>=0) {
					reason = &(_nogood[_reason[x]+1]);
					rsize = _nogood[_reason[x]];
					for(int j=1; j<rsize; ++j) {
						if(_rank[reason[j-1]] <= _rank[reason[j]]) {
							std::cout << "ERROR: unsorted nogood (" << msg << ")"<< std::endl;
							std::cout << "literal " << x << " clause " << _reason[x] << ":";
							for(int k=0; k<rsize; ++k) {
								std::cout << " " << reason[k];	
							}
							std::cout << std::endl;
						}
						
						assert(_rank[reason[j-1]] > _rank[reason[j]]);
					}	
				}
			}
			
			for(int i=0; i<_graph.node.size; ++i) {
				int x = _graph.node[i];
				assert(!_cover.contain(x));
				for(int j=0; j<_graph.neighbor[x].size; ++j) {
					int y = _graph.neighbor[x][j];
					
					if(_cover.contain(y)) {
						std::cout << msg << " inconsistency graph/cover: " << y << " cannot be both in the graph and in the cover " << _cover 
							<< "\n" << _graph << std::endl;
						
					}
					
					assert(!_cover.contain(y));
				}
				//for(int j=_graph.neighbor[x].size; j<_graph.original_size[x]; ++j) {
				for(int j=0; j<_graph.codegree(x); ++j) {
					int y = _graph.co_neighbor[x][j];
					
					if(!_cover.contain(y)) {
						std::cout << msg << " inconsistency: since " << x << " is in the graph, " << y << " must be either in the graph or in the cover " << _cover 
							<< "\nco_neighbor["<< x << "] = " << _graph.co_neighbor[x] << "\n" 
								<< _graph << std::endl;
						
					}
					
					assert(_cover.contain(y));
				}
			}
			
		}
		
		void print_trail(const int full=false) { 
			if(nb_nodes>1000) exit(1);
			
			std::cout << std::setw(5) << nb_nodes << ":"; 
			
			int kthdec = _decision.size-1;
			int kthcov = _cover.size-1;
			int kthnod;
			for(kthnod=_trail.size-1; kthnod>=0; --kthnod) {
				if(kthdec>=0 &&_trail[kthnod] == _decision[kthdec]) {
					--kthdec;
					if(kthcov>=0 && _trail[kthnod] == _cover[kthcov]) {
						std::cout << " [" << _trail[kthnod] << "]";
						--kthcov;
					} else  {
						std::cout << " [-" << _trail[kthnod] << "]";
					}
						
				} else if(kthcov>=0 && _trail[kthnod] == _cover[kthcov]) {
					std::cout << " +" << _trail[kthnod];
					--kthcov;
				} else {
					std::cout << " -" << _trail[kthnod];
				}
			}
			std::cout << std::endl;
			
			if(full){
				
				std::cout << "\n                  Decisions: " << _decision << std::endl;
				std::cout << " Analyze conflict with trail " ; //_cover << std::endl;
				for(int i=_trail.size-1; i>=0; --i)
					std::cout << " " << std::setw(2) << _trail[i];
				std::cout << std::endl << "                      Level: " ;
				for(int i=_trail.size-1; i>=0; --i)
					std::cout << " " << std::setw(2) << _level[_trail[i]];
				std::cout << std::endl << "                       Rank: " ;
				for(int i=_trail.size-1; i>=0; --i)				
					std::cout << " " << std::setw(2) << _rank[_trail[i]];
				std::cout << std::endl << "                     Reason: " ;
				for(int i=_trail.size-1; i>=0; --i)								
#ifdef _STORE_NOGOOD	
					std::cout << " " << (_reason[_trail[i]] ? "  " : " *");
#else
				std::cout << " " << std::setw(2) << _reason[_trail[i]];
#endif
				std::cout << std::endl;
				
			}
			
		}
		
		void verify_cover() {
			util_set.clear();
			for(int i=0; i<_cover.size; ++i)
				util_set.add(_cover[i]);
			
			int n = _graph.capacity;
			for(int x=0; x<n; ++x) {
				for(int i=0; i<_graph.original_size[x]; ++i) {
					int y = _graph.neighbor[x][i];
					if(!util_set.contain(x) and !util_set.contain(y)) {
						
						std::cout << "residual graph:\n" << _graph << std::endl;
						std::cout << "The edge (" << x << ","<< y << ") is not covered by " << util_set << std::endl;
						exit(1);
					}
				}
			}	
		}
		
		int explain_with(int* beg, int* end) {
			int reason = _nogood.size;
			int sz;
			if(beg<end) {
				sz = (int)(end-beg); 
				_nogood.add(sz);
				for(int *cur=beg; cur!=end; ++cur) {
					_nogood.add(*cur);
				}
			} else {
				sz = (int)(beg-end); 
				_nogood.add(sz);
				for(int *cur=beg; cur!=end;) {
					--cur;
					_nogood.add(*cur);
				}
			}
			return reason;
		}
	
	
		int greedy_lower_bound(const int limit=-1)
		{
			_lower_bound = 0; 
			int i, k=0, covered=0, m=_graph.num_edges, l, c, stop=limit;
			
			for(int d=_graph.max_degree; d>=_graph.min_degree; --d) {
				for(int j=_graph.node_of_degree[d].size; --j>=0;) {		
					degree_count[k++] = d;
				}
			}
			
			if(stop<0) stop = k;
			--k;
		
			while(covered<m) {
				c = degree_count[_lower_bound];
				if(++_lower_bound>=stop) break;
				
				covered += c;
				
				l = k-c;	
				while(k>l && --degree_count[k]==0) --k;
			 	for(i=k; --i>l;) --degree_count[i];
			}

			_lower_bound += _cover.size;
			return _lower_bound;
		}
		
		
    // push on the trail			
		inline void save(const int x, const int reason=NO_REASON) {
			_rank[x] = _trail.size;
			_trail.add(x);
			_reason[x] = reason;
			_level[x] = _decision.size;
		}
		
		// pop from the trail
		inline int restore() {
			return _trail.pop();						
		}
		
			
		inline void add_to_cover(const int x, const int reason=NO_REASON) {
			save(x, reason);
			_cover.add(x);
		}
			
		inline void exclude_from_cover(const int x, const int reason=NO_REASON) {
			save(x, reason);
			int y;
			
#ifdef _UNOGOOD_
			int reason_for_y = _nogood.size;
			_nogood.add(1);
			_nogood.add(x);
#else
			int reason_for_y = reason;
			if(reason==NO_REASON && use_backjump) {
				reason_for_y = _nogood.size;
				_nogood.add(1);
				_nogood.add(x);
			}
#endif
			
			for(int i=_graph.neighbor[x].size; --i>=0;) {	
				y = _graph.neighbor[x][i];
				
#ifdef _DEBUG_REDUCTION			
				if(_DEBUG_REDUCTION) {
					std::cout << "add " << y << " b/c " << x << std::endl;
				}
#endif
				
// #ifdef _UNOGOOD_
// 				add_to_cover(y, reason_for_y);
// #else
// 				add_to_cover(y, reason);
// #endif
				add_to_cover(y, reason_for_y);
				//
				_graph.rem_and_update(y);
			}
		}


		inline int undo() {
			
			// print_trail();
			// std::cout << _cover << std::endl;
			
			_graph.clear_cliques();

			int x = restore(), y = _decision.pop();
			int l = _level[y];
			
			
#ifdef _TRACE_
				std::cout << "undo: level " << l << std::endl;
#endif
			
			
			while(x!=y) {
				_graph.add_and_update(x);
				x = restore();
			}

			do {
				x = _cover.pop();	
			} while(!_cover.empty() && _level[_cover.back()]>=l);

			return (x==y ? IN_COVER(y) : IN_STABLE(y));
		}
		
		
		int kernelize_cliques() {
			int reduction = false;
			
#ifdef _TRACE_
			int sizebf = _graph.size();
			int startt, endt;
#endif
			
			if(first_propag || !incremental) {
				
				
				reduction = _graph.kernelize_cliques(*this);
				
#ifdef _TRACE_
				std::cout << "first propag: " << (sizebf - _graph.size()) << std::endl;
#endif
				
				first_propag = false;
			} else {
				
#ifdef _TRACE_
				startt = _last_checked;	
				endt = 	_trail.size;
#endif
				
				reduction = _graph.kernelize_cliques_incremental(*this);
				
#ifdef _TRACE_
				std::cout << "propag: " << (sizebf - _graph.size()) << " |";
				for(int i=startt; i<endt; ++i) {
					std::cout << " " << _trail[i];
				}
				std::cout << std::endl;
#endif
				
			}
			return reduction;
		}
		
		
		inline bool suboptimal() {
			_lower_bound = _cover.size;
			if(_graph.num_edges) {

				int handled = (2*use_buss + use_cliques);
				int potential_reduction = handled;

				do {
					if(use_cliques && potential_reduction&1) {
						potential_reduction |= 2*(kernelize_cliques());
						potential_reduction ^= 1;
					}
					if(use_buss && potential_reduction&2) {
						potential_reduction |= _graph.kernelize_buss(*this);
						potential_reduction ^= 2;
					}		
				} while(potential_reduction&handled);

				if(_upper_bound<_graph.capacity) {
					_lower_bound = _cover.size + (_graph.size() - _graph.min_clique_cover_bitset());
				}
			}
			
#ifdef _TRACE_
			std::cout << "bounds = [" << _lower_bound << ".." << _upper_bound << "]" << std::endl;
#endif
			
			return _lower_bound >= _upper_bound;
		}

		
		
		bool backjump() {
			
#ifdef _VERIF_TRAIL
			verify_trail("before backjump");
#endif			

#ifdef _DEBUG_BACKJUMPS
			if(_DEBUG_BACKJUMPS) {
				print_trail(true);
			}
#endif			
			
			++nb_conflicts;	
			if(_decision.size<=init_level) {
				return false;
			}

			// first, compute a minimal explanation for the last fail
			int i, j, x, y, z, u;
			_explanation.quick_clear();
			
			
			
			// std::cout << std::endl << "cover: " << _cover << std::endl;
			// std::cout << "nodes: " << _graph.node << std::endl << "dummy: ";
			// for(i=_cover.size; i<_graph.capacity; ++i) {
			// 	//x = _cover[i];
			// 	if(!_graph.node.contain(_cover[i])) {
			//
			// 		//std::cout << "i: " << i << ", x: " << x << std::endl;
			//
			// 		std::cout << " " << _cover[i];
			// 	}
			// }
			// std::cout << std::endl;
			// std::cout << "nc: " << _graph.num_cliques << " / " << _graph.capacity << std::endl;
			
			for(i=_cover.size; i<_graph.capacity; ++i) {
				x = _cover[i];
				if(!_graph.node.contain(x)) {
					
					//std::cout << "i: " << i << ", x: " << x << "nc: " << _graph.num_cliques << " / " << _graph.capacity << std::endl;
					
					_graph.candidates[_graph.num_cliques].copy(_graph.matrix[x]);
					++_graph.num_cliques;	
				}	
			}
			
			//std::cout << "nc: " << _graph.num_cliques << " / " << _graph.capacity << std::endl << std::endl;

			int slack = _lower_bound-_upper_bound;
			//if(slack>0) std::cout << "slack = " << slack << std::endl;
			
			for(i=_cover.size-1; i>=0; --i) {
				x = _cover[i];
					
				for(j=0; j<_graph.num_cliques; ++j) {
					if(_graph.candidates[j].contain(x)) {
						_graph.candidates[j].intersect_with(_graph.matrix[x]);
						break;
					}
				}
				if(j>=_graph.num_cliques) {
					_explanation.add_end(x);		
				} else if(slack>0) {
					_graph.candidates[_graph.num_cliques].copy(_graph.matrix[x]);
					++_graph.num_cliques;
					--slack;
				}
			}
			

#ifdef _DEBUG_BACKJUMPS
			if(_DEBUG_BACKJUMPS) {
				if(_cover.empty()) {
					std::cout << "Extract contradiction!" << std::endl;
				} else {
					std::cout << " Extract a minimal conflict " << _explanation << std::endl
										<< "                 from cover [" << _cover[_cover.size-1]; //_cover << std::endl;
					for(int k=_cover.size-2; k>=0; --k) {
						std::cout << " " << _cover[k];
					}
					std::cout << "]" << std::endl;
				}
			}
#endif
			
			// now use resolution to remove all but one of this level's literals
			if(_explanation.size == 0) {
				
#ifdef _NOGOOD_STAT
				++backjump_sizes[_decision.size];
				++nogood_sizes[0];
#endif
								
				backjump_size += _decision.size;
				return false;
			} else {
				// resolution
	
				x = _explanation.get_first();
				y = _explanation.next[x];
				

				int new_nogood = NO_REASON;
				int lst_reason = NO_REASON;
				
				int nbacktracks = _decision.size;
				
				int min_reason = _nogood.size;
				
				if(y != _explanation.capacity) {
					while(_level[y]==_level[x]) {
						// resolve x since there are at least one other literal from the same level
							
						//std::cout << "           current conflict " << _explanation << std::endl;
							
						z = y;				
						_explanation.remove(x);
						

#ifdef _DEBUG_BACKJUMPS							
						if(_reason[x] < 0) {
							std::cout << "node " << x << " level " << _level[x] << " rank " << _rank[x] << std::endl;
							std::cout << nb_nodes << std::endl;
						}
						
						assert(_reason[x]>=0);
						assert(_reason[x]<_nogood.size);
#endif
						
						// the previous resolved literal had the same reason, no need to resolve twice		
						if(_reason[x] != lst_reason) {
							lst_reason = _reason[x];
							
							int* reason = _nogood.stack_+_reason[x]+1;
							int   rsize	= _nogood[_reason[x]];
						
							min_reason = _reason[x];
						
						
#ifdef _DEBUG_BACKJUMPS
						
							for(i=1; i<rsize; ++i) {
								assert(_rank[reason[i-1]] > _rank[reason[i]]);
							}
						
						
							if(_DEBUG_BACKJUMPS) {
								std::cout << " resolve " << x << " with its explanation" ;
								for(int k=0; k<rsize; ++k)
									std::cout << " " << reason[k];
								std::cout << std::endl;
							}
#endif	

							for(i=0; i<rsize; ++i) {
								u = reason[i];
								while(_rank[z]>_rank[u]) {
									z = _explanation.next[z];
								}
								if(_rank[z]<_rank[u]) { // otherwise z and u are equal
									_explanation.add_before(u,z);
								} 
							}
						
						}
						// saved resolution
	
						x = _explanation.get_first();
						y = _explanation.next[x];
						
						if(_level[x] == 0) {
							
#ifdef _NOGOOD_STAT
							++backjump_sizes[_decision.size];
							++nogood_sizes[0];
#endif
							
							backjump_size += _decision.size;
							return false;
						} 
						
						//else std::cout << "level["<<x<<"] = " << _level[x] << std::endl;
						
					}

				
					// x is the last literal from this level and y is the literal of highest level otherwise
					_explanation.remove(x);

					_nogood.size = min_reason;
					
#ifdef _DEBUG_BACKJUMPS
					if(_DEBUG_BACKJUMPS) {	
						std::cout << " new backjump nogood:" ;
					}
#endif				
															
					new_nogood = _nogood.size;
					_nogood.add(_explanation.size);
					int l = _explanation.get_first();
					do {
						_nogood.push_back(l);
						l = _explanation.next[l];
					} while(l < _explanation.capacity);
					
					if(y<_explanation.capacity)
						nbacktracks -= _level[y];
					
				}	
				
						if(_level[x] == 0) {
							
#ifdef _NOGOOD_STAT
							++backjump_sizes[_decision.size];
							++nogood_sizes[0];
#endif
							
							backjump_size += _decision.size;
							return false;
						}
				
#ifdef _DEBUG_BACKJUMPS
				if(_DEBUG_BACKJUMPS) {
					std::cout << "                            " << _explanation << std::endl;
					std::cout << " backjump " << nbacktracks << " levels (of " << _decision.size << ") and exclude " << x << std::endl << std::endl;
				}
#endif
				
				
				nogood_size += _explanation.size;
				backjump_size += nbacktracks;
#ifdef _NOGOOD_STAT
				++backjump_sizes[nbacktracks];
				++nogood_sizes[_explanation.size];
#endif
				
				bool was_in_cover = _cover.contain(x);
	
				y = undo();
				
				while(--nbacktracks>0) {

					z = NODE(y);							

					_graph.add_and_update(z);		

					y = undo();	

				}	
				
				_last_checked = _trail.size;
				
				z = NODE(y);

				if(z!=x) {
			
					_graph.add_and_update(z);
	
					_graph.rem_and_update(x);
				}
				if(was_in_cover){
					exclude_from_cover(x, new_nogood);
				} else {
					add_to_cover(x, new_nogood);
				}
				
			}
			
#ifdef _TRACE_
			print_trail();
#endif
			
#ifdef _VERIF_TRAIL					
				verify_trail("after backjump");
#endif			

			return true;
		}
		
		
		bool backtrack() {
			++nb_conflicts;
			if(_decision.size>init_level) {
				++backjump_size;
				int e = undo();
				_last_checked = _trail.size;
				
				if(STATUS(e)==STA) {
					add_to_cover(NODE(e), DEDUCTION);
				} else {
					exclude_from_cover(NODE(e), DEDUCTION);
				}
#ifdef _TRACE_
			print_trail();
#endif
			
#ifdef _VERIF_TRAIL
			verify_trail("after branch");
#endif
			
				return true;
			}
			return false;
		}
			
		inline void branch_on_maxdegree() {
			int x = _graph.node_of_degree[_graph.max_degree].back();				
			_decision.add(x);
			_graph.rem_and_update(x);
			add_to_cover(x, NO_REASON);
				
#ifdef _TRACE_
			print_trail();
#endif
			
#ifdef _VERIF_TRAIL
			verify_trail("after branch");
#endif
				
		}
					
		inline void reset() {
			int e=0;
			while(_decision.size>init_level) {
				e = undo();
				_graph.add_and_update(NODE(e));
				//_graph.add_and_update((x < 0 ? -x-1 : x));
		  }
		}
			
		void print_head() {
			std::cout << std::setw(3) << "#" << " | "  
								<< std::setw(6) << "cover" << " " 
								<< std::setw(6) << "IS" << " " 
								<< std::setw(8) << "time" << " "
								<< std::setw(5) << "rst" << " "
								<< std::setw(10) << "fails" << " "
								<< std::setw(10) << "nodes" << " "
								<< std::setw(6) << "|ngd|" << " "
								<< std::setw(6) << "bkjp" << " "
								<< std::setw(12) << "memory" << " "
								<< std::endl
								<< std::setfill('=') << std::setw(83) << "="
							  << std::endl << std::setfill(' ');
		}


		void print_summary() {
			double time_now = cpu_time();
			std::cout << std::setfill('=') << std::setw(83) << "=" << std::endl << std::setfill(' ');
			print_head();
			print_statistics();
			std::cout << std::setfill('=') << std::setw(83) << "=" << std::endl << std::setfill(' ')
								<< std::setw(12) << (optimal ? "optimal " : "timeout ")		
								<< std::setw(9) << (int)(nb_conflicts/(time_now-time_start)) << " fails/s  "
								<< std::setw(9) << (int)(nb_nodes/(time_now-time_start)) << " nodes/s  "
								<< "       memory: " << std::setw(9) << std::left << nogood_size << std::right << std::endl;
								//<< std::setfill('=') << std::setw(62) << "=" << std::endl << std::setfill(' ');
	
		#ifdef _NOGOOD_STAT
			for(int i=0; i<_graph.capacity; ++i) {
				if(backjump_sizes[i] || nogood_sizes[i]) {
					std::cout << std::right << std::setw(8) << nogood_sizes[i] << " nogoods," 
						<< std::right << std::setw(8) << backjump_sizes[i] << " backjumps of size " 
							<< std::left << std::setw(4) << i << std::endl;
				}
			}
		#endif
	
		}


		void print_statistics() {
			std::cout << std::setw(3) << nb_solutions << " | "  
								<< std::setw(6) << _upper_bound << " " 
								<< std::setw(6) << (_graph.capacity-_upper_bound) << " " 
								<< std::setw(8) << (int)((cpu_time()-time_start)*1000) << " "
								<< std::setw(5) << nb_restarts << " "
								<< std::setw(10) << nb_conflicts << " "
								<< std::setw(10) << nb_nodes << " "
								<< std::setw(6) << (int)(nb_conflicts ? (double)nogood_size/(double)nb_conflicts : 0) << " "
								<< std::setw(6) << (double)((int)(((nb_conflicts ? (double)backjump_size/(double)nb_conflicts : 0)*10)))/10.0 << " "
								<< std::setw(12) << 
		#ifdef _STORE_NOGOOD
									nogood_size
		#else
									_nogood.capacity
		#endif
								<< " "
								<< std::endl;
		}


		void new_cover() {
			++nb_solutions;
			_upper_bound = _cover.size;		
			
			
			//if(_cover.size < _algo.caapcity/2)
			_best_cover.clear();
			for(int i=0; i<_cover.size; ++i) {
				_best_cover.add(_cover[i]);
			}
			
			
			
#ifdef _TRACE_
			std::cout << " new cover of size " << _upper_bound << std::endl;
#endif
			
			if(verbose) print_statistics();
			verify_cover();
		}
		
		bool run(const int limit, const int ub_limit)
		{			
			init_level = _decision.size;
	
			int nb_restart_conflicts = 0, feasible=true;			
			while(feasible && (limit<0 || nb_restart_conflicts<limit) && _upper_bound>ub_limit) {
				++nb_nodes;
				if(suboptimal()) {
					++nb_restart_conflicts;
					feasible = (use_backjump ? backjump() : backtrack());
				} else {
					if(_graph.num_edges==0) {
						new_cover();
						++nb_restart_conflicts;
					} else {
						branch_on_maxdegree();
					}	
				}
			}
	
			reset();
			
			if(!feasible) optimal = true;
			return feasible;
		}
			
		void solve()
		{		
			time_start = cpu_time();

			if(verbose) print_head();
	
			run(fail_limit, -1);
			//run(300, -1);
		}

	
	};


	class ReversibleCompactGraph : public CompactGraph {
	
	public:
		int min_degree;
		int max_degree;
		
		Vector<int>    *co_neighbor;
		
		Vector<int> *node_of_degree;
		int           *degree_index;
		
		int *original_size;
		
		
		int num_cliques;
		
#ifdef _NOT_SO_COMPACT		
		BitSet  *candidates;
#endif		
		
		
		// // for each vertex the list of vertices for which it is a "watcher"
		// Vector<int> *watched_by;
		// 	  // for each watched vertex, ist index in the other watcher list,
		// // so that when we remove it from one list, we can also remove it from the second
		// Vector<int> *watched_index;
		// Vector<int> *watched_other;
		

		TwoWatchedThing<int> watchers;

		
		
		int* saturation;
		//Vector<int> max_satur;
		
		
		IntStack util_stack;
		Vector<int> util_vec;
		BitSet util_set;
		
#ifdef _COUNT_OP
		int count;
#endif		
		
		
		ReversibleCompactGraph() {
	
#ifdef _COUNT_OP
		count = 0;
#endif			
			
		}
	
		ReversibleCompactGraph(const int n) {
			
#ifdef _COUNT_OP
		count = 0;
#endif	
		
			CompactGraph::initialise(n);
			initialise_degree();
			//initialise_watchers();
		}
	
		ReversibleCompactGraph(const CompactGraph& g) : CompactGraph(g) {	
			
#ifdef _COUNT_OP
		count = 0;
#endif	
		
			initialise_degree();
			//initialise_watchers();
		}
		
		ReversibleCompactGraph(const ReversibleCompactGraph& g) : CompactGraph(g) {	
			initialise_memory();
			// co_neighbor = new Vector<int>[capacity];
			// node_of_degree = new Vector<int>[capacity];
			// watched_by = new Vector<int>[capacity];
			// watched_index = new Vector<int>[capacity];
			// degree_index = new int[capacity];
			min_degree = g.min_degree;
			max_degree = g.max_degree;
			original_size = new int[capacity];
			for(int i=0; i<capacity; ++i) {
				original_size[i] = g.original_size[i];
				degree_index[i]	= g.degree_index[i];	
				node_of_degree[degree(i)].add(i);
			}
			// util_stack.initialise(0, capacity-1, capacity, false);
			// saturation = new int[capacity];
			
			//verify("init");
		}
		
		virtual ~ReversibleCompactGraph() {
			delete []    co_neighbor;
			delete [] node_of_degree;
			delete []  original_size;
			delete []   degree_index;
			// delete []     watched_by;
			// delete []  watched_index;
			// delete []  watched_other;
			delete []     saturation;
			
#ifdef _NOT_SO_COMPACT
			delete []     candidates;
#endif
			
		}
	
		void initialise_memory() {

			num_cliques = 0;
			// watched_by = new Vector<int>[capacity];
			// watched_index = new Vector<int>[capacity];
			// watched_other = new Vector<int>[capacity];
			co_neighbor = new Vector<int>[capacity];
			node_of_degree = new Vector<int>[capacity];
			degree_index = new int[capacity];
			util_stack.initialise(0, capacity-1, capacity, false);
			util_set.initialise(0, capacity-1, BitSet::empt);
			saturation = new int[capacity];
			watchers.initialise(capacity);
					
#ifdef _NOT_SO_COMPACT
			candidates = new BitSet[capacity];
			for(int i=0; i<capacity; ++i) {
				candidates[i].initialise(0, capacity-1, BitSet::full);
			}
#endif
			
		}
	
		void initialise_degree() {
			initialise_memory();
			min_degree = capacity-1;
			max_degree = 0;
			original_size = new int[capacity];
			for(int i=0; i<capacity; ++i) {
				int d = degree(i);
				original_size[i] = d;
								
				assert(d<=capacity-1);

				if(d < min_degree) min_degree = d;
				if(d > max_degree) max_degree = d;
				
				degree_index[i] = node_of_degree[d].size;
				node_of_degree[d].add(i);
			}
			
			verify("init");
		}
		
		int real_size() { return node.size-node_of_degree[0].size; }
		
		
		void verify_watchers(const char* msg) {
			for(size_t i=0; i<node.size; ++i) {
				int w1 = node[i];
				
				if(watchers.watched_by[w1].size != watchers.watched_index[w1].size || watchers.watched_by[w1].size != watchers.watched_other[w1].size) {
					std::cout << msg << ": inconsistency in size of the watched lists of " << w1 << ": " << std::endl;
					std::cout << watchers.watched_by[w1] << std::endl;
					std::cout << watchers.watched_index[w1] << std::endl;
					std::cout << watchers.watched_other[w1] << std::endl;
					
					exit(1);
				}
				
				for(size_t ixw1=0; ixw1<watchers.watched_by[w1].size; ++ixw1) {
					int x = watchers.watched_by[w1][ixw1];
					int w2 = watchers.watched_other[w1][ixw1];
					int ixw2 = watchers.watched_index[w1][ixw1];
					
					if(watchers.watched_by[w2][ixw2] != x) {
						std::cout << msg << ": inconsistency in the watched lists (element) of " << x << " (" << w1 << " & " << w2 << "): " << std::endl;
						std::cout << watchers.watched_by[w1] << " / " << watchers.watched_by[w2] << std::endl;
						std::cout << watchers.watched_index[w1] << " / " << watchers.watched_index[w2] << std::endl;
						std::cout << watchers.watched_other[w1] << " / " << watchers.watched_other[w2] << std::endl;
						
						exit(1);
					}
					
					if(watchers.watched_other[w2][ixw2] != w1) {
						std::cout << msg << ":inconsistency in the watched lists (other) of " << x << " (" << w1 << " & " << w2 << "): " << std::endl;
						std::cout << watchers.watched_by[w1] << " / " << watchers.watched_by[w2] << std::endl;
						std::cout << watchers.watched_index[w1] << " / " << watchers.watched_index[w2] << std::endl;
						std::cout << watchers.watched_other[w1] << " / " << watchers.watched_other[w2] << std::endl;
						
						exit(1);
					}
					
					if(watchers.watched_index[w2][ixw2] != (int)ixw1) {
						std::cout << msg << ": inconsistency in the watched lists (index) of " << x << " (" << w1 << " & " << w2 << "): " << std::endl;
						std::cout << watchers.watched_by[w1] << " / " << watchers.watched_by[w2] << std::endl;
						std::cout << watchers.watched_index[w1] << " / " << watchers.watched_index[w2] << std::endl;
						std::cout << watchers.watched_other[w1] << " / " << watchers.watched_other[w2] << std::endl;
						
						exit(1);
					}
					
				}
			}
			
		}
		
		
		// void initialise_watchers(VCAlgo<ReversibleCompactGraph>& algo) {
		// 	watched_by = new Vector<int>[capacity];
		// 	bool reduction = clique_dominance(algo);
		// 	int start = 0;
		// 	if(reduction)
		// 		clique_dominance_incremental(algo);
		// 	// // while(reduction) {
		// 	// // 	reduction = clique_dominance_incremental();
		// 	// // }
		// 	//
		// 	// std::cout << reduction << std::endl;
		// 	//
		// 	// for(int i=0; i<node.size; ++i) {
		// 	// 	std::cout << watched_by[node[i]] << std::endl;
		// 	// }
		// }
		
		void dominance_brute_force(VCAlgo<ReversibleCompactGraph>& algo) {
			
			bool print = false;
			if(num_edges<4) {
				print =true;
			}
			
			if(print)
				std::cout << "clique dom on\n" << this << std::endl << node << std::endl;
			
			bool reduction = true;
			while(reduction)
				reduction = clique_dominance(algo);
			
			if(print)
				std::cout << "end" << std::endl;
			
		}
		
	
		
#ifdef _NOT_SO_COMPACT		
		void clear_cliques() {
			for(int i=0; i<num_cliques; ++i) {
				candidates[i].fill();
			}
			num_cliques=0;
		}
		
		
		int min_clique_cover_bitset()
		{
			int d, j, i, x;
			clear_cliques();
			for(d=min_degree; d<=max_degree; ++d) {
				for(j=node_of_degree[d].size; --j>=0;) {
					x = node_of_degree[d][j];
					
					//std::cout << "color " << x << " with" ;
					
					
					for(i=0; i<num_cliques; ++i) {
						//std::cout << " " << i << "?";
						if(candidates[i].fast_contain(x)) break;
					}
					if(i==num_cliques) {
						++num_cliques;
					}
					candidates[i].intersect_with(matrix[x]);
					
					//std::cout << std::endl;
					
				}
			}			
			return num_cliques;
		}
		
		
		int min_clique_cover_bitset2()
		{
			int d, j, i, x, n;
			for(i=0; i<num_cliques; ++i) {
				candidates[i].fill();
			}
			num_cliques=0;
			for(d=min_degree; d<=max_degree; ++d) {
				n = node_of_degree[d].size;
				for(j=0; j<n; ++j) {
					x = node_of_degree[d][j];
					for(i=0; i<num_cliques; ++i) {
						if(candidates[i].fast_contain(x)) break;
					}
					if(i==num_cliques) ++num_cliques;
					candidates[i].intersect_with(matrix[x]);
				}
			}			
			return num_cliques;
		}
		
		
		// stupi implementation of dsatur, just as a witness that it does not seem to help
		int min_clique_cover_bitset_sat(const int limit=-1)
		{
			size_t i;
			int x, y;
	
			for(int c=0; c<num_cliques; ++c) {
				candidates[c].fill();
			}
			num_cliques=0;
			
			util_stack.clear();
			for(i=0; i<node.size; ++i) {
				util_stack.add(node[i]);
			}
			
			int *saturation = new int[capacity];
			std::fill(saturation, saturation+capacity, 0);
			

			while(!util_stack.empty()) {
				x = util_stack[0];
				for(i=1; i<util_stack.size; ++i) {
					y = util_stack[i];
					//if(saturation[x]<saturation[y] || (saturation[x]==saturation[y] && degree(x)>degree(y))) {
					if(degree(x)>degree(y) || (degree(x)==degree(y) && saturation[x]<saturation[y])) {
						x = y;
					}
				}
				
				util_stack.remove(x);

				for(i=0; (int)i<num_cliques; ++i) {
					if(candidates[i].fast_contain(x)) break;
				}
				if((int)i==num_cliques) ++num_cliques;
				candidates[i].intersect_with(matrix[x]);
		
				for(i=0; (int)i<degree(x); ++i) {
					--saturation[neighbor[x][i]];
				}
			}			
			
			
			return num_cliques;
		}
		
		
		bool neighboring_clique(const int x, int& w1, int& w2) {
			int i, y;

			w1 = -1;
			w2 = -1;
			
			
			
#ifdef _DEBUG_NEIGHCLIQUE			
			if(_DEBUG_NEIGHCLIQUE) {
				std::cout << "check if " << x << "'s neighbors = " << neighbor[x] << " is a clique" << std::endl;
			}
#endif
			
			if(neighbor[x].size>1) {		
				// if(neighbor[x].size<4) {
				// 	for(int i=1; i<neighbor[x].size; ++i) {
				// 		for(int j=0; j<i; ++j) {
				// 			a = neighbor[x][i];
				// 			b = neighbor[x][j];
				// 			if(!neighbor[a].contain(b)) {
				// 				w1 = a;
				// 				w2 = b;
				// 				return false;
				// 			}
				// 		}
				// 	}
				// } else {
				util_set.copy(matrix[x]);
				util_set.intersect_with(node_set);
					
				for(i=neighbor[x].size; --i>=0;) {
					y = neighbor[x][i];
					util_set.remove(y);
						
#ifdef _DEBUG_NEIGHCLIQUE			
					if(_DEBUG_NEIGHCLIQUE) {
						std::cout << "      -> " << y << ": " << matrix[y] << " is a super set of " << util_set << "?" << std::endl;
					}
#endif
						
						
					// std::cout << util_set.empty() << std::endl;
					// std::cout << util_set.size() << std::endl;
						
					//std::cout << util_set << " \\ " << matrix[y] << std::endl;
					if(!util_set.included(matrix[y])) {
						util_set.setminus_with(matrix[y]);
						
						//std::cout << " = " << util_set << std::endl;
						
						w1 = util_set.min();
						
						//std::cout << w1 << std::endl;
						
						w2 = y;
						break;
					} // else {
 // 						std::cout << " = {}" << std::endl;
 // 					}
					util_set.add(y);
				}
				// }
			}
			
#ifdef _DEBUG_NEIGHCLIQUE			
			if(_DEBUG_NEIGHCLIQUE) {
				if(w1<0) {
					std::cout << "      -> clique!" << std::endl;
				} else {
					std::cout << "      -> not a clique b/c " << w1 << " and " << w2 << " are not neighbors" << std::endl;
				}
			}
#endif
			
			return (w1<0);
		}
		
		
		inline bool kernelize_buss(VCAlgo<ReversibleCompactGraph>& algo) {
			bool pruning = false;
			if(num_edges && max_degree >= (algo._upper_bound - (int)(algo._cover.size))) {
				

				
				int* begr = algo._cover.begin();
				int offset = std::min((int)(algo._cover.size), (max_degree - (algo._upper_bound - (int)(algo._cover.size))));
				int* endr = algo._cover.end()-offset;
				
				int x, reason; 
				do { 
					
					x = node_of_degree[max_degree].back();
					
					util_vec.clear();
					//for(int i=algo._cover.size-1; i>=0; --i) {
					for(int i=0; i<(int)(algo._cover.size)-offset; ++i) {
						if(!matrix[x].contain(algo._cover[i])) {
							util_vec.add(algo._cover[i]);
						}
					}
					begr = algo._cover.begin();
					endr = algo._cover.end();

					//std::cout << util_vec.size << " " << algo._cover.size-offset << std::endl;
					
					
#ifdef _DEBUG_REDUCTION			
					if(_DEBUG_REDUCTION) {
						std::cout << " add " << x << " (Buss rule): " << algo._cover << std::endl;
					}
#endif
					
					rem_and_update(x);
					reason = (algo.use_backjump ? algo.explain_with(endr, begr) : NO_REASON);
					algo.add_to_cover(x, reason);
				}			
				while(num_edges && max_degree >= (algo._upper_bound - (int)(algo._cover.size)));
				pruning = true;
			}
			return pruning;
		}
		
		bool kernelize_cliques(VCAlgo<ReversibleCompactGraph>& algo)
		{
			
#ifdef _DEBUG_KERNCLIQUE			
			if(_DEBUG_KERNCLIQUE) {
				std::cout << "kernelize for the clique dominance: " << node << std::endl;
			}
#endif
			
			
			int d, j, x, reason, w1, w2;
			bool fix_point = false;
			bool reduction = false;
			while(num_edges && !fix_point) {
				fix_point = true;
				//d = std::max(min_degree,2);
				d = min_degree;
				
#ifdef _DEBUG_KERNCLIQUE			
				if(_DEBUG_KERNCLIQUE) {
					std::cout << " fixed point iteration" << std::endl;
				}
#endif
				
				//while(num_edges && d<=max_degree) {
				do {
					
					while(num_edges && min_degree<=1) {
						x = node_of_degree[min_degree][0];
					// while(num_edges && node_of_degree[1].size>0) {
					// 	x = node_of_degree[1][0];
						
#ifdef _DEBUG_REDUCTION			
					if(_DEBUG_REDUCTION) {				
// #ifdef _DEBUG_KERNCLIQUE
// 						if(_DEBUG_KERNCLIQUE) {
							std::cout << "  remove " << x << " b/c low degree (" << degree(x) << "): " << co_neighbor[x] << std::endl;
						}
#endif
						
						rem_and_update(x);
						reason = (algo.use_backjump ? algo.explain_with(co_neighbor[x].end(), co_neighbor[x].begin()) : NO_REASON);
						algo.exclude_from_cover(x, reason);			
						if(!size()) return true;
						fix_point = false;
						reduction = true;
					}
					
					//if(d<2) d=2;
					
#ifdef _DEBUG_KERNCLIQUE			
			if(_DEBUG_KERNCLIQUE) {
				std::cout << "  iteration for nodes of degree " << d << ":"; //<< std::endl;
			}
#endif
					
					for(j=0; j<(int)(node_of_degree[d].size);) {
						x = node_of_degree[d][j];
						
#ifdef _DEBUG_KERNCLIQUE			
			if(_DEBUG_KERNCLIQUE) {
				std::cout << " " << x ; //<< std::endl;
			}
#endif
						
						if(neighboring_clique(x, w1, w2)) {
							
#ifdef _DEBUG_REDUCTION			
					if(_DEBUG_REDUCTION) {							
// #ifdef _DEBUG_KERNCLIQUE
// 						if(_DEBUG_KERNCLIQUE) {
							std::cout << "\n  remove " << x << " b/c its neighborhood N-" << co_neighbor[x] << "is a clique" << std::endl;
						}
#endif
							
							rem_and_update(x);
							reason = (algo.use_backjump ? algo.explain_with(co_neighbor[x].end(), co_neighbor[x].begin()) : NO_REASON);
							algo.exclude_from_cover(x, reason);							
							//if(!size()) return true;
							fix_point = (!size());//false;
							reduction = true;
						} else {
							// // store witnesses
							// watched_index[w1].add(watched_by[w2].size);
							// watched_index[w2].add(watched_by[w1].size);
							// watched_other[w1].add(w2);
							// watched_other[w2].add(w1);
							// watched_by[w1].add(x);
							// watched_by[w2].add(x);
							

#ifdef _DEBUG_KERNCLIQUE			
						if(_DEBUG_KERNCLIQUE) {
							std::cout << "\n  the neighborhood of " << x << " is not a clique b/c " << w1 << " and " << w2 << " are not neighbors" << std::endl;
						}
#endif
							
							watchers.watch(x,w1,w2);
							
#ifdef _DEBUG_KERNCLIQUE			
						if(_DEBUG_KERNCLIQUE) {
							std::cout << " watch ok "<< std::endl;
						}
#endif
							
							
#ifdef _VERIF_WATCHERS
							verify_watchers("after add");
#endif
							
							// move on
							++j;
						}
					}
					
					
#ifdef _DEBUG_KERNCLIQUE			
			if(_DEBUG_KERNCLIQUE) {
				std::cout << std::endl;
			}
#endif
					
					if(min_degree<=1) {
						d = min_degree;
					} else {
						++d;
					}
				} while(num_edges && d<=max_degree);
			}
			
			
#ifdef _VERIF_WATCHERS
			std::cout << "VERIFY AT THE END OF THE INIT PHASE!" << std::endl;
					verify_watchers("end init");
#endif
			
			return reduction;
		}
		
		
		bool kernelize_cliques_incremental(VCAlgo<ReversibleCompactGraph>& algo) {
			int x, y, w1, w2, reason, reduction=false; 
			
#ifdef _VERIF_WATCHERS
					verify_watchers("start incr");
#endif
			
			
#ifdef _DEBUG_KERNCLIQUEINCR			
			if(_DEBUG_KERNCLIQUEINCR) {
				std::cout << "\nkernelize for the clique dominance (incremental, checking the removed witnesses:";
				for(i=algo._last_checked; i<algo._trail.size; ++i) {
					std::cout << " " << algo._trail[i];
				}
				std::cout << ")" << std::endl;
			}
#endif
			
			
			while(num_edges && algo._last_checked < (int)(algo._trail.size)) {
				
				// first apply the low degree rule
				while(num_edges && min_degree<=1) {
					x = node_of_degree[min_degree][0];
				// while(num_edges && node_of_degree[1].size>0) {
				// 	x = node_of_degree[1][0];
	
#ifdef _DEBUG_REDUCTION			
					if(_DEBUG_REDUCTION) {							
// #ifdef _DEBUG_KERNCLIQUEINCR
// 					if(_DEBUG_KERNCLIQUEINCR) {
						std::cout << "  remove " << x << " b/c low degree (" << degree(x) << ") " << co_neighbor[x] << std::endl;
					}
#endif
					
					reason = (algo.use_backjump ? algo.explain_with(co_neighbor[x].end(), co_neighbor[x].begin()) : NO_REASON);
					rem_and_update(x);
					algo.exclude_from_cover(x, reason);	
					reduction = true;		
				}
				if(!size()) break;
				
				
				// now take the first removed node y and check the nodes for which it was witness 
				y = algo._trail[algo._last_checked];
				
#ifdef _DEBUG_KERNCLIQUEINCR			
				if(_DEBUG_KERNCLIQUEINCR) {
					std::cout << "  * check non cliques witnessed by " << y << ": " << watchers.watched_by[y] << " " << num_edges << " " << (watchers.watched_by[y].size) << std::endl;					
				}
#endif
				
				for(int q=watchers.watched_by[y].size-1; num_edges && q>=0; --q) {
					x = watchers.watched_by[y][q];
					// x   = watched_by[y].pop();
					// ixz = watched_index[y].pop();
					// z   = watched_other[y].pop();
					//
					// // we also need to remove x from its other witness' list
					// t = watched_other[z].pop();
					// iut = watched_index[z].pop();
					// //int u = watched_by[z].pop();
					//
					// watched_index[t][iut] = ixz;
					//
					// watched_index[z][ixz] = iut;
					// watched_other[z][ixz] = t;
					// watched_by[z][ixz] = watched_by[z].pop();
					
// #ifdef _VERIF_WATCHERS
// 					verify_watchers("after rem (incr)");
// #endif
					
#ifdef _DEBUG_KERNCLIQUEINCR			
					if(_DEBUG_KERNCLIQUEINCR) {
						std::cout << "   -> " << x; //<< std::endl;
						std::cout.flush();
					}
#endif				
					
					// if x has itself been removed, there's no point
					if(node.contain(x)) {
						
						
						if(neighboring_clique(x, w1, w2)) {
							
#ifdef _DEBUG_REDUCTION			
				if(_DEBUG_REDUCTION) {							
// #ifdef _DEBUG_KERNCLIQUEINCR
// 							if(_DEBUG_KERNCLIQUEINCR) {
								std::cout << "  remove " << x << " b/c its neighborhood N-" << co_neighbor[x] << " is a clique" << std::endl;
							}
#endif
							
							reason = (algo.use_backjump ? algo.explain_with(co_neighbor[x].end(), co_neighbor[x].begin()) : NO_REASON);
							rem_and_update(x);
							algo.exclude_from_cover(x, reason);	
							reduction = true;			
							
										
							//if(!size()) return true;
							//fix_point = (!size());//false;
							//reduction = true;
						} else {
							
#ifdef _DEBUG_KERNCLIQUEINCR		
							if(_DEBUG_KERNCLIQUEINCR) {
								std::cout << ":(" << w1 << "," << w2 << ")" << std::endl;
									//std::cout.flush();
							}
#endif
														

#ifdef _VERIF_WATCHERS	
							verify_watchers("before rem (incr)");
#endif							
							
							watchers.unwatch(q,y);
							
							
#ifdef _DEBUG_KERNCLIQUEINCR		
							if(_DEBUG_KERNCLIQUEINCR) {
								std::cout << " unwatch ok" << std::endl;
									//std::cout.flush();
							}
#endif						
							
#ifdef _VERIF_WATCHERS	
							verify_watchers("after rem (incr)");
#endif
							
							watchers.watch(x,w1,w2);
							
#ifdef _DEBUG_KERNCLIQUEINCR		
							if(_DEBUG_KERNCLIQUEINCR) {
								std::cout << " watch ok" << std::endl;
									//std::cout.flush();
							}
#endif		

							
#ifdef _VERIF_WATCHERS
							verify_watchers("after add (incr)");
#endif

						}
						
					}
#ifdef _DEBUG_KERNCLIQUEINCR	
					else if(_DEBUG_KERNCLIQUEINCR) {
						std::cout << " not in the graph" << std::endl ;
					}
#endif					
				}		
						
				
	
				
				++algo._last_checked;
			}
			
#ifdef _VERIF_WATCHERS
					verify_watchers("end incr");
#endif
			
			
			return reduction;
		}
		
		
		bool kernelize_crowns(VCAlgo<ReversibleCompactGraph>& algo)
		{
			find_maximal_matching();
			if (nmatch > algo._upper_bound-(int)(algo._cover.size)) return false;

			BipartiteGraph *H = new BipartiteGraph(capacity);
			int u,v;
			bool empty_crown = true;
			
			for (size_t i=0; i<node.size; i++) {
				u = node[i];
				if (!is_matched(u)) {
					empty_crown = false;
					H->add_node(u);
					H->set_right(u);
					for (size_t j=0; j<neighbor[u].size; j++) {
						v = neighbor[u][j];
						if (!H->node.contain(v)) {
							H->add_node(v);
							H->set_left(v);
						}
						H->add_undirected(u,v);
					}
				}
			}
			
			if (empty_crown) {
				delete H;
				return true;
			}
			
			H->maximum_matching();
			
			if (H->nmatch > algo._upper_bound-(int)(algo._cover.size)) {
				delete H;
				return false;
			}
			
			if (H->nmatch == (int)(H->left.size)) {
				for (size_t i=0; i<H->left.size; i++) {
					u = H->left[i];
					algo.add_to_cover(u, NO_REASON);
				}
				for (size_t i=0; i<H->right.size;i++) {
					u = H->right[i];
					algo.save(u, NO_REASON);
					rem_and_update(u);
					//remove_and_save(u, NO_REASON);
					//exclude_from_cover(u, NO_REASON);
				}
			} else {
				IntStack *I = new IntStack(0,H->capacity-1,H->capacity,false);
				for (size_t i=0; i<H->right.size; i++) {
					u = H->right[i];
					if (!H->is_matched(u)) I->add(u);
				}
				
				if (I->size == 0) {
					delete I;
					delete H;
					return true;
				}
				
				bool updated = true;
				while (updated) {
					updated=false;
					for (size_t i=0; i<I->size; i++) {
						u = (*I)[i];
						for (size_t j=0; j<H->neighbor[u].size;j++) {
							v = H->neighbor[u][j];
							if (H->is_matched(v) && (!I->contain(H->matching[v]))) {
								I->add(H->matching[v]);
								updated=true;
							}
						}
					}
				}
				int ns = 0;
				int nf = 0;
				for (size_t i=0; i<I->size;i++) {
					u = (*I)[i];
					algo.save(u, NO_REASON);
					rem_and_update(u);
					//remove_and_save(u, NO_REASON);
					//exclude_from_cover(u, NO_REASON);
					ns++;
					for (size_t j=0; j<H->neighbor[u].size; j++) {
						v = H->neighbor[u][j];
						if (node.contain(v)) {
							algo.add_to_cover(v, NO_REASON);
							nf++;
						}
					}
				}
				assert(ns > nf);
				delete I;
			}
	
			delete H;
			return true;
		}
		
		bool kernelize_crowns_multi(VCAlgo<ReversibleCompactGraph>& algo)
		{
			int prevsize = size()+1;
			while (prevsize != size()) {
				prevsize = size();
				for (int i=0; i<algo.nbtry; i++) {
					if (!kernelize_crowns(algo)) return false;
				}
			}
			return true;
		}
		
		bool kernelize_NT(VCAlgo<ReversibleCompactGraph>& algo)
		{
			int v;
			
		  	// Generate the bipartite graph
		  	BipartiteGraph *H = new BipartiteGraph(capacity*2);
		  	
		  	for (int i=0; i<size(); i++) {
		  		H->add_node(node[i]);
		  		H->add_node(node[i]+capacity);
		  		H->set_left(node[i]);
		  		H->set_right(node[i]+capacity);
		  	}
		  	
		  	for (int i=0; i<size(); i++) {
		  		for (size_t j=0; j<neighbor[node[i]].size; j++) {
		  			H->add_undirected(node[i], neighbor[node[i]][j]+capacity);
				}
			}
			
			H->vertex_cover();
			
			for (size_t i=0; i<H->left.size; i++) {
				v = H->left[i];
				if ((H->vcover.contain(v)) && (H->vcover.contain(v+capacity))) {
					//add_to_cover(H->node[i], NO_REASON);
					algo.add_to_cover(v, NO_REASON);
				} else if ((!H->vcover.contain(v)) && (!H->vcover.contain(v+capacity))) {
					algo.save(v, NO_REASON);
					rem_and_update(v);
					//remove_and_save(v, NO_REASON);
				}
			}
			
			delete H;
			
			// check the number of vertices
			if (size() > 2*(algo._upper_bound-(int)(algo._cover.size))) return false;
			
			return true;	
		  		
		}
		
		bool kernelize_NT_strong(VCAlgo<ReversibleCompactGraph>& algo)
		{
		
			int v;
			
		  	// Generate the bipartite graph
		  	BipartiteGraph *H = new BipartiteGraph(capacity*2);
		  	
		  	for (int i=0; i<size(); i++) {
		  		H->add_node(node[i]);
		  		H->add_node(node[i]+capacity);
		  		H->set_left(node[i]);
		  		H->set_right(node[i]+capacity);
		  	}
		  	
		  	for (int i=0; i<size(); i++) {
		  		for (size_t j=0; j<neighbor[node[i]].size; j++) {
		  			H->add_undirected(node[i], neighbor[node[i]][j]+capacity);
				}
			}
			
			H->persistently_matched_left_vertices();
		
			for (size_t i=0; i<H->left.size; i++) {
				v = H->left[i];
				if (!H->lpmatched.contain(v)) {
					while (neighbor[v].size != 0) {
						algo.add_to_cover(neighbor[v][0], NO_REASON);
					}
					algo.save(v, NO_REASON);
					rem_and_update(v);
				}
			}
			
			// check the number of vertices
			if (size() > 2*(algo._upper_bound-(int)(algo._cover.size))) return false;
			
			delete H;
			
			return true;
		}
		
		
		
		bool kernelize_cliques_incremental_old(VCAlgo<ReversibleCompactGraph>& algo) {
			int x, y, z, t, w1, w2, ixz, iut, reason, reduction=false; 
			
#ifdef _DEBUG_KERNCLIQUEINCR			
			if(_DEBUG_KERNCLIQUEINCR) {
				std::cout << "\nkernelize for the clique dominance (incremental, checking the removed witnesses:";
				for(i=algo._last_checked; i<algo._trail.size; ++i) {
					std::cout << " " << algo._trail[i];
				}
				std::cout << ")" << std::endl;
			}
#endif
			
			
			while(num_edges && algo._last_checked < (int)(algo._trail.size)) {
				
				// first apply the low degree rule
				while(num_edges && min_degree<=1) {
					x = node_of_degree[min_degree][0];
				// while(num_edges && node_of_degree[1].size>0) {
				// 	x = node_of_degree[1][0];
					
#ifdef _DEBUG_KERNCLIQUEINCR		
					if(_DEBUG_KERNCLIQUEINCR) {
						std::cout << "  remove " << x << " b/c low degree (" << degree(x) << ")" << std::endl;
					}
#endif
					
					reason = (algo.use_backjump ? algo.explain_with(co_neighbor[x].begin(), co_neighbor[x].end()) : NO_REASON);
					rem_and_update(x);
					algo.exclude_from_cover(x, reason);	
					reduction = true;		
				}
				if(!size()) break;
				
				
				// now take the first removed node y and check the nodes for which it was witness 
				y = algo._trail[algo._last_checked];
				
#ifdef _DEBUG_KERNCLIQUEINCR			
				if(_DEBUG_KERNCLIQUEINCR) {
					std::cout << "  NOT HERE check non cliques witnessed by " << y << ": " << watchers.watched_by[y] << std::endl;
				}
#endif
				
				while(num_edges && !watchers.watched_by[y].empty()) {
					x   = watchers.watched_by[y].pop();
					ixz = watchers.watched_index[y].pop();
					z   = watchers.watched_other[y].pop();
					
					// we also need to remove x from its other witness' list
					t = watchers.watched_other[z].pop();
					iut = watchers.watched_index[z].pop();
					//int u = watchers.watched_by[z].pop();
					
					watchers.watched_index[t][iut] = ixz;
					
					watchers.watched_index[z][ixz] = iut;
					watchers.watched_other[z][ixz] = t;
					watchers.watched_by[z][ixz] = watchers.watched_by[z].pop();
					
#ifdef _VERIF_WATCHERS	
					verify_watchers("after rem (incr)");
#endif
					
#ifdef _DEBUG_KERNCLIQUEINCR			
					if(_DEBUG_KERNCLIQUEINCR) {
						std::cout << "   -> " << x; //<< std::endl;
					}
#endif				
					
					// if x has itself been removed, there's no point
					if(node.contain(x)) {
						
						
						if(neighboring_clique(x, w1, w2)) {
							
#ifdef _DEBUG_KERNCLIQUEINCR		
							if(_DEBUG_KERNCLIQUEINCR) {
								std::cout << "  remove " << x << " b/c its neighborhood is a clique" << std::endl;
							}
#endif
							
							reason = (algo.use_backjump ? algo.explain_with(co_neighbor[x].begin(), co_neighbor[x].end()) : NO_REASON);
							rem_and_update(x);
							algo.exclude_from_cover(x, reason);	
							reduction = true;						
							//if(!size()) return true;
							//fix_point = (!size());//false;
							//reduction = true;
						} else {
							// store witnesses
							watchers.watched_index[w1].add(watchers.watched_by[w2].size);
							watchers.watched_index[w2].add(watchers.watched_by[w1].size);
							watchers.watched_other[w1].add(w2);
							watchers.watched_other[w2].add(w1);
							watchers.watched_by[w1].add(x);
							watchers.watched_by[w2].add(x);
							// move on
							
							
#ifdef _DEBUG_KERNCLIQUEINCR			
							if(_DEBUG_KERNCLIQUEINCR) {
								std::cout << "  watched" << std::endl;
							}
#endif
							
#ifdef _VERIF_WATCHERS
							verify_watchers("after add (incr)");
#endif

						}
						
					}
#ifdef _DEBUG_KERNCLIQUEINCR	
					else if(_DEBUG_KERNCLIQUEINCR) {
						std::cout << " not in the graph" << std::endl ;
					}
#endif					
				}		
						
	
				
				++algo._last_checked;
			}
			
			return reduction;
		}
		
		
#endif
				
		
		int my_fucking_clique_cover() {
			
			int size_before = node.size;
			Vector<int> clique;
			
			
			
			
			
			int num_cliques = 0;
			while(!node.empty()) {
				
				//display(std::cout);
				
				++num_cliques;
				clique.clear();
				
				int saved_size = node.size;
				while(!node.empty()) {
					// start a clique with the vertex of lowest degree, then on next steps select the node of highest degree
					int d = (clique.empty() ? min_degree : max_degree);
					int x = node_of_degree[d].back();
					for(int i=node_of_degree[d].size-2; i>=0; --i) {
						int y = node_of_degree[d][i];
						if(codegree(y)<codegree(x)) {
							x = y;
						}
					}
					
					// if(clique.empty()) {
					// 	std::cout << std::endl;
					// }
					//
					// std::cout << std::endl << "choose " << x << std::endl;
					
					// 1/ remove x from the graph
					rem_and_update(x);
					clique.add(x);
					// 2/ remove its non-neighbors 
					int num_candidates = 0;
					//    a/ put its neighbors at the front
					for(size_t i=0; i<neighbor[x].size; ++i) {
						int y = neighbor[x][i];
						if(node.contain(y)) {
							node.move(y, num_candidates++);
						}
					}
					//    b/ remove the other
					while((int)(node.size)>num_candidates) {
						rem_and_update(node.back());
					}
					
					
					//display(std::cout);
					
				}
				
				// at this point the graph is empty, undo the removals from the last clique, and remove only the vertices of the clique
				while((int)(node.size)<saved_size) {
					add_and_update(node.next());
				}
				while(!clique.empty()) {
					rem_and_update(clique.pop());
				}
				
			}
			
			while((int)(node.size)<size_before) {
				add_and_update(node.next());
			}
	
			return num_cliques;
			
		}
		
		
		
		int min_clique_cover_list() {
			
#ifdef _DEBUG_CLIQUECOV			
			std::cout << "compute CC " << std::endl;
#endif
			
			if(node.empty()) return 0;
			
			int i, x = node_of_degree[min_degree].back(), y;
			int start_clique;
			int start_candidates;
			int end_candidates;
			
			int saved_size = node.size;

#ifdef _DEBUG_CLIQUECOV			
			std::cout << "start with " << x << std::endl;
#endif
			
			for(i=0; i<(int)(node.size); ++i) { saturation[node[i]] = 0; }
			
			int num_cliques = 0;
			while(node.size) {
#ifdef _DEBUG_CLIQUECOV
				std::cout << "build a clique around " << x << std::endl;
#endif
				
				++num_cliques;
				neighbor_clique(x, node, start_clique, start_candidates, end_candidates, true);
				
#ifdef _DEBUG_CLIQUECOV
				std::cout << "clique:";
				for(size_t j=start_clique; j<node.size; ++j) {
					std::cout << " " << node[j];
				}
				std::cout << std::endl;
#endif
				
				node.size = start_clique;				
				
				x = node[0];
				
#ifdef _DEBUG_CLIQUECOV
				std::cout << "pick next node: " << x << "(" << saturation[x] << "/" << degree(x) << ")";
#endif

				for(i=1; i<(int)(node.size); ++i) { 
					y = node[i];
					
#ifdef _DEGREE_					
					if(degree(y) < degree(x)) { 
#elif _SAT_DEG_
					if(saturation[y] > saturation[x] || (saturation[y] == saturation[x] && degree(y) < degree(x))) { 
#elif _DEG_SAT_
					if(degree(y) < degree(x) || (degree(y) == degree(x) && saturation[y] > saturation[x])) { 
#elif _DEG_PLUS_SAT_						
					if(degree(y)-saturation[y] < degree(x)-saturation[x]) {
#endif
						
						x=y; 
					
#ifdef _DEBUG_CLIQUECOV
						std::cout << " " << x << "(" << saturation[x] << "/" << degree(x) << ")";
#endif
				
					} 
				}
				
#ifdef _DEBUG_CLIQUECOV				
				std::cout << std::endl;
#endif
				
			}

			node.size = saved_size;

			return num_cliques;
		}
		
		
		
		
		// returns true if the neighbors of x in V form a clique
		// if true, the clique is at the back of V starting at "start_clique"
		// if false, a maximal clique with a subset of neighbors of x is at the back of V and all vertices at the front
		// (until "end_candidates") are not neighbors of the node at index "start_clique"
		bool neighbor_clique(const int x, IntStack& V,
		int& start_clique,
		int& start_candidates,
		int& end_candidates,
		const bool maximize) { //, int& witness_a, int& witness_b) {

#ifdef _DEBUG_ADJCC
			if(_DEBUG_ADJCC) {
				if(maximize)
					std::cout << std::endl << "  find a maximal clique in " << V  << " around " << x << std::endl;
				else
					std::cout << std::endl << "  check if " << x << "'s neighborhood in " << V << " is a clique" << std::endl;
			}
#endif

			int a, b, c, d=capacity+1;
#ifndef _DEGREE_		
			int s=-d;
#endif
			bool full_clique=true; 


			start_clique = V.size;
			V.move(x, --start_clique);
			start_candidates = 0;
			end_candidates = V.size;
			for(int i=degree(x)-1; start_candidates<end_candidates && i>=0; --i) {

#ifdef _COUNT_OP
				++count;
#endif

				b = neighbor[x][i];
				--saturation[b];
				if(V.contain(b)) {
					V.move(b, start_candidates++);
#ifdef _DEGREE_						
					if(degree(b) < d) {
#elif _SAT_DEG_
					if(saturation[b] > s || (saturation[b] == s && degree(b) < d)) {
#elif _DEG_SAT_
					if(degree(b) < d || (degree(b) == d && saturation[b] > s)) {
#elif _DEG_PLUS_SAT_						
					if(degree(b)-saturation[b] < d-s) {
#endif
						c = b;
						d = degree(b);
#ifndef _DEGREE_		
						s = saturation[b];
#endif
					}
				}
			}

			//end_candidates = start_candidates;

#ifdef _DEBUG_ADJCC
			if(_DEBUG_ADJCC) {
				std::cout << "  (1) [";
				for(int i=0; i<start_candidates; ++i) {
					if(i) std::cout << " ";
					std::cout << V[i];
				}
				std::cout << "] - {...} - [" ;
				for(int i=start_clique; i<V.size; ++i) {
					if(i>start_clique) std::cout << " ";
					std::cout << V[i];
				}
				std::cout << "]" << std::endl;
			}
#endif
			
#ifdef _DEBUG_ADJCC
			if(_DEBUG_ADJCC) {
				std::cout << "here " << start_candidates << std::endl;
			}
#endif
			

			end_candidates = start_candidates;
			while(start_candidates>0 && (maximize || full_clique)) {
				end_candidates = start_candidates-1;

				// set a as the candidate of lowest degree
				a = c;
				// remove it from [candidates]
				V.move(a, --start_candidates);
				// move it to [clique]
				V.move(a, --start_clique);

				// at this point we have [candidates]-end_candidate-[not-clique]-start_clique-[clique]
				// with "a" in clique and V empty

#ifdef _DEBUG_ADJCC
				if(_DEBUG_ADJCC) {
					std::cout << "  move " << a << " to the clique and intersect with its neighbors = " << neighbor[a] << std::endl
						<< "  (2) [";
					for(int i=0; i<start_candidates; ++i) {
						if(i) std::cout << " ";
						std::cout << V[i];
					}
					std::cout << "] - {...} - [" ;
					for(int i=start_clique; i<V.size; ++i) {
						if(i>start_clique) std::cout << " ";
						std::cout << V[i];
					}
					std::cout << "]" << std::endl;
				}
#endif

				// we go through a's neighbors and add them into V (i.e. put them first in [candidates])
				start_candidates = 0;
				d=capacity+1;
				
#ifndef _DEGREE_		
				s=-d;
#endif
						
			
				
				for(int i=degree(a)-1; start_candidates<end_candidates && i>=0; --i) {

#ifdef _COUNT_OP
					++count;
#endif

					b = neighbor[a][i];
					--saturation[b];
					if(V.contain(b, 0, end_candidates)) {
						V.move(b,start_candidates++);
#ifdef _DEGREE_						
						if(degree(b) < d) {
#elif _SAT_DEG_
						if(saturation[b] > s || (saturation[b] == s && degree(b) < d)) {
#elif _DEG_SAT_
						if(degree(b) < d || (degree(b) == d && saturation[b] > s)) {
#elif _DEG_PLUS_SAT_						
						if(degree(b)-saturation[b] < d-s) {
#endif
							c = b;
							d = degree(b);

#ifndef _DEGREE_		
							s = saturation[b];
#endif
						}
					}
				}


				if(end_candidates > start_candidates) {

					full_clique = false;

				} 

#ifdef _DEBUG_ADJCC
				if(_DEBUG_ADJCC) {
					std::cout << "  (3) [";
					for(int i=0; i<start_candidates; ++i) {
						if(i) std::cout << " ";
						std::cout << V[i];
					}
					std::cout << "] - {...} - [" ;
					for(int i=start_clique; i<V.size; ++i) {
						if(i>start_clique) std::cout << " ";
						std::cout << V[i];
					}
					std::cout << "]" << std::endl;
				}
#endif

			}
			
			
			//std::cout << "return " << full_clique << std::endl;

			return full_clique;
		}
		
		
		void clique_dominance_incremental(VCAlgo<ReversibleCompactGraph>& algo) {
			int i, x, y, w1, w2; 
			int start_clique;
			int start_candidates;
			int end_candidates;
			
			while(algo._last_checked < (int)(algo._trail.size)) {
				y = algo._trail[algo._last_checked];
				
				//for(i=0; i<watched_by[y].size; ++i) {
				for(i=watchers.watched_by[y].size-1; i>=0; ++i) {
					x = watchers.watched_by[y][i];
						
					if(neighbor_clique(x, node, start_clique, start_candidates, end_candidates, false)) {
					
#ifdef _DEBUG_CLIQUEDOM			
						if(_DEBUG_CLIQUEDOM) {		
						std::cout << "-> clique:" ;
						for(int j=start_clique; j<node.size; ++j) {
							std::cout << " " << node[j];
						}
						std::cout << std::endl;
					}
#endif
					

					
						for(size_t j=start_clique; j<node.size; ++j) {						
							util_stack.remove(node[j]);
						}
										 
						
						rem_and_update(x);
						algo.exclude_from_cover(x, algo.explain_with(co_neighbor[x].begin(), co_neighbor[x].end()));			
					
					
					
					} else {
					
#ifdef _DEBUG_CLIQUEDOM			
						if(_DEBUG_CLIQUEDOM) {				
						std::cout << "-> not a full clique: {" ;
						for(int j=start_clique; j<node.size; ++j) {
							std::cout << " " << node[j];
						}
						std::cout << " } / {";
						for(int j=0; j<start_candidates; ++j) {
							std::cout << " " << node[j];
						}	
						std::cout << " |";
						for(int j=start_candidates; j<end_candidates; ++j) {
							std::cout << " " << node[j];
						}							
						std::cout << " }" << std::endl;
					}
#endif				

						// x's neighborhood is not a clique, and node[start_clique],node[end_candidates] is a witness
						w1 = node[start_clique];
						w2 = node[start_candidates];
					
						//std::cout << "witnesses = " << w1 << " " << w2 << std::endl;
					
						watchers.watched_by[w1].add(x);
						watchers.watched_by[w2].add(x);
						for(size_t j=start_clique+2; j<node.size; ++j) {
							if(util_stack.contain(node[j])) {
								util_stack.remove(node[j]);
								watchers.watched_by[w1].add(node[j]);
								watchers.watched_by[w2].add(node[j]);
							}
						}
					}
				}
				
				++algo._last_checked;
			}
		}
		
		bool clique_dominance(VCAlgo<ReversibleCompactGraph>& algo) {
			int i, x, w1, w2; 
			int start_clique;
			int start_candidates;
			int end_candidates;
			bool reduction = false;

			if(num_edges) {
				util_stack.clear();
				for(i=node.size-1; i>=0; --i) {
					util_stack.add(node[i]);
				}
			
				while(!util_stack.empty()) {
					
					x = util_stack.pop();
				
	#ifdef _DEBUG_CLIQUEDOM
					if(_DEBUG_CLIQUEDOM) {		
					std::cout << std::endl << util_stack ; //<< std::endl << "check " << x ;
				}
	#endif
					
					
				
					if(neighbor_clique(x, node, start_clique, start_candidates, end_candidates, false)) {
					
	#ifdef _DEBUG_CLIQUEDOM					
						if(_DEBUG_CLIQUEDOM) {		
							
								//std::cout << "start clique = " << start_clique << " size " << node.size << std::endl;
							
						std::cout << "-> clique:" ;
						for(int j=start_clique; j<node.size; ++j) {
							std::cout << " " << node[j];
						}
						std::cout << std::endl;
					}
	#endif
					
						reduction = true;
						for(size_t j=start_clique; j<node.size; ++j) {
							
							if(node[j] != x) {
#ifdef _DEBUG_CLIQUEDOM			
						if(_DEBUG_CLIQUEDOM) {		
							std::cout << "remove " << node[j] << " from " << util_stack << std::endl;
						}
	#endif
							
							util_stack.remove(node[j]);
						}
						}
						
						
#ifdef _DEBUG_CLIQUEDOM			
						if(_DEBUG_CLIQUEDOM) {		
							std::cout << "remove " << x << " from cover" << std::endl;
						}
	#endif		
						
						rem_and_update(x);
						algo.exclude_from_cover(x, algo.explain_with(co_neighbor[x].begin(), co_neighbor[x].end()));			
						
						
#ifdef _DEBUG_CLIQUEDOM			
						if(_DEBUG_CLIQUEDOM) {		
							std::cout << "ok" << std::endl;
						}
	#endif		
					
					} else {
					
	#ifdef _DEBUG_CLIQUEDOM			
						if(_DEBUG_CLIQUEDOM) {				
						std::cout << "-> not a full clique: {" ;
						for(int j=start_clique; j<node.size; ++j) {
							std::cout << " " << node[j];
						}
						std::cout << " } / {";
						for(int j=0; j<start_candidates; ++j) {
							std::cout << " " << node[j];
						}	
						std::cout << " |";
						for(int j=start_candidates; j<end_candidates; ++j) {
							std::cout << " " << node[j];
						}							
						std::cout << " }" << std::endl;
					}
	#endif				

						// x's neighborhood is not a clique, and node[start_clique],node[end_candidates] is a witness
						w1 = node[start_clique];
						w2 = node[start_candidates];
					
						//std::cout << "witnesses = " << w1 << " " << w2 << std::endl;
					
						watchers.watched_by[w1].add(x);
						watchers.watched_by[w2].add(x);
						for(size_t j=start_clique+2; j<node.size; ++j) {
							if(util_stack.contain(node[j])) {
								util_stack.remove(node[j]);
								watchers.watched_by[w1].add(node[j]);
								watchers.watched_by[w2].add(node[j]);
							}
						}
					}
				}
			}
			
			
			return reduction;
		}
		
		
		inline int codegree(const int x) { return co_neighbor[x].size; }
		
		
		void verify(const char* msg)
		{
			
			for(int i=0; i<capacity; ++i) {
				assert(original_size[i]==degree(i)+codegree(i));
			}
			
			// check that edges are in both directions
			for(size_t i=0; i<node.size; ++i) {
				int x = node[i];
				for(size_t j=0; j<neighbor[x].size; ++j) {
					int y = neighbor[x][j];			
					if(neighbor[y][nb_index[x][j]] != x) {
						// std::cout << std::endl;
						// display(std::cout);
						std::cout << msg << ": inconsistency for edge " << x << "," << y << std::endl;
						exit(1);
					}
					//assert(neighbor[y][nb_index[x][j]] == x);
				}
			}
			
			// check min/max degree
			int lb=capacity;
			int ub=0;
			int ne=0;
			for(size_t i=0; i<node.size; ++i) {
				int x = node[i];
				int d = degree(x);
				if(d < lb) lb = d;
				if(d > ub) ub = d;
				
				ne += d;
				
				//assert(node_of_degree[d].contain(x));
				if(node_of_degree[d][degree_index[x]] != x) {
					std::cout << *this << std::endl << msg << ": node_of_degree[" << d << "] does not contain " 
						<< x << " (" << node_of_degree[d] << ")" << std::endl;
					exit(1);
				}
			}
			
			if(node.size>0) {
			if(lb!=min_degree) {
				std::cout << *this << std::endl << msg << ": min_degree=" << min_degree << " (should be " << lb << ")" << std::endl;
				exit(1);
			}
			if(ub!=max_degree) {
				std::cout << *this << std::endl << msg << ": max_degree=" << max_degree << " (should be " << ub << ")" << std::endl;
				exit(1);
			}
			}
			ne /= 2;
			
			
			//std::cout << "NUM EDGES: " << num_edges << " == " << ne << std::endl;
			
			if(ne!=num_edges) {
				std::cout << *this << std::endl << msg << ": num_edges=" << num_edges << " (should be " << ne << ")" << std::endl;
				exit(1);
			}
			
			// check degree ordering
			for(int d=min_degree; d<=max_degree; ++d) {
				for(size_t i=0; i<node_of_degree[d].size; ++i) {
					if(degree(node_of_degree[d][i]) != d) {
						std::cout << *this << std::endl << msg << ": degree of " << node_of_degree[d][i] << " is not " << d << std::endl;
						exit(1);
					}
				}
			}
					
			
		}
		
		
		//#define _DEBUG_UPDATE

		void add_and_update(const int x) {
			
#ifdef _VERIFIED_RCG
			verify("before add");
#endif			
			

#ifdef _DEBUG_UPDATE
			assert(x>=0);
			assert(x<capacity);
			assert(!node.contain(x));
			std::cout << std::endl << *this << "\nadd node " << x << " " ;
			neighbor[x].display(std::cout);
			std::cout << std::endl;
#endif			
			
			node.add(x);
			node_set.add(x);
			int i = degree(x), y, d, l, idx;
			degree_index[x] = node_of_degree[i].size;
			node_of_degree[i].add(x);
			num_edges += i;
			
			while(i--) {
				y = neighbor[x][i];
				d = degree(y);
				
				//std::cout << " add in " << y << "'s neighbor list" << std::endl; 
				
				idx = degree_index[y], l=node_of_degree[d].pop();
				node_of_degree[d][idx] = l;
				degree_index[l] = idx;
				
				degree_index[y] = node_of_degree[d+1].size;				
				node_of_degree[d+1].add(y);
				
				
				if(max_degree==d)
					max_degree = d+1;
				if(min_degree==d && node_of_degree[d].empty()) {
					min_degree = d+1;
					//std::cout << "increase min degree (" << d << " => " << (d+1) <<") b/c " << "node_of_degree[" << d << "] = " << node_of_degree[d] << " is empty" << std::endl; 
				}
				
				// the new rank of x in y's neighbors
				nb_index[x][i] = neighbor[y].size;
				
				// the rank of y in x's neighbors
				nb_index[y].add(i);
				
				// add x in y's neighbors
				neighbor[y].add(x);
				
				// pop x from y's co-neighbor
				assert(co_neighbor[y].back()==x);
				co_neighbor[y].pop();
			}
			
			d = degree(x);
			if(max_degree < d) max_degree = d;
			if(min_degree > d) min_degree = d;
	
#ifdef _VERIFIED_RCG_
			verify("after add");
#endif
				
		}

		
		void rem_and_update(const int x) {
			
#ifdef _VERIFIED_RCG
			verify("before rem");
#endif	

			
#ifdef _DEBUG_UPDATE			
			assert(x>=0);
			assert(x<capacity);
			assert(node.contain(x));
			std::cout << std::endl << *this << "\nremove node " << x << " " ;
			neighbor[x].display(std::cout);
			std::cout << std::endl;
#endif			

			node.remove(x);
			node_set.remove(x);
			int i = degree(x), y, d, l=node_of_degree[i].pop(), idx=degree_index[x];
			node_of_degree[i][idx] = l;
			degree_index[l] = idx;
			num_edges -= i;
			
			while(i--) {
				y = neighbor[x][i];
				d = degree(y);
			
				idx = degree_index[y], l=node_of_degree[d].pop();
				node_of_degree[d][idx] = l;
				degree_index[l] = idx;
				
				degree_index[y] = node_of_degree[d-1].size;				
				node_of_degree[d-1].add(y);

				if(min_degree==d)
					min_degree = d-1;
				
				if(max_degree==d && node_of_degree[d].empty())
					max_degree = d-1;
			
				// x is at rank idx in neighbor[y]
				idx = nb_index[x][i];
				// replace x by z
				int z = neighbor[y].pop();
				neighbor[y][idx] = z;
				// y is at rank idy in neighbor[z]
				int idy = nb_index[y].pop();
				// change the rank of z in neighbor[y]
				nb_index[z][idy] = idx;
				// copy the rank of y in neighbor[z] to its new location
				nb_index[y][idx] = idy;				
							
				// add x to y's co-neighbor
				co_neighbor[y].add(x);
			}
			
			d = degree(x);
			if(max_degree==d) {
				while(d>=min_degree && node_of_degree[d].empty()) --d;
				max_degree = d;
			}
			if(min_degree==d) {
				while(d<=max_degree && node_of_degree[d].empty()) ++d;
				min_degree = d;
			}

#ifdef _DEBUG_UPDATE					
			verify("rem");
#endif
			
#ifdef _VERIFIED_RCG
			verify("after rem");
#endif	
					
		}
		
		

		
		
// 		// returns true if the neighbors of x in V form a clique
// 		// if true, the clique is at the back of V starting at "start_clique"
// 		// if false, a maximal clique with a subset of neighbors of x is at the back of V and all vertices at the front
// 		// (until "end_candidates") are not neighbors of the node at index "start_clique"
// 		bool neighbor_clique(const int x, IntStack& V,
// 		int& start_clique,
// 		int& start_candidates,
// 		int& end_candidates,
// 		const bool maximize) {
//
// #ifdef _DEBUG_ADJCC
// 			std::cout << std::endl << "  check if " << x << "'s neighborhood in " << V << " is a clique" << std::endl;
// #endif
//
// 			int a=x, b, c, d=capacity+1, full_clique=true;
// 			start_clique = V.size;
//
// 			V.move(x, --start_clique);
// 			start_candidates = 0;
// 			end_candidates = start_clique;
// 			//end_candidates = -1;
//
// 			do {
// 				for(int i=degree(a)-1; i>=0; --i) {
// #ifdef _COUNT_OP
// 					++count;
// #endif
// 					b = neighbor[a][i];
// 					if(V.contain(b, 0, end_candidates)) {
// 						V.move(b, start_candidates++);
// 						if(degree(b) < d) {
// 							c = b;
// 							d = degree(b);
// 						}
// 					}
// 				}
//
// 				if(a!=x && end_candidates > start_candidates) {
//
// 					full_clique = false;
//
// 				} else {
//
// 					end_candidates = start_candidates-1;
//
// 					// set a as the candidate of lowest degree
// 					a = c;
// 					// remove it from [candidates]
// 					V.move(a, --start_candidates);
// 					// move it to [clique]
// 					V.move(a, --start_clique);
//
// 					// at this point we have [candidates]-end_candidate-[not-clique]-start_clique-[clique]
// 					// with "a" in clique and V empty
//
// #ifdef _DEBUG_ADJCC
// 					std::cout << "move " << a << " to the clique and intersect with its neighbors = " << neighbor[a] << std::endl
// 						<< "  [";
// 					for(int i=0; i<start_candidates; ++i) {
// 						if(i) std::cout << " ";
// 						std::cout << V[i];
// 					}
// 					std::cout << "] - {...} - [" ;
// 					for(int i=start_clique; i<capacity; ++i) {
// 						if(i>start_clique) std::cout << " ";
// 						std::cout << V[i];
// 					}
// 					std::cout << "]" << std::endl;
// #endif
//
// 					// we go through a's neighbors and add them into V (i.e. put them first in [candidates])
// 					start_candidates = 0;
// 					d=capacity+1;
// 				}
// 			} while(start_candidates>0 && (maximize || full_clique));
//
//
// 			return full_clique;
// 		}
		
	
		
// 		void clique_dominance(VCAlgo<ReversibleCompactGraph>& algo) {
// 			int i, x; //, size_before = node.size;
// 			int start_clique;
// 			int start_candidates;
// 			int end_candidates;
//
// 			bool maximize = false;
//
// 			util_stack.clear();
// 			for(i=node.size-1; i>=0; --i) {
// 				util_stack.add(node[i]);
// 			}
//
// 			//for(i=util_stack.size-1; i>=0; --i) {
// 			while(!util_stack.empty()) {
// 				x = util_stack.pop();
//
// #ifdef _DEBUG_CLIQUEDOM
// 				std::cout << std::endl << util_stack ; //<< std::endl << "check " << x ;
// #endif
//
// 				//size_before = node.size;
//
// 				if(neighbor_clique(x, node, start_clique, start_candidates, end_candidates, maximize)) {
//
// #ifdef _DEBUG_CLIQUEDOM
// 					std::cout << "-> clique:" ;
// 					for(int j=start_clique; j<node.size; ++j) {
// 						std::cout << " " << node[j];
// 					}
// 					std::cout << std::endl;
// #endif
//
// 					for(int j=start_clique; j<node.size; ++j) {
// 						util_stack.remove(node[j]);
// 						//rem_and_update(node[j]);
// 					}
// 					rem_and_update(x);
// 					exclude_from_cover(x);
//
// 				} else {
//
// #ifdef _DEBUG_CLIQUEDOM
// 					std::cout << "-> not a full clique: {" ;
// 					for(int j=start_clique; j<node.size; ++j) {
// 						std::cout << " " << node[j];
// 					}
// 					std::cout << " } / {";
// 					for(int j=0; j<start_candidates; ++j) {
// 						std::cout << " " << node[j];
// 					}
// 					std::cout << " |";
// 					for(int j=start_candidates; j<end_candidates; ++j) {
// 						std::cout << " " << node[j];
// 					}
// 					std::cout << " }" << std::endl;
// #endif
//
// 					//if(maximize) {
// 						for(int j=start_clique; j<node.size-2; ++j) {
// 							if(util_stack.contain(node[j]))
// 								util_stack.remove(node[j]);
// 						}
// 						//}
//
// 				}
//
//
// 				//node.size = size_before;
// 			}
//
//
//
// //
// // 			util_stack.clear();
// //
// //
// // 			for(i=node.size-1; i>=0; --i) {
// // 				x = node[i];
// //
// // #ifdef _DEBUG_CLIQUEDOM
// // 				std::cout << " check " << x << std::endl;
// // #endif
// //
// // 				if(!util_stack.contain(x)) {
// // 					if(neighbor_clique(x, node, start_clique, end_candidates, false)) {
// // 						// the neighbors of x form a clique, we can remove it
// //
// // 						// first restore the set of nodes
// // 						node.size = size_before;
// //
// // 						rem_and_update(x);
// // 						algo.exclude_from_cover(x);
// //
// // 						// re-save the node size
// // 						size_before = node.size;
// // 					}
// // 				}
// //
// // #ifdef _DEBUG_CLIQUEDOM
// // 				else std::cout << "";
// // #endif
// //
// //
// // 			}
//
//
// #ifdef _COUNT_OP
// 			std::cout << "Num operations = " << count << std::endl;
// #endif
// 		}
//
//
		
	
		std::ostream& display(std::ostream& os) const {
			CompactGraph::display(os);
			
			std::cout << min_degree << ".." << max_degree << std::endl;
			
			int d = max_degree;
			
			while(d >= min_degree) {
				if(!node_of_degree[d].empty()) {
					os << d << ": [" << node_of_degree[d][0];
					for(size_t i=1; i<node_of_degree[d].size; ++i) {
						os << " " << node_of_degree[d][i];
					}
					os << "] ";
				}
				--d;
			}
			return os;
		}
	
	};
	
	
	std::ostream& operator<< (std::ostream& os, const ReversibleCompactGraph& x);

	std::ostream& operator<< (std::ostream& os, const ReversibleCompactGraph* x);
  






	class VertexCoverAlgo {


	public:

		int init_level;

		ReversibleGraph             _graph;
		int                   _lower_bound;
		int                   _upper_bound;
		Vector<int>                 _trail;
		Vector<int>              _decision;
		
#ifdef _COVER_SET
		IntStack                    _cover;
#else		
		Vector<int>                 _cover;
#endif
		
		// coloring bound [not used]
		int* _coloring;
		IntList  color_list;
		
		// greedy bound [not used]
		int * degree_count;
		
		// clique cover bound 
		int num_cliques;
#ifdef _USE_CLIQUES
		BitSet   * clique_set;
	  IntStack * clique_lst;
#endif
		BitSet   * candidates;
	
		// explanation stuff
		int*                    _rank;
		int*                   _level;
		IntList          _explanation;
		
#ifdef _STORE_NOGOOD		
		Vector< Array<int>* > _nogood;
		Array<int>**          _reason;
#else
		Vector< int >    _nogood;
		int*             _reason;
#endif 
		
		// used for dominance detection
		BitSet util_set;
		
		
		int suggested_node;

	public:		
		
		// stats
		double time_start;
		int nb_nodes;
		int nb_restarts;
		int nb_solutions;
		int nb_conflicts;
		long int nogood_size;
		long int backjump_size;
		bool optimal;
		
#ifdef _NOGOOD_STAT
		int* nogood_sizes;
		int* backjump_sizes;
#endif
				
		//params
		bool verbose;
		int nbtry;
		bool use_backjump;
		int fail_limit;
				
				
		int get_upper_bound() const { return _upper_bound; }

		VertexCoverAlgo(const Graph& g, const int k)
			: _graph(g), _upper_bound(k)
		{
			
			
			
			
			int n = _graph.capacity;
			
			init_level = 0;
			
			nb_nodes = 0;
			nb_restarts = -1;
			nb_solutions = 0;
			nb_conflicts = 0;
			nogood_size = 0;
			backjump_size = 0;
			optimal = false;
			
#ifdef _NOGOOD_STAT
			nogood_sizes   = new int[n];
			backjump_sizes = new int[n];
			std::fill(nogood_sizes,   nogood_sizes+n,   0);
			std::fill(backjump_sizes, backjump_sizes+n, 0);
#endif
			
			suggested_node = 0;
			//use_backjump = false;
			verbose = true;
			nbtry = 10;
			fail_limit = -1;
			
#ifdef _COVER_SET			
			_cover.initialise(0, n-1, n, false);
#endif
			
			degree_count = new int[n];
			

			num_cliques = 0;
#ifdef _USE_CLIQUES
			clique_set= new   BitSet[n];
			clique_lst= new IntStack[n];
#endif
			candidates= new   BitSet[n];
			for(int i=0; i<n; ++i) {
#ifdef _USE_CLIQUES
				clique_set [i].initialise(0, n-1, BitSet::empt);
				clique_lst [i].initialise(0, n-1, n, false);
#endif
				candidates [i].initialise(0, n-1, BitSet::full);
			}
			
			
			color_list.initialise(n);			
			_coloring = new int[n];
			std::fill(_coloring, _coloring+n, -1);
			
			
			_rank = new int[n+1];
			_rank[n] = -1;
			_level = new int[n+1];
			_level[n] = -1;
			_explanation.initialise(n);

#ifdef _STORE_NOGOOD	
			_reason = new Array<int>*[n];
#else
			_reason = new int[n];
#endif
			
			
			util_set.initialise(0, n-1, BitSet::empt);
		}

		virtual ~VertexCoverAlgo()
		{
			delete [] degree_count;
			
			delete [] candidates;
#ifdef _USE_CLIQUES
			delete [] clique_set;
			delete [] clique_lst;
#endif			
			
			delete [] _coloring;
			
			delete [] _rank;
			delete [] _level;
			delete [] _reason;
			
#ifdef _NOGOOD_STAT
			delete [] nogood_sizes;
			delete [] backjump_sizes;
#endif
		
#ifdef _STORE_NOGOOD	
			while(!_nogood.empty()) {
				Array<int>* cl = _nogood.pop();
				free(cl);
			}
#endif
			
		}
		
		void verify_cover() {
			util_set.clear();
			for(size_t i=0; i<_cover.size; ++i)
				util_set.add(_cover[i]);
			
			int n = _graph.capacity;
			for(int x=0; x<n; ++x) {
				for(int y=x+1; y<n; ++y) {
					if(_graph.matrix[x].contain(y)) {
						if(!util_set.contain(x) and !util_set.contain(y)) {
							
							std::cout << "residual graph:\n" << _graph << std::endl;
							std::cout << "The edge (" << x << ","<< y << ") is not covered by " << util_set << std::endl;
							exit(1);
						}
					}
				}
			}	
		}
		
		
		int explain_with(int* beg, int* end) {
			int reason = _nogood.size;
			int sz = (int)(end-beg); 
			_nogood.add(sz);
			for(int *cur=beg; cur!=end; ++cur) {
				_nogood.add(*cur);
			}
			return reason;
		}
		
		
		void verify_trail(const char* msg) {
		
			int x;
			int *reason;
			int rsize;
			int dec_idx = init_level;
			for(size_t i=0; i<_trail.size; ++i) {
				x = _trail[i];
				
				assert(_rank[x]==(int)i);
				if(i) assert(_level[x]>=_level[_trail[i-1]]);
				
				if(dec_idx<(int)(_decision.size) && x==_decision[dec_idx]) {
					assert(_reason[x]==NO_REASON);
					if(i) assert(_level[x]==_level[_trail[i-1]]+1);
					else assert(_level[x]<=init_level+1);
					++dec_idx; 
				} else {
					
					if(_reason[x]==NO_REASON) {
						std::cout << "ERROR: reason/trail of " << x << " (" << msg << ")"<< std::endl;
						print_trail(true);
					}
						
					assert(_reason[x]!=NO_REASON);
					if(_reason[x]!=DEDUCTION) {
						
						//std::cout << " (" << msg << ") reason[" << x << "] = " << _reason[x] << std::endl;
						
						if(_reason[x]<0) {
							std::cout << " (" << msg << ") reason[" << x << "] = " << _reason[x] << std::endl;
						}
						
						assert(_reason[x]>=0);
						
						if(_reason[x]>=(int)(_nogood.size)) {
							std::cout << " (" << msg << ") reason[" << x << "] = " << _reason[x] << " nogoods = " << _nogood << std::endl;
						}
						
						assert(_reason[x]<(int)(_nogood.size));
						
						if(_reason[x]+_nogood[_reason[x]]>=(int)(_nogood.size)) {
							std::cout << " (" << msg << ") reason[" << x << "] = " << _reason[x] << ", |reason[" 
								<< x << "]| = " << _nogood[_reason[x]] << " nogoods = " << _nogood << std::endl;
						}
						
						
						assert(_reason[x]+_nogood[_reason[x]]<(int)(_nogood.size));
					}
				}
				
				if(_reason[x]>=0) {
					reason = &(_nogood[_reason[x]+1]);
					rsize = _nogood[_reason[x]];
					for(int j=1; j<rsize; ++j) {
						if(_rank[reason[j-1]] <= _rank[reason[j]]) {
							std::cout << "ERROR: unsorted nogood (" << msg << ")"<< std::endl;
							std::cout << "literal " << x << " clause " << _reason[x] << ":";
							for(int k=0; k<rsize; ++k) {
								std::cout << " " << reason[k];	
							}
							std::cout << std::endl;
						}
						
						assert(_rank[reason[j-1]] > _rank[reason[j]]);
					}	
				}
			}
			
			for(size_t i=0; i<_graph.node.size; ++i) {
				int x = _graph.node[i];
				assert(!_cover.contain(x));
				for(size_t j=0; j<_graph.neighbor[x].size; ++j) {
					int y = _graph.neighbor[x][j];
					
					if(_cover.contain(y)) {
						std::cout << msg << " inconsistency graph/cover: " << y << " cannot be both in the graph and in the cover " << _cover 
							<< "\n" << _graph << std::endl;
						
					}
					
					assert(!_cover.contain(y));
				}
				for(int j=_graph.neighbor[x].size; j<_graph.original_size[x]; ++j) {
					int y = _graph.neighbor[x][j];
					
					if(!_cover.contain(y)) {
						std::cout << msg << " inconsistency: since " << x << " is in the graph, " << y << " must be either in the graph or in the cover " << _cover 
							<< "\n" << _graph << std::endl;
						
					}
					
					assert(_cover.contain(y));
				}
			}
			
		}
		
		
		bool kernelize_buss()
		{
			int high;
			while (_graph.max_degree >= (int)(_upper_bound-_cover.size)) {
				high = _graph.node_of_degree[_graph.max_degree].back();
				add_to_cover(high, NO_REASON);
			}
			
			// check the number of edges
			if (_graph.num_edges > _graph.max_degree*(_upper_bound-(int)(_cover.size))) {
				return false;
			}
			
			return true;
		}
		
		
		bool kernelize_cliques()
		{
			int d, i, j, x, y, reason, ds, k;
			bool fix_point = false;
			bool reduction = false;
			while(_graph.num_edges && !fix_point) {
				fix_point = true;
				//std::cout << _graph.min_degree << " <> " << _graph.max_degree << std::endl;
				d = _graph.min_degree;
				while(_graph.num_edges && d<=_graph.max_degree) {
					while(_graph.num_edges && _graph.min_degree<=1) {
						//std::cout << "there are nodes of degree <= 1 (" << _graph.min_degree << ")" << std::endl;
						
						//std::cout << _graph.node_of_degree[_graph.min_degree].size << std::endl;
						
						x = _graph.node_of_degree[_graph.min_degree][0];
						
						//std::cout << "remove " << x << " and add " << _graph.neighbor[x] << std::endl;
						//remove_and_save(x, NO_REASON);
						//save(x, NO_REASON);
						
						
						
						reason = NO_REASON;
						if(use_backjump) {
								
#ifdef _DEBUG_BACKJUMPS
							if(_DEBUG_BACKJUMPS) {	
								std::cout << "build nogood for clique kernel:";
							}
#endif								
								
							reason = _nogood.size;
							ds = _graph.original_size[x];
							_nogood.add(ds-_graph.degree(x));
							for(k=_graph.degree(x); k<ds; ++k) {
									
#ifdef _DEBUG_BACKJUMPS
								if(_DEBUG_BACKJUMPS) {	
									std::cout << " " << _graph.neighbor[x][k] ;
								}
#endif										
									
								_nogood.add(_graph.neighbor[x][k]);
							}
								
#ifdef _DEBUG_BACKJUMPS		
							if(_DEBUG_BACKJUMPS) {						
								std::cout << std::endl << "   -> explains " << x << " and its neighbors " << _graph.neighbor[x] << std::endl;
							}
#endif
								
						}
							
	
						_graph.rem_and_update(x);
						exclude_from_cover(x, reason);			
						if(!_graph.size()) return true;
						fix_point = false;
						reduction = true;
					}
					for(j=0; j<(int)(_graph.node_of_degree[d].size);) {
						x = _graph.node_of_degree[d][j];
						util_set.copy(_graph.matrix[x]);
						util_set.intersect_with(_graph.node_set);
						//std::cout << " examine " << x << " current neighbors = " << util_set << " " << _graph.neighbor[x] << std::endl;
						for(i=_graph.neighbor[x].size; --i>=0;) {
							y = _graph.neighbor[x][i];
							util_set.remove(y);
							if(!util_set.included(_graph.matrix[y])) {
								//std::cout << " not a clique b/c of " <<  y << " neighbors: " << (_graph.matrix[y]) << std::endl;
								break;
							}
							util_set.add(y);
						}
						if(i<0) {
							//std::cout << " -> it's a clique! remove " << x << " and add " << _graph.neighbor[x] << std::endl;
							// N(x) is a clique
							//remove_and_save(x, NO_REASON);
							//save(x, NO_REASON);
							
							
							
							reason = NO_REASON;
							if(use_backjump) {
								
	#ifdef _DEBUG_BACKJUMPS
								if(_DEBUG_BACKJUMPS) {	
									std::cout << "build nogood for clique kernel:";
								}
	#endif								
								
								reason = _nogood.size;
								ds = _graph.original_size[x];
								_nogood.add(ds-_graph.degree(x));
								for(k=_graph.degree(x); k<ds; ++k) {
									
	#ifdef _DEBUG_BACKJUMPS
									if(_DEBUG_BACKJUMPS) {	
										std::cout << " " << _graph.neighbor[x][k] ;
									}
	#endif										
									
									_nogood.add(_graph.neighbor[x][k]);
								}
								
	#ifdef _DEBUG_BACKJUMPS		
								if(_DEBUG_BACKJUMPS) {						
									std::cout << std::endl << "   -> explains " << x << " and its neighbors " << _graph.neighbor[x] << std::endl;
								}
	#endif
								
							}
							
							
							_graph.rem_and_update(x);
							exclude_from_cover(x, reason);							
							if(!_graph.size()) return true;
							fix_point = false;
							reduction = true;
						} else {
							// move on
							++j;
						}
					}
					if(_graph.min_degree<=1) {
						d = _graph.min_degree;
					} else {
						++d;
					}
				}
			}
			//std::cout << _graph << std::endl << "n=" << _graph.size() << std::endl;
			return reduction;
		}
		
		
		bool kernelize_crowns()
		{
			_graph.find_maximal_matching();
			if (_graph.nmatch > (int)(_upper_bound-_cover.size)) return false;

			BipartiteGraph *H = new BipartiteGraph(_graph.capacity);
			int u,v;
			bool empty_crown = true;
			
			for (size_t i=0; i<_graph.node.size; i++) {
				u = _graph.node[i];
				if (!_graph.is_matched(u)) {
					empty_crown = false;
					H->add_node(u);
					H->set_right(u);
					for (size_t j=0; j<_graph.neighbor[u].size; j++) {
						v = _graph.neighbor[u][j];
						if (!H->node.contain(v)) {
							H->add_node(v);
							H->set_left(v);
						}
						H->add_undirected(u,v);
					}
				}
			}
			
			if (empty_crown) {
				delete H;
				return true;
			}
			
			H->maximum_matching();
			
			if (H->nmatch > (int)(_upper_bound-_cover.size)) {
				delete H;
				return false;
			}
			
			if (H->nmatch == (int)(H->left.size)) {
				for (size_t i=0; i<H->left.size; i++) {
					u = H->left[i];
					add_to_cover(u, NO_REASON);
				}
				for (size_t i=0; i<H->right.size;i++) {
					u = H->right[i];
					save(u, NO_REASON);
					_graph.rem_and_update(u);
					//remove_and_save(u, NO_REASON);
					//exclude_from_cover(u, NO_REASON);
				}
			} else {
				IntStack *I = new IntStack(0,H->capacity-1,H->capacity,false);
				for (size_t i=0; i<H->right.size; i++) {
					u = H->right[i];
					if (!H->is_matched(u)) I->add(u);
				}
				
				if (I->size == 0) {
					delete I;
					delete H;
					return true;
				}
				
				bool updated = true;
				while (updated) {
					updated=false;
					for (size_t i=0; i<I->size; i++) {
						u = (*I)[i];
						for (size_t j=0; j<H->neighbor[u].size;j++) {
							v = H->neighbor[u][j];
							if (H->is_matched(v) && (!I->contain(H->matching[v]))) {
								I->add(H->matching[v]);
								updated=true;
							}
						}
					}
				}
				int ns = 0;
				int nf = 0;
				for (size_t i=0; i<I->size;i++) {
					u = (*I)[i];
					save(u, NO_REASON);
					_graph.rem_and_update(u);
					//remove_and_save(u, NO_REASON);
					//exclude_from_cover(u, NO_REASON);
					ns++;
					for (size_t j=0; j<H->neighbor[u].size; j++) {
						v = H->neighbor[u][j];
						if (_graph.node.contain(v)) {
							add_to_cover(v, NO_REASON);
							nf++;
						}
					}
				}
				assert(ns > nf);
				delete I;
			}
	
			delete H;
			return true;
		}
		
		bool kernelize_crowns_multi()
		{
			int prevsize = _graph.size()+1;
			while (prevsize != _graph.size()) {
				prevsize = _graph.size();
				for (int i=0; i<nbtry; i++) {
					if (!kernelize_crowns()) return false;
				}
			}
			return true;
		}
		
		bool kernelize_NT()
		{
			int v;
			
		  	// Generate the bipartite graph
		  	BipartiteGraph *H = new BipartiteGraph(_graph.capacity*2);
		  	
		  	for (int i=0; i<_graph.size(); i++) {
		  		H->add_node(_graph.node[i]);
		  		H->add_node(_graph.node[i]+_graph.capacity);
		  		H->set_left(_graph.node[i]);
		  		H->set_right(_graph.node[i]+_graph.capacity);
		  	}
		  	
		  	for (int i=0; i<_graph.size(); i++) {
		  		for (size_t j=0; j<_graph.neighbor[_graph.node[i]].size; j++) {
		  			H->add_undirected(_graph.node[i], _graph.neighbor[_graph.node[i]][j]+_graph.capacity);
				}
			}
			
			H->vertex_cover();
			
			for (size_t i=0; i<H->left.size; i++) {
				v = H->left[i];
				if ((H->vcover.contain(v)) && (H->vcover.contain(v+_graph.capacity))) {
					//add_to_cover(H->node[i], NO_REASON);
					add_to_cover(v, NO_REASON);
				} else if ((!H->vcover.contain(v)) && (!H->vcover.contain(v+_graph.capacity))) {
					save(v, NO_REASON);
					_graph.rem_and_update(v);
					//remove_and_save(v, NO_REASON);
				}
			}
			
			delete H;
			
			// check the number of vertices
			if (_graph.size() > 2*(_upper_bound-(int)(_cover.size))) return false;
			
			return true;	
		  		
		}
		
		bool kernelize_NT_strong()
		{
		
			int v;
			
		  	// Generate the bipartite graph
		  	BipartiteGraph *H = new BipartiteGraph(_graph.capacity*2);
		  	
		  	for (int i=0; i<_graph.size(); i++) {
		  		H->add_node(_graph.node[i]);
		  		H->add_node(_graph.node[i]+_graph.capacity);
		  		H->set_left(_graph.node[i]);
		  		H->set_right(_graph.node[i]+_graph.capacity);
		  	}
		  	
		  	for (int i=0; i<_graph.size(); i++) {
		  		for (size_t j=0; j<_graph.neighbor[_graph.node[i]].size; j++) {
		  			H->add_undirected(_graph.node[i], _graph.neighbor[_graph.node[i]][j]+_graph.capacity);
				}
			}
			
			H->persistently_matched_left_vertices();
		
			for (size_t i=0; i<H->left.size; i++) {
				v = H->left[i];
				if (!H->lpmatched.contain(v)) {
					while (_graph.neighbor[v].size != 0) {
						add_to_cover(_graph.neighbor[v][0], NO_REASON);
					}
					save(v, NO_REASON);
					_graph.rem_and_update(v);
				}
			}
			
			// check the number of vertices
			if (_graph.size() > 2*(_upper_bound-(int)(_cover.size))) return false;
			
			delete H;
			
			return true;
		}

		
		int greedy_lower_bound(const int limit=-1)
		{
			int lb=0, i, k=0, covered=0, m=_graph.num_edges, l, c, stop=limit;
			
			for(int d=_graph.max_degree; d>=_graph.min_degree; --d) {
				for(int j=_graph.node_of_degree[d].size; --j>=0;) {		
					degree_count[k++] = d;
				}
			}
			
			if(stop<0) stop = k;
			--k;
		
			while(covered<m) {
				c = degree_count[lb];
				if(++lb>=stop) break;
				
				covered += c;
				
				l = k-c;	
				while(k>l && --degree_count[k]==0) --k;
			 	for(i=k; --i>l;) --degree_count[i];
			}

			
			return _cover.size+lb;
		}
		

		
		
		int cliques_lower_bound(const int limit=-1)
		{
			int lb=0;
			int d, j, i, x, best;
	

			for(i=0; i<num_cliques; ++i) {
				candidates[i].fill();
#ifdef _USE_CLIQUES
				clique_set[i].clear();
				clique_lst[i].clear();
#endif
			}
			num_cliques=0;
			int smallest_clique=0;

			for(d=_graph.min_degree; d<=_graph.max_degree; ++d) {
				for(j=_graph.node_of_degree[d].size; --j>=0;) {
					x = _graph.node_of_degree[d][j];
					best = -1;										
					for(i=smallest_clique; i<num_cliques; ++i) {
						if(candidates[i].fast_contain(x)) {
							best=i;
							break;
						}
					}

					if(best<0) {
						best = num_cliques;
						++num_cliques;
						//while(smallest_clique<num_cliques && candidates[smallest_clique].empty()) ++smallest_clique;
					} else {
							//if(++lb==limit) break;
						++lb;
					}

					candidates[best].intersect_with(_graph.matrix[x]);

					
#ifdef _USE_CLIQUES
					clique_set[best].add(x);
					clique_lst[best].add(x);
#endif
				}
			}			

			return _cover.size+lb;
		}
		
		
		
		int cliques_dual_bound(const int limit=-1)
		{
	
			int i, x;
	
			for(i=0; i<num_cliques; ++i) {
				candidates[i].fill();
#ifdef _USE_CLIQUES
				clique_set[i].clear();
				clique_lst[i].clear();
#endif
			}
			
			num_cliques=0;
			util_set.copy(_graph.node_set);
		
			//while(!util_set.empty()) {
			for(i=_graph.size(); i>0;) {
				candidates[num_cliques].copy(util_set);
				while(!candidates[num_cliques].empty()) {
					x = candidates[num_cliques].min();
					util_set.remove(x);
#ifdef _USE_CLIQUES
					clique_set[num_cliques].add(x);
					clique_lst[num_cliques].add(x);
#endif
					candidates[num_cliques].intersect_with(_graph.matrix[x]);
					
					--i;
				}
				++num_cliques;
			}
			
			suggested_node = x;
	
			return _cover.size+(_graph.size()-num_cliques);
		}
		
		
		bool clique_dominance() {
#ifdef _VERIF_TRAIL	
			verify_trail("before dominance");
#endif			
			int i, j, k, x, ds, reason = DEDUCTION;
			bool dominance = false;
			
			for(i=0; i<num_cliques; ++i) {
				//util_set.copy(candidates[i]);
				//util_set.intersect_with(_graph.node_set);
				if(!candidates[i].intersect(_graph.node_set)) {
					// possibility that a node has all its neighbors in the clique
					for(j=0; j<(int)(clique_lst[i].size); ++j) {
						x = clique_lst[i][j];
						util_set.copy(_graph.matrix[x]);
						util_set.intersect_with(_graph.node_set);
						if(util_set.included(clique_set[i])) {
							// x has all its current neighbors in the same clique, its dominated
							
							if(use_backjump) {
								
#ifdef _DEBUG_BACKJUMPS
								if(_DEBUG_BACKJUMPS) {	
									std::cout << "build nogood for clique dominance:";
								}
#endif								
								
								reason = _nogood.size;
								ds = _graph.original_size[x];
								_nogood.add(ds-_graph.degree(x));
								for(k=_graph.degree(x); k<ds; ++k) {
									
#ifdef _DEBUG_BACKJUMPS
									if(_DEBUG_BACKJUMPS) {	
										std::cout << " " << _graph.neighbor[x][k] ;
									}
#endif										
									
									_nogood.add(_graph.neighbor[x][k]);
								}
								
#ifdef _DEBUG_BACKJUMPS		
								if(_DEBUG_BACKJUMPS) {						
									std::cout << std::endl << "   -> explains " << x << " and its neighbors " << _graph.neighbor[x] << std::endl;
								}
#endif
								
							}
							
							_graph.rem_and_update(x);
							//save(x, reason);
							//remove_and_save(x, reason);
							exclude_from_cover(x, reason);
							dominance = true;
							
							break;
						}
					}
				}
			}
#ifdef _VERIF_TRAIL			
			verify_trail("after dominance");
#endif			
			
			return dominance;
		}
		
		
		
		
#ifdef _USE_CLIQUES		
		void print_cliques() 
		{
			for(int i=0; i<num_cliques; ++i) {
				std::cout << clique_set[i] << "\t\t\t" << candidates[i] << std::endl;
			}
		}
#endif		
		
		void print_trail(const int full=false) { 
			std::cout << std::setw(5) << nb_nodes << ":"; 
						
			int kthdec = _decision.size-1;
			int kthcov = _cover.size-1;
			int kthnod;
			for(kthnod=_trail.size-1; kthnod>=0; --kthnod) {
				if(kthdec>=0 &&_trail[kthnod] == _decision[kthdec]) {
					--kthdec;
					if(kthcov>=0 && _trail[kthnod] == _cover[kthcov]) {
						std::cout << " [" << _trail[kthnod] << "]";
						--kthcov;
					} else  {
						std::cout << " [-" << _trail[kthnod] << "]";
					}
						
				} else if(kthcov>=0 && _trail[kthnod] == _cover[kthcov]) {
					std::cout << " +" << _trail[kthnod];
					--kthcov;
				} else {
					std::cout << " -" << _trail[kthnod];
				}
			}
			std::cout << std::endl;
			
			if(full){
				
				std::cout << "\n                  Decisions: " << _decision << std::endl;
				std::cout << " Analyze conflict with trail " ; //_cover << std::endl;
				for(int i=_trail.size-1; i>=0; --i)
					std::cout << " " << std::setw(2) << _trail[i];
				std::cout << std::endl << "                      Level: " ;
				for(int i=_trail.size-1; i>=0; --i)
					std::cout << " " << std::setw(2) << _level[_trail[i]];
				std::cout << std::endl << "                       Rank: " ;
				for(int i=_trail.size-1; i>=0; --i)				
					std::cout << " " << std::setw(2) << _rank[_trail[i]];
				std::cout << std::endl << "                     Reason: " ;
				for(int i=_trail.size-1; i>=0; --i)								
#ifdef _STORE_NOGOOD	
					std::cout << " " << (_reason[_trail[i]] ? "  " : " *");
#else
				std::cout << " " << std::setw(2) << _reason[_trail[i]];
#endif
				std::cout << std::endl;
				
			}
			
		}
		
		
		bool backtrack() {
				
			++nb_conflicts;

			if((int)(_decision.size)>init_level) {
				
				++backjump_size;
				
				int e = undo();
				if(STATUS(e)==COV) {
					exclude_from_cover(NODE(e));
				} else {
					add_to_cover(NODE(e));
				}
				
#ifdef _TRACE_
			print_trail();
#endif
					
				return true;
			}
			return false;
		}
		
		
		bool backjump() {
			
#ifdef _VERIF_TRAIL
			verify_trail("before backjump");
#endif			

#ifdef _DEBUG_BACKJUMPS
			if(_DEBUG_BACKJUMPS) {
				print_trail(true);
			}
#endif			
			
			
			++nb_conflicts;	
			if((int)(_decision.size)<=init_level) {
				return false;
			}

			// first, compute a minimal explanation for the last fail
			int i, j, x, y, z, u;
			_explanation.quick_clear();
			
			
			for(i=_cover.size; i<_graph.capacity; ++i) {
				x = _cover[i];
				if(!_graph.node.contain(x)) {
					candidates[num_cliques].copy(_graph.matrix[x]);
					++num_cliques;	
				}	
			}

			int slack = _lower_bound-_upper_bound;
			//if(slack>0) std::cout << "slack = " << slack << std::endl;
			
			for(i=_cover.size-1; i>=0; --i) {
				x = _cover[i];
					
				for(j=0; j<num_cliques; ++j) {
					if(candidates[j].contain(x)) {
						candidates[j].intersect_with(_graph.matrix[x]);
						break;
					}
				}
				if(j>=num_cliques) {
					_explanation.add_end(x);		
				} else if(slack>0) {
					candidates[num_cliques].copy(_graph.matrix[x]);
					++num_cliques;
					--slack;
				}
			}
			

#ifdef _DEBUG_BACKJUMPS
			if(_DEBUG_BACKJUMPS) {
				if(_cover.empty()) {
					std::cout << "Extract contradiction!" << std::endl;
				} else {
					std::cout << " Extract a minimal conflict " << _explanation << std::endl
										<< "                 from cover [" << _cover[_cover.size-1]; //_cover << std::endl;
					for(int k=_cover.size-2; k>=0; --k) {
						std::cout << " " << _cover[k];
					}
					std::cout << "]" << std::endl;
				}
			}
#endif
			
			// now use resolution to remove all but one of this level's literals
			if(_explanation.size == 0) {
				
#ifdef _NOGOOD_STAT
				++backjump_sizes[_decision.size];
				++nogood_sizes[0];
#endif
								
				backjump_size += _decision.size;
				return false;
			} else {
				// resolution
	
				x = _explanation.get_first();
				y = _explanation.next[x];
				
#ifdef _STORE_NOGOOD
				Array<int>* new_nogood = NO_REASON;
#else
				int new_nogood = NO_REASON;
#endif
				
				int nbacktracks = _decision.size;
				
#ifdef _FORGET_NOGOOD
				int min_reason = _nogood.size;
#endif				
				
				if(y != _explanation.capacity) {
					while(_level[y]==_level[x]) {
						// resolve x since there are at least one other literal from the same level
							
						z = y;				
						_explanation.remove(x);
						
#ifdef _STORE_NOGOOD
						Array<int>& reason(*(_reason[x]));
						int   rsize = reason.size;
#else

#ifdef _DEBUG_BACKJUMPS							
						if(_reason[x] < 0) {
							std::cout << "node " << x << " level " << _level[x] << " rank " << _rank[x] << std::endl;
							std::cout << nb_nodes << std::endl;
						}
						
						assert(_reason[x]>=0);
						assert(_reason[x]<_nogood.size);
#endif
												
						int* reason = _nogood.stack_+_reason[x]+1;
						int   rsize	= _nogood[_reason[x]];
						
#ifdef _FORGET_NOGOOD
						min_reason = _reason[x];
#endif

#endif
						
						
#ifdef _DEBUG_BACKJUMPS
						if(_DEBUG_BACKJUMPS) {
#ifdef _STORE_NOGOOD
						std::cout << " resolve " << x << " with its explanation " << reason << std::endl;
#else
						std::cout << " resolve " << x << " with its explanation" ;
						for(int k=0; k<rsize; ++k)
							std::cout << " " << reason[k];
						std::cout << std::endl;
#endif
					}
#endif	

						for(i=0; i<rsize; ++i) {
							u = reason[i];
							while(_rank[z]>_rank[u]) {
								z = _explanation.next[z];
							}
							if(_rank[z]<_rank[u]) { // otherwise z and u are equal
								_explanation.add_before(u,z);
							} 
						
						}
	
						x = _explanation.get_first();
						y = _explanation.next[x];
						
						if(_level[x] == 0) {
							
#ifdef _NOGOOD_STAT
							++backjump_sizes[_decision.size];
							++nogood_sizes[0];
#endif
							
							backjump_size += _decision.size;
							return false;
						}
					}
				
					// x is the last literal from this level and y is the literal of highest level otherwise
					_explanation.remove(x);

#ifdef _STORE_NOGOOD
					new_nogood = Array<int>::Array_new(_explanation);
					_nogood.add(new_nogood);
#else
					
#ifdef _FORGET_NOGOOD
					_nogood.size = min_reason;
#endif
					
#ifdef _DEBUG_BACKJUMPS
					if(_DEBUG_BACKJUMPS) {	
						std::cout << " new backjump nogood:" ;
					}
#endif				
															
					new_nogood = _nogood.size;
					_nogood.add(_explanation.size);
					int l = _explanation.get_first();
					do {
						_nogood.push_back(l);
						l = _explanation.next[l];
					} while(l < _explanation.capacity);

#endif
					
					if(y<_explanation.capacity)
						nbacktracks -= _level[y];
					
				}	
				
#ifdef _DEBUG_BACKJUMPS
				if(_DEBUG_BACKJUMPS) {
					std::cout << "                            " << _explanation << std::endl;
					std::cout << " backjump " << nbacktracks << " levels (of " << _decision.size << ") and exclude " << x << std::endl << std::endl;
				}
#endif
				
				
				nogood_size += _explanation.size;
				backjump_size += nbacktracks;
#ifdef _NOGOOD_STAT
				++backjump_sizes[nbacktracks];
				++nogood_sizes[_explanation.size];
#endif
				
				bool was_in_cover = _cover.contain(x);
				

					//int z;

				//print_trail();

				y = undo();

				
				//std::cout << "undone until decision " << y << std::endl;
				
				while(--nbacktracks>0) {

					//print_trail();

					z = NODE(y);							
					//z = (y<0 ? -y-1 : y);
					
					//std::cout << "add " << z << " into " << std::endl << _graph << std::endl;
					//std::cout << "add " << z << " into G" << std::endl;
					
					_graph.add_and_update(z);		
					
					//std::cout << _graph << std::endl;
					
					y = undo();	
					
					//std::cout << "undone until decision " << y << std::endl;
				}	
				
				z = NODE(y);
				//z = (y<0 ? -y-1 : y);		
				
				//print_trail();
				
				//std::cout << "last decision = " << z << " / nogood = " << x << std::endl; 
				
				if(z!=x) {
					
					//std::cout << "add " << z << " into G" << std::endl;
					
					_graph.add_and_update(z);
					
					//std::cout << "rem " << x << " from G" << std::endl;
					
					_graph.rem_and_update(x);
				}
				if(was_in_cover){
					exclude_from_cover(x, new_nogood);
				} else {
					add_to_cover(x, new_nogood);
				}
				
				// if(y!=x) {
				// 	restore();
				// 	remove_and_save(x, new_nogood);
				// } else {
				// 	// trick to avoid popping x out of the trail and back in again
				// 	_reason[x] = new_nogood;
				// 	_level[x] = _decision.size;
				// }
				//
				// exclude_from_cover(x, new_nogood);	

			}
			
#ifdef _VERIF_TRAIL					
				verify_trail("after backjump");
#endif			

			return true;
		}
	

		
// #define _DEBUG_COLORING
// #define _VERIFY_COLORING

		int complement_coloring() {
				
			int x, y, n=_graph.capacity-1, max_col=-1;
			//std::fill(_coloring, _coloring+n+1, -1);
				
					
#ifdef _DEBUG_COLORING				
			std::cout << "\nCompute coloring in complement of\n" << _graph << "\nwith colors " << color_list << std::endl;
#endif
				
			for(int d=_graph.min_degree; d<=_graph.max_degree; ++d) {
				for(size_t i=0; i<_graph.node_of_degree[d].size; ++i) {
					x = _graph.node_of_degree[d][i];
						
#ifdef _DEBUG_COLORING								
					std::cout << "color " << x << " with";
#endif
						
					for(int j=_graph.original_size[x]; j<n; ++j)
					{
						y = _graph.neighbor[x][j];
							
						if(_coloring[y]>=0 && color_list.contain(_coloring[y])) {
							color_list.remove(_coloring[y]);

#ifdef _DEBUG_COLORING										
							std::cout << std::endl << y << "=[" << _coloring[y] << "]"<< " -> " << color_list;
#endif
																
						}
					}
						
					_coloring[x] = color_list.get_first();

#ifdef _DEBUG_COLORING								
					std::cout << " => " << _coloring[x] << " smallest in " << color_list << std::endl;
#endif			
									
					color_list.fill();
						
					if(max_col<_coloring[x]) max_col = _coloring[x];
						
				}
			}
				
#ifdef _DEBUG_COLORING						
			std::cout << "=> coloring of size " << (max_col+1) << std::endl;
#endif
							

#ifdef _VERIFY_COLORING
			for(size_t i=0; i<_graph.node.size; ++i) {
				x = _graph.node[i];
				for(int j=_graph.neighbor[x].size; j<_graph.capacity-1; ++j) {
					y = _graph.neighbor[x][j];
						
					if(_graph.node.contain(y) && _coloring[x] == _coloring[y]) {
							
							
						std::cout << " " << x << " and " << y << " cannot be in the same clique since they aren't neighbors!" << std::endl;
						exit(1);
					}
				}
					
			}
#endif					
							
			for(size_t i=0; i<_graph.node.size; ++i) {
				_coloring[_graph.node[i]] = -1;
			}			
								
			return max_col+1;
		}
			
			
			
    // push on the trail			
		inline void save(const int x, 
#ifdef _STORE_NOGOOD				
		Array<int>* reason=NO_REASON
#else
		const int reason=NO_REASON
#endif			
		) {
			
// #ifdef _VERIF_TRAIL
// 			verify_trail("before rem and save");
// #endif
			
			_rank[x] = _trail.size;
			_trail.add(x);
			_reason[x] = reason;
			_level[x] = _decision.size;
			//_graph.rem_and_update(x);
			
// #ifdef _VERIF_TRAIL
// 			verify_trail("after rem and save");
// #endif
			
		}
		
		// pop from the trail
		inline int restore() {
			
			int x = _trail.pop();
			//if(_cover.back()==x) _cover.pop();
				
			// _graph.add_and_update(x);
			return x; //_trail.pop();
						
		}
		
		

			
// 		inline void remove_and_save(const int x,
// #ifdef _STORE_NOGOOD
// 		Array<int>* reason=NO_REASON
// #else
// 		const int reason=NO_REASON
// #endif
// 		) {
//
// #ifdef _VERIF_TRAIL
// 			verify_trail("before rem and save");
// #endif
//
// 			// _rank[x] = _trail.size;
// 			// _trail.add(x);
// 			// _reason[x] = reason;
// 			// _level[x] = _decision.size;
// 			save(x,reason);
// 			_graph.rem_and_update(x);
//
// #ifdef _VERIF_TRAIL
// 			verify_trail("after rem and save");
// #endif
// 		}
			
		inline void add_to_cover(const int x, 
#ifdef _STORE_NOGOOD				
		Array<int>* reason=NO_REASON
#else
		const int reason=NO_REASON
#endif				
		) {
			
// #ifdef _VERIF_TRAIL
// 			verify_trail("before add to cover");
// #endif		
			
			//std::cout << "add " << x << " b/c " << reason << std::endl;
			
			//remove_and_save(x, reason);
			save(x, reason);
			_cover.add(x);
			
// #ifdef _VERIF_TRAIL
// 			verify_trail("after add to cover");
// #endif	
			
		}
			
		inline void exclude_from_cover(const int x, 
#ifdef _STORE_NOGOOD				
		Array<int>* reason=NO_REASON
#else
		const int reason=NO_REASON
#endif				
		) {
		
// #ifdef _VERIF_TRAIL
// 			verify_trail("before exclude from cover");
// #endif	
			
			//std::cout << "exclude " << x << " b/c " << reason << std::endl;
			
			save(x, reason);
			int y;
			for(int i=_graph.neighbor[x].size; --i>=0;) {	
				y = _graph.neighbor[x][i];
				add_to_cover(y, reason);
				_graph.rem_and_update(y);
			}
			
// #ifdef _VERIF_TRAIL
// 			verify_trail("after exclude from cover");
// #endif		
			
		}
		
		// inline void save() {
		// 	_decision.add(-1);
		// }
					


		inline int undo() {
				
// #ifdef _VERIF_TRAIL
// 			verify_trail("before undo");
// #endif		
					

				
				// << _decision << std::endl
			// 		<< _cover << std::endl
			// 			<< _trail << std::endl;
					
					
			// int x, y = _decision.pop();
			// do {
			// 	x = _cover.pop();
			// 	std::cout << " " << x ;
			// } while(x != y);
			// std::cout << std::endl;
			
						//	std::cout << _decision.size << std::endl;
			
			
			
				int x = restore(), y = _decision.pop();
				
				// std::cout << "restore " << x << std::endl;
				//
				// std::cout << "decision = " << y << std::endl;
				
				int l = _level[y];
				
				// std::cout << l << std::endl;
				
				while(x!=y) {
					
					// std::cout << " add " << x << " to" << std::endl << _graph << std::endl;
					
					_graph.add_and_update(x);
					
					// std::cout << "restore" << std::endl;

					x = restore();
					
					// std::cout << "restored " << x << std::endl;
					//if(x==_cover.back()) _cover.pop();
				}
				
				// std::cout << "restore cover " << _cover << std::endl;

				do {
					x = _cover.pop();					
				} while(!_cover.empty() && _level[_cover.back()]>=l);
				
			
			
			
			//
			//
			// int x, y = _decision.pop();
			// //std::cout << y << std::endl;
			// int l = _level[y];
			// std::cout << l << std::endl;
			//
			// while(!_cover.empty() && _level[_cover.back()] >= l) {
			// 	std::cout << "pop " << _cover.back() << " " << _level[_cover.back()] << std::endl;
			// 	//std::cout.flush();
			// 	x = _cover.pop();
			// }
			//
			//
			// while(_trail.back() != y)
			// {
			// 	std::cout << "restore " << _trail.back() << std::endl;
			//
			// 	restore();
			// }
			//
			// //std::cout << _trail << std::endl; //std::endl;
	
// #ifdef _VERIF_TRAIL
// 			verify_trail("after undo");
// #endif
			
				
			//return (x==y ? y : -(y+1));
				return (x==y ? IN_COVER(y) : IN_STABLE(y));
			//return ()
		}
			
// 		inline bool suboptimal(const bool dominance) {
//
// #ifdef _VERIF_TRAIL
// 			verify_trail("before consistency");
// #endif
//
// 			int lb = _cover.size;
//
// 			if(_graph.num_edges) {
// 				apply_buss();
//
// 				lb = cliques_lower_bound(_upper_bound-_cover.size);
//
// 				if(lb<_upper_bound && dominance) {
// 					clique_dominance();
// 					lb = cliques_lower_bound(_upper_bound-_cover.size);
// 				}
// 			}
//
// #ifdef _VERIF_TRAIL
// 			verify_trail("after consistency");
// #endif
//
// 			return lb >= _upper_bound;
// 		}
		
		inline bool suboptimal(const bool dominance) {

#ifdef _VERIF_TRAIL
			verify_trail("before consistency");
#endif

			_lower_bound = _cover.size;

			if(_graph.num_edges) {
				apply_buss();
				kernelize_cliques();


				_lower_bound = cliques_lower_bound(_upper_bound-_cover.size);
				
				// if(_lower_bound>_upper_bound) {
				// 	std::cout << "overhead= " << (_lower_bound-_upper_bound) << std::endl;
				// }
				
			}

#ifdef _VERIF_TRAIL
			verify_trail("after consistency");
#endif
			
#ifdef _TRACE_
			std::cout << "bounds = [" << _lower_bound << ".." << _upper_bound << "]" << std::endl;
#endif

			return _lower_bound >= _upper_bound;
		}
		
// 		inline bool suboptimal(const bool dominance) {
//
// #ifdef _VERIF_TRAIL
// 			verify_trail("before consistency");
// #endif
//
// 			int lb = _cover.size;
//
// 			if(_graph.num_edges) {
// 				apply_buss();
//
// 				lb = cliques_lower_bound(_upper_bound-_cover.size);
//
// 				while(lb<_upper_bound && dominance && (clique_dominance() || apply_buss())) {
// 					lb = cliques_lower_bound(_upper_bound-_cover.size);
// 				}
// 			}
//
// #ifdef _VERIF_TRAIL
// 			verify_trail("after consistency");
// #endif
//
// 			return lb >= _upper_bound;
// 		}
			
			
		inline bool apply_buss() {
				
#ifdef _VERIF_TRAIL					
			verify_trail("before buss");
#endif				
			
			bool pruning = false;
			if(_graph.max_degree >= (_upper_bound - (int)(_cover.size))) {
				int x, reason = DEDUCTION;
					
				if(use_backjump) {
						
#ifdef _DEBUG_BACKJUMPS
					if(_DEBUG_BACKJUMPS) {	
						std::cout << "build nogood from buss rule:";
					}
#endif							
						
					reason = _nogood.size;
					_nogood.add(_cover.size);
					for(int i=_cover.size-1; i>=0; --i) {
					// _nogood.add(_decision.size);
					// for(int i=_decision.size-1; i>=0; --i) {
						_nogood.add(_cover[i]);
							
#ifdef _DEBUG_BACKJUMPS
						if(_DEBUG_BACKJUMPS) {	
							std::cout << " " << _cover[i];
						}
#endif	
							
					}	
						
#ifdef _DEBUG_BACKJUMPS
					if(_DEBUG_BACKJUMPS) {	
						std::cout << std::endl << "   -> explains";
					}
#endif					
						
				}
				do { 
					
					//std::cout << "buss " << _graph.min_degree << " .. " << _graph.max_degree << std::endl;
					
					
					x = _graph.node_of_degree[_graph.max_degree].back();
					_graph.rem_and_update(x);
					add_to_cover(x, reason);
					
#ifdef _DEBUG_BACKJUMPS
					if(_DEBUG_BACKJUMPS) {	
						std::cout << " " << x;
					}
#endif	
				}			
				while(_graph.num_edges && _graph.max_degree >= (_upper_bound - (int)(_cover.size)));
				
#ifdef _DEBUG_BACKJUMPS
				if(_DEBUG_BACKJUMPS) {	
					std::cout << std::endl ;
				}
#endif
					
				pruning = true;
			}

#ifdef _VERIF_TRAIL					
			verify_trail("after buss");
#endif

			return pruning;
		}
			
		inline void branch_on_maxdegree() {
			
#ifdef _VERIF_TRAIL
			verify_trail("before branch");
#endif			
				
			int x = _graph.node_of_degree[_graph.max_degree].back();				
			_decision.add(x);
			_graph.rem_and_update(x);
			add_to_cover(x, NO_REASON);
				
				
#ifdef _TRACE_
			print_trail();
#endif
			
#ifdef _VERIF_TRAIL
			verify_trail("after branch");
#endif
				
		}
		
		inline void branch_on_mindegree() {
			
#ifdef _VERIF_TRAIL
			verify_trail("before branch");
#endif			
				
			int x = _graph.node_of_degree[_graph.min_degree].back();				
			_decision.add(x);
			//add_to_cover(x, NO_REASON);
			//remove_and_save(x, NO_REASON);
			_graph.rem_and_update(x);
			int reason = DEDUCTION;
			if(use_backjump) {
				reason = _nogood.size;
				_nogood.add(1);
				_nogood.add(x);
			}
			exclude_from_cover(x, reason);
			///HACK!!
			_reason[x] = NO_REASON;
				
				
#ifdef _TRACE_
			
			
			// for(int i=_cover.size-1; i>=0; --i) {
			// 	std::cout << " " << _cover[i];
			// }
			// std::cout << std::endl;
			
			int kthdec = _decision.size-1;
			int kthcov = _cover.size-1;
			int kthnod;
			for(kthnod=_trail.size-1; kthnod>=0; --kthnod) {
				if(kthdec>=0 &&_trail[kthnod] == _decision[kthdec]) {
					--kthdec;
					if(kthcov>=0 && _trail[kthnod] == _cover[kthcov]) {
						std::cout << " [" << _trail[kthnod] << "]";
						--kthcov;
					} else  {
						std::cout << " [-" << _trail[kthnod] << "]";
					}
							
				} else if(kthcov>=0 && _trail[kthnod] == _cover[kthcov]) {
					std::cout << " +" << _trail[kthnod];
					--kthcov;
				} else {
					std::cout << " -" << _trail[kthnod];
				}
				std::cout.flush();
			}
			std::cout << std::endl;
			
			// int kthdec = _decision.size-1;
			// int kthcov = _cover.size-1;
			// int kthnod;
			// for(kthnod=_trail.size-1; kthnod>=0; --kthnod) {
			// 	if(kthdec>=0 &&_trail[kthnod] == _decision[kthdec]) {
			// 		std::cout << " [" << _trail[kthnod] << "]";
			// 		--kthdec;
			// 		--kthcov;
			// 	} else if(kthcov>=0 &&_trail[kthnod] == _cover[kthcov]) {
			// 		std::cout << " +" << _trail[kthnod];
			// 		--kthcov;
			// 	} else {
			// 		std::cout << " -" << _trail[kthnod];
			// 	}
			// }
			// std::cout << std::endl;
#endif
			
#ifdef _VERIF_TRAIL
			verify_trail("after branch");
#endif
				
		}
			
		inline void reset() {
			// if(_decision.size>init_level) {
			// 	int x;
			// 	do {
			// 		x = undo();
			// 	} while(_decision.size>init_level);
			// 	_graph.add_and_update((x < 0 ? -x-1 : x));
			// }
			
			//std::cout << "reset" << std::endl;
			
			int x=0;
			
			//print_trail();
			
			while((int)(_decision.size)>init_level) {
				x = undo();
				
				// print_trail();
				//
				// std::cout << "add " << (x < 0 ? -x-1 : x) << " to" << std::endl << _graph << std::endl;
				//
				_graph.add_and_update((x < 0 ? -x-1 : x));
				
				//std::cout << _graph << std::endl;
				
		  }
			
			// std::cout << _graph << std::endl;
			//
			//
			//
			// std::cout << _graph << std::endl;
			//
			// //if(_trail.size && _level[_trail.back()]>init_level) restore();
		}
			
		void branch_on_maxdegree_rand();
			
		void print_statistics();
			
		void print_summary();
			
		void print_head();

		void new_cover();
			
		bool run(const bool bmax=true, const bool dominance=true, const int limit=-1, const int ub_limit=-1, const bool randomized=false);
			
		void solve(const bool bmax=true, const bool use_bkjp=true, const bool dominance=true, const bool kernelize=false, const bool randomized=false, const int base=-1, const double factor=1.1);

		void lns(const int limit=100, const double factor=.7);

		


	};

}



		








// class IntervalList {

// public:

//   unsigned int size;
//   Vector<int> min;
//   Vector<int> max;

//   bool contain(const int v) {
//     int lmin=0;
//     int lmax=0;

//     bool is_in = true;

//     while(is_in) {
//       if(v < min[lmin] || v > max[lmax]) {
// 	is_in = false;
//       }
//     }
//   }

// };

// // [0,10][23,25][40,100][120,150]
// //    min   max
// /*   -inf   inf
// 	0   150
//        23    10
//       120   100
//        40    25
//       inf  -inf
//  */




#endif // __STRUCTURE_HPP
