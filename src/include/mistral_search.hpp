
/*
  Mistral 2.0 is a constraint satisfaction and optimisation library
  Copyright (C) 2009  Emmanuel Hebrard
  
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


#ifndef __SEARCH_HPP
#define __SEARCH_HPP


#include <mistral_global.hpp>
#include <mistral_structure.hpp>
#include <mistral_constraint.hpp>


namespace Mistral {


  class Solver;
  //class Search;
  class RestartPolicy {
    
  public:

    unsigned int base;
    
    RestartPolicy(const unsigned int b=256);
    virtual void reset(unsigned int& limit) = 0;
    
  };


  class NoRestart : public RestartPolicy {
    
  public:
    
    NoRestart();
    virtual ~NoRestart();
    
    void reset(unsigned int& limit) {
      limit = base;
    }
    
  };


  class Geometric : public RestartPolicy {
    
  public:
    
    unsigned int increment;
    double factor;

    Geometric(const unsigned int b=256, const double f=1.333);
    virtual ~Geometric();
    
    void reset(unsigned int& limit) {
      limit += increment;
      increment = (unsigned int)((double)increment * factor);
    }
    
  };

  class Luby : public RestartPolicy {

  private:
    
    unsigned int luby_seq(const int iter) {
      unsigned int thelog = log2(iter);
      if( iter == (1 << (thelog + 1))-1 )
	return (1 << thelog);
      return luby_seq(iter - (1 << thelog) + 1);
    }
    
    
  public:
    
    unsigned int iteration;

    Luby(const unsigned int b=100);
    virtual ~Luby();
    
    void reset(unsigned int& limit) {
      limit += (base * luby_seq(++iteration));
    }
    
  };

  class VarOrdering {

  public: 

    Solver *solver;
    unsigned int& length;
    IntVar *variables;

    VarOrdering(Solver *s);
    virtual ~VarOrdering();
    
    virtual IntVar select() = 0;

  };

  class NoOrder : public VarOrdering {

  public: 

    NoOrder(Solver *s);
    virtual ~NoOrder();
    
    virtual IntVar select();

  };


  /**********************************************
   * Generic heuristic
   **********************************************/
  
  template < class Selector >
  class GenericDVO : public VarOrdering
  {
  public: 
    
    /**@name Parameters*/
    //@{ 
    Selector best;
    Selector current;
    //@}

    /**@name Constructors*/
    //@{
    GenericDVO(Solver* s) : VarOrdering(s) {}
    //@}
    
    /**@name Utils*/
    //@{ 
    inline IntVar select()
    {    
      IntVar var = variables[0];
      best = var;
      for(unsigned int i=1; i<length; ++i) {
	current = variables[i];
	if( current < best ) {
	  best = current;
	  var = variables[i];
	} 
      }
      return var;
    }
    //@}
  };


  /**********************************************
   * GenericRandom heuristic
   **********************************************/
  
  template < class Selector >
  class GenericRandomDVO : public VarOrdering
  {
  public: 

    /**@name Parameters*/
    //@{ 
    Selector *bests;
    Selector current;
    IntVar *bestvars;
    int size;
    //@}

    /**@name Constructors*/
    //@{
    GenericRandomDVO(Solver* s, const int sz)  : VarOrdering(s) 
    {
      size = sz;
      bests = new Selector[size+1];
      bestvars = new IntVar[size+1];
    }

    virtual ~GenericRandomDVO() 
    {
      delete [] bests;
      delete [] bestvars;
    }
    //@}

    /**@name Utils*/
    //@{ 
    inline IntVar select()
    {
      unsigned int realsize=1, i, j;
      bests[0] = bestvars[0] = variables[0];
      for(j=1; j<length; ++j)
	{  
	  current = variables[j];
	  i = realsize;
	  while( i && current < bests[i-1] ) {
	    bests[i] = bests[i-1];
	    bestvars[i] = bestvars[i-1];
	    --i;
	  }
	  bests[i] = current;
	  bestvars[i] = variables[j];
	  
	  if(realsize<size) ++realsize;
	}
      return bestvars[randint(realsize)];
    }
    //@}
  };



  
//   class Search {

//   public:
//     Stack< IntVar > sequence;

//     Vector< unsigned int > trail_;
//     Vector< IntVar > decision;

//     VarOrdering *heuristic;   
//     RestartPolicy *policy;   
    
//     SolverParameters parameters;
//     SolverStatistics statistics;

//     Search() {
//       heuristic = NULL;
//       policy = NULL;
//     }

//     void initialise(Solver *s);
//     IntVar backtrack(); 
    
//     IntVar make_node() {
//       ++statistics.num_nodes;
//       trail_.add(sequence.size);
//       decision.add(heuristic->select());
//       return decision.back();
//     }

//     void assign(IntVar x);
//     Outcome satisfied();
//     bool limitsExpired();

//     void init_search(Vector< IntVar >& seq, VarOrdering *h, RestartPolicy *p);
//   };


//   inline Outcome Search::satisfied() {    
// #ifdef _DEBUG_SEARCH
//     std::cout << "c";
//     for(unsigned int k=0; k<=decision.size; ++k) std::cout << " ";
//     std::cout << " SAT!" << std::endl; 
// #endif
    
//     ++statistics.num_solutions;
    
//     return SAT;
//   }
  
// //   inline bool Search::limitsExpired() {

// //     return (parameters.limit && 
// // 	    ((parameters.time_limit > 0.0 && (getRunTime() - statistics.start_time) > parameters.time_limit) ||
// // 	     (parameters.node_limit > 0 && (statistics.num_nodes > parameters.node_limit)) ||
// // 	     (parameters.fail_limit > 0 && (statistics.num_failures > parameters.fail_limit)) ||
// // 	     (parameters.restart_limit > 0 && (statistics.num_restarts > parameters.restart_limit)) ||
// // 	     (parameters.backtrack_limit > 0 && (statistics.num_backtracks > parameters.backtrack_limit))
// // 	     ));
// //   }

}

#endif // __SEARCH_HPP
