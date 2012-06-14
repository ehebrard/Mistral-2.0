
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


/*! \file mistral_variable.hpp
  \brief Header for the variables
*/

#include <vector>

#include <mistral_backtrack.hpp>


#ifndef _MISTRAL_VARIABLE_HPP
#define _MISTRAL_VARIABLE_HPP


//#define _DEBUG_HISTORY true


//
namespace Mistral {

  /*!
    The 'Variable' class involves:
    1/ a pointer to an implementation;
    2/ a type flag;
    3/ a list of accessors;
    
    When calling an accessor, the type flag is checked,
    the pointer is cast accordingly and called with the same accessor.
  */  

  class VarArray;
  class Variable;
  class Solver;
  class VariableImplementation {
  public:

    /// Attributes common to all variables:
    /////////////////////////////////////
    /// The linked solver
    Environment *solver;
    /// unique identifier, corresponds to its rank in the variables list of the solver
    int id;
    /////////////////////////////////////

    VariableImplementation() {solver=NULL; id=-1;}
    virtual ~VariableImplementation() {}

    bool is_initialised() const { 
      return id>-1;
    }

    virtual int get_solution_int_value() const ;//{ return solver->last_solution_lb[id] } ; 
    virtual std::string get_solution_str_value() const ;//{ return solver->last_solution_lb[id] } ; 
    virtual int get_solution_min() const ;//{ return solver->last_solution_lb[id] } ; 
    virtual int get_solution_max() const ;//{ return solver->last_solution_ub[id] } ; 

    //void trigger_value_event_and_save(Variable *x);
    void trigger_value_event_and_save();
    void trigger_event_and_save(const Event evt);
  };


  class BitsetDomain {
  
  public:
  
    BitsetDomain() {min = NOVAL; max = NOVAL; size = 0;}
    BitsetDomain(const int lb, const int ub);
    BitsetDomain(Vector< int >& vals);
    virtual ~BitsetDomain() {}
    void initialise(const int lb, const int ub, const bool vals=true);
    void initialise(Vector< int >& vals);
  
    int min;
    int max;
    int size;
    BitSet values;
  
    virtual std::ostream& display(std::ostream& os) const {
      
      //os << "<" << min << "|" << values << "|" << max << "> " ;

      if(min == max-1) {
	os << "[" << min << "," << max << "]";
      } else if(values.table && (max-min) >= size
		) {
	os << values ; //<< " [" << min << ".." << max << "] (" << size << ")";
      } else if(min==max){
	os << min;
      } else {
	os << "[" << min << ".." << max << "]";
      }
      return os;
    }
  
  };

  std::ostream& operator<< (std::ostream& os, const BitsetDomain& x);
  std::ostream& operator<< (std::ostream& os, const BitsetDomain* x);


  template < class WORD_TYPE >
  class VariableBitset : public VariableImplementation
  {
    
  public:

    /// Attributes specific to bitset variables:
    /////////////////////////////////////
    /// trail, used for backtracks
    Vector<int> trail_;
    
    /// Includes a bitset of values, a lower and an upper bound
    BitsetDomain domain;

    /// trail for the bitset representation
    WORD_TYPE **delta_;
    int **level_;
    WORD_TYPE **delta_abs;
    int **level_abs;

    /*!@name Constructors*/
    //@{
    VariableBitset() {
      id = -1;
      solver = NULL;
    };

    VariableBitset(const int lb, const int ub) {
      id = -1;
      solver = NULL;
      initialise(lb, ub);
    };

    virtual void initialise(const int lb, const int ub) {
      domain.initialise(lb, ub);
      initialise_trail();
    }

    VariableBitset(Vector< int >& values) {
      id = -1;
      solver = NULL;
      initialise(values);
    };

    virtual void initialise(Vector< int >& values) {
      domain.initialise(values);
      initialise_trail();
    }

    std::string get_history() {
      std::ostringstream buf;
      int k = trail_.size-1, n=domain.values.pos_words-domain.values.neg_words;
      BitSet dom(trail_[0], trail_[1], BitSet::empt);
      int counter[n];
      for(int i=0; i<n; ++i) counter[i] = 0;

      while(k>0) {
	buf << " " << trail_[k] << ":";
	if(trail_[k-1] == 1) buf << trail_[k-2];
	else if(trail_[k-1] == trail_[k-2]-trail_[k-3]+1) buf << "[" << trail_[k-3] << "," << trail_[k-2] << "]";
	//else {
	  int lvl = trail_[k], j=domain.values.neg_words;
	  dom.clear();
	  for(int i=domain.values.pos_words-1; i>=j; --i) {
	    dom.table[i] = *(delta_[i]-counter[i-j]);
	    if(*(level_[i]-counter[i-j]) == lvl) ++counter[i-j];
	  }
	  buf << dom ;
	  //}
	k -= 4;
      }
      return buf.str();
    }

    void initialise_trail() {
      int i, k;

      trail_.initialise(0,8);
      trail_.add(domain.min);
      trail_.add(domain.max);
      trail_.add(domain.size);
      trail_.add(-1);

      delta_abs = new WORD_TYPE*[domain.values.pos_words-domain.values.neg_words];
      delta_abs -= domain.values.neg_words;
      delta_ = new WORD_TYPE*[domain.values.pos_words-domain.values.neg_words];
      delta_ -= domain.values.neg_words;
      level_abs = new int*[domain.values.pos_words-domain.values.neg_words];
      level_abs -= domain.values.neg_words;
      level_ = new int*[domain.values.pos_words-domain.values.neg_words];
      level_ -= domain.values.neg_words;
      for(i=domain.values.neg_words; i<domain.values.pos_words; ++i) {
	k = domain.values.size(i)+1;
	delta_[i] = new WORD_TYPE[k];
	delta_abs[i] = delta_[i];
	delta_[i][0] = domain.values.table[i];
	level_[i] = new int[k];
	level_abs[i] = level_[i];
	level_[i][0] = -1;
      }
    }

    virtual ~VariableBitset() {

      for(int i=domain.values.neg_words; i<domain.values.pos_words; ++i) {
	delete [] delta_abs[i];
	delete [] level_abs[i];
      }
      delta_abs += domain.values.neg_words;
      level_abs += domain.values.neg_words;
      delta_ += domain.values.neg_words;
      level_ += domain.values.neg_words;

      delete [] delta_abs;
      delta_abs = NULL;
      delete [] level_abs;
      level_abs = NULL;
      delete [] delta_;
      delta_ = NULL;
      delete [] level_;
      level_ = NULL;

      if(domain.values.table) {
	domain.values.table += domain.values.neg_words;
	delete [] domain.values.table;
	domain.values.table = NULL;
      } //else {
      domain.values.neg_words = 0;
      //}
    }

    void set_bound_history(const int lb, const int ub, const int level) {

#ifdef _DEBUG_HISTORY
      std::cout << domain << std::endl;
      std::cout << "set bound history " << lb << ".." << ub << " @ " << level << std::endl;
#endif

      // first find out which words have changed
      int prev_lb = domain.min; //trail_.back(-4);
      int prev_ub = domain.max; //trail_.back(-3);
      int i, j, k, l;
      WORD_TYPE w;

#ifdef _DEBUG_HISTORY
      std::cout << prev_lb << " -> " << lb << ".." << ub << " <- " << prev_ub << std::endl;
#endif

      if(prev_lb < lb || prev_ub > ub) {
	j = (prev_lb >> BitSet::EXP); //BitSet::word_index(prev_lb);
	i = (lb >> BitSet::EXP);
	l = (ub >> BitSet::EXP);
	k = (prev_ub >> BitSet::EXP); //BitSet::word_index(prev_ub);
	
#ifdef _DEBUG_HISTORY
	std::cout << "change table[" << j << ":" << i << "]" << std::endl;
	std::cout << "change table[" << l << ":" << k << "]" << std::endl;
#endif

	while( j < i ) {

#ifdef _DEBUG_HISTORY
	  std::cout << "set table[" << (j) << "] = ";
	  print_bitset(domain.values.table[j], (j), std::cout);
	  std::cout << " to void" << std::endl ;
#endif

	  *(++delta_[j]) = BitSet::empt;
	  domain.values.table[j] = BitSet::empt;
	  *(++level_[j]) = level;

	  ++j;
	}
	
	while( k > l ) {

#ifdef _DEBUG_HISTORY
	  std::cout << "set table[" << k << "] = ";
	  print_bitset(domain.values.table[k], k, std::cout);
	  std::cout << " to void" << std::endl ;
#endif

	  *(++delta_[k]) = BitSet::empt;
	  domain.values.table[k] = BitSet::empt;
	  *(++level_[k]) = level;

	  --k;
	}

	if(i == l) {
	  w = (BitSet::full << (lb & BitSet::CACHE)) & (BitSet::full >> (BitSet::CACHE - (ub & BitSet::CACHE)));

#ifdef _DEBUG_HISTORY
	  std::cout << "change " ;
	  print_bitset(domain.values.table[i], i, std::cout);
	  std::cout << " to " ;
	  print_bitset(w, i, std::cout);
	  std::cout << std::endl;
#endif

	  domain.values.table[i] = w;
	  *(++delta_[i]) = w;
	  *(++level_[i]) = level;
	} else {
	  if(prev_lb < lb) {
	    w = (BitSet::full << (lb & BitSet::CACHE));

#ifdef _DEBUG_HISTORY
	  std::cout << "change " ;
	  print_bitset(domain.values.table[i], i, std::cout);
	  std::cout << " to " ;
	  print_bitset(w, i, std::cout);
	  std::cout << std::endl;
#endif

	    domain.values.table[i] = w;
	    *(++delta_[i]) = w;
	    *(++level_[i]) = level;
	  }
	  if(prev_ub > ub) {
	    w = (BitSet::full >> (BitSet::CACHE - (ub & BitSet::CACHE)));

#ifdef _DEBUG_HISTORY
	    std::cout << "change "  ;
	  print_bitset(domain.values.table[l], l, std::cout);
	  std::cout << " to " ;
	  print_bitset(w, l, std::cout);
	  std::cout << std::endl;
#endif

	    domain.values.table[l] = w;
	    *(++delta_[l]) = w;
	    *(++level_[l]) = level;
	  }
	}
      }


      domain.min = lb;
      domain.max = ub;

      trail_.add(lb);
      trail_.add(ub);
      trail_.add(ub-lb+1);
      trail_.add(level);

#ifdef _DEBUG_HISTORY
      std::cout << domain << " " << domain.values << std::endl << std::endl;
#endif

    }

    /*!@name Static Accessors and Iterators*/
    //@{
    /// Returns the assigned value if it exists
    inline int get_value() const { return domain.min; }
    /// Returns the domain size
    inline unsigned int get_size() const { return domain.size; }
    /// Returns the minimum value in the domain
    inline int get_min() const { return domain.min; }
    /// Returns the maximum value in the domain
    inline int get_max() const { return domain.max; }
    /// Returns the minimum value that could belong to the domain
    inline int get_initial_min() const { return trail_[0]; }
    /// Returns 1 + the maximum value that could belong to the domain
    inline int get_initial_max() const { return trail_[1]; }
    /// Returns the minimum value in [1..infty] \\inter domain, min if there are none
    inline int get_min_pos() const {
      if( domain.min > 0 ) return domain.min;
      if( domain.max < 1 ) return domain.min;
      if( domain.values.fast_contain(1) ) return 1;
      return domain.values.next( 1 );
    }
    /// Returns the maximum value in [-infty..-1] \\inter domain, max if there are none
    inline int get_max_neg() const  {
      if( domain.max < 0 ) return domain.max;
      if( domain.min > -1 ) return domain.max;
      if( domain.values.fast_contain(-1) ) return -1;
      return domain.values.prev( -1 );
    } 
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int next(const int v) const { return (((domain.size != domain.max - domain.min + 1)) ? domain.values.next(v) : (v < domain.max ? v+1 : v)); }
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int prev(const int v) const { return (((domain.size != domain.max - domain.min + 1)) ? domain.values.prev(v) : (v > domain.min ? v-1 : v)); }


    /// Whether or not the Variable is currently an interval
    inline bool is_range() const { return (domain.size == domain.max - domain.min + 1); }
    /// Whether or not the Variable is bound to a ground value
    inline bool is_ground() const { return (domain.size == 1); }
    /// Whether or not the Variable is bound to a given ground value
    inline bool equal(const int v) const { return (domain.min == v && domain.max == v); }
    /// Whether the value "v" is currently contained in the domain
    inline bool contain(const int v) const { return (domain.min <= v && domain.max >= v && domain.values.fast_contain(v)); }

    /// Whether the domain has a nonempty intersection with the interval [l..u]
    inline bool intersect(const int lo, const int up) const { return (domain.min <= up && domain.max >= lo && domain.values.intersect(lo, up)); }
    /// Whether the domain is included in the interval [l..u]
    inline bool included(const int lo, const int up) const { return (domain.min >= lo && domain.max <= up); }
    /// Whether the domain is included in the interval [l..u]
    inline bool includes(const int lo, const int up) const { return (domain.min <= lo && domain.max >= up && domain.values.includes(lo, up)); }

    /// Whether the domain has a nonempty intersection with the set s 
    inline bool intersect(const BitSet& s) const { return domain.values.intersect(s); }
    /// Whether the domain is included in the set s 
    inline bool included(const BitSet& s) const { return domain.values.included(s); }
    /// Whether the domain is included in the set s 
    inline bool includes(const BitSet& s) const { return domain.values.includes(s); }

    /// Whether the domain has a nonempty intersection with the Variable x
    inline bool intersect(const VariableBitset* x) const { return x->intersect(domain.values); }
    /// Whether the domain is included in the Variable x 
    inline bool included(const VariableBitset* x) const { return x->includes(domain.values); }
    /// Whether the domain is included in the Variable x 
    inline bool includes(const VariableBitset* x) const { return x->included(domain.values); }

    //     /// Whether the domain has a nonempty intersection with the Variable x
    //     bool intersect(Variable x) const ;
    //     /// Whether the domain is included in the Variable x 
    //     bool included(Variable x) const ;
    //     /// Whether the domain is included in the Variable x 
    //     bool includes(Variable x) const ;

    /// Intersect its domain with a set s
    inline void intersect_to( BitSet& s ) const { s.intersect_with(domain.values); }
    /// Do the union of its domain with a set s
    inline void union_to( BitSet& s ) const { s.union_with(domain.values); }
    //@}

    /*!@name Domain handling methods*/
    //@{
    


    inline Event remove(const int v) {
      Event removal = DOMAIN_EVENT;

      // first check if we can abort early
      if(!contain(v)) return NO_EVENT;
      if(domain.size == 1) return FAIL_EVENT;

      // if(id==1) {
      // 	std::cout << "remove " << v << " from " << this << " in " << domain << std::endl;
      // }


      save();

      // then change the static domain
      if(--domain.size == 1) removal |= VALUE_C; 
      domain.values.fast_remove(v);
	
      if(removal & VALUE_C) {
	if(domain.min == v) {
	  removal |= LB_EVENT;
	  domain.min = domain.max;
	} else {
	  removal |= UB_EVENT;
	  domain.max = domain.min;
	}
      } else {
	if(domain.max == v) {
	  removal |= UB_EVENT;
	  domain.max = domain.values.prev(v);
	} else if(domain.min == v) {
	  removal |= LB_EVENT;
	  domain.min = domain.values.next(v);
	}       
      }  
      

      // if(id==1) {
      // 	std::cout << event2str(removal) << " event ==> " << this << " = " << domain << std::endl;
      // }

      solver->trigger_event(id, removal);
      return removal; 
    }

    inline Event remove_wo_trigger(const int v) {
      Event removal = DOMAIN_EVENT;

      // first check if we can abort early
      if(!contain(v)) return NO_EVENT;
      if(domain.size == 1) return FAIL_EVENT;

      save();

      // then change the static domain
      if(--domain.size == 1) removal = VALUE_EVENT; 
      domain.values.fast_remove(v);
	
      if(removal == VALUE_EVENT)
	if(domain.min == v) domain.min = domain.max;
	else domain.max = domain.min;
      else {
	if(domain.max == v) {
	  removal |= UB_EVENT;
	  domain.max = domain.values.max();
	} else if(domain.min == v) {
	  removal |= LB_EVENT;
	  domain.min = domain.values.min();
	}       
      }  
      
      //solver->trigger_event(id, removal);
      return removal; 
    }

    /// Remove all values but "v"
    inline Event set_domain(const int v) {
      Event setdomain = VALUE_C;

      // this->display(std::cout);
      // std::cout << " <- " << v << " (was: " << domain << " " << domain.values << ") " 
      //  		<< contain(v) << "/" << domain.size << std::endl;

      // first check if we can abort early
      if(!contain(v)) return FAIL_EVENT;
      else if(domain.size == 1) return NO_EVENT;

      //setdomain = VALUE_EVENT;
      save();


      // std::cout << domain.values.table << std::endl;

      if(domain.values.table) 
	domain.values.set_to(v);
      domain.size = 1;
      if(domain.min != v) {
	domain.min = v;
	setdomain |= LB_EVENT;
      }
      if(domain.max != v) {
	domain.max = v;
	setdomain |= UB_EVENT;
      }
      
      solver->trigger_event(id, setdomain);	
      //trigger_event(setdomain);

      return setdomain; 
    }

    /// Remove all values strictly lower than l
    inline Event set_min(const int lo) {
      Event lower_bound = LB_EVENT;

      // first check if we can abort early
      if(domain.max <  lo) return FAIL_EVENT;
      if(domain.min >= lo) return NO_EVENT;
      
      // if(id==1) {
      // 	std::cout << "set_min " << lo << " of " << this << " in " << domain << std::endl;
      // }



      save();
      
      // then change the static domain
      domain.values.set_min(lo);
      if(lo == domain.max) {
	domain.min = lo;
	domain.size = 1;
	lower_bound |= VALUE_C;
      } else {
	domain.size = domain.values.size();
	if(domain.values.contain(lo)) domain.min = lo;
	else domain.min = domain.values.next(lo-1);
	if(domain.size == 1) lower_bound |= VALUE_C;
      }


      // if(id==1) {
      // 	std::cout << event2str(lower_bound) << " event ==> " << this << " = " << domain << std::endl;
      // }
      
      solver->trigger_event(id, lower_bound);
      //trigger_event(lower_bound);
      return lower_bound; 
    }

    /// Remove all values strictly greater than u
    inline Event set_max(const int up) {
      Event upper_bound = UB_EVENT;

      //std::cout << "x_set_max(" << up << ") " << domain << " " << domain.values << std::endl;



      // first check if we can abort early
      if(domain.min >  up) return FAIL_EVENT;
      if(domain.max <= up) return NO_EVENT;
      
      save();

      // then change the static domain
      domain.values.set_max(up);
      if(up == domain.min) {
	domain.max = up;
	domain.size = 1;
	upper_bound |= VALUE_C;
      } else {
	domain.size = domain.values.size();
	if(domain.values.contain(up)) domain.max = up;
	else domain.max = domain.values.prev(up);
	if(domain.size == 1) upper_bound |= VALUE_C;
      }

      //std::cout << domain << std::endl;

      solver->trigger_event(id, upper_bound);
      //trigger_event(upper_bound);
      return upper_bound; 
    }


    /// Remove all values that do not appear in the set "s"
    inline Event set_domain(const BitSet& s) {
      Event intersection = DOMAIN_EVENT;

      if( !domain.values.intersect(s) ) return FAIL_EVENT;
      if( domain.values.included(s) ) return NO_EVENT;

      save();

      // then change the static domain
      domain.values.intersect_with(s);
      domain.size = domain.values.size();
      if(!s.contain(domain.min)) { intersection |= LB_EVENT; domain.min = domain.values.next(domain.min); }
      if(!s.contain(domain.max)) { intersection |= UB_EVENT; domain.max = domain.values.prev(domain.max); }
      if(domain.min == domain.max) intersection |= VALUE_EVENT;
    
      solver->trigger_event(id, intersection);
      //trigger_event(intersection);
      return intersection; 
    }


    //   Event set_domain(Variable x) ;
    //     /// Remove all values that do not appear in the current domain of the Variable "x"
    //     inline Event set_domain(Variable x) {
    //       if(domain.size == 1) {
    // 	return(x.contain(domain.min) ? NO_EVENT : FAIL_EVENT);
    //       } else {
    // 	int xsize = x.get_size();
    // 	if(xsize == 1) {
    // 	  return set_domain(x.get_value());
    // 	} else {
    // 	  int xmin = x.get_min();
    // 	  int xmax = x.get_max();
    // 	  if(xsize == (xmax-xmin+1)) {
    // 	    return(set_min(xmin) | set_max(xmax));
    // 	  } else if(x.includes(domain.values)) {
    // 	    return NO_EVENT;
    // 	  } else if(x.intersect(domain.values)) {
    // 	    Event intersection = DOMAIN_EVENT;
    // 	    save();
    // 	    x.intersect_to(domain.values);
    // 	    domain.size = domain.values.size();
    // 	    if(!domain.values.contain(domain.min)) { intersection |= LB_EVENT; domain.min = domain.values.next(domain.min); }
    // 	    if(!domain.values.contain(domain.max)) { intersection |= UB_EVENT; domain.max = domain.values.prev(domain.max); }
    // 	    if(domain.min == domain.max) intersection |= VALUE_EVENT;
    
    // 	    solver->trigger_event(id, intersection);
	    
    // 	  } else return FAIL_EVENT;
    // 	}
    //       }
    //       return NO_EVENT;
    //     }

    /// Remove all values that belong to the set "s"
    inline Event removeSet(const BitSet& s) {
      Event setdifference = NO_EVENT;

      //       if(type != VIRTUAL_VAR) {
      // 	// first check if we can abort early
      // 	if(domain.values.table) {
      // 	  if( !domain.values.intersect(s) ) return NO_EVENT;
      // 	  if( domain.values.included(s) ) return FAIL_EVENT;

      // 	  _save_();
	  
      // 	  // then change the static domain
      // 	  setdifference = DOMAIN_EVENT;
      // 	  domain.values.setminus_with(s);
      // 	  domain.size = domain.values.size();
      // 	  if(s.contain(domain.min)) { setdifference |= LB_EVENT; domain.min = domain.values.next(domain.min); }
      // 	  if(s.contain(domain.max)) { setdifference |= UB_EVENT; domain.max = domain.values.prev(domain.max); }
      // 	  if(domain.min == domain.max) setdifference |= VALUE_EVENT;
      // 	} else {
      // 	  std::cerr << "not supported" << std::endl;
      // 	  exit(0);
      // 	}
      //       } 

      //       trigger_event(setdifference);
      return setdifference; // do nothing more for BOOL or CONSTANT or BITSET or RANGE
    }

    /// Remove all values in the interval [l..u]
    inline Event remove_interval(const int lo, const int up) {

      //std::cout << "REMOVE [" << lo << ".." << up << "] from " << domain << std::endl;

      if(lo <= domain.min) {
	//std::cout << "EQUIVALENT TO >= " << (up+1) << std::endl;
	return set_min(up+1);
      }
      if(up >= domain.max) {
	//std::cout << "EQUIVALENT TO <= " << (lo-1) << std::endl;
	return set_max(lo-1);
      }
      if(domain.size != domain.max-domain.min+1) {
	if(!domain.values.intersect(lo, up)) {
	  //std::cout << "NO INTERSECTION -> NO EVENT" << std::endl;
	  return NO_EVENT;
	}
	save();
	domain.values.remove_interval(lo, up);
	domain.size = domain.values.size();
      } else {
	save();
	domain.values.remove_interval(lo, up);
	domain.size -= (up-lo+1);
      }

      solver->trigger_event(id, DOMAIN_EVENT);
      //trigger_event(DOMAIN_EVENT);
      return DOMAIN_EVENT; 
    }
    //@}


    /*!@name Virtual Accessors and Iterators*/
    //@{

    //@}   

    inline void save() {

      // std::cout << "SAVE " << this << " at level " << solver->level << std::endl
      //  		<< trail_  ;
      // // for(int i=domain.values.pos_words; i-->domain.values.neg_words;)
      // // 	std::cout << " " << delta_[i]-delta_abs[i] ;
      // std::cout << std::endl;

      int i, j;//, lvl = get_current_level();
      if(trail_.back() != solver->level) {

	//if(trail_.back() != lvl) {
	//store();
	//Variable self(this, BITSET_VAR);
	//solver->save(self);
	//solver->save(this, BITSET_VAR);
	solver->save(id);
	//store(this, BITSET_VAR);

	//int i = domain.values.pos_words, j = domain.values.neg_words;
	i = domain.values.pos_words;
	j = domain.values.neg_words;

	trail_.add(domain.min);
	trail_.add(domain.max);
	trail_.add(domain.size);
	
	if(delta_) 
	  while( i --> j ) {
	    //WORD_TYPE buf = domain.values.table[i];
	    if(*(delta_[i]) != domain.values.table[i]) {
	      *(++delta_[i]) = domain.values.table[i];
	      *(++level_[i]) = solver->level; //lvl;
	    }
	  }
	
	trail_.add(solver->level);
      }


    // std::cout << "===> " << this << std::endl
    //    		<< trail_  ;
    //   std::cout << std::endl;
      
    }

    inline Event restore() {

      // backtrack from solver->level to solver->level-1;
      // the domain has changed at solver->level
      int i = domain.values.pos_words;
      int j = domain.values.neg_words;
      //int lvl = get_current_level();

      trail_.pop();
      trail_.pop(domain.size);
      trail_.pop(domain.max);
      trail_.pop(domain.min);

      if(delta_) 
	while( i --> j ) {
	  domain.values.table[i] = *(delta_[i]);
	  if(*(level_[i]) == solver->level) {
	    //if(*(level_[i]) == lvl) {
	    --level_[i];
	    --delta_[i];
	  }
	}

      return NO_EVENT;
    }

    virtual std::ostream& display(std::ostream& os) const
    {
      os << "x" << id;
      return os;
    }
    
    virtual void debug_print() const
    {

      //       //Variable self((VariableImplementation*)this, BITSET_VAR);

      std::cout << "x" << id << " in " 
		<< domain << " "
     		<< trail_ << " " << std::endl;

      for(int i=domain.values.neg_words; i<domain.values.pos_words; ++i) {
     	std::cout << "word[" << i << "] " ;
     	//for(WORD_TYPE* it = delta_abs[i]; it <= delta_[i]; ++it) {
	for(int j=0; j<=(delta_[i]-delta_abs[i]); ++j) {
	  std::cout << level_abs[i][j] << ":";
     	  print_bitset(delta_abs[i][j], i, std::cout);
     	}
     	std::cout << std::endl;
      }
    }

  };


  template < class WORD_TYPE >
  std::ostream& operator<< (std::ostream& os, const VariableBitset< WORD_TYPE >& x) {
    return x.display(os);
  }

  template < class WORD_TYPE >
  std::ostream& operator<< (std::ostream& os, const VariableBitset< WORD_TYPE > * x) {
    return x->display(os);
  }



  typedef VariableBitset< unsigned long long int > VariableBitset64;
  typedef VariableBitset< unsigned int > VariableBitset32;

#ifdef _BIT64
  
  typedef VariableBitset64 VariableBitmap;

#else

  typedef VariableBitset32 VariableBitmap;

#endif


  template < class WORD_TYPE, int NUM_WORDS >
  class VariableWord : public VariableBitset< WORD_TYPE > {

  public: 

    WORD_TYPE word[NUM_WORDS];

    /*!@name Constructors*/
    //@{
    VariableWord(const int lb, const int ub) 
    //: VariableBitset< WORD_TYPE >(lb, ub) 
    {
      //VariableBitset< WORD_TYPE >::id = -1;
      //VariableBitset< WORD_TYPE >::solver = NULL;
      initialise(lb, ub);
    }

    virtual void initialise(const int lb, const int ub) {



      int nw = (lb >> BitSet::EXP);
      int pw = (ub >> BitSet::EXP)+1;


      // std::cout << " initialise " << lb << " " << ub 
      // 		<< " -> " << nw << " " << pw 
      // 		<< std::endl;

      //assert((pw-nw) <= NUM_WORDS);
      VariableBitset< WORD_TYPE >::domain.initialise(lb, ub, false);
      VariableBitset< WORD_TYPE >::domain.values.neg_words = nw;
      VariableBitset< WORD_TYPE >::domain.values.pos_words = pw;
      VariableBitset< WORD_TYPE >::domain.values.table = ((&(word[0]))-nw);

      std::fill(VariableBitset< WORD_TYPE >::domain.values.table+nw, 
		VariableBitset< WORD_TYPE >::domain.values.table+pw,
		0);

      VariableBitset< WORD_TYPE >::domain.values.add_interval(lb,ub);
      VariableBitset< WORD_TYPE >::initialise_trail();
      
      //  this->VariableBitset< WORD_TYPE >::display(std::cout);
      //  std::cout  << " in " << VariableBitset< WORD_TYPE >::domain << std::endl;
      // std::cout  << " in " << VariableBitset< WORD_TYPE >::domain.values << std::endl;
      //  std::cout  << " in " << VariableBitset< WORD_TYPE >::domain.values.table[0] << std::endl;
    }


    VariableWord(Vector< int >& values) 
    {
      initialise(values);
    }

    virtual void initialise(Vector< int >& values) {

      int nw = (values.front() >> BitSet::EXP);
      int pw = (values.back() >> BitSet::EXP)+1;


      // std::cout << " initialise " << lb << " " << ub 
      // 		<< " -> " << nw << " " << pw 
      // 		<< std::endl;

      //assert((pw-nw) <= NUM_WORDS);
      VariableBitset< WORD_TYPE >::domain.initialise(values.front(), values.back(), false);
      VariableBitset< WORD_TYPE >::domain.values.neg_words = nw;
      VariableBitset< WORD_TYPE >::domain.values.pos_words = pw;
      VariableBitset< WORD_TYPE >::domain.values.table = ((&(word[0]))-nw);

      std::fill(VariableBitset< WORD_TYPE >::domain.values.table+nw, 
		VariableBitset< WORD_TYPE >::domain.values.table+pw,
		0);

      //VariableBitset< WORD_TYPE >::domain.values.add_interval(lb,ub);

      for(unsigned int i=0; i<values.size; ++i) 
	VariableBitset< WORD_TYPE >::domain.values.add(values[i]);

      VariableBitset< WORD_TYPE >::initialise_trail();
      
      //  this->VariableBitset< WORD_TYPE >::display(std::cout);
      //  std::cout  << " in " << VariableBitset< WORD_TYPE >::domain << std::endl;
      // std::cout  << " in " << VariableBitset< WORD_TYPE >::domain.values << std::endl;
      //  std::cout  << " in " << VariableBitset< WORD_TYPE >::domain.values.table[0] << std::endl;
    }

    virtual ~VariableWord() {

      //std::cout << VariableBitset< WORD_TYPE >::id << " delete subclass" << std::endl;
      
      VariableBitset< WORD_TYPE >::domain.values.table = NULL;
      //VariableBitset< WORD_TYPE >::domain.values.neg_words = 0;
    }

  };



  class VariableList : public VariableBitmap {};
  class VariableRange : public VariableImplementation {

  public:

    int min;
    int max;

    Vector<int> trail_;

    /*!@name Constructors*/
    //@{
    VariableRange(const int lb, const int ub) : VariableImplementation() {
      initialise(lb, ub);
    };

    virtual void initialise(const int lb, const int ub) {
      min = lb;
      max = ub;
      trail_.add(min);
      trail_.add(max);
      trail_.add(-1);
    }

    virtual ~VariableRange() {
    }


   std::string get_history() {
      std::ostringstream buf;
      int k = trail_.size-1;
      while(k>0) {
	buf << " " << trail_[k] << ":";
	if(trail_[k-1] == trail_[k-2]) buf << trail_[k-2];
	else buf << "[" << trail_[k-2] << "," << trail_[k-1] << "]";
	k -= 3;
      }
      return buf.str();
    }

    // set the history of X to match self
    void set_history(VariableBitmap *X) {

#ifdef _DEBUG_HISTORY
      std::cout << "setup a bitset trail for ";
      display(std::cout);
      std::cout << " in [" << min << "," << max << "] at level " << solver->level << std::endl;
      std::cout << "current trail: " << trail_ << std::endl;
#endif

      for(unsigned int i = 3; i<trail_.size; i+=3) {
	//if(trail_[i] != trail_[i-3] || trail_[i+1] != trail_[i-2])
	X->set_bound_history(trail_[i], trail_[i+1], trail_[i+2]);
      }

      // if(min > X->domain.min || max < X->domain.max) {
      // 	X->set_bound_history(min, max, solver->level);
      // }

      //X->set_min(min);
      //X->set_max(max);

      X->domain.values.set_min(min);
      X->domain.values.set_max(max);
      X->domain.min = min;
      X->domain.max = max;
      X->domain.size = (max-min+1);

      // if(id == 35) {

#ifdef _DEBUG_HISTORY
      std::cout << "new trail: " << X->get_history() << std::endl
		<< "domain: " << X->domain << std::endl;
#endif
      // 	exit(1);
      // }
    }

    // void set_bound_history(const int lb, const int ub, const int level) {
    //   // first find out which words have changed
    //   min = lb;
    //   max = ub;
      
    //   trail_.add(lb);
    //   trail_.add(ub);
    //   trail_.add(level);
    // }


    /*!@name Static Accessors and Iterators*/
    //@{
    /// Returns the assigned value if it exists
    inline int get_value() const { return min; }
    /// Returns the domain size
    inline unsigned int get_size() const { return max-min+1; }
    /// Returns the minimum value in the domain
    inline int get_min() const { return min; }
    /// Returns the maximum value in the domain
    inline int get_max() const { return max; }
    /// Returns the minimum value that could belong to the domain
    inline int get_initial_min() const { return trail_[0]; }
    /// Returns 1 + the maximum value that could belong to the domain
    inline int get_initial_max() const { return trail_[1]; }
    /// Returns the minimum value in [1..infty] \\inter domain, min if there are none
    inline int get_min_pos() const {
      if( min > 0 ) return min;
      if( max < 1 ) return min;
      return 1;
    }
    /// Returns the maximum value in [-infty..-1] \\inter domain, max if there are none
    inline int get_max_neg() const  {
      if( max < 0 ) return max;
      if( min > -1 ) return max;
      return -1;
    } 
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int next(const int v) const { return (v < max ? (v>=min ? v+1 : min) : v); }
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int prev(const int v) const { return (v > min ? (v<=max ? v-1 : max) : v); }


    /// Whether or not the Variable is currently an interval
    inline bool is_range() const { return true; }
    /// Whether or not the Variable is bound to a ground value
    inline bool is_ground() const { return min == max; }
    /// Whether or not the Variable is bound to a given ground value
    inline bool equal(const int v) const { return (min == v && max == v); }
    /// Whether the value "v" is currently contained in the domain
    inline bool contain(const int v) const { return (min <= v && max >= v); }

    /// Whether the domain has a nonempty intersection with the interval [l..u]
    inline bool intersect(const int lo, const int up) const { return (min <= up && max >= lo); }
    /// Whether the domain is included in the interval [l..u]
    inline bool included(const int lo, const int up) const { return (min >= lo && max <= up); }
    /// Whether the domain is included in the interval [l..u]
    inline bool includes(const int lo, const int up) const { return (min <= lo && max >= up); }

    /// Whether the domain has a nonempty intersection with the set s 
    //inline bool intersect(const BitSet& s) const { return (min <= s.max() && max >= s.min()); }
    inline bool intersect(const BitSet& s) const { return s.intersect(min, max); }
    /// Whether the domain is included in the set s 
    inline bool included(const BitSet& s) const { return s.includes(min, max); }
    /// Whether the domain is included in the set s 
    inline bool includes(const BitSet& s) const { return s.included(min, max); }

    //    /// Whether the domain has a nonempty intersection with the Variable x
    //     inline bool intersect(const Variable x) const { return x.intersect(min,max); }
    //     /// Whether the domain is included in the Variable x 
    //     inline bool included(const Variable x) const { return x.includes(min,max); }
    //     /// Whether the domain is included in the Variable x 
    //     inline bool includes(const Variable x) const { return x.included(min,max); }

    //     /// Whether the domain has a nonempty intersection with the Variable x
    //     bool intersect(const Variable x) const ;
    //     /// Whether the domain is included in the Variable x 
    //     bool included(const Variable x) const ;
    //     /// Whether the domain is included in the Variable x 
    //     bool includes(const Variable x) const ;

    /// Intersect its domain with a set s
    inline void intersect_to( BitSet& s ) const { s.set_min(min); s.set_max(max); }
    /// Do the union of its domain with a set s
    inline void union_to( BitSet& s ) const { // s.fill(min,max);
      s.fill(min,max);
    }
    //@}

    /*!@name Domain handling methods*/
    //@{
    


     Event remove(const int v);//  {
    //   Event removal = DOMAIN_EVENT;

    //   // first check if we can abort early
    //   if(min!=v && max!=v) {
    // 	solver->make_non_convex(id);
    // 	return solver->variables[id].remove(v);
    // 	//return NO_EVENT;
    //   }
    //   if(min==max) return FAIL_EVENT;

    //   save();

    //   if(min==v) {
    // 	++min;
    // 	removal |= LB_EVENT;
    //   } else {
    // 	--max;
    // 	removal |= UB_EVENT;
    //   }

    //   if(min == max) removal |= VALUE_EVENT; 

    //   //       display(std::cout);
    //   //       std::cout << " in [" << min << ".." << max << "]: ";
    //   //       std::cout << "TRIGGER " << removal << std::endl;
    //   //       std::cout << (ASSIGNED(removal)) << std::endl;

    //   solver->trigger_event(id, removal);
    //   return removal; 
    // }

    /// Remove all values but "v"
    inline Event set_domain(const int v) {
      Event setdomain = VALUE_EVENT;

      // first check if we can abort early
      if(!contain(v)) return FAIL_EVENT;
      if(min == max) return NO_EVENT;

      save();

      if(v==min) setdomain^=LB_C;
      if(v==max) setdomain^=UB_C;

      min = max = v;
      
      solver->trigger_event(id, setdomain);	

      return setdomain; 
    }

    /// Remove all values strictly lower than l
    inline Event set_min(const int lo) {
      Event lower_bound = LB_EVENT;

      // first check if we can abort early
      if(max <  lo) return FAIL_EVENT;
      if(min >= lo) return NO_EVENT;
      
      save();

      min = lo;
      if(min == max) lower_bound |= VALUE_EVENT;

      //      display(std::cout);
      //       std::cout << " in [" << min << ".." << max << "]: ";
      //       std::cout << "TRIGGER " << lower_bound << std::endl;
      //       std::cout << (ASSIGNED(lower_bound)) << std::endl;
      
      solver->trigger_event(id, lower_bound);
      return lower_bound; 
    }

    /// Remove all values strictly greater than u
    inline Event set_max(const int up) {

      //std::cout << "set_max(" << up << ")" << std::endl;

      Event upper_bound = UB_EVENT;

      // first check if we can abort early
      if(min >  up) return FAIL_EVENT;
      if(max <= up) return NO_EVENT;
      
      save();

      max = up;
      if(max == min) upper_bound |= VALUE_EVENT;
      
      //       display(std::cout);
      //       std::cout  << " in [" << min << ".." << max << "]: ";
      //       std::cout << "TRIGGER " << upper_bound << std::endl;
      //       std::cout << (ASSIGNED(upper_bound)) << std::endl;

      solver->trigger_event(id, upper_bound);
      return upper_bound; 
    }


    /// Remove all values that do not appear in the set "s"
     Event set_domain(const BitSet& s);//  {
    //   //return set_domain(s.next(min-1), s.prev(max+1));
    //   return NO_EVENT;
    // }


    Event set_domain(Variable x) ;
    /// Remove all values that do not appear in the current domain of the Variable "x"
    //     inline Event set_domain(Variable x) {
    //       //return set_domain(x.next(min-1), x.prev(max+1));
    //       return NO_EVENT;
    //     }

    /// Remove all values that belong to the set "s"
    inline Event removeSet(const BitSet& s) {
      Event setdifference = NO_EVENT;
      return setdifference; // do nothing more for BOOL or CONSTANT or BITSET or RANGE
    }

    /// Remove all values in the interval [l..u]
    inline Event remove_interval(const int lo, const int up) {
      if(lo <= min) return set_min(up+1);
      if(up >= max) return set_max(lo-1);
      return NO_EVENT; 
    }
    //@}


    /*!@name Virtual Accessors and Iterators*/
    //@{

    //@}   

    inline void save() {
      if(trail_.back() != solver->level) {
	//Variable self(this, RANGE_VAR);
	//solver->save(self);
	//solver->save(this, RANGE_VAR);
	solver->save(id);
	trail_.add(min);
	trail_.add(max);
	trail_.add(solver->level);
      }
    }

    inline Event restore() {
      trail_.pop();
      trail_.pop(max);
      trail_.pop(min);
      return NO_EVENT;
    }

    virtual std::ostream& display(std::ostream& os) const
    {
      os << "x" << id;
      return os;
    }


  };
  class VariableVirtual : public VariableBitmap {};


  //   class VariableStruct {
  //     union {
  //       int domain_type;
  //       int* bool_domain;
  //     };
  //     union {
  //       int constant_value;
  //       VariableImplementation* variable;
  //       VariableVirtual* virtual_domain;
  //       VariableBitmap* bitset_domain;
  //       VariableRange* range_domain;
  //       VariableList* list_domain;
  //       //Expression* expression;
  //       Constraint* constraint;
  //     };
  //   };

  //   class VariableBitmap;
  //   class VariableRange;
  //   class VariableList;


  class Expression;


  class Variable  {

  public:

    union {
      int domain_type;
      int* bool_domain;
    };
    union {
      int constant_value;
      VariableImplementation* variable;
      VariableVirtual* virtual_domain;
      VariableBitmap* bitset_domain;
      VariableRange* range_domain;
      VariableList* list_domain;
      Expression* expression;
      ConstraintImplementation* constraint;
    };

    Variable();
    Variable(const int value);
    Variable(VariableImplementation* impl, const int type=DYN_VAR);
    Variable(Expression* impl);
    //     Variable(const int lo, const int up, const int type=DYN_VAR);
    //     Variable(Vector< int >& values, const int type=DYN_VAR);
    Variable(const int lo, const int up, const int type=EXPRESSION);
    Variable(Vector< int >& values, const int type=EXPRESSION);
    Variable(Variable X, bool h);
    //Variable(const int lo, const int up, BitSet& values, const int type=EXPRESSION);


    ~Variable() { } //std::cout << "del "; this->display(std::cout); std::cout << std::endl; }

    // class iterator {
      
    // private:
    //   int *values;
    //   int size;
    //   int current;

    // public:
    //   iterator(Variable X) {
	
    //   }

    //   virtual ~iterator()

    // };

    bool is_void() { return (domain_type != CONST_VAR && variable == NULL); }
    void initialise_domain(const int lo, const int up, const int type);
    void initialise_domain(Vector< int >& values, const int type);
    //void set_history(Variable X);
    //void set_bound_history(const int lb, const int ub, const int level);


    //Variable get_children();
    Variable get_var();
    const Variable get_var() const;
    bool is_expression() { return domain_type == EXPRESSION; }

    Variable operator+(Variable);
    Variable operator+(const int);
    Variable operator-(Variable);
    Variable operator-(const int);
    Variable operator*(Variable);
    Variable operator*(const int);
    Variable operator/(Variable);
    //    Variable operator/(const int);
    Variable operator&&(Variable);
    Variable operator||(Variable);
    //    Variable operator^(Variable);
    //    Variable operator->(Variable);
    Variable operator==(Variable);
    Variable operator==(const int);
    Variable operator!=(Variable);
    Variable operator!=(const int);
    Variable operator<(Variable);
    Variable operator<=(Variable);
    Variable operator>(Variable);
    Variable operator>=(Variable);
    Variable operator<(const int);
    Variable operator<=(const int);
    Variable operator>(const int);
    Variable operator>=(const int);
    Variable operator-();
    Variable operator!();

    void add_to(Solver *s) {
    
      //std::cout << "before init: " << domain_type << std::endl;

      initialise(s);
    
      //std::cout << "after init: " << domain_type << std::endl;

      //Variable x = get_var();

      //std::cout << "pointed var: " << x.domain_type << std::endl;

      *this = get_var();

      //std::cout << "after copy: " << domain_type << std::endl;
    }
    void initialise(Solver *s, const bool top=true);
    //    void initialise(Solver *s);

    Event setValue( const int val );    
    Event setState( const int vals );    
    inline int id() const {return (domain_type == CONST_VAR ? -1 : variable->id);}
    //inline 
    Solver* get_solver();//  {
    //   std::cout << (int*)(variable) << std::endl;

    //   return (Solver*)(variable->solver);
    // }

    /*!@name Constant Accessors and Iterators*/
    //@{
    /// Returns the assigned value if it exists
    int get_value() const ; 
    std::string get_solution_str_value() const ; 
    int get_solution_int_value() const ; 
    int get_solution_min() const ; 
    int get_solution_max() const ; 
    /// Returns the domain 
    //Bitset_domain get_domain() const ; 
    std::string get_domain() const ; 
    std::string get_history() const ; 
    /// Returns the domain size
    unsigned int get_size() const ; 
    /// Returns the degree (number of constraints)
    unsigned int get_degree() const ; 
    /// Returns the minimum value in the domain
    int get_min() const ;
    /// Returns the maximum value in the domain
    int get_max() const ; 
    /// Returns the minimum value that could belong to the domain
    int get_initial_min() const ;
    /// Returns 1 + the maximum value that could belong to the domain
    int get_initial_max() const ;
    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    int get_min_pos() const ;
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    int get_max_neg() const ;
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    int next(const int v) const ;
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    int prev(const int v) const ;
    /// Whether or not the Variable is currently an interval
    bool is_range() const ;
    /// Whether or not the Variable is bound to a ground value
    bool is_ground() const ;
    /// Whether or not the Variable is bound to a given ground value
    bool equal(const int v) const ;
    /// Whether the value "v" is currently contained in the domain
    bool contain(const int v) const ;

    /// Whether the domain has a nonempty intersection with the interval [l..u]
    bool intersect(const int lo, const int up) const ;
    /// Whether the domain is included in the interval [l..u]
    bool included(const int lo, const int up) const ;
    /// Whether the domain is included in the interval [l..u]
    bool includes(const int lo, const int up) const ;

    /// Whether the domain has a nonempty intersection with the set s 
    bool intersect(const BitSet& s) const ;
    /// Whether the domain is included in the set s 
    bool included(const BitSet& s) const ;
    /// Whether the domain is included in the set s 
    bool includes(const BitSet& s) const ;

    /// Whether the domain has a nonempty intersection with the Variable x
    bool intersect(const Variable& x) const ;
    /// Whether the domain is included in the Variable x 
    bool included(const Variable& x) const ;
    /// Whether the domain is included in the Variable x 
    bool includes(const Variable& x) const ;

    /// Intersect its domain with a set s
    void intersect_to( BitSet& s ) const ;
    /// Do the union of its domain with a set s
    void union_to( BitSet& s ) const ;
    //@}

    /*!@name Domain changing methods*/
    //@{
    /// Remove value "v"
    Event remove(const int v) ;
    Event remove_wo_trigger(const int v) ;
    /// Remove all values but "v"
    Event set_domain(const int v) ;
    /// Remove all values strictly lower than lo
    Event set_min(const int lo) ;
    /// Remove all values strictly greater than up
    Event set_max(const int up) ;
    /// Remove all values that do not appear in the set "s"
    Event set_domain(const BitSet& s) ;
    /// Remove all values that do not appear in the current domain of the Variable "x"
    Event set_domain(Variable& x) ;
    /// Remove all values that belong to the set "s"
    Event removeSet(const BitSet& s) ;
    /// Remove all values in the interval [l..u]
    Event remove_interval(const int lo, const int up) ;
    //@}

    /*!@name Backtracking methods*/
    //@{
    /// Restore the last solved state
    Event restore();
    //@}

    /*!@name printing methods*/
    //@{
    //void full_print() const ;
    std::ostream& display(std::ostream& os) const ;
    //void debug_print() const ;
    //@}
  };

  std::ostream& operator<< (std::ostream& os, const Variable& x);
  std::ostream& operator<< (std::ostream& os, const Variable* x);


  /**********************************************
   * Decision
   **********************************************/
  /*! \class Decision
    \brief Representation of simplification decision

    A Decision is an object with a method 
    propagate() that somehow simplifies the problem
    and another method invert() that returns the 
    complementary (logically inverse) decision 
  */
  class Decision {
  public:
    // 0 -> removal 'n'
    static const int    REMOVAL = 0;
    // 1 -> assignment 'e'
    static const int ASSIGNMENT = 1;
    // 2 -> lower bound 'g'
    static const int LOWERBOUND = 2;
    // 3 -> upper bound 'l'
    static const int UPPERBOUND = 3;

    /**
       2 bits for the type 
       30 bits for the value
    */
    int   _data_;
    Variable var;

    inline int type() const {return _data_&3; }
    inline int value() const {return _data_>>2; }
    
    Decision() {
      init_data(ASSIGNMENT,NOVAL);
    }

    Decision(Variable x, const int t, const int v) {
      init_data(t,v);
      var = x;
    }
    
    void init_data(const int t, const int v) {
      _data_ = (((v - (t == LOWERBOUND)) << 2) | t);
    }

    Decision(ConstraintImplementation *c) {
      _data_ = -1;
      var = Variable();
      var.constraint = c;
    }

    virtual ~Decision() {
    }

    inline void invert() { _data_^=1; }
    
    inline bool make() {

      //std::cout << _data_ << std::endl;

      if(_data_ == -1) return propagateRelation();
      switch(type()) {
      case REMOVAL:    return !IS_FAIL(var.remove(value()));
      case ASSIGNMENT: return !IS_FAIL(var.set_domain(value()));
      case LOWERBOUND: return !IS_FAIL(var.set_min(value()+1));
      case UPPERBOUND: return !IS_FAIL(var.set_max(value()));
      }
      return true;
    }

    bool propagateRelation();

    std::ostream& display(std::ostream& os) const {
      os << var;
      switch(type()) {
      case REMOVAL:    { os << " =/= " ; } break;
      case ASSIGNMENT: { os << " == "  ; } break;
      case LOWERBOUND: { os << " > "   ; } break;
      case UPPERBOUND: { os << " <= "  ; } break;
      }
      os << value();
      return os;
    }

    bool operator==(const Decision& d) {
      return(var.id() == d.var.id() && _data_ == d._data_);
    }

  };

  std::ostream& operator<< (std::ostream& os, const Decision& x);
  std::ostream& operator<< (std::ostream& os, const Decision* x);


  /*
    Variables are created as expressions, with a pointer to a var object
    When expressions are added to the solver, the pointed vars are built
    either as range or boolean variable. Also, when posting a constraint
    expression variables might be marked as non convex. 

    At any time, expressions can be replaced by their variable-object
    in the constraints' scopes. 
  */

  // 1/ creation: A Variable object with domain type "EXPRESSION" and
  //              a pointer to an expression. 
  //              The expression itself points to a variable, 
  //               - either a Bitset or a Boolean variable 
  //              operations on Expression variables are as follows:
  //                 if(self_initialised()) {
  //                   self.operation()
  //                 } else {
  //                   self.unsaved_operation()
  //                 }
  //              
  // 2/ add to a Model: the reification variable is given a solver, 
  //                    an id, and added to the solver variables list.
  //                    Its trail is also initialised.
  //                    Moreover, if it was 

  //              These objects contain the domain information which can be read.


  // 2/ use in expressions
  // 3/ add to solver | here we need to 



  // Variable objects can be expression (domain_type = EXPRESSION)
  // In this case, their implementation is an object predicate that:
  // 1/ keeps the tree structure (children)
  // 2/ can be 'extracted' to produce 
  //     a/ an extra variable and a constraint (when nested)
  //     b/ a constraint (when at the top level)
  // when adding a variable, the solver
  // 1/ check what type it is
  //   a/ if it is a variable, it initilise it and add it to the stack
  //   b/ if it is an expression, it first recursively add the children
  //      then:
  //      i/ at top level: extract a constraint and post it
  //      ii/ otherwise: extract a variable 
  //          then extract a constraint and post it
  class Expression : public VariableImplementation {

  public :
    
    //int id;
    Variable self;
    Vector< Variable > children;
    
    
    Expression() : VariableImplementation() {  id=-1; }
    Expression(Variable X);
    Expression(Variable X, Variable Y);
    Expression(const int lo, const int up);
    Expression(Vector< Variable >& args);
    Expression(Vector< Variable >& args, const int lo, const int up);
    Expression(std::vector< Variable >& args);
    Expression(std::vector< Variable >& args, const int lo, const int up);
    Expression(Vector< int >& values);
    virtual ~Expression();
    
    void add(Variable X) { children.add(X); }

    virtual void extract_constraint(Solver*) {}
    virtual void extract_predicate(Solver*) {}
    //virtual void reify(Solver *s, Variable X);
    virtual void extract_variable(Solver*);
    virtual const char* get_name() const { return "var"; }
    bool is_set() const {
      std::string name = get_name();
      return (name[0] == 's' && name[1] == 'e' && name[2] == 't');
    }

    virtual std::ostream& display(std::ostream& os) const;    
  };






  // class BinaryExpression : public Expression {

  // public :

  //   BinaryExpression(Variable X, Variable Y);
  //   ~BinaryExpression();

  // };

  class AddExpression : public Expression {

  public:

    AddExpression(Variable X, Variable Y);
    virtual ~AddExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  class MulExpression : public Expression {

  public:

    MulExpression(Variable X, Variable Y);
    virtual ~MulExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  // class DivExpression : public Expression {

  // public:

  //   DivExpression(Variable X, Variable Y);
  //   virtual ~DivExpression();

  //   virtual void extract_constraint(Solver*);
  //   virtual void extract_variable(Solver*);
  //   virtual void extract_predicate(Solver*);
  //   virtual const char* get_name() const;

  // };

  class OffsetExpression : public Expression {

  public:

    int offset;

    OffsetExpression(Variable X, const int ofs);
    virtual ~OffsetExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  class FactorExpression : public Expression {

  public:

    int factor;

    FactorExpression(Variable X, const int ofs);
    virtual ~FactorExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  class SubExpression : public Expression {

  public:

    SubExpression(Variable X, Variable Y);
    virtual ~SubExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  class AndExpression : public Expression {

  public:

    AndExpression(Variable X, Variable Y);
    virtual ~AndExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  class OrExpression : public Expression {

  public:

    OrExpression(Variable X, Variable Y);
    virtual ~OrExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  class NotExpression : public Expression {

  public:

    NotExpression(Variable X);
    virtual ~NotExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };


  class NegExpression : public Expression {

  public:

    NegExpression(Variable X);
    virtual ~NegExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  // class NeqExpression : public Expression {

  // public:

  //   NeqExpression(Variable X, Variable Y);
  //   virtual ~NeqExpression();

  //   virtual void extract_constraint(Solver*);
  //   virtual void extract_variable(Solver*);
  //   virtual void extract_predicate(Solver*);
  //   virtual const char* get_name() const;

  // };

  // class EqualExpression : public Expression {

  // public:

  //   EqualExpression(Variable X, Variable Y);
  //   virtual ~EqualExpression();

  //   virtual void extract_constraint(Solver*);
  //   virtual void extract_variable(Solver*);
  //   virtual void extract_predicate(Solver*);
  //   virtual const char* get_name() const;

  // };

  class EqualExpression : public Expression {

  public:
  
    int spin;
    int value;

    EqualExpression(Variable X, Variable Y, const int sp=1);
    EqualExpression(Variable X, const int y, const int sp=1);
    virtual ~EqualExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  class EqualSetExpression : public Expression {

  public:
  
    int spin;
    int value;

    EqualSetExpression(Variable X, Variable Y, const int sp=1);
    EqualSetExpression(Variable X, const int y, const int sp=1);
    virtual ~EqualSetExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  class PrecedenceExpression : public Expression {

  public:

    int spin;
    int offset;

    PrecedenceExpression(Variable X, 
			 const int of=0, const int sp=1);
    PrecedenceExpression(Variable X, Variable Y, 
			 const int of=0, const int sp=1);
    virtual ~PrecedenceExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable Precedence(Variable X, const int d, Variable Y);

  class DisjunctiveExpression : public Expression {

  public:

    int processing_time[2];

    DisjunctiveExpression(Variable X, 
  			  Variable Y,
  			  const int p0=1, 
  			  const int p1=1);

    virtual ~DisjunctiveExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable Disjunctive(Variable X, Variable Y, const int px, const int py);
  

  class ReifiedDisjunctiveExpression : public Expression {

  public:

    int processing_time[2];

    ReifiedDisjunctiveExpression(Variable X, 
  				 Variable Y,
  				 const int p0=1, 
  				 const int p1=1);
    
    virtual ~ReifiedDisjunctiveExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable ReifiedDisjunctive(Variable X, Variable Y, const int px, const int py);
  

  class FreeExpression : public Expression {

  public:
    FreeExpression(Variable X);
    
    virtual ~FreeExpression();

    virtual const char* get_name() const;

  };

  Variable Free(Variable X);
  

  class AllDiffExpression : public Expression {

  public:

    int consistency_level;

    AllDiffExpression(Vector< Variable >& args, const int ct);
    virtual ~AllDiffExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable AllDiff(Vector< Variable >& args, const int ct=BOUND_CONSISTENCY);


  class LexExpression : public Expression {

  public:

    int strict;
    //int row_size;

    LexExpression(Vector< Variable >& r1, Vector< Variable >& r2, const int st_);
    virtual ~LexExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable LexLess(VarArray& r1, VarArray& r2);
  Variable LexLeq(VarArray& r1, VarArray& r2);


  class BoolSumExpression : public Expression {

  public:

    int lb;
    int ub;

    BoolSumExpression() : Expression() {}
    BoolSumExpression(const int l, const int u);
    BoolSumExpression(Vector< Variable >& args, const int l, const int u);
    virtual ~BoolSumExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable BoolSum(Vector< Variable >& args);
  Variable BoolSum(Vector< Variable >& args, const int l, const int u);
  Variable BoolSum(std::vector< Variable >& args, const int l, const int u);



  class LinearExpression : public Expression {

  public:

    int lower_bound;
    int upper_bound;
    Vector< int > weight;


    LinearExpression(Vector< Variable >& args, Vector< int >& wgt, const int l, const int u);
    LinearExpression(std::vector< Variable >& args, std::vector< int >& wgt, const int l, const int u);
    virtual ~LinearExpression();
    void initialise_bounds();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  //Variable Sum(Vector< Variable >& args);
  Variable Sum(Vector< Variable >& args, Vector< int >& wgts, Variable T);
  Variable Sum(std::vector< Variable >& args, std::vector< int >& wgts, Variable T);
  Variable Sum(Vector< Variable >& args, Vector< int >& wgts, const int l=-INFTY, const int u=INFTY);
  Variable Sum(std::vector< Variable >& args, std::vector< int >& wgts, const int l=-INFTY, const int u=INFTY);



  class ElementExpression : public Expression {

  public:

    int lower_bound;
    int upper_bound;
    int offset;
    //BitSet domain;
    Vector< int > values;

    ElementExpression(Vector< Variable >& args, Variable X, int ofs);
    virtual ~ElementExpression();
    void initialise_domain();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable Element(Vector<Variable>& X, Variable selector, int offset=0);
  Variable Element(VarArray& X, Variable selector, int offset=0);


  class MinExpression : public Expression {

  public:

    MinExpression() : Expression() {};
    MinExpression(Vector< Variable >& args);
    virtual ~MinExpression();
    //void initialise_domain();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable Min(Vector<Variable>& X);
  Variable Min(VarArray& X);
  Variable Min(Variable X, Variable Y);

  class MaxExpression : public Expression {

  public:

    MaxExpression() : Expression() {};
    MaxExpression(Vector< Variable >& args);
    virtual ~MaxExpression();
    //void initialise_domain();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable Max(Vector<Variable>& X);
  Variable Max(VarArray& X);
  Variable Max(Variable X, Variable Y);


  class SetExpression : public BoolSumExpression {
    
  public:
    
    //VarArray elements;
    Vector< int > elts_ub;
    Vector< int > elts_lb;
    int num_args;
    // BitSet required;
    // BitSet allowed;
 
    SetExpression() : BoolSumExpression() {}
    SetExpression(const int lelt, const int uelt, const int clb, const int cub);
    // SetExpression(const BitSet& lb, const BitSet& ub, const int clb, const int cub);
    SetExpression(const Vector<int>& lb, const Vector<int>& ub, const int clb, const int cub);
    void initialise_elements();
    virtual ~SetExpression();
    
    Variable get_index_var(const int idx) {return children[num_args+idx];}
    Variable get_elt_var(const int vali);
    int get_element_index(const int vali);
    
    // bool can_be_in(const int vali) { return allowed.contain(vali); }
    // bool must_be_in(const int vali) { return required.contain(vali); }
    // int get_smallest_elt() { return elements.front(); }
    // int get_highest_elt() { return elements.back(); }

    virtual const char* get_name() const;
    virtual std::ostream& display(std::ostream& os) const;
    virtual int get_solution_int_value() const ;
    virtual std::string get_solution_str_value() const ;
  };

  Variable SetVariable(const int lelt, const int uelt, 
  		       const int clb=0, const int cub=INFTY);
  Variable SetVariable(const Vector<int>& lb, const Vector<int>& ub, 
  		       const int clb, const int cub);
  // Variable SetVariable(std::vector<int> lb, std::vector<int> ub, const int clb, const int cub);
  // Variable SetVariable(BitSet lb, Bitset ub, const int clb, const int cub);

  Variable Card(Variable S); //{ return S; }


  class ElementSetExpression : public SetExpression {

  public:

    int offset;
    //VarArray elements;
    //Vector< int > values;

    ElementSetExpression(Vector< Variable >& args, Variable X, int ofs);
    virtual ~ElementSetExpression();
    void initialise_domain();

    //virtual Variable get_index_var(const int idx);

    virtual void extract_constraint(Solver*);
    //virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable ElementSet(Vector<Variable>& X, Variable selector, int offset=0);
  Variable ElementSet(VarArray& X, Variable selector, int offset=0);



  class IntersectionExpression : public SetExpression {

  public:

    IntersectionExpression(Variable X, Variable Y);
    virtual ~IntersectionExpression();
    void initialise_domain();

    virtual void extract_constraint(Solver*);
    //virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable Intersection(Variable X, Variable Y);



  class SubsetExpression : public Expression {

  public:

    SubsetExpression(Variable X, Variable Y);
    virtual ~SubsetExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable Subset(Variable X, Variable Y);


  class MemberExpression : public Expression {

  public:

    int lb;
    int ub;
    int size;
    BitSet values;

    MemberExpression(Variable X, Variable Y);
    MemberExpression(Variable X, const int lo, const int up);
    MemberExpression(Variable X, const BitSet& s);
    MemberExpression(Variable X, const Vector<int>& s);
    virtual ~MemberExpression();

    virtual void extract_constraint(Solver*);
    virtual void extract_variable(Solver*);
    virtual void extract_predicate(Solver*);
    virtual const char* get_name() const;

  };

  Variable Member(Variable X, Variable Y);
  Variable Member(Variable X, const int lo, const int up);
  Variable Member(Variable X, const BitSet& s);
  Variable Member(Variable X, const Vector<int>& s);


  class VarArray : public Vector< Variable > {
  public:
    
    VarArray() : Vector< Variable >() {}
    VarArray(const int n, int lb=NOVAL, int ub=NOVAL, int type=EXPRESSION) 
      : Vector< Variable >() 
    {
      if(lb==NOVAL) { lb=0; ub=1; }
      else if(ub==NOVAL) { lb=0; ub=lb-1; }
      
      for(int i=0; i<n; ++i) {
	Variable x(lb, ub, type);
	add(x);
      }
    }

    virtual ~VarArray() {}
    Variable operator[](Variable X);
    Variable operator[](const int X);
    //Variable& operator[](const int X);
    void set(const int i, Variable x);


    Variable operator<(VarArray& X);
    Variable operator<=(VarArray& X);
    Variable operator>(VarArray& X);
    Variable operator>=(VarArray& X);
  };


  class Goal {

  public:
    
    enum method { MINIMIZATION, MAXIMIZATION, SATISFACTION, ALLSOLUTIONS };
    method            type;
    int        lower_bound;
    int        upper_bound;
    Variable     objective;

    Goal(method t); 
    Goal(method t, Variable X);

    virtual ~Goal();

    virtual std::ostream& display(std::ostream& os) const;

    int value() const;
    bool enforce();
    Outcome notify_solution(Solver *solver);
    Outcome notify_exhausted();
     
  };


  std::ostream& operator<< (std::ostream& os, const Goal& x);
  std::ostream& operator<< (std::ostream& os, const Goal* x);

  //   std::ostream& operator<< (std::ostream& os, const VarArray& x);
  //   std::ostream& operator<< (std::ostream& os, const VarArray* x);


  //   template < class WORD_TYPE >
  //   /// Remove all values that do not appear in the current domain of the Variable "x"
  //   Event VariableBitset< WORD_TYPE >::set_domain(Variable x) {
  //     if(domain.size == 1) {
  //       return(x.contain(domain.min) ? NO_EVENT : FAIL_EVENT);
  //     } else {
  //       int xsize = x.get_size();
  //       if(xsize == 1) {
  // 	return set_domain(x.get_value());
  //       } else {
  // 	int xmin = x.get_min();
  // 	int xmax = x.get_max();
  // 	if(xsize == (xmax-xmin+1)) {
  // 	  return(set_min(xmin) | set_max(xmax));
  // 	} else if(x.includes(domain.values)) {
  // 	  return NO_EVENT;
  // 	} else if(x.intersect(domain.values)) {
  // 	  Event intersection = DOMAIN_EVENT;
  // 	  save();
  // 	  x.intersect_to(domain.values);
  // 	  domain.size = domain.values.size();
  // 	  if(!domain.values.contain(domain.min)) { intersection |= LB_EVENT; domain.min = domain.values.next(domain.min); }
  // 	  if(!domain.values.contain(domain.max)) { intersection |= UB_EVENT; domain.max = domain.values.prev(domain.max); }
  // 	  if(domain.min == domain.max) intersection |= VALUE_EVENT;
	  
  // 	  solver->trigger_event(id, intersection);
	  
  // 	} else return FAIL_EVENT;
  //       }
  //     }
  //     return NO_EVENT;
  //   }


  //     /// Whether the domain has a nonempty intersection with the Variable x
  //   template < class WORD_TYPE >
  //   bool VariableBitset< WORD_TYPE >::intersect(Variable x) const { return x.intersect(domain.values); }
  //     /// Whether the domain is included in the Variable x 
  //   template < class WORD_TYPE >
  //   bool VariableBitset< WORD_TYPE >::included(Variable x) const { return x.includes(domain.values); }
  //     /// Whether the domain is included in the Variable x 
  //   template < class WORD_TYPE >
  //   bool VariableBitset< WORD_TYPE >::includes(Variable x) const { return x.included(domain.values); }


  //Event VariableRange::remove(const int v);
  // inline Event VariableRange::remove(const int v) {
  //   Event removal = DOMAIN_EVENT;
    
  //   // first check if we can abort early
  //   if(min>v || max<v) {
  //     return NO_EVENT;
  //   }
  //   if(min!=v && max!=v) {
  //     solver->make_non_convex(id);
  //     removal = solver->variables[id].remove(v);
  //     return removal;
  //     //return NO_EVENT;
  //   }
  //   if(min==max) return FAIL_EVENT;

  //     save();

  //     if(min==v) {
  // 	++min;
  // 	removal |= LB_EVENT;
  //     } else {
  // 	--max;
  // 	removal |= UB_EVENT;
  //     }

  //     if(min == max) removal |= VALUE_EVENT; 
  //     solver->trigger_event(id, removal);
  //     return removal; 
  //   }


  /// Remove all values that do not appear in the set "s"
  //Event VariableRange::set_domain(const BitSet& s) ; 
   // /// Remove all values that do not appear in the set "s"
   //  inline Event VariableRange::set_domain(const BitSet& s) {
   //    Event setdomain = NO_EVENT;

   //    // std::cout << "include " << s.includes(min, max) << std::endl;
   //    // std::cout << "intersect " << s.intersect(min, max) << std::endl;
   //    // std::cout << "[" << s.next(min-1) << ".." << s.prev(max+1) << "]" << std::endl;
      

   //    if(s.includes(min, max)) return NO_EVENT;
   //    if(!s.intersect(min, max)) return FAIL_EVENT;
   //    int lb = s.next(min-1);
   //    int ub = s.prev(max+1);
   //    if(s.includes(lb, ub)) {
   // 	if(lb>min) {
   // 	  setdomain |= set_min(lb);
   // 	} 
   // 	if(ub<max) {
   // 	  setdomain |= set_max(ub);
   // 	}
   //    } else {

   // 	// std::cout << "the intersection is not convex" << std::endl;

   // 	solver->make_non_convex(id);


   // 	// std::cout << solver->variables[id] << " in " 
   // 	// 	  << solver->variables[id].get_domain() << std::endl;

   // 	return solver->variables[id].set_domain(s);
   //    }

   //    //return set_domain(s.next(min-1), s.prev(max+1));
   //    return setdomain;
   //  }


}

#endif // _VARIABLE_2_HPP
