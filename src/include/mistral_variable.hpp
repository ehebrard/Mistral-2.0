
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


#include <mistral_global.hpp>
#include <mistral_solver.hpp>
#include <mistral_structure.hpp>


#ifndef _MISTRAL_VARIABLE_HPP
#define _MISTRAL_VARIABLE_HPP

//class Solver;
namespace Mistral {

  /*!
    The 'Variable' class involves:
    1/ a pointer to an implementation;
    2/ a type flag;
    3/ a list of accessors;
    
    When calling an accessor, the type flag is checked,
    the pointer is cast accordingly and called with the same accessor.
  */  


  class VariableImplementation {
  public:

    /// Attributes common to all variables:
    /////////////////////////////////////
    /// The linked solver
    Solver *solver;
    /// unique identifier, corresponds to its rank in the variables list of the solver
    int id;
    /////////////////////////////////////

    VariableImplementation() {solver=NULL; id=-1;}
    virtual ~VariableImplementation() {}

    bool is_initialised() { 
      return id>-1 // && solver
	;
      //(self.domain_type != DYN_VAR && self.implementation != NULL);
    }

    int get_solution_value() const ;//{ return solver->last_solution_lb[id] } ; 
    int get_solution_lb() const ;//{ return solver->last_solution_lb[id] } ; 
    int get_solution_ub() const ;//{ return solver->last_solution_ub[id] } ; 

    //inline void triggerValueEvent() {solver->trigger_event(id, VALUE_EVENT);}
    void trigger_value_event_and_save(Variable *x);
  };


  class BitsetDomain {
  
  public:
  
    BitsetDomain() {min = NOVAL; max = NOVAL; size = 0;}
    BitsetDomain(const int lb, const int ub);
    //BitsetDomain(const int lb, const int ub, const int size, BitSet vals);
    void initialise(const int lb, const int ub, const bool vals=true);
  
    int min;
    int max;
    int size;
    BitSet values;
  
    virtual std::ostream& display(std::ostream& os) const {
      if(min == max-1) {
	os << "{" << min << "," << max << "}";
      } else if(values.table) {
	os << values;
      } else {
	os << "[" << min << ".." << max << "]";
      }
      return os;
    }
  
  };

  std::ostream& operator<< (std::ostream& os, const BitsetDomain& x);
  std::ostream& operator<< (std::ostream& os, const BitsetDomain* x);

  //class Variable;

  class Expression;
  class Variable {

  public:

    int domain_type;
    VariableImplementation *implementation;

    Variable();
    Variable(const int value);
    Variable(VariableImplementation* impl, const int type=DYN_VAR);
    Variable(Expression* impl);
    Variable(const int lo, const int up, const int type=DYN_VAR);

    //Variable get_children();
    Variable get_var();

    Variable operator+(Variable);
    Variable operator-(Variable);
//     //Variable operator*(Variable);
//     //Variable operator/(Variable);
//     //Variable operator&&(Variable);
//     //Variable operator||(Variable);
//     //Variable operator^(Variable);
//     //Variable operator->(Variable);
//     Variable operator==(Variable);
    Variable operator!=(Variable);
    Variable operator<(Variable);
    Variable operator<=(Variable);
    Variable operator>(Variable);
    Variable operator>=(Variable);
//     //Variable operator-();
//     //Variable operator!();


    void initialise(Solver *s, const bool top=true);
//    void initialise(Solver *s);

    Event setValue( const int val );    
    inline int id() const {return implementation->id;}
    inline Solver* get_solver() {return implementation->solver;}

    /*!@name Constant Accessors and Iterators*/
    //@{
    /// Returns the assigned value if it exists
    int get_value() const ; 
    int get_solution_value() const ; 
    /// Returns the domain 
    //BitsetDomain get_domain() const ; 
    std::string get_domain() const ; 
    /// Returns the domain size
    unsigned int get_size() const ; 
    /// Returns the minimum value in the domain
    int get_min() const ;
    /// Returns the maximum value in the domain
    int get_max() const ; 
    /// Returns the minimum value that could belong to the domain
    int get_minCapacity() const ;
    /// Returns 1 + the maximum value that could belong to the domain
    int get_maxCapacity() const ;
    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    int get_minPosAbs() const ;
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    int get_minNegAbs() const ;
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    int next(const int v) const ;
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    int prev(const int v) const ;
    /// Whether or not the Variable is currently an interval
    bool is_range() const ;
    /// Whether or not the Variable is bound to a ground value
    bool is_ground() const ;
    /// Whether or not the Variable is bound to a given ground value
    bool equal(const int v) ;
    /// Whether the value "v" is currently contained in the domain
    bool contain(const int v) ;

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
    void intersectTo( BitSet& s ) const ;
    /// Do the union of its domain with a set s
    void unionTo( BitSet& s ) const ;
    //@}

    /*!@name Domain changing methods*/
    //@{
    /// Remove value "v"
    Event remove(const int v) ;
    /// Remove all values but "v"
    Event setDomain(const int v) ;
    /// Remove all values strictly lower than lo
    Event setMin(const int lo) ;
    /// Remove all values strictly greater than up
    Event setMax(const int up) ;
    /// Remove all values that do not appear in the set "s"
    Event setDomain(const BitSet& s) ;
    /// Remove all values that do not appear in the current domain of the Variable "x"
    Event setDomain(Variable& x) ;
    /// Remove all values that belong to the set "s"
    Event removeSet(const BitSet& s) ;
    /// Remove all values in the interval [l..u]
    Event removeRange(const int lo, const int up) ;
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


  //  std::string toString(const BitsetDomain& x);


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

    //     /// vector of variables to restore
    //     Vector<Variable> *store;

    /*!@name Constructors*/
    //@{
    VariableBitset(const int lb, const int ub) {
      id = -1;
      solver = NULL;
      initialise(lb, ub);
    };

    virtual void initialise(const int lb, const int ub) {
 
      domain.initialise(lb, ub);
      trail_.initialise(0,8);
      trail_.add(domain.min);
      trail_.add(domain.max);
      trail_.add(domain.size);
      trail_.add(-1);
      initialise_trail();
    }

    void initialise_trail() {
      int i, k;
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
      //delete [] delta_abs;
      //delete [] level_abs;
      //delete [] delta_;
      //delete [] level_;
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
    inline int get_minCapacity() const { return trail_[0]; }
    /// Returns 1 + the maximum value that could belong to the domain
    inline int get_maxCapacity() const { return trail_[1]; }
    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    inline int get_minPosAbs() const {
      if( domain.max < 0 ) return NOVAL;
      if( domain.min > 0 ) return domain.min;
      if( domain.values.fastContain(0) ) return 0;
      return domain.values.next( 0 );
    }
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    inline int get_minNegAbs() const  {
      if( domain.min > 0 ) return NOVAL;
      if( domain.max < 0 ) return domain.max;
      if( domain.values.fastContain(0) ) return 0;
      return domain.values.prev( 0 );
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
    inline bool contain(const int v) const { return (domain.min <= v && domain.max >= v && domain.values.fastContain(v)); }

    /// Whether the domain has a nonempty intersection with the interval [l..u]
    inline bool intersect(const int lo, const int up) const { return (domain.min <= up && domain.max >= lo); }
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
    inline void intersectTo( BitSet& s ) const { s.intersectWith(domain.values); }
    /// Do the union of its domain with a set s
    inline void unionTo( BitSet& s ) const { s.unionWith(domain.values); }
    //@}

    /*!@name Domain handling methods*/
    //@{
    


    inline Event remove(const int v) {
      Event removal = DOMAIN_EVENT;

      // first check if we can abort early
      if(!contain(v)) return NO_EVENT;
      if(domain.size == 1) return FAIL_EVENT;

      save();

      // then change the static domain
      if(--domain.size == 1) removal = VALUE_EVENT; 
      domain.values.fastErase(v);
	
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

      solver->trigger_event(id, removal);
      return removal; 
    }

    /// Remove all values but "v"
    inline Event setDomain(const int v) {
      Event setdomain = VALUE_EVENT;

      // first check if we can abort early
      if(!contain(v)) return FAIL_EVENT;
      else if(domain.size == 1) return NO_EVENT;

      setdomain = VALUE_EVENT;
      save();

      if(domain.values.table) 
	domain.values.setTo(v);
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

      return setdomain; 
    }

    /// Remove all values strictly lower than l
    inline Event setMin(const int lo) {
      Event lower_bound = LB_EVENT;

      // first check if we can abort early
      if(domain.max <  lo) return FAIL_EVENT;
      if(domain.min >= lo) return NO_EVENT;
      
      save();
      
      // then change the static domain
      domain.values.setMin(lo);
      if(lo == domain.max) {
	domain.min = lo;
	domain.size = 1;
	lower_bound |= VALUE_EVENT;
      } else {
	domain.size = domain.values.size();
	if(domain.values.contain(lo)) domain.min = lo;
	else domain.min = domain.values.next(lo-1);
	if(domain.size == 1) lower_bound |= VALUE_EVENT;
      }
      
      solver->trigger_event(id, lower_bound);
      return lower_bound; 
    }

    /// Remove all values strictly greater than u
    inline Event setMax(const int up) {
      Event upper_bound = UB_EVENT;

      // first check if we can abort early
      if(domain.min >  up) return FAIL_EVENT;
      if(domain.max <= up) return NO_EVENT;
      
      save();

      // then change the static domain
      domain.values.setMax(up);
      if(up == domain.min) {
	domain.max = up;
	domain.size = 1;
	upper_bound |= VALUE_EVENT;
      } else {
	domain.size = domain.values.size();
	if(domain.values.contain(up)) domain.max = up;
	else domain.max = domain.values.prev(up);
	if(domain.size == 1) upper_bound |= VALUE_EVENT;
      }

      solver->trigger_event(id, upper_bound);
      return upper_bound; 
    }


    /// Remove all values that do not appear in the set "s"
    inline Event setDomain(const BitSet& s) {
      Event intersection = DOMAIN_EVENT;

      if( !domain.values.intersect(s) ) return FAIL_EVENT;
      if( domain.values.included(s) ) return NO_EVENT;

      save();
	  
      // then change the static domain
      domain.values.intersectWith(s);
      domain.size = domain.values.size();
      if(!s.contain(domain.min)) { intersection |= LB_EVENT; domain.min = domain.values.next(domain.min); }
      if(!s.contain(domain.max)) { intersection |= UB_EVENT; domain.max = domain.values.prev(domain.max); }
      if(domain.min == domain.max) intersection |= VALUE_EVENT;
    
      solver->trigger_event(id, intersection);
      return intersection; 
    }


    //    Event setDomain(Variable x) ;
    /// Remove all values that do not appear in the current domain of the Variable "x"
    inline Event setDomain(Variable x) {
      if(domain.size == 1) {
	return(x.contain(domain.min) ? NO_EVENT : FAIL_EVENT);
      } else {
	int xsize = x.get_size();
	if(xsize == 1) {
	  return setDomain(x.get_value());
	} else {
	  int xmin = x.get_min();
	  int xmax = x.get_max();
	  if(xsize == (xmax-xmin+1)) {
	    return(setMin(xmin) | setMax(xmax));
	  } else if(x.includes(domain.values)) {
	    return NO_EVENT;
	  } else if(x.intersect(domain.values)) {
	    Event intersection = DOMAIN_EVENT;
	    save();
	    x.intersectTo(domain.values);
	    domain.size = domain.values.size();
	    if(!domain.values.contain(domain.min)) { intersection |= LB_EVENT; domain.min = domain.values.next(domain.min); }
	    if(!domain.values.contain(domain.max)) { intersection |= UB_EVENT; domain.max = domain.values.prev(domain.max); }
	    if(domain.min == domain.max) intersection |= VALUE_EVENT;
    
	    solver->trigger_event(id, intersection);
	    
	  } else return FAIL_EVENT;
	}
      }
      return NO_EVENT;
    }

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
      // 	  domain.values.setminusWith(s);
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
    inline Event removeRange(const int lo, const int up) {
      if(lo <= domain.min) return setMin(up+1);
      if(up >= domain.max) return setMax(lo-1);
      if(domain.size != domain.max-domain.min+1) {
	if(!domain.values.intersect(lo, up)) return NO_EVENT;
	save();
	domain.values.removeInterval(lo, up);
	domain.size = domain.values.size();
      } else {
	save();
	domain.values.removeInterval(lo, up);
	domain.size -= (up-lo+1);
      }

      solver->trigger_event(id, DOMAIN_EVENT);

      return DOMAIN_EVENT; 
    }
    //@}


    /*!@name Virtual Accessors and Iterators*/
    //@{

    //@}   

    inline void save() {
      if(trail_.back() != solver->level) {
	//store();
	Variable self(this, BITSET_VAR);
	solver->save(self);
	//solver->save(this, BITSET_VAR);

	int i = domain.values.pos_words, j = domain.values.neg_words;

	trail_.add(domain.min);
	trail_.add(domain.max);
	trail_.add(domain.size);
	
	if(delta_) 
	  while( i --> j ) {
	    //WORD_TYPE buf = domain.values.table[i];
	    if(*(delta_[i]) != domain.values.table[i]) {
	      *(++delta_[i]) = domain.values.table[i];
	      *(++level_[i]) = solver->level;
	    }
	  }
	
	trail_.add(solver->level);
      }
    }

    inline Event restore() {
      // backtrack from solver->level to solver->level-1;
      // the domain has changed at solver->level
      int i = domain.values.pos_words, j = domain.values.neg_words;

      trail_.pop();
      trail_.pop(domain.size);
      trail_.pop(domain.max);
      trail_.pop(domain.min);

      if(delta_) 
	while( i --> j ) {
	  domain.values.table[i] = *(delta_[i]);
	  if(*(level_[i]) == solver->level) {
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
    
    //     virtual void debug_print() const; 
    //     {

    //       //Variable self((VariableImplementation*)this, BITSET_VAR);

    //       std::cout << this << " " 
    // 		<< domain << " "
    // 		<< trail_ << " " << std::endl;

    //       for(int i=domain.values.neg_words; i<domain.values.pos_words; ++i) {
    // 	std::cout << "word[" << i << "] " ;
    // 	for(WORD_TYPE* it = delta_abs[i]; it <= delta_[i]; ++it) {
    // 	  printBitset(*it, i, std::cout);
    // 	}
    // 	std::cout << std::endl;
    //       }
    //     }

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
      : VariableBitset< WORD_TYPE >(lb, ub) {}

    virtual void initialise(const int lb, const int ub) {
      int nw = (lb >> BitSet::EXP);
      int pw = (ub >> BitSet::EXP);
      assert((pw-nw) <= NUM_WORDS);
      VariableBitset< WORD_TYPE >::domain.initialise(lb, ub, false);
      VariableBitset< WORD_TYPE >::domain.values.neg_words = nw;
      VariableBitset< WORD_TYPE >::domain.values.pos_words = pw;
      VariableBitset< WORD_TYPE >::domain.values.table = ((&(word[0]))-nw);
      VariableBitset< WORD_TYPE >::initialise_trail();
    }

    virtual ~VariableWord() {
      VariableBitset< WORD_TYPE >::domain.values.table = NULL;
      VariableBitset< WORD_TYPE >::domain.values.neg_words = 0;
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
    inline int get_minCapacity() const { return trail_[0]; }
    /// Returns 1 + the maximum value that could belong to the domain
    inline int get_maxCapacity() const { return trail_[1]; }
    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    inline int get_minPosAbs() const {
      if( max < 0 ) return NOVAL;
      if( min > 0 ) return min;
      return 0;
    }
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    inline int get_minNegAbs() const  {
      if( min > 0 ) return NOVAL;
      if( max < 0 ) return max;
      return 0;
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
    inline bool intersect(const BitSet& s) const { return (min <= s.max() && max >= s.min()); }
    /// Whether the domain is included in the set s 
    inline bool included(const BitSet& s) const { return s.includes(min, max); }
    /// Whether the domain is included in the set s 
    inline bool includes(const BitSet& s) const { return s.included(min, max); }

    /// Whether the domain has a nonempty intersection with the Variable x
    inline bool intersect(const Variable x) const { return x.intersect(min,max); }
    /// Whether the domain is included in the Variable x 
    inline bool included(const Variable x) const { return x.includes(min,max); }
    /// Whether the domain is included in the Variable x 
    inline bool includes(const Variable x) const { return x.included(min,max); }

//     /// Whether the domain has a nonempty intersection with the Variable x
//     bool intersect(const Variable x) const ;
//     /// Whether the domain is included in the Variable x 
//     bool included(const Variable x) const ;
//     /// Whether the domain is included in the Variable x 
//     bool includes(const Variable x) const ;

    /// Intersect its domain with a set s
    inline void intersectTo( BitSet& s ) const { s.setMin(min); s.setMax(max); }
    /// Do the union of its domain with a set s
    inline void unionTo( BitSet& s ) const { // s.fill(min,max);
    }
    //@}

    /*!@name Domain handling methods*/
    //@{
    


    inline Event remove(const int v) {
      Event removal = DOMAIN_EVENT;

      // first check if we can abort early
      if(min!=v && max!=v) return NO_EVENT;
      if(min==max) return FAIL_EVENT;

      save();

      if(min==v) {
	++min;
	removal |= LB_EVENT;
      } else {
	--max;
	removal |= UB_EVENT;
      }

      if(min == max) removal |= VALUE_EVENT; 

//       display(std::cout);
//       std::cout << " in [" << min << ".." << max << "]: ";
//       std::cout << "TRIGGER " << removal << std::endl;
//       std::cout << (ASSIGNED(removal)) << std::endl;

      solver->trigger_event(id, removal);
      return removal; 
    }

    /// Remove all values but "v"
    inline Event setDomain(const int v) {
      Event setdomain = VALUE_EVENT;

      // first check if we can abort early
      if(!contain(v)) return FAIL_EVENT;
      if(min == max) return NO_EVENT;

      save();

      min = max = v;
      
      solver->trigger_event(id, setdomain);	

      return setdomain; 
    }

    /// Remove all values strictly lower than l
    inline Event setMin(const int lo) {
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
    inline Event setMax(const int up) {
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
    inline Event setDomain(const BitSet& s) {
      //return setDomain(s.next(min-1), s.prev(max+1));
      return NO_EVENT;
    }


    //    Event setDomain(Variable x) ;
    /// Remove all values that do not appear in the current domain of the Variable "x"
    inline Event setDomain(Variable x) {
      //return setDomain(x.next(min-1), x.prev(max+1));
      return NO_EVENT;
    }

    /// Remove all values that belong to the set "s"
    inline Event removeSet(const BitSet& s) {
      Event setdifference = NO_EVENT;
      return setdifference; // do nothing more for BOOL or CONSTANT or BITSET or RANGE
    }

    /// Remove all values in the interval [l..u]
    inline Event removeRange(const int lo, const int up) {
      if(lo <= min) return setMin(up+1);
      if(up >= max) return setMax(lo-1);
      return NO_EVENT; 
    }
    //@}


    /*!@name Virtual Accessors and Iterators*/
    //@{

    //@}   

    inline void save() {
      if(trail_.back() != solver->level) {
	Variable self(this, RANGE_VAR);
	solver->save(self);
	//solver->save(this, RANGE_VAR);
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
      os << "r" << id;
      return os;
    }


  };
  class VariableVirtual : public VariableBitmap {};




//   class Expression;
//   class Variable {

//   public:

//     int domain_type;
//     VariableImplementation *implementation;

//     Variable();
//     Variable(const int value);
//     Variable(VariableImplementation* impl, const int type=DYN_VAR);
//     Variable(Expression* impl);
//     Variable(const int lo, const int up, const int type=DYN_VAR);

//     //Variable get_children();
//     Variable get_var();

//     Variable operator+(Variable);
//     Variable operator-(Variable);
// //     //Variable operator*(Variable);
// //     //Variable operator/(Variable);
// //     //Variable operator&&(Variable);
// //     //Variable operator||(Variable);
// //     //Variable operator^(Variable);
// //     //Variable operator->(Variable);
// //     Variable operator==(Variable);
//     Variable operator!=(Variable);
//     Variable operator<(Variable);
//     Variable operator<=(Variable);
//     Variable operator>(Variable);
//     Variable operator>=(Variable);
// //     //Variable operator-();
// //     //Variable operator!();


//     void initialise(Solver *s, const bool top=true);
// //    void initialise(Solver *s);

//     Event setValue( const int val );    
//     inline int id() const {return implementation->id;}
//     inline Solver* get_solver() {return implementation->solver;}

//     /*!@name Constant Accessors and Iterators*/
//     //@{
//     /// Returns the assigned value if it exists
//     int get_value() const ; 
//     int get_solution_value() const ; 
//     /// Returns the domain 
//     //BitsetDomain get_domain() const ; 
//     std::string get_domain() const ; 
//     /// Returns the domain size
//     unsigned int get_size() const ; 
//     /// Returns the minimum value in the domain
//     inline int get_min() const {
//     //int Mistral::Variable::get_min() const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->get_min();
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->get_min();
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->get_min();
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->get_min();
//     else if(domain_type ==   CONST_VAR) return (int)implementation;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.get_min();
//     else  return (!(*((int*)domain_type) & 1));
//   }

//     /// Returns the maximum value in the domain
//     int get_max() const ; 
//     /// Returns the minimum value that could belong to the domain
//     int get_minCapacity() const ;
//     /// Returns 1 + the maximum value that could belong to the domain
//     int get_maxCapacity() const ;
//     /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
//     int get_minPosAbs() const ;
//     /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
//     int get_minNegAbs() const ;
//     /// Return the smallest value currently in the domain that is strictly greater than "v"
//     int next(const int v) const ;
//     /// Return the smallest value currently in the domain that is strictly greater than "v"
//     int prev(const int v) const ;
//     /// Whether or not the Variable is currently an interval
//     bool is_range() const ;
//     /// Whether or not the Variable is bound to a ground value
//     bool is_ground() const ;
//     /// Whether or not the Variable is bound to a given ground value
//     bool equal(const int v) ;
//     /// Whether the value "v" is currently contained in the domain
//     bool contain(const int v) ;

//     /// Whether the domain has a nonempty intersection with the interval [l..u]
//     bool intersect(const int lo, const int up) const ;
//     /// Whether the domain is included in the interval [l..u]
//     bool included(const int lo, const int up) const ;
//     /// Whether the domain is included in the interval [l..u]
//     bool includes(const int lo, const int up) const ;

//     /// Whether the domain has a nonempty intersection with the set s 
//     bool intersect(const BitSet& s) const ;
//     /// Whether the domain is included in the set s 
//     bool included(const BitSet& s) const ;
//     /// Whether the domain is included in the set s 
//     bool includes(const BitSet& s) const ;

//     /// Whether the domain has a nonempty intersection with the Variable x
//     bool intersect(const Variable& x) const ;
//     /// Whether the domain is included in the Variable x 
//     bool included(const Variable& x) const ;
//     /// Whether the domain is included in the Variable x 
//     bool includes(const Variable& x) const ;

//     /// Intersect its domain with a set s
//     void intersectTo( BitSet& s ) const ;
//     /// Do the union of its domain with a set s
//     void unionTo( BitSet& s ) const ;
//     //@}

//     /*!@name Domain changing methods*/
//     //@{
//     /// Remove value "v"
//     Event remove(const int v) ;
//     /// Remove all values but "v"
//     Event setDomain(const int v) ;
//     /// Remove all values strictly lower than lo
//     Event setMin(const int lo) ;
//     /// Remove all values strictly greater than up
//     Event setMax(const int up) ;
//     /// Remove all values that do not appear in the set "s"
//     Event setDomain(const BitSet& s) ;
//     /// Remove all values that do not appear in the current domain of the Variable "x"
//     Event setDomain(Variable& x) ;
//     /// Remove all values that belong to the set "s"
//     Event removeSet(const BitSet& s) ;
//     /// Remove all values in the interval [l..u]
//     Event removeRange(const int lo, const int up) ;
//     //@}

//     /*!@name Backtracking methods*/
//     //@{
//     /// Restore the last solved state
//     Event restore();
//     //@}

//     /*!@name printing methods*/
//     //@{
//     //void full_print() const ;
//     std::ostream& display(std::ostream& os) const ;
//     //void debug_print() const ;
//     //@}
//   };

//   std::ostream& operator<< (std::ostream& os, const Variable& x);
//   std::ostream& operator<< (std::ostream& os, const Variable* x);


//   //  std::string toString(const BitsetDomain& x);




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

    //     inline bool satisfied() {
    //       switch(type()) {
    //       case REMOVAL: return !var.contain(value());
    //       case ASSIGNMENT: return var.equal(value());
    //       case LOWERBOUND: return var.get_min() >  value(); // > v
    //       case UPPERBOUND: return var.get_max() <= value();   // <= v
    //       }
    //       return true;
    //     }

    //     inline bool violated() {
    //       switch(type()) {
    //       case REMOVAL: return var.equal(value());
    //       case ASSIGNMENT: return !var.contain(value());
    //       case LOWERBOUND: return var.get_max() <= value(); // > v
    //       case UPPERBOUND: return var.get_min() >  value();   // <= v
    //       }
    //       return true;
    //     }
    
    Decision(Variable x, const int t, const int v) {
      init_data(t,v);
      var = x;
    }
    
    void init_data(const int t, const int v) {
      _data_ = (((v - (t == LOWERBOUND)) << 2) | t);
    }

    Decision(Constraint *c) {
      _data_ = -1;
      var = Variable();
      var.implementation = (VariableImplementation*)c;
    }

    virtual ~Decision() {
    }

    inline void invert() { _data_^=1; }
    
    inline bool make() {

      if(_data_ == -1) return propagateRelation();
      switch(type()) {
      case REMOVAL:    return !IS_FAIL(var.remove(value()));
      case ASSIGNMENT: return !IS_FAIL(var.setDomain(value()));
      case LOWERBOUND: return !IS_FAIL(var.setMin(value()+1));
      case UPPERBOUND: return !IS_FAIL(var.setMax(value()));
	//       case REMOVAL: return (var.remove(value()) == FAIL_EVENT ? var : ok);
	//       case ASSIGNMENT: return (var.setDomain(value()) == FAIL_EVENT ? var : ok);
	//       case LOWERBOUND: return (var.setMin(value()+1) == FAIL_EVENT ? var : ok); // > v
	//       case UPPERBOUND: return (var.setMax(value()) == FAIL_EVENT ? var : ok);   // <= v
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


  // 1/ creation
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

    Expression() : VariableImplementation() {}
  Expression(Vector< Variable >& args);
  ~Expression();

  
  virtual void extract_constraint(Solver*)=0;
  virtual void extract_predicate(Solver*)=0;
  virtual void extract_variable(Solver*)=0;
  virtual const char* get_name()=0;

//   /// Returns the assigned value if it exists
//   int get_value() const { return self.get_value(); } 
//   int get_solution_value() const { return self.get_solution_value(); } 
//   /// Returns the domain 
//   //BitsetDomain get_domain() const ; 
//   std::string get_domain() const { return self.get_domain(); }
//   /// Returns the domain size
//     unsigned int get_size() const { return self.get_size(); }
//     /// Returns the minimum value in the domain
//     int get_min() const { return self.get_min(); }
//     /// Returns the maximum value in the domain
//     int get_max() const { return self.get_max(); }
//     /// Returns the minimum value that could belong to the domain
//     int get_minCapacity() const { return self.get_minCapacity(); }
//     /// Returns 1 + the maximum value that could belong to the domain
//     int get_maxCapacity() const { return self.get_minCapacity(); }
//     /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
//     int get_minPosAbs() const { return self.get_minPosAbs(); }
//     /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
//     int get_minNegAbs() const { return self.get_minNegAbs(); }
//     /// Return the smallest value currently in the domain that is strictly greater than "v"
//     int next(const int v) const { return self.next(v); }
//     /// Return the smallest value currently in the domain that is strictly greater than "v"
//     int prev(const int v) const { return self.prev(v); }
//     /// Whether or not the Variable is currently an interval
//     bool is_range() const { return self.is_range(); }
//     /// Whether or not the Variable is bound to a ground value
//     bool is_ground() const { return self.is_ground(); }
//     /// Whether or not the Variable is bound to a given ground value
//     bool equal(const int v) { return self.equal(v); }
//     /// Whether the value "v" is currently contained in the domain
//     bool contain(const int v) { return self.contain(v); }

//     /// Whether the domain has a nonempty intersection with the interval [l..u]
//   bool intersect(const int lo, const int up) const { return self.intersect(lo, up); }
//     /// Whether the domain is included in the interval [l..u]
//     bool included(const int lo, const int up) const ;
//     /// Whether the domain is included in the interval [l..u]
//     bool includes(const int lo, const int up) const ;

//     /// Whether the domain has a nonempty intersection with the set s 
//     bool intersect(const BitSet& s) const ;
//     /// Whether the domain is included in the set s 
//     bool included(const BitSet& s) const ;
//     /// Whether the domain is included in the set s 
//     bool includes(const BitSet& s) const ;

//     /// Whether the domain has a nonempty intersection with the Variable x
//     bool intersect(const Variable& x) const ;
//     /// Whether the domain is included in the Variable x 
//     bool included(const Variable& x) const ;
//     /// Whether the domain is included in the Variable x 
//     bool includes(const Variable& x) const ;

//     /// Intersect its domain with a set s
//     void intersectTo( BitSet& s ) const ;
//     /// Do the union of its domain with a set s
//     void unionTo( BitSet& s ) const ;
//     //@}

//     /*!@name Domain changing methods*/
//     //@{
//     /// Remove value "v"
//     Event remove(const int v) ;
//     /// Remove all values but "v"
//     Event setDomain(const int v) ;
//     /// Remove all values strictly lower than lo
//     Event setMin(const int lo) ;
//     /// Remove all values strictly greater than up
//     Event setMax(const int up) ;
//     /// Remove all values that do not appear in the set "s"
//     Event setDomain(const BitSet& s) ;
//     /// Remove all values that do not appear in the current domain of the Variable "x"
//     Event setDomain(Variable& x) ;
//     /// Remove all values that belong to the set "s"
//     Event removeSet(const BitSet& s) ;
//     /// Remove all values in the interval [l..u]
//     Event removeRange(const int lo, const int up) ;
//     //@}

//     /*!@name Backtracking methods*/
//     //@{
//     /// Restore the last solved state
//     Event restore();
//     //@}

//     /*!@name printing methods*/
//     //@{
//     //void full_print() const ;
//     std::ostream& display(std::ostream& os) const ;
//     //void debug_print() const ;
//     //@}



};

class BinaryExpression : public Expression {

public :

  BinaryExpression(Variable X, Variable Y);
  ~BinaryExpression();

};

class AddExpression : public BinaryExpression {

public:

  AddExpression(Variable X, Variable Y);
  ~AddExpression();

  virtual void extract_constraint(Solver*);
  virtual void extract_variable(Solver*);
  virtual void extract_predicate(Solver*);
  virtual const char* get_name();

};

class SubExpression : public BinaryExpression {

public:

  SubExpression(Variable X, Variable Y);
  ~SubExpression();

  virtual void extract_constraint(Solver*);
  virtual void extract_variable(Solver*);
  virtual void extract_predicate(Solver*);
  virtual const char* get_name();

};


class NeqExpression : public BinaryExpression {

public:

  NeqExpression(Variable X, Variable Y);
  ~NeqExpression();

  virtual void extract_constraint(Solver*);
  virtual void extract_variable(Solver*);
  virtual void extract_predicate(Solver*);
  virtual const char* get_name();

};

class PrecedenceExpression : public BinaryExpression {

public:
  int offset;

  PrecedenceExpression(Variable X, Variable Y, const int of=0);
  ~PrecedenceExpression();

  virtual void extract_constraint(Solver*);
  virtual void extract_variable(Solver*);
  virtual void extract_predicate(Solver*);
  virtual const char* get_name();

};


class AllDiffExpression : public Expression {

public:

  int consistency_level;

  AllDiffExpression(Vector< Variable >& args, const int ct);
  ~AllDiffExpression();

  virtual void extract_constraint(Solver*);
  virtual void extract_variable(Solver*);
  virtual void extract_predicate(Solver*);
  virtual const char* get_name();

};

  Variable AllDiff(Vector< Variable >& args, const int ct=BOUND_CONSISTENCY);



  class VarArray : public Vector< Variable > {
  public:
    
    VarArray() : Vector< Variable >() {}
    VarArray(const int n, int lb=NOVAL, int ub=NOVAL, int type=DYN_VAR) : Vector< Variable >() 
    {
      if(lb==NOVAL) { lb=0; ub=1; }
      else if(ub==NOVAL) { lb=0; ub=lb-1; }

      for(int i=0; i<n; ++i) {
	Variable x(lb, ub, type);
	add(x);
      }
    }

    ~VarArray() {}

    // std::ostream& display(std::ostream& os) const;

  };

//   std::ostream& operator<< (std::ostream& os, const VarArray& x);
//   std::ostream& operator<< (std::ostream& os, const VarArray* x);


//   template < class WORD_TYPE >
//   /// Remove all values that do not appear in the current domain of the Variable "x"
//   Event VariableBitset< WORD_TYPE >::setDomain(Variable x) {
//     if(domain.size == 1) {
//       return(x.contain(domain.min) ? NO_EVENT : FAIL_EVENT);
//     } else {
//       int xsize = x.get_size();
//       if(xsize == 1) {
// 	return setDomain(x.get_value());
//       } else {
// 	int xmin = x.get_min();
// 	int xmax = x.get_max();
// 	if(xsize == (xmax-xmin+1)) {
// 	  return(setMin(xmin) | setMax(xmax));
// 	} else if(x.includes(domain.values)) {
// 	  return NO_EVENT;
// 	} else if(x.intersect(domain.values)) {
// 	  Event intersection = DOMAIN_EVENT;
// 	  save();
// 	  x.intersectTo(domain.values);
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


  

}

#endif // _VARIABLE_2_HPP
