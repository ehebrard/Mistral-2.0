
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
#include <mistral_backtrack.hpp>
#include <mistral_structure.hpp>


#ifndef __VARIABLE_HPP
#define __VARIABLE_HPP


/*!
  A variable is composed (mainly) of
  a/ a domain
  c/ a list of constraints

  Editing a domain should work as follows:
  1/ check if it really is an edit (non-virtual)
  2/ check if it will not cause a wipe-out (non-virtual)
  3/ (as a rev struct) check if this is the first modif (non-virtual)
  4/ modify the domain (part non-virtual / part virtual) maybe check the type?
 */

namespace Mistral {

  class FiniteDomain {
    
  public:
    
    void initialise(const int lb, const int ub, const bool vals=true);

    int min;
    int max;
    int size;
    BitSet values;

    void _assert_();

//     std::string getString() const {
//       std::string return_str = "";
//       if(values.table)
// 	return_str += (toString(values));
//       else
// 	return_str += "["+toString(min)+".."+toString(max)+"]";
//       return return_str;
//     }
    virtual std::ostream& display(std::ostream& os) const {
      if(values.table)
	os << values;
      else
	os << "[" << min << ".." << max << "]";
      return os;
    }

  };


  //  std::string toString(const FiniteDomain& x);
  std::ostream& operator<< (std::ostream& os, const FiniteDomain& x);
  std::ostream& operator<< (std::ostream& os, const FiniteDomain* x);


  template < class WORD_TYPE >
  class VariableInt 
  {
    
  public:

    /// trail, used for backtracks
    Vector<int> trail_;

    /// Includes a bitset of values, a lower and an upper bound
    FiniteDomain domain;

    /// The list of all constraints whose scope includes this variable
    ConstraintList constraints;

    /// trail for the bitset representation
    //Vector< WORD_TYPE > *values_trail;
    //Vector< int > *level_trail;
    WORD_TYPE **delta_;
    int **level_;
    WORD_TYPE **delta_abs;
    int **level_abs;
    //int *num_changes;


    /// The linked solver
    Solver *solver;

    /// The stack of variables it belongs to
    Stack< IntVar > *var_list;

    /// type of the variable (NOTYPE, BITSET, BOOL, RANGE, LIST, VIRTUAL, CONSTANT)
    const int type;

    /// unique identifier, corresponds to its rank in the variables list of the solver
    int id;

    /// weight value, used for branching 
    float weight;

    void _assert_constraints_() {

      // check active constraints 
      if(domain.size > 1) {
	ConstraintNode nd = constraints.first(_value_);
	while( constraints.next(nd) )
	  {
	    bool is_active = false;
	    int i=0;
	    for(; !is_active && i<nd.elt.constraint->arity; ++i) {
	      is_active = (nd.elt.constraint->scope[i] != this && nd.elt.constraint->scope[i]->domain.size > 1);
	    }
	    if(!is_active) {
	      std::cerr << "Unactive constraint " << (nd.elt.constraint) << " => " << (this) << std::endl;

	      
	      std::cerr << solver << std::endl;

	      exit(1);
	    }
	  }
      }

      // check scope 
      Constraint *cons;
      for(unsigned int i=0; i<constraints.data.size; ++i) {
	cons = ((ConstraintTrigger)(constraints.data[i])).constraint;
	if(cons) {
	  bool is_in = false;
	  for(int i=0; !is_in && i<cons->arity; ++i) {
	    is_in = (cons->scope[i] == this);
	  }
	  if(!is_in) {
	    std::cerr << "Wrong scope" << std::endl;
	    exit(1);
	  }
	}
      }
      
    }
    

    /*!@name Constructors*/
    //@{
    VariableInt(const int lb, const int ub, const int TYPE=BITSET_VAR) : type(TYPE) {
      id = -1;
      weight = 0.0;
      solver = NULL;
      var_list = NULL;
      initialise(lb, ub, type);
    };

    virtual void initialise(const int lb, const int ub, const int var_type) {
 
      domain.initialise(lb, ub);
      trail_.initialise(0,8);
      trail_.add(domain.min);
      trail_.add(domain.max);
      trail_.add(domain.size);
      trail_.add(-1);


//       values_trail = new Vector< WORD_TYPE >[domain.values.pos_words-domain.values.neg_words];
//       values_trail -= domain.values.neg_words;
//       level_trail = new Vector< int >[domain.values.pos_words-domain.values.neg_words];
//       level_trail -= domain.values.neg_words;
//       for(int i=domain.values.neg_words; i<domain.values.pos_words; ++i) {
// 	std::cout << "init value trail[" << i << "]: " 
// 		  << (domain.values.size(i)) << " " 
// 		  << (domain.values.table[i]) << std::endl;
// 	values_trail[i].initialise(0, domain.values.size(i));
// 	values_trail[i].add(domain.values.table[i]);
// 	level_trail[i].initialise(0, domain.values.size(i));
// 	level_trail[i].add(0);
//       }
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
	k = domain.values.size(i);
	delta_[i] = new WORD_TYPE[k];
	delta_abs[i] = delta_[i];
	delta_[i][0] = domain.values.table[i];
	level_[i] = new int[k];
	level_abs[i] = level_[i];
	level_[i][0] = -1;
      }
    }

    virtual ~VariableInt() {
      for(int i=domain.values.neg_words; i<domain.values.pos_words; ++i) {
	delete [] delta_abs[i];
	delete [] level_abs[i];
      }
      delta_abs += domain.values.neg_words;
      level_abs += domain.values.neg_words;
      delta_ += domain.values.neg_words;
      level_ += domain.values.neg_words;
      delete [] delta_abs;
      delete [] level_abs;
      delete [] delta_;
      delete [] level_;
      //values_trail += domain.values.neg_words;
      //delete [] values_trail;
    }

    void _assert_() { domain._assert_(); }

    inline void unassign() {

      //if(this->id == 8 || this->id == 9 || this->id == 23) std::cout << "unassign " << this << std::endl;

      ConstraintNode nd = constraints.first(_value_);
      while( constraints.next(nd) )
	nd.elt.constraint->unassign(nd.elt.index);
      
      //std::cout << std::endl;
    }

    inline void triggerEvent(const Event evt) {
      ConstraintNode nd;

      if(is_value(evt)) {
	//solver->assign(this);

	//std::cout << "remove " << this << " from " << (*var_list) << std::endl;
	//if(this->id == 8 || this->id == 9 || this->id == 23) std::cout << "assign " << this << std::endl;

	  //<<  ((*var_list)[0]->id >= 12 ? " from aux" : " from seq") << std::endl;

	//assert(var_list->member(this));
	var_list->erase(this);
	nd = constraints.first(_value_);
      } else if(is_range(evt)) {
	nd = constraints.first(_range_);
      } else if(is_domain(evt)) {
	nd = constraints.first(_domain_);
      }

      if(is_value(evt)) while( constraints.next(nd) ) {
	  solver->trigger(nd.elt.constraint, nd.elt.index, evt);
	  nd.elt.constraint->assign(nd.elt.index);
	}
      else 
	while( constraints.next(nd) )
	  solver->trigger(nd.elt.constraint, nd.elt.index, evt);

      //if(is_value(evt)) 
      //std::cout << std::endl;

    }

    /*!@name Static Accessors and Iterators*/
    //@{
    /// Returns the assigned value if it exists
    inline int get_value() const { return solver->solution[id]; }
    /// Returns the minimum value in the domain
    inline int get_min() const { return domain.min; }
    /// Returns the maximum value in the domain
    inline int get_max() const { return domain.max; }
    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    inline int minPosAbs() const {
      if( domain.max < 0 ) return NOVAL;
      if( domain.min > 0 ) return domain.min;
      if( !domain.values.table || domain.values.fastMember(0) ) return 0;
      return domain.values.next( 0 );
    }
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    inline int minNegAbs() const  {
      if( domain.min > 0 ) return NOVAL;
      if( domain.max < 0 ) return domain.max;
      if( !domain.values.table || domain.values.fastMember(0) ) return 0;
      return domain.values.prev( 0 );
    } 
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int next(const int v) const { return ((domain.values.table && (domain.size != domain.max - domain.min + 1)) ? domain.values.next(v) : (v < domain.max ? v+1 : NOVAL)); }
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int prev(const int v) const { return ((domain.values.table && (domain.size != domain.max - domain.min + 1)) ? domain.values.prev(v) : (v > domain.min ? v-1 : NOVAL)); }


    /// Whether or not the Variable is currently an interval
    inline bool isRange() const { return (!domain.values.table || (domain.size == domain.max - domain.min + 1)); }
    /// Whether or not the Variable is bound to a ground value
    inline bool isGround() const { return (domain.min == domain.max); }
    /// Whether or not the Variable is bound to a given ground value
    inline bool equal(const int v) const { return (domain.min == v && domain.max == v); }
    /// Whether the value "v" is currently contained in the domain
    inline bool contain(const int v) const { return (domain.min <= v && domain.max >= v && (!domain.values.table || domain.values.fastMember(v))); }

    /// Whether the domain has a nonempty intersection with the interval [l..u]
    inline bool intersect(const int lo, const int up) const { return (domain.min <= up && domain.max >= lo); }
    /// Whether the domain is included in the interval [l..u]
    inline bool included(const int lo, const int up) const { return (domain.min >= lo && domain.max <= up); }
    /// Whether the domain is included in the interval [l..u]
    inline bool includes(const int lo, const int up) const { return (domain.min <= lo && domain.max >= up && (!domain.values.table || domain.values.includes(lo, up))); }

    /// Whether the domain has a nonempty intersection with the set s 
    inline bool intersect(const BitSet& s) const { return (domain.values.table ? domain.values.intersect(s) : s.intersect(domain.min, domain.max)); }
    /// Whether the domain is included in the set s 
    inline bool included(const BitSet& s) const { return (domain.values.table ? domain.values.included(s) : s.includes(domain.min, domain.max)); }
    /// Whether the domain is included in the set s 
    inline bool includes(const BitSet& s) const { return (domain.values.table ? domain.values.includes(s) : s.included(domain.min, domain.max)); }

    /// Whether the domain has a nonempty intersection with the Variable x
    inline bool intersect(const VariableInt* x) const { return (domain.values.table ? x->intersect(domain.values) : x->intersect(domain.min, domain.max)); }
    /// Whether the domain is included in the Variable x 
    inline bool included(const VariableInt* x) const { return (domain.values.table ? x->includes(domain.values) : x->includes(domain.min, domain.max)); }
    /// Whether the domain is included in the Variable x 
    inline bool includes(const VariableInt* x) const { return (domain.values.table ? x->included(domain.values) : x->included(domain.min, domain.max)); }

    /// Intersect its domain with a set s
    inline void intersectTo( BitSet& s ) const { if(domain.values.table) s.intersectWith(domain.values); else {s.setMin(domain.min); s.setMax(domain.max);} }
    /// Do the union of its domain with a set s
    inline void unionTo( BitSet& s ) const { (domain.values.table ? s.unionWith(domain.values) : s.addInterval(domain.min, domain.max)); }
    //@}

    /*!@name Domain handling methods*/
    //@{
    inline void _save_() {
      if(trail_.capacity) {
	if(trail_.back() != solver->level) {
	  solver->save(this);
	  save();
	}
      } else {
	solver->save(this);
      }
    }

    inline Event remove(const int v) {
      Event removal = DOMAIN_EVENT;

      if(type != VIRTUAL_VAR) {
	// first check if we can abort early
	if(!contain(v)) return NO_EVENT;
	if(domain.size == 1) return FAIL_EVENT;

	if(trail_.capacity) {
	  if(trail_.back() != solver->level) {
	    solver->save(this);
	    save();
	  }
	} else {
	  solver->save(this);
	}

	//_save_();

	// then change the static domain
	if(--domain.size == 1) removal = VALUE_EVENT; // : DOMAIN_EVENT);
	
	if(domain.values.table) 
	  domain.values.fastErase(v);
	
	if(removal == VALUE_EVENT)
	  if(domain.min == v) domain.min = domain.max;
	  else domain.max = domain.min;
	else if(domain.values.table) {
	  if(domain.max == v) {
	    removal |= UB_EVENT;
	    domain.max = domain.values.max();
	  } else if(domain.min == v) {
	    removal |= LB_EVENT;
	    domain.min = domain.values.min();
	  }       
	} else {
	  if(domain.max == v) { removal |= UB_EVENT; --domain.max; }
	  else if(domain.max == v) { removal |= LB_EVENT; ++domain.min; }
	}
      } 

      triggerEvent(removal);
      return removal; // do nothing more for BOOL or CONSTANT or BITSET or RANGE
    }

    /// Remove all values but "v"
    inline Event setDomain(const int v) {
      Event setdomain = VALUE_EVENT;

      if(type != VIRTUAL_VAR) {

	// first check if we can abort early
	if(!contain(v)) return FAIL_EVENT;

// 	if(trail_.capacity) {
// 	  if(trail_.back() != solver->level) {
// 	    //trail_.add(solver->level);
// 	    solver->save(this);
// 	    save();
// 	  }
// 	} else {
// 	  solver->save(this);
// 	}

	_save_();

	if(domain.size > 1) {
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
	}
      } 

      triggerEvent(setdomain);
      return setdomain; // do nothing more for BOOL or CONSTANT or BITSET or RANGE
    }

    /// Remove all values strictly lower than l
    inline Event setMin(const int low) {
      Event lower_bound = LB_EVENT;

      if(type != VIRTUAL_VAR) {
	// first check if we can abort early
	if(domain.max <  low) return FAIL_EVENT;
	if(domain.min >= low) return NO_EVENT;

// 	if(trail_.capacity) {
// 	  if(trail_.back() != solver->level) {
// 	    //trail_.add(solver->level);
// 	    solver->save(this);
// 	    save();
// 	  }
// 	} else {
// 	  solver->save(this);
// 	}

	_save_();

	// then change the static domain
	if(domain.values.table) {
	  domain.values.setMin(low);
	  if(low == domain.max) {
	    domain.min = low;
	    domain.size = 1;
	    lower_bound |= VALUE_EVENT;
	  } else {
	    domain.size = domain.values.size();
	    if(domain.values.member(low)) domain.min = low;
	    else domain.min = domain.values.next(low-1);
	    if(domain.size == 1) lower_bound |= VALUE_EVENT;
	  }
	} else {
	  domain.max = low;
	  domain.size = domain.max - low + 1;
	}
      } 

      triggerEvent(lower_bound);
      return lower_bound; // do nothing more for BOOL or CONSTANT or BITSET or RANGE
    }

    /// Remove all values strictly greater than u
    inline Event setMax(const int up) {
      Event upper_bound = UB_EVENT;

      if(type != VIRTUAL_VAR) {
	// first check if we can abort early
	if(domain.min >  up) return FAIL_EVENT;
	if(domain.max <= up) return NO_EVENT;

	_save_();

	// then change the static domain
	if(domain.values.table) {
	  domain.values.setMax(up);
	  if(up == domain.min) {
	    domain.max = up;
	    domain.size = 1;
	    upper_bound |= VALUE_EVENT;
	  } else {
	    domain.size = domain.values.size();
	    if(domain.values.member(up)) domain.max = up;
	    else domain.max = domain.values.prev(up);
	    if(domain.size == 1) upper_bound |= VALUE_EVENT;
	  }
	} else {
	  domain.max = up;
	  domain.size = up - domain.min + 1;
	}
      } 

      triggerEvent(upper_bound);
      return upper_bound; // do nothing more for BOOL or CONSTANT or BITSET or RANGE
    }


    /// Remove all values that do not appear in the set "s"
    inline Event setDomain(const BitSet& s) {
      Event intersection = DOMAIN_EVENT;

      if(type != VIRTUAL_VAR) {
	// first check if we can abort early
	if(domain.values.table) {
	  if( !domain.values.intersect(s) ) return FAIL_EVENT;
	  if( domain.values.included(s) ) return NO_EVENT;

	  _save_();
	  
	  // then change the static domain
	  domain.values.intersectWith(s);
	  domain.size = domain.values.size();
	  if(!s.member(domain.min)) { intersection |= LB_EVENT; domain.min = domain.values.next(domain.min); }
	  if(!s.member(domain.max)) { intersection |= UB_EVENT; domain.max = domain.values.prev(domain.max); }
	  if(domain.min == domain.max) intersection |= VALUE_EVENT;
	} else return setMin(s.min()) && setMax(s.max());
      } 

      triggerEvent(intersection);
      return intersection; // do nothing more for BOOL or CONSTANT or BITSET or RANGE
    }

    /// Remove all values that do not appear in the current domain of the Variable "x"
    inline Event setDomain(VariableInt* x) {
      FiniteDomain& xdom = x->domain;
      if(xdom.size == 1) return setDomain(xdom.min);
      else if(xdom.size == xdom.max - xdom.min + 1)
	return (setMin(xdom.min) & setMax(xdom.max));
      else return setDomain(xdom.values);
    }

    /// Remove all values that belong to the set "s"
    inline Event removeSet(const BitSet& s) {
      Event setdifference = NO_EVENT;

      if(type != VIRTUAL_VAR) {
	// first check if we can abort early
	if(domain.values.table) {
	  if( !domain.values.intersect(s) ) return NO_EVENT;
	  if( domain.values.included(s) ) return FAIL_EVENT;

	  _save_();
	  
	  // then change the static domain
	  setdifference = DOMAIN_EVENT;
	  domain.values.setminusWith(s);
	  domain.size = domain.values.size();
	  if(s.member(domain.min)) { setdifference |= LB_EVENT; domain.min = domain.values.next(domain.min); }
	  if(s.member(domain.max)) { setdifference |= UB_EVENT; domain.max = domain.values.prev(domain.max); }
	  if(domain.min == domain.max) setdifference |= VALUE_EVENT;
	} else {
	  std::cerr << "not supported" << std::endl;
	  exit(0);
	}
      } 

      triggerEvent(setdifference);
      return setdifference; // do nothing more for BOOL or CONSTANT or BITSET or RANGE
    }

    /// Remove all values in the interval [l..u]
    inline Event removeRange(const int low, const int up) {
      Event removal = DOMAIN_EVENT;

      if(type != VIRTUAL_VAR) {
	if(low <= domain.min) return setMin(up +1);
	if(up  >= domain.max) return setMax(low-1);
	// first check if we can abort early
	if(domain.values.table) {
	  if(!domain.values.intersect(low, up)) return NO_EVENT;

	  _save_();

	  domain.values.removeInterval(low, up);
	  domain.size = domain.values.size();
	}
      }
      
      triggerEvent(removal);
      return removal; // do nothing more for BOOL or CONSTANT or BITSET or RANGE
    }
    //@}


    /*!@name Virtual Accessors and Iterators*/
    //@{
    /// Returns the minimum value that could belong to the domain
    inline int minCapacity() const { return (trail_.size ? trail_[0] : 0); }
    /// Returns 1 + the maximum value that could belong to the domain
    inline int maxCapacity() const { return (trail_.size ? trail_[1] : 1); }
    //@}   

    virtual void full_print() {
      std::cout << this << " " 
		<< domain << " ";
      ConstraintNode nd = constraints.first(_value_);
      while( constraints.next(nd) )
	std::cout << "(" << (nd.elt.constraint) << ") ";
    }

    inline void save() {
      int i = domain.values.pos_words, j = domain.values.neg_words;

      trail_.add(domain.min);
      trail_.add(domain.max);
      trail_.add(domain.size);

      //      if(values_trail) 
      if(delta_) 
	while( i --> j ) {
// 	  std::cout << this << std::endl;
// 	  std::cout << i << std::endl;
// 	  std::cout << values_trail[i] << std::endl;

// 	  if(//domain.values.table[i] && 
// 	     values_trail[i].back() != domain.values.table[i]) {
// 	    values_trail[i].add(domain.values.table[i]);
// 	    level_trail[i].add(solver->level);
// 	  }
	  if(*(delta_[i]) != domain.values.table[i]) {
	    *(++delta_[i]) = domain.values.table[i];
	    *(++level_[i]) = solver->level;
	  }
	}

      trail_.add(solver->level);
    }

    virtual void restore() {

      // backtrack from solver->level to solver->level-1;
      // the domain has changed at solver->level
      int i = domain.values.pos_words, j = domain.values.neg_words;

      trail_.pop();
      trail_.pop(domain.size);
      trail_.pop(domain.max);
      trail_.pop(domain.min);

//       if(values_trail) 
// 	while( i --> j ) {
// 	  //if(values_trail[i].back() != domain.values.table[i]) 
// 	  domain.values.table[i] = values_trail[i].back();
// 	  if(level_trail[i].back() == solver->level) {
// 	    level_trail[i].pop();
// 	    values_trail[i].pop();
// 	    //values_trail[i].pop(domain.values.table[i]);
// 	  }
// 	}

      //if(values_trail) 
      if(delta_) 
	while( i --> j ) {
	  domain.values.table[i] = *(delta_[i]);
	  if(*(level_[i]) == solver->level) {
	    --level_[i];
	    --delta_[i];
	  }
	}

    }

//     virtual std::string getString() const {  
//       std::string return_str = 
// 	"x_"+toString(id); 
//       return return_str;
//     }

    virtual std::ostream& display(std::ostream& os) const {
      os << "x_" << id;
      return os;
    }
    
    virtual void debug_print() const {
      std::cout << this << " " 
		<< domain << " "
		<< trail_ << " " << std::endl;

      for(int i=domain.values.neg_words; i<domain.values.pos_words; ++i) {
	std::cout << "word[" << i << "] " ;
	for(WORD_TYPE* it = delta_abs[i]; it <= delta_[i]; ++it) {
	  printBitset(*it, i, std::cout);
	}
	std::cout << std::endl;
      }
    }

  };



  typedef VariableInt< unsigned long long int > VariableInt64;
  typedef VariableInt< unsigned int > VariableInt32;


#ifdef _BIT64
  
  typedef VariableInt64 IntegerVar;
  typedef VariableInt64* IntVar;

#else

  typedef VariableInt32 IntegerVar;
  typedef VariableInt32* IntVar;

#endif



  template < class WORD_TYPE, int NUM_WORDS >
  class VariableWord : public VariableInt< WORD_TYPE > {

  public: 

    WORD_TYPE word[NUM_WORDS];


    /*!@name Constructors*/
    //@{
    VariableWord(const int lb, const int ub, const int TYPE=BITSET_VAR) 
      : VariableInt< WORD_TYPE >(lb, ub, TYPE) {}

    virtual void initialise(const int lb, const int ub, const int var_type) {
      int nw = (lb >> BitSet::EXP);
      int pw = (ub >> BitSet::EXP);
      assert((pw-nw) <= NUM_WORDS);
      VariableInt< WORD_TYPE >::domain.initialise(lb, ub, false);

      VariableInt< WORD_TYPE >::domain.values.neg_words = nw;
      VariableInt< WORD_TYPE >::domain.values.pos_words = pw;
      VariableInt< WORD_TYPE >::domain.values.table = ((&(word[0]))-nw);

//       VariableInt< WORD_TYPE >::trail_.initialise(0,8);
//       VariableInt< WORD_TYPE >::trail_.add(VariableInt< WORD_TYPE >::domain.min);
//       VariableInt< WORD_TYPE >::trail_.add(VariableInt< WORD_TYPE >::domain.max);
//       VariableInt< WORD_TYPE >::trail_.add(VariableInt< WORD_TYPE >::domain.size);
//       VariableInt< WORD_TYPE >::trail_.add(-1);
//       VariableInt< WORD_TYPE >::values_trail = new Vector< WORD_TYPE >[pw-nw];
//       VariableInt< WORD_TYPE >::values_trail -= nw;
//       for(int i=nw; i<pw; ++i) {
// 	VariableInt< WORD_TYPE >::values_trail[i].initialise(0, VariableInt< WORD_TYPE >::domain.values.size(i));
// 	VariableInt< WORD_TYPE >::values_trail[i].add(VariableInt< WORD_TYPE >::domain.values.table[i]);
//       }
      VariableInt< WORD_TYPE >::initialise_trail();
    }

    virtual ~VariableWord() {
      VariableInt< WORD_TYPE >::domain.values.table = NULL;
      VariableInt< WORD_TYPE >::domain.values.neg_words = 0;
    }

  };


  class VarSelectorDomain 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorDomain() {d_ = NOVAL;}
    //@}

    /**@name Parameters*/
    //@{ 
    int d_;
    //@}  

    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorDomain& x ) const { return d_ < x.d_; }
    inline void operator=( VarSelectorDomain& x ) { d_ = x.d_; }
    inline void operator=( IntVar x ) { d_ = x->domain.size; }
    //@}  
  };

  class VarSelectorDomainMinVal 
  {
  public: 
    
    /**@name Constructors*/
    //@{
    VarSelectorDomainMinVal() {d_ = NOVAL; v_ = NOVAL;}
    //@}
    
    /**@name Parameters*/
    //@{ 
    int d_;
    int v_;
    //@}  

    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorDomainMinVal& x ) const { return d_ < x.d_  || ( d_ <= x.d_ && v_ < x.v_) ;}
    inline void operator=( VarSelectorDomainMinVal& x ) { d_ = x.d_; v_ = x.v_; }
    inline void operator=( IntVar x ) { d_ = x->domain.size; v_ = x->domain.min; }
    //@}  
  };

  class VarSelectorDegree 
  {
  public: 
    
    /**@name Constructors*/
    //@{
    VarSelectorDegree() {d_ = 0;}
    //@}
    
    /**@name Parameters*/
    //@{ 
    int d_;
    //@}    
    
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorDegree& x ) const { return d_ > x.d_; }
    inline void operator=( VarSelectorDegree& x ) { d_ = x.d_; }
    inline void operator=( IntVar x ) { d_ = x->constraints.degree; }
    //@}  
  };

  class VarSelectorDomainDegree 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorDomainDegree() {dom_ = NOVAL; deg_ = 0;}
    //@}
    
    /**@name Parameters*/
    //@{ 
    int deg_;
    int dom_;
    //@}  
    
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorDomainDegree& x ) const { 
      return dom_ < x.dom_  || ( dom_ <= x.dom_ && deg_ > x.deg_) ; }   
    inline void operator=( VarSelectorDomainDegree& x ) { deg_ = x.deg_; dom_ = x.dom_; }
    inline void operator=( IntVar x ) { deg_ = x->constraints.degree; dom_ = x->domain.size; }
    //@}  
  };

  class VarSelectorDomainOverDegree 
  {
  public: 
    
    /**@name Constructors*/
    //@{
    VarSelectorDomainOverDegree() {dom_ = NOVAL; deg_ = 0;}
    //@}
    
    /**@name Parameters*/
    //@{ 
    int deg_;
    int dom_;
    //@}  
    
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorDomainOverDegree& x ) const { return dom_ * x.deg_ < x.dom_ * deg_ ; }
    inline void operator=( VarSelectorDomainOverDegree& x ) { deg_ = x.deg_; dom_ = x.dom_; }
    inline void operator=( IntVar x ) { deg_ = x->constraints.degree; dom_ = x->domain.size; }
    //@}  
  };

  class VarSelectorWeight 
  {
  public: 
    
    /**@name Constructors*/
    //@{
    VarSelectorWeight() {w_ = 0.0;}
    //@}
    
    /**@name Parameters*/
    //@{ 
    float w_;
    //@}  
    
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorWeight& x ) const { return w_ > x.w_; }
    inline void operator=( VarSelectorWeight& x ) { w_ = x.w_; }
    inline void operator=( IntVar x ) { w_ = x->weight; }
    //@}  
  };

  class VarSelectorDomainOverWeight 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorDomainOverWeight() {dom_ = NOVAL; wgt_ = 0;}
    //@}
    
    /**@name Parameters*/
    //@{ 
    float wgt_;
    int dom_;
    //@}  
    
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorDomainOverWeight& x ) const { return dom_ * x.wgt_ < x.dom_ * wgt_ ; }
    inline void operator=( VarSelectorDomainOverWeight& x ) { wgt_ = x.wgt_; dom_ = x.dom_; }
    inline void operator=( IntVar x ) { wgt_ = x->weight; dom_ = x->domain.size; }
    //@}
  };


//   //class VariableInt;
//   template < class WORD_TYPE >
//   std::string toString(const VariableInt< WORD_TYPE >* x) {
//     return x->getString();
//   }



  template < class WORD_TYPE >
  std::ostream& operator<< (std::ostream& os, const VariableInt< WORD_TYPE >* x) {
    //return os << x->getString();
    return x->display(os);
  }

  template < class WORD_TYPE >
  std::ostream& operator<< (std::ostream& os, const VariableInt< WORD_TYPE >& x) {
    //return os << x.getString();
    return x.display(os);
  }



}

#endif // __VARIABLE_HPP
