
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


namespace Mistral {

//   class VarWrapper {
//   public:
//     int domain_type;
//     int index;


//   };
  


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
    }

//     VariableInt(VariableInt& x) : type(x.type) {
//       id = x.id;
//       weight = x.weight;
//       solver = x.solver;
//       var_list = x.var_list;
//       domain = x.domain;
//       constraints = x.constraints;

//       trail_ = x.trail_;
//       level_ = x.level_;
//       delta_ = x.delta_;
//       level_ = x.level_abs;
//       delta_ = x.delta_abs;
//     };

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
	var_list->remove(this);
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

    inline void triggerDomainEvent() {
      ConstraintNode nd;
      nd = constraints.first(_range_);
      while( constraints.next(nd) )
	solver->trigger(nd.elt.constraint, nd.elt.index, DOMAIN_EVENT);
    }

    inline void triggerRangeEvent() {
      ConstraintNode nd;
      nd = constraints.first(_range_);
      while( constraints.next(nd) )
	solver->trigger(nd.elt.constraint, nd.elt.index, RANGE_EVENT);
    }

    inline void triggerValueEvent() {
      ConstraintNode nd;
      var_list->remove(this);
      nd = constraints.first(_value_);
      while( constraints.next(nd) ) {
	solver->trigger(nd.elt.constraint, nd.elt.index, VALUE_EVENT);
	nd.elt.constraint->assign(nd.elt.index);
      }
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
    inline Event setMin(const int lo) {
      Event lower_bound = LB_EVENT;

      if(type != VIRTUAL_VAR) {
	// first check if we can abort early
	if(domain.max <  lo) return FAIL_EVENT;
	if(domain.min >= lo) return NO_EVENT;

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
	  domain.values.setMin(lo);
	  if(lo == domain.max) {
	    domain.min = lo;
	    domain.size = 1;
	    lower_bound |= VALUE_EVENT;
	  } else {
	    domain.size = domain.values.size();
	    if(domain.values.member(lo)) domain.min = lo;
	    else domain.min = domain.values.next(lo-1);
	    if(domain.size == 1) lower_bound |= VALUE_EVENT;
	  }
	} else {
	  domain.max = lo;
	  domain.size = domain.max - lo + 1;
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
    inline Event removeRange(const int lo, const int up) {
      Event removal = DOMAIN_EVENT;

      if(type != VIRTUAL_VAR) {
	if(lo <= domain.min) return setMin(up +1);
	if(up  >= domain.max) return setMax(lo-1);
	// first check if we can abort early
	if(domain.values.table) {
	  if(!domain.values.intersect(lo, up)) return NO_EVENT;

	  _save_();

	  domain.values.removeInterval(lo, up);
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

//     inline void operator=( VariableInt< WORD_TYPE > & x ) {  
    
//     }

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

      VariableInt< WORD_TYPE >::initialise_trail();
    }

    virtual ~VariableWord() {
      VariableInt< WORD_TYPE >::domain.values.table = NULL;
      VariableInt< WORD_TYPE >::domain.values.neg_words = 0;
    }

  };







#ifdef _STATIC_CAST  
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
    /// The list of all constraints whose scope includes this variable
    ConstraintList constraints;
    /// The linked solver
    Solver *solver;
    /// The stack of variables it belongs to
    Stack< IntVar > *var_list;
    /// corresponds to its rank in the variable stack
    int stack_id;
    /// unique identifier, corresponds to its rank in the variables list of the solver
    int id;
    /////////////////////////////////////

    inline void unassign() {
      ConstraintNode nd = constraints.first(_value_);
      while( constraints.next(nd) )
	nd.elt.constraint->unassign(nd.elt.index);
    }

    inline void triggerEvent(const Event evt) {

      assert(evt != NO_EVENT);

      ConstraintNode nd;
      if(is_value(evt)) {
	var_list->remove(this);
	nd = constraints.first(_value_);
	while( constraints.next(nd) ) {
	  solver->trigger(nd.elt.constraint, nd.elt.index, evt);
	  nd.elt.constraint->assign(nd.elt.index);
	}
      } else {
	if(is_range(evt)) {
	  nd = constraints.first(_range_);
	} else //if(is_domain(evt)) 
	  {
	    nd = constraints.first(_domain_);
	  }
	while( constraints.next(nd) )
	  solver->trigger(nd.elt.constraint, nd.elt.index, evt);
      }

//       if(is_value(evt)) while( constraints.next(nd) ) {
// 	  solver->trigger(nd.elt.constraint, nd.elt.index, evt);
// 	  nd.elt.constraint->assign(nd.elt.index);
// 	}
//       else 
// 	while( constraints.next(nd) )
// 	  solver->trigger(nd.elt.constraint, nd.elt.index, evt);
    }

    virtual std::ostream& display(std::ostream& os) const = 0;
  };


  class BitsetDomain {
  
  public:
  
    void initialise(const int lb, const int ub, const bool vals=true);
  
    int min;
    int max;
    int size;
    BitSet values;
  
    //void _assert_();
  
    virtual std::ostream& display(std::ostream& os) const {
      if(values.table)
	os << values;
      else
	os << "[" << min << ".." << max << "]";
      return os;
    }
  
  };


  //  std::string toString(const BitsetDomain& x);
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

    /// weight value, used for branching 
    float weight;

    void _assert_constraints_() {

//       // check active constraints 
//       if(domain.size > 1) {
// 	ConstraintNode nd = constraints.first(_value_);
// 	while( constraints.next(nd) )
// 	  {
// 	    bool is_active = false;
// 	    int i=0;
// 	    for(; !is_active && i<nd.elt.constraint->arity; ++i) {
// 	      is_active = (nd.elt.constraint->scope[i].implementation != this && nd.elt.constraint->scope[i]->domain.size > 1);
// 	    }
// 	    if(!is_active) {
// 	      std::cerr << "Unactive constraint " << (nd.elt.constraint) << " => " << (this) << std::endl;

	      
// 	      std::cerr << solver << std::endl;

// 	      exit(1);
// 	    }
// 	  }
//       }

//       // check scope 
//       Constraint *cons;
//       for(unsigned int i=0; i<constraints.data.size; ++i) {
// 	cons = ((ConstraintTrigger)(constraints.data[i])).constraint;
// 	if(cons) {
// 	  bool is_in = false;
// 	  for(int i=0; !is_in && i<cons->arity; ++i) {
// 	    is_in = (cons->scope[i].implementation == this);
// 	  }
// 	  if(!is_in) {
// 	    std::cerr << "Wrong scope" << std::endl;
// 	    exit(1);
// 	  }
// 	}
//       }
      
    }
    

    /*!@name Constructors*/
    //@{
    VariableBitset(const int lb, const int ub) {
      id = -1;
      weight = 0.0;
      solver = NULL;
      var_list = NULL;
      initialise(lb, ub, type);
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
	k = domain.values.size(i);
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
      delete [] level_abs;
      delete [] delta_;
      delete [] level_;
    }

    void _assert_() { domain._assert_(); }

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
      if( domain.values.fastMember(0) ) return 0;
      return domain.values.next( 0 );
    }
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    inline int minNegAbs() const  {
      if( domain.min > 0 ) return NOVAL;
      if( domain.max < 0 ) return domain.max;
      if( domain.values.fastMember(0) ) return 0;
      return domain.values.prev( 0 );
    } 
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int next(const int v) const { return (((domain.size != domain.max - domain.min + 1)) ? domain.values.next(v) : (v < domain.max ? v+1 : v)); }
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int prev(const int v) const { return (((domain.size != domain.max - domain.min + 1)) ? domain.values.prev(v) : (v > domain.min ? v-1 : v)); }


    /// Whether or not the Variable is currently an interval
    inline bool isRange() const { return (domain.size == domain.max - domain.min + 1); }
    /// Whether or not the Variable is bound to a ground value
    inline bool isGround() const { return (domain.size == 1); }
    /// Whether or not the Variable is bound to a given ground value
    inline bool equal(const int v) const { return (domain.min == v && domain.max == v); }
    /// Whether the value "v" is currently contained in the domain
    inline bool contain(const int v) const { return (domain.min <= v && domain.max >= v && domain.values.fastMember(v)); }

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

    /// Intersect its domain with a set s
    inline void intersectTo( BitSet& s ) const { s.intersectWith(domain.values); }
    /// Do the union of its domain with a set s
    inline void unionTo( BitSet& s ) const { s.unionWith(domain.values); }
    //@}

    /*!@name Domain handling methods*/
    //@{
    inline void _save_() {
      if(trail_.back() != solver->level) {
	solver->save(this);
	save();
      }
    }

    inline Event remove(const int v) {
      Event removal = DOMAIN_EVENT;

      if(type != VIRTUAL_VAR) {
	// first check if we can abort early
	if(!contain(v)) return NO_EVENT;
	if(domain.size == 1) return FAIL_EVENT;

	if(trail_.back() != solver->level) {
	  solver->save(this, BITSET_VAR);
	  save();
	}

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
      } 

      triggerEvent(removal);
      return removal; 
    }

    /// Remove all values but "v"
    inline Event setDomain(const int v) {
      Event setdomain = VALUE_EVENT;

      if(type != VIRTUAL_VAR) {

	// first check if we can abort early
	if(!contain(v)) return FAIL_EVENT;
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
    inline Event setMin(const int lo) {
      Event lower_bound = LB_EVENT;

      if(type != VIRTUAL_VAR) {
	// first check if we can abort early
	if(domain.max <  lo) return FAIL_EVENT;
	if(domain.min >= lo) return NO_EVENT;

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
	  domain.values.setMin(lo);
	  if(lo == domain.max) {
	    domain.min = lo;
	    domain.size = 1;
	    lower_bound |= VALUE_EVENT;
	  } else {
	    domain.size = domain.values.size();
	    if(domain.values.member(lo)) domain.min = lo;
	    else domain.min = domain.values.next(lo-1);
	    if(domain.size == 1) lower_bound |= VALUE_EVENT;
	  }
	} else {
	  domain.max = lo;
	  domain.size = domain.max - lo + 1;
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
    inline Event setDomain(VariableBitset* x) {
      BitsetDomain& xdom = x->domain;
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
    inline Event removeRange(const int lo, const int up) {
      Event removal = DOMAIN_EVENT;

      if(type != VIRTUAL_VAR) {
	if(lo <= domain.min) return setMin(up +1);
	if(up  >= domain.max) return setMax(lo-1);
	// first check if we can abort early
	if(domain.values.table) {
	  if(!domain.values.intersect(lo, up)) return NO_EVENT;

	  _save_();

	  domain.values.removeInterval(lo, up);
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




  class Variable {
  public:
//     enum domain_type { BITSET_VAR = 0xe0000000, 
// 		       LIST_VAR = 0xc0000000, 
// 		       RANGE_VAR = 0xa0000000, 
// 		       VIRTUAL_VAR = 0x60000000, 
// 		       CONST_VAR = 0x40000000 
//     };

    
  
  private:
    VariableImplementation *implementation;

  public:

    Variable(const int value);
    Variable(const int lo, const int up, const int type=BITSET_VAR);


    /*!@name Constant Accessors and Iterators*/
    //@{
    /// Returns the assigned value if it exists
    inline int get_value() const ; 
    /// Returns the minimum value in the domain
    inline int get_min() const ; 
    /// Returns the maximum value in the domain
    inline int get_max() const ; 
    /// Returns the minimum value that could belong to the domain
    inline int minCapacity() const ;
    /// Returns 1 + the maximum value that could belong to the domain
    inline int maxCapacity() const ;
    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    inline int minPosAbs() const ;
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    inline int minNegAbs() const ;
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int next(const int v) const ;
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int prev(const int v) const ;
    /// Whether or not the Variable is currently an interval
    inline bool isRange() const ;
    /// Whether or not the Variable is bound to a ground value
    inline bool isGround() const ;
    /// Whether or not the Variable is bound to a given ground value
    inline bool equal(const int v) ;
    /// Whether the value "v" is currently contained in the domain
    inline bool contain(const int v) ;

    /// Whether the domain has a nonempty intersection with the interval [l..u]
    inline bool intersect(const int lo, const int up) const ;
    /// Whether the domain is included in the interval [l..u]
    inline bool included(const int lo, const int up) const ;
    /// Whether the domain is included in the interval [l..u]
    inline bool includes(const int lo, const int up) const ;

    /// Whether the domain has a nonempty intersection with the set s 
    inline bool intersect(const BitSet& s) const ;
    /// Whether the domain is included in the set s 
    inline bool included(const BitSet& s) const ;
    /// Whether the domain is included in the set s 
    inline bool includes(const BitSet& s) const ;

    /// Whether the domain has a nonempty intersection with the Variable x
    inline bool intersect(const Variable& x) const ;
    /// Whether the domain is included in the Variable x 
    inline bool included(const Variable& x) const ;
    /// Whether the domain is included in the Variable x 
    inline bool includes(const Variable& x) const ;

    /// Intersect its domain with a set s
    inline void intersectTo( BitSet& s ) const ;
    /// Do the union of its domain with a set s
    inline void unionTo( BitSet& s ) const ;
    //@}


    /*!@name Domain changin methods*/
    //@{
    /// Remove value "v"
    inline Event remove(const int v) ;
    /// Remove all values but "v"
    inline Event setDomain(const int v) ;
    /// Remove all values strictly lower than lo
    inline Event setMin(const int lo) ;
    /// Remove all values strictly greater than up
    inline Event setMax(const int up) ;
    /// Remove all values that do not appear in the set "s"
    inline Event setDomain(const BitSet& s) ;
    /// Remove all values that do not appear in the current domain of the Variable "x"
    inline Event setDomain(Variable& x) ;
    /// Remove all values that belong to the set "s"
    inline Event removeSet(const BitSet& s) ;
    /// Remove all values in the interval [l..u]
    inline Event removeRange(const int lo, const int up) ;
    //@}


    /*!@name Virtual Accessors and Iterators*/
    //@{
    void full_print() const ;
    std::ostream& display(std::ostream& os) const ;
    void debug_print() const ;

  };


  Variable::Variable(const int value) {
    domain_type = CONST_VAR;
    implementation = (VariableImplementation*)value;
  }
  
  Variable::Variable(Solver *s, const int lo, const int up, const int type=DYN_VAR) {
    if(type == CONST_VAR || (type == DYN_VAR && lo == up)) {
      domain_type = CONST_VAR;
      implementation = (VariableImplementation*)lo;
    } else if(type == BOOL_VAR || (type == DYN_VAR && lo==0 && up==1)) {
      domain_type = (int)(s->getNextBooleanSlot());
      implementation = new VariableImplementation(s);
    } else if(type == RANGE_VAR) {
      domain_type = type;
      implementation = new VariableRange(s, lo, up);
    } else if(type == LIST_VAR) {
      domain_type = type;
      implementation = new VariableList(s, lo, up);
    } else {
      domain_type = type;
      implementation = new VariableBitset(s, lo, up);
    }
  }


  inline Event Variable::setValue( int val ) 
  {
    int nstat = (*((int*)domain_type) & val);
    if( !nstat ) return FAIL_EVENT;

    *((int*)domain_type) = nstat;

    implementation->triggerValueEvent();
    return VALUE_EVENT;  
  }

  inline int Variable::get_value() const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->get_value();
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->get_value();
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->get_value();
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->get_value();
    else if(domain__type ==   CONST_VAR) return (int)implementation;
    else  return (*((int*)domain_type)-1);
  }

  inline int Variable::get_min() const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->get_min();
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->get_min();
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->get_min();
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->get_min();
    else if(domain__type ==   CONST_VAR) return (int)implementation;
    else  return (!(*((int*)domain_type) & 1));
  }

  inline int Variable::get_max() const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->get_max();
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->get_max();
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->get_max();
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->get_max();
    else if(domain__type ==   CONST_VAR) return (int)implementation;
    else  return (*((int*)domain_type) >> 1);
  }

  inline int Variable::get_minCapacity() const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->get_minCapacity();
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->get_minCapacity();
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->get_minCapacity();
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->get_minCapacity();
    else if(domain__type ==   CONST_VAR) return (int)implementation;
    else  return 0;
  }

  inline int Variable::get_maxCapacity() const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->get_maxCapacity();
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->get_maxCapacity();
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->get_maxCapacity();
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->get_maxCapacity();
    else if(domain__type ==   CONST_VAR) return (int)implementation;
    else  return 1;
  }

  inline int Variable::get_minPosAbs() const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->get_minPosAbs();
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->get_minPosAbs();
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->get_minPosAbs();
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->get_minPosAbs();
    else if(domain__type ==   CONST_VAR) return (int)implementation;
    else  return (!(*((int*)domain_type) & 1));
  }

  inline int Variable::get_minNegAbs() const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->get_minNegAbs();
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->get_minNegAbs();
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->get_minNegAbs();
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->get_minNegAbs();
    else if(domain__type ==   CONST_VAR) return (int)implementation;
    else  return (!(*((int*)domain_type) & 1));
  }

  inline int Variable::next(const int v) const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->next(v);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->next(v);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->next(v);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->next(v);
    else if(domain__type ==   CONST_VAR) return (int)implementation;
    else  return (*((int*)domain_type) >> 1);
  }

  inline int Variable::prev(const int v) const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->prev(v);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->prev(v);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->prev(v);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->prev(v);
    else if(domain__type ==   CONST_VAR) return (int)implementation;
    else  return (!(*((int*)domain_type) & 1));
  }

  inline bool Variable::isRange() const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->isRange();
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->isRange();
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->isRange();
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->isRange();
    //else if(domain__type ==   CONST_VAR) return true;
    else return true;
  }

  inline bool Variable::isGround() const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->isGround();
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->isGround();
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->isGround();
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->isGround();
    else if(domain__type ==   CONST_VAR) return true;
    else  return (*((int*)domain_type) != 3);
  }

  inline bool Variable::equal(const int v) {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->equal(v);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->equal(v);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->equal(v);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->equal(v);
    else if(domain__type ==   CONST_VAR) return ((int)implementation == v);
    else  return (*((int*)domain_type)-1 == v);
  }

  inline bool Variable::contain(const int v) {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->contain(v);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->contain(v);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->contain(v);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->contain(v);
    else if(domain__type ==   CONST_VAR) return ((int)implementation == v);
    else  return (!(v >> 1) && (*((int*)domain_type) & (v+1)));
  }

  inline bool Variable::intersect(const int lo, const int up) const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->intersect(lo, up);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->intersect(lo, up);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->intersect(lo, up);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->intersect(lo, up);
    else if(domain__type ==   CONST_VAR) return ((int)implementation >= lo && (int)implementation <= up);
    else  return ((lo<=0 | 2*(up>0)) & *((int*)domain_type));
  }

  inline bool Variable::included(const int lo, const int up) const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->included(lo, up);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->included(lo, up);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->included(lo, up);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->included(lo, up);
    else if(domain__type ==   CONST_VAR) return ((int)implementation >= lo && (int)implementation <= up);
    else  {
      int state = *((int*)domain_type);
      return ( up >= (state >> 1) && (lo <= !(state & 1)) );
    }
  }

  inline bool Variable::includes(const int lo, const int up) const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->includes(lo, up);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->includes(lo, up);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->includes(lo, up);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->includes(lo, up);
    else if(domain__type ==   CONST_VAR) return ((int)implementation == lo && (int)implementation == up);
    else  {
      int state = *((int*)domain_type);
      return ( up <= (state >> 1) && (lo >= !(state & 1)) );
    }
  }

  inline bool Variable::intersect(const BitSet& s) const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->intersect(s);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->intersect(s);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->intersect(s);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->intersect(s);
    else if(domain__type ==   CONST_VAR) return (s.member((int)implementation));
    else  return s.intersect(*((int*)domain_type));
  }

  inline bool Variable::included(const BitSet& s) const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->included(s);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->included(s);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->included(s);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->included(s);
    else if(domain__type ==   CONST_VAR) return (s.member((int)implementation));
    else  return s.includes(*((int*)domain_type));
  }

  inline bool Variable::includes(const BitSet& s) const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->includes(s);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->includes(s);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->includes(s);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->includes(s);
    else if(domain__type ==   CONST_VAR) return (s.size() == 1 && s.member((int)implementation));
    else  return s.included(*((int*)domain_type));
  }

  inline bool Variable::intersect(const Variable& x) const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->intersect(x);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->intersect(x);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->intersect(x);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->intersect(x);
    else if(domain__type ==   CONST_VAR) return x.contain((int)implementation);
    else  {
      int state = *((int*)domain_type);
      if(state == 3) return x.intersect(0,1); 
      else if(state == 2) return x.contain(1); 
      else return x.contain(0);
    }
  }

  inline bool Variable::included(const Variable& x) const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->included(x);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->included(x);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->included(x);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->included(x);
    else if(domain__type ==   CONST_VAR) return x.contain((int)implementation);
    else  {
      int state = *((int*)domain_type);
      if(state == 3) return x.includes(0,1); 
      else if(state == 2) return x.contain(1); 
      else return x.contain(0);
    }
  }

  inline bool Variable::includes(const Variable& x) const {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->includes(x);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->includes(x);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->includes(x);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->includes(x);
    else if(domain__type ==   CONST_VAR) return (x.size() == 1 && x.contain((int)implementation));
    else  {
      int state = *((int*)domain_type);
      if(state == 3) return x.included(0,1); 
      else if(state == 2) return x.equal(1); 
      else return x.equal(0);
    }
  }

  inline void Variable::intersectTo( BitSet& s ) const {
    if     (domain__type ==  BITSET_VAR) ((VariableBitset  *)implementation)->intersectTo(s);
    else if(domain__type ==    LIST_VAR) ((VariableList    *)implementation)->intersectTo(s);
    else if(domain__type ==   RANGE_VAR) ((VariableRange   *)implementation)->intersectTo(s);
    else if(domain__type == VIRTUAL_VAR) ((VariableVirtual *)implementation)->intersectTo(s);
    else if(domain__type ==   CONST_VAR) {
      if(s.member((int)implementation)) {
	s.clear();
	s.insert((int)implementation);
      } else s.clear();
    }
    else  return s.intersectWith(*((int*)domain_type));
  }

  inline void Variable::unionTo( BitSet& s ) const {
    if     (domain__type ==  BITSET_VAR) ((VariableBitset  *)implementation)->unionTo(s);
    else if(domain__type ==    LIST_VAR) ((VariableList    *)implementation)->unionTo(s);
    else if(domain__type ==   RANGE_VAR) ((VariableRange   *)implementation)->unionTo(s);
    else if(domain__type == VIRTUAL_VAR) ((VariableVirtual *)implementation)->unionTo(s);
    else if(domain__type ==   CONST_VAR) s.insert((int)implementation);
    else  return s.unionWith(*((int*)domain_type));
  }

  inline Event Variable::remove(const int v) {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->remove(v);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->remove(v);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->remove(v);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->remove(v);
    else if(domain__type ==   CONST_VAR) return ((int)implementation == v ? FAIL_EVENT : NO_EVENT);
    else  return ((v>1 || v<0) ? NO_EVENT : setValue(2-v));
  }

  inline Event Variable::setDomain(const int v) {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->setDomain(v);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->setDomain(v);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->setDomain(v);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->setDomain(v);
    else if(domain__type ==   CONST_VAR) return ((int)implementation != v ? FAIL_EVENT : NO_EVENT);
    else  return ((v>1 || v<0) ? NO_EVENT : setValue(1+v));
  }

  inline Event Variable::setMin(const int lo) {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->setMin(lo);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->setMin(lo);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->setMin(lo);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->setMin(lo);
    else if(domain__type ==   CONST_VAR) return ((int)implementation < lo ? FAIL_EVENT : NO_EVENT);
    else  return (lo<1 ? NO_EVENT : (lo>1 ? FAIL_EVENT : setValue(2)));
  }

  inline Event Variable::setMax(const int up) {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->setMax(up);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->setMax(up);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->setMax(up);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->setMax(up);
    else if(domain__type ==   CONST_VAR) return ((int)implementation > up ? FAIL_EVENT : NO_EVENT);
    else  return (up>0 ? NO_EVENT : (up<0 ? FAIL_EVENT : setValue(1)));
  }

  inline Event Variable::setDomain(const BitSet& s) {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->setDomain(s);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->setDomain(s);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->setDomain(s);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->setDomain(s);
    else if(domain__type ==   CONST_VAR) return (s.member((int)implementation) ? NO_EVENT : FAIL_EVENT);
    else  return ((s.pos_words<1 || s.neg_words>0) ? FAIL_EVENT : ((s.table[0]&3)==3 ? NO_EVENT : setValue(s.table[0])));
  }

  inline Event Variable::setDomain(Variable& x) {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->setDomain(x);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->setDomain(x);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->setDomain(x);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->setDomain(x);
    else if(domain__type ==   CONST_VAR) return (x.contain((int)implementation) ? NO_EVENT : FAIL_EVENT);
    else  {
      if(!x.contain(0)) {
	if(!x.contain(1)) return FAIL_EVENT;
	else return setValue(2);
      } else {
	if(!x.contain(1)) return setValue(1)
	else return NO_EVENT;
      }
    }
  }

  inline Event Variable::removeSet(const BitSet& s) {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->removeSet(s);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->removeSet(s);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->removeSet(s);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->removeSet(s);
    else if(domain__type ==   CONST_VAR) return (s.member((int)implementation) ? FAIL_EVENT : NO_EVENT);
    else  return ((s.pos_words<1 || s.neg_words>0 || (s.table[0]^3)==3) ? NO_EVENT : setValue(s.table[0]^3));
  }

  inline Event Variable::removeRange(const int lo, const int up) {
    if     (domain__type ==  BITSET_VAR) return ((VariableBitset  *)implementation)->removeRange(lo, up);
    else if(domain__type ==    LIST_VAR) return ((VariableList    *)implementation)->removeRange(lo, up);
    else if(domain__type ==   RANGE_VAR) return ((VariableRange   *)implementation)->removeRange(lo, up);
    else if(domain__type == VIRTUAL_VAR) return ((VariableVirtual *)implementation)->removeRange(lo, up);
    else if(domain__type ==   CONST_VAR) return (((int)implementation < lo || (int)implementation > up) ? NO_EVENT : FAIL_EVENT);
    else  return (lo==1 ? setValue(1) : (up==0 ? setValue(2) : ((lo>1 || up<0) ? NO_EVENT : FAIL_EVENT)));
  }

  std::ostream& Variable::display(std::ostream& os) const {
    if(*((int*)domain_type)  ==   CONST_VAR) os << (int)implementation;
    else os = implementation->display(os);
    return os;
  }

#endif




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
