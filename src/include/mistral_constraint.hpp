
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


/** \file mistral_constraint.hpp
    \brief Header for the Constraint library.
*/


#include <string>
#include <vector>

#include <mistral_global.hpp>
#include <mistral_backtrack.hpp>
#include <mistral_structure.hpp>
//#include <mistral_sat.hpp>


#ifndef __CONSTRAINT_HPP
#define __CONSTRAINT_HPP

#define FILTER( var, method )				    \
  event_type[(var)] = scope[(var)].method ;		    \
  if(event_type[(var)] == FAIL_EVENT) wiped = FAILURE(var); \
  else if(event_type[(var)] != NO_EVENT && !changes.contain(var)) changes.add(var);


namespace Mistral {

 
  //class ConstraintList;
  class Solver;
  class Variable;
  class Decision;
  /**********************************************
   * Constraint
   **********************************************/
  /*! \class Constraint
    \brief Representation of Constraints.

    The constrained variables are stored in the array _scope[]. 
    The method propagate() is called during preprocessing.
    Generic propagate methods, using AC3
    with residual supports are implemented in
    the class Constraint.
  */
  class Constraint {

    //protected:
  public:

    /*!@name Parameters*/
    //@{
    /// The trail, used for backtracking
    Vector< int > trail_;

    /// 
    int   *solution;
    int ***supports;
    //@}
    

  public:

    /*!@name Parameters*/
    //@{
    ///
    Solver *solver;

    /// Whether the propagation is delayed
    int priority;
    /// An unique id
    int id;

    /// stress >= active.size means that the constraint is entailed
    /// stress <  active.size means that the constraint is active
    /// If the propagator enforces at least FC then stress should be >= 1
    unsigned int stress;

    /// whether the constraint is posted/relaxed
    bool is_posted;

    /// The indices of unassigned variables
    IntStack active;

    /// We use two lists so that the active constraint can add events to its own 
    /// list (events) without changing the one used during propagation (changes)
    IntStack changes; // this is the list that is accessible from a propagator
    IntStack events; // this is the list that collects the events
    /// The type of event for each modified variable
    Event *event_type;

    // used to post/relax the constraints on the right triggers
    int *trigger;
    int *self;

    /// The constrained variables.
    Vector< Variable > scope;
    //@}


    /*!@name Constructors*/
    //@{
    /// The _scope is build by copying an array "scp" containing "l" variables
    virtual void mark_domain();
    virtual void initialise();
    virtual Constraint *clone() = 0;
    Constraint();
    Constraint(Vector< Variable >& scp);
    Constraint(std::vector< Variable >& scp);
    virtual ~Constraint();
    
    void add(Variable X);

    /// An idempotent constraint should not be called on events triggered by itself.
    /// To forbid that, the lists 'events' and 'changes' are merged
    inline void set_idempotent(const bool idp) {
      if(idp) {
	events.size = 0;
	events.capacity = scope.size;
	events.list_ = changes.list_;
	events.index_ = changes.index_;
	events.start_ = NULL;
      } else {
	events.initialise(0, scope.size-1, false);
      }
    }

    // called before propagation, the events strored in the list 'events'
    // are copied onto the list 'changes'.
    inline Constraint* freeze() {
      // if the constraint is idempotent, the two lists point to the same
      // elements, we just set the size up to what it should be.
      changes.size = events.size;

      // otherwise,
      // before each propagation, the lists events and changes are swapped 
      if(changes.list_ != events.list_) {

	// make the changes list points to the events list, and clear the other
	events.size = 0;
	
	int *iaux = events.list_;
	events.list_ = changes.list_;
	changes.list_ = iaux;
	
	unsigned int *uaux = events.index_;
	events.index_ = changes.index_;
	changes.index_ = uaux;
	
	return NULL;
      }
      return this;
    }

    inline void defrost() {
      if(is_posted && active.size <= stress) {

#ifdef _DEBUG_CGRAPH
	std::cout << " ---> relax " ;
	display(std::cout);
	std::cout << std::endl;
#endif

	relax();
      }
      if(changes.list_ == events.list_) 
	// if idempotent, clear the events list
	events.size = 0;      
    }

    inline void trigger_on(const int t, const int x) {
      trigger[x] = t;
    }


    //Constraint* notify_assignment(const int var, const int level) ;
    inline Constraint* notify_assignment(const int var, const int level) {

      //   std::cout << "remove " << scope[var] << " from " ;
      //   //<< this //<< std::endl;
      //   display(std::cout);
      //   std::cout << " " << active << " -> "; 
      
      Constraint *r = NULL;

      if(active.contain(var)) {
	if(trail_.back() != level) {
	  trail_.add(active.size);
	  trail_.add(level);
	  r = this;
	}
      
	active.remove(var);
      }
      //   std::cout << active << std::endl;
      
      return r;
    }

    void post(Solver *s);
    void relax();
    void relax_from(const int var);   



    // void monitor(const int k, PropagationOutcome& wiped) {
    //   if(event_type[k] == FAIL_EVENT) wiped = FAILURE(k);
    // }


    //inline
    void restore();//  {

    //   if(active.size <= stress) {
    // 	solver->notify_post(this);
    //   }

    //   trail_.pop();
    //   active.size = trail_.pop();
    // }


    /*!@name Propagators*/
    //@{
    /*!
     * This methods returns 0 if the constraint is satisfied
     * by the combination and any positive value otherwise.  
     */
    virtual int check(const int*) const = 0;
    /*!
     *  This method is called when the domain of at least one   
     *  variable in _scope has been modified. The list of 
     *  changed variables is in 'modified', and the type
     *  of event in 'event_type'.
     *  Some domains in _scope may be reduced, and NULL is
     *  returned on success (there is still at least one consistent
     *  assignment) or a ptr to the wiped-out variable otherwise. 
     */
    PropagationOutcome checker_propagate() { // std::cout << "checker propagate\n";
      return Constraint::propagate(); 
    }
    PropagationOutcome bound_propagate();
    PropagationOutcome bound_wordy_propagate();
    virtual PropagationOutcome propagate(); // { return NULL; }
    virtual PropagationOutcome rewrite();
    virtual void consolidate();

    virtual int get_backtrack_level();// {return solver->level-1;}
    virtual Decision get_decision();// { return solver->decisions.back(0); }
    /*!
     *  Check if the cached support is valid.
     *  Otherwise initialize an iterator on supports
     */
    bool first_support(const int, const int);
    /*!
     *  Find the first valid support. 
     *  Return true if a support has been found, 
     *  or false if no such support exists.
     */
    bool find_support(const int, const int);
    bool find_bound_support(const int, const int);
    //@}


    /*!@name Miscellanous methods*/
    //@{
    /// Print the constraint
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "c"; }
    //@}
  };



  /**********************************************
   * == Constraint
   **********************************************/ 
  /*! \class ConstraintEqual
    \brief  Binary equal Constraint (x0 == x1).
  */
  class ConstraintEqual : public Constraint {

  public:  
    /**@name Constructors*/
    //@{
    ConstraintEqual() : Constraint() {}
    ConstraintEqual(Vector< Variable >& scp) 
      : Constraint(scp) {}
    ConstraintEqual(std::vector< Variable >& scp) 
      : Constraint(scp) {}
    virtual Constraint *clone() { return new ConstraintEqual(scope); }
    virtual void initialise();
    virtual ~ConstraintEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check(const int* sol) const { return (sol[0] != sol[1]); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "not_equal"; }
    //@}
  };

  /**********************************************
   * =/= Constraint
   **********************************************/ 
  /*! \class ConstraintNotEqual
    \brief  Binary not-equal Constraint (x0 =/= x1).
  */
  class ConstraintNotEqual : public Constraint {

  public:  
    /**@name Constructors*/
    //@{
    ConstraintNotEqual() : Constraint() {}
    ConstraintNotEqual(Vector< Variable >& scp) 
      : Constraint(scp) {}
    ConstraintNotEqual(std::vector< Variable >& scp) 
      : Constraint(scp) {}
    virtual Constraint *clone() { return new ConstraintNotEqual(scope); }
    virtual void initialise();
    virtual void mark_domain();
    virtual ~ConstraintNotEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check(const int* sol) const { return (sol[0] == sol[1]); }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "not_equal"; }
    //@}
  };

  /**********************************************
   * <= Constraint
   **********************************************/ 
  /*! \class ConstraintLess
    \brief  Binary Less Than Constraint (x0 + k <= x1).
  */
  class ConstraintLess : public Constraint {

  public: 
    /**@name Parameters*/
    //@{  
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintLess() : Constraint() {}
    ConstraintLess(const int ofs=0) 
      : Constraint() { offset = ofs; }
    ConstraintLess(Vector< Variable >& scp, const int ofs=0) 
      : Constraint(scp) { offset = ofs; }
    ConstraintLess(std::vector< Variable >& scp, const int ofs=0) 
      : Constraint(scp) { offset = ofs; }
    virtual Constraint *clone() { return new ConstraintLess(scope, offset); }
    virtual void initialise();
    virtual ~ConstraintLess() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return (sol[0]+offset > sol[1]); }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "less_than"; }
    //@}
  };


  /**********************************************
   * Disjunctive Constraint (implemented with two precedence constraints)
   **********************************************/ 
  /*! \class ConstraintDisjunctive
    \brief  Binary Disjunctive Constraint (x0 + p0 <= x1 || x1 + p1 <= x0).
  */
  class ConstraintDisjunctive : public Constraint {

  public: 
    /**@name Parameters*/
    //@{
    int processing_time[2];
    Constraint *precedence[2];
    //@}

    /**@name Constructors*/
    //@{
    ConstraintDisjunctive() : Constraint() {}
    ConstraintDisjunctive(Vector< Variable >& scp, const int p0, const int p1); 
    ConstraintDisjunctive(std::vector< Variable >& scp, const int p0, const int p1); 
//       : Constraint(scp) { 
//       processing_time[0] = p0; 
//       processing_time[1] = p1;

//       precedence[0] = new ConstraintLess(scope, processing_time[0]);
//       precedence[1] = new ConstraintLess(processing_time[1]);
//       precedence[1]->add(scope[1]);
//       precedence[1]->add(scope[0]);
//     }
    virtual Constraint *clone() { return new ConstraintDisjunctive(scope, processing_time[0], processing_time[1]); }
    virtual void initialise();
    virtual ~ConstraintDisjunctive() {}
    //@}

    /**@name Solving*/
    //@{
    void decide(const int choice);
    virtual int check( const int* sol ) const { 
      return ((sol[0]+processing_time[0] > sol[1])
	      &&
	      (sol[1]+processing_time[1] > sol[0])); 
    }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    virtual void consolidate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "disjunctive"; }
    //@}
  };


  /**********************************************
   * <= Predicate
   **********************************************/ 
  /*! \class PredicateLess

    \brief  Truth value of a precedence ((x0 + k <= x1) <-> y)
  */
  class PredicateLess : public ConstraintLess {

  public: 

    /**@name Constructors*/
    //@{
    PredicateLess() : ConstraintLess(0) {}
    PredicateLess(Vector< Variable >& scp, const int ofs=0) 
      : ConstraintLess(scp, ofs) {}
    PredicateLess(std::vector< Variable >& scp, const int ofs=0) 
      : ConstraintLess(scp, ofs) {}
    virtual Constraint *clone() { return new PredicateLess(scope, offset); }
    virtual void initialise();
    virtual ~PredicateLess() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0]+offset <= sol[1]) != sol[2]); 
    }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "less_than"; }
    //@}
  };


  /**********************************************
   * <=K Predicate
   **********************************************/ 
  /*! \class PredicateUpperBound

    \brief  Truth value of a precedence ((x0 <= k) <-> y)
  */
  class PredicateUpperBound : public Constraint {

  public: 

    int bound;

    /**@name Constructors*/
    //@{
    PredicateUpperBound() : Constraint() {}
    PredicateUpperBound(Vector< Variable >& scp, const int b=0) 
      : Constraint(scp) { bound = b; }
    PredicateUpperBound(std::vector< Variable >& scp, const int b=0) 
      : Constraint(scp) { bound = b; }
    virtual Constraint *clone() { return new PredicateUpperBound(scope, bound); }
    virtual void initialise();
    virtual ~PredicateUpperBound() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0] <= bound) != sol[1]); 
    }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "upper_bound"; }
    //@}
  };


  /**********************************************
   * >=K Predicate
   **********************************************/ 
  /*! \class PredicateLowerBound

    \brief  Truth value of a precedence ((x0 >= k) <-> y)
  */
  class PredicateLowerBound : public Constraint {

  public: 

    int bound;

    /**@name Constructors*/
    //@{
    PredicateLowerBound() : Constraint() {}
    PredicateLowerBound(Vector< Variable >& scp, const int b=0) 
      : Constraint(scp) { bound = b; }
    PredicateLowerBound(std::vector< Variable >& scp, const int b=0) 
      : Constraint(scp) { bound = b; }
    virtual Constraint *clone() { return new PredicateLowerBound(scope, bound); }
    virtual void initialise();
    virtual ~PredicateLowerBound() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0] >= bound) != sol[1]); 
    }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "lower_bound"; }
    //@}
  };


  /**********************************************
   * Equality Predicate
   **********************************************/
  /*! \class PredicateEqual
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateEqual : public Constraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    //@}

    /**@name Constructors*/
    //@{
    PredicateEqual() : Constraint() {}
    PredicateEqual(Vector< Variable >& scp, const int sp=1) 
      : Constraint(scp) { spin = sp; }
    PredicateEqual(std::vector< Variable >& scp, const int sp=1) 
      : Constraint(scp) { spin = sp; }
    virtual Constraint *clone() { return new PredicateEqual(scope, spin); }
    virtual void initialise();
    virtual void mark_domain();
    virtual ~PredicateEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((sol[0] == sol[1]) == (sol[2] ^ spin)); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "are_equal"; }
    //@}
  };


  /**********************************************
   * Interval Membership Predicate
   **********************************************/
  /*! \class PredicateIntervalMember
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateIntervalMember : public Constraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    int lower_bound;
    int upper_bound;
    //@}

    /**@name Constructors*/
    //@{
    PredicateIntervalMember() : Constraint() {}
    PredicateIntervalMember(Vector< Variable >& scp, const int lb=-INFTY, const int ub=+INFTY, const int sp=1) 
      : Constraint(scp) { spin = sp; lower_bound=lb; upper_bound=ub; }
    PredicateIntervalMember(std::vector< Variable >& scp, const int lb, const int ub, const int sp=1) 
      : Constraint(scp) { spin = sp; lower_bound=lb; upper_bound=ub; }
    virtual Constraint *clone() { return new PredicateIntervalMember(scope, lower_bound, upper_bound, spin); }
    virtual void initialise();
    virtual void mark_domain();
    virtual ~PredicateIntervalMember() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return(((sol[0]>=lower_bound) && (sol[0]<=upper_bound)) 
						       == (sol[1] ^ spin)); }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "are_equal"; }
    //@}
  };


  /**********************************************
   * Set Membership Predicate
   **********************************************/
  /*! \class PredicateSetMember
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateSetMember : public Constraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    BitSet values;
    BitSet non_values;
    //@}

    /**@name Constructors*/
    //@{
    PredicateSetMember() : Constraint() {}
    PredicateSetMember(Vector< Variable >& scp, const int sp=1) 
      : Constraint(scp) { spin = sp; }
    PredicateSetMember(Vector< Variable >& scp, const BitSet& vals, const int sp=1) 
      : Constraint(scp) { spin = sp; values=vals; }
    PredicateSetMember(std::vector< Variable >& scp, const BitSet& vals, const int sp=1) 
      : Constraint(scp) { spin = sp; values=vals; }
    virtual Constraint *clone() { return new PredicateSetMember(scope, values, spin); }
    virtual void initialise();
    virtual void mark_domain();
    virtual ~PredicateSetMember() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((values.contain(sol[0]) == (sol[1] ^ spin))); }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "are_equal"; }
    //@}
  };


  /**********************************************
   * ConstantEquality Predicate
   **********************************************/
  /*! \class PredicateConstantEqual
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateConstantEqual : public Constraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    int value;
    //@}

    /**@name Constructors*/
    //@{
    PredicateConstantEqual() : Constraint() {}
    PredicateConstantEqual(Vector< Variable >& scp, const int val=0, const int sp=1) 
      : Constraint(scp) { value = val; spin = sp; }
    PredicateConstantEqual(std::vector< Variable >& scp, const int val=0, const int sp=1) 
      : Constraint(scp) { value = val; spin = sp; }
    virtual Constraint *clone() { return new PredicateConstantEqual(scope, value, spin); }
    virtual void initialise();
    virtual void mark_domain();
    virtual ~PredicateConstantEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((sol[0] == value) == (sol[1] ^ spin)); }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "equal_to"; }
    //@}
  };


  /**********************************************
   * Offset Predicate
   **********************************************/
  /*! \class PredicateOffset
    \brief  Truth value of a conjunction (!x0 <-> y)
  */
  class PredicateOffset : public Constraint
  {

  public:

    /**@name Parameters*/
    //@{  
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    PredicateOffset() : Constraint() {}
    PredicateOffset(Vector< Variable >& scp, const int ofs) 
      : Constraint(scp) { offset=ofs; }
    PredicateOffset(std::vector< Variable >& scp, const int ofs) 
      : Constraint(scp) { offset=ofs; }
    virtual Constraint *clone() { return new PredicateOffset(scope, offset); }
    virtual void initialise();
    virtual ~PredicateOffset() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0]+offset) != sol[1]);
    }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "offset"; }
    //@}
  };


  /**********************************************
   * Factor Predicate
   **********************************************/
  /*! \class PredicateFactor
    \brief  Truth value of a conjunction (!x0 <-> y)
  */
  class PredicateFactor : public Constraint
  {

  public:

    /**@name Parameters*/
    //@{  
    int factor;
    //@}

    /**@name Constructors*/
    //@{
    PredicateFactor() : Constraint() {}
    PredicateFactor(Vector< Variable >& scp, const int fct) 
      : Constraint(scp) { factor=fct; }
    PredicateFactor(std::vector< Variable >& scp, const int fct) 
      : Constraint(scp) { factor=fct; }
    virtual Constraint *clone() { return new PredicateFactor(scope, factor); }
    virtual void initialise();
    virtual ~PredicateFactor() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0]*factor) != sol[1]);
    }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "factor"; }
    //@}
  };


  /**********************************************
   * Not Predicate
   **********************************************/
  /*! \class PredicateNot
    \brief  Truth value of a conjunction (!x0 <-> y)
  */
  class PredicateNot : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateNot() : Constraint() {}
    PredicateNot(Vector< Variable >& scp) 
      : Constraint(scp) {}
    PredicateNot(std::vector< Variable >& scp) 
      : Constraint(scp) {}
    virtual Constraint *clone() { return new PredicateNot(scope); }
    virtual void initialise();
    virtual ~PredicateNot() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0]>0) == sol[1]);
    }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "not"; }
    //@}
  };


  /**********************************************
   * And Predicate
   **********************************************/
  /*! \class PredicateAnd
    \brief  Truth value of a conjunction ((x0 and x1) <-> y)
  */
  class PredicateAnd : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateAnd() : Constraint() {}
    PredicateAnd(Vector< Variable >& scp) 
      : Constraint(scp) {}
    PredicateAnd(std::vector< Variable >& scp) 
      : Constraint(scp) {}
    virtual Constraint *clone() { return new PredicateAnd(scope); }
    virtual void initialise();
    virtual ~PredicateAnd() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0] && sol[1]) != (sol[2])); 
    }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "and"; }
    //@}
  };

  /**********************************************
   * And Constraint
   **********************************************/
  /*! \class ConstraintAnd
    \brief  Conjunction (x0 and x1)
  */
  class ConstraintAnd : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    ConstraintAnd() : Constraint() {}
    ConstraintAnd(Vector< Variable >& scp) 
      : Constraint(scp) {}
    ConstraintAnd(std::vector< Variable >& scp) 
      : Constraint(scp) {}
    virtual Constraint *clone() { return new ConstraintAnd(scope); }
    virtual void initialise();
    virtual ~ConstraintAnd() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return(!(sol[0] && sol[1])); 
    }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "and"; }
    //@}
  };


  /**********************************************
   * Or Predicate
   **********************************************/
  /*! \class PredicateOr
    \brief  Truth value of a disjunction ((x0 or x1) <-> y)
  */
  class PredicateOr : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateOr() : Constraint() {}
    PredicateOr(Vector< Variable >& scp) 
      : Constraint(scp) {}
    PredicateOr(std::vector< Variable >& scp) 
      : Constraint(scp) {}
    virtual Constraint *clone() { return new PredicateOr(scope); }
    virtual void initialise();
    virtual ~PredicateOr() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0] || sol[1]) != (sol[2])); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "and"; }
    //@}
  };

  /**********************************************
   * Or Constraint
   **********************************************/
  /*! \class ConstraintOr
    \brief  Disjunction (x0 or x1)
  */
  class ConstraintOr : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    ConstraintOr(Vector< Variable >& scp) 
      : Constraint(scp) { }
    ConstraintOr(std::vector< Variable >& scp) 
      : Constraint(scp) { }
    virtual Constraint *clone() { return new ConstraintOr(scope); }
    virtual void initialise();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return(!(sol[0] || sol[1])); 
    }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "and"; }
    //@}
  };


  /**********************************************
   * Lex Constraint
   **********************************************/
  /*! \class ConstraintLex
    \brief  Basic bloc of a lex lt/leq constraint

    let x0 and x1 being two cells with same rank on two rows/columns 
    and b0, b1 be two Boolean variables.
    This constraint ensures that 
       - x0 =/= x1 => b1
       - b0 < b1 => x0 < x1
       - b0 <= b1 
  */
  class ConstraintLex : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    ConstraintLex() : Constraint() {}
    ConstraintLex(Vector< Variable >& scp) 
      : Constraint(scp) {}
    ConstraintLex(std::vector< Variable >& scp) 
      : Constraint(scp) {}
    virtual Constraint *clone() { return new ConstraintLex(scope); }
    virtual void initialise();
    virtual ~ConstraintLex() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return( ((!sol[2] && !sol[3]) > (sol[0] == sol[1])
	       || (sol[2] < sol[3]) > (sol[0] <  sol[1])
	       || sol[2] > sol[3])
	      ); 
    }
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "lex"; }
    //@}
  };

  // /* leading constraint in the decomposition */
  // class ConstraintLexf : public Constraint
  // {

  // public:
  //   /**@name Constructors*/
  //   //@{
  //   ConstraintLexf() : Constraint() {}
  //   ConstraintLexf(Vector< Variable >& scp) 
  //     : Constraint(scp) {}
  //   ConstraintLexf(std::vector< Variable >& scp) 
  //     : Constraint(scp) {}
  //   virtual Constraint *clone() { return new ConstraintLex(scope); }
  //   virtual void initialise();
  //   virtual ~ConstraintLexf() {}
  //   //@}

  //   /**@name Solving*/
  //   //@{
  //   virtual int check( const int* sol ) const { 
  //     return( (sol[0]>sol[1]) || (sol[2]!=(sol[0]<sol[1])) );

  //   }
  //   virtual PropagationOutcome propagate();
  //   //virtual PropagationOutcome rewrite();
  //   //@}

  //   /**@name Miscellaneous*/
  //   //@{  
  //   virtual std::ostream& display(std::ostream&) const ;
  //   virtual std::string name() const { return "lexf"; }
  //   //@}
  // };

  /**********************************************
   * Addition Predicate
   **********************************************/
  /*! \class PredicateAdd
    \brief  Binary addition predicate (x0 + x1 = y)
  */
  class PredicateAdd : public Constraint
  {
    
  public:
    /**@name Constructors*/
    //@{
    PredicateAdd() : Constraint() {}
    PredicateAdd(Vector< Variable >& scp) 
      : Constraint(scp) {}
    PredicateAdd(std::vector< Variable >& scp) 
      : Constraint(scp) {}
    virtual Constraint *clone() { return new PredicateAdd(scope); }
    virtual void initialise();
    virtual ~PredicateAdd() {}
    //@}
    
    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return (sol[2] != (sol[0]+sol[1])); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome rewrite();
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "plus"; }
    //@}
  };


  /**********************************************
   * Multiplication Predicate
   **********************************************/
  /*! \class PredicateMul
    \brief  Binary mulition predicate (x0 + x1 = y)
  */
  class PredicateMul : public Constraint
  {
    
  public:

    // class SplitDomain {
    // private:
    //   int min_neg; // minimum negative value, if there are none, then this is the minimum value
    //   int max_neg; // maximum negative value, if there are none, then this is the maximum value
    //   int min_pos; // minimum positive value, if there are none, then this is the minimum value
    //   int max_pos; // maximum positive value, if there are none, then this is the maximum value
    //   int polarity; // 1->can be =0 // 2->can be >0 // 4->can be <0

    // public:
    //   inline bool con_be_zero() { return polarity&1; }
    //   inline bool con_be_positive() { return polarity&2; }
    //   inline bool con_be_negative() { return polarity&4; }
      
    //   // set the bounds to
    //   inline bool set_pos_bounds(const int l, const int u) {
    // 	bool change = false;
    // 	if(l>u) {
    // 	  change = (max_pos)
    // 	}

    // 	  return false
    //   }

    // };



    int 
    min_pos[3],
      max_pos[3],
      min_neg[3],
      max_neg[3],
      zero[3];

    /**@name Constructors*/
    //@{
    PredicateMul() : Constraint() {}
    PredicateMul(Vector< Variable >& scp) 
      : Constraint(scp) {}
    PredicateMul(std::vector< Variable >& scp) 
      : Constraint(scp) {}
    virtual Constraint *clone() { return new PredicateMul(scope); }
    virtual void initialise();
    virtual ~PredicateMul() {}
    //@}
    
    /**@name Solving*/
    //@{
    PropagationOutcome pruneZeros(const int changedIdx);
    PropagationOutcome pruneUnary(const int last); 
    PropagationOutcome pruneBinary(const int otherIdx, const int reviseIdx, const int v);
    PropagationOutcome pruneTernary(const int reviseIdx);

    void compute_division(Variable X, Variable Y, int& lb, int& ub);
    void compute_multiplication(Variable X, Variable Y, int& lb, int& ub);
    void refine_bounds(const int evt_idx);

    PropagationOutcome revise_division(const int X, const int Y, const int Z);
    PropagationOutcome revise_multiplication(const int X, const int Y, const int Z);
    PropagationOutcome prune(const int lb_neg, 
			     const int ub_neg, 
			     const int lb_pos, 
			     const int ub_pos,
			     const bool pzero,
			     const int Z);

    virtual int check( const int* sol ) const { return (sol[2] != (sol[0]*sol[1])); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome rewrite();
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "plus"; }
    //@}
  };

  #define NE 0
  #define EQ 1
  #define GT 2
  #define LE 3

  class Literal {

  public:
    
    unsigned int __data__;
    int __value__;

    Literal() {
      __data__ = 0;
      __value__ = 0;
    }
    Literal(const int id, const int tp=EQ, const int va=1) {
      __value__ = va;
      __data__ = tp + (id*4);
    }

    virtual ~Literal() {
    }
    
    inline unsigned int get_type() const {
      // 0 -> inequality
      // 1 -> equality
      // 2 -> lower bound
      // 3 -> upper bound
      return (__data__&3);
    }

    inline unsigned int get_atom() const {
      return __data__/4;
    }

    inline int get_value() const {
      return __value__;
    }

    inline void invert() {
      __data__^=1;
    }

    bool check(const int x) const;

    bool is_true(Variable x) const;

    bool is_false(Variable x) const;

    bool apply(Variable x) const;

    std::ostream& display(std::ostream&) const ;

  };


  // /***********************************************
  //  * NogoodBase Constraint (forward checking).
  //  ***********************************************/
  // /*! \class ConstraintNogoodBase
  //   \brief   Constraint.
  // */
  // class ConstraintNogoodBase : public Constraint {

  // public:

  //   /**@name Parameters*/
  //   //@{ 
  //   // minimum values, used as an offset when accessing the base
  //   //Vector< int > minimums;
  //   int* minimums;
  //   // list of nogoods
  //   Vector< Array < Literal >* > nogood;
  //   // the watched literals data structure
  //   Vector< Vector< Array < Literal >* >* > watch_structure;
  //   //@}
    
  //   /**@name Constructors*/
  //   //@{
  //   ConstraintNogoodBase() : Constraint() {}
  //   ConstraintNogoodBase(Vector< Variable >& scp);
  //   virtual void mark_domain();
  //   virtual Constraint *clone() { return new ConstraintNogoodBase(scope); }
  //   virtual void initialise();
  //   virtual ~ConstraintNogoodBase();
    
  //   void add( Variable x );
  //   void add( Vector < Literal >& clause );
  //   //@}

  //   /**@name Solving*/
  //   //@{
  //   virtual int check( const int* sol ) const ;
  //   virtual PropagationOutcome propagate();
  //   //virtual PropagationOutcome rewrite();
  //   //@}

  //   /**@name Miscellaneous*/
  //   //@{  
  //   virtual std::ostream& display(std::ostream&) const ;
  //   virtual std::string name() const { return "nogood_base"; }
  //   //@}
    
  // };



  /***********************************************
   * AllDifferent Constraint (forward checking).
   ***********************************************/
  /*! \class ConstraintCliqueNotEqual
    \brief  Clique of NotEqual Constraint.
  */
  class ConstraintCliqueNotEqual : public Constraint {

  public:
    
    /**@name Constructors*/
    //@{
    ConstraintCliqueNotEqual() : Constraint() {}
    ConstraintCliqueNotEqual(Vector< Variable >& scp);
    ConstraintCliqueNotEqual(std::vector< Variable >& scp);
    virtual void mark_domain();
    virtual Constraint *clone() { return new ConstraintCliqueNotEqual(scope); }
    virtual void initialise();
    virtual ~ConstraintCliqueNotEqual();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    //virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "clique_ne"; }
    //@}
    
  };


  /**********************************************
   * BoolSum Equal Constraint
   **********************************************/
  //  lb <= x1 + ... + xn <= ub
  /// Constraint on the value of the sum of a set of variables.
  class ConstraintBoolSumEqual : public Constraint {
 
  public:
    /**@name Parameters*/
    //@{
    int total;
    ReversibleNum<int> min_;
    ReversibleNum<int> max_;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintBoolSumEqual() : Constraint() {}
    ConstraintBoolSumEqual(Vector< Variable >& scp, const int t);
    ConstraintBoolSumEqual(std::vector< Variable >& scp, const int t);
    virtual Constraint *clone() { return new ConstraintBoolSumEqual(scope, total); }
    virtual void initialise();
    virtual ~ConstraintBoolSumEqual();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "bool_sum"; }
    //@}
  };


  /**********************************************
   * BoolSum Interval Constraint
   **********************************************/
  //  lb <= x1 + ... + xn <= ub
  /// Constraint on the value of the sum of a set of variables.
  class ConstraintBoolSumInterval : public Constraint {
 
  public:
    /**@name Parameters*/
    //@{
    int lb;
    int ub;
    ReversibleNum<int> min_;
    ReversibleNum<int> max_;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintBoolSumInterval() : Constraint() {}
    ConstraintBoolSumInterval(Vector< Variable >& scp, const int l, const int u);
    ConstraintBoolSumInterval(std::vector< Variable >& scp, const int l, const int u);
    virtual Constraint *clone() { return new ConstraintBoolSumInterval(scope, lb, ub); }
    virtual void initialise();
    virtual ~ConstraintBoolSumInterval();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "bool_sum"; }
    //@}
  };


  /**********************************************
   * BoolSum  Predicate
   **********************************************/
  //  x1 + ... + xn-1 = xn
  /// predicate on the value of the sum of a set of variables.
  class PredicateBoolSum : public Constraint {
 
  public:
    /**@name Parameters*/
    //@{
    int lb;
    int ub;
    ReversibleNum<int> min_;
    ReversibleNum<int> max_;
    //@}

    /**@name Constructors*/
    //@{
    PredicateBoolSum() : Constraint() {}
    PredicateBoolSum(Vector< Variable >& scp);
    PredicateBoolSum(std::vector< Variable >& scp);
    PredicateBoolSum(Vector< Variable >& scp, Variable tot);
    PredicateBoolSum(std::vector< Variable >& scp, Variable tot);
    virtual Constraint *clone() { return new PredicateBoolSum(scope); }
    virtual void initialise();
    virtual ~PredicateBoolSum();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "bool_sum"; }
    //@}
  };


  /**********************************************
   * Element Predicate
   **********************************************/
  /*! \class PredicateElement
    \brief  
  */
  class PredicateElement : public Constraint {

  public:
    /**@name Parameters*/
    //@{ 
    BitSet aux_dom;
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    PredicateElement(Vector< Variable >& scp, const int o=0);
    PredicateElement(std::vector< Variable >& scp, const int o=0);
    virtual ~PredicateElement();
    virtual Constraint *clone() { return new PredicateElement(scope, offset); }
    virtual void initialise();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "sum"; }
    //@}
  };


  // /**********************************************
  //  * BoolElement Predicate
  //  **********************************************/
  // /*! \class PredicateBoolElement
  //   \brief  
  // */
  // class PredicateBoolElement : public Constraint {

  // public:
  //   /**@name Parameters*/
  //   //@{ 
  //   BitSet aux_dom;
  //   int offset;
  //   //@}

  //   /**@name Constructors*/
  //   //@{
  //   PredicateBoolElement(Vector< Variable >& scp, const int o=0);
  //   PredicateBoolElement(std::vector< Variable >& scp, const int o=0);
  //   virtual ~PredicateBoolElement();
  //   virtual Constraint *clone() { return new PredicateBoolElement(scope, offset); }
  //   virtual void initialise();
  //   //@}

  //   /**@name Solving*/
  //   //@{
  //   virtual int check( const int* sol ) const ;
  //   virtual PropagationOutcome propagate();
  //   //virtual PropagationOutcome rewrite();
  //   //@}

  //   /**@name Miscellaneous*/
  //   //@{  
  //   virtual std::ostream& display(std::ostream&) const ;
  //   virtual std::string name() const { return "sum"; }
  //   //@}
  // };


  /**********************************************
   * WeightedSum Constraint
   **********************************************/
  /*! \class ConstraintWeightedSum
    \brief  Constraint on a sum of variables (a1 * x1 + ... + an * xn = total)
  */
  class PredicateWeightedSum : public Constraint {

  public:
    /**@name Parameters*/
    //@{ 
    // Lower bound of the linear experssion
    int lower_bound;

    // Upper bound of the linear experssion
    int upper_bound;

    Vector< int > weight;
    // from index 0 to wpos (not included), the coefficients are all 1s
    int wpos;
    // from index wpos to wneg (not included), the coefficients are all >0
    int wneg;
    // from index wneg to size (not included), the coefficients are all <0

    // utils for the propagation
    int *lo_bound;
    int *up_bound;
    int *span;
    ReversibleNum<int> parity;
    ReversibleIntStack unknown_parity;
    
    //@}

    /**@name Constructors*/
    //@{
    PredicateWeightedSum(Vector< Variable >& scp,
			 const int L=0, const int U=0);
    PredicateWeightedSum(Vector< Variable >& scp,
			 Vector< int >& coefs,
			 const int L=0, const int U=0);
    PredicateWeightedSum(std::vector< Variable >& scp,
			 std::vector< int >& coefs,
			 const int L=0, const int U=0);
    // PredicateWeightedSum(Vector< Variable >& scp,
    // 			 Vector< int >& coefs,
    // 			 Variable tot);
    // PredicateWeightedSum(std::vector< Variable >& scp,
    // 			 std::vector< int >& coefs,
    // 			 Variable tot);
    virtual ~PredicateWeightedSum();
    virtual Constraint *clone() { return new PredicateWeightedSum(scope, weight); }
    virtual void initialise();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "sum"; }
    //@}
  };


  /**********************************************
   * Min Predicate
   **********************************************/
  /*! \class PredicateMin
    \brief  Predicate on the mininum of a set of variables z = (x1 , ... , xn)
  */
  class PredicateMin : public Constraint {

  public:
    /**@name Parameters*/
    //@{ 
    int last_min;
    ReversibleIntStack candidates;
    //@}

    /**@name Constructors*/
    //@{
    PredicateMin(Vector< Variable >& scp);
    virtual ~PredicateMin();
    virtual Constraint *clone() { return new PredicateMin(scope); }
    virtual void initialise();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "min"; }
    //@}
  };


  /**********************************************
   * Max Predicate
   **********************************************/
  /*! \class PredicateMax
    \brief  Predicate on the maxinum of a set of variables z = (x1 , ... , xn)
  */
  class PredicateMax : public Constraint {

  public:
    /**@name Parameters*/
    //@{ 
    int last_max;
    ReversibleIntStack candidates;
    //@}

    /**@name Constructors*/
    //@{
    PredicateMax(Vector< Variable >& scp);
    virtual ~PredicateMax();
    virtual Constraint *clone() { return new PredicateMax(scope); }
    virtual void initialise();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "max"; }
    //@}
  };


  /***********************************************
   * All Different Constraint (bounds consistency).
   ***********************************************/
  typedef struct {
    int min, max;		// start, end of interval
    int minrank, maxrank; // rank of min & max in bounds[] of an adcsp
  } Interval;
  
  /*! \class ConstraintAllDiff
    \brief  AllDifferent Constraint.

    Constraint of difference on a set of variables.
    Ensures that a set of variables are assigned distinct 
    values (only Bounds Consistency is implemented)
    The code is from Lopez-Ortiz, Quimper, Tromp and van Beek
  */
  class ConstraintAllDiff : public Constraint {
    
  private:
    /**@name Parameters*/
    //@{  
    //int *level;
    int lastLevel;
    int *t;		// tree links
    int *d;		// diffs between critical capacities
    int *h;		// hall interval links
    Interval *iv;
    Interval **minsorted;
    Interval **maxsorted;
    int *bounds;  // bounds[1..nb] hold set of min & max in the niv intervals
                  // while bounds[0] and bounds[nb+1] allow sentinels
    int nb;
    void sortit();
    int filterlower();
    int filterupper();
    void propagateValue();
    //@}

  public:
    /**@name Constructors*/
    //@{
    ConstraintAllDiff() : Constraint() {}
    ConstraintAllDiff(Vector< Variable >& scp);
    ConstraintAllDiff(std::vector< Variable >& scp);
    virtual void mark_domain();
    virtual Constraint *clone() { return new ConstraintAllDiff(scope); }
    virtual void initialise();
    virtual ~ConstraintAllDiff();
    void init();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    //virtual std::string getString() const ;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "alldiff"; }
    //@}
  };


  // /***********************************************
  //  * Global Cardinality Constraint (bounds consistency).
  //  ***********************************************/
  // typedef struct {
  //   int firstValue;
  //   int lastValue;
  //   int* sum;
  //   int* ds;
  // } partialSum;
   

  // /*! \class ConstraintGlobalCardinality
  //  \brief  Global Cardinality Constraint.

  //  User defined propagator for enforcing bounds consistency
  //  on the restricted gcc constraint when bounds on
  //  occurrences are [a_i,b_i].
  //  A value "v" must be assigned to at least
  //  minOccurrences[v - firstDomainValue] variables and at most
  //  maxOccurrences[v - firstDomainValue] variables
  //  The code is from Lopez-Ortiz, Quimper, Tromp and van Beek
  // */   
  // class ConstraintGlobalCardinality : public Constraint
  // {

  //  private:
  //   /**@name Parameters*/
  //   //@{  
  //   int lastLevel;
  //   int *t;			// tree links
  //   int *d;			// diffs between critical capacities
  //   int *h;			// hall interval links
  //   int *stableInterval;	// stable sets
  //   int *potentialStableSets;	// links elements that potentialy belong to same stable set
  //   int *newMin;
  //   Interval *iv;
  //   Interval **minsorted;
  //   Interval **maxsorted;
  //   int *bounds;  // bounds[1..nb] hold set of min & max of the n intervals
  //   // while bounds[0] and bounds[nb+1] allow sentinels
  //   int nb;
  
  //   partialSum* l; 
  //   partialSum* u;
  //   partialSum* initializePartialSum(const int firstValue, 
  // 				   int count, const int* elements);
  //   void destroyPartialSum(partialSum *p);
  //   int  sum(partialSum *p, int from, int to);
  //   int  searchValue(partialSum *p, int value);
  //   int  minValue(partialSum *p);
  //   int  maxValue(partialSum *p);
  //   int  skipNonNullElementsRight(partialSum *p, int value);
  //   int  skipNonNullElementsLeft(partialSum *p, int value);
  
  //   void sortit();
  //   int  filterLowerMax();
  //   int  filterUpperMax();
  //   int  filterLowerMin(int *tl, int *c,
  // 		      int* stableAndUnstableSets,
  // 		      int* stableInterval,
  // 		      int* potentialStableSets,
  // 		      int* newMin);
  //   int  filterUpperMin(int *tl, int *c,
  // 		      int* stableAndUnstableSets,
  // 		      int* stableInterval,
  // 		      int* newMax);
  //   //@}  

  //  public:
  //   /**@name Constructors*/
  //   //@{
  //   ConstraintGlobalCardinality( Solver *s,
  // 			       VariableInt **v,
  // 			       const int n,
  // 			       const int firstDomainValue,
  // 			       const int lastDomainValue,
  // 			       const int* minOccurrences,
  // 			       const int* maxOccurrences);
  //   ~ConstraintGlobalCardinality();
  //   //@}

  //   /**@name Solving*/
  //   //@{  
  //   inline int check( const int* ) const ;
  //   inline bool propagate();
  //   inline bool propagate(const int changedIdx, const int e); 
  //   //@}

  //   /**@name Miscellaneous*/
  //   //@{    
  //   virtual void print(std::ostream& o) const ;
  //   //@}
  // };





  /**********************************************
   * MultiQueue
   **********************************************/
  
  /*! \class MultiQueue
    \brief MultiQueue Class
  */
  template < int NUM_PRIORITY >
  class MultiQueue {

  public:

    /**@name Parameters*/
    //@{  
    Queue triggers[NUM_PRIORITY];
    int num_actives;
    BitSet _set_;
    //@}

    /**@name Constructors*/
    //@{    
    MultiQueue()
    {
      num_actives = 0;
      //active.initialise(0, NUM_PRIORITY-1, false);
    }
    void initialise(const int n)
    {
      _set_.initialise(0, n-1, BitSet::empt);
      for(int i=0; i<NUM_PRIORITY; ++i) triggers[i].initialise(n);
    }
    virtual ~MultiQueue() {}
    //@}

    /**@name Accessors*/
    //@{    
    inline bool empty() { return !num_actives; }

    inline void trigger(Constraint* cons, const int var, const Event evt)//;
    {
      int priority = cons->priority, cons_id = cons->id;
      if(_set_.fast_contain(cons_id)) {
	if(cons->events.contain(var)) {
	  cons->event_type[var] |= evt;
	} else {
	  cons->events.add(var);
	  cons->event_type[var] = evt;
	}
      } else {
	_set_.fast_add(cons_id);
	triggers[priority].add(cons_id);
	++num_actives;
	cons->events.clear();
	cons->events.add(var);
	cons->event_type[var] = evt;
      }
    }

    inline int pop()//;
    {
      int priority = NUM_PRIORITY;
      while(--priority) if(!triggers[priority].empty()) break;
      int cons = triggers[priority].pop();
      _set_.fast_remove(cons);
      --num_actives;
      return cons;
    }

    inline void clear() {
      if(!empty()) {
	for(int i=0; i<NUM_PRIORITY; ++i) triggers[i].clear();
	num_actives = 0;
	_set_.clear();
      }
    }
    //@}

  };

  std::ostream& operator<< (std::ostream& os, const Constraint& x);
  std::ostream& operator<< (std::ostream& os, const Constraint* x);
  std::ostream& operator<< (std::ostream& os, const Literal& x);
  std::ostream& operator<< (std::ostream& os, const Literal* x);


}

#endif //__CONSTRAINT_HPP




