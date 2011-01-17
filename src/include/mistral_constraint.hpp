
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

#include <mistral_global.hpp>
#include <mistral_backtrack.hpp>
#include <mistral_structure.hpp>


#ifndef __CONSTRAINT_HPP
#define __CONSTRAINT_HPP


namespace Mistral {

 
  //class ConstraintList;
  class Solver;
  class Variable;
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

  protected:

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
    Vector<Variable> scope;
    //@}


    /*!@name Constructors*/
    //@{
    /// The _scope is build by copying an array "scp" containing "l" variables
    virtual void initialise();
    Constraint();
    Constraint(Vector< Variable >& scp);
    virtual ~Constraint();

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

    // called 
    inline Constraint* freeze() {
      changes.size = events.size;

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

      if(active.size == 1) relax();

      if(changes.list_ == events.list_) 
	// if idempotent, clear the events list
	events.size = 0;      
    }

    inline void trigger_on(const int t, const int x) {
      trigger[x] = t;
    }

    Constraint* notify_assignment(const int var, const int level);
//     inline Constraint* notify_assignment(const int var, const int level) {

//       //display(std::cout);
//       //std::cout << std::endl;
//       std::cout << "remove " << scope[var] << " from " ;
//       display(std::cout);
//       std::cout << active << " => " ;
		

//       Constraint *r = NULL;
//       if(trail_.back() != level) {
// 	trail_.add(active.size);
// 	trail_.add(level);
// 	r = this;
//       }
//       active.erase(var);

//       std::cout << active << std""endl;

//       return r;
//     }

    void post(Solver *s);
    void relax();
    
    inline void restore() {
      trail_.pop();
      active.size = trail_.pop();
    }


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
    virtual PropagationOutcome propagate() { return NULL; }
    /*!
     *  Check if the cached support is valid.
     *  Otherwise initialize an iterator on supports
     */
    bool firstSupport(const int, const int);
    /*!
     *  Find the first valid support. 
     *  Return true if a support has been found, 
     *  or false if no such support exists.
     */
    bool findSupport(const int, const int);
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
    virtual void initialise() {
      Constraint::initialise();
      trigger_on(_domain_, 0);
      trigger_on(_domain_, 1);
      set_idempotent(true);
    }
    virtual ~ConstraintEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check(const int* sol) const { return (sol[0] != sol[1]); }
    virtual PropagationOutcome propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    //virtual std::string getString() const ;
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
    virtual void initialise() {
      Constraint::initialise();
      trigger_on(_value_, 0);
      trigger_on(_value_, 1);
      set_idempotent(true);
    }
    virtual ~ConstraintNotEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check(const int* sol) const { return (sol[0] == sol[1]); }
    virtual PropagationOutcome propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    //virtual std::string getString() const ;
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
    ConstraintLess(Vector< Variable >& scp, const int ofs=0) 
      : Constraint(scp) { offset = ofs; }

    virtual void initialise() {
      Constraint::initialise();
      trigger_on(_range_, 0);
      trigger_on(_range_, 1);
      set_idempotent(true);
    }
    virtual ~ConstraintLess() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return (sol[0]+offset > sol[1]); }
    virtual PropagationOutcome propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "less_than"; }
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
    PredicateLess() : ConstraintLess() {}
    PredicateLess(Vector< Variable >& scp, const int ofs=0) 
      : ConstraintLess(scp) { offset = ofs; }

    virtual void initialise() {
      Constraint::initialise();
      trigger_on(_range_, 0);
      trigger_on(_range_, 1);
      trigger_on(_value_, 2);
      set_idempotent(false);
    }
    virtual ~PredicateLess() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0]+offset <= sol[1]) != sol[2]); 
    }
    virtual PropagationOutcome propagate();
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

    virtual void initialise() {
      Constraint::initialise();
      trigger_on(_range_, 0);
      trigger_on(_value_, 1);
      set_idempotent(false);
    }
    virtual ~PredicateUpperBound() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0] <= bound) != sol[1]); 
    }
    virtual PropagationOutcome propagate();
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

    virtual void initialise() {
      Constraint::initialise();
      trigger_on(_range_, 0);
      trigger_on(_value_, 1);
      set_idempotent(false);
    }
    virtual ~PredicateLowerBound() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0] <= bound) != sol[1]); 
    }
    virtual PropagationOutcome propagate();
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

    virtual void initialise() {
      Constraint::initialise();
      trigger_on(_domain_, 0);
      trigger_on(_domain_, 1);
      trigger_on(_value_, 2);
      set_idempotent(false);
    }
    virtual ~PredicateEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((sol[0] == sol[1]) == (sol[2] ^ spin)); }
    virtual PropagationOutcome propagate();
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
    PredicateConstantEqual(Vector< Variable >& scp, const int val, const int sp=1) 
      : Constraint(scp) { value = val; spin = sp; }

    virtual void initialise() {
      Constraint::initialise();
      trigger_on(_domain_, 0);
      trigger_on(_value_, 1);
      set_idempotent(false);
    }
    virtual ~PredicateConstantEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((sol[0] == value) == (sol[1] ^ spin)); }
    virtual PropagationOutcome propagate();
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

    virtual void initialise() {
      Constraint::initialise();
      for(int i=0; i<2; ++i)
	trigger_on(_range_, i);
      set_idempotent(false);
    }
    virtual ~PredicateOffset() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0]+offset) != sol[1]);
    }
    virtual PropagationOutcome propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "offset"; }
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

    virtual void initialise() {
      Constraint::initialise();
      for(int i=0; i<2; ++i)
	trigger_on(_value_, i);
      set_idempotent(false);
    }
    virtual ~PredicateNot() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0]>0) == sol[1]);
    }
    virtual PropagationOutcome propagate();
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

    virtual void initialise() {
      Constraint::initialise();
      for(int i=0; i<3; ++i)
	trigger_on(_value_, i);
      set_idempotent(false);
    }
    virtual ~PredicateAnd() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0] && sol[1]) != (sol[2])); 
    }
    virtual PropagationOutcome propagate();
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

    virtual void initialise() {
      Constraint::initialise();
      for(int i=0; i<2; ++i)
	trigger_on(_range_, i);
      set_idempotent(false);
    }
    virtual ~ConstraintAnd() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return(!(sol[0] && sol[1])); 
    }
    virtual PropagationOutcome propagate();
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

    virtual void initialise() {
      Constraint::initialise();
      for(int i=0; i<3; ++i)
	trigger_on(_value_, i);
      set_idempotent(false);
    }
    virtual ~PredicateOr() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0] || sol[1]) != (sol[2])); 
    }
    virtual PropagationOutcome propagate();
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

//   private: 
//     int* data[0];

  public:
    /**@name Constructors*/
    //@{
    ConstraintOr(Vector< Variable >& scp) 
      : Constraint(scp) { }

    virtual void initialise() {
      Constraint::initialise();
      for(int i=0; i<2; ++i)
	trigger_on(_range_, i);
      set_idempotent(true);
    }

 //    static ConstraintOr* ConstraintOr_new(Vector< Variable >& X);
//     //ConstraintOr() : InlineConstraint() {}
//     ConstraintOr(Vector< Variable >& scp);

//     virtual void initialise() {
//       Constraint::initialise();
//       for(int i=0; i<2; ++i)
// 	trigger_on(_range_, i);
//       set_idempotent(true);
//     }
//     virtual ~ConstraintOr() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return(!(sol[0] || sol[1])); 
    }
    virtual PropagationOutcome propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "and"; }
    //@}
  };


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

    virtual void initialise() {
      Constraint::initialise();
      trigger_on(_range_, 0);
      trigger_on(_range_, 1);
      trigger_on(_range_, 2);
      set_idempotent(true);
    }
    virtual ~PredicateAdd() {}
    //@}
    
    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return (sol[2] != (sol[0]+sol[1])); }
    virtual PropagationOutcome propagate();
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "plus"; }
    //@}
  };



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
    virtual void initialise();
    virtual ~ConstraintCliqueNotEqual();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
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
    virtual void initialise();
    virtual ~ConstraintBoolSumEqual();
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
    int *level;
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
    virtual void initialise();
    virtual ~ConstraintAllDiff();
    void init();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
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
      if(_set_.fastContain(cons_id)) {
	if(cons->events.contain(var)) {
	  cons->event_type[var] |= evt;
	} else {
	  cons->events.add(var);
	  cons->event_type[var] = evt;
	}
      } else {
	_set_.fastAdd(cons_id);
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
      _set_.fastErase(cons);
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


}

#endif //__CONSTRAINT_HPP




