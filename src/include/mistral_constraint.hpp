
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


#include <mistral_structure.hpp>
#include <mistral_global.hpp>


#ifndef __CONSTRAINT_HPP
#define __CONSTRAINT_HPP


namespace Mistral {

  template < class WORD_TYPE >
  class VariableInt;

#ifdef _BIT64
  
  typedef VariableInt< unsigned long long int >* IntVar;

#else

  typedef VariableInt< unsigned int >* IntVar;

#endif
  

  //typedef VariableInt* IntVar;

  /**********************************************
   * Constraint
   **********************************************/
  /*! \class Constraint
    \brief Representation of Constraints.

    The constrained variables are stored in the array _scope[]. 
    The method propagate() is called during preprocessing.
    The method propagate(const Event, const Vector) is
    called when computing the GAC closure. 
    Generic propagate methods, using AC3
    with residual supports are implemented in
    the class Constraint.
  */
  class Constraint {

  public:

    /*!@name Parameters*/
    //@{
    /// Whether the propagation is delayed
    int priority;
    /// An unique id
    int id;
    /// The number of constrained variables
    int arity;
    /// The indices of unassigned variables
    IntStack active;
    IntStack trail_act;
    /// Weight, used in some heuristics
    float weight;

    /// The constrained 'real' variables.
    IntVar* _scope;
    /// The variables whose domains are filtered (may be references)
    IntVar* scope;
    /// The list of modified variables' indices
    /// We use two lists so that the active constraint can add events to its own 
    /// list (events) without changing the one used during propagation (changes)
    IntStack changes; // this is the list that is accessible from a propagator
    IntStack events; // this is the list that collects the events
    /// The type of event for each modified variable
    Event *event_type;

    /// For each variable in the scope, the constraint list containing self
    //ConstraintList** self_list;
    /// The element of the list corresponding to this constraint
    int* self;
    /// The type of trigger (domain/range/value)
    int* trigger;

    /// residual supports used for generic ac propagation
    int ***supports;
    int *solution;
    //@}


    /*!@name Constructors*/
    //@{
    /// The _scope is build by copying an array "scp" containing "l" variables
    void initialise(Vector< IntVar >& scp, const float w=1.0);
    void triggerOn(const int t, const int x);
    void notifyAssignment(const int var);
    void entail();
    void restore();
    Constraint();
    Constraint(Vector< IntVar >& scp, const float w=1.0);
    virtual ~Constraint();

    /// An idempotent constraint should not be called on events triggered by itself.
    /// To forbid that, the lists 'events' and 'changes' are merged
    inline void set_idempotent(const bool idp) {
      if(idp) {
	events.size = 0;
	events.capacity = arity;
	events.list_ = changes.list_;
	events.index_ = changes.index_;
	events.start_ = NULL;
      } else {
	events.initialise(0, arity-1, false);
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
      if(changes.list_ == events.list_) 
	// if idempotent, clear the events list
	events.size = 0;      
    }


    void post();
    void relax();
    //@}

    void assign(const int var);
    void unassign(const int var);

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
    virtual IntVar propagate() { return NULL; }
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
    //virtual std::string getString() const ;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "c"; }
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
    ConstraintNotEqual(Vector< IntVar >& scp) : Constraint(scp) {
      triggerOn(_value_, 0);
      triggerOn(_value_, 1);
      set_idempotent(true);
    }
    virtual ~ConstraintNotEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check(const int* sol) const { return (sol[0] == sol[1]); }
    virtual IntVar propagate();
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
    ConstraintLess(Vector< IntVar >& scp, const int ofs=0) : Constraint(scp) {
      offset = ofs;
      triggerOn(_range_, 0);
      triggerOn(_range_, 1);
      set_idempotent(true);
    }
    virtual ~ConstraintLess() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return (sol[0]+offset > sol[1]); }
    virtual IntVar propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    //virtual std::string getString() const ;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "less_than"; }
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
    PredicateEqual(Vector< IntVar >& scp, const int sp=1) : Constraint(scp) {
      spin = sp;
      triggerOn(_domain_, 0);
      triggerOn(_domain_, 1);
      triggerOn(_value_, 2);
      set_idempotent(false);
    }
    virtual ~PredicateEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((sol[0] == sol[1]) == (sol[2] ^ spin)); }
    virtual IntVar propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    //virtual std::string getString() const ;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "are_equal"; }
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
    PredicateAdd(Vector< IntVar >& scp) : Constraint(scp) {
      triggerOn(_range_, 0);
      triggerOn(_range_, 1);
      triggerOn(_range_, 2);
      set_idempotent(true);
    }
    virtual ~PredicateAdd() {}
    //@}
    
    /**@name Solving*/
    //@{
    inline int check( const int* sol ) const { return (sol[2] != (sol[0]+sol[1])); }
    inline IntVar propagate();
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    //virtual std::string getString() const ;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "plus"; }
    //@}
  };


// /**********************************************
//  * Substraction Predicate
//  **********************************************/
// /*! \class PredicateSub
//  \brief  Binary Substraction predicate (x0 - x1 = y)
// */
// class PredicateSub : public Constraint
//   {
    
//   public:
//     /**@name Constructors*/
//     //@{
//     PredicateSub(Vector< IntVar >& scp) : Constraint(scp) {
//       triggerOn(_range_, 0);
//       triggerOn(_range_, 1);
//       triggerOn(_range_, 2);
//       set_idempotent(true);
//     }
//     virtual ~PredicateSub() {}
//     //@}
    
//     /**@name Solving*/
//     //@{
//     inline int check( const int* sol ) const { return (sol[2] != (sol[0]-sol[1])); }
//     inline IntVar propagate();
//     //@}
    
//     /**@name Miscellaneous*/
//     //@{  
//     virtual std::string getString() const ;
//     virtual std::string name() const { return "minus"; }
//     //@}
//   };


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
    ConstraintAllDiff(Vector< IntVar >& scp);
    virtual ~ConstraintAllDiff();
    void init();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual IntVar propagate();
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

    Queue triggers[NUM_PRIORITY];
    //IntStack active;
    //bool active[NUM_PRIORITY+1];
    int num_actives;
    BitSet _set_;
    
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
    
    inline bool empty() { return !num_actives; }
      //return !active.size; }//;
    inline void trigger(Constraint* cons, const int var, const Event evt)//;
    {
      int priority = cons->priority, cons_id = cons->id;
      if(_set_.fastMember(cons_id)) {
	if(cons->events.member(var)) {
	  cons->event_type[var] |= evt;
	} else {
	  cons->events.insert(var);
	  cons->event_type[var] = evt;
	}
      } else {
	_set_.fastInsert(cons_id);
	//if(triggers[priority].empty()) active.ordered_insert(priority);
	//if(triggers[priority].empty()) active.insert(priority);
	triggers[priority].add(cons_id);
	++num_actives;
	cons->events.clear();
	cons->events.insert(var);
	cons->event_type[var] = evt;
      }
    }
    inline int pop()//;
    {
      //int priority = active.back();
      //Constraint *cons = Solver::constraints[triggers[priority].pop()];
      int priority = NUM_PRIORITY;
      while(--priority) if(!triggers[priority].empty()) break;
      int cons = triggers[priority].pop();
      _set_.fastErase(cons);
      --num_actives;
      //if(triggers[priority].empty()) active.erase(priority);
      return cons;
    }
    inline void clear() {
      if(!empty()) {
	for(int i=0; i<NUM_PRIORITY; ++i) triggers[i].clear();
	//active.clear();
	num_actives = 0;
	_set_.clear();
      }
    }
  };



}

#endif //__CONSTRAINT_HPP
