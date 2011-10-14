
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
#include <mistral_variable.hpp>
#include <mistral_solver.hpp>
//#include <mistral_sat.hpp>



namespace Mistral {

  
  /**********************************************
   * Listener
   **********************************************/

  /*! \class RestartListener
    \brief RestartListener Class

    * Called whenever the solver restarts *
    
    This is used to implement procedures triggered by restarts
  */
  class RestartListener {
  public:
    int rid;
    virtual void notify_restart() = 0;
    
    std::ostream& display(std::ostream& os) { os << "restart-L"; return os; }    
  };

  /*! \class DecisionListener
    \brief DecisionListener Class

    * Called whenever the solver restarts *
    
    This is used to implement procedures triggered by decisions
  */
  class DecisionListener {
  public:
    int did;
    virtual void notify_decision() = 0;

    std::ostream& display(std::ostream& os) { os << "decision-L"; return os; }    
  };

  /*! \class SuccessListener
    \brief SuccessListener Class

    * Called whenever the solver restarts *
    
    This is used to implement procedures triggered on successful propagation steps
  */
  class SuccessListener {
  public:
    int sid;
    virtual void notify_success() = 0;

    std::ostream& display(std::ostream& os) { os << "success-L"; return os; }    
  };

  /*! \class FailureListener
    \brief FailureListener Class

    * Called whenever the solver restarts *
    
    This is used to implement procedures triggered on failed propagation steps
  */
  class FailureListener {
  public:
    int fid;
    virtual void notify_failure() = 0;

    std::ostream& display(std::ostream& os) { os << "failure-L"; return os; }    
  };

  /*! \class ConstraintListener
    \brief ConstraintListener Class

    * Called whenever the solver restarts *
    
    This is used to implement procedures triggered wheneever a new
    constraint is added to the solver, or an old constraint is posted/relaxed
  */
  class ConstraintListener {
  public:
    int cid;
    virtual void notify_post(Constraint *c) = 0;
    virtual void notify_relax(Constraint *c) = 0;
    virtual void notify_add(Constraint *c) = 0;

    std::ostream& display(std::ostream& os) { os << "constraint-L"; return os; }    
  };

  /*! \class VariableListener
    \brief VariableListener Class

    * Called whenever the solver restarts *
    
    This is used to implement procedures triggered wheneever a new
    variable is added to the solver, or when its domain representation changes
  */
  class VariableListener {
  public:
    int vid;
    virtual void notify_add() = 0;
    virtual void notify_change(const int idx) = 0;

    std::ostream& display(std::ostream& os) { os << "variable-L"; return os; }    
  };

  std::ostream& operator<<(std::ostream& os, DecisionListener& x);
  std::ostream& operator<<(std::ostream& os, RestartListener& x);
  std::ostream& operator<<(std::ostream& os, SuccessListener& x);
  std::ostream& operator<<(std::ostream& os, FailureListener& x);
  std::ostream& operator<<(std::ostream& os, ConstraintListener& x);
  std::ostream& operator<<(std::ostream& os, VariableListener& x);

  // std::ostream& operator<<(std::ostream& os, DecisionListener* x);
  // std::ostream& operator<<(std::ostream& os, RestartListener* x);
  // std::ostream& operator<<(std::ostream& os, SuccessListener* x);
  // std::ostream& operator<<(std::ostream& os, FailureListener* x);
  // std::ostream& operator<<(std::ostream& os, ConstraintListener* x);
  // std::ostream& operator<<(std::ostream& os, VariableListener* x);


  /*! \class LiteralActivityManager
    \brief LiteralActivityManager Class

    * Listener interface for literal activity *

    This structure is shared with the clause base constraint.
    The constraint updates the activity when doing conflict analysis
    while this listener implements the decay.
  */
  class LiteralActivityManager : public DecisionListener {

  public:

    /*\ TODO: change double* to vectors* and make it a variable listener \*/
    Solver *solver;
    double *lit_activity; 
    double *var_activity;
    int n_vars;
 
    double decay;

    //LiteralActivityManager(Solver *s, void *a=NULL) ;
    LiteralActivityManager(Solver *s) ;
    virtual ~LiteralActivityManager() ;

    virtual void notify_decision() ;    
    double *get_weight() ;
    
  };


  /*! \class FailureCountManager
    \brief FailureCountManager Class

    * Listener interface for weighted degree *
  */
  class FailureCountManager : public FailureListener, public ConstraintListener {

  public:

    Solver *solver;

    /*\ TODO: make it a variable listener \*/
    Vector<double> constraint_weight;
    Vector<double> variable_weight;

    //FailureCountManager(Solver *s, void *a=NULL) : solver(s) {// }
    FailureCountManager(Solver *s) : solver(s) {// }

      for(unsigned int i=0; i<solver->variables.size; ++i) {
	variable_weight.add(solver->variables[i].get_degree());
      }
      for(unsigned int i=0; i<solver->constraints.size; ++i) {
	constraint_weight.add(1);
      }

      solver->add((FailureListener*)this);
      solver->add((ConstraintListener*)this);
    }

    virtual ~FailureCountManager() {// }

      solver->remove((ConstraintListener*)this);
      solver->remove((FailureListener*)this);
    }

    double *get_weight() { return variable_weight.stack_; }     

    virtual void notify_failure() {
      Constraint *con = solver->culprit;
      if(con) {
	int i = con->scope.size, idx;
	++constraint_weight[con->id];
	while(i--) {
	  idx = con->scope[i].id();
	  if(idx>=0) ++variable_weight[idx];
	}
      }
    }

    virtual void notify_post(Constraint *con) {
      int i = con->stress, idx;
      while(i--) {
	idx = con->scope[con->active[i]].id();
	if(idx>=0) variable_weight[idx] += constraint_weight[con->id];
      }
    }

    virtual void notify_relax(Constraint *con) {
      int i = con->active.size, idx;
      while(i--) {
	idx = con->scope[con->active[i]].id();
	if(idx>=0) variable_weight[idx] -= constraint_weight[con->id];
      }
    }

    virtual void notify_add(Constraint *con) {
      unsigned int idx = con->id;
      while(constraint_weight.size <= idx) {
	constraint_weight.add(1.0);
      }
    }
  };


  /*! \class PruningCountManager
    \brief PruningCountManager Class

    * Listener interface for ABS *
  */
  class PruningCountManager : public SuccessListener {

  public:

    Solver *solver;

    Vector<double> variable_weight;

    //PruningCountManager(Solver *s, void *a=NULL) : solver(s) {
    PruningCountManager(Solver *s) : solver(s) {
      for(unsigned int i=0; i<solver->variables.size; ++i) {
	variable_weight.add(solver->variables[i].get_degree());
      }
      solver->add((SuccessListener*)this);
    }

    virtual ~PruningCountManager() {
      solver->remove((SuccessListener*)this);
    }

    double *get_weight() { return variable_weight.stack_; }     

    virtual void notify_success() {
      int // last_decision = solver->decisions.back().var.id(),
	id;
      unsigned int i = solver->trail_.back(), n=solver->saved_vars.size;
      while(++i<n) {
	id = solver->saved_vars[i]; //.id();
	//if(id != last_decision) 
	++variable_weight[id];
      }
    }
  };


  /*! \class RestartPolicy
    \brief  Interface RestartPolicy

    super class for restart-cutoff sequence generators
  */
  class RestartPolicy {
    
  public:

    unsigned int base;
    
    RestartPolicy(const unsigned int b=256);
    virtual void reset(unsigned int& limit) = 0;
    virtual void initialise(unsigned int& limit) = 0;
    
  };


  class NoRestart : public RestartPolicy {
    
  public:
    
    NoRestart();
    virtual ~NoRestart();
    
    void reset(unsigned int& limit) {
      limit = base;
    }

    void initialise(unsigned int& limit) {
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

    void initialise(unsigned int& limit) {
      limit = 0;
      increment = base;
      reset(limit);
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

    void initialise(unsigned int& limit) {
      iteration = 0;
      reset(limit);
    }
    
    
  };


  /**********************************************
   * Search Strategies
   **********************************************/

  /*! \class BranchingHeuristic
    \brief  Interface BranchingHeuristic

    The branching heuristic is queried by the
    solver on each new node in order to take
    a branching decision.

    It implements one method: branch(), which return a 
    Decision object representing the branching decision
    to make.
  */
  class BranchingHeuristic {

  public:
    
    Solver *solver;

    BranchingHeuristic() {}
    BranchingHeuristic(Solver *s) {solver = s;}
    virtual ~BranchingHeuristic() {}

    //virtual void initialise() {}
    virtual void initialise(VarStack< Variable >& seq) {}
    virtual void close() {}

    virtual Decision branch() = 0;

  };


  /**********************************************
   * Generic Variable/Value Ordering heuristics
   **********************************************/
  /*! \class GenericHeuristic
    \brief  Class GenericHeuristic

    Generic branching heuristic: 
    1/ a method to select the next variable to branch on (var)
    2/ a method to reduce the domain of this variable (choice)
  */
  template < class VarSelector, class ValSelector >
  class GenericHeuristic : public BranchingHeuristic {
  public:

    VarSelector var;
    ValSelector choice;

    GenericHeuristic(Solver *s) 
      : BranchingHeuristic(s) {
      var.initialise(s);
      choice = ValSelector(s);
    }


    // GenericHeuristic(Solver *s, void *a) 
    //   : BranchingHeuristic(s) {
    //   var.initialise(s,a);
    //   choice = ValSelector(s,a);
    // }

    virtual void initialise(VarStack< Variable >& seq) {var.initialise(seq);}

    virtual Decision branch() {
      return choice.make(var.select());
    }

  };


  /**********************************************
   * Generic variable heuristic
   **********************************************/
  
  /*! \class GenericDVO
    \brief  Class GenericDVO

    Generic (dynamic) variable ordering heuristic. 
    - A parameterized comparison method is used to select the 'best' variable
  */
  template < class VarComparator >
  class GenericDVO 
  {
  public: 
    
    /**@name Parameters*/
    //@{ 
    Solver   *solver;
    VarComparator best;
    VarComparator current;
    //@}

    /**@name Constructors*/
    //@{
    GenericDVO() { solver = NULL; }
    GenericDVO(Solver* s) : solver(s) {}
    //virtual void initialise(Solver *s, void *a=NULL) { solver = s; }
    virtual void initialise(Solver *s) { solver = s; }
    virtual void initialise(VarStack< Variable >& seq) {}
    //@}
    
    /**@name Utils*/
    //@{ 
    Variable select()
    {    
      Variable *variables = solver->sequence.list_;
      unsigned int length = solver->sequence.size-1;
      Variable var = variables[length];
      best = var;
      for(unsigned int i=length; i--;) {
	current = variables[i];
	if( current < best ) {
	  best = current;
	  var = variables[i];
	  
	  //std::cout << "*" << std::endl;
	}
      }
      return var;
    }
    //@}
  };


  /*! \class GenericRandomDVO
    \brief  Class GenericRandomDVO

    Randomized Generic (dynamic) variable ordering heuristic 
    - A parameterized comparison method is used to select 
      the k 'best' variables, then one of them is randomly selected
  */
  template < class VarComparator, int RAND >
  class GenericRandomDVO 
  {
  public: 

    /**@name Parameters*/
    //@{ 
    Solver      *solver;
    VarComparator   bests[RAND+1];
    VarComparator current;
    Variable   bestvars[RAND+1];
    //@}

    /**@name Constructors*/
    //@{
    GenericRandomDVO()
    {
      solver = NULL;
    }
    GenericRandomDVO(Solver* s)  : solver(s) 
    {
    }
    //virtual void initialise(Solver *s, void *a=NULL) { solver = s; }
    virtual void initialise(Solver *s) { solver = s; }
    virtual void initialise(VarStack< Variable >& seq) {}

    virtual ~GenericRandomDVO() 
    {
    }
    //@}

    /**@name Utils*/
    //@{ 
    Variable select()
    {
      Variable *variables = solver->sequence.list_;
      unsigned int length = solver->sequence.size-1;
      unsigned int realsize=1, i, j;
      bests[0] = bestvars[0] = variables[length];
      for(j=length; j--;)
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
	  
	  if(realsize<RAND) ++realsize;
	}
      return bestvars[randint(realsize)];
    }
    //@}
  };


  /*! \class GenericWeightedDVO
    \brief  Class GenericWeightedDVO

    Generic (dynamic) variable ordering heuristic with a weight function. 
    A 'WeightManager' is used to compute each variable's weight
  */
  template < class WeightManager, class VarComparator >
  class GenericWeightedDVO : public GenericDVO< VarComparator >
  {
  public: 
    
    /**@name Parameters*/
    //@{ 
    WeightManager *manager;
    //@}

    /**@name Constructors*/
    //@{
    GenericWeightedDVO() : GenericDVO< VarComparator >() { manager=NULL; }
    GenericWeightedDVO(Solver* s) : GenericDVO< VarComparator >(s) {
      initialise_manager();
    }
    virtual ~GenericWeightedDVO() { // manager->close();
      delete manager; }
    //virtual void initialise(Solver *s, void *a=NULL) { 
    virtual void initialise(Solver *s) { 
      GenericDVO< VarComparator >::initialise(s);
      initialise_manager();
    }
    //virtual void close() { manager->close(); }
    //virtual void initialise() { manager->initialise(); }
    virtual void initialise(VarStack< Variable >& seq) {}
    virtual void initialise_manager() {
      manager = new WeightManager(GenericDVO< VarComparator >::solver);
      //GenericDVO< VarComparator >::solver->add(manager);
      GenericDVO< VarComparator >::best.weight = manager->get_weight();
      GenericDVO< VarComparator >::current.weight = manager->get_weight();
    }
    //@}
    
  };

  
  /*! \class GenericWeightedRandomDVO
    \brief  Class GenericWeightedRandomDVO

    Generic (dynamic) variable ordering heuristic with a weight function. 
    A 'WeightManager' is used to compute each variable's weight.
    This is the randomized version. i.e., the k best variables are computed
    then one of them is randomly selected
  */
  template < class WeightManager, class VarComparator, int RAND >
  class GenericWeightedRandomDVO : public GenericRandomDVO< VarComparator, RAND >
  {
  public: 
    
    /**@name Parameters*/
    //@{ 
    WeightManager *manager;
    //@}

    /**@name Constructors*/
    //@{
    GenericWeightedRandomDVO() : GenericRandomDVO< VarComparator, RAND >() { manager=NULL; }
    GenericWeightedRandomDVO(Solver* s) : GenericRandomDVO< VarComparator, RAND >(s) {
      initialise_manager();
    }
    virtual ~GenericWeightedRandomDVO() { manager->close(); delete manager; }
    //virtual void initialise(Solver *s, void *a=NULL) { 
    virtual void initialise(Solver *s) { 
      GenericRandomDVO< VarComparator, RAND >::initialise(s);
      initialise_manager();
    }
    //virtual void initialise(Solver *s) { manager->initialise(s); }
    virtual void initialise(VarStack< Variable >& seq) {}
    virtual void initialise_manager() {
      manager = new WeightManager(GenericRandomDVO< VarComparator, RAND >::solver);
      //GenericRandomDVO< VarComparator >::solver->add(manager);
      GenericRandomDVO< VarComparator, RAND >::current.weight = manager->get_weight();
      for(int i=0; i<=RAND; ++i)
	GenericRandomDVO< VarComparator, RAND >::bests[i].weight = manager->get_weight();
    }
    //@}
    
  };


  /*! \class NoOrder
    \brief  Class NoOrder

    This heuristic selects the first variable in the sequence, 
    there is no way of knowing which it will be (in particular
    this is NOT a lexicographic heuristic, and this is not
    even a static heuristic)
  */
  class NoOrder {

  public: 
    Solver *solver;

    NoOrder() { solver = NULL; }
    NoOrder(Solver *s);
    virtual ~NoOrder();
    void initialise(Solver *s) { solver = s; }
    //void initialise(Solver *s, void *a) { solver = s; }
    void initialise(VarStack< Variable >& seq) {}
    
    Variable select();

  };


  /*! \class Lexicographic
    \brief  Class Lexicographic

    This heuristic selects the variable with lowest rank
    in the initial sequence of search variables.
  */
  class Lexicographic : public VariableListener {

  public: 
    
    Solver *solver;
    Vector< Variable >           order;
    Vector< int >                index;
    ReversibleNum< unsigned int > last;

    Lexicographic() { solver = NULL; }
    Lexicographic(Solver *s);
    void initialise(Solver *s);
    //void initialise(Solver *s, void *a=NULL);
    void initialise(VarStack< Variable >& seq);
    virtual ~Lexicographic();

    virtual void notify_add() {};
    virtual void notify_change(const int idx);
    
    Variable select();

  };


  /**********************************************
   * Variable Comparators
   **********************************************/

  /*! \class MinDomain
    \brief  Class MinDomain

    Order two variables by their domain sizes
  */
  class MinDomain 
  {
  public: 

    /**@name Constructors*/
    //@{
    MinDomain() {dom_ = LARGE_VALUE;}
    //@}

    /**@name Parameters*/
    //@{ 
    int dom_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (double)dom_; } 
    inline bool operator<( MinDomain& x ) const { return dom_ < x.dom_; }
    inline void operator=( MinDomain& x ) { dom_ = x.dom_; }
    inline void operator=( Variable x ) { dom_ = x.get_size(); }
    //@}  

    std::ostream& display(std::ostream& os) {
      os << dom_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MinDomain& x);


  /*! \class MinDomainOverDegree
    \brief  Class MinDomainOverDegree

    Order two variables by the ratio of their domain sizes and degree
  */
  class MinDomainOverDegree 
  {
  public: 

    /**@name Constructors*/
    //@{
    MinDomainOverDegree() {dom_ = LARGE_VALUE; deg_ = 0;}
    //@}

    /**@name Parameters*/
    //@{ 
    int dom_;
    int deg_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (deg_ ? (double)dom_/(double)deg_ : (double)INFTY); } 
    inline bool operator<( MinDomainOverDegree& x ) const { return dom_*x.deg_ < x.dom_*deg_; }
    inline void operator=( MinDomainOverDegree& x ) { dom_ = x.dom_; deg_ = x.deg_; }
    inline void operator=( Variable x ) { dom_ = x.get_size(); deg_ = x.get_degree(); }
    //@}  

    std::ostream& display(std::ostream& os) {
      os << dom_ << "/" << deg_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MinDomainOverDegree& x);


  /*! \class MinDomainOverWeight
    \brief  Class MinDomainOverWeight

    Order two variables by the ratio of their domain sizes and weight
  */
  class MinDomainOverWeight
  {
  public: 

    /**@name Constructors*/
    //@{
    MinDomainOverWeight() {dom_ = LARGE_VALUE; wei_ = 0;}
    //MinDomainOverWeight(void *w) : weight((double*)w) {dom_ = LARGE_VALUE; wei_ = 0;}
    //@}

    /**@name Parameters*/
    //@{ 
    double *weight;
    int dom_;
    double wei_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (wei_ ? (double)dom_/wei_ : (double)INFTY); } 
    inline bool operator<( MinDomainOverWeight& x ) const { return dom_*x.wei_ < x.dom_*wei_; }
    inline void operator=( MinDomainOverWeight& x ) { dom_ = x.dom_; wei_ = x.wei_; }
    inline void operator=( Variable x ) { dom_ = x.get_size(); wei_ = weight[x.id()]; }
    //@}  

    std::ostream& display(std::ostream& os) {
      os << dom_ << "/" << wei_;
      return os;
    }
  };


  /*! \class MinNeighborDomainOverNeighborWeight
    \brief  Class MinNeighborDomainOverNeighborWeight

    Order two variables by the ratio of the domain sizes and weight
    of a vector of other variables
  */
  class MinNeighborDomainOverNeighborWeight
  {
  public: 

    /**@name Constructors*/
    //@{
    MinNeighborDomainOverNeighborWeight() {dom_ = LARGE_VALUE; wei_ = 0;}
    //@}

    /**@name Parameters*/
    //@{ 
    double *weight;
    Vector< Variable > *map;
    int dom_;
    double wei_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (wei_ ? (double)dom_/wei_ : (double)INFTY); } 
    inline bool operator<( MinNeighborDomainOverNeighborWeight& x ) const { return dom_*x.wei_ < x.dom_*wei_; }
    inline void operator=( MinNeighborDomainOverNeighborWeight& x ) { dom_ = x.dom_; wei_ = x.wei_; }
    inline void operator=( Variable x ) { 
      int idx = x.id();
      int i = map[idx].size;
      Variable y;
      dom_ = 0;
      wei_ = 0;
      while(i--) {
  	y = map[idx][i];
  	dom_ += y.get_size(); 
  	wei_ += weight[y.id()];
      } 

      //std::cout << x << ": " << dom_ << "/" << wei_ << std::endl;

    }
    //@}  

    std::ostream& display(std::ostream& os) {
      os << dom_ << "/" << wei_;
      return os;
    }
  };


  /*! \class MinNeighborDomainOverWeight
    \brief  Class MinNeighborDomainOverWeight

    Order two variables by the ratio of the domain sizes and weight
    of a vector of other variables
  */
  class MinNeighborDomainOverWeight
  {
  public: 

    /**@name Constructors*/
    //@{
    MinNeighborDomainOverWeight() {dom_ = LARGE_VALUE; wei_ = 0;}
    //@}

    /**@name Parameters*/
    //@{ 
    double *weight;
    Vector< Variable > *map;
    int dom_;
    double wei_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (wei_ ? (double)dom_/wei_ : (double)INFTY); } 
    inline bool operator<( MinNeighborDomainOverWeight& x ) const { return dom_*x.wei_ < x.dom_*wei_; }
    inline void operator=( MinNeighborDomainOverWeight& x ) { dom_ = x.dom_; wei_ = x.wei_; }
    inline void operator=( Variable x ) { 
      int idx = x.id();
      int i = map[idx].size;
      Variable y;
      wei_ = weight[idx];
      dom_ = 0;
      while(i--) {
  	y = map[idx][i];
  	dom_ += y.get_size(); 
      } 
    }
    //@}  

    std::ostream& display(std::ostream& os) {
      os << dom_ << "/" << wei_;
      return os;
    }
  };


  /*! \class MaxWeight
    \brief  Class MaxWeight

    Order two variables by their weights
  */
  class MaxWeight
  {
  public: 

    /**@name Constructors*/
    //@{
    MaxWeight() {wei_ = 0;}
    //MaxWeight(int *w) : weight(w) {wei_ = 0;}
    //@}

    /**@name Parameters*/
    //@{ 
    double *weight;
    double wei_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return wei_ ; } 
    inline bool operator<( MaxWeight& x ) const { return x.wei_ < wei_; }
    inline void operator=( MaxWeight& x ) { wei_ = x.wei_; }
    inline void operator=( Variable x ) { 
      //std::cout << "activity of " << x << " = " << weight[x.id()] << std::endl;
      wei_ = weight[x.id()]; 
    }
    //@}  

    std::ostream& display(std::ostream& os) {
      os << wei_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MaxWeight& x);


  /**********************************************
   * Branching Decisions
   **********************************************/

  /*! \class MinValue
    \brief  Class MinValue

    Assigns the variable to its minimum value.
  */
   class MinValue {

  public: 
    
    MinValue() {}
    MinValue(Solver *s) {}
     //MinValue(Solver *s, void *a) {}
    virtual ~MinValue() {};
    
    inline Decision make(Variable x) {
      Decision d(x, Decision::ASSIGNMENT, x.get_min());
      return d;
    }

  };

  /*! \class MaxValue
    \brief  Class MaxValue

    Assigns the variable to its maximum value.
  */
  class MaxValue {

  public: 
    
    MaxValue() {}
    MaxValue(Solver *s) {}
    //MaxValue(Solver *s, void *a) {}
    virtual ~MaxValue() {};
    
    inline Decision make(Variable x) {
      Decision d(x, Decision::ASSIGNMENT, x.get_max());
      return d;
    }

  };

  /*! \class HalfSplit
    \brief  Class HalfSplit

    Set the upper bound of the variable to (ub+lb)/2.
  */
  class HalfSplit {

  public: 
    
    HalfSplit() {}
    HalfSplit(Solver *s) {}
    //HalfSplit(Solver *s, void *a) {}
    virtual ~HalfSplit() {};
    
    inline Decision make(Variable x) {
      Decision d(x, Decision::UPPERBOUND, (x.get_min()+x.get_max())/2);
      return d;
    }

  };

  /*! \class RandomMinMax
    \brief  Class RandomMinMax

    Assigns the variable randomly its minimum or maximum value.
  */
  class RandomMinMax {

  public: 
    
    RandomMinMax() {}
    RandomMinMax(Solver *s) {}
    //RandomMinMax(Solver *s, void *a) {}
    virtual ~RandomMinMax() {};
    
    inline Decision make(Variable x) {
      Decision d(x, Decision::ASSIGNMENT, 
		 (randint(2) ? x.get_min() : x.get_max()));
      return d;
    }

  };


  /*! \class MinWeightValue
    \brief  Class MinWeightValue

    Assigns the variable to its value with minimum weight for some weight matrix.
  */
  class MinWeightValue {

  public: 

    double **weight;
    
    MinWeightValue() {}
    MinWeightValue(Solver *s) {}
    //MinWeightValue(Solver *s, void *a) { weight = (double**)a; }
    virtual ~MinWeightValue() {};
    
    inline Decision make(Variable x) {
      int // id_x = x.id(),
	best_val = x.get_min();
      double *wgt = weight[x.id()];
      //double min_weight = weight[id_x][best_val]// [best_val][id_x]
      double min_weight = wgt[best_val]// [best_val][id_x]
	, aux_weight;
      int vali, vnxt=x.next(best_val);
      do {
	vali = vnxt;
	vnxt = x.next(vali);
	aux_weight = wgt[vali]; //weight[id_x][vali]; //weight[vali][id_x];
	if(aux_weight < min_weight) {
	  min_weight = aux_weight;
	  best_val = vali;
	}
      } while(vali<vnxt);
      
      Decision d(x, Decision::ASSIGNMENT, best_val);
      return d;
    }

  };


#define NEG(a) ((2*a))
#define POS(a) ((2*a+1))
  /*! \class BoolMinWeightValue
    \brief  Class BoolMinWeightValue

    Assigns the variable to its value with minimum weight for some weight matrix.
    Assumes that the variable is Boolean
  */
  class BoolMinWeightValue {

  public: 

    //double **weight;
    double *weight;
    
    BoolMinWeightValue() {}
    BoolMinWeightValue(Solver *s) {}
    // BoolMinWeightValue(Solver *s, void *a) {
    //   weight = (double*)a;
    // }
    virtual ~BoolMinWeightValue() {};
    
    inline Decision make(Variable x) {

      // this thing is tricky: the value 'lit_activity' of a literal gets incremented when we post a new clause 
      // involving this literal. Now, this means that the opposite literal is found to be more constrained, 
      // hence when weighting clause one should use the opposite literal, and when branching one should use the literal
      // with highest activity

      Atom a = x.id();
      Decision d(x, Decision::ASSIGNMENT, (weight[NEG(a)]<weight[POS(a)] ? 1 : 0));
      return d;
    }
    
  };
  
  
  
  class VSIDS : public GenericHeuristic< GenericWeightedDVO< LiteralActivityManager, MaxWeight >, 
					 BoolMinWeightValue > {
  public:
    
    VSIDS(Solver *s) : GenericHeuristic< GenericWeightedDVO< LiteralActivityManager, MaxWeight >, 
					 BoolMinWeightValue >(s) {
      choice.weight = s->get_literal_activity();
  }
};


// template < Branching > 
// class SchedulingWeightedDegree 
//   : public GenericHeuristic< GenericWeightedDVO< FailureCountManager,  VarComparator >,
// 			     Branching > {
// public:
  
//   SchedulingWeightedDegree(Solver *s, Vector< Vector< Variable > >& graph) 
//     : GenericHeuristic< GenericWeightedDVO< FailureCountManager, VarComparator >, 
// 			Branching >(s) {
//     var.best.map = graph.stack_;
//     var.current.map = graph.stack_;
//   }  
// };

}

#endif // __SEARCH_HPP
