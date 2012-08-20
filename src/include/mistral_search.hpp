
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

#include <algorithm>

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

  /*! \class SolutionListener
    \brief SolutionListener Class

    * Called whenever the solver solutions *
    
    This is used to implement procedures triggered by solutions
  */
  class SolutionListener {
  public:
    int mid;
    virtual void notify_solution() = 0;
    
    std::ostream& display(std::ostream& os) { os << "solution-L"; return os; }    
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
    virtual void notify_post   (Constraint c) = 0;
    virtual void notify_relax  (Constraint c) = 0;
    virtual void notify_add_con(Constraint c) = 0;

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
    virtual void notify_add_var() = 0;
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

    virtual std::ostream& display(std::ostream& os, const bool all) const {return os;}
    
  };


  //double *weight_sorting_array;
  int decreasing_weight(const void *x, const void *y);


  /*! \class NoManager
    \brief NoManager Class
    
    * Manager doing nothing, used for genericity *
    */
  class NoManager : public FailureListener, public ConstraintListener {
    
  public:

    NoManager(Solver *s) {}
    virtual ~NoManager() {}

    double *get_weight() { return NULL; }  

    std::ostream& display(std::ostream& os) { return os; }    
    // void initialise(Solver* s);
    // void initialise_structures(VarComparator &v);
  };



  /*! \class FailureCountManager
    \brief FailureCountManager Class

    * Listener interface for weighted degree *
  */
  //template< float DECAY > 
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

    virtual void check_consistency() {
      
      double xweight;

      //std::cout << solver << std::endl;
      
      solver->display(std::cout, 1);


      for(unsigned int i=0; i<variable_weight.size; ++i) {
	
	if(!(solver->domain_types[i] & REMOVED_VAR) && solver->sequence.contain(i)) {

	  xweight = 0;
	  for(Event trig = 0; trig<3; ++trig) 
	    for(int cons = solver->constraint_graph[i].on[trig].size; --cons>=0;) {
	      xweight += constraint_weight[solver->constraint_graph[i].on[trig][cons].id()];
	    }

	  if(xweight != variable_weight[i]) {

	    std::cout << "WARNING! inconsistency: on " << solver->variables[i] << ": " 
		      << variable_weight[i] << " should be " << xweight << std::endl;

	  } else {
	    
	    std::cout << "OK!" << std::endl;

	  }
	}
      }

    }
  

    virtual void notify_failure() {
      int i;

      // if(DECAY >= 0) {
      // 	for(i=variable_weight.size; --i>=0;)
      // 	  variable_weight[i] *= DECAY;
      // 	for(i=constraint_weight.size; --i>=0;)
      // 	  constraint_weight[i] *= DECAY;
      // }


      //std::cout << "increment weight of ";

      Constraint con = solver->culprit;




      if(!con.empty()) {


	// std::cout << "--failure on [" << con.id() << "] " << con << "--\n";

	Variable *scope = con.get_scope();
	int idx;
	i = con.arity();
	++constraint_weight[con.id()];
	while(i--) {
	  idx = scope[i].id();
	  if(idx>=0) {

	    // std::cout << scope[i] << " (" << variable_weight[idx] << " => ";

	    ++variable_weight[idx];

	    // std::cout << variable_weight[idx] << ") ";
	    
	  }
	}
      } 

      // else {

      // 	std::cout << "--failure on [nill]?? " << solver->wiped_idx << "--\n";
	
      // }

      // std::cout << std::endl;


      // display(std::cout);

      // std::cout << std::endl;

      //check_consistency();

    }

    virtual void notify_post(Constraint con) {
      int i = con.num_active(), idx;
      Variable *scope = con.get_scope();
      while(i--) {
	idx = scope[con.get_active(i)].id();
	if(idx>=0) variable_weight[idx] += constraint_weight[con.id()];
      }
    }

    virtual void notify_relax(Constraint con) {
      int i = con.num_active(), idx;
      Variable *scope = con.get_scope();
      while(i--) {
	idx = scope[con.get_active(i)].id();
	if(idx>=0) variable_weight[idx] -= constraint_weight[con.id()];
      }
    }

    virtual void notify_add_con(Constraint con) {
      unsigned int idx = con.id();
      //Variable *scope = con.get_scope();
      while(constraint_weight.size <= idx) {
	constraint_weight.add(1.0);
      }
    }

    virtual std::ostream& display(std::ostream& os, const bool all) const ;
// {
      
//       int *all_variables = new int[variable_weight.size];
//       int *all_constraints = new int[constraint_weight.size];

//       for(unsigned int i=0; i<variable_weight.size; ++i) {
// 	all_variables[i] = i;
//       }

//       for(unsigned int i=0; i<constraint_weight.size; ++i) {
// 	all_constraints[i] = i;
//       }

//       weight_sorting_array = variable_weight.stack_;
//       qsort(all_variables, variable_weight.size, sizeof(int), decreasing_weight);

//       weight_sorting_array = constraint_weight.stack_;
//       qsort(all_constraints, constraint_weight.size, sizeof(int), decreasing_weight);

//       for(unsigned int i=0; i<variable_weight.size; ++i) {
// 	os << std::setw(5) << solver->variables[all_variables[i]] << " ";
//       }
//       os << std::endl;
//       for(unsigned int i=0; i<variable_weight.size; ++i) {
// 	os << std::setw(5) << variable_weight[i] << " ";
//       }
//       os << std::endl;

//      for(unsigned int i=0; i<constraint_weight.size; ++i) {
// 	os << std::setw(5) << solver->constraints[all_constraints[i]] << " ";
//       }
//       os << std::endl;
//       for(unsigned int i=0; i<constraint_weight.size; ++i) {
// 	os << std::setw(5) << constraint_weight[i] << " ";
//       }
//       os << std::endl;

//     }    

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
      int i = solver->trail_.back(), n=solver->saved_vars.size;


      // if(DECAY >= 0) {
      // 	for(i=variable_weight.size; --i>=0;)
      // 	  variable_weight[i] *= DECAY;
      // }


      //std::cout << "increment weight of ";
      while(++i<n) {
	
	id = solver->saved_vars[i]; //.id();

	//std::cout << solver->variables[solver->saved_vars[i]] << " ";
	
	//if(id != last_decision) 
	++variable_weight[id];
      }
      
      //std::cout << std::endl;
    }

    virtual std::ostream& display(std::ostream& os, const bool all) const ;
  };



  // /*! \class ProgressSavingManager
  //   \brief ProgressSavingManager Class

  //   * Listener interface for progress saving *
  // */
  // class ProgressSavingManager : public FailureListener {

  // public:

  //   Solver *solver;

  //   int max_solution_length;
  //   int best_objective;

  //   Vector<int> progress;
    

  //   //ProgressSavingManager(Solver *s, void *a=NULL) : solver(s) {
  //   ProgressSavingManager(Solver *s) : solver(s) {
  //     best_objective = objective->value();
  //     max_solution_length = 0;
  //     for(unsigned int i=0; i<solver->variables.size; ++i) {
  // 	progress.add(solver->variables[i].get_min());
  //     }
  //     solver->add((FailureListener*)this);
  //   }

  //   virtual ~ProgressSavingManager() {
  //     solver->remove((FailureListener*)this);
  //   }

  //    virtual void notify_failure() {

  //   }
  // };


  // /*! \class GuidedSearchManager
  //   \brief GuidedSearchManager Class

  //   * Listener interface for progress saving *
  // */
  // class GuidedSearchManager : public FailureListener {

  // public:

  //   Solver *solver;

  //   int best_objective;

  //   Vector<int> best_solution;
    

  //   //GuidedSearchManager(Solver *s, void *a=NULL) : solver(s) {
  //   GuidedSearchManager(Solver *s) : solver(s) {
  //     best_objective = objective->value();
  //     max_solution_length = 0;
  //     for(unsigned int i=0; i<solver->variables.size; ++i) {
  // 	progress.add(solver->variables[i].get_min());
  //     }
  //     solver->add((FailureListener*)this);
  //   }

  //   virtual ~GuidedSearchManager() {
  //     solver->remove((FailureListener*)this);
  //   }

  //    virtual void notify_failure() {

  //   }
  // };


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
      // unsigned int increment = (base * luby_seq(++iteration));
      // std::cout << "restart for " << increment << std::endl;
      // limit += increment;

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
    virtual void initialise(VarStack< Variable, ReversibleNum<int> >& seq) {}
    virtual void close() {}

    virtual Decision branch() = 0;

    virtual std::ostream& display(std::ostream& os) = 0;
  };


// std::ostream& operator<<(std::ostream& os, BranchingHeuristic& x) {
//   return x.display(os);
// }

// std::ostream& operator<<(std::ostream& os, BranchingHeuristic* x) {
//   return x->display(os);
// }


  // std::ostream& operator<<(std::ostream& os, BranchingHeuristic& x);
  // std::ostream& operator<<(std::ostream& os, BranchingHeuristic* x);


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

    virtual void initialise(VarStack< Variable, ReversibleNum<int> >& seq) {var.initialise(seq);}

    virtual Decision branch() {
      return choice.make(var.select());
    }

    virtual std::ostream& display(std::ostream& os) {
      //os << "Branch on the variable " << std::endl;
      var.display(os);
      os << " c ";
      choice.display(os);
      os << std::endl ;
      return os;
    }
  };





  template< class VarComparator >
  class Identifiable {
    
  public :
    
    VarComparator criterion;
    int id;
    
    /**@name Utils*/
    //@{
    inline double value() { return criterion.value(); } 
    inline bool operator<( const Identifiable<VarComparator>& x ) const { 
      if(criterion < x.criterion) return true; 
      else if(x.criterion < criterion) return false;
      return (id > x.id);
    }
    inline void operator=( const Identifiable<VarComparator>& x ) { criterion = x.criterion; id = x.id; }
    inline void operator=( const Variable x ) { criterion = x; id = x.id(); }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      return criterion.display_criterion(os);
    }

    std::ostream& display(std::ostream& os) const {
      return criterion.display(os);
    }

  };



  /**********************************************
   * Generic variable heuristic
   **********************************************/
  
  template < class VarComparator >
  class GenericDVOI 
  {
  public: 
    
    Solver   *solver;

    GenericDVOI() { solver = NULL; }
    GenericDVOI(Solver* s) : solver(s) {}
    virtual void initialise(Solver *s) { solver = s; }

    virtual ~GenericDVOI() {}
    
    virtual std::ostream& display(std::ostream& os,  const int n, double* weights) {

      os << " c Select the " ;
      if(n>1) os << n << " ";
      os << "variable" << (n > 1 ? "s " : " ") << "with minimal value of ";
      VarComparator v;
      v.display_criterion(os);
      if(n>1) os << " and pick one uniformly at random" ;
      os << std::endl;

      Variable *variables = solver->sequence.list_;
      unsigned int length = solver->sequence.size-1;
      Variable var = variables[length];


      os << "--> branch in [";

      std::vector< Identifiable< VarComparator > > all_vars;
      for(unsigned int i=0; i<=length; ++i) {
      	Identifiable<VarComparator> vc;
	//if(weights) vc.criterion.weight = weights;
      	vc = variables[i];
      	vc.id = i;
      	all_vars.push_back(vc);

	os << variables[i] << " in " << variables[i].get_domain() << " ";

      }
      
      os << std::endl;


      sort(all_vars.begin(), all_vars.end());

      os << " c [" << variables[all_vars[0].id].id() << ":";
      all_vars[0].display(os);

      for(int i=1; i<n; ++i) {
	os << " " << variables[all_vars[i].id].id() << ":";
      	all_vars[i].display(os);
      }

      os << "]";

      for(unsigned int i=n; i<all_vars.size(); ++i) {
	os << " " << variables[all_vars[i].id].id() << ":";
      	all_vars[i].display(os);
      }
      os << std::endl;
      
      return os;
    }

  };


  /*! \class GenericDVO
    \brief  Class GenericDVO

    Generic (dynamic) variable ordering heuristic. 
    - A parameterized comparison method is used to select the 'best' variable
  */
  template < class VarComparator >
  class GenericDVO : public GenericDVOI< VarComparator >
  {
  public: 
    
    /**@name Parameters*/
    //@{ 
    VarComparator best;
    VarComparator current;
    //@}

    /**@name Constructors*/
    //@{
    GenericDVO() : GenericDVOI< VarComparator >() {}
    GenericDVO(Solver* s) : GenericDVOI< VarComparator >(s) {}
    virtual void initialise(Solver *s) { GenericDVOI< VarComparator >::initialise(s); }
    //virtual void initialise(Solver *s, void *a=NULL) { solver = s; }
    virtual void initialise(VarStack< Variable, ReversibleNum<int> >& seq) {}
   
    virtual ~GenericDVO() {}
    //@}
    
    /**@name Utils*/
    //@{ 
    Variable select()
    {    
      Variable *variables = GenericDVOI< VarComparator >::solver->sequence.list_;
      unsigned int length = GenericDVOI< VarComparator >::solver->sequence.size-1;
      Variable var = variables[length];
      best = var;

      //std::cout << variables[length] << "*" << std::endl ;

      for(unsigned int i=length; i--;) {
	current = variables[i];
	//std::cout << variables[i] ;

	if( current < best ) {
	  best = current;
	  var = variables[i];
	  
	  //std::cout << "*" ;
	}
	//std::cout << std::endl;
      }
      return var;
    }
    //@}


    virtual std::ostream& display(std::ostream& os) {
      return GenericDVOI< VarComparator >::display(os, 1, NULL);
    }


  };


  /*! \class GenericRandomDVO
    \brief  Class GenericRandomDVO

    Randomized Generic (dynamic) variable ordering heuristic 
    - A parameterized comparison method is used to select 
      the k 'best' variables, then one of them is randomly selected
  */
  template < class VarComparator, int RAND >
  class GenericRandomDVO : public GenericDVOI< VarComparator >
  {
  public: 

    /**@name Parameters*/
    //@{ 
    //Solver      *solver;
    VarComparator   bests[RAND+1];
    VarComparator current;
    Variable   bestvars[RAND+1];
    //@}

    /**@name Constructors*/
    //@{
    GenericRandomDVO() : GenericDVOI< VarComparator >() {}
    GenericRandomDVO(Solver* s) : GenericDVOI< VarComparator >(s) {}
   
    //virtual void initialise(Solver *s, void *a=NULL) { solver = s; }
    //virtual void initialise(Solver *s) { solver = s; }
    virtual void initialise(VarStack< Variable, ReversibleNum<int> >& seq) {}

    virtual ~GenericRandomDVO() {}
    //@}

    /**@name Utils*/
    //@{ 
    Variable select()
    {
      Variable *variables = GenericDVOI< VarComparator >::solver->sequence.list_;
      unsigned int length = GenericDVOI< VarComparator >::solver->sequence.size-1;
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
      return bestvars[(realsize>1 ? randint(realsize) : 0)];
    }
    //@}


    virtual std::ostream& display(std::ostream& os) {
      return GenericDVOI< VarComparator >::display(os, RAND, NULL);
    }

    // virtual std::ostream& display(std::ostream& os) const {
    //   bests[0].display_criterion(os);
    //   os << std::endl;
    //   Variable *variables = solver->sequence.list_;
    //   unsigned int length = solver->sequence.size-1;
    //   unsigned int realsize=1, i, j;

    //   // VarComparator   _bests[RAND+1];
    //   // VarComparator _current;
    //   // Variable   _bestvars[RAND+1];

    //   current.initialise();
    //   for(i=0; i<RAND; ++i) bests[i].initialise();

    //   bests[0] = bestvars[0] = variables[length];
    //   for(j=length; j--;)
    // 	{  
    // 	  current = variables[j];
    // 	  i = realsize;
    // 	  while( i && current < bests[i-1] ) {
    // 	    bests[i] = bests[i-1];
    // 	    bestvars[i] = bestvars[i-1];
    // 	    --i;
    // 	  }
    // 	  bests[i] = current;
    // 	  bestvars[i] = variables[j];
	  
    // 	  if(realsize<RAND) ++realsize;
    // 	}
      
    //   for(i=length; i--;) {
    // 	current = variables[i];
    // 	os << variables[i] << ": " ;
    // 	current.display(os) ;

    // 	for(j=0; j<realsize; ++j) {
    // 	  if(bestvars[j].id() == variables[i].id()) os << " <<==[CHOICE " << (j+1) << "]" ;
    // 	}

    // 	os << std::endl;
    //   }

    //   return os;
    // }

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
      GenericDVOI< VarComparator >::initialise(s);
      initialise_manager();
    }
    //virtual void close() { manager->close(); }
    //virtual void initialise() { manager->initialise(); }
    virtual void initialise(VarStack< Variable, ReversibleNum<int> >& seq) {}
    virtual void initialise_manager() {
      manager = new WeightManager(GenericDVO< VarComparator >::solver);
      //GenericDVO< VarComparator >::solver->add(manager);
      GenericDVO< VarComparator >::best.weight = manager->get_weight();
      GenericDVO< VarComparator >::current.weight = manager->get_weight();
    }
    //@}

    virtual std::ostream& display(std::ostream& os) {
      manager->display(os, false);
      return GenericDVOI< VarComparator >::display(os, 1, manager->get_weight());
    }
    
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
    virtual ~GenericWeightedRandomDVO() { //manager->close(); 
      delete manager; }
    //virtual void initialise(Solver *s, void *a=NULL) { 
    virtual void initialise(Solver *s) { 
      GenericDVOI< VarComparator >::initialise(s);
      initialise_manager();
    }
    //virtual void initialise(Solver *s) { manager->initialise(s); }
    virtual void initialise(VarStack< Variable, ReversibleNum<int> >& seq) {}
    virtual void initialise_manager() {
      manager = new WeightManager(GenericRandomDVO< VarComparator, RAND >::solver);
      //GenericRandomDVO< VarComparator >::solver->add(manager);
      if(manager->get_weight()) {
	GenericRandomDVO< VarComparator, RAND >::current.weight = manager->get_weight();
	for(int i=0; i<=RAND; ++i)
	  GenericRandomDVO< VarComparator, RAND >::bests[i].weight = manager->get_weight();
      }
    }

    virtual std::ostream& display(std::ostream& os) {
      manager->display(os, false);
      return GenericDVOI< VarComparator >::display(os, RAND, manager->get_weight());
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
  template <class WeightManager, template< class T > class Aggregator, class VarComparator, int RAND >
  class GenericNeighborDVO : public GenericWeightedRandomDVO< WeightManager, Aggregator< VarComparator >, RAND >
  {
  public: 
    
    /**@name Parameters*/
    //@{ 
    Solver *solver;
    Vector< Variable > *neighborhood;
    //@}

    /**@name Constructors*/
    //@{
    GenericNeighborDVO() : GenericWeightedRandomDVO< WeightManager, Aggregator< VarComparator >, RAND >() { neighborhood=NULL; }
    GenericNeighborDVO(Solver* s) : GenericWeightedRandomDVO< WeightManager, Aggregator< VarComparator >, RAND >(s) {
      neighborhood = NULL;
    }
    virtual ~GenericNeighborDVO() {
      delete [] neighborhood;
    }
    //virtual void initialise(Solver *s, void *a=NULL) { 
    virtual void initialise(Solver *s) { 
      GenericWeightedRandomDVO< WeightManager, Aggregator< VarComparator >, RAND >::initialise(s);
    }
    //virtual void initialise(Solver *s) { manager->initialise(s); }
    virtual void initialise(VarStack< Variable, ReversibleNum<int> >& seq) 
    {
      int i, j, k, cons, self_idx;
      Constraint constraint;
      Variable *scope;
      Event trig;
      bool is_in;

      //std::cout << "NEIGHBORHOOD " << neighborhood << std::endl;

      if(!neighborhood) {

	//std::cout << "NEIGHBORHOOD" << std::endl;

	neighborhood = new Vector< Variable >[GenericDVOI< Aggregator< VarComparator > >::solver->variables.size];
	std::fill(neighborhood, neighborhood+GenericDVOI< Aggregator< VarComparator > >::solver->variables.size, NULL);
	for(i=seq.size; --i>=0;) {

	  self_idx = seq[i].id();
	  for(trig = 0; trig<3; ++trig) {
	    
	    for(cons = GenericDVOI< Aggregator< VarComparator > >::solver->constraint_graph[self_idx].on[trig].size; --cons>=0;) {
	      
	      constraint = GenericDVOI< Aggregator< VarComparator > >::solver->constraint_graph[self_idx].on[trig][cons];

	      for(scope=constraint.get_scope(), j=constraint.arity(); --j>=0;) if(scope[j].id() != seq[i].id()) {
		  is_in = false;
		  for(k = neighborhood[self_idx].size; --k>=0 && !is_in;)
		    if(neighborhood[self_idx][k].id() == scope[j].id()) is_in = true;
		  if(!is_in) neighborhood[self_idx].add(scope[j]);
		}
	    }
	  }
	  //std::cout << seq[i] << ": " << neighborhood[self_idx] << std::endl;
	}
      }

      GenericRandomDVO< Aggregator< VarComparator >, RAND >::current.map = neighborhood;
      for(int i=0; i<=RAND; ++i)
	GenericRandomDVO< Aggregator< VarComparator >, RAND >::bests[i].map = neighborhood;
      
    }

    // virtual void initialise_manager() {
    //   manager = new WeightManager(GenericRandomDVO< Aggregator< VarComparator >, RAND >::solver);
    //   //GenericRandomDVO< VarComparator >::solver->add(manager);
    //   if(manager->get_weight()) {
    // 	GenericRandomDVO< Aggregator< VarComparator >, RAND >::current.weight = manager->get_weight();
    // 	for(int i=0; i<=RAND; ++i)
    // 	  GenericRandomDVO< Aggregator< VarComparator >, RAND >::bests[i].weight = manager->get_weight();
    //   }
    // }

    virtual std::ostream& display(std::ostream& os) {
      GenericWeightedRandomDVO< WeightManager, Aggregator< VarComparator >, RAND >::manager->display(os, false);
      return GenericDVOI< Aggregator< VarComparator > >::display(os, RAND, GenericWeightedRandomDVO< WeightManager, Aggregator< VarComparator >, RAND >::manager->get_weight());
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
    void initialise(VarStack< Variable, ReversibleNum<int> >& seq) {}
    
    Variable select();

    virtual std::ostream& display(std::ostream& os) const {
      os << "Go by the default sequence: " << solver->sequence.back(); 
      return os;
    }
  };


  /*! \class Lexicographic
    \brief  Class Lexicographic

    This heuristic selects the variable with lowest rank
    in the initial sequence of search variables.
  */
  class Lexicographic : public VariableListener {

  public: 
    
    Solver             *solver;
    Vector< Variable >   order;
    Vector< int >        index;
    ReversibleNum< int >  last;

    Lexicographic() { solver = NULL; }
    Lexicographic(Solver *s);
    void initialise(Solver *s);
    //void initialise(Solver *s, void *a=NULL);
    void initialise(VarStack< Variable, ReversibleNum<int> >& seq);
    virtual ~Lexicographic();

    virtual void notify_add_var() {};
    virtual void notify_change(const int idx);
    
    Variable select();

    virtual std::ostream& display(std::ostream& os) const {
      os << "Go by lexicographic order: " ;

      int i = last;
      while(i<(int)(order.size) && order[i].is_ground()) { 
	++i;
      }
      os << order[i];
      return os;
    }

  };


  /*! \class Lexicographic
    \brief  Class Lexicographic

    This heuristic selects the variable with lowest rank
    in the initial sequence of search variables.
  */
  class ConsolidateListener : public VariableListener, public ConstraintListener {

  public: 
    
    Solver                                      *solver;
    VarStack < Variable, ReversibleNum<int> > *sequence;
    Vector< Vector< Constraint > >          constraints;


    ConsolidateListener(Solver *s);
    virtual ~ConsolidateListener();

    virtual void notify_post   (Constraint c);
    virtual void notify_relax  (Constraint c);
    virtual void notify_add_con(Constraint c);
    virtual void notify_add_var();
    virtual void notify_change (const int idx);

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
    void initialise() {dom_ = LARGE_VALUE;}
    //@}

    /**@name Parameters*/
    //@{ 
    int dom_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (double)dom_; } 
    inline bool operator<( const MinDomain& x ) const { return dom_ < x.dom_; }

    // inline MinDomain& operator*( const int x ) const { MinDomain d; d.dom_=(dom_ * x); return d; }
    // inline MinDomain& operator+( const int x ) const { MinDomain d; d.dom_=(dom_ + x); return d; }
    // inline MinDomain& operator-( const int x ) const { MinDomain d; d.dom_=(dom_ - x); return d; }
    // inline MinDomain& operator/( const int x ) const { MinDomain d; d.dom_=(dom_ / x); return d; }

    inline MinDomain& operator*=( const int x ) { dom_ *= x; return *this; }
    inline MinDomain& operator+=( const int x ) { dom_ += x; return *this; }
    inline MinDomain& operator-=( const int x ) { dom_ -= x; return *this; }
    inline MinDomain& operator/=( const int x ) { dom_ /= x; return *this; }

    // inline MinDomain& operator*( const MinDomain& x ) const { MinDomain d; d.dom_=(dom_ * x.dom_); return d; }
    // inline MinDomain& operator+( const MinDomain& x ) const { MinDomain d; d.dom_=(dom_ + x.dom_); return d; }
    // inline MinDomain& operator-( const MinDomain& x ) const { MinDomain d; d.dom_=(dom_ - x.dom_); return d; }
    // inline MinDomain& operator/( const MinDomain& x ) const { MinDomain d; d.dom_=(dom_ / x.dom_); return d; }

    inline MinDomain& operator*=( const MinDomain& x ) { dom_ *= x.dom_; return *this; }
    inline MinDomain& operator+=( const MinDomain& x ) { dom_ += x.dom_; return *this; }
    inline MinDomain& operator-=( const MinDomain& x ) { dom_ -= x.dom_; return *this; }
    inline MinDomain& operator/=( const MinDomain& x ) { dom_ /= x.dom_; return *this; }

    inline void operator=( const MinDomain& x ) { dom_ = x.dom_; }
    inline void operator=( const Variable x ) { dom_ = x.get_size(); }

    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << " with minimum domain size";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
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
    void initialise() {dom_ = LARGE_VALUE; deg_ = 0;}
    //@}

    /**@name Parameters*/
    //@{ 
    int dom_;
    int deg_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (deg_ ? (double)dom_/(double)deg_ : (double)INFTY); } 
    inline bool operator<( const MinDomainOverDegree& x ) const { return dom_*x.deg_ < x.dom_*deg_; }

    // inline MinDomainOverDegree& operator*( const MinDomainOverDegree& x ) const { MinDomainOverDegree d; d.dom_=(dom_ * x.dom_); d.deg_=(deg_ * x.deg_); return d; }
    // inline MinDomainOverDegree& operator+( const MinDomainOverDegree& x ) const { MinDomainOverDegree d; d.dom_=(dom_ + x.dom_); d.deg_=(deg_ + x.deg_); return d; }
    // inline MinDomainOverDegree& operator-( const MinDomainOverDegree& x ) const { MinDomainOverDegree d; d.dom_=(dom_ - x.dom_); d.deg_=(deg_ - x.deg_); return d; }
    // inline MinDomainOverDegree& operator/( const MinDomainOverDegree& x ) const { MinDomainOverDegree d; d.dom_=(dom_ / x.dom_); d.deg_=(deg_ / x.deg_); return d; }

  

    inline MinDomainOverDegree& operator*=( const int x ) { dom_ *= x; return *this; }
    inline MinDomainOverDegree& operator+=( const int x ) { dom_ += x * deg_; return *this; }
    inline MinDomainOverDegree& operator-=( const int x ) { dom_ -= x * deg_; return *this; }
    inline MinDomainOverDegree& operator/=( const int x ) { deg_ *= x; return *this; }


    inline MinDomainOverDegree& operator*=( const MinDomainOverDegree& x ) { dom_ *= x.dom_; deg_ *= x.deg_; return *this; }
    inline MinDomainOverDegree& operator+=( const MinDomainOverDegree& x ) { 
      dom_ = (dom_ * x.deg_ + x.dom_ * deg_); 
      deg_ *= x.deg_; return *this; }
    inline MinDomainOverDegree& operator-=( const MinDomainOverDegree& x ) { 
      dom_ = (dom_ * x.deg_ - x.dom_ * deg_); 
      deg_ *= x.deg_; return *this; }
    inline MinDomainOverDegree& operator/=( const MinDomainOverDegree& x ) { dom_ /= x.dom_; deg_ /= x.deg_; return *this; }

    // inline MinDomainOverDegree& operator*=( const int x ) { dom_ *= x; deg_ *= x; return *this; }
    // inline MinDomainOverDegree& operator+=( const int x ) { dom_ += x; deg_ += x; return *this; }
    // inline MinDomainOverDegree& operator-=( const int x ) { dom_ -= x; deg_ -= x; return *this; }
    // inline MinDomainOverDegree& operator/=( const int x ) { dom_ /= x; deg_ /= x; return *this; }


    // inline MinDomainOverDegree& operator*=( const MinDomainOverDegree& x ) { dom_ *= x.dom_; deg_ *= x.deg_; return *this; }
    // inline MinDomainOverDegree& operator+=( const MinDomainOverDegree& x ) { dom_ += x.dom_; deg_ += x.deg_; return *this; }
    // inline MinDomainOverDegree& operator-=( const MinDomainOverDegree& x ) { dom_ -= x.dom_; deg_ -= x.deg_; return *this; }
    // inline MinDomainOverDegree& operator/=( const MinDomainOverDegree& x ) { dom_ /= x.dom_; deg_ /= x.deg_; return *this; }

    inline void operator=( const MinDomainOverDegree& x ) { dom_ = x.dom_; deg_ = x.deg_; }
    inline void operator=( const Variable x ) { dom_ = x.get_size(); deg_ = x.get_degree(); }
    //@}  


    std::ostream& display_criterion(std::ostream& os) const {
      os << " with minimum (domain size / degree)";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
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
    void initialise() {dom_ = LARGE_VALUE; wei_ = 0;}
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
    inline bool operator<( const MinDomainOverWeight& x ) const { return dom_*x.wei_ < x.dom_*wei_; }


    // inline MinDomainOverWeight& operator*( const MinDomainOverWeight& x ) const { MinDomainOverWeight d; d.dom_=(dom_ * x.dom_); d.wei_=(wei_ * x.wei_); return d; }
    // inline MinDomainOverWeight& operator+( const MinDomainOverWeight& x ) const { MinDomainOverWeight d; d.dom_=(dom_ + x.dom_); d.wei_=(wei_ + x.wei_); return d; }
    // inline MinDomainOverWeight& operator-( const MinDomainOverWeight& x ) const { MinDomainOverWeight d; d.dom_=(dom_ - x.dom_); d.wei_=(wei_ - x.wei_); return d; }
    // inline MinDomainOverWeight& operator/( const MinDomainOverWeight& x ) const { MinDomainOverWeight d; d.dom_=(dom_ / x.dom_); d.wei_=(wei_ / x.wei_); return d; }


    inline void reduce() {
      while(dom_>100000 || !(dom_&1)) {
	dom_ /=2;
	wei_ /=2;
      }
    }


    inline MinDomainOverWeight& operator*=( const int x ) { dom_ *= x; return *this; }
    inline MinDomainOverWeight& operator+=( const int x ) { dom_ += x * wei_; return *this; }
    inline MinDomainOverWeight& operator-=( const int x ) { dom_ -= x * wei_; return *this; }
    inline MinDomainOverWeight& operator/=( const int x ) { wei_ *= x; return *this; }


    inline MinDomainOverWeight& operator*=( const MinDomainOverWeight& x ) { dom_ *= x.dom_; wei_ *= x.wei_; return *this; }
    inline MinDomainOverWeight& operator+=( const MinDomainOverWeight& x ) { 
      dom_ = (dom_ * x.wei_ + x.dom_ * wei_); 
      wei_ *= x.wei_; 
      reduce();
      return *this; }
    inline MinDomainOverWeight& operator-=( const MinDomainOverWeight& x ) { 
      dom_ = (dom_ * x.wei_ - x.dom_ * wei_); 
      wei_ *= x.wei_; return *this; }
    inline MinDomainOverWeight& operator/=( const MinDomainOverWeight& x ) { dom_ /= x.dom_; wei_ /= x.wei_; return *this; }


    inline void operator=( const MinDomainOverWeight& x ) { dom_ = x.dom_; wei_ = x.wei_; }
    inline void operator=( const Variable x ) { 
      dom_ = x.get_size(); wei_ = weight[x.id()]; 
      //std::cout << x << ": " << dom_ << "/" << wei_ << std::endl;
    }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << " with minimum (domain size / weight)";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << dom_ << "/" << wei_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MinDomainOverWeight& x);

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
    void initialise() {dom_ = LARGE_VALUE; wei_ = 0;}
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
    inline bool operator<( const MinNeighborDomainOverNeighborWeight& x ) const { return dom_*x.wei_ < x.dom_*wei_; }
    inline void operator=( const MinNeighborDomainOverNeighborWeight& x ) { dom_ = x.dom_; wei_ = x.wei_; }
    inline void operator=( const Variable x ) { 
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

    std::ostream& display_criterion(std::ostream& os) const {
      os << " with minimum (Sum of neighrbors' domain sizes / Sum of neighrbors' weights)";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << dom_ << "/" << wei_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MinNeighborDomainOverNeighborWeight& x);


  /*! \class MinNeighborDomainOverNeighborWeight
    \brief  Class MinNeighborDomainOverNeighborWeight

    Order two variables by the ratio of the domain sizes and weight
    of a vector of other variables
  */
  template< class VarComparator >
  class SelfPlusAverage
  {
  public: 

    /**@name Constructors*/
    //@{
    SelfPlusAverage() { map = NULL; }
    void initialise() { crit.initialise(); map = NULL; }
    //@}

    /**@name Parameters*/
    //@{ 
    double *weight;
    VarComparator crit;
    Vector< Variable > *map;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return 1; } 
    inline bool operator<( const SelfPlusAverage<VarComparator> & x ) const { return crit < x.crit; }
    inline void operator=( const SelfPlusAverage<VarComparator>& x ) { crit = x.crit; }
    inline void operator=( const Variable x ) { 



      crit.weight = weight;
      crit = x;


      // std::cout << "\n => ";
      // crit.display(std::cout);
      // std::cout << std::endl;
      // std::cout << "scan neighborhood" << std::endl;



      VarComparator crit_neighbor, aux;
      crit_neighbor.weight = weight;
      aux.weight = weight;

      int idx = x.id();

      int n = map[idx].size;

      int i = n-1;
      Variable y = map[idx][i];
      crit_neighbor = y;

      // std::cout << "   ---> ";
      // crit_neighbor.display(std::cout);
      // std::cout << std::endl;

      while(--i>=0) {
  	y = map[idx][i];
	aux = y;
	crit_neighbor += aux;

	// std::cout << "      + ";
	// aux.display(std::cout);
	// std::cout << " = " ;
	// crit_neighbor.display(std::cout);
	// std::cout << std::endl;
      } 
      crit_neighbor /= n;

      // std::cout << "   ---> ";
      //  crit_neighbor.display(std::cout);
      // std::cout << std::endl;
      
      crit += crit_neighbor;


      // std::cout << "   ---> " ;
      //  crit.display(std::cout);
      // std::cout << std::endl;


      //std::cout << x << ": " << dom_ << "/" << wei_ << std::endl;

    }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << " sum of neighbors' " ;
      crit.display_criterion(os);
      os << " plus self" ;
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      //os << crit;
      crit.display(os);
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MinNeighborDomainOverNeighborWeight& x);


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
    void initialise() {dom_ = LARGE_VALUE; wei_ = 0;}
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
    inline bool operator<( const MinNeighborDomainOverWeight& x ) const { return dom_*x.wei_ < x.dom_*wei_; }
    inline void operator=( const MinNeighborDomainOverWeight& x ) { dom_ = x.dom_; wei_ = x.wei_; }
    inline void operator=( const Variable x ) { 
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

    std::ostream& display_criterion(std::ostream& os) const {
      os << " with minimum (Sum of neighrbors' domain sizes / weight)";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << dom_ << "/" << wei_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MinNeighborDomainOverWeight& x);


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
    void initialise() {wei_ = 0;}
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
    inline bool operator<( const MaxWeight& x ) const { return x.wei_ < wei_; }
    inline void operator=( const MaxWeight& x ) { wei_ = x.wei_; }
    inline void operator=( const Variable x ) { 
      //std::cout << "activity of " << x << " = " << weight[x.id()] << std::endl;
      wei_ = weight[x.id()]; 
    }
    //@}  


    std::ostream& display_criterion(std::ostream& os) const {
      os << " with maximum weight";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
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

     std::ostream& display(std::ostream& os) const {
       os << "assign it to the minimum value in its domain";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, MinValue& x);

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

     std::ostream& display(std::ostream& os) const {
       os << "assign it to the maximum value in its domain";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, MaxValue& x);

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

     std::ostream& display(std::ostream& os) const {
       os << "set its upper bound to (min+max)/2";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, HalfSplit& x);

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

     std::ostream& display(std::ostream& os) const {
       os << "assign it to its minimum or maximum value at random with equal probability";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, RandomMinMax& x);

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

     std::ostream& display(std::ostream& os) const {
       os << "assign it to the value with minimum weight in its domain";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, MinWeightValue& x);


  /*! \class Guided
    \brief  Class Guided

    Assigns the variable to its minimum value.
  */
  class Guided {
    
  public: 
    
    Solver *solver;
    
    Guided() {}
    Guided(Solver *s) {solver=s;}
    virtual ~Guided() {};
    
    inline Decision make(Variable x) {
      int val = solver->last_solution_lb[x.id()];
      Decision d(x, Decision::ASSIGNMENT, val);
      if(!x.contain(val)) d.set_value(x.get_min());
      return d;
    }


     std::ostream& display(std::ostream& os) const {
       os << "assign it to the value of this variable in the last solution";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, Guided& x);


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

     std::ostream& display(std::ostream& os) const {
       os << "assign it to the value with minimum weight in its (Boolean) domain";
       return os;
     }
        
  };
  
  std::ostream& operator<<(std::ostream& os, BoolMinWeightValue& x);






  template< class VarComparator >
  std::ostream& operator<<(std::ostream& os, Identifiable<VarComparator>& x) {
    x.criterion.display(os);
    return os;
  }

  
  
  class VSIDS : public GenericHeuristic< GenericWeightedDVO< LiteralActivityManager, MaxWeight >, 
					 //MinValue > {
					 BoolMinWeightValue > {
  public:
    
    VSIDS(Solver *s) : GenericHeuristic< GenericWeightedDVO< LiteralActivityManager, MaxWeight >, 
					 //MinValue >(s) {
					 BoolMinWeightValue >(s) {
      choice.weight = s->get_literal_activity();
    }
  };



}

#endif // __SEARCH_HPP
