
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


/*! \file mistral_solver.hpp
    \brief Header file for the solver.
*/


#ifndef __SOLVER_HPP
#define __SOLVER_HPP


//#include <mistral_search.hpp>
#include <mistral_global.hpp>
#include <mistral_backtrack.hpp>
#include <mistral_structure.hpp>
#include <mistral_constraint.hpp>


namespace Mistral {

  class SolverParameters {

  public:

    SolverParameters();
    SolverParameters(const SolverParameters&);
    virtual ~SolverParameters();
    void initialise();
    void copy(const SolverParameters&);

    /// level of verbosity
    int verbosity;
    /// Number of solutions to find, if equal to -1, then all solutions are listed
    int find_all;

    /// Limit on the number of nodes
    unsigned int node_limit;
    /// Limit on the number of backtracks
    unsigned int backtrack_limit;
    /// Limit on the number of failures
    unsigned int fail_limit;
    /// fail limit used for restarts
    unsigned int restart_limit;
    /// flag to know if there is a limit
    unsigned int limit;

    /// Limit on cpu time
    double time_limit;  
  };

  class SolverStatistics {
  
  public:

    SolverStatistics();
    SolverStatistics(const SolverStatistics&);
    virtual ~SolverStatistics();
    void initialise();
    void copy(const SolverStatistics&);
    void update(const SolverStatistics&);

    /// Number of nodes, that is recursive calls to  the dfs algo
    unsigned long int num_nodes; 
    /// Number of backtracks, that is unsuccesful recursive calls 
    unsigned long int num_backtracks;
    /// Number of constraint failures
    unsigned long int num_failures; 
    /// Number of constraint failures
    unsigned long int num_restarts; 
    /// Number of calls to a constraint propagator
    unsigned long int num_propagations;
    /// Number of solutions found so far
    unsigned long int num_solutions;
    /// Number of inference steps
    unsigned long int num_filterings;
    /// Search outcome
    Outcome outcome;
    
    /// Number of nodes, that is recursive calls to  the dfs algo
    unsigned long int num_variables; 
    /// Number of nodes, that is recursive calls to  the dfs algo
    unsigned long int num_values; 
    /// Number of nodes, that is recursive calls to  the dfs algo
    unsigned long int num_constraints; 

    /// timestamp
    double start_time;
    /// timestamp
    double end_time;

    //virtual std::string getString() const;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::ostream& print_full(std::ostream&) const ;
  };
  


  class Solver;
  /**********************************************
   * MultiQueue
   **********************************************/
  
  /*! \class MultiQueue
    \brief MultiQueue Class
  */
  class ConstraintQueue // : public Vector<Constraint*>
  {

  private:
    Solver *solver;
    int min_priority;
    /// The constraint being processed
    Constraint *taboo_constraint;
    //Constraint **stack_;
    //int *min_index;
    //int *max_index;

  public:

    int      cardinality;
    Queue      *triggers;
    int  higher_priority;

    BitSet _set_;
    
    ConstraintQueue();
    void initialise(Solver *s);
    void initialise(const int min_p, const int max_p, const int size);
    void declare(Constraint *c, Solver *s);
    virtual ~ConstraintQueue();
    inline bool empty() { return higher_priority<min_priority; }
    void trigger(Constraint* cons);
    void trigger(Constraint* cons, const int var, const Event evt);
//     inline void trigger(Constraint* cons, const int var, const Event evt)
//     //void Mistral::ConstraintQueue::trigger(Constraint* cons, const int var, const Event evt)//;
// {

// #ifdef _DEBUG_AC
//   std::cout << "  triggers " << cons << "(" << (cons->id) << ") by " 
// 	    << cons->scope[var] << std::endl;
// #endif


//   if(cons != taboo_constraint) {
//     int priority = cons->priority, cons_id = cons->id;
//     if(_set_.fastContain(cons_id)) {

// // #ifdef _DEBUG_AC
// //       std::cout << "update " << std::endl;
// // #endif

//       if(cons->events.contain(var)) {
// 	cons->event_type[var] |= evt;
//       } else {
// 	cons->events.add(var);
// 	cons->event_type[var] = evt;
//       }
//     } else {

// // #ifdef _DEBUG_AC
// //       std::cout << "add " << std::endl;
// // #endif

//       _set_.fastAdd(cons_id);
//       if(priority > higher_priority) higher_priority = priority;
//       triggers[priority].add(cons_id);
//       cons->events.setTo(var);
//       //cons->events.clear();
//       //cons->events.add(var);
//       cons->event_type[var] = evt;
//     }
//   } 

//  #ifdef _DEBUG_AC
//   else
//     std::cout << cons << "(" << (cons->id) 
// 	      << ") is currently being processed " 
// 	      << std::endl;
// #endif

//  #ifdef _DEBUG_AC
//   std::cout << " => " << cons->events << std::endl;
// //   else
// //     std::cout << cons << "(" << (cons->id) 
// // 	      << ") is currently being processed " 
// // 	      << std::endl;
// #endif

// }

    inline void reset_higher_priority() {
      while(--higher_priority>=min_priority && triggers[higher_priority].empty());
    }
    inline Constraint* select(Vector<Constraint*>& constraints) {
      int cons_id = triggers[higher_priority].pop();
      Constraint *cons = constraints[cons_id];
      _set_.fastErase(cons_id);
      if(triggers[higher_priority].empty()) reset_higher_priority();


//     std::cout << "events of " << cons << " before freeze: " 
// 	      << cons->events << std::endl;
//     std::cout << "changes of " << cons << " before freeze: " 
// 	      << cons->changes << std::endl;


      taboo_constraint = cons->freeze();
      return cons;
    }
    inline void clear() {
      while(higher_priority>=min_priority) triggers[higher_priority--].clear();
      _set_.clear();
      taboo_constraint = NULL;
    }
    virtual std::ostream& display(std::ostream&) ;
  };



  /**********************************************
   * Solver
   **********************************************/

  /*! \class Solver
    \brief Solver Class

    The solver class do the following:
    * it keeps a list of the constraints and variables 
    * it controls 
    - the propagation queue
    - the bracktracking structures
    - the search stack
  */

  class BranchingHeuristic;
  class RestartPolicy;
  class Reversible;
  class Decision;
  class Variable;
  class VariableImplementation;
  class Solver : public Environment {

  public:
    /*!@name Parameters*/
    //@{
    /// The current level in the search tree (number of left branches from the root)
    /*! The current level in the search tree (number of left branches from the root)
      level -1: state of the domains as stated when modelling
      level 0: level at wich constraints are added (+ initial propagation)
      level k: after k decisions
    */
    //int level;

    /// The set of variables, in the initial order, that is as loaded from the model
    Vector< Variable > variables;

    /// The set of constraints, with special accessors for triggers
    Vector< Constraint* >         constraints;
    ConstraintQueue        active_constraints;
    
    Vector< ConstraintList* > constraint_graph;
    Vector< Vector <ConstraintTrigger> > triggers;

    /// For each level, the list of reversible objects that changed at this level, 
    /// and will need to be restored
    //Vector< Reversible* > saved_objs;
    /// For each level, the list of reversible objects that changed at this level, 
    /// and will need to be restored
    Vector< Variable > saved_vars;
    /// For each level, the list of reversible objects that changed at this level, 
    /// and will need to be restored
    Vector< Constraint* > saved_cons;

    /// Stores the last solution
    Vector< int > last_solution_lb;
    Vector< int > last_solution_ub;

    /// Set of parameters for the solver
    SolverParameters parameters;
    /// Set of statistics collected by the solver
    SolverStatistics statistics;
    //@}

    /// The search part
    /// These are the search variables.
    VarStack < Variable > sequence;
    //IntStack sequence;
    Vector< Decision > decisions;

    /// The delimitation between different levels is kept by this vector of integers
    Vector< int > trail_;

    BranchingHeuristic *heuristic;
    RestartPolicy *policy; 

    /// the constraint responsible for the last fail (NULL otherwise) 
    Constraint *culprit;
    /// the constraint-index of the last wiped_out variable in the culprit constraint (-1 otherwise)
    int wiped_idx;
    /// the overall-index of the variable responsible for the last trigger before a failure
    int wiper_idx;

    Outcome satisfied();

    /*!@name Constructor*/
    //@{
    Solver();
    //void initialise();


    class BooleanMemoryManager {
    public:
      Vector<int*> slots;

      BooleanMemoryManager() {}
      virtual ~BooleanMemoryManager() {
	while(!slots.empty()) 
	  delete [] slots.pop();
      }

      int* get_next() { 
	int *dom = new int[1];
	*dom = 3;
	slots.add(dom);
	return dom;
      }
    };

    //int *booleans;
    BooleanMemoryManager booleans;

    int* getNextBooleanSlot() { return booleans.get_next(); }
    //@}

    /*!@name Destructor*/
    //@{
    virtual ~Solver();
    //@}
  
    /*!@name Declarations*/
    //@{
    /// add a variable (prior to search!!!)
    void add(Variable x);
    int declare(Variable x);
    /// add a constraint
    void add(Constraint* x); 
    //@}

    /*!@name Trail accessors*/
    //@{
    inline void save() {
      trail_.add(sequence.size);
      trail_.add(saved_objs.size);
      trail_.add(saved_cons.size);
      trail_.add(saved_vars.size);

      ++statistics.num_nodes;
      ++level;
    }
    /// called by reversible objects when they want to be restored when backtrack to this level
    inline void save(Reversible* object) { saved_objs.add(object); }
    inline void save(Constraint* c) { saved_cons.add(c); }
    void save(Variable x); 
    void save(VariableImplementation *x, int dtype); 
    //@}

    /*!@name Propagation accessors*/
    //@{
    void notify_failure();
    void notify_success();
    void notify_decision();
    /// called when the var'th variable of constraint cons changes (with event type evt)
    //void trigger(Constraint* cons, const int var, const Event evt);
    //void activate_constraint(Constraint* cons, const int var, const Event evt);
    

    void trigger_event(const int var, const Event evt);
//     inline void trigger_event(const int var, const Event evt) 
// {
//  #ifdef _DEBUG_AC
//   std::cout << (ASSIGNED(evt) ? "value" : (BOUND_CHANGED(evt) ? "range" : "domain"))
// 	    << " event on " << variables[var] << std::endl;
// #endif

//   Constraint *c;
//   ConstraintNode nd;

//   nd = constraint_graph[var]->first(EVENT_TYPE(evt));

//   if(ASSIGNED(evt)) {
//     while(constraint_graph[var]->next(nd)) {
//       active_constraints.trigger(nd.elt.constraint, nd.elt.index, evt);
//       c = nd.elt.constraint->notify_assignment(nd.elt.index, level);
//       if(c) saved_cons.add(c);
//     }
//     if(sequence.contain(variables[var])) sequence.remove(variables[var]);
//   } else while(constraint_graph[var]->next(nd)) {
//       active_constraints.trigger(nd.elt.constraint, nd.elt.index, evt);
//     }

//   wiper_idx = var;
  
// }

    /// achieve propagation closure
    bool propagate(); 
    //@}

    /*!@name Search accesors*/
    //@{
    /// make a branching decision
    void branch_left();
    /// undo the previous branching decision and enforce the complementary constraint
    void branch_right();


    /// backtrack one level
    void restore();
    /// backtrack to a given level
    void restore(const int);


    void initialise_search(Vector< Variable >& seq, 
			   BranchingHeuristic *heu=NULL, 
			   RestartPolicy *pol=NULL);
    Outcome get_next_solution();

    /*!
      Launches a depth first search on the sequence 'seq' of variables
      with variable ordering 'heu' and restart policy 'pol'.
      Returns an outcome (SAT/UNSAT/OPT/UNKNOWN) and undo all decisions before returning.
    */
    Outcome depth_first_search(Vector< Variable >& seq, 
			       BranchingHeuristic *heu=NULL, 
			       RestartPolicy *pol=NULL);
    /*!
      Launches a depth first search on all variables
      with variable ordering 'heu' and restart policy 'pol'.
      Returns an outcome (SAT/UNSAT/OPT/UNKNOWN) and undo all decisions before returning.
    */
    Outcome depth_first_search(BranchingHeuristic *heu=NULL, 
			       RestartPolicy *pol=NULL);
    ///
    bool limits_expired();

    /// depth first search algorithm
    Outcome iterative_dfs();
//     //Mistral::Outcome Mistral::Solver::iterative_dfs() 
// {
//   int status = UNKNOWN;
//   while(status == UNKNOWN) {
//     if(propagate()) {
//       if( sequence.empty() ) status = satisfied();
//       else branch_left();
//     } else {
//       if( decisions.empty() ) status = UNSAT;
//       else if( limits_expired() ) status = LIMITOUT;
//       else branch_right();
//     }
//   }
//   return status;
// }
    //@}


    /*!@name Search accesors*/
    //@{
    void debug_print();
    void full_print();
    virtual std::ostream& display(std::ostream&) ;
    //@}
  };


  std::ostream& operator<< (std::ostream& os, Solver& x);
  std::ostream& operator<< (std::ostream& os, Solver* x);

  std::ostream& operator<< (std::ostream& os, ConstraintQueue& x);
  std::ostream& operator<< (std::ostream& os, ConstraintQueue* x);

  std::ostream& operator<< (std::ostream& os, const SolverStatistics& x);
  std::ostream& operator<< (std::ostream& os, const SolverStatistics* x);

}

#endif // __SOLVER_HPP
