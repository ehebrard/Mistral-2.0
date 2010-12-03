
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


#include <mistral_search.hpp>
#include <mistral_global.hpp>
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

  class Reversible;
  template< int N >
  class MultiQueue;
  typedef MultiQueue< 3 > ACQueue;
  class Variable;
  class Solver {

  public:
    /*!@name Parameters*/
    //@{
    /// The current level in the search tree (number of left branches from the root)
    /*! The current level in the search tree (number of left branches from the root)
      level -1: state of the domains as stated when modelling
      level 0: level at wich constraints are added (+ initial propagation)
      level k: after k decisions
    */
    int level;

    /// The set of variables, in the initial order, that is as loaded from the model
    //#ifdef _STATIC_CAST
    Vector< Variable > variables_2;
    Vector< Variable > boolean_variables;
    //#else
    Vector< IntVar > variables;
    //#endif

    /// The set of constraints
    Vector< Constraint* > constraints;

    /// The constraints queue, used for achieving the propagation closure
    ACQueue triggered_constraints;
    /// The constraint being processed
    Constraint *taboo_constraint;


    /// For each level, the list of reversible objects that changed at this level, and will need to be restored
    Vector< Variable > saved_variables;
    Vector< Reversible* > saved_objects;
    /// The delimitation between different levels is kept by this vector of integers
    Vector< int > obj_trail_size;

    /// For each level, the list of reversible objects that changed at this level, and will need to be restored
    Vector< IntVar > saved_vars;
    /// The delimitation between different levels is kept by this vector of integers
    Vector< int > var_trail_size;


    /// For each level, the list of reversible objects that changed at this level, and will need to be restored
    Vector< Constraint* > saved_cons;
    /// The delimitation between different levels is kept by this vector of integers
    Vector< int > con_trail_size;

    /// There is only one goal at a time, but one can solve several goals back to back
    /// (or force backtracking to start back from a previous state) externally
    ///Search search;

    /// Stores the last solution
    Vector< int > solution;

    /// Set of parameters for the solver
    SolverParameters parameters;
    /// Set of statistics collected by the solver
    SolverStatistics statistics;
    //@}

    /// The search part
    /// These are the search variables.
    Stack< IntVar > sequence;
    /// These are the variables we do not want to branch on.
    Stack< IntVar > auxilliary;




    Vector< unsigned int > trail_seq;
    Vector< unsigned int > trail_aux;
    Vector< IntVar > decision;

    VarOrdering *heuristic;   
    RestartPolicy *policy;   

    //void assign(IntVar x);
    Outcome satisfied();

    /*!@name Constructor*/
    //@{
    Solver();
    void initialise();
    //@}

    /*!@name Destructor*/
    //@{
    virtual ~Solver();
    //@}
  
    /*!@name Declarations*/
    //@{
    /// add a variable (prior to search!!!)
    void add(IntVar x);
    /// add a constraint
    void add(Constraint* x); 



    void add(Variable x);
    //@}

    /*!@name Trail accessors*/
    //@{
    /// called by reversible objects when they want to be restored when backtrack to this level
    inline void save(Reversible* object) { saved_objects.add(object); }
    inline void save(IntVar x) { saved_vars.add(x); }
    inline void save(Constraint* c) { saved_cons.add(c); }
//     inline void save(VariableImplementation* v, int type) { 
//       Variable x(v,type);
//       std::cout << "saving " << x << std::endl;
//     }
    void save(Variable x); // { saved_variables.add(x); }
    //@}

    /*!@name Propagation accessors*/
    //@{
    /// called when the var'th variable of constraint cons changes (with event type evt)
    void trigger(Constraint* cons, const int var, const Event evt);
    void triggerEvent(const int var, const Event evt);
    /// achieve propagation closure
    bool propagate(); 
    //@}

    /*!@name Search accesors*/
    //@{
    /// make a branching decision
    void branchLeft();
    /// undo the previous branching decision and enforce the complementary constraint
    void branchRight();
    /// do the necessary setup of the trailling structures in order to advance in search
    void make_node();
    /// backtrack one level
    void backtrack();
    /// backtrack to a given level
    void backtrack(const int&);
    /*!
      Launches a depth first search on the sequence 'seq' of variables
      with variable ordering 'heu' and restart policy 'pol'.
      Returns an outcome (SAT/UNSAT/OPT/UNKNOWN) and undo all decisions before returning.
    */
    Outcome depth_first_search(Vector< IntVar >& seq, 
			       VarOrdering *heu=NULL, 
			       RestartPolicy *pol=NULL);
    /*!
      Launches a depth first search on all variables
      with variable ordering 'heu' and restart policy 'pol'.
      Returns an outcome (SAT/UNSAT/OPT/UNKNOWN) and undo all decisions before returning.
    */
    Outcome depth_first_search(VarOrdering *heu=NULL, 
			       RestartPolicy *pol=NULL);
    ///
    bool limitsExpired();
    /// depth first search algorithm
    
    Outcome iterative_dfs();
    
//     inline Outcome iterative_dfs() 
//     {
//       int status = UNKNOWN;
//       while(status == UNKNOWN) {
// 	if(propagate()) {

// 	  std::cout << (*this) << std::endl;

// 	  if( sequence.empty() ) status = satisfied();
// 	  else branchLeft();
// 	} else {
// 	  if( decision.empty() ) status = UNSAT;
// 	  else if( limitsExpired() ) status = LIMITOUT;
// 	  else branchRight();
// 	}
//       }
//       return status;
//     }
    //@}


    /*!@name Search accesors*/
    //@{
    void debug_print();
    void full_print();
    virtual std::ostream& display(std::ostream&) const ;
    //@}
  };


  std::ostream& operator<< (std::ostream& os, const Solver& x);
  std::ostream& operator<< (std::ostream& os, const Solver* x);

}

#endif // __SOLVER_HPP
