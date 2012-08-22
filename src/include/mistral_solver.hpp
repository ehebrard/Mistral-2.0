
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


#include <mistral_constraint.hpp>




namespace Mistral {




  class Solution {
    
  public: 
    
    Vector< Variable > variables;
    int min_id;
    int max_id;
    int *values;

    Solution( Vector< Variable >& vars );
    virtual ~Solution();

    virtual std::ostream& display(std::ostream&) const ;
    int& operator[](Variable x) ;

  };


  class SolverParameters {

  public:

    SolverParameters();
    SolverParameters(const SolverParameters&);
    virtual ~SolverParameters();
    void initialise();
    void copy(const SolverParameters&);

    /// random seed
    int seed;

    /// degree of verbosity
    int verbosity;
    /// Number of solutions to find, if equal to -1, then all solutions are listed
    int find_all;

    /// Limit on the number of nodes
    unsigned int node_limit;
    /// Limit on the number of backtracks
    unsigned int backtrack_limit;
    /// Limit on the number of failures
    unsigned int fail_limit;
    /// Limit on the number of failures
    unsigned int propagation_limit;
    /// fail limit used for restarts
    unsigned int restart_limit;
    /// flag to know if there is a limit
    unsigned int limit;

    /// Limit on cpu time
    double time_limit;  

    /// restart policy settings 
    int restart_policy;
    /// restart policy settings (base limit)
    unsigned int restart_base;
    /// restart policy settings (geometric increment)
    double restart_factor;

    /// type of randomization
    unsigned int randomization;
    /// variables sequence shuffle between restarts
    bool shuffle;


    int backjump;

    /// whether solutions are checked
    // 0 -> not checked
    // 1 -> check constraints which variables are all assigned
    // 2 -> 1+check that there exist a support for other constraints
    // 3 -> 2+check that all values are consistent for other constraints
    // NOT IMPLEMENTED! 4 -> check that the solution can be extended to all variables
    int checked; 
     


    /////// PARAMETERS FOR SAT SOLVING ///////
    double            activity_increment;
    double            normalize_activity;
    int               init_activity;
    double            forgetfulness;
    double            activity_decay;

    int               value_selection;
    int               dynamic_value;

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
    /// Objective value (ub for minimization, lb for maximization, -1 otherwise)
    long int objective_value;

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

    //unsigned int memory;

    /////// STATISTICS FOR SAT SOLVING ///////
    unsigned int      literals;
    unsigned int      small;
    double            base_avg_size;
    double            learnt_avg_size;


#ifdef _PROFILING

    void init_prof();

    double total_propag_time;
    //std::vector< double > constraint_propag_time;
    
    double total_branching_time;

    double total_restore_time;

    //int vartype_index[17];

#ifdef _PROFILING_PRIMITIVE

    double prof_time[NUM_METHODS][NUM_VARTYPES];
    unsigned long long int prof_num[NUM_METHODS][NUM_VARTYPES];
    virtual std::ostream& print_profile(std::ostream&) const ;

#endif

#endif



    //virtual std::string getString() const;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::ostream& print_full(std::ostream&) const ;
    virtual std::ostream& print_short(std::ostream&) const ;
  };
  

  class Solver;
  class SearchMonitor {
    
  public:
    
    SearchMonitor(Solver* s) { solver = s; }
    virtual ~SearchMonitor() {}

    Solver *solver;

    Vector<int> sequence;
    
    std::vector<const char*> strs;

    void add(Variable);
    void add(Constraint);
    void add(const char*);

    virtual std::ostream& display( std::ostream& os ) const;

  };


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
    //ConstraintImplementation* taboo_constraint;

  public:

    int      cardinality;
    Queue      *triggers;
    int  higher_priority;

    BitSet _set_;
    
    ConstraintQueue();
    void initialise(Solver *s);
    void initialise(const int min_p, const int max_p, const int size);
    void declare(Constraint c, Solver *s);
    virtual ~ConstraintQueue();
    inline bool empty() { return higher_priority<min_priority; }
    
    void add( Constraint cons );
    void add( ConstraintImplementation *cons);
    void trigger( BinaryConstraint *cons);
    void trigger(TernaryConstraint *cons);
    void trigger( GlobalConstraint *cons);
    void trigger( GlobalConstraint *cons, const int var, const Event evt);
    /*
    void trigger(Constraint cons);
    void trigger(Constraint cons, const Event evt);
    */

    void reset_higher_priority();
    // inline void reset_higher_priority() {
    //   while(--higher_priority>=min_priority && triggers[higher_priority].empty());
    // }

    Constraint select(Vector<Constraint>& constraints);
    // inline Constraint select(Vector<Constraint>& constraints) {
    //   int cons_id = triggers[higher_priority].pop();
    //   Constraint cons = constraints[cons_id];
    //   _set_.fast_remove(cons_id);
    //   if(triggers[higher_priority].empty()) reset_higher_priority();
    //   //taboo_constraint = cons.freeze();
    //   return cons;
    // }
    // // inline void select(Constraint cons) {
    // //   //_set_.remove(cons->id);
    // //   taboo_constraint = cons.freeze();
    // // }

    void clear();
    // inline void clear() {
    //   while(higher_priority>=min_priority) triggers[higher_priority--].clear();
    //   _set_.clear();
    //   //taboo_constraint = NULL;
    // }
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
  //typedef Varstack< ConstraintElement > ConstraintStack;

  class ConsolidateListener;
  class SolutionListener;
  class RestartListener;
  class SuccessListener;
  class FailureListener;
  class DecisionListener;
  class VariableListener;
  class ConstraintListener;
  class BranchingHeuristic;
  class RestartPolicy;
  class Reversible;
  class Expression;
  class Decision;
  class Variable;
  class VarArray;
  class Goal;
  class VariableImplementation;
  class ConstraintClauseBase;
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
    bool     search_started;
    int      backtrack_level;
    //Decision deduction;

    /// The set of variables, in the initial order, that is as loaded from the model
    Vector< Variable >   variables;
    /// The set of variables, as originally declared (i.e., as expressions)
    Vector< Variable >   declared_variables;
    Vector< int >     domain_types;
    Vector< int > assignment_level;

#ifdef _MONITOR
    SearchMonitor monitor_list;

    Vector< unsigned int > monitored;
    Vector< unsigned int > monitored_index;
    void monitor(Vector<Variable>& X);
    void monitor(Variable X);
#endif


    /// The set of constraints, with special accessors for triggers
    Vector< Constraint >         constraints;
    IntStack              posted_constraints;
    ConstraintQueue       active_constraints;
    //VariableQueue           active_variables;
    
    
    Vector< ConstraintTriggerArray > constraint_graph;


    /// For each level, the list of reversible objects that changed at this level, 
    // /// and will need to be restored
    // //Vector< Variable > saved_vars;
    // Vector< int > saved_vars;
    // /// For each level, the list of reversible objects that changed at this level, 
    // /// and will need to be restored
    // Vector< Constraint > saved_post;
    // Vector< Constraint > saved_relax;

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
    VarStack < Variable, ReversibleNum<int> >   sequence;
    //ReversibleNum<int> sequence_size;

    Vector< Decision > decisions;
    Vector< Clause* > reason;
    Vector< Lit > learnt_clause;
    BitSet visited;
    int num_search_variables;
    ConstraintClauseBase *base;


    /// Variable selection and branching
    BranchingHeuristic *heuristic;
    /// Restart policy
    RestartPolicy *policy; 
    /// Goal 
    Goal *objective;

    Vector<SolutionListener*>     solution_triggers;
    Vector<RestartListener*>       restart_triggers;
    Vector<SuccessListener*>       success_triggers;
    Vector<FailureListener*>       failure_triggers;
    Vector<DecisionListener*>     decision_triggers;
    Vector<VariableListener*>     variable_triggers;
    Vector<ConstraintListener*> constraint_triggers;

    /// the constraint responsible for the last fail (NULL otherwise) 
    Constraint culprit;
    /// the constraint-index of the last wiped_out variable in the culprit constraint (-1 otherwise)
    int wiped_idx;
    /// the overall-index of the variable responsible for the last trigger before a failure
    int wiper_idx;

    Vector< Expression* > expression_store;

    ConsolidateListener *consolidate_manager;

    Outcome satisfied();
    Outcome exhausted();

    /*!@name Constructor*/
    //@{
    Solver();

    void parse_dimacs(const char* filename);
    void set_parameters(SolverParameters& p);

    class BooleanMemoryManager {
    public:
      Vector< unsigned int > size;
      Vector<int*> slots;

      BooleanMemoryManager() {
	size.add(0);
	int *nslot = new int[1024];
	std::fill(nslot, nslot+1024, 3);
	slots.add(nslot);
      }
      virtual ~BooleanMemoryManager() {
	while(!slots.empty()) 
	  delete [] slots.pop();
      }

      void add(Variable *x);
      void add(Vector< Variable >& bool_vars);
    };

    unsigned int initialised_vars;
    unsigned int initialised_cons;
    BooleanMemoryManager booleans;
    //@}

    /*!@name Destructor*/
    //@{
    virtual ~Solver();
    bool is_initialised() {return variables.size <= initialised_vars;}
    //@}
  
    /*!@name Declarations*/
    //@{
    int declare(Variable x);
    /// add a variable (prior to search!!!)
    void add(Variable x);
    void remove(Variable x);
    void add(VarArray& x);
    void add(Constraint x); 
    void add(Vector< Lit >& clause); 
    //void add(ConstraintW x); 
    //void add(BranchingHeuristic* h);

    void add(SolutionListener* l);
    void add(RestartListener* l);
    void add(SuccessListener* l);
    void add(FailureListener* l);
    void add(DecisionListener* l);
    void add(VariableListener* l);
    void add(ConstraintListener* l);
    
    void remove(SolutionListener* l);
    void remove(RestartListener* l);
    void remove(SuccessListener* l);
    void remove(FailureListener* l);
    void remove(DecisionListener* l);
    void remove(VariableListener* l);
    void remove(ConstraintListener* l);

    void mark_non_convex(const int i) { 
      domain_types[i] &= (~RANGE_VAR); 
    }
    //void add(Vector<Literal>& clause);
    //@}

    /*!@name Trail accessors*/
    //@{
    // inline void save() {
    //   trail_.add(sequence.size);
    //   trail_.add(saved_objs.size);
    //   trail_.add(saved_cons.size);
    //   trail_.add(saved_vars.size);
    //   trail_.add(saved_post.size);
    //   trail_.add(saved_relax.size);

    //   ++statistics.num_nodes;
    //   ++level;
    // }
    /// called by reversible objects when they want to be restored when backtrack to this level
    // inline void save(Reversible* object) { saved_objs.add(object); }
    // inline void save(Constraint c) { saved_cons.add(c); }
    // //void save(Variable x); 
    // void save(const int idx); 
    // //void save(VariableImplementation *x, int dtype); 
    //@}

    /*!@name Propagation accessors*/
    //@{
    void notify_failure();
    void notify_success();
    void notify_decision();
    void notify_restart();
    void notify_relax(Constraint c);
    void notify_post(Constraint c);
    void notify_add_constraint(Constraint c);
    void notify_change_variable(const int idx);
    void notify_add_variable();
    /// called when the var'th variable of constraint cons changes (with event type evt)
    //void trigger_event(const int var, const Event evt);

    /// achieve propagation closure
    PropagationOutcome propagate(Constraint c, const bool t=true);
    PropagationOutcome checker_propagate(Constraint c, const bool t=true);
    PropagationOutcome bound_checker_propagate(Constraint c, const bool t=true);
    bool propagate(); 
    
    void fail();

    bool rewrite(); 
    void consolidate(); 
    void make_non_convex(const int idx);
    //void specialise(); 
    //void initialise_domains(); 
    //@}

    /*!@name Search accesors*/
    //@{
    /// make a branching decision
    void branch_left();
    /// undo the previous branching decision and enforce the complementary constraint
    void branch_right();
    void backjump();

    /// backtrack one level
    void restore();
    /// backtrack to a given level
    void restore(const int);
    

    //Solution *get_solution(Vector< Variable >& seq);
    void initialise_search(Vector< Variable >& seq, 
			   BranchingHeuristic *heu=NULL, 
			   RestartPolicy *pol=NULL,
			   Goal *goal=NULL);
    Outcome get_next_solution();

    /*!
      Launches a depth first search on the sequence 'seq' of variables
      with variable ordering 'heu' and restart policy 'pol'.
      Returns an outcome (SAT/UNSAT/OPT/UNKNOWN) and undo all decisions before returning.
    */
    Outcome depth_first_search(Vector< Variable >& seq, 
			       BranchingHeuristic *heu=NULL, 
			       RestartPolicy *pol=NULL,
			       Goal *goal=NULL);
    /*!
      Launches a depth first search on all variables
      with variable ordering 'heu' and restart policy 'pol'.
      Returns an outcome (SAT/UNSAT/OPT/UNKNOWN) and undo all decisions before returning.
    */
    Outcome depth_first_search(BranchingHeuristic *heu=NULL, 
			       RestartPolicy *pol=NULL,
			       Goal *goal=NULL);

    Outcome search();

    /*!
      Black box search.
    */
    Outcome solve();
    Outcome minimize(Variable X);
    Outcome maximize(Variable X);

    ///
    bool limits_expired();

    /// depth first search algorithm
    Outcome chronological_dfs();

    // /// sat search algorithm
    Outcome conflict_directed_backjump();
    void learn_nogood();
    void forget();
    // //@}


    BranchingHeuristic *heuristic_factory(std::string var_ordering, std::string branching);
    // {
    //   BranchingHeuristic *heu = NULL;
    //   if(var_ordering == "dom/wdeg") {
    // 	if(branching == "minval") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >, MinValue > (this); 
    // 	} else if(branching == "maxval") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >, MaxValue > (this); 
    // 	} else if(branching == "halfsplit") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >, HalfSplit > (this); 
    // 	} else if(branching == "randminmax") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >, RandMinMax > (this); 
    // 	} else if(branching == "minweight") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >, MinWeightValue > (this); 
    // 	} 
    //   } else if(var_ordering == "dom/activity") {
    // 	if(branching == "minval") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < PruningCountManager, MinDomainOverWeight >, MinValue > (this); 
    // 	} else if(branching == "maxval") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < PruningCountManager, MinDomainOverWeight >, MaxValue > (this); 
    // 	} else if(branching == "halfsplit") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < PruningCountManager, MinDomainOverWeight >, HalfSplit > (this); 
    // 	} else if(branching == "randminmax") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < PruningCountManager, MinDomainOverWeight >, RandMinMax > (this); 
    // 	} else if(branching == "minweight") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < PruningCountManager, MinDomainOverWeight >, MinWeightValue > (this); 
    // 	} 
    //   }						
    //   return heu;
    // }

    RestartPolicy *restart_factory(std::string rpolicy);
    // {
    //   RestartPolicy pol = new NoRestart();
    //   if(policy == "luby") pol = new Luby(); 
    //   else if(policy == "geom") pol = new Geometric(); 
    //   return pol;
    // }

    double *get_literal_activity();

    //void check_constraints();

    void check_constraint_graph_integrity();

    void extract_instance_statistics();

    void initialise_random_seed(const int seed);
    void set_time_limit(const double limit);

    /*!@name Printing*/
    //@{
    void debug_print();
    void full_print();
    void print_clist(int k) ;
    virtual std::ostream& display(std::ostream&, const int current=0)  ;
    //@}
  };

  std::ostream& operator<< (std::ostream& os, Solution& x);
  std::ostream& operator<< (std::ostream& os, Solution* x);

  std::ostream& operator<< (std::ostream& os, Solver& x);
  std::ostream& operator<< (std::ostream& os, Solver* x);

  std::ostream& operator<< (std::ostream& os, ConstraintQueue& x);
  std::ostream& operator<< (std::ostream& os, ConstraintQueue* x);

  std::ostream& operator<< (std::ostream& os, const SolverStatistics& x);
  std::ostream& operator<< (std::ostream& os, const SolverStatistics* x);

  std::ostream& operator<< (std::ostream& os, ConstraintTriggerArray& x);
  std::ostream& operator<< (std::ostream& os, ConstraintTriggerArray* x);



  SearchMonitor& operator<< (SearchMonitor& os, Variable& x);
  SearchMonitor& operator<< (SearchMonitor& os, VarArray& x);
  SearchMonitor& operator<< (SearchMonitor& os, Constraint& x);
  //SearchMonitor& operator<< (SearchMonitor& os, std::string& x);
  SearchMonitor& operator<< (SearchMonitor& os, const char* x);
  SearchMonitor& operator<< (SearchMonitor& os, const int x);


}

#endif // __SOLVER_HPP
