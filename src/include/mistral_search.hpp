
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



namespace Mistral {


  class RestartListener {
  public:
    int rid;
    virtual void notify_restart() = 0;
  };

  class DecisionListener {
  public:
    int did;
    virtual void notify_decision() = 0;
  };

  class SuccessListener {
  public:
    int sid;
    virtual void notify_success() = 0;
  };

  class FailureListener {
  public:
    int fid;
    virtual void notify_failure() = 0;
  };

  class ConstraintListener {
  public:
    int cid;
    virtual void notify_post(Constraint *c) = 0;
    virtual void notify_relax(Constraint *c) = 0;
  };

  class VariableListener {
  public:
    int vid;
    virtual void notify_add() = 0;
    virtual void notify_change(const int idx) = 0;
  };


  class LiteralActivityManager : public DecisionListener {

  public:

    Solver *solver;
    double *lit_activity; 
    double *var_activity;
    int n_vars;
 
    double decay;

    LiteralActivityManager(Solver *s) ;
    virtual ~LiteralActivityManager() ;

    virtual void notify_decision() ;    
    double *get_weight() ;
    
  };


  class FailureCountManager : public FailureListener, public ConstraintListener {

  public:

    Solver *solver;

    Vector<double> constraint_weight;
    Vector<double> variable_weight;

    FailureCountManager(Solver *s) : solver(s) {
      for(unsigned int i=0; i<solver->variables.size; ++i) {
	variable_weight.add(solver->variables[i].get_degree());
      }
      for(unsigned int i=0; i<solver->constraints.size; ++i) {
	constraint_weight.add(1);
      }

      solver->add((FailureListener*)this);
      solver->add((ConstraintListener*)this);
    }

    virtual ~FailureCountManager() {}

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
  };


  class PruningCountManager : public SuccessListener {

  public:

    Solver *solver;

    Vector<double> variable_weight;

    PruningCountManager(Solver *s) : solver(s) {
      for(unsigned int i=0; i<solver->variables.size; ++i) {
	variable_weight.add(solver->variables[i].get_degree());
      }
      solver->add((SuccessListener*)this);
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





  class RestartPolicy {
    
  public:

    unsigned int base;
    
    RestartPolicy(const unsigned int b=256);
    virtual void reset(unsigned int& limit) = 0;
    
  };


  class NoRestart : public RestartPolicy {
    
  public:
    
    NoRestart();
    virtual ~NoRestart();
    
    void reset(unsigned int& limit) {
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
    
  };


  /**********************************************
   * Variable Selection
   **********************************************/

  class BranchingHeuristic {

  public:
    
    Solver *solver;

    BranchingHeuristic(Solver *s) {solver = s;}
    virtual ~BranchingHeuristic() {}

    virtual void initialise(VarStack< Variable >& seq) {}
    
    virtual Decision branch() = 0;

  };


  // class VarOrdering {

  // public: 
  //   Solver *solver;

  //   VarOrdering() {}
  //   VarOrdering(Solver *s);
  //   void initialise(Solver *s) { solver = s; }
  //   virtual ~VarOrdering();
    
  // };

  class NoOrder {

  public: 
    Solver *solver;

    NoOrder() { solver = NULL; }
    NoOrder(Solver *s);
    virtual ~NoOrder();
    void initialise(Solver *s) { solver = s; }
    void initialise(VarStack< Variable >& seq) {}
    
    Variable select();

  };

  class Lexicographic : public VariableListener {

  public: 
    
    Solver *solver;
    Vector< Variable >           order;
    Vector< int >                index;
    ReversibleNum< unsigned int > last;

    Lexicographic() { solver = NULL; }
    Lexicographic(Solver *s);
    void initialise(Solver *s);
    void initialise(VarStack< Variable >& seq);
    virtual ~Lexicographic();

    virtual void notify_add() {};
    virtual void notify_change(const int idx);
    
    Variable select();

  };


  /**********************************************
   * Variable Selection
   **********************************************/

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


  class MinDomainOverWeight
  {
  public: 

    /**@name Constructors*/
    //@{
    MinDomainOverWeight() {dom_ = LARGE_VALUE; wei_ = 0;}
    //MinDomainOverWeight(int *w) : weight(w) {dom_ = LARGE_VALUE; wei_ = 0;}
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
   * Value Selection
   **********************************************/

  class MinValue {

  public: 
    
    MinValue() {}
    MinValue(Solver *s) {}
    virtual ~MinValue() {};
    
    inline Decision make(Variable x) {
      Decision d(x, Decision::ASSIGNMENT, x.get_min());
      return d;
    }

  };

  class MaxValue {

  public: 
    
    MaxValue() {}
    MaxValue(Solver *s) {}
    virtual ~MaxValue() {};
    
    inline Decision make(Variable x) {
      Decision d(x, Decision::ASSIGNMENT, x.get_max());
      return d;
    }

  };

  class HalfSplit {

  public: 
    
    HalfSplit() {}
    HalfSplit(Solver *s) {}
    virtual ~HalfSplit() {};
    
    inline Decision make(Variable x) {
      Decision d(x, Decision::UPPERBOUND, (x.get_min()+x.get_max())/2);
      return d;
    }

  };

  class RandomMinMax {

  public: 
    
    RandomMinMax() {}
    RandomMinMax(Solver *s) {}
    virtual ~RandomMinMax() {};
    
    inline Decision make(Variable x) {
      Decision d(x, Decision::ASSIGNMENT, 
		 (randint(2) ? x.get_min() : x.get_max()));
      return d;
    }

  };


  class MinWeightValue {

  public: 

    double **weight;
    
    MinWeightValue() {}
    MinWeightValue(Solver *s) {}
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


  class BoolMinWeightValue {

  public: 

    double **weight;
    
    BoolMinWeightValue() {}
    BoolMinWeightValue(Solver *s) {}
    virtual ~BoolMinWeightValue() {};
    
    inline Decision make(Variable x) {
      double *wgt = weight[x.id()];
      Decision d(x, Decision::ASSIGNMENT, (wgt[0]<wgt[1] ? 0 : 1));
      return d;
    }

  };


  /**********************************************
   * Generic heuristic
   **********************************************/
  
  template < class VarSelector >
  class GenericDVO 
  {
  public: 
    
    /**@name Parameters*/
    //@{ 
    Solver   *solver;
    VarSelector best;
    VarSelector current;

    //ValSelector choice;
    //@}

    /**@name Constructors*/
    //@{
    GenericDVO() { solver = NULL; }
    GenericDVO(Solver* s) : solver(s) {}
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
      //std::cout << var << ":" ;
      //best.display(std::cout);
      //for(unsigned int i=1; i<length; ++i) {
      for(unsigned int i=length; i--;) {
	current = variables[i];
	if( current < best ) {
	  best = current;
	  var = variables[i];
	  //std::cout << " < " << var << ":" ;
	  //best.display(std::cout);
	}
      }
      //std::cout << std::endl;
      return var;
    }
    //@}
  };


  /**********************************************
   * GenericRandom heuristic
   **********************************************/
  
  template < class VarSelector, int RAND >
  class GenericRandomDVO 
  {
  public: 

    /**@name Parameters*/
    //@{ 
    Solver      *solver;
    VarSelector   bests[RAND+1];
    VarSelector current;
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



  template < class WeightManager, class VarSelector >
  class GenericWeightedDVO : public GenericDVO< VarSelector >
  {
  public: 
    
    /**@name Parameters*/
    //@{ 
    WeightManager *manager;
    //@}

    /**@name Constructors*/
    //@{
    GenericWeightedDVO() : GenericDVO< VarSelector >() { manager=NULL; }
    GenericWeightedDVO(Solver* s) : GenericDVO< VarSelector >(s) {
      initialise_manager();
    }
    virtual ~GenericWeightedDVO() { delete manager; }
    virtual void initialise(Solver *s) { 
      GenericDVO< VarSelector >::initialise(s);
      initialise_manager();
    }
    virtual void initialise(VarStack< Variable >& seq) {}
    virtual void initialise_manager() {
      manager = new WeightManager(GenericDVO< VarSelector >::solver);
      //GenericDVO< VarSelector >::solver->add(manager);
      GenericDVO< VarSelector >::best.weight = manager->get_weight();
      GenericDVO< VarSelector >::current.weight = manager->get_weight();
    }
    //@}
    
  };


  /**********************************************
   * Generic Random Weighted heuristic
   **********************************************/
  
  template < class WeightManager, class VarSelector, int RAND >
  class GenericWeightedRandomDVO : public GenericRandomDVO< VarSelector, RAND >
  {
  public: 
    
    /**@name Parameters*/
    //@{ 
    WeightManager *manager;
    //@}

    /**@name Constructors*/
    //@{
    GenericWeightedRandomDVO() : GenericRandomDVO< VarSelector, RAND >() { manager=NULL; }
    GenericWeightedRandomDVO(Solver* s) : GenericRandomDVO< VarSelector, RAND >(s) {
      initialise_manager();
    }
    virtual ~GenericWeightedRandomDVO() { delete manager; }
    virtual void initialise(Solver *s) { 
      GenericRandomDVO< VarSelector, RAND >::initialise(s);
      initialise_manager();
    }
    virtual void initialise(VarStack< Variable >& seq) {}
    virtual void initialise_manager() {
      manager = new WeightManager(GenericRandomDVO< VarSelector, RAND >::solver);
      //GenericRandomDVO< VarSelector >::solver->add(manager);
      GenericRandomDVO< VarSelector, RAND >::current.weight = manager->get_weight();
      for(int i=0; i<=RAND; ++i)
	GenericRandomDVO< VarSelector, RAND >::bests[i].weight = manager->get_weight();
    }
    //@}
    
  };


  template < class VarSelector, class ValSelector >
  class GenericHeuristic : public BranchingHeuristic {
  public:

    VarSelector var;
    ValSelector choice;

    GenericHeuristic(Solver *s) 
      : BranchingHeuristic(s) 
    {
      var.initialise(s);
      choice = ValSelector(s);
    }

    virtual void initialise(VarStack< Variable >& seq) {var.initialise(seq);}

    virtual Decision branch() {
      return choice.make(var.select());
    }

  };


   // typedef GenericHeuristic< 
   //  GenericWeightedDVO< LiteralActivityManager, MaxWeight >, 
   //  MinValue >                                                          VSIDS;

  class VSIDS : public GenericHeuristic< GenericWeightedDVO< LiteralActivityManager, MaxWeight >, 
					 BoolMinWeightValue > {

  public:
    
    VSIDS(Solver *s) : GenericHeuristic< GenericWeightedDVO< LiteralActivityManager, MaxWeight >, 
					 BoolMinWeightValue >(s) {
      int j, i, n=var.manager->n_vars;
      choice.weight = new double*[n];
      for(i=0; i<n; ++i) {
	choice.weight[i] = new double[2];
      }
      // choice.weight[0] = new double[n];
      // choice.weight[1] = new double[n];

      
	for(i=0; i<n; ++i) 
	  for(j=0; j<2; ++j)
	    choice.weight[i][j] = var.manager->lit_activity[2*i+j];
    }

    virtual ~VSIDS() {
      int  i, n=var.manager->n_vars;
      for(i=0; i<n; ++i) 
	delete [] choice.weight[i];
  
      // delete [] choice.weight[0];
      // delete [] choice.weight[1];
      delete [] choice.weight;
    }
    
  };

//   template < class VarSelector, class ValSelector >
//   class GenericPtrHeuristic : public BranchingHeuristic {
//   public:

//     VarSelector *var;
//     ValSelector choice;

//     GenericPtrHeuristic(Solver *s, VarSelector *var_ord) 
//       : BranchingHeuristic(s) 
//     {
//       var = var_ord; 
//       choice = ValSelector(s);
//     }

//     virtual Decision branch() {
//       return choice.make(var->select());
//     }

//   };


//   class PCP : public BranchingHeuristic {

//   public:
    
//     Solver *solver;

//     int nJobs;
//     int nMachines;
    
//     Vector< Variable > *conflict;

//     Constraint ***precedences;


//     BranchingHeuristic(Solver *s, VarArray& tasks, int nj, int nm) {
//       solver = s;
//       nJobs = nj;
//       nMachines = nm;

//       conflict = new Vector< Variable >[nJobs*nMachines];
//       for(int i=0; i<nJobs; ++i) {
// 	for(int j=0; j<nMachines; ++j) {
// 	  // record all variables that may be in conflict with 
// 	}
//       }

//       precedences = new Constraint**[nJobs];
//       for(int i=0; i<nJobs; ++i)
// 	precedences[i] = new Constraint*[nMachines];


//     }
//     virtual ~BranchingHeuristic() {}
    
//     virtual Decision branch() = 0;

//   };
  
//   class Search {

//   public:
//     Stack< IntVar > sequence;

//     Vector< unsigned int > trail_;
//     Vector< IntVar > decision;

//     VarOrdering *heuristic;   
//     RestartPolicy *policy;   
    
//     SolverParameters parameters;
//     SolverStatistics statistics;

//     Search() {
//       heuristic = NULL;
//       policy = NULL;
//     }

//     void initialise(Solver *s);
//     IntVar backtrack(); 
    
//     IntVar make_node() {
//       ++statistics.num_nodes;
//       trail_.add(sequence.size);
//       decision.add(heuristic->select());
//       return decision.back();
//     }

//     void assign(IntVar x);
//     Outcome satisfied();
//     bool limitsExpired();

//     void init_search(Vector< IntVar >& seq, VarOrdering *h, RestartPolicy *p);
//   };


//   inline Outcome Search::satisfied() {    
// #ifdef _DEBUG_SEARCH
//     std::cout << "c";
//     for(unsigned int k=0; k<=decision.size; ++k) std::cout << " ";
//     std::cout << " SAT!" << std::endl; 
// #endif
    
//     ++statistics.num_solutions;
    
//     return SAT;
//   }
  
// //   inline bool Search::limitsExpired() {

// //     return (parameters.limit && 
// // 	    ((parameters.time_limit > 0.0 && (getRunTime() - statistics.start_time) > parameters.time_limit) ||
// // 	     (parameters.node_limit > 0 && (statistics.num_nodes > parameters.node_limit)) ||
// // 	     (parameters.fail_limit > 0 && (statistics.num_failures > parameters.fail_limit)) ||
// // 	     (parameters.restart_limit > 0 && (statistics.num_restarts > parameters.restart_limit)) ||
// // 	     (parameters.backtrack_limit > 0 && (statistics.num_backtracks > parameters.backtrack_limit))
// // 	     ));
// //   }

}

#endif // __SEARCH_HPP
