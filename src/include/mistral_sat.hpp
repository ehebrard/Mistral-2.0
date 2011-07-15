
/*
  Mistral is a constraint satisfaction and optimisation library
  Copyright (C) 2003-2005  Emmanuel Hebrard
  
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

  The author can be contacted electronically at 
  ehebrard@cse.unsw.edu.au.
*/

/** \file sat.h
    \brief Header for the SAT Solver.
*/


#ifndef _SAT_H
#define _SAT_H

#include <iostream>
#include <iomanip>
#include <assert.h>

#include <mistral_search.hpp>

// -n..-1,1..n -> 0..2n
// dimacs to literal x>0 -> 2(x-1)+1
// dimacs to literal x<0 -> -2(x+1)

// 0..2n -> 0..n
// literal to atom l -> l/2


namespace Mistral {

  // typedef unsigned int Lit;
  // typedef unsigned int Atom;
  // typedef unsigned int Value;

  // typedef Array<Lit> Clause;


  std::ostream& operator<< (std::ostream& os, const Clause& x);
  std::ostream& operator<< (std::ostream& os, const Clause* x);
  

#define SIGN(l) ((l)&1)
#define NEG(l) ((l)^1)
#define UNSIGNED(l) ((l)/2)
#define ATOM(l) (state[(l)/2])
#define LEVEL(a) ((a)/2)

#define V_TRUE    1
#define V_FALSE   0
#define V_UNKNOWN 2

  void print_clause(std::ostream& o, Clause* cl) ;
  void print_literal(std::ostream& o, Lit l, bool dir=true) ;

  class RestartPolicy;

  class SatSolver 
  {

  public:
    /**@name Parameters*/
    //@{
    /// UNKNOWN - SAT - UNSAT - LIMITOUT
    int status;

    // (first bit stands for the truth value, 
    // the rest stands for the index in the vector 'assumptions')
    Vector< unsigned int > state;

    // The list of literals that have been decided/inferred
    Vector< Atom > assumptions;


    inline int truth(Atom a) {
      return SIGN(state[a]);
    }

    inline int assigned(Atom a) {
      return LEVEL(state[a])<=assumptions.size;
    }
    
    inline bool satisfied(Lit p) {
      unsigned int st = state[UNSIGNED(p)];
      return LEVEL(st) <= assumptions.size && SIGN(st) == SIGN(p);
    }
    
    inline void make_assumption(Lit p) {
      
      // swap p and assumptions[assumptions.size]
      Atom q_atom = assumptions.back(0);
      Atom p_atom = UNSIGNED(p);
      unsigned int p_level = LEVEL(state[p_atom]);
      state[p_atom] = assumptions.size*2 + SIGN(p);
      state[q_atom] = p_level*2 + SIGN(state[q_atom]);

      assumptions[p_level] = q_atom;
      assumptions.add(p_atom);
    }

    // size of the vector 'assumptions' at each level
    Vector< unsigned int > decisions;

    ///// List of decision "levels", used to backtrack
    //Vector< int > decision_level;

    /// the clause that entailed this atom
    Vector< Clause* > reason;

    /// pointers to the scope of the (original) clauses,
    Vector< Clause* > original;

    /// pointers to the scope of the current state of the (original) clauses,
    Vector< Clause* > base;

    /// pointers to the scope of the (learnt) clauses,
    Vector< Clause* > learnt;

    /// for each Lit, the list of clauses it is watched by
    Vector< Vector< Clause*> > is_watched_by;

    /// Lit Activity 
    Vector< double > activity;

  
    /// utils
    BitSet visited;
    Vector< Lit > learnt_clause;
    int next_deduction;

    /// search statistics
    SolverStatistics stats;
    /// search parameters
    SolverParameters params;
    /// Restart
    RestartPolicy *restart_policy;
    //@}


    /**@name Constructors*/
    //@{
    SatSolver();
    SatSolver(const char* filename);
    //SatSolver(CSP& model);
    virtual ~SatSolver();
    void set_parameters(SolverParameters& p);
    void set_policy(const int p);
    virtual void init_vars(const int n, const int m);
    virtual void init_watchers();
    void parse_dimacs(const char* filename);
    //@}

    /**@name Miscellanous*/
    //@{
    std::ostream& display(std::ostream& o) const;
    void print_all(std::ostream& o) const;
    void print_watchers(std::ostream& o, int beg=NOVAL, int end=NOVAL) const;
    void print_decisions(std::ostream& o, bool m=true) const;
    void print_clauses(std::ostream& o) const;
    //@}

    /**@name Solving Methods*/
    //@{
    /// Solves the problem and set the status accordingly
    virtual int solve();
    /// Minisat style conflict directed search
    int iterative_search();
    /// Finds an explanation and uses it to backjump
    int analyze( Clause *conflict, Lit& lit, const bool learn );
    /// Chooses the next Lit to branch on
    Lit choice();
    /// Create a choice point and add l to the clause base
    void make_decision(const Lit l);
    /// Backjump to the choice point where l was entailed
    void backtrack_to( const int backtrackLevel );
    /// Add l to the clause base
    //void add_lit(const Atom a, const Lit l);
    void add_lit(const Lit l);
    /// Returns the conflict index, or -1 if there was no conflict
    Clause* unit_propagate();
    /// Returns the conflict index, or -1 if there was no conflict
    Clause* update_watcher(const int cw, const Lit p);
    //@}

    /**@name Clause Base Methods*/
    //@{
    /// normalize the activity values in the interval [0,1]
    void normalize_activity(const double M);
    /// Reduces every literal's activity by a constant factor
    void decay_activity();
    /// Returns NULL if the clause is always true, and a reduced version otherwise
    Clause* reduce(Clause* clause);
    /// Reduces all clauses that can be
    void simplify_data_base();
    /// Add a clause to the base/learnt
    void add_clause( Vector<Lit>& conflict );
    void add_clause( Vector<Clause*>& clauseList, 
		    Vector<Lit>& conflict,
		    double& avgsize );
    /// Remove a clause from the base/learnt
    void remove_clause( Vector<Clause*>& clauseList, 
		       const int cidx,
		       double& avgsize );
    /// Forget learnt clauses that do not meet a given criterion
    void forget();
    /// Add a clause to the original base
    void add_original_clause( Vector<Lit>& conflict );
    //@}

    /**@name Utils*/
    //@{
    int check_solution();
    //int atom( const Lit l ) const;
    void shuffle();
    bool rlimit_expired();
    bool limit_expired();
    //@}
  };

  std::ostream& operator<< (std::ostream& os, SatSolver& x);
  std::ostream& operator<< (std::ostream& os, SatSolver* x);


  /***********************************************
   * NogoodBase Constraint (forward checking).
   ***********************************************/
  /*! \class ConstraintNogoodBase
    \brief   Constraint.
  */
  class ConstraintClauseBase : public Constraint {

  public:

    /**@name Parameters*/
    //@{ 
    // if there was a conflict it is stored there:
    Clause* conflict;
    //Vector< Clause* > reason;
    Clause** reason;
    Vector< double > lit_activity;
    Vector< double > var_activity;
    // list of clauses
    Vector< Clause* > clauses;
    Vector< Clause* > learnt;
    // the watched literals data structure
    Vector< Vector< Clause* > > is_watched_by;
    //@}
    
    /**@name Constructors*/
    //@{
    ConstraintClauseBase() : Constraint() { conflict = NULL; }
    ConstraintClauseBase(Vector< Variable >& scp);
    virtual void mark_domain();
    virtual Constraint *clone() { return new ConstraintClauseBase(scope); }
    virtual void initialise();
    virtual ~ConstraintClauseBase();
    
    void add( Variable x );
    void add( Vector < Lit >& clause, double init_activity=0.0 );
    void learn( Vector < Lit >& clause );
    void remove( const int cidx );
    void forget( double forgetfulness );
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    Clause* update_watcher(const int cw, const Lit p, PropagationOutcome& o);
    //virtual PropagationOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "clause_base"; }
    //@}
    
  };

};

using namespace Mistral;




int compar(const void *a, const void *b);
void initSort(double *sa);


inline void SatSolver::shuffle()
{
  unsigned int i, j, k;
  for(i=assumptions.size; i<state.size; ++i)
    {
      j = i+randint(state.size-i);
      k = assumptions[j];
      assumptions[j] = assumptions[i];
      assumptions[i] = k;
      state[k] = 2*i+SIGN(state[k]);
    }
}


inline Clause* SatSolver::reduce(Clause* clause)
{
  Atom x;
  Lit* data = clause->data;
  int j, sz;

  // check the watchs first
  for(j=2; j;) {
    x = ATOM(data[--j]);
    if( LEVEL(x)<assumptions.size && SIGN(x) == SIGN(data[j]) ) {
      return NULL;
    }
  }
  sz = clause->size;
  for(j=sz; j>2;) {
    x = ATOM(data[--j]);
    if( LEVEL(x)>=assumptions.size ) continue;
    if( SIGN(x) == SIGN(data[j]) ) return NULL;
    data[j] = data[--sz];
  }

  clause->size = sz;
  return clause;
}

inline void SatSolver::simplify_data_base()
{
  unsigned int i;
  double oldsize;
  for(i=base.size; i;) {
    oldsize = base[i-1]->size;
    if(!reduce(base[--i])) {
      remove_clause(base, i, stats.base_avg_size);
    } else if(oldsize != base[i]->size) {
      oldsize -= base[i]->size;
      stats.base_avg_size -= oldsize/(double)(base.size);
    }
  }
  for(i=learnt.size; i;) {
    oldsize = learnt[i-1]->size;
    if(!reduce(learnt[--i])) {
      remove_clause(learnt, i, stats.learnt_avg_size );
    } else if(oldsize != learnt[i]->size) {
      oldsize -= learnt[i]->size;
      stats.learnt_avg_size -= oldsize/(double)(learnt.size);
    }
  }
}

inline void SatSolver::add_original_clause( Vector<Lit>& conflict )
{
  Clause *cl = Clause::Array_new(conflict);
  original.add( cl );
}

inline void SatSolver::add_clause( Vector<Lit>& conflict )
{
  if(conflict.size > 1) {
    Clause *cl = Clause::Array_new(conflict);
    base.add( cl );
    double size = base.size;
    stats.base_avg_size = (stats.base_avg_size*(size-1) + double(conflict.size))/size;
    if( conflict.size < 4 ) ++stats.small;
  } else {
    Lit p = conflict[0];
    Atom x = ATOM(p);
    if(LEVEL(x) > decisions.size)
      add_lit(p);
    else if(SIGN(p) != SIGN(x)) 
      status = UNSAT;
  }
}

inline void SatSolver::add_clause( Vector<Clause*>& clauseList, 
				  Vector<Lit>& conflict,
				  double& avgsize )
{
  if(conflict.size > 1) {
    Clause *cl = Clause::Array_new(conflict);
    clauseList.add( cl );
    is_watched_by[conflict[0]].add(cl);
    is_watched_by[conflict[1]].add(cl);
    double size = clauseList.size;
    avgsize = (avgsize*(size-1) + double(conflict.size))/size;
    if( conflict.size < 4 ) ++stats.small;

  } else {
    Lit p = conflict[0];
    Atom x = ATOM(p);
    if(LEVEL(x) > decisions.size)
      add_lit(p);
    else if(SIGN(p) != SIGN(x)) 
      status = UNSAT;
  }
}

inline void SatSolver::remove_clause( Vector<Clause*>& clauseList, 
				     const int cidx,
				     double& avgsize )
{
  Clause *clause = clauseList[cidx];

  is_watched_by[clause->data[0]].remove_elt( clause );
  is_watched_by[clause->data[1]].remove_elt( clause );
  clauseList.remove( cidx );

  double size = clauseList.size;
  if(size > 0) avgsize = (avgsize*(size+1) - double(clause->size))/size;
  else         avgsize = 0;

  free(clause);
}

inline void SatSolver::forget()
{

  if( params.forgetfulness > 0.0 ) {
    int nlearnt = learnt.size;
    double sa[nlearnt];
    Clause *tmp[nlearnt];
    int i, j, order[nlearnt];
    initSort(&(sa[0]));
    for(i=0; i<nlearnt; ++i)
      {
	order[i] = i;
	sa[i] = 0.0;
	Clause& clause = *(learnt[i]);
	j=clause.size;
	while(j--)
	  sa[i] += activity[clause[j]];
	sa[i] /= (clause.size * clause.size);
      }
    qsort(order, nlearnt, sizeof(int), compar);
    for(i=0; i<nlearnt; ++i)
      tmp[i] = learnt[order[i]]; 
    for(i=0; i<nlearnt; ++i) {
      learnt[i] = tmp[i];
    }
    
    int keep = (int)((double)nlearnt * (1.0-params.forgetfulness));
    
    for(i=nlearnt; i>keep;)
      remove_clause( learnt, --i, stats.learnt_avg_size );
    while(i>1) {
      --i;
      if(sa[order[i]] == sa[order[i-1]])
	remove_clause( learnt, i, stats.learnt_avg_size );
    }
  }
}

inline bool SatSolver::rlimit_expired()
{
  return( params.restart_limit && (stats.num_failures > params.restart_limit) );
}

inline bool SatSolver::limit_expired()
{
  return( params.time_limit>0 && ((getRunTime() - stats.start_time) > params.time_limit) );
}

inline int SatSolver::iterative_search()
{

  Lit p = 0;
  Atom a;
  while(status == UNKNOWN) {
    
    Clause *conflict = unit_propagate();
    
// #ifdef _DEBUG_SEARCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
//   print_decisions(std::cout);
// #endif

  if(decisions.size == 0 && assumptions.size) simplify_data_base();
  
  if( conflict ) {

#ifdef _DEBUG_SEARCH
      for(unsigned int d=0; d<decisions.size; ++d) std::cout << " ";
      std::cout << "conflict: ";
      print_clause(std::cout, conflict);
#endif

      if( decisions.size == 0 ) status = UNSAT;
      else if( limit_expired() ) {

#ifdef _DEBUG_SEARCH
	std::cout << " limit reached" << std::endl;
#endif

	return UNKNOWN;
      } else {

	int bt_level = 	analyze( conflict, p, true );
	
	if( rlimit_expired() ) {
	  
// #ifdef _DEBUG_SEARCH
// 	  std::cout << " restart limit reached" << std::endl;
// #endif
	  
	  backtrack_to(0);
	  status = LIMITOUT;
	} 
	else backtrack_to( bt_level ); //analyze( conflict, p, true ) );
      }
      if(status != LIMITOUT) {
	a = ATOM(p);

	if( !decisions.size && LEVEL(a)<assumptions.size )
	  status = UNSAT;
	else {	
	  next_deduction = assumptions.size;

#ifdef _DEBUG_SEARCH
	  std::cout << "c " ;
	  for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
	  //std::cout << "deduce: " ;
	  print_literal(std::cout, p, false);
	  std::cout << std::endl;
#endif

	  add_lit(p);
	}
      }
    } else {
      if( assumptions.size < state.size ) make_decision(choice());
      else status = check_solution();
    }
  }

  return status;
}

inline int SatSolver::analyze( Clause *conflict, Lit& lit, const bool learn )
{
#ifdef _DEBUG_NOGOOD
#ifdef _DEBUG_SEARCH
  std::cout << std::endl;
#endif
  print_decisions( std::cout , 0 );
#endif

#ifdef _CHRONOLOGICAL

  Atom a = assumptions[decisions.back()];
  lit = (2*a) | (SIGN(state[a]) ^ 1);
  unsigned int backtrackLevel = decisions.size-1;

#else
  
  learnt_clause.clear();

  unsigned int j, backtrackLevel = 0;
  int pathC = 0, index = assumptions.size;
  Lit p=0, q;
  Atom a;
  unsigned int lvl;
  
  learnt_clause.add(p);
  
  do {
    
    // add the parents of the conflict to the current set of visited atoms
    Clause& con = *conflict;
    
#ifdef _DEBUG_NOGOOD
    print_clause( std::cout, conflict );
    std::cout << std::endl;
#endif

    for(j=0; j<con.size; ++j) {
      q = con[j];
      a = UNSIGNED(q);
      lvl = LEVEL(state[a]);

#ifdef _DEBUG_NOGOOD
      std::cout << "\t" ;
      print_literal(std::cout, q); 
      std::cout << ": ";
#endif

      if( !visited.fast_contain(a) ) {

	activity[q] += params.activity_increment;
	visited.fast_add(a);
	// we'll need to replace 'a' by its parents since its level is too high
	if(lvl >= decisions.back()) {

#ifdef _DEBUG_NOGOOD
	  std::cout << "expend" << std::endl;
#endif

	  ++pathC;
	} else {
	  // q's level is below the current level, hence we are not expending it further

	  learnt_clause.add(q);

#ifdef _DEBUG_NOGOOD
	  std::cout << "add to the clause" ;
	  for(unsigned int k=0; k<learnt_clause.size; ++k) {
	    std::cout << " ";
	    print_literal(std::cout, learnt_clause[k]);
	  }
	  std::cout << std::endl;
#endif

	  if(lvl > backtrackLevel)
	    backtrackLevel = lvl;
	}
      }

#ifdef _DEBUG_NOGOOD
      else {
	std::cout << "visited" << std::endl;
      }
#endif

    }
    // jump to the next visited atom that need be further expended
    while(!visited.fast_contain(assumptions[--index]));
    a = assumptions[index];
    p = (2*a) | SIGN(state[a]); //polarity[a];
    lvl = LEVEL(state[a]);

#ifdef _DEBUG_NOGOOD
    std::cout << "explore ";
    print_literal(std::cout, p); 
    std::cout << " ";
    //std::cout.flush();
#endif

    if( pathC > 1 ) {
      // there are still atoms to expand, we start with 'a'
      conflict = reason[a];
      visited.fast_add(a);
    } 
#ifdef _DEBUG_NOGOOD
    else {
      std::cout << std::endl;
    }
#endif

  } while( --pathC );
  // p is the last decision, since all atoms above it in the
  // assumption stack have been skipped or expended.
  learnt_clause[0] = NEG(p);    

#ifdef _DEBUG_SEARCH
  std::cout << " (";
  for(unsigned int i=0; i<learnt_clause.size; ++i) {
    std::cout << " " ;//<< learnt_clause[i];
    print_literal(std::cout, learnt_clause[i]);
  }
  std::cout << " )" << std::endl;
#endif
        
  if( learn && learnt_clause.size != 1 ) {
    add_clause( learnt, learnt_clause, stats.learnt_avg_size );
    reason[UNSIGNED(p)] = learnt.back();
  }
  visited.clear();
  lit = NEG(p); 

  // get the backtrack level from the index
  lvl = decisions.size-1;
  while(lvl && decisions[lvl-1]>backtrackLevel) {
    --lvl;
  }
  backtrackLevel = lvl;

#ifdef _DEBUG_NOGOOD
  std::cout << "backtrackLevel = " << backtrackLevel << "/" << (decisions.size) << std::endl;
#endif

#endif

  return backtrackLevel;
}

inline void SatSolver::normalize_activity( const double M ) 
{
  double max_act = activity[0];
  double min_act = activity[0];
  for(unsigned int i=1; i<2*state.size; ++i) {
    if(activity[i] > max_act) max_act = activity[i];
    if(activity[i] < min_act) min_act = activity[i];
  }

  max_act -= min_act;

  if(max_act == 0) max_act = 1;

  for(unsigned int i=0; i<2*state.size; ++i) {
    //    std::cout << activity[i] << " -> ";
    activity[i] -= min_act;
    activity[i] /= max_act;
    activity[i] *= M;
    //    std::cout << activity[i] << std::endl;
  }
  
  //  exit(1);
}

inline void SatSolver::decay_activity()
{
  if(params.activity_decay >= 0) {
    unsigned int i = 2*state.size;
    while(i--) {
      activity[i] *= params.activity_decay;
    }
  }
}

#define FALSE_FIRST 0
#define TRUE_FIRST 1
#define LEAST_ACTIVE 2
#define MOST_ACTIVE 3
#define RANDOM 4
#define PERSISTENT 5

inline Lit SatSolver::choice()
{
  decay_activity();

  Lit p;

  if(params.randomization) {
    unsigned int crd = 0;
    double best[params.randomization], cur;
    unsigned int i=assumptions.size, j, k;
    Atom x[params.randomization], y;
    
    while(i < state.size)
      {	       
	y = 2*assumptions[i];
	
	cur = activity[y];
	cur += activity[NEG(y)]; 
	
	//std::cout << "activity of b" << UNSIGNED(y) << " = " << cur << std::endl;
	
	for(j=crd; j && cur>best[j-1]; --j);
	for(k=crd; k>j; --k) {
	  x[k] = x[k-1];
	  best[k] = best[k-1];
	}
	best[j] = cur;
	x[j] = y;
	
	if(crd<params.randomization) ++crd;
	++i;
      }
    
    j = randint(crd);
    
    //case FALSE_FIRST : 
    p = x[j]; //break;
  } else {
    p = 2*assumptions.back(0);
  }


  // switch( params.value_selection ) 
  //   {
  //   case TRUE_FIRST  : p = NEG(x[j]); break;
  //   case LEAST_ACTIVE: p = (activity[NEG(x[j])] > activity[x[j]] ? NEG(x[j]) : x[j]); break;
  //   case MOST_ACTIVE : p = (activity[NEG(x[j])] > activity[x[j]] ? x[j] : NEG(x[j])); break;
  //   case RANDOM      : p = (x[j] | (randint(2))); break;
  //   case PERSISTENT  : p = (x[j] | SIGN(state[UNSIGNED(x[j])])); break;
  //   }


  switch( params.value_selection ) 
    {
    case TRUE_FIRST  : p = NEG(p); break;
    case LEAST_ACTIVE: p = (activity[NEG(p)] > activity[p] ? NEG(p) : p); break;
    case MOST_ACTIVE : p = (activity[NEG(p)] > activity[p] ? p : NEG(p)); break;
    case RANDOM      : p |= (randint(2)); break;
    case PERSISTENT  : p |= SIGN(state[UNSIGNED(p)]); break;
    default: ; //{std::cout << "default choice" << std::endl;} ;
    }


  return p;    
}

inline void SatSolver::make_decision(const Lit l)
{

#ifdef _DEBUG_SEARCH
  std::cout << "c  ";
  for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
  //std::cout << "decide: " ;
  print_literal(std::cout, l);
  std::cout << std::endl;
#endif

  ++stats.num_nodes;
  next_deduction = assumptions.size;
  decisions.add( next_deduction );
  Atom a = (UNSIGNED(l));
  reason[a] = NULL;
  add_lit(l);

}

inline void SatSolver::backtrack_to( const int backtrackLevel )
{

  assumptions.pop_until( decisions.pop_until(backtrackLevel) );

// #ifdef _DEBUG_SEARCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
//   std::cout << "backtrack to level " << backtrackLevel << std::endl;
// #endif
  
}

inline void SatSolver::add_lit(const Lit l)
{    
  if(!decisions.size) ++stats.literals;  
  make_assumption(l);
}

inline Clause* SatSolver::unit_propagate()
{
  Clause* conflict = NULL;
  int cw;
  Lit p;
  Atom a;
  while( next_deduction < (int)(assumptions.size) ) {
    ++stats.num_propagations;
    a = assumptions[next_deduction];
    p = (2*a)|NEG((SIGN(state[a])));

#ifdef _DEBUG_UNITPROP
    for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
    std::cout << "propagate " ;
    print_literal(std::cout, NEG(p));
    std::cout << " " << is_watched_by[p].size << std::endl;
#endif

    cw = is_watched_by[p].size;
    while(cw-- && !conflict) {
      conflict = update_watcher(cw, p);
    }
    ++next_deduction;

    if(conflict) {
      ++stats.num_failures;
      break;
    }
  }


  return conflict;
}
  
inline Clause* SatSolver::update_watcher(const int cw, const Lit p)
{
  Clause *cl = is_watched_by[p][cw];
  Clause& clause = *cl;
  unsigned int j;
  Lit q, r;
  Atom v, w;

#ifdef _DEBUG_WATCH
  std::cout << "update watchers for " << clause 
	    << " because " << (SIGN(p) ? "" : "~") << UNSIGNED(p)
	    << " <-> b" << UNSIGNED(p) << " in " 
	    << (SIGN(p) ? "{0}" : "{1}") << std::endl;
#endif
// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
//   std::cout << " ";
//   print_clause( std::cout, is_watched_by[p][cw] );
//   std::cout << std::endl;
// #endif

  //ensure that p is the second watched lit
  if( clause[1] != p ) {
    q = clause[1];
    clause[0] = q;
    clause[1] = p;
  } else q = clause[0];
  v=state[UNSIGNED(q)];

// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
//    std::cout << " second watched: ";
//    print_literal(std::cout, q);
//    std::cout << std::endl;

//    if( LEVEL(v) < assumptions.size && SIGN(v) == SIGN(q) ) {
//      for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
//      std::cout << " satisfied by second watched" << std::endl;
//    }
// #endif

  //check if the other watched lit is assigned
  if( LEVEL(v) >= assumptions.size || SIGN(v) != SIGN(q) ) {

#ifdef _DEBUG_WATCH    
    std::cout << "  the second watcher does not satisfy the clause, we need a replacement" << std::endl;
#endif

    for(j=2; j<clause.size; ++j) {
      // for each literal q of the clause,
      r = clause[j];
      w = state[UNSIGNED(r)];

#ifdef _DEBUG_WATCH
      std::cout << "    what about " << (SIGN(r) ? "" : "~") << UNSIGNED(r)
		<< " <-> b" << UNSIGNED(r) << " in " << (LEVEL(w) >= assumptions.size ? "{0,1}" : 
							 (SIGN(w) ? "{1}" : "{0}")) << std::endl; 
#endif
// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
//        std::cout << "\tcheck another lit: ";
//        print_literal(std::cout, r);
//        std::cout << std::endl;
// #endif

      if( LEVEL(w) >= assumptions.size ) { // this literal is not set
	// then it is a good candidate to replace p
	clause[1] = r;
	clause[j] = p;
	is_watched_by[p].remove(cw);
	is_watched_by[r].add(cl);

#ifdef _DEBUG_WATCH
	std::cout << "    ok!" // << clause << " " << (cl)
		  << std::endl;
#endif
// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// 	std::cout << " new watched" << std::endl;
// #endif

	break;	
      }
      // if it is set true, then the clause is satisfied
      else if( SIGN(w) == SIGN(r) ) {

#ifdef _DEBUG_WATCH
	std::cout << "    ok! (satisfied)" << std::endl;
#endif
// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// 	std::cout << " satisfied" << std::endl;
// #endif

	break;
      }
    }
      
    if( j == clause.size ) // no replacement could be found
      { 

#ifdef _DEBUG_WATCH
	std::cout << "  couldn't find a replacement!" << std::endl;
#endif

// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// 	std::cout << " no more watched" << std::endl;
// #endif

	if( LEVEL(v) >= assumptions.size ) {
	  // the last literal (other watched lit) is not set yet, we set it
	  add_lit(q);
	  reason[UNSIGNED(q)] = cl;

#ifdef _DEBUG_WATCH
	  std::cout << "    -> b" << UNSIGNED(q) << " in " << (SIGN(q) ? "{1}" : "{0}") << std::endl;
#endif
// 	  //#ifdef _DEBUGNOGOOD
// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// 	  std::cout //<< "unit prune disjunct b" << y 
// 	    << " because ";
// 	  print_clause( std::cout, cl );
// 	  std::cout << std::endl;
// #endif

	  //	  std::cout << " prop " << q << std::endl;
	  
	} else 
	  // it is set to false already, we fail
	  if( // polarity[y] == -q
	     SIGN(v) != SIGN(q)
	      ) {

#ifdef _DEBUG_WATCH
	    std::cout << "    -> fail!" << std::endl;
#endif
// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// 	    std::cout << " fail" << std::endl;
// #endif

	    return cl;
	  }
      }
  }

  return NULL;
}






#endif
