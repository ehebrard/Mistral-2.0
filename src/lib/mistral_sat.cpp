
#include <mistral_sat.hpp>


#include <fstream>
#include <stdio.h>
#include <stdlib.h>
 

double *sorting_array;
int compar(const void *a, const void *b)
{
  if (sorting_array[*(int*)a]==sorting_array[*(int*)b])
    return 0;
  else
    if (sorting_array[*(int*)a] < sorting_array[*(int*)b])
      return 1;
    else
      return -1;
}
void initSort(double *sa)
{
  sorting_array = sa;
}

using namespace Mistral;
using namespace std;

void SatSolver::set_parameters(SolverParameters& p) {
  params = p;
  set_policy(params.restart_policy);
}


void SatSolver::set_policy( const int policy )
{
  if( policy == LUBY )
    restart_policy = new Luby(params.restart_base);
  else
    restart_policy = new Geometric(params.restart_base, params.restart_factor);
}

SatSolver::SatSolver() 
{
  restart_policy = NULL;
}

SatSolver::SatSolver(const char* filename) 
{
  restart_policy = NULL;
  parse_dimacs(filename);
}

SatSolver::~SatSolver() 
{
  for(unsigned int i=0; i<learnt.size; ++i)
    free(learnt[i]);
  for(unsigned int i=0; i<base.size; ++i)
    free(base[i]);
  delete restart_policy;
}

int SatSolver::solve()
{
  usrand(params.seed);
  stats.start_time = getRunTime();
  if(!restart_policy)
    set_policy( params.restart_policy );
  if(params.verbosity>1)
    cout << endl 
	 << "c  ==================================[ Mistral (Sat module) ]===================================" << endl
	 << "c  |      SEARCH  STATISTICS           |                  PROBLEM STATISTICS                   |" << endl
	 << "c  |  Conflicts |       Nodes  CPUTime |  Atoms  Clauses (Size)   Learnt (Size)   Added (Size) |" << endl
	 << "c  =============================================================================================" << endl;
  int result = (unsigned int)(stats.base_avg_size)+1;
  while(1) {
    //print_decisions(std::cout);
    if(params.shuffle) shuffle();
    if(params.dynamic_value) params.value_selection = randint(6);

    result = iterative_search();
    if(restart_policy) restart_policy->reset(params.restart_limit);
    if(params.verbosity>1)
      cout << "c  | " << setw(10) << stats.num_failures
	   << " | " << setw(10) << stats.num_nodes
	   << " " << setw(9) << (getRunTime() - stats.start_time)
	   << " | " << setw(6) << (state.size - stats.literals - 1)
	   << " " << setw(8) << base.size
	   << " " << setw(5) << setprecision(3) << stats.base_avg_size     
	   << "  " << setw(8) << learnt.size
	   << " " << setw(5) << setprecision(3) << stats.learnt_avg_size ;
    forget();
    if(params.verbosity>1)
      cout << "  " << setw(7) << learnt.size
	   << " " << setw(5) << setprecision(3) << stats.learnt_avg_size 
	   << "  |" << endl;
    if(result == LIMITOUT)
      status = UNKNOWN;
    else break;
  } 
  
  cout << "c  =============================================================================================" << endl;
  if( result == SAT ) 
    cout << "c  |                                    SATISFIABLE!                                           |" << endl;
  else if( result == UNSAT ) 
    cout << "c  |                                     UNSATISFIABLE!                                        |" << endl;
  else if( result == UNKNOWN ) 
    cout << "c  |                                       UNKNOWN!                                            |" << endl;
  double cputime = (getRunTime() - stats.start_time);
  if( cputime < 0.001 ) cputime = 0.001;
  cout << "c  =============================================================================================" << endl ;
  if(params.verbosity>0)    
    cout << "c  |      Conflicts  |          Nodes  |    CPU Time  | Conflicts/s |     Nodes/s | node/confl |" << endl
	 << "c  | " << setw(14) << stats.num_failures 
	 << "  | " << setw(14) << stats.num_nodes 
	 << "  | " << setw(11) << cputime
	 << "  | " << setw(11) << int(double(stats.num_failures)/cputime) 
	 << " | " << setw(11) << int(double(stats.num_nodes)/cputime) 
	 << " | " << setw(10) << (double(stats.num_nodes)/double(stats.num_failures))
	 << " |" << endl
	 << "c  =============================================================================================" << endl << endl;
  if(params.verbosity > 3) {
    if(result == SAT)
      {
	cout << "s SATISFIABLE\nv";
	for(unsigned int i=1; i<state.size; ++i) 
	  cout << " " << SIGN(state[i]) ;
	cout << endl;
      } else cout << "s UNSATISFIABLE" << endl;
    cout << "d ASSIGNMENTS " << stats.num_nodes << endl
	 << "d FAILS " << stats.num_failures << endl
      //<< "d BACKTRACKS " << stats.num_failures << endl
	 << "d NODES/s " << int(double(stats.num_nodes)/cputime) << endl;
  }


//   if( result == SAT )
//     print_decisions(cout);

  return status;
}

std::ostream& SatSolver::display(std::ostream& o) const {
  print_clauses(o);
  return o;
}

void SatSolver::print_all(ostream& o) const
{
  print_decisions( o );
  print_clauses( o );
  print_watchers( o );
}

void SatSolver::print_watchers(ostream& o, int beg, int end) const
{
  if( beg == NOVAL ) beg = -state.size;
  if( end == NOVAL ) end = state.size;
  
  Lit l;
  for(int i=beg; i<=end; ++i)
    {
      if(i) {
	l = (i<0 ? -2*(i+1) : (2*(i-1)+1));
	print_literal(o, l);
	o << " is watched by";
	o.flush();
	for(unsigned int j=0; j<is_watched_by[l].size; ++j) {
	  assert(is_watched_by[l][j]);
	  print_clause(cout, is_watched_by[l][j]);
	}
	o << endl;
      }
    }
}

void SatSolver::print_decisions(ostream& o, bool mode) const
{
  if( mode ) {
    cout << "assumptions: ";
    for(unsigned int i=0; i<assumptions.size; ++i) {
      o << " " ;
      Lit l = ((2*assumptions[i]) | SIGN(state[assumptions[i]]));
      print_literal(o, l);
      
      if(LEVEL(state[assumptions[i]]) != i)
	o << " warning - inconsistent level: " << LEVEL(state[assumptions[i]]);
    }
    o << endl;
  } else {
    unsigned int tlvl = 0;
    unsigned int clvl = decisions[0];
    for(unsigned int i=0; i<assumptions.size; ++i) {
      if(tlvl<decisions.size && i >= clvl) clvl = decisions[++tlvl];
	
      o << LEVEL(state[assumptions[i]]) // lvl[assumptions[i]]
	<< "\t" << tlvl
	<< "\t";
      print_literal(o, (2*assumptions[i]) | SIGN(state[assumptions[i]]));
      //o << "\t" << SIGN(state[assumptions[i]]) // polarity[assumptions[i]]
      o << "\t";
      o.flush();
      if( tlvl && reason[assumptions[i]] )
	print_clause(o, reason[assumptions[i]]);
      else if(tlvl)
	cout << "decision";
      else
	cout << "data";
      cout << endl;
    }
  }
}

void SatSolver::print_literal(ostream& o, Lit l) const
{
  o << (l%2 ? "+" : "-") << UNSIGNED(l)+1;
}

void SatSolver::print_clause(ostream& o, Clause *cl) const
{
  Clause& clause = *cl;
  o //<< " " << cl 
    << "(";
  for(unsigned int i=0; i<clause.size-1; ++i) {
    print_literal(o,clause[i]);
    o << " ";
  }
  print_literal(o,clause[clause.size-1]);
  o << ")";
}

void SatSolver::print_clauses(ostream& o) const
{
  o << "base (" << base.size << "):" << endl;
  for(unsigned int i=0; i<base.size; ++i)
    {
      o << "c" << i ;
      print_clause( o, base[i] );
      o << endl;
    }
  o << endl;
  o << "learnt (" << learnt.size << "):" << endl;
  for(unsigned int i=0; i<learnt.size; ++i)
    {
      o << "c" << i ;
      print_clause( o, learnt[i] );
      o << endl;
    }
  o << endl;
}

void SatSolver::init_vars(const int n, const int m)
{
  usrand(12345);

  status = UNKNOWN;
  decisions.initialise(0,n);
  assumptions.initialise(0,n);

  state.initialise(0,n);
  reason.initialise(0,n);
  //decisions_level.init(0,n+1);
  //assumptions_size.init(0,n+1);

  is_watched_by.initialise(0,2*n);
  activity.initialise(0,2*n);

  base.initialise(0,m);
  learnt.initialise(0,m);
  
//   init_num_atoms = n;
//   init_num_base = m;

  next_deduction = 0;

  visited.initialise(0, n-1, BitSet::empt);

  for(int i=0; i<n; ++i) {
    state.add(2*i + randint(2));
    reason[i] = NULL;
    
    assumptions[i] = i;
    decisions[i] = 0;
  }

  for(int i=0; i<2*n; ++i) {
    activity[i] = 0;
  }
}

void SatSolver::init_watchers()
{
  unsigned int i, j=base.size;

  for(i=0; i<2*state.size; ++i) {
    is_watched_by[i].initialise(0,(int)(2*activity[i]));
  }

  for(i=0; i<j; ++i) {
    Clause& clause = *(base[i]);
    is_watched_by[clause[0]].add(base[i]);
    is_watched_by[clause[1]].add(base[i]);
  }
} 

void SatSolver::parse_dimacs(const char* filename) 
{
  unsigned int LARGENUMBER = 131072;
  ifstream infile( filename );
  char c=' ';
  string word;
  int N, M, l=0, lit;

  // skip comments
  infile >> c;
  while( c != 'p' ) {
    infile.ignore( LARGENUMBER, '\n' );
    infile >> c;
  }

  infile >> word;
  assert( word == "cnf" );
  
  // get number of atoms and clauses
  infile >> N;
  infile >> M;

  //init(N, M);
  init_vars(N, M);
  learnt_clause.initialise(0, N);

  for(int i=0; i<M; ++i)
    {
      learnt_clause.clear();
      do {
	infile >> l;
	if(l!=0) {
	  if(l>0) lit = (l-1)*2+1;
	  else lit = (l+1)*-2;
	  learnt_clause.add(lit);
	  
	  if(params.init_activity == 1)
	    activity[lit] += params.activity_increment;
	}
      } while(l && infile.good());
      add_clause( learnt_clause );
      if(params.checked) add_original_clause( learnt_clause );
    }

  init_watchers();

  if(params.normalize_activity != 0)
    normalize_activity(params.normalize_activity);
}

int SatSolver::check_solution()
{
  if(params.checked) {
    Atom x;
    bool correct = true;
    
    for(unsigned int i=0; correct && i<original.size; ++i) {
      bool satisfied = false;
      Clause& clause = *(original[i]);
      
      for(unsigned int j=0; !satisfied && j<clause.size; ++j) {
	x = ATOM(clause[j]);
	satisfied = (SIGN(x) == SIGN(clause[j]));
      }
      if(!satisfied) {
	print_clause(cerr, original[i]);
	cerr << endl;
      }
      
      correct = satisfied;
    }
    if(!correct) {
      cerr << "/!\\ The solution is not correct /!\\" << endl; 
      exit(1);
      return UNKNOWN;
    }
  }
  return SAT;
}


Mistral::ConstraintClauseBase::ConstraintClauseBase(Vector< Variable >& scp) 
  : Constraint(scp) { }

void Mistral::ConstraintClauseBase::mark_domain() {
  for(unsigned int i=0; i<scope.size; ++i) {
    solver->mark_non_convex(scope[i].id());
  }
}

void Mistral::ConstraintClauseBase::initialise() {
  Constraint::initialise();
  for(unsigned int i=0; i<scope.size; ++i) {
    trigger_on(_VALUE_, i);
  }
  set_idempotent(true);

  is_watched_by.initialise(0,2*scope.size);
  reason.initialise(0,scope.size);
  

  for(unsigned int i=0; i<scope.size; ++i) {
    reason.add(NULL);
  //   is_watched_by[i] = new Vector< Clause* >;
  //   is_watched_by[i]-=scope[i].get_min();
  }
}

Mistral::ConstraintClauseBase::~ConstraintClauseBase() {
}

void Mistral::ConstraintClauseBase::add(Variable x) {
  unsigned int idx = x.id();
  if(idx == scope.size) {
    scope.add(x);
    reason.add(NULL);
    while(is_watched_by.capacity <= 2*idx)
      is_watched_by.extendStack();
  } else if(idx > scope.size) {
    while(scope.capacity <= idx)
      scope.extendStack();
    scope[idx] = x;

    while(reason.capacity <= idx)
      reason.extendStack();
    reason[idx] = NULL;

    while(is_watched_by.capacity <= idx)
      is_watched_by.extendStack();
  }
}

void Mistral::ConstraintClauseBase::add( Vector < Lit >& clause ) {
 if(clause.size > 1) {
   Clause *cl = Clause::Array_new(clause);
   clauses.add( cl );
   is_watched_by[clause[0]].add(cl);
   is_watched_by[clause[1]].add(cl);
 } else {
   scope[UNSIGNED(clause[0])].set_domain(SIGN(clause[0]));
 }
}

int Mistral::ConstraintClauseBase::check( const int* sol ) const {
  unsigned int i, j;
  bool falsified=false;

  for(i=0; !falsified && i<clauses.size; ++i) {
    falsified = true;
    Clause& clause = *(clauses[i]);

    std::cout << clause << " " ;
    for(j=0; j<clause.size; ++j) {
      std::cout << " " << (sol[UNSIGNED(clause[j])]) ;
    }
    std::cout << std::endl;

    for(j=0; falsified && j<clause.size; ++j) {
      falsified = //clause[j].check(sol[j]);
	(sol[UNSIGNED(clause[j])] != (int)SIGN(clause[j]));
    }
  }
  
  return falsified;
}

Mistral::PropagationOutcome Mistral::ConstraintClauseBase::propagate() {
  Clause *conflict=NULL;
  PropagationOutcome wiped = CONSISTENT;

  int x, v, cw;
  Lit p;
  //unsigned int i;
  while( !changes.empty() ) {
    //std::cout << changes << std::endl;
    //std::cout << scope << std::endl;
    x = changes.pop();
    //std::cout << scope[x] << " in " << scope[x].get_domain() << std::endl;
    v = scope[x].get_min();

    p = NEG(2*x+v);

    //std::cout << x << "=" << v << ": " << is_watched_by[p] << std::endl;

    cw = is_watched_by[p].size;
    while(cw-- && !conflict) {
      conflict = update_watcher(cw, p, wiped);
    }
  }


  // if(scope[112].is_ground() && scope[243].is_ground() && scope[77].is_ground())
  //   {
  //     exit(1);
  //   }

  // if(conflict) {

  //     // here we should do some conflict analysis
  // }

  return wiped;
}

#define _DEBUG_WATCH true

inline Clause* ConstraintClauseBase::update_watcher(const int cw, 
						    const Lit p,
						    PropagationOutcome& po)
{
  Clause *cl = is_watched_by[p][cw];
  Clause& clause = *cl;
  unsigned int j;


  // BitSet e(0, scope.size, BitSet::empt);
  // for(j=0; j<clause.size; ++j) 
  //   {
  //     e.add(UNSIGNED(clause[j]));
  //   }
  // std::cout << e << std::endl;


  Lit q, r;
  //Atom v, w;
  Variable v, w;

#ifdef _DEBUG_WATCH
  std::cout << "update watchers for " << clause 
	    << " because " << (SIGN(p) ? "" : "~") << UNSIGNED(p)
	    << " <-> " << scope[UNSIGNED(p)] << " in " 
	    << scope[UNSIGNED(p)].get_domain() << std::endl;
#endif

  //ensure that p is the second watched lit
  if( clause[1] != p ) {
    q = clause[1];
    clause[0] = q;
    clause[1] = p;
  } else q = clause[0];
  v=scope[UNSIGNED(q)];

  //check if the other watched lit is assigned
  if( !v.is_ground() || v.get_min() != (int)SIGN(q) ) {

#ifdef _DEBUG_WATCH    
    std::cout << "  the second watcher does not satisfy the clause, we need a replacement" << std::endl;
#endif

    for(j=2; j<clause.size; ++j) {
      // for each literal r of the clause,
      r = clause[j];
      w = scope[UNSIGNED(r)];

#ifdef _DEBUG_WATCH
      std::cout << "    what about " << (SIGN(r) ? "" : "~") << UNSIGNED(r)
		<< " <-> " << w << " in " << w.get_domain() << std::endl; 
#endif

      if( !w.is_ground() ) { // this literal is not set
	// then it is a good candidate to replace p

	clause[1] = r;
	clause[j] = p;
	is_watched_by[p].remove(cw);
	is_watched_by[r].add(cl);

#ifdef _DEBUG_WATCH
	std::cout << "    ok!" // << clause << " " << (cl)
		  << std::endl;
#endif

	break;	
      }
      // if it is set true, then the clause is satisfied
      else if( w.get_min() == (int)SIGN(r) ) {

#ifdef _DEBUG_WATCH
	std::cout << "    ok! (satisfied)" << std::endl;
#endif

	break;
      }
    }
      
    if( j == clause.size ) // no replacement could be found
      { 

#ifdef _DEBUG_WATCH
	std::cout << "  couldn't find a replacement!" << std::endl;
#endif

	if( !v.is_ground() ) {
	  // the last literal (other watched lit) is not set yet, we set it
	  //add_lit(q);
	  v.set_domain(SIGN(q));
	  reason[UNSIGNED(q)] = cl;
	  changes.add(UNSIGNED(q));

#ifdef _DEBUG_WATCH
	  std::cout << "    -> " << v << " in " << v.get_domain() << std::endl;
#endif

	} else 
	  // it is set to false already, we fail
	  if( v.get_min() != (int)SIGN(q) ) {

#ifdef _DEBUG_WATCH
	    std::cout << "    -> fail!" << std::endl;
#endif
	    po = FAILURE(UNSIGNED(q));

	    return cl;
	  }
      }
  }

  return NULL;
}

std::ostream& Mistral::ConstraintClauseBase::display(std::ostream& os) const {
  // os << " (";
  // if(clauses.size>0) {
  //   os << clauses[0];
  //   for(unsigned int i=1; i<clauses.size; ++i)
  //     os << " " << clauses[i]  ;
  // }
  // os << ")";
  os << "nogoods";
  return os;
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Clause& x) {
  os << "(" << (SIGN(x[0]) ? "" : "~") << UNSIGNED(x[0]) ;
  for(unsigned int i=1; i<x.size; ++i) {
    os << " " << (SIGN(x[i]) ? "" : "~") << UNSIGNED(x[i]) ;
  }
  os << ")";
  return os;
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Clause* x) {
  return os << (*x);
}

std::ostream& Mistral::operator<< (std::ostream& os, Mistral::SatSolver& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, Mistral::SatSolver* x) {
  return x->display(os);
}

