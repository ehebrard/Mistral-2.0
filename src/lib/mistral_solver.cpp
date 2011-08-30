
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

#include <sstream>
#include <fstream>
#include <signal.h>
#include <assert.h>

#include <mistral_sat.hpp>
#include <mistral_solver.hpp>
#include <mistral_variable.hpp>
#include <mistral_constraint.hpp>


Mistral::Solver* active_solver;
static void Mistral_SIGINT_handler(int signum) {
  std::cout << std::endl 
	    << " c ********************************* INTERRUPTED *********************************" 
	    << std::endl;
  std::cout << active_solver->statistics << std::endl
	    << " c ********************************* INTERRUPTED *********************************" 
	    << std::endl;
  exit(1);
}

#ifdef _DEBUG_PROPAG 
std::ostringstream *o_propag = NULL;
#endif

Mistral::Solution::Solution( Vector< Variable >& vars ) {
  min_id = INFTY;
  max_id = -INFTY;
  int aux;
  for(unsigned int i=0; i<vars.size; ++i) {
    variables.add(vars[i].get_var());
    aux = variables[i].id();
    if(aux<min_id) min_id = aux;
    if(aux>max_id) max_id = aux;
  }

  values = new int[max_id-min_id+1];
  values -= min_id;

  for(unsigned int i=0; i<variables.size; ++i) {
    values[variables[i].id()] = variables[i].get_solution_int_value();
  }
}
Mistral::Solution::~Solution() {
  values += min_id;
  delete [] values;
}

std::ostream& Mistral::Solution::display(std::ostream& os) const {
  if(!variables.empty()) {
    os << variables[0] << ":" << values[variables[0].id()];
    for(unsigned int i=1; i<variables.size; ++i) {
      os << " " << variables[i] << ":" << values[variables[i].id()];
    }
  }
  return os;
}
   
int& Mistral::Solution::operator[](Variable x) {
  return values[x.id()];
}


Mistral::SolverParameters::SolverParameters() {initialise();}
Mistral::SolverParameters::~SolverParameters() {}
void Mistral::SolverParameters::initialise() {
  find_all = 0;
  node_limit = 0;
  backtrack_limit = 0;
  fail_limit = 0;
  //restart_limit = 0;
  limit = 0;
  seed = 11041979;

  verbosity = 0;
  time_limit = -1;
  seed = 11041979;
  restart_policy = GEOMETRIC;
  restart_base = 256;
  restart_limit = 256;
  restart_factor = 1.33;
  activity_increment = 1e-2;
  normalize_activity = 0;
  init_activity = 1;
  forgetfulness = .75;
  randomization = 1; //2;
  shuffle = false; //true;
  activity_decay = 0.96;
  checked = 1;
  backjump = 0;
  value_selection = 0; //6;
  dynamic_value = 0; //1;
}
Mistral::SolverParameters::SolverParameters(const SolverParameters& sp) {
  copy(sp);
}
void Mistral::SolverParameters::copy(const SolverParameters& sp) {
  verbosity = sp.verbosity;
  find_all = sp.find_all;
  node_limit = sp.node_limit;
  backtrack_limit = sp.backtrack_limit;
  fail_limit = sp.fail_limit;
  restart_limit = sp.restart_limit;
  limit = sp.limit;
  time_limit = sp.time_limit;
}

Mistral::SolverStatistics::SolverStatistics() { initialise(); }
Mistral::SolverStatistics::~SolverStatistics() {}
void Mistral::SolverStatistics::initialise() {
  num_nodes = -1; 
  num_restarts = 0; 
  num_backtracks = 0;
  num_failures = 0; 
  num_propagations = 0;
  num_solutions = 0;
  num_filterings = 0;
  start_time = 0.0;
  end_time = -1.0;

  base_avg_size = 0;
  learnt_avg_size = 0;
  literals = 0;
  small = 0;

}
std::ostream& Mistral::SolverStatistics::print_full(std::ostream& os) const {
  os << " c +" << std::setw(78) << std::setfill('=')
    //"=============================================================================
     << "+" << std::endl << std::setfill(' ')
     << std::left << std::setw(40) << " s  ";

  switch(outcome) {
  case SAT: 
    os << std::right << std::setw(42) << "SAT |" ;
    break;
  case OPT: 
    os << std::right << std::setw(42) << "OPT |" ;
    break;
  case UNSAT: 
    os << std::right << std::setw(42) << "UNSAT |" ;
    break;
  case UNKNOWN: 
    os << std::right << std::setw(42) << "UNKNOWN |" ;
    break;
  case LIMITOUT: 
    os << std::right << std::setw(42) << "LIMITOUT |" ;
    break;
  }
  os << std::endl
     << std::left << std::setw(40) << " d  TIME"
     << std::right << std::setw(40) << (end_time - start_time) << " |" << std::endl
     << std::left << std::setw(40) << " d  NODES"
     << std::right << std::setw(40) << num_nodes << " |" << std::endl
     << std::left << std::setw(40) << " d  RESTARTS"
     << std::right << std::setw(40) << num_restarts << " |" << std::endl
     << std::left << std::setw(40) << " d  FAILURES"
     << std::right << std::setw(40) << num_failures << " |" << std::endl
     << std::left << std::setw(40) << " d  BACKTRACKS"
     << std::right << std::setw(40) << num_backtracks << " |" << std::endl
     << std::left << std::setw(40) << " d  PROPAGATIONS"
     << std::right << std::setw(40) << num_propagations << " |" << std::endl
     << std::left << std::setw(40) << " d  FILTERINGS"
     << std::right << std::setw(40) << num_filterings << " |" << std::endl
     << " c +" << std::setw(78) << std::setfill('=') << "+" << std::endl << std::setfill(' ');
  //<< " c +=============================================================================+" << std::endl;
  return os;
}
std::ostream& Mistral::SolverStatistics::print_short(std::ostream& os) const {
  os << " c |"
     << std::right << std::setw(7) << num_variables << " |"
     << std::right << std::setw(8) << num_values << " |"
     << std::right << std::setw(7) << num_constraints << " |"
     << std::right << std::setw(9) << num_nodes << " |" 
     << std::right << std::setw(11) << num_filterings << " |" 
     << std::right << std::setw(13) << num_propagations << " |" 
     << std::right << std::setw(9) << std::setprecision(5) 
     << (get_run_time() - start_time) << " |" ;
  return os;
}
std::ostream& Mistral::SolverStatistics::display(std::ostream& os) const {
  if(end_time >= 0.0) print_full(os);
  else print_short(os);
  return os;
}
Mistral::SolverStatistics::SolverStatistics(const SolverStatistics& sp) {
  copy(sp);
}
void Mistral::SolverStatistics::copy(const SolverStatistics& sp) {
  num_nodes = sp.num_nodes; 
  num_restarts = sp.num_restarts; 
  num_backtracks = sp.num_backtracks;
  num_failures = sp.num_failures; 
  num_propagations = sp.num_propagations;
  num_solutions = sp.num_solutions;
  num_filterings = sp.num_filterings;
  start_time = sp.start_time;
  end_time = sp.end_time;
}
void Mistral::SolverStatistics::update(const SolverStatistics& sp) {
  num_nodes += sp.num_nodes; 
  num_restarts += sp.num_restarts; 
  num_backtracks += sp.num_backtracks;
  num_failures += sp.num_failures; 
  num_propagations += sp.num_propagations;
  num_solutions += sp.num_solutions;
  num_filterings += sp.num_filterings;
  if(end_time < sp.end_time) end_time = sp.end_time;
}

Mistral::ConstraintQueue::ConstraintQueue()
{
  min_priority = 0;
  cardinality = -1;
  higher_priority  = -1;
  triggers = NULL;
  taboo_constraint = NULL;
  _set_.initialise();
}

void Mistral::ConstraintQueue::declare(Constraint *c, Solver *s) {

  // std::cout << "declare " << c << "(" << c->priority << ") to the constraint queue" << std::endl;
  //  std::cout << "was: [" << min_priority << ","
  // 	    << min_priority+cardinality-1 << "]" << std::endl;

  int cons_idx = c->id;
  int cons_priority = c->priority;

  int new_min_p = min_priority;
  int new_max_p = min_priority+cardinality-1;

  if(cons_idx == 0) solver = s;
  if(cons_priority < new_min_p || cons_priority > new_max_p) {
    if(cardinality > 0) {
      if(cons_priority < new_min_p) new_min_p = cons_priority;
      if(cons_priority > new_max_p) new_max_p = cons_priority;
    } else {
      new_min_p = cons_priority;
      new_max_p = cons_priority;
    }

//     std::cout << "now: [" << new_min_p << ","
// 	    << new_max_p << "]" << std::endl;


    int new_cardinality = (new_max_p-new_min_p+1);
    Queue *aux_triggers = triggers;
    triggers = new Queue[new_cardinality];
    triggers -= new_min_p;
    for(int i=min_priority; i<min_priority+cardinality; ++i) {
      triggers[i] = aux_triggers[i];
      aux_triggers[i].cancel();
    }

    triggers[cons_priority].initialise(cons_idx, cons_idx+7);

    aux_triggers += min_priority;
    delete [] aux_triggers;

    min_priority = new_min_p;
    cardinality = new_cardinality;


    // std::cout << " ==> init " << triggers[cons_priority] << std::endl;


  } else {
    
     // std::cout << "no need to create a new trigger list" << std::endl;
     // std::cout << "extend " << triggers[cons_priority] << " with " << cons_idx << std::endl;
     if(!triggers[cons_priority].is_initialised()) {
       triggers[cons_priority].initialise(cons_idx, cons_idx+7);
     } else {
       triggers[cons_priority].extend(cons_idx);
     }
  }
  
  if(_set_.table)
    _set_.extend(cons_idx);
  else 
    _set_.initialise(cons_idx, cons_idx, BitSet::empt);

}

void Mistral::ConstraintQueue::initialise(Solver *s)
{
  solver = s;
  _set_.initialise(0,2*s->constraints.capacity,BitSet::empt);
  
  int min_p = INFTY;
  int max_p = -INFTY;
  for(unsigned int i=0; i<s->constraints.size; ++i) {
    if(s->constraints[i]->priority < min_p) min_p = s->constraints[i]->priority;
    if(s->constraints[i]->priority > max_p) max_p = s->constraints[i]->priority;
  }
  initialise(min_p, max_p, s->constraints.size);
}

void Mistral::ConstraintQueue::initialise(const int min_p, 
					  const int max_p,
					  const int size)
{
  cardinality = max_p-min_p+1;
  min_priority = min_p;
  triggers = new Queue[cardinality];
  for(int i=0; i<cardinality; ++i) {
    triggers[i].initialise(size);
  }
  triggers -= min_priority;
}

Mistral::ConstraintQueue::~ConstraintQueue() {
  triggers += min_priority;
  delete [] triggers;
}

void Mistral::ConstraintQueue::trigger(Constraint* cons)//;
{
#ifdef _DEBUG_AC
  std::cout << " initial trigger for " << cons << "(" << (cons->id) << ")" << std::endl;
#endif

  int priority = cons->priority, cons_id = cons->id, triggered=false;
  Event evt;
  
  if(!_set_.fast_contain(cons_id)) {
    cons->events.clear();
    for(unsigned int i=0; i<cons->scope.size; ++i) {
      evt = (// cons->scope[i].domain_type != BOOL_VAR &&
	     cons->scope[i].is_ground() ? VALUE_EVENT : (LB_EVENT|UB_EVENT));

      //std::cout << EVENT_TYPE(evt) << " <=? " << cons->trigger[i] << std::endl;

      if(EVENT_TYPE(evt) <= cons->trigger[i]) {
	cons->events.add(i);
	cons->event_type[i] = evt;
	triggered=true;
      }
    }
    if(triggered) {

      //std::cout << "TRIGGER!" << std::endl;

      _set_.fast_add(cons_id);
      if(priority > higher_priority) higher_priority = priority;
      triggers[priority].add(cons_id);
    }  
  } else {


    //std::cout << "ALREADY IN!" << std::endl;

    for(unsigned int var=0; var<cons->scope.size; ++var) {
      evt = (// cons->scope[i].domain_type != BOOL_VAR &&
	     cons->scope[var].is_ground() ? VALUE_EVENT : (LB_EVENT|UB_EVENT));
      if(cons->events.contain(var)) {
	cons->event_type[var] |= evt;
      } else {
	cons->events.add(var);
	cons->event_type[var] = evt;
      }
    }
  }

  //std::cout << "EVENTS: " << cons->events << std::endl;

  //std::cout << _set_ << std::endl;
}

void Mistral::ConstraintQueue::trigger(Constraint* cons, const int var, const Event evt)//;
{

#ifdef _DEBUG_AC
  std::cout << "  triggers " << cons << " after a " 
	    << (ASSIGNED(evt) ? "value" : (RANGE_CHANGED(evt) ? "range" : "domain"))
	    << " event on "
	    << cons->scope[var] << " in " << cons->scope[var].get_domain() 
	    << std::endl;
#endif


  if(cons != taboo_constraint) {
    int priority = cons->priority, cons_id = cons->id;
    if(_set_.fast_contain(cons_id)) {
      if(cons->events.contain(var)) {
	cons->event_type[var] |= evt;
      } else {
	cons->events.add(var);
	cons->event_type[var] = evt;
      }
    } else {
      _set_.fast_add(cons_id);
      if(priority > higher_priority) higher_priority = priority;
      triggers[priority].add(cons_id);
      cons->events.set_to(var);
      cons->event_type[var] = evt;
    }
  } 

#ifdef _DEBUG_QUEUE

  else
    std::cout << cons << "(" << (cons->id) 
	      << ") is currently being processed " 
	      << std::endl;
      
#endif


#ifdef _DEBUG_QUEUE

  std::cout << "=> " << cons->events << std::endl;

#endif

}


std::ostream& Mistral::ConstraintQueue::display(std::ostream& os) {
  int elt;
  int end;

  for(int i=0; i<cardinality; ++i) {
    if(!triggers[i+min_priority].empty()) {
      elt = triggers[i+min_priority].first();
      end = triggers[i+min_priority]._head;
      os << "P" << (i+min_priority) << ": " 
	//<< solver->constraints[elt];
	 << "["<< solver->constraints[elt]->id << "]";

      if(!_set_.contain(solver->constraints[elt]->id)) {
	std::cout << "inconsistent constraint queue" <<std::endl;
	exit(1);
      }

      while(triggers[i+min_priority].next[elt] != end) {
	elt = triggers[i+min_priority].next[elt];
	os << ", " //<< solver->constraints[elt];
	   << " [" << solver->constraints[elt]->id << "]";

	if(!_set_.contain(solver->constraints[elt]->id)) {
	  std::cout << "inconsistent constraint queue" <<std::endl;
	  exit(1);
	}
      }
      os << std::endl;
    }
  }
  return os;
}

Mistral::Solver::Solver() { 
  search_started = false;

  // search stuf
  heuristic = NULL;
  policy = NULL;
  objective = NULL; //new Goal(Goal::SATISFACTION);
  //backjump_policy = NULL;

  // variables & constraints
  domain_types.initialise(0,128);
  variables.initialise(0,128);
  assignment_level.initialise(0,128);
  visited.initialise(0,1023);
  reason.initialise(0,128);
  constraints.initialise(0,256);
  //constraints.initialise(256);
  posted_constraints.initialise(0,256,false);
  sequence.initialise(128);
  initialised_vars = 0;
  initialised_cons = 0;
  num_search_variables = 0;
  base = NULL;

  // trail stuff
  level = -1;
  saved_objs.initialise(0,4096); 
  saved_vars.initialise(0,4096); 
  saved_cons.initialise(0,4096); 
  trail_.initialise    (0,4096);
  decisions.initialise (0,4096);
  con_trail_.initialise(0,512);
  con_trail_.add(0);
  con_trail_.add(-1);

  // params
  parameters.initialise();

  // statistics
  statistics.initialise();

  heuristic = NULL; //new GenericHeuristic< GenericDVO< MinDomain >, MinValue >(this);
  policy = NULL; //new Geometric();

  usrand(parameters.seed);

  save();
}

void Mistral::Solver::parse_dimacs(const char* filename) {
  unsigned int LARGENUMBER = 131072;
  std::ifstream infile( filename );
  char c=' ';
  std::string word;
  int N, M, l=0;
  Lit lit;
  Vector< Lit > new_clause;

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

  for(int i=0; i<N; ++i) {
    Variable x(0,1);
    add(x);
  }

  new_clause.initialise(0,N);
  //ConstraintClauseBase *cbase = new ConstraintClauseBase(variables);
  //add(cbase);

  for(int i=0; i<M; ++i)
    {
      new_clause.clear();
      do {
	infile >> l;
	if(l!=0) {
	  if(l>0) lit = (l-1)*2+1;
	  else lit = (l+1)*-2;
	  new_clause.add(lit);
	  // if(parameters.init_activity == 1)
	  //   base->activity[lit] += parameters.activity_increment;
	}
      } while(l && infile.good());
      //cbase->add( new_clause );
      add( new_clause );

      //std::cout << new_clause << std::endl;
      
      //if(params.checked) add_original_clause( new_clause );
    }
  
  //init_watchers();
  
  //if(params.normalize_activity != 0)
  //normalize_activity(params.normalize_activity);
  
  //  std::cout << base << std::endl;

  //cbase->reason = reason.stack_;
}


void Mistral::Solver::add(Vector< Lit >& clause) {
  if(!base) {
    base = new ConstraintClauseBase(variables);
    add(base);
    base->reason = reason.stack_;
  }
  base->add( clause, (parameters.init_activity ? parameters.activity_increment : 0.0) );
}

void Mistral::Solver::set_parameters(SolverParameters& p) {
  parameters = p;
}

void Mistral::Solver::add(VarArray& x) {
  for(unsigned int i=0; i<x.size; ++i)
    x[i].initialise(this);
}

void Mistral::Solver::add(Variable x) { 

  //std::cout << "add a variable " << x << std::endl;

  x.initialise(this); 
}

int Mistral::Solver::declare(Variable x) {
  if(x.domain_type > DYN_VAR) booleans.add(&x);

  // add the variables to the set of vars
  x.variable->id = variables.size;
  x.variable->solver = this;
  visited.extend(variables.size);
  variables.add(x);
  assignment_level.add(INFTY);
  reason.add(NULL);
  domain_types.add(DYN_VAR|(x.is_range() ? RANGE_VAR : 0));

  last_solution_lb.add(-INFTY);
  last_solution_ub.add( INFTY);
  
  // create a new empty list of constraints
  ConstraintList *l = new ConstraintList(this);
  constraint_graph.add(l);

  notify_add_variable();

  return variables.size-1;
}


// void Mistral::Solver::add(Constraint* c) { 
//   for()
// }

void Mistral::Solver::add(Constraint* c) { 
  if(c->id < 0) {
    // //if(parameters.verbosity>0) {
    // //std::cout << "c add a constraint: " << c << std::endl; 
    //   //}

    for(unsigned int i=0; i<c->scope.size; ++i) {
      c->scope[i].initialise(this, false);
    }

    c->initialise();

    // get a name for the constraint and add it to the list
    c->id = constraints.size; 
    //constraints.declare(c);
    constraints.add(c); 

    //
    active_constraints.declare(c, this);

  }

  // then post it on the constraint graph
  c->post(this);
  active_constraints.trigger(c);

  // if(!active_constraints._set_.contain(c->id)) {

  //   std::cout << "problem when triggering " << c << std::endl;
  //   //exit(1);
  // }
  

  if(con_trail_.back() != level) {
    con_trail_.add(posted_constraints.size);
    con_trail_.add(level);
  }

  if(!posted_constraints.safe_contain(c->id)) {
    posted_constraints.extend(c->id);
    posted_constraints.add(c->id);
  }
}

Mistral::Outcome Mistral::Solver::solve() {
  BranchingHeuristic *heu = new GenericHeuristic <
    GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >,
    RandomMinMax 
    > (this); 
  RestartPolicy *pol = new Geometric();
  return depth_first_search(variables, heu, pol); 
  //return (search_outcome == SAT || search_outcome == OPT);
}


Mistral::Outcome Mistral::Solver::minimize(Variable X) {
  BranchingHeuristic *heu = new GenericHeuristic <
    GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >,
    RandomMinMax 
    > (this); 
  RestartPolicy *pol = new Geometric();
  Goal *goal = new Goal(Goal::MINIMIZATION, X.get_var());
  return depth_first_search(variables, heu, pol, goal); 
  //return (search_outcome == OPT);
}


Mistral::Outcome Mistral::Solver::maximize(Variable X) {
  BranchingHeuristic *heu = new GenericHeuristic <
    GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >,
    RandomMinMax 
    > (this); 
  RestartPolicy *pol = new Geometric();
  Goal *goal = new Goal(Goal::MAXIMIZATION, X.get_var());
  return depth_first_search(variables, heu, pol, goal); 
  //return (search_outcome == OPT);
}

Mistral::Outcome Mistral::Solver::depth_first_search(BranchingHeuristic *heu, 
						     RestartPolicy *pol,
						     Goal *goal) {
  return depth_first_search(variables, heu, pol, goal);
}

Mistral::Outcome Mistral::Solver::depth_first_search(Vector< Variable >& seq, 
						     BranchingHeuristic *heu, 
						     RestartPolicy *pol,
						     Goal *goal) 
{
  initialise_search(seq, heu, pol, goal);

  Outcome satisfiability = UNKNOWN;
  while(satisfiability == UNKNOWN) {

    statistics.num_variables = sequence.size;
    statistics.num_values = 0;
    for(unsigned int i=0; i<sequence.size; ++i)
      statistics.num_values += sequence[i].get_size();

    if(parameters.verbosity) {
      statistics.print_short(std::cout);
      //if(base) std::cout << " " << base->learnt.size ;
      std::cout << std::endl;
    }

    ++statistics.num_restarts;
    satisfiability = chronological_dfs(); 
    //conflict_directed_backjump(); //

    if(satisfiability == LIMITOUT) {
      policy->reset(parameters.restart_limit);    
      if(!limits_expired()) satisfiability = UNKNOWN;
      forget();
    }

    restore(0);
    decisions.clear();
  }

  statistics.outcome = satisfiability;
  statistics.end_time = get_run_time();

  if(parameters.verbosity)  {
    // switch(satisfiability) {
    // case UNKNOWN: std::cout << " s  UNKNOWN" << std::endl; break;
    // case SAT: std::cout << " s  SATISFIABLE" << std::endl; break;
    // case UNSAT: std::cout << " s  UNSATISFIABLE" << std::endl; break;
    // case OPT: std::cout << " s  OPTIMAL" << std::endl; break;
    // }
    std::cout << statistics << std::endl;
  }

  return satisfiability;
}

Mistral::Outcome Mistral::Solver::get_next_solution()  
{
  Outcome satisfiability = UNSAT;

//   std::cout << "get next solution " << (decisions.size) << " "
// 	    << search_started << std::endl;

  if(search_started) {
    if(decisions.size) 
      branch_right();
    else return satisfiability;
  }
   
  search_started = true;
  
  statistics.num_variables = sequence.size;
  statistics.num_values = 0;
  for(unsigned int i=0; i<sequence.size; ++i)
    statistics.num_values += sequence[i].get_size();
  
  //display(std::cout);
  satisfiability = chronological_dfs();
  
  if(parameters.verbosity) {
    statistics.print_short(std::cout);
  }

  return satisfiability;
}

void Mistral::Solver::BooleanMemoryManager::add(Variable *x) {

  if(size.back() < 1024) {
    x->bool_domain = slots.back()+size.back();
    ++size.back();
  } else {
    int *nslot = new int[1024];
    std::fill(nslot, nslot+1024, 3);
    size.add(1);
    slots.add(nslot);
    x->bool_domain = nslot;
  }


  //std::cout << "zzz " << *x << ": " << x->domain_type << std::endl;
}


void Mistral::Solver::initialise_search(Vector< Variable >& seq, 
					BranchingHeuristic *heu, 
					RestartPolicy *pol,
					Goal *goal) 
{
  consolidate();

  if(level < 0) save();

  active_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  if(statistics.start_time == 0.0) statistics.start_time = get_run_time();

  sequence.clear();
  decisions.clear();
  for(unsigned int i=seq.size; i;) {
    Variable x = seq[--i].get_var();
    if(!x.is_ground() && !sequence.contain(x) && !(domain_types[x.id()]&REMOVED_VAR)) sequence.add(x);
  }
  num_search_variables = sequence.size;


  // for(int i=0; i<sequence.size; ++i) {
  //   std::cout << sequence[i] << " in " << sequence[i].get_domain() << std::endl;
  // }
  // std::cout << std::endl;
  // //exit(1);

  if(heu) { delete heuristic; heuristic = heu; }
  else if(!heuristic) heuristic = new GenericHeuristic< Lexicographic, MinValue >(this);
  if(pol) { delete policy;    policy    = pol; }
  else if(!policy)    policy    = new NoRestart();
  if(goal){ delete objective; objective = goal;}
  else if(!objective) objective = new Goal(Goal::SATISFACTION);

  heuristic->initialise(sequence);

  parameters.restart_limit = policy->base;
  parameters.limit = (policy->base > 0);

  statistics.num_constraints = constraints.size;

  if(parameters.verbosity)  std::cout << " c +" << std::setw(78) << std::setfill('=')
			      //=============================================================================
				      << "+" << std::endl << std::setfill(' ') 
				      << " c |      INSTANCE STATS       |                    SEARCH STATS                 |" << std::endl 
				      << " c |   vars |    vals |   cons |    nodes | filterings | propagations | cpu time |" << std::endl;
}



Mistral::Solver::~Solver() {

  //std::cout << "delete solver " << std::endl;  

  delete heuristic;
  delete policy;
  delete objective;

  //std::cout << "delete constraints" << std::endl;
  for(unsigned int i=0; i<constraints.size; ++i) {

    //std::cout << "  delete " << constraints[i] << std::endl;

    delete constraints[i];
  }

  //std::cout << "delete expressions" << std::endl;
  for(unsigned int i=0; i<expression_store.size; ++i) {

    //Variable x(expression_store[i]);
    //std::cout << "  delete " << expression_store[i] << std::endl;

    delete expression_store[i];
  }

  //std::cout << "delete variables" << std::endl;
  for(unsigned int i=0; i<variables.size; ++i) {

    //std::cout << "  delete " << variables[i] << " in " << variables[i].get_domain() << std::endl;

    int domain_type = variables[i].domain_type;
    if     (domain_type ==  BITSET_VAR) delete variables[i].bitset_domain;
    else if(domain_type ==    LIST_VAR) delete variables[i].list_domain;
    else if(domain_type ==   RANGE_VAR) delete variables[i].range_domain;
    else if(domain_type == VIRTUAL_VAR) delete variables[i].virtual_domain;
    else if(domain_type ==  EXPRESSION) delete variables[i].expression;
    else if(domain_type !=   CONST_VAR) delete variables[i].variable;
    delete constraint_graph[i];
  }

}


void Mistral::Solver::trigger_event(const int var, 
				    const Mistral::Event evt) {

  
#ifdef _DEBUG_QUEUE
  std::cout << (ASSIGNED(evt) ? "value" : (RANGE_CHANGED(evt) ? "range" : "domain"))
	    << " event on " << variables[var] 
	    << " in " << variables[var].get_domain() << std::endl;
  // if(var==1287)
  //   std::cout << sequence << std::endl;
#endif

  Constraint *c;
  ConstraintNode nd;

  nd = constraint_graph[var]->first(EVENT_TYPE(evt));
  if(ASSIGNED(evt)) {

    //std::cout << variables[var] << " is assigned " << std::endl;

    while(constraint_graph[var]->next(nd)) {

      //std::cout << "\t" << nd.elt.constraint << std::endl;

#ifdef _DEBUG_CGRAPH
      std::cout << "value event on " << variables[var] << ", notify " << nd.elt.constraint 
		<< " " << nd.elt.constraint->active ;
#endif

      active_constraints.trigger(nd.elt.constraint, nd.elt.index, evt);
      c = nd.elt.constraint->notify_assignment(nd.elt.index, level);

#ifdef _DEBUG_CGRAPH
      std::cout << " -> " << nd.elt.constraint->active << std::endl;
#endif

      if(c) saved_cons.add(c);

    }
    if(sequence.contain(variables[var])) {
      sequence.remove(variables[var]);
      assignment_level[var] = level;
    }
  } else while(constraint_graph[var]->next(nd)) {
      active_constraints.trigger(nd.elt.constraint, nd.elt.index, evt);
    }

  wiper_idx = var;

  // if(var==1287)
  // std::cout << sequence << std::endl;
  
}


// void Mistral::Solver::save(Variable x) { 
//   saved_vars.add(x); 
// }

// void Mistral::Solver::save(VariableImplementation *impl, int dtype) { 
//   Variable x(impl, dtype);
//   saved_vars.add(x); 
// }

void Mistral::Solver::save(const int idx) {

    // if(idx == 17) {
    //   for(int i=0; i<level; ++i)
    // 	std::cout << " ";
    //   std::cout << "save " << variables[idx] 
    // 		<< " in " << variables[idx].get_domain() 
    // 		<< " "
    // 		<< ((VariableBitmap*)(variables[idx].variable))->trail_
    // 		<< std::endl;
    // }

 
  //saved_vars.add(variables[idx]);
  saved_vars.add(idx);
}

void Mistral::Solver::restore() {
  unsigned int previous_level;

  previous_level = trail_.pop();
  while( saved_vars.size > previous_level ) {
    // if(saved_vars.back() == 4) {
    //   std::cout << level << " RESTORE X4: " << variables[4] << " in " << variables[4].get_domain() 
    // 		<< " " << (variables[4].domain_type == RANGE_VAR ? 
    // 			   ((VariableRange*)(variables[4].variable))->trail_ :
    // 			   ((VariableBitmap*)(variables[4].variable))->trail_)
    // 		<< std::endl;
    // }

    variables[saved_vars.pop()].restore();

    // if(saved_vars.back(0) == 4) {
    //   std::cout << level << " ======> X4: " << variables[4] << " in " << variables[4].get_domain() 
    // 		<< " " << (variables[4].domain_type == RANGE_VAR ? 
    // 			   ((VariableRange*)(variables[4].variable))->trail_ :
    // 			   ((VariableBitmap*)(variables[4].variable))->trail_)
    // 		<< std::endl << std::endl;
    // }
  }
  previous_level = trail_.pop();
  while( saved_cons.size > previous_level ) 
    saved_cons.pop()->restore();

  previous_level = trail_.pop();
  while( saved_objs.size > previous_level )
    saved_objs.pop()->restore();

  if(con_trail_.back() == level) {
    con_trail_.pop();
    posted_constraints.size = con_trail_.pop();
  }

  //  decisions.pop(); 
  sequence.size = trail_.pop();

  ++statistics.num_backtracks;
  --level;
}

void Mistral::Solver::restore(const int lvl) {
  while(lvl < level) restore();
}

void Mistral::Solver::add(Mistral::RestartListener* l) {
  l->rid = restart_triggers.size;
  restart_triggers.add(l);
}
void Mistral::Solver::add(Mistral::SuccessListener* l) {
  l->sid = success_triggers.size;
  success_triggers.add(l);
}
void Mistral::Solver::add(Mistral::FailureListener* l) {
  l->fid = failure_triggers.size;
  failure_triggers.add(l);
}
void Mistral::Solver::add(Mistral::DecisionListener* l) {
  l->did = decision_triggers.size;
  decision_triggers.add(l);
}
void Mistral::Solver::add(Mistral::VariableListener* l) {
  l->vid = variable_triggers.size;
  variable_triggers.add(l);
}
void Mistral::Solver::add(Mistral::ConstraintListener* l) {
  l->cid = constraint_triggers.size;
  constraint_triggers.add(l);
}


void Mistral::Solver::notify_failure() { //Constraint *con, const int idx) {
  for(unsigned int i=0; i<failure_triggers.size; ++i) {
    failure_triggers[i]->notify_failure();
  }
} 

void Mistral::Solver::notify_success() { //Variable* changes, const int n) {
  for(unsigned int i=0; i<success_triggers.size; ++i) {
    success_triggers[i]->notify_success();
  }
} 

void Mistral::Solver::notify_decision() { //Decision d) {
  for(unsigned int i=0; i<decision_triggers.size; ++i) {
    decision_triggers[i]->notify_decision();
  }
} 

void Mistral::Solver::notify_restart() { 
  for(unsigned int i=0; i<restart_triggers.size; ++i) {
    restart_triggers[i]->notify_restart();
  }
} 

void Mistral::Solver::notify_relax(Constraint *c) { 
  for(unsigned int i=0; i<constraint_triggers.size; ++i) {
    constraint_triggers[i]->notify_relax(c);
  }
} 

void Mistral::Solver::notify_post(Constraint *c) { 
  for(unsigned int i=0; i<constraint_triggers.size; ++i) {
    constraint_triggers[i]->notify_post(c);
  }
} 

void Mistral::Solver::notify_change_variable(const int idx) { 
  for(unsigned int i=0; i<variable_triggers.size; ++i) {
    variable_triggers[i]->notify_change(idx);
  }
} 

void Mistral::Solver::notify_add_variable() { 
  for(unsigned int i=0; i<variable_triggers.size; ++i) {
    variable_triggers[i]->notify_add();
  }
} 

void Mistral::Solver::consolidate() 
{

  for(; initialised_vars<variables.size; ++initialised_vars) {
    variables[initialised_vars] = variables[initialised_vars].get_var();
    if(!(domain_types[initialised_vars]&RANGE_VAR) 
       && variables[initialised_vars].domain_type == RANGE_VAR
       && !variables[initialised_vars].is_ground() 
       && variables[initialised_vars].get_degree()>0) {

      int minval = variables[initialised_vars].get_min();
      int maxval = variables[initialised_vars].get_max();
      Variable X(minval, 
		 maxval,
		 domain_types[initialised_vars]);

      X.variable->solver = this;
      X.variable->id = initialised_vars;
      variables[initialised_vars] = X;

      if(X.domain_type > DYN_VAR) {
	booleans.add(variables.stack_+initialised_vars);
      // if(minval==0 && maxval==1) {
      // 	variables+
      // }
      }
    }
  }
//   if(!is_initialised()) {
//     Vector< int > bool_vars;
//     while(initialised_vars < variables.size) {
//       if(variables[initialised_vars].domain_type == BOOL_VAR)
// 	bool_vars.add(initialised_vars);
//       ++initialised_vars;
//     }
//     //booleans.add(bool_vars);
//     //void Mistral::Solver::BooleanMemoryManager::add(Vector< int >& bool_vars) {
//     int *pool = new int[bool_vars.size];
//     std::fill(pool, pool+bool_vars.size, 3);
//     booleans.slots.add(pool);
//     for(unsigned int i=0; i<bool_vars.size; ++i) {
//       variables[bool_vars[i]].bool_domain = pool+i;
//     }
    
//     //     for(int i=0; i<non_convex.size; ++i) {
//     //       non_convex->initialise_bitset_domain();
//     //     }
//     //     boolean_pool.add(new int[booleans.size]);
//     //     for(int i=0; i<booleans.size; ++i) {
//     //       booleans->initialise_boolean_domain(boolean_pool.back()+i);
//     //     }
  while(initialised_cons < constraints.size)
    constraints[initialised_cons++]->consolidate();
//   for(unsigned int i; i<constraints.size; ++i) 
//     constraints[i]->consolidate();
//     {
//     for(unsigned int j=0; j<constraints[i]->scope.size; ++j) {
//       constraints[i]->scope[j] = constraints[i]->scope[j].get_var();
//     }
//}
//   }
}


void Mistral::Solver::make_non_convex(const int idx) 
{

  // if(idx == 426) {
  //    std::cout << std::endl << std::endl << idx << std::endl;
  //    std::cout << "turn " << variables[idx] << " in " ;
  //   std::cout<< variables[idx].get_domain() << " into a bitset variable" << std::endl; 
  // }

  if(variables[idx].domain_type == RANGE_VAR) {
    ConstraintNode nd;

    // std::cout << "SWAP" << std::endl;

    Variable X(variables[idx]// , DYN_VAR&(~RANGE_VAR)
	       , true);
    variables[idx] = X;

    // if(idx == 426) {
    //   std::cout << " => " << X << " in " << X.get_domain() << std::endl; 
    //   std::cout << X.bool_domain << std::endl;
    //   exit(1);
    // }

    
    // std::cout << sequence << std::endl;
    // std::cout << sequence.capacity << std::endl;
    // std::cout << sequence.size << std::endl;
    // std::cout << sequence.offset << std::endl;
    // std::cout << sequence.index_ << std::endl;
    // std::cout << sequence.index_[idx] << std::endl;
    // std::cout << (int*)(sequence.list_) << std::endl;
    

    int ids = sequence.index(idx);

    //std::cout << "ids: " << ids << std::endl;

    if(ids>=0) sequence.list_[ids] = X;

     //std::cout << 11 << std::endl;

    
    nd = constraint_graph[idx]->first(_VALUE_);
    while(constraint_graph[idx]->next(nd)) {
      nd.elt.constraint->consolidate();
    }

    notify_change_variable(idx);

    //std::cout << 22 << std::endl;
  }
}

// void Mistral::Solver::specialise() 
// {
//   ConstraintNode nd;
//   unsigned int i;
//   for(i=0; i<variables.size; ++i) {
//     if((domain_types[i]&RANGE_VAR) && variables[i].domain_type != RANGE_VAR) {

//       //std::cout << "spec " << variables[i] << " => ";

//       Variable R(variables[i].get_min(), variables[i].get_max(), RANGE_VAR);
//       R.variable->solver = this;
//       R.variable->id = i;
//       variables[i] = R;

//       //std::cout << variables[i] << std::endl;

//       nd = constraint_graph[i]->first(_VALUE_);
//       while(constraint_graph[i]->next(nd)) {
	
// 	nd.elt.constraint->consolidate();

// // 	std::cout << "\t" << nd.elt.constraint << " (" 
// // 		  << nd.elt.constraint->scope[nd.elt.index] 
// // 		  << ") => ";

// // 	nd.elt.constraint->scope[nd.elt.index] = variables[i];

// // 	std::cout << "(" 
// // 		  << nd.elt.constraint->scope[nd.elt.index] 
// // 		  << ")" << nd.elt.constraint << std::endl;
//       }
//     }
//   }
// }

bool Mistral::Solver::rewrite() 
{
  
#ifdef _DEBUG_REWRITE
  std::cout << "start rewrite loop: " << (statistics.num_filterings) 
	    << std::endl << active_constraints << std::endl;
#endif

  Constraint *transformed = NULL;

  wiped_idx = -1;
  culprit = NULL;

  ++statistics.num_filterings;
  while( IS_OK(wiped_idx) && !active_constraints.empty() ) {

    // std::cout << "rewrite loop: " << (statistics.num_filterings) 
    //  	      << std::endl << active_constraints << std::endl;

    // std::cout << constraints[35] << std::endl;

    do {
      culprit = active_constraints.select(constraints);
    } while (!culprit->is_posted);

    // std::cout << 11 << std::endl;

    // std::cout << (int*)culprit << std::endl;


    // std::cout << culprit->id << std::endl;

    // std::cout << culprit->scope << std::endl;

    
#ifdef _DEBUG_REWRITE
    int size_before = 0;
    delete o_propag;
    o_propag = new std::ostringstream();
    (*o_propag) << "rewrite " << (culprit) << " b/c" ;
    for(unsigned int i=0; i<culprit->changes.size; ++i) 
      (*o_propag) << " " << culprit->scope[culprit->changes[i]];
    (*o_propag) << std::endl;
    for(unsigned int i=0; i<culprit->scope.size; ++i) {
      size_before += culprit->scope[i].get_size();
      (*o_propag) << culprit->scope[i] << ": " << (culprit->scope[i].get_domain()) << " ";
    }
#endif

    ++statistics.num_propagations;

    
    transformed = NULL;

    //  std::cout << "rewrite **" << culprit->id << "** " << culprit 
    // // 	      // << " " << active_constraints._set_
    //  	      << std::endl;

    wiped_idx = culprit->rewrite();

    //std::cout << "DONE\n\n";

    // std::cout << "rewrite **" << culprit->id << "** " << culprit  
    // 	      << " " << active_constraints._set_ << std::endl;

//     if(IS_OK(wiped_idx) && transformed) 
//       discarded_constraints.add(culprit);
    


//     std::cout << "events of " << culprit << " after defrost: " 
// 	      << culprit->events << std::endl;
//     std::cout << "changes of " << culprits << " after defrost: " 
// 	      << culprit->changes << std::endl;
   
#ifdef _DEBUG_REWRITE
    if(!IS_OK(wiped_idx)) {
      std::cout << (o_propag->str()) << std::endl << culprit->scope[wiped_idx] 
		<< " was wiped out" << std::endl;
    } else {
      for(unsigned int i=0; i<culprit->scope.size; ++i)
	size_before -= culprit->scope[i].get_size();
      if(size_before) {
	std::cout << (o_propag->str()) << std::endl;
	for(unsigned int i=0; i<culprit->scope.size; ++i)
	  std::cout << culprit->scope[i] << ": " 
		    << (culprit->scope[i].get_domain()) << " ";
	std::cout << std::endl << std::endl;
      } 
//       else {
// 	std::cout << (o_propag.str()) << std::endl;
// 	std::cout << "no pruning" << std::endl;
//       }
    }
    if(active_constraints.empty()) 
      std::cout << "Get out of the AC loop because the closure is achieved" << std::endl;
    //else
    //std::cout << active_constraints << std::endl;

#endif 

    //std::cout << "END LOOP" << std::endl;

  }

  //  taboo_constraint = NULL;
  active_constraints.clear();


#ifdef _DEBUG_AC
  if(!IS_OK(wiped_idx)) {
    std::cout << "inconsistency found!" << std::endl;
  } else {
    std::cout << "done" << std::endl;
  }
#endif 

  // std::cout << statistics.num_variables << " x " << statistics.num_constraints << std::endl;

  if(IS_OK(wiped_idx)) {
    statistics.num_constraints = 0;
    for(unsigned int i=0; i<constraints.size; ++i)
      if(constraints[i]->is_posted) ++statistics.num_constraints;

    statistics.num_variables = 0;
    for(unsigned int i=0; i<variables.size; ++i)
      if(variables[i].get_degree() && !variables[i].is_ground()) ++statistics.num_variables;

    // std::cout << statistics.num_variables << " x " << statistics.num_constraints << std::endl;

    return true;
  } else {

    std::cout << "FAIL DURING PREPROCESSING" << std::endl; 

    ++statistics.num_failures;
    notify_failure();
    return false;
  }
  //return !wiped_out;


  return true;

}

Mistral::PropagationOutcome Mistral::Solver::propagate(Constraint *cons) {
  culprit = cons;
  //active_constraints.taboo_constraint = culprit->freeze();
  active_constraints.select(culprit);
  wiped_idx = culprit->propagate();
  culprit->defrost();
  if(IS_OK(wiped_idx)) culprit = NULL;
  return wiped_idx;
}

bool Mistral::Solver::propagate() 
{

#ifdef _DEBUG_QUEUE
  std::cout << "start propagation loop: " << (statistics.num_filterings) 
	    << std::endl << active_constraints << std::endl;
#endif

  wiped_idx = -1;
  culprit = NULL;

  if(objective && objective->enforce())
    wiped_idx = objective->objective.id();

  ++statistics.num_filterings;
  while( IS_OK(wiped_idx) && !active_constraints.empty() ) {
    culprit = active_constraints.select(constraints);
 
#ifdef _DEBUG_PROPAG
    int size_before = 0;
    delete o_propag;
    o_propag = new std::ostringstream();
    (*o_propag) << "propagate " << (culprit) << " b/c" ;
    for(unsigned int i=0; i<culprit->changes.size; ++i) 
      (*o_propag) << " " << culprit->scope[culprit->changes[i]];
    (*o_propag) << std::endl;
    for(unsigned int i=0; i<culprit->scope.size; ++i) {
      size_before += culprit->scope[i].get_size();
      (*o_propag) << culprit->scope[i] << ": " << (culprit->scope[i].get_domain()) << " ";
    }
#endif

    ++statistics.num_propagations;

    wiped_idx = culprit->propagate();

    culprit->defrost();

   
#ifdef _DEBUG_PROPAG
    // bool concerned = false;
    // for(unsigned int zz=0; zz<culprit->scope.size; ++zz) {
    //   if(culprit->scope[zz].id() == 1287) concerned = true;
    // }
    // if(concerned) {
    if(!IS_OK(wiped_idx)) {
      std::cout << (o_propag->str()) << std::endl << culprit->scope[wiped_idx] 
		<< " was wiped out" << std::endl;
    } else {
      for(unsigned int i=0; i<culprit->scope.size; ++i)
	size_before -= culprit->scope[i].get_size();
      if(size_before) {
	std::cout << (o_propag->str()) << std::endl;
	for(unsigned int i=0; i<culprit->scope.size; ++i)
	  std::cout << culprit->scope[i] << ": " 
		    << (culprit->scope[i].get_domain()) << " ";
	std::cout << std::endl << std::endl;
      } 
    }
    if(active_constraints.empty()) 
      std::cout << "Get out of the AC loop because the closure is achieved" << std::endl;
    //}
#endif 


  }

  //  taboo_constraint = NULL;
  active_constraints.clear();


#ifdef _DEBUG_AC
  if(!IS_OK(wiped_idx)) {
    std::cout << "inconsistency found!" << std::endl;
  } else {
    std::cout << "done" << std::endl;
  }
#endif 

  //std::cout << "end propagate: " << culprit << std::endl;

  if(IS_OK(wiped_idx)) {
    notify_success();

#ifdef _DEBUG_PRUNING
    for(unsigned int k=0; k<variables.size; ++k) {
      std::cout << variables[k].get_domain() << " ";
    }
    std::cout << std::endl;
#endif

    return true;
  } else {
    ++statistics.num_failures;
    notify_failure();

#ifdef _DEBUG_PRUNING
    std::cout << "{}" << std::endl;
#endif

    return false;
  }
  //return !wiped_out;
}


void Mistral::Solver::full_print() {
  for(int j=0; j<level; ++j) std::cout << "  ";
  std::cout << sequence << std::endl;
  for(unsigned int i=0; i<variables.size; ++i) {
    for(int j=0; j<level; ++j) std::cout << "  ";
    //variables[i].full_print();
    std::cout << std::endl;
  }
}

void Mistral::Solver::debug_print() {

//   std::cout << std::endl << this << std::endl;

// //   std::cout << "    currently changed: " << std::endl; 
// //   int k = changed_objs.size;
// //   while(k --> 0) {
// //     changed_objs[k]->debug_print();
// //     //std::cout << std::endl;
// //   }
// //   std::cout << std::endl;
//   int i, k;
//   for(i=var_trail_size.size-1; i>0; --i) {
//     std::cout << "    changed at level " << (i-1) << ": " << std::endl; 
//     k = var_trail_size[i];
//     while(k --> var_trail_size[i-1])
//       saved_vars[k]->debug_print();
//     //std::cout << std::endl;
//   }


//   std::cout << "    AC Queue: " << std::endl; 
// //   std::cout << active_constraints.active << std::endl;
// //   std::cout << active_constraints.triggers[0] << std::endl;
// //   std::cout << active_constraints.triggers[1] << std::endl;
// //   std::cout << active_constraints.triggers[2] << std::endl;

//   for(unsigned int i=0; i<3; ++i) {
//     if(!active_constraints.triggers[i].empty()) {
//       std::cout << "priority " << i << std::endl;
//       int elt = active_constraints.triggers[i].first();
//       while(elt != active_constraints.triggers[i]._head) {
// 	Constraint *cons = constraints[elt];
// 	std::cout << "\t" << cons << ": ";
// 	for(unsigned int j=0; j<cons->changes.size; ++j) {
// 	  int var = cons->changes[j];
// 	  int evt = cons->event_type[var];
// 	  std::cout << cons->scope[var] << "/" 
// 		    << (is_value(evt) ? "v" : (is_range(evt) ? "r" : "d") )
// 		    << (is_upper_bound(evt) ? "u" : "")
// 		    << (is_lower_bound(evt) ? "l" : "")
// 		    << " ";
// 	}
// 	std::cout << std::endl;
// 	elt = active_constraints.triggers[i].next[elt];
//       } 
//     }
//   }

// //   for(unsigned int i=0; i<active_constraints.size; ++i) {
// //     Constraint *cons = constraints[active_constraints[i]];
// //     std::cout << cons << ": ";
// //     for(unsigned int j=0; j<cons->changes.size; ++j) {
// //       int var = cons->changes[j];
// //       int evt = cons->event_type[var];
// //       std::cout << cons->scope[var] << "/" 
// // 		<< (is_value(evt) ? "v" : (is_range(evt) ? "r" : "d") )
// // 		<< (is_upper_bound(evt) ? "u" : "")
// // 		<< (is_lower_bound(evt) ? "l" : "")
// // 		<< " ";
// //     }
// //     std::cout << std::endl;
// //   }
}

// std::string Mistral::Solver::getString() const {
//   std::string return_str = "Variables:\n";
//   for(unsigned int i=0; i<variables.size; ++i)
//     return_str += ("\t"+toString(variables[i])+" in "+toString(variables[i]->domain)+"\n");

//   return_str += "\nConstraints:\n";
//   for(unsigned int i=0; i<constraints.size; ++i)
//     return_str += ("\t"+toString(constraints[i])+"\n");

//   return_str += ("\nSearch on "+toString(search.sequence)+"\n");

//   return return_str;
// }

 void Mistral::Solver::print_clist(int k) {
   ConstraintNode nd = constraint_graph[k]->first(_VALUE_);
   if(constraint_graph[k]->next(nd)) {
     std::cout << " [ " << nd.elt.constraint->id;
     while( constraint_graph[k]->next(nd) ) 
       std::cout << " ][ " << nd.elt.constraint->id;
     std::cout << " ] \n";
   }
 }

std::ostream& Mistral::Solver::display(std::ostream& os) {
  os << "Variables:\n";
  ConstraintNode nd;
  //Vector< int > con_id;

  for(unsigned int i=0; i<variables.size; ++i) {
    //if(variables[i].get_degree()>0) {
    os << "  " << variables[i] << " in " << variables[i].get_domain() ; //<< "\n";
    
    //if(!variables[i].is_ground()) {
    nd = constraint_graph[i]->first(_VALUE_);
    if(constraint_graph[i]->next(nd)) {
      os << ": [" << nd.elt.constraint->id;
      while( constraint_graph[i]->next(nd) ) {
	os << "] [" << nd.elt.constraint->id;
      }
      os << "]";
    }
    //}
    os << "\n";
    //}
  }
  
  os << "\nConstraints:\n";
  for(unsigned int i=0; i<constraints.size; ++i)
    if(constraints[i]->is_posted)
      os << "  [" << constraints[i]->id << "]: " << constraints[i] << "\n";
  
  return os;
}


//Mistral::Decision 
void Mistral::Solver::learn_nogood() {

#ifdef _DEBUG_NOGOOD
#ifdef _DEBUG_SEARCH
  for(int i=0; i<level; ++i) std::cout << " "; 
  std::cout << "conflict: " ; 
  //std::cout << (int*)base << std::endl;
  //std::cout << (int*)(base->conflict) << std::endl ;
  print_clause(std::cout, base->conflict);
  std::cout << std::endl;
#endif

  //#ifdef _DEBUG_NOGOOD
  for(unsigned int i=num_search_variables-1; i>=sequence.size; --i) {
    std::cout << num_search_variables-1-i << "\t"
	      << assignment_level[sequence[i].id()] << "\t"
	      << sequence[i] << " == " << sequence[i].get_min() << "\t";
    if(reason[sequence[i].id()]) {
      print_clause(std::cout, reason[sequence[i].id()]);
    } else {
      std::cout << "decision";
    }
    std::cout << std::endl;
  }
#endif

  //if(base->conflict) {

  unsigned int j;
  int pathC = 0, index = sequence.size-1;//, index = assumptions.size;
  Lit p=0, q;
  Atom     a;
  Variable x;
  int lvl;
  Clause *current_clause = base->conflict;
  double *lit_activity = base->lit_activity.stack_;
  double *var_activity = base->var_activity.stack_;


  backtrack_level = 0;
  //std::cout << current_clause->size << std::endl;

  
  learnt_clause.clear();
  learnt_clause.add(p);
  //visited.clear();

  // for(int i=num_search_variables-1; i>=sequence.size; --i) {
  //   std::cout << " " << sequence[i];
  // }

  do {
    // add the parents of the conflict to the current set of visited atoms
    Clause& con = *(current_clause);


    //std::cout << con.size << std::endl;

#ifdef _DEBUG_NOGOOD
    print_clause( std::cout, current_clause );
    std::cout << std::endl;
#endif
    for(j=0; j<con.size; ++j) {
      q = con[j];
      a = UNSIGNED(q);
      x = variables[a];
      lvl = assignment_level[a];
      //lvl = LEVEL(state[a]);
#ifdef _DEBUG_NOGOOD
      std::cout << "\t" ;
      print_literal(std::cout, q); 
      std::cout << ": ";
#endif
      if( !visited.fast_contain(a) ) {
	lit_activity[q] += parameters.activity_increment;
	var_activity[a] += parameters.activity_increment;
	visited.fast_add(a);
	// we'll need to replace 'a' by its parents since its level is too high
	//if(lvl >= decisions.back()) {
	if(lvl >= level) {
#ifdef _DEBUG_NOGOOD
	  std::cout << "expend" << std::endl;
#endif
	  ++pathC;
	} else {
	  // q's level is below the current level, we are not expending it further
	  learnt_clause.add(q);
#ifdef _DEBUG_NOGOOD
	  std::cout << "add to the clause" ;
	  for(unsigned int k=0; k<learnt_clause.size; ++k) {
	    std::cout << " ";
	    print_literal(std::cout, learnt_clause[k]);
	  }
	  std::cout << std::endl;
#endif
	  if(lvl > backtrack_level)
	    backtrack_level = lvl;
	}
      }
#ifdef _DEBUG_NOGOOD
      else {
	std::cout << "visited" << std::endl;
      }
#endif
    }
    // jump to the next visited atom that need be further expended
    //if(index<variables.size-1) {
    while(!visited.fast_contain(sequence[++index].id()))//  {
    //   if(index==variables.size-1) { 
    // 	  pathC = 1;
    // 	  break;
    // 	}
    //   }
    // }
      ;
    x = sequence[index];
    a = x.id();
    p = ((2*a) | (x.get_min()));
    lvl = assignment_level[a];
    //p = (2*a) | SIGN(state[a]); //polarity[a];
    //lvl = LEVEL(state[a]);

#ifdef _DEBUG_NOGOOD
    std::cout << "explore ";
    print_literal(std::cout, p); 
    std::cout << " ";
    //std::cout.flush();
#endif

    if( pathC > 1 ) {
      // there are still atoms to expand, we start with 'a'
      current_clause = reason[a];
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
        
  //exit(1);




  if( learnt_clause.size != 1 ) {
    base->learn(learnt_clause);
    //add_clause( learnt, learnt_clause, stats.learnt_avg_size );
    reason[UNSIGNED(p)] = base->learnt.back();
  }
  visited.clear();


  //backjump_decision = decision(variables[UNSIGNED(p)], Decision::REMOVAL, SIGN(p));

#ifdef _DEBUG_NOGOOD
  std::cout << "backtrackLevel = " << backtrack_level << "/" << (decisions.size) << std::endl;
#endif

//   while(level>backtrack_level) {
//     restore();
//     decisions.pop();
//   }

  //return decision;

}

void Mistral::Solver::forget() {
  if(base) base->forget(parameters.forgetfulness);
}

void Mistral::Solver::branch_right() {
  Mistral::Decision deduction;
  //backtrack_level = level-1;
  // if(parameters.backjump) {
  //   //int backtrack_level=level-1;
  //   decision = learn_nogood_and_backjump(); //backtrack_level);
  //   else {
  //     while(level>backtrack_level) {
  // 	restore();
  // 	decisions.pop();
  //     }
  //   }
  // } else {


  if(parameters.backjump) {
    decisions.size += (backtrack_level-level);
    Lit p = learnt_clause[0];
    deduction = Decision(variables[UNSIGNED(p)], Decision::REMOVAL, NEG(SIGN(p)));
  } else {
    backtrack_level = level-1;
    deduction = decisions.pop();
    deduction.invert();
  }

  restore(backtrack_level);  

  //decision = decisions.pop(); 
  //decision.invert();
  //}

  // if( limits_expired() ) status = LIMITOUT;
  // else {
  
#ifdef _DEBUG_SEARCH
  std::cout << "c";
  for(unsigned int k=0; k<=decisions.size; ++k) std::cout << " ";
  std::cout << deduction << std::endl;
#endif
  
  //decisions.back(-1).make();
  //decision.make();
  deduction.make();
  //}



}


void Mistral::Solver::backjump() {
  int backtrack_level = culprit->get_backtrack_level();
  decisions.size -= (level - backtrack_level);
  restore(backtrack_level);
  Decision decision = culprit->get_decision();
  
#ifdef _DEBUG_SEARCH
  std::cout << "c";
  for(unsigned int k=0; k<=decisions.size; ++k) std::cout << " ";
  std::cout << decision << std::endl;
#endif

  decision.make();
}

void Mistral::Solver::branch_left() {

  save();
  Mistral::Decision decision = heuristic->branch();

  //std::cout << decision << " " << decision.var.get_domain() << std::endl;
  
  //if(decision.var.is_ground()) exit(1);


  // if(decisions.size && (decisions.back() == decision)) {
  //   exit(1);
  // }

  reason[decision.var.id()] = NULL;
  decisions.add(decision);

#ifdef _DEBUG_SEARCH
  std::cout << "c";
  for(unsigned int k=0; k<=decisions.size; ++k) std::cout << " ";
  std::cout << decision << std::endl;
#endif

  //if(decision.var.is_ground()) exit(1);
  
  decision.make();
  notify_decision();


  // //if(level>=2)
  // std::cout << "X4's trail: " << (variables[4].domain_type == RANGE_VAR ?
  // 				  ((VariableRange*)(variables[4].variable))->trail_ :
  // 				  ((VariableBitmap*)(variables[4].variable))->trail_)
  // 	    << " " << variables[4].get_domain() << std::endl;
}


 Mistral::Outcome Mistral::Solver::satisfied() {    
#ifdef _DEBUG_SEARCH
  std::cout << "c";
  for(unsigned int k=0; k<=decisions.size; ++k) std::cout << " ";
  std::cout << " SAT!" << std::endl; 
#endif

    unsigned int i, j, k;

  if(parameters.checked) {
    /// check the current solution
    Vector< int > tmp_sol;
    Constraint *C;
    //for(i=0; i<constraints.size; ++i) {
    for(i=0; i<posted_constraints.size; ++i) {
      C = constraints[posted_constraints[i]];
      
      k=C->scope.size;
      for(j=0; j<k; ++j) {
	if(C->scope[j].is_ground())
	  tmp_sol.add(C->scope[j].get_value());
	else break;
      }

      bool consistent = true;
      if(tmp_sol.size < k) {
	if(parameters.checked==2) {
	  j = C->scope[0].get_min();
	  consistent = (C->first_support(0,j) ||
			C->find_support(0,j));
	} else if(parameters.checked>2) {
	  active_constraints.trigger(C);
	  culprit = active_constraints.select(constraints);
	  if(culprit != C) {
	    std::cout << "c Triggered constraints during checking!" << std::endl;
	  }
	  if(IS_FAIL(C->Constraint::propagate())) {
	    consistent = false;
	  } else {
	    if(!active_constraints.empty()) {
	      consistent = propagate();
	      if(consistent) {
		std::cout << "c Warning, values are not all consistent!" << std::endl;
	      }
	    }
	  }
	}
      } else {
	consistent = !C->check(tmp_sol.stack_);
      }
      if(!consistent)
      {
	
	if(tmp_sol.size < k) {
	  std::cerr << "\nError: solution does not satisfies c" << C->id << ": " << C << tmp_sol << " (backtracking)"<< std::endl;
	  exit(0);
	} else {
	  std::cerr << "\nError: solution does not satisfies c" << C->id << ": " << C ;
	  for(j=0; j<k; ++j) {
	    std::cerr << " " << C->scope[j].get_domain();
	  }
	  std::cerr << " (backtracking)"<< std::endl;
	  exit(0);
	}
	if( decisions.empty() ) return UNSAT;
	else if( limits_expired() ) return LIMITOUT;
	else {
	  branch_right();
	  return UNKNOWN;
	}
      }
      tmp_sol.clear();
    }
  }

  /// store the solution 
  for(i=0; i<variables.size; ++i) {
    last_solution_lb[i] = variables[i].get_min();
    last_solution_ub[i] = variables[i].get_max();

    //std::cout << variables[i] << " := " << last_solution_lb[i] << " ";

  }
  //std::cout << std::endl;
  ++statistics.num_solutions;

  /// notify the objective and return the outcome
  return objective->notify_solution(this);
  
  //return SAT;
}


 Mistral::Outcome Mistral::Solver::exhausted() {    
#ifdef _DEBUG_SEARCH
  std::cout << "c UNSAT!" << std::endl; 
#endif

  return objective->notify_exhausted();
}

 bool Mistral::Solver::limits_expired() {
  
  return (parameters.limit && 
	  ((parameters.time_limit > 0.0 && (get_run_time() - statistics.start_time) > parameters.time_limit) ||
	   (parameters.node_limit > 0 && (statistics.num_nodes > parameters.node_limit)) ||
	   (parameters.fail_limit > 0 && (statistics.num_failures > parameters.fail_limit)) ||
	   (parameters.restart_limit > 0 && (statistics.num_failures > parameters.restart_limit)) ||
	   (parameters.backtrack_limit > 0 && (statistics.num_backtracks > parameters.backtrack_limit))
	   ));
}

// void Mistral::Search::init_search(Vector< Variable >& seq, VarOrdering *h, RestartPolicy *p) {
//   for(unsigned int i=0; i<seq.size; ++i) {
//     std::cout << "insert " << seq[i] << std::endl;
//     sequence.insert(seq[i]);    
//   }
//   if(heuristic) delete heuristic;
//   if(policy) delete policy;
//   heuristic = (h ? h : new NoOrder(sequence));
//   policy = p;
// }


// std::string Mistral::toString(const Mistral::Solver& x) {
//   return x.getString();
// }

Mistral::Outcome Mistral::Solver::chronological_dfs() 
{
  int status = UNKNOWN;
  while(status == UNKNOWN) {
    if(propagate()) {
      if( sequence.empty()  ) status = satisfied();
      else branch_left();
    } else {
      if( parameters.backjump ) learn_nogood();
      if( decisions.empty() ) status = exhausted();
      else if( limits_expired() ) status = LIMITOUT;
      else branch_right();
    }
  }
  return status;
}


Mistral::Outcome Mistral::Solver::conflict_directed_backjump()
{
  int status = UNKNOWN;
  while(status == UNKNOWN) {
    if(propagate()) {
      if( sequence.empty() ) status = satisfied();
      else branch_left();
    } else {
      if( decisions.empty() ) status = UNSAT;
      else if( limits_expired() ) status = LIMITOUT;
      else backjump();
    }
  }

  return status;
}

std::ostream& Mistral::operator<< (std::ostream& os, Mistral::Solution& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, Mistral::Solution* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, Mistral::Solver& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, Mistral::Solver* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, Mistral::ConstraintQueue& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, Mistral::ConstraintQueue* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::SolverStatistics& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::SolverStatistics* x) {
  return x->display(os);
}

