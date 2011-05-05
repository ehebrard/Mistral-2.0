
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
  std::cout << active_solver->statistics //<< std::endl
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
    values[variables[i].id()] = variables[i].get_solution_value();
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
  restart_limit = 0;
  limit = 0;

  verbosity = 0;
  time_limit = -1;
  seed = 11041979;
  restart_policy = GEOMETRIC;
  restart_base = 200;
  restart_limit = 200;
  restart_factor = 1.05;
  activity_increment = 1e-2;
  normalize_activity = 0;
  init_activity = 1;
  forgetfulness = .75;
  randomization = 2;
  shuffle = true;
  activity_decay = 0.96;
  checked = true;
  value_selection = 6;
  dynamic_value = 1;
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
  os << " c +=============================================================================+" << std::endl
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
     << " c +=============================================================================+" << std::endl;
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
     << (getRunTime() - start_time) << " |" ;
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

//   std::cout << "declare " << c << " to the constraint queue" << std::endl;
//   std::cout << "was: [" << min_priority << ","
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


//     std::cout << "==> " << triggers[cons_priority] << std::endl;


  } else {
    
//     std::cout << "no need to create a new trigger list" << std::endl;
//     std::cout << "extend " << triggers[cons_priority] << " with " << cons_idx << std::endl;

    triggers[cons_priority].extend(cons_idx);
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
  
  cons->events.clear();
  for(unsigned int i=0; i<cons->scope.size; ++i) {
    evt = (// cons->scope[i].domain_type != BOOL_VAR &&
	   cons->scope[i].is_ground() ? VALUE_EVENT : (LB_EVENT|UB_EVENT));
    if(EVENT_TYPE(evt) <= cons->trigger[i]) {
      cons->events.add(i);
      cons->event_type[i] = evt;
      triggered=true;
    }
  }
  if(triggered) {
    _set_.fastAdd(cons_id);
    if(priority > higher_priority) higher_priority = priority;
    triggers[priority].add(cons_id);
  }  
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

      _set_.fastAdd(cons_id);
      if(priority > higher_priority) higher_priority = priority;
      triggers[priority].add(cons_id);
      cons->events.setTo(var);
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
	 << solver->constraints[elt];
      while(triggers[i+min_priority].next[elt] != end) {
	elt = triggers[i+min_priority].next[elt];
	os << ", " << solver->constraints[elt];
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

  // variables & constraints
  domain_types.initialise(0,128);
  variables.initialise(0,128);
  constraints.initialise(0,256);
  //constraints.initialise(256);
  posted_constraints.initialise(0,256,false);
  sequence.initialise(128);
  initialised_vars = 0;
  initialised_cons = 0;

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
  ConstraintClauseBase *base = new ConstraintClauseBase(variables);
  add(base);

  for(int i=0; i<M; ++i)
    {
      new_clause.clear();
      do {
	infile >> l;
	if(l!=0) {
	  if(l>0) lit = (l-1)*2+1;
	  else lit = (l+1)*-2;
	  new_clause.add(lit);
	  //if(params.init_activity == 1)
	  //activity[lit] += params.activity_increment;
	}
      } while(l && infile.good());
      base->add( new_clause );

      //std::cout << new_clause << std::endl;
      
      //if(params.checked) add_original_clause( new_clause );
    }
  
  //init_watchers();
  
  //   if(params.normalize_activity != 0)
  //     normalize_activity(params.normalize_activity);
  
  //  std::cout << base << std::endl;
}

void Mistral::Solver::set_parameters(SolverParameters& p) {
  parameters = p;
}

void Mistral::Solver::add(VarArray& x) {
  for(unsigned int i=0; i<x.size; ++i)
    x[i].initialise(this);
}

void Mistral::Solver::add(Variable x) { 
  x.initialise(this); 
}

int Mistral::Solver::declare(Variable x) {
  if(x.domain_type > DYN_VAR) booleans.add(&x);

  // add the variables to the set of vars
  x.variable->id = variables.size;
  x.variable->solver = this;
  variables.add(x);
  domain_types.add(DYN_VAR|(x.is_range() ? RANGE_VAR : 0));

  last_solution_lb.add(-INFTY);
  last_solution_ub.add( INFTY);
  
  // create a new empty list of constraints
  ConstraintList *l = new ConstraintList(this);
  constraint_graph.add(l);

  return variables.size-1;
}


void Mistral::Solver::add(Constraint* c) { 
  if(c->id < 0) {
    if(parameters.verbosity>0) {
      std::cout << "c add a constraint: " << c << std::endl; 
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

  if(con_trail_.back() != level) {
    con_trail_.add(posted_constraints.size);
    con_trail_.add(level);
  }

  if(!posted_constraints.contain(c->id)) {
    posted_constraints.extend(c->id);
    posted_constraints.add(c->id);
  }
}

Mistral::Outcome Mistral::Solver::solve() {
  BranchingHeuristic *heu = NULL; 
  RestartPolicy *pol = NULL;
  return depth_first_search(variables, heu, pol);
}

Mistral::Outcome Mistral::Solver::depth_first_search(BranchingHeuristic *heu, 
						     RestartPolicy *pol) {
  return depth_first_search(variables, heu, pol);
}

Mistral::Outcome Mistral::Solver::depth_first_search(Vector< Variable >& seq, 
						     BranchingHeuristic *heu, 
						     RestartPolicy *pol) 
{
  initialise_search(seq, heu, pol);

  Outcome satisfiability = UNKNOWN;
  while(satisfiability == UNKNOWN) {

    statistics.num_variables = sequence.size;
    statistics.num_values = 0;
    for(unsigned int i=0; i<sequence.size; ++i)
      statistics.num_values += sequence[i].get_size();

    if(parameters.verbosity) {
      statistics.print_short(std::cout);
      std::cout << std::endl;
    }

    ++statistics.num_restarts;
    satisfiability = iterative_dfs();

    if(satisfiability == LIMITOUT) {
      policy->reset(parameters.fail_limit);    
      if(!limits_expired()) satisfiability = UNKNOWN;
    }

    restore(0);
    decisions.clear();
  }

  statistics.outcome = satisfiability;
  statistics.end_time = getRunTime();
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
  satisfiability = iterative_dfs();
  
  if(parameters.verbosity) {
    statistics.print_short(std::cout);
  }

  return satisfiability;
}

void Mistral::Solver::BooleanMemoryManager::add(Variable *x) {
  if(size.back() < 1000) {
    x->bool_domain = slots.back()+size.back();
    ++size.back();
  } else {
    int *nslot = new int[1000];
    std::fill(nslot, nslot+1000, 3);
    size.add(0);
    slots.add(nslot);
    x->bool_domain = nslot;
  }
}


void Mistral::Solver::initialise_search(Vector< Variable >& seq, 
					BranchingHeuristic *heu, 
					RestartPolicy *pol) 
{
  consolidate();

  if(level < 0) save();

  active_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  if(statistics.start_time == 0.0) statistics.start_time = getRunTime();

  sequence.clear();
  decisions.clear();
  for(unsigned int i=seq.size; i;) {
    Variable x = seq[--i].get_var();
    if(!x.is_ground() && !sequence.contain(x)) sequence.add(x);
  }

  std::cout <<  sequence << std::endl;

  if(heu) { delete heuristic; heuristic = heu; }
  else if(!heuristic) heuristic = new GenericHeuristic< NoOrder, MinValue >(this);
  if(pol) { delete policy;    policy    = pol; }
  else if(!policy) policy = new NoRestart();

  heuristic->initialise(sequence);

  parameters.fail_limit = policy->base;
  parameters.limit = (policy->base > 0);

  statistics.num_constraints = constraints.size;

  if(parameters.verbosity)  std::cout << " c +=============================================================================+" << std::endl 
				  << " c |      INSTANCE STATS       |                    SEARCH STATS                 |" << std::endl 
				  << " c |   vars |    vals |   cons |    nodes | filterings | propagations | cpu time |" << std::endl;
}



Mistral::Solver::~Solver() {
  delete heuristic;
  delete policy;

  for(unsigned int i=0; i<variables.size; ++i) {
    int domain_type = variables[i].domain_type;
    if     (domain_type ==  BITSET_VAR) delete variables[i].bitset_domain;
    else if(domain_type ==    LIST_VAR) delete variables[i].list_domain;
    else if(domain_type ==   RANGE_VAR) delete variables[i].range_domain;
    else if(domain_type == VIRTUAL_VAR) delete variables[i].virtual_domain;
    else if(domain_type ==  EXPRESSION) delete variables[i].expression;
    delete constraint_graph[i];
  }
}


void Mistral::Solver::trigger_event(const int var, 
				    const Mistral::Event evt) {


#ifdef _DEBUG_QUEUE
  std::cout << (ASSIGNED(evt) ? "value" : (BOUND_CHANGED(evt) ? "range" : "domain"))
	    << " event on " << variables[var] 
	    << " in " << variables[var].get_domain() << std::endl;
#endif

  Constraint *c;
  ConstraintNode nd;

  nd = constraint_graph[var]->first(EVENT_TYPE(evt));
  if(ASSIGNED(evt)) {
    while(constraint_graph[var]->next(nd)) {
      active_constraints.trigger(nd.elt.constraint, nd.elt.index, evt);
      c = nd.elt.constraint->notify_assignment(nd.elt.index, level);
      if(c) saved_cons.add(c);
    }
    if(sequence.contain(variables[var])) sequence.remove(variables[var]);
  } else while(constraint_graph[var]->next(nd)) {
      active_constraints.trigger(nd.elt.constraint, nd.elt.index, evt);
    }

  wiper_idx = var;
  
}


// void Mistral::Solver::trigger(Constraint* cons,
// 			      const int var, 
// 			      const Mistral::Event evt) {

//  #ifdef _DEBUG_AC
//   if(statistics.num_filterings > 10748) {
//     std::cout << "trigger " << cons << "(" << (cons->id) << ") by a " << (is_value(evt) ? "value" : (is_range(evt) ? "range" : "domain")) << " event on " << cons->scope[var] << std::endl;
//   }
// #endif

// //   if(is_value(evt)) {
// //     //assign(cons->_scope[var]);
// //     cons->assign(var);
// //   }

//   if(cons != taboo_constraint) 
//     active_constraints.trigger(cons, var, evt);

// #ifdef _DEBUG_AC
//   if(statistics.num_filterings > 10748) {
//     std::cout << "after trigger: " << active_constraints.active ;
//     for(unsigned int i=0; i<active_constraints.active.size; ++i)
//       std::cout << " [" << active_constraints.triggers[active_constraints.active[i]] << "]";
//     std::cout << std::endl;  
//   }
// #endif
// }

// void Mistral::Solver::assign(Variable x) { 
 
//   std::cout << "assign " << x << " | " << sequence << std::endl;

//    //if(sequence.member(x)) 
//    assert(sequence.member(x));

//    sequence.erase(x); 
// }

// void Mistral::Solver::make_node() {
//   obj_trail_size.add(saved_objs.size);
//   var_trail_size.add(saved_vars.size);
//   con_trail_size.add(saved_cons.size);

//   trail_seq.add(sequence.size);
//   trail_aux.add(auxilliary.size);
//   //decision.add(heuristic->select());

//   ++statistics.num_nodes;
//   ++level;
// }

void Mistral::Solver::save(Variable x) { 
  saved_vars.add(x); 
}

void Mistral::Solver::save(VariableImplementation *impl, int dtype) { 
  Variable x(impl, dtype);
  saved_vars.add(x); 
}

void Mistral::Solver::save(const int idx) { 
  saved_vars.add(variables[idx]);
}

void Mistral::Solver::restore() {
  unsigned int previous_level;

  previous_level = trail_.pop();
  while( saved_vars.size > previous_level ) 
    saved_vars.pop().restore();

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

//   std::cout << "X5CLIST AFTER RESTORE ";
//     //<< constraint_graph[5] << std::endl;
//   print_clist(5);
//   std::cout << "C10ALIST AFTER RESTORE "
// 	    << constraints[10]->active << " " << constraints[10]->trail_ << std::endl;
  

}

void Mistral::Solver::restore(const int lvl) {
  //std::cout << lvl << " " << level << std::endl;
  while(lvl < level) restore();
}

void Mistral::Solver::add(Mistral::RestartListener* l) {
  l->rid = restart_triggers.size;
  restart_triggers.add(l);
}
void Mistral::Solver::add(Mistral::DecisionListener* l) {
  l->did = decision_triggers.size;
  decision_triggers.add(l);
}
void Mistral::Solver::add(Mistral::SuccessListener* l) {
  l->sid = success_triggers.size;
  success_triggers.add(l);
}
void Mistral::Solver::add(Mistral::FailureListener* l) {
  l->fid = failure_triggers.size;
  failure_triggers.add(l);
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

void Mistral::Solver::consolidate() 
{

  for(unsigned int i=0; i<variables.size; ++i)
    variables[i] = variables[i].get_var();

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

void Mistral::Solver::specialise() 
{
  ConstraintNode nd;
  unsigned int i;
  for(i=0; i<variables.size; ++i) {
    if((domain_types[i]&RANGE_VAR) && variables[i].domain_type != RANGE_VAR) {

      //std::cout << "spec " << variables[i] << " => ";

      Variable R(variables[i].get_min(), variables[i].get_max(), RANGE_VAR);
      R.variable->solver = this;
      R.variable->id = i;
      variables[i] = R;

      //std::cout << variables[i] << std::endl;

      nd = constraint_graph[i]->first(_VALUE_);
      while(constraint_graph[i]->next(nd)) {
	
	nd.elt.constraint->consolidate();

// 	std::cout << "\t" << nd.elt.constraint << " (" 
// 		  << nd.elt.constraint->scope[nd.elt.index] 
// 		  << ") => ";

// 	nd.elt.constraint->scope[nd.elt.index] = variables[i];

// 	std::cout << "(" 
// 		  << nd.elt.constraint->scope[nd.elt.index] 
// 		  << ")" << nd.elt.constraint << std::endl;
      }
    }
  }
}

bool Mistral::Solver::rewrite() 
{
  
#ifdef _DEBUG_REWRITE
  std::cout << "start rewrite loop: " << (statistics.num_filterings) 
	    << std::endl << active_constraints << std::endl;
#endif

  int wiped_idx = -1;
  Constraint *culprit = NULL, *transformed = NULL;

  ++statistics.num_filterings;
  while( IS_OK(wiped_idx) && !active_constraints.empty() ) {
   
    culprit = active_constraints.select(constraints);
 
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

    //std::cout << "rewrite " << culprit << std::endl;
    wiped_idx = culprit->rewrite();
//     if(IS_OK(wiped_idx) && transformed) 
//       discarded_constraints.add(culprit);
    


//     std::cout << "events of " << culprit << " after defrost: " 
// 	      << culprit->events << std::endl;
//     std::cout << "changes of " << culprit << " after defrost: " 
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

  if(IS_OK(wiped_idx)) {
    return true;
  } else {
    ++statistics.num_failures;
    notify_failure();
    return false;
  }
  //return !wiped_out;

  return true;

}

bool Mistral::Solver::propagate() 
{


#ifdef _DEBUG_QUEUE
  std::cout << "start propagation loop: " << (statistics.num_filterings) 
	    << std::endl << active_constraints << std::endl;
// active_constraints.active << " ";
//   for(unsigned int i=0; i<active_constraints.active.size; ++i) {
//     std::cout << active_constraints.triggers[active_constraints.active[i]] << " ";
//   }
//   std::cout << std::endl;
#endif

  //Variable wiped_out = NULL;
  //Constraint *cons;
  //int 
  int wiped_idx = -1;
  Constraint *culprit = NULL;

  ++statistics.num_filterings;
  while( IS_OK(wiped_idx) && !active_constraints.empty() ) {
//     std::cout << std::endl << "propagation step: " 
// 	      << active_constraints ; //<< std::endl;
   
    culprit = active_constraints.select(constraints);


    //display(std::cout);
 
#ifdef _DEBUG_PROPAG
    int size_before = 0;
    //std::ostringstream o_propag;
    delete o_propag;
    o_propag = new std::ostringstream();

 //    o_propag << std::endl << "Propagation step" << std::endl;
//     for(unsigned int i=0; i<variables.size; ++i) {
//       o_propag << variables[i] << " in " << variables[i].get_domain() << ": ";
//       ConstraintNode nd;
//       nd = constraint_graph[i]->first(_VALUE_);
//       o_propag << "[" ;
//       while( constraint_graph[i]->next(nd) ) {
// 	for(unsigned int j=0; j<nd.elt.constraint->scope.size; ++j)
// 	  if(nd.elt.constraint->scope[j].id() == variables[i].id())
// 	    o_propag << nd.elt.constraint->scope[j].get_domain();
//       }
//       o_propag << "]" << std::endl;
//     }
//     o_propag << "\n";

    

    (*o_propag) << "propagate " << (culprit) << " b/c" ;
    for(unsigned int i=0; i<culprit->changes.size; ++i) 
      (*o_propag) << " " << culprit->scope[culprit->changes[i]];
    (*o_propag) << std::endl;
    for(unsigned int i=0; i<culprit->scope.size; ++i) {
      size_before += culprit->scope[i].get_size();
      (*o_propag) << culprit->scope[i] << ": " << (culprit->scope[i].get_domain()) << " ";
    }
    //(*o_propag) << std::endl;
#endif

    ++statistics.num_propagations;


    //culprit->freeze();

//     std::cout << "events of " << culprit << " after freeze: " 
// 	      << culprit->events << std::endl;
//     std::cout << "changes of " << culprit << " after freeze: " 
// 	      << culprit->changes << std::endl;
    //std::cout << "propagate " << (culprit) << std::endl;
    wiped_idx = culprit->propagate();

    culprit->defrost();

//     std::cout << "events of " << culprit << " after defrost: " 
// 	      << culprit->events << std::endl;
//     std::cout << "changes of " << culprit << " after defrost: " 
// 	      << culprit->changes << std::endl;
   
#ifdef _DEBUG_PROPAG
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

  if(IS_OK(wiped_idx)) {
    return true;
  } else {
    ++statistics.num_failures;
    notify_failure();
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
  for(unsigned int i=0; i<variables.size; ++i) {
    os << "\t" << variables[i] << " in " << variables[i].get_domain() ; //<< "\n";

    //os << domain_types[i] << "&" << RANGE_VAR << " ";
    if(domain_types[i] & RANGE_VAR)
      os << "(r/" << constraint_graph[i]->degree << ")";
    else
      os << "(d/" << constraint_graph[i]->degree << ")";

    if(!variables[i].is_ground()) {
      nd = constraint_graph[i]->first(_VALUE_);
      if(constraint_graph[i]->next(nd)) {
	os << ": [ " << nd.elt.constraint;

// 	if(!constraint_graph[i]->contain(nd.elt.constraint->self[nd.elt.index])) {
// 	  std::cout << std::endl << " is not in " << std::endl;
// 	  exit(1);
// 	}
	
	while( constraint_graph[i]->next(nd) ) {
	  os << " ][ " << nd.elt.constraint;

// 	  if(!constraint_graph[i]->contain(nd.elt.constraint->self[nd.elt.index])) {
// 	    std::cout << std::endl << " is not in " << std::endl;
// 	    exit(1);
// 	  }
	}
	os << " ]";
      }
    }
    os << "\n";
  }

//   os << "\nConstraints:\n";
//   for(unsigned int i=0; i<constraints.size; ++i)
//     os << "\t" << constraints[i] << "\n";

  //return_str += ("\nSearch on "+toString(search.sequence)+"\n");
  return os;
}

void Mistral::Solver::branch_right() {
  restore();  
  Mistral::Decision decision = decisions.pop(); //search.backtrack();
  decision.invert();
  //decisions.back(-1).invert();

#ifdef _DEBUG_SEARCH
  std::cout << "c";
  for(unsigned int k=0; k<=decisions.size; ++k) std::cout << " ";
  std::cout << decision << std::endl;
#endif

  //decisions.back(-1).make();
  decision.make();
}

void Mistral::Solver::branch_left() {
  save();
  Mistral::Decision decision = heuristic->branch(); 
  decisions.add(decision);

#ifdef _DEBUG_SEARCH
  std::cout << "c";
  for(unsigned int k=0; k<=decisions.size; ++k) std::cout << " ";
  std::cout << decision << std::endl;
#endif

  if(decision.var.is_ground()) exit(1);
  
  decision.make();
}


 Mistral::Outcome Mistral::Solver::satisfied() {    
#ifdef _DEBUG_SEARCH
  std::cout << "c";
  for(unsigned int k=0; k<=decisions.size; ++k) std::cout << " ";
  std::cout << " SAT!" << std::endl; 
#endif


  Vector< int > tmp_sol;
  unsigned int i, j, k;
  Constraint *C;
  //for(i=0; i<constraints.size; ++i) {
  for(i=0; i<posted_constraints.size; ++i) {
    C = constraints[posted_constraints[i]];
    
    k=C->scope.size;
    for(j=0; j<k; ++j) 
      tmp_sol.add(C->scope[j].get_value());


    if(C->check(tmp_sol.stack_))
      {
	std::cerr << "\nError: solution does not satisfies " << C << tmp_sol << " (backtracking)"<< std::endl;
	exit(0);
	if( decisions.empty() ) return UNSAT;
	else if( limits_expired() ) return LIMITOUT;
	else {
	  branch_right();
	  return UNKNOWN;
	}
      }
    tmp_sol.clear();
  }
  
  for(i=0; i<variables.size; ++i) {
    last_solution_lb[i] = variables[i].get_min();
    last_solution_ub[i] = variables[i].get_min();
  }
  ++statistics.num_solutions;
  
  return SAT;
}

 bool Mistral::Solver::limits_expired() {
  
  return (parameters.limit && 
	  ((parameters.time_limit > 0.0 && (getRunTime() - statistics.start_time) > parameters.time_limit) ||
	   (parameters.node_limit > 0 && (statistics.num_nodes > parameters.node_limit)) ||
	   (parameters.fail_limit > 0 && (statistics.num_failures > parameters.fail_limit)) ||
	   (parameters.restart_limit > 0 && (statistics.num_restarts > parameters.restart_limit)) ||
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

Mistral::Outcome Mistral::Solver::iterative_dfs() 
{
  int status = UNKNOWN;
  while(status == UNKNOWN) {
    if(propagate()) {

      //std::cout << "new node " << sequence << std::endl;
      //        display(std::cout);

      if( sequence.empty() ) status = satisfied();
      else branch_left();
    } else {
      if( decisions.empty() ) status = UNSAT;
      else if( limits_expired() ) status = LIMITOUT;
      else branch_right();

//        std::cout << "restore:\n";
//        display(std::cout);

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

