
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


#include <mistral_solver.hpp>
#include <mistral_variable.hpp>
#include <mistral_constraint.hpp>
#include <sstream>




Mistral::SolverParameters::SolverParameters() {}
Mistral::SolverParameters::~SolverParameters() {}
void Mistral::SolverParameters::initialise() {
  verbosity = 0;
  find_all = 0;
  node_limit = 0;
  backtrack_limit = 0;
  fail_limit = 0;
  restart_limit = 0;
  limit = 0;
  time_limit = 0.0;
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

Mistral::SolverStatistics::SolverStatistics() {}
Mistral::SolverStatistics::~SolverStatistics() {}
void Mistral::SolverStatistics::initialise() {
  num_nodes = 0; 
  num_restarts = 0; 
  num_backtracks = 0;
  num_failures = 0; 
  num_propagations = 0;
  num_solutions = 0;
  num_filterings = 0;
  start_time = 0.0;
  end_time = -1.0;
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
std::ostream& Mistral::SolverStatistics::display(std::ostream& os) const {
  if(end_time >= 0.0) print_full(os);
  else {
    os << " c |"
       << std::right << std::setw(7) << num_variables << " |"
       << std::right << std::setw(8) << num_values << " |"
       << std::right << std::setw(7) << num_constraints << " |"
       << std::right << std::setw(9) << num_nodes << " |" 
       << std::right << std::setw(11) << num_filterings << " |" 
       << std::right << std::setw(13) << num_propagations << " |" 
       << std::right << std::setw(9) << std::setprecision(5) 
       << (getRunTime() - start_time) << " |" ;
  }
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


Mistral::Solver::Solver() { 
  // search stuf
  heuristic = NULL;
  policy = NULL;

  // variables & constraints
  variables.initialise(0,128);
  constraints.initialise(0,256);

  // trail stuff
  level = 0;
  saved_objects.initialise(0,4096); 
  saved_vars.initialise(0,4096); 
  obj_trail_size.initialise(0,512);
  obj_trail_size.add(0);
  var_trail_size.initialise(0,512);
  var_trail_size.add(0);
  con_trail_size.initialise(0,512);
  con_trail_size.add(0);

  // params
  parameters.initialise();

  // statistics
  statistics.initialise();

  // other stuff
  taboo_constraint = NULL;
}

void Mistral::Solver::initialise() {
  //Mistral::Constraint *cons;
  Mistral::IntVar var;
  unsigned int i;

  //search.initialise(this);
  sequence.initialise(variables, false);
  auxilliary.initialise(variables, false);
  trail_seq.initialise(0, variables.size);
  trail_aux.initialise(0, variables.size);
  decision.initialise(0, variables.size);

  triggered_constraints.initialise(constraints.size);

  // trigger evrything in order to achieve a full propagation round 

  for(i=0; i<variables.size; ++i) {
    var = variables[i];
    if(var->isGround()) var->triggerEvent(VALUE_EVENT);
    else                var->triggerEvent(RANGE_EVENT);
  }    

//   for(i=0; i<constraints.size; ++i) {
//     cons = constraints[i];

//     triggered_constraints.trigger(cons, 0, (cons->scope[0]->isGround() ? VALUE_EVENT : RANGE_EVENT));
//     cons->events.fill();
//     j = cons->arity;
//     while(j --> 1) cons->event_type[j] = (cons->scope[j]->isGround() ? VALUE_EVENT : RANGE_EVENT);
//   }

  heuristic = new NoOrder(this);
  policy = new NoRestart();
}

void Mistral::Solver::add(Mistral::IntVar x) { 
  if(x->id < 0) {
    x->id = variables.size; 
    x->solver = this;
    x->constraints.solver = this;
    variables.add(x);
    solution.add(x->domain.min);
  }
}

void Mistral::Solver::add(Mistral::Constraint* c) { 
  if(c->id < 0) {
    for(int i=0; i<c->arity; ++i) add(c->scope[i]);
    c->id = constraints.size; 
    constraints.add(c); 
  }
}



Mistral::Outcome Mistral::Solver::depth_first_search(VarOrdering *heu, 
						     RestartPolicy *pol) {
  return depth_first_search(variables, heu, pol);
}


Mistral::Outcome Mistral::Solver::depth_first_search(Vector< IntVar >& seq, 
						     VarOrdering *heu, 
						     RestartPolicy *pol) 
{

  if(statistics.start_time == 0.0) statistics.start_time = getRunTime();

  auxilliary.fill();
  for(unsigned int i=variables.size; i;)
    variables[--i]->var_list = &auxilliary;
  sequence.clear();
  for(unsigned int i=seq.size; i;) {
    sequence.insert(seq[--i]);
    sequence.back()->var_list = &sequence;
  }

  if(heuristic) delete heuristic;
  if(policy) delete policy;

  heuristic = (heu ? heu : new NoOrder(this));
  policy = (pol ? pol : new NoRestart());

  parameters.fail_limit = policy->base;
  parameters.limit = (policy->base > 0);

  statistics.num_constraints = constraints.size;

  Outcome s = UNKNOWN;

  std::cout << " c +=============================================================================+" << std::endl 
	    << " c |      INSTANCE STATS       |                    SEARCH STATS                 |" << std::endl 
	    << " c |   vars |    vals |   cons |    nodes | filterings | propagations | cpu time |" << std::endl;

  while(s == UNKNOWN) {

    statistics.num_variables = sequence.size;
    statistics.num_values = 0;
    for(unsigned int i=0; i<sequence.size; ++i)
      statistics.num_values += sequence[i]->domain.size;

    std::cout << statistics << std::endl;

    s = iterative_dfs();

    if(s == LIMITOUT) {
      policy->reset(parameters.fail_limit);    
      if(!limitsExpired()) s = UNKNOWN;
    }

    backtrack(0);
  }

  statistics.outcome = s;
  statistics.end_time = getRunTime();
  return s;
}

Mistral::Solver::~Solver() {
  delete heuristic;
  delete policy;
}

void Mistral::Solver::trigger(Constraint* cons,
			      const int var, 
			      const Mistral::Event evt) {

 #ifdef _DEBUG_AC
  if(statistics.num_filterings > 10748) {
    std::cout << "trigger " << cons << "(" << (cons->id) << ") by a " << (is_value(evt) ? "value" : (is_range(evt) ? "range" : "domain")) << " event on " << cons->scope[var] << std::endl;
  }
#endif

//   if(is_value(evt)) {
//     //assign(cons->_scope[var]);
//     cons->assign(var);
//   }

  if(cons != taboo_constraint) 
    triggered_constraints.trigger(cons, var, evt);

#ifdef _DEBUG_AC
  if(statistics.num_filterings > 10748) {
    std::cout << "after trigger: " << triggered_constraints.active ;
    for(unsigned int i=0; i<triggered_constraints.active.size; ++i)
      std::cout << " [" << triggered_constraints.triggers[triggered_constraints.active[i]] << "]";
    std::cout << std::endl;  
  }
#endif
}

// void Mistral::Solver::assign(IntVar x) { 
 
//   std::cout << "assign " << x << " | " << sequence << std::endl;

//    //if(sequence.member(x)) 
//    assert(sequence.member(x));

//    sequence.erase(x); 
// }

void Mistral::Solver::make_node() {
  obj_trail_size.add(saved_objects.size);
  var_trail_size.add(saved_vars.size);
  con_trail_size.add(saved_cons.size);

  trail_seq.add(sequence.size);
  trail_aux.add(auxilliary.size);
  decision.add(heuristic->select());

  ++statistics.num_nodes;
  ++level;
}


void Mistral::Solver::backtrack() {
  unsigned int previous_level;
  IntVar var;
  Reversible *object;
  Constraint *con;
  
  previous_level = var_trail_size.pop();
  while( saved_vars.size > previous_level ) {
    saved_vars.pop(var);
    var->restore();
  }

  previous_level = con_trail_size.pop();
  while( saved_cons.size > previous_level ) {
    saved_cons.pop(con);
    con->restore();
  }

  previous_level = obj_trail_size.pop();
  while( saved_objects.size > previous_level ) {
    saved_objects.pop(object);
    object->restore();
  }

  previous_level = sequence.size;
  trail_seq.pop(sequence.size);
  while(previous_level < sequence.size) {
    //std::cout << "unassign " << (sequence[previous_level]) << std::endl;
    sequence[previous_level++]->unassign();
  }
  //return decision.pop();

  previous_level = auxilliary.size;
  trail_aux.pop(auxilliary.size);

  while(previous_level < auxilliary.size) {
    //std::cout << "unassign " << (auxilliary[previous_level]) << std::endl;
    auxilliary[previous_level++]->unassign();
  }
  //return decision.pop();

  ++statistics.num_backtracks;
  --level;
}

void Mistral::Solver::backtrack(const int& lvl) {
  while(lvl < level) backtrack();
}

bool Mistral::Solver::propagate() 
{

#ifdef _DEBUG_AC
  std::cout << "start propagation loop: " << (statistics.num_filterings) << " " << triggered_constraints.active << " ";
  for(unsigned int i=0; i<triggered_constraints.active.size; ++i) {
    std::cout << triggered_constraints.triggers[triggered_constraints.active[i]] << " ";
  }
  std::cout << std::endl;
#endif

  IntVar wiped_out = NULL;
  Constraint *cons;
  

  ++statistics.num_filterings;
  while( !wiped_out && !triggered_constraints.empty() ) {

    cons = constraints[triggered_constraints.pop()];
    taboo_constraint = cons->freeze();
    
#ifdef _DEBUG_PROPAG
    int size_before = 0;
    std::ostringstream o_propag;

    o_propag << "propagate " << (cons) << " b/c" ;
    for(unsigned int i=0; i<cons->changes.size; ++i) 
      o_propag << " " << cons->scope[cons->changes[i]];
    o_propag << std::endl;
    for(int i=0; i<cons->arity; ++i) {
      size_before += cons->scope[i]->domain.size;
      o_propag << cons->scope[i] << ": " << (cons->scope[i]->domain) << " ";
    }
    //o_propag << std::endl;
#endif

    ++statistics.num_propagations;
    wiped_out = cons->propagate();
    cons->defrost();

#ifdef _DEBUG_PROPAG
    if(wiped_out) {
      std::cout << (o_propag.str()) << std::endl << wiped_out 
		<< " was wiped out" << std::endl;
    } else {
      for(int i=0; i<cons->arity; ++i)
	size_before -= cons->scope[i]->domain.size;
      if(size_before) {
	std::cout << (o_propag.str()) << std::endl;
	for(int i=0; i<cons->arity; ++i)
	  std::cout << cons->scope[i] << ": " 
		    << (cons->scope[i]->domain) << " ";
	std::cout << std::endl << std::endl;
      } 
//       else {
// 	std::cout << (o_propag.str()) << std::endl;
// 	std::cout << "no pruning" << std::endl;
//       }
    }
#endif 

  }

  taboo_constraint = NULL;
  triggered_constraints.clear();


#ifdef _DEBUG_AC
  if(wiped_out) {
    std::cout << "inconsistency found!" << std::endl;
  } else {
    std::cout << "done" << std::endl;
  }
#endif 

  if(wiped_out) {
    ++statistics.num_failures;
    return false;
  } else return true;
  //return !wiped_out;
}


void Mistral::Solver::full_print() {
  for(int j=0; j<level; ++j) std::cout << "  ";
  std::cout << sequence << std::endl;
  for(unsigned int i=0; i<variables.size; ++i) {
    for(int j=0; j<level; ++j) std::cout << "  ";
    variables[i]->full_print();
    std::cout << std::endl;
  }
}

void Mistral::Solver::debug_print() {

  std::cout << std::endl << this << std::endl;

//   std::cout << "    currently changed: " << std::endl; 
//   int k = changed_objects.size;
//   while(k --> 0) {
//     changed_objects[k]->debug_print();
//     //std::cout << std::endl;
//   }
//   std::cout << std::endl;
  int i, k;
  for(i=var_trail_size.size-1; i>0; --i) {
    std::cout << "    changed at level " << (i-1) << ": " << std::endl; 
    k = var_trail_size[i];
    while(k --> var_trail_size[i-1])
      saved_vars[k]->debug_print();
    //std::cout << std::endl;
  }


  std::cout << "    AC Queue: " << std::endl; 
//   std::cout << triggered_constraints.active << std::endl;
//   std::cout << triggered_constraints.triggers[0] << std::endl;
//   std::cout << triggered_constraints.triggers[1] << std::endl;
//   std::cout << triggered_constraints.triggers[2] << std::endl;

  for(unsigned int i=0; i<3; ++i) {
    if(triggered_constraints.active.member(i)) {
      std::cout << "priority " << i << std::endl;
      int elt = triggered_constraints.triggers[i].first();
      while(elt != triggered_constraints.triggers[i]._head) {
	Constraint *cons = constraints[elt];
	std::cout << "\t" << cons << ": ";
	for(unsigned int j=0; j<cons->changes.size; ++j) {
	  int var = cons->changes[j];
	  int evt = cons->event_type[var];
	  std::cout << cons->scope[var] << "/" 
		    << (is_value(evt) ? "v" : (is_range(evt) ? "r" : "d") )
		    << (is_upper_bound(evt) ? "u" : "")
		    << (is_lower_bound(evt) ? "l" : "")
		    << " ";
	}
	std::cout << std::endl;
	elt = triggered_constraints.triggers[i].next[elt];
      } 
    }
  }

//   for(unsigned int i=0; i<triggered_constraints.size; ++i) {
//     Constraint *cons = constraints[triggered_constraints[i]];
//     std::cout << cons << ": ";
//     for(unsigned int j=0; j<cons->changes.size; ++j) {
//       int var = cons->changes[j];
//       int evt = cons->event_type[var];
//       std::cout << cons->scope[var] << "/" 
// 		<< (is_value(evt) ? "v" : (is_range(evt) ? "r" : "d") )
// 		<< (is_upper_bound(evt) ? "u" : "")
// 		<< (is_lower_bound(evt) ? "l" : "")
// 		<< " ";
//     }
//     std::cout << std::endl;
//   }
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

std::ostream& Mistral::Solver::display(std::ostream& os) const {
  os << "Variables:\n";
  ConstraintNode nd;
  for(unsigned int i=0; i<variables.size; ++i) {
    os << "\t" << variables[i] << " in " << variables[i]->domain ; //<< "\n";
    nd = variables[i]->constraints.first(_value_);
    while( variables[i]->constraints.next(nd) ) 
      os << " " << nd.elt.constraint;
    os << "\n";
  }

//   os << "\nConstraints:\n";
//   for(unsigned int i=0; i<constraints.size; ++i)
//     os << "\t" << constraints[i] << "\n";

  //return_str += ("\nSearch on "+toString(search.sequence)+"\n");
  return os;
}

void Mistral::Solver::branchRight() {
  backtrack();  
  Mistral::IntVar x = decision.pop(); //search.backtrack();

#ifdef _DEBUG_SEARCH
  std::cout << "c";
  for(unsigned int k=0; k<=decision.size; ++k) std::cout << " ";
  std::cout << " " << x << " in  " << (x->domain) << " =/= " << (x->domain.min) << std::endl;
#endif

  //x->branch->right();
  x->remove(x->domain.min);
}

void Mistral::Solver::branchLeft() {
  make_node();
  Mistral::IntVar x = decision.back(); //search.make_node();

#ifdef _DEBUG_SEARCH
  std::cout << "c";
  for(unsigned int k=0; k<=decision.size; ++k) std::cout << " ";
  std::cout << " " << x << " := " << (x->domain.min) << std::endl;
#endif

  //x->branch->left();
  x->setDomain(x->domain.min);
  
  //assert(sequence.member(x));
  //if(sequence.member(x)) 
  //sequence.erase(x); 
  //search.assign(x);
}


inline Mistral::Outcome Mistral::Solver::satisfied() {    
#ifdef _DEBUG_SEARCH
  std::cout << "c";
  for(unsigned int k=0; k<=decision.size; ++k) std::cout << " ";
  std::cout << " SAT!" << std::endl; 
#endif
  
  for(unsigned int i=0; i<sequence.capacity; ++i)
    solution[sequence[i]->id] = sequence[i]->domain.min;
  ++statistics.num_solutions;
  
  return SAT;
}

inline bool Mistral::Solver::limitsExpired() {
  
  return (parameters.limit && 
	  ((parameters.time_limit > 0.0 && (getRunTime() - statistics.start_time) > parameters.time_limit) ||
	   (parameters.node_limit > 0 && (statistics.num_nodes > parameters.node_limit)) ||
	   (parameters.fail_limit > 0 && (statistics.num_failures > parameters.fail_limit)) ||
	   (parameters.restart_limit > 0 && (statistics.num_restarts > parameters.restart_limit)) ||
	   (parameters.backtrack_limit > 0 && (statistics.num_backtracks > parameters.backtrack_limit))
	   ));
}

// void Mistral::Search::init_search(Vector< IntVar >& seq, VarOrdering *h, RestartPolicy *p) {
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
      
      //std::cout << (*this) << std::endl;
      for(unsigned int i=0; i<variables.size; ++i) {
	variables[i]->_assert_constraints_();
      }

      
      if( sequence.empty() ) status = satisfied();
      else branchLeft();
    } else {
      if( decision.empty() ) status = UNSAT;
      else if( limitsExpired() ) status = LIMITOUT;
      else branchRight();
    }
  }
  return status;
}


std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Solver& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Solver* x) {
  return x->display(os);
}
