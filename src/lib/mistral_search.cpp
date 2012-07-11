 
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


#include <mistral_search.hpp>

#include <mistral_sat.hpp>
#include <mistral_variable.hpp>



Mistral::ConsolidateListener::ConsolidateListener(Mistral::Solver *s) 
  : VariableListener(), ConstraintListener(), solver(s)
{
  sequence = &(solver->sequence);
  constraints.initialise(solver->variables.size);
  int n = solver->variables.size;
  for(int i=0; i<n; ++i) {
    Vector< Constraint > neighborhood;
    neighborhood.initialise(solver->variables[i].get_degree());
    for(int k=0; k<3; ++k) {
      for(int j=solver->constraint_graph[i].on[k].size-1; j>=0; --j)
     	neighborhood.add(solver->constraint_graph[i].on[k][j]);
    }
    constraints.add(neighborhood);
  }
}


Mistral::ConsolidateListener::~ConsolidateListener() {}

void Mistral::ConsolidateListener::notify_add_var() {
  Variable x = solver->variables.back();
  while((int)(constraints.size)<x.id()) {
    Vector< Constraint > neighborhood;
    constraints.add(neighborhood);
  }
  Vector< Constraint > neighborhood;
  neighborhood.initialise(x.get_degree());
  for(int k=0; k<3; ++k) {
    for(int j=solver->constraint_graph[x.id()].on[k].size-1; j>=0; --j)
      neighborhood.add(solver->constraint_graph[x.id()].on[k][j]);
  }
  constraints.add(neighborhood);
}

void Mistral::ConsolidateListener::notify_post (Constraint c) {};

void Mistral::ConsolidateListener::notify_relax(Constraint c) {};

void Mistral::ConsolidateListener::notify_add_con(Constraint c) {
  Variable *scope = c.get_scope();
  int arity = c.arity();
  for(int i=0; i<arity; ++i) {
    constraints[scope[i].id()].add(c);
    constraints[scope[i].id()].back().set_index(i);
  }
}

void Mistral::ConsolidateListener::notify_change(const int idx) {

  // std::cout << "BEG REACT TO CHANGE ON " << solver->variables[idx] << std::endl;

  Variable X = solver->variables[idx];
  int ids = sequence->index(idx);
  if(ids>=0) sequence->list_[ids] = X;

  for(int i=constraints[idx].size-1; i>=0; --i) {
    constraints[idx][i].consolidate_var();
  }

  // std::cout << "END REACT TO CHANGE ON " << solver->variables[idx] << std::endl;
}



//Mistral::LiteralActivityManager::LiteralActivityManager(Solver *s, void *a) 
Mistral::LiteralActivityManager::LiteralActivityManager(Solver *s) 
  : solver(s) {
  lit_activity = s->base->lit_activity.stack_;
  var_activity = s->base->var_activity.stack_;
  decay = solver->parameters.activity_decay;
  n_vars = solver->base->scope.size;
  solver->add((DecisionListener*)this);
}

Mistral::LiteralActivityManager::~LiteralActivityManager() {
  solver->remove((DecisionListener*)this);
}

double *Mistral::LiteralActivityManager::get_weight() { return var_activity; }     

void Mistral::LiteralActivityManager::notify_decision() {
  int i=n_vars;
  while(i--) {
    //std::cout << i << " " << var_activity[i] << " -> ";
    var_activity[i] *= decay;
    //std::cout << var_activity[i] << std::endl;
  }    
  i=2*n_vars;
  while(i--) lit_activity[i] *= decay;
}    

Mistral::RestartPolicy::RestartPolicy(const unsigned int b) {
  base = b;
}

Mistral::NoRestart::NoRestart() 
  : RestartPolicy(-1)
{
}

Mistral::NoRestart::~NoRestart() {}

Mistral::Geometric::Geometric(const unsigned int b, const double f) 
  : RestartPolicy(b)
{
  increment = b;
  factor = f;
}

Mistral::Geometric::~Geometric() {}

Mistral::Luby::Luby(const unsigned int b) 
  : RestartPolicy(b)
{
  iteration = 0;
}

Mistral::Luby::~Luby() {}

Mistral::NoOrder::NoOrder(Solver *s) 
  : solver(s) {}

Mistral::NoOrder::~NoOrder() {}

Mistral::Variable Mistral::NoOrder::select() {
  return solver->sequence.back();
}

// Mistral::Lexicographic::Lexicographic(Solver *s) 
//   : solver(s) {
//   index.initialise(s->variables.size);
//   solver->add(this);
//   last.initialise(0,s);
// }

Mistral::Lexicographic::Lexicographic(Solver *s) 
  : solver(s) {
  index.initialise(s->variables.size);
  solver->add(this);
  last.initialise(s,0);
}

void Mistral::Lexicographic::initialise(VarStack< Variable, ReversibleNum<int> >& seq) {
  int n = solver->variables.size;
  std::fill(index.stack_, index.stack_+n, -1);
  for(int i=0; i<seq.size; ++i) {
    index[seq[i].id()] = order.size;
    order.add(seq[i]);
  }
}

//void Mistral::Lexicographic::initialise(Solver *s, void *a) {
void Mistral::Lexicographic::initialise(Solver *s) {
  solver = s;
  index.initialise(s->variables.size);
  solver->add(this);
  last.initialise(s,0);
}

Mistral::Lexicographic::~Lexicographic() {}

Mistral::Variable Mistral::Lexicographic::select() {
  while(last<(int)(order.size) && order[last].is_ground()) { 
    ++last;
  }
  return order[last];
}

void Mistral::Lexicographic::notify_change(const int idx) {
  int ido = index[idx];
  if(ido>=0) {
    order[ido] = solver->variables[idx];
  }
}

// std::ostream& operator<<(std::ostream& os, Mistral::BranchingHeuristic& x) {
//   return x.display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::BranchingHeuristic* x) {
//   return x->display(os);
// }


std::ostream& operator<<(std::ostream& os, Mistral::DecisionListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RestartListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::SuccessListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::FailureListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::ConstraintListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::VariableListener& x) {
  return x.display(os);
}

// std::ostream& operator<<(std::ostream& os, Mistral::DecisionListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::RestartListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::SuccessListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::FailureListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::ConstraintListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::VariableListener* x) {
//   return x->display(os);
// }

double *weight_sorting_array;

int Mistral::decreasing_weight(const void *x, const void *y) {
    int _x = *(int*)x;
    int _y = *(int*)y;

    int r_value = 0;
    if(weight_sorting_array[_x] > weight_sorting_array[_y]) {
      r_value = -1;
    } else if(weight_sorting_array[_x] < weight_sorting_array[_y]) {
      r_value = 1;
    }

    return r_value;
  }

int log10(const double x){
  long int y = (long int)x;
  int log = 0;
  while(y>0) {
    ++log;
    y/=10;
  }
  return log;
}

std::ostream& Mistral::FailureCountManager::display(std::ostream& os, const bool all) const {
      
      int *all_variables = new int[variable_weight.size];
      int *all_constraints = new int[constraint_weight.size];


      int w, 
	xwidth; //, // = log10(solver->variables[variable_weight.size-1].id()),
	//cwidth; // = log10(solver->constraints[constraint_weight.size-1].id());
    
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	all_variables[i] = i;

	// w = log10(variable_weight[i]);
	// if(w>xwidth) xwidth = w;

      }

      for(unsigned int i=0; i<constraint_weight.size; ++i) {
	all_constraints[i] = i;

	// w = log10(constraint_weight[i]);
	// if(w>cwidth) cwidth = w;

      }


      weight_sorting_array = variable_weight.stack_;
      qsort(all_variables, variable_weight.size, sizeof(int), decreasing_weight);

      weight_sorting_array = constraint_weight.stack_;
      qsort(all_constraints, constraint_weight.size, sizeof(int), decreasing_weight);

      os << " c variable weight: \n c id: ";
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(variable_weight[all_variables[i]]);
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;

	  os << std::setw(xwidth) << all_variables[i] << " ";
	}
      }
      os << "\n c va: ";
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(variable_weight[all_variables[i]]);
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  os << std::setw(xwidth) << variable_weight[all_variables[i]] << " ";
	}
      }
      os << "\n c constraint weight: \n c id: ";
      for(unsigned int i=0; i<constraint_weight.size; ++i) {
	xwidth = log10(constraint_weight[all_constraints[i]]);
	w = log10(all_constraints[i]);
	if(w>xwidth) xwidth = w;

	os << std::setw(xwidth) << all_constraints[i] << " ";
      }
      os << "\n c va: ";
      for(unsigned int i=0; i<constraint_weight.size; ++i) {
	xwidth = log10(constraint_weight[all_constraints[i]]);
	w = log10(all_constraints[i]);
	if(w>xwidth) xwidth = w;

	os << std::setw(xwidth) << constraint_weight[all_constraints[i]] << " ";
      }
      os << std::endl;

      delete [] all_constraints;
      delete [] all_variables;

      return os;
    }    

std::ostream& Mistral::PruningCountManager::display(std::ostream& os, const bool all) const {
      
      int *all_variables = new int[variable_weight.size];


      int w, 
	xwidth; //, // = log10(solver->variables[variable_weight.size-1].id()),
	//cwidth; // = log10(solver->constraints[constraint_weight.size-1].id());
    
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	all_variables[i] = i;

	// w = log10(variable_weight[i]);
	// if(w>xwidth) xwidth = w;

      }


      weight_sorting_array = variable_weight.stack_;
      qsort(all_variables, variable_weight.size, sizeof(int), decreasing_weight);

      os << " c variable weight: \n c id: ";
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(variable_weight[all_variables[i]]);
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  os << std::setw(xwidth) << all_variables[i] << " ";
	}
      }
      os << "\n c va: ";
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(variable_weight[all_variables[i]]);
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  os << std::setw(xwidth) << variable_weight[all_variables[i]] << " ";
	}
      }
      os << std::endl;

    delete [] all_variables;

      return os;
    }    

std::ostream& operator<<(std::ostream& os, Mistral::MinDomain& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainOverDegree& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainOverWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinNeighborDomainOverWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinNeighborDomainOverNeighborWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::HalfSplit& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomMinMax& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeightValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::Guided& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::BoolMinWeightValue& x) {
  return x.display(os);
}
