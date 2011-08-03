
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


Mistral::LiteralActivityManager::LiteralActivityManager(Solver *s) 
  : solver(s) {
  lit_activity = s->base->lit_activity.stack_;
  var_activity = s->base->var_activity.stack_;
  decay = solver->parameters.activity_decay;
  n_vars = solver->base->scope.size;
  solver->add((DecisionListener*)this);
}

Mistral::LiteralActivityManager::~LiteralActivityManager() {}
    
double *Mistral::LiteralActivityManager::get_weight() { return var_activity; }     

void Mistral::LiteralActivityManager::notify_decision() {
  // int i=solver->sequence.size, a;
  // Variable *seq = solver->sequence.list_;
  // while(i--) {
  //   a = seq[i].id();
  //   var_activity[a] *= decay;
  //   lit_activity[2*a] *= decay;
  //   lit_activity[2*a+1] *= decay;
  // }
  int i=n_vars;
  while(i--) var_activity[i] *= decay;    
  i=n_vars;
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

Mistral::Lexicographic::Lexicographic(Solver *s) 
  : solver(s) {
  index.initialise(s->variables.size);
  solver->add(this);
  last.initialise(0,s);
}

void Mistral::Lexicographic::initialise(VarStack< Variable >& seq) {
  int n = solver->variables.size;
  std::fill(index.stack_, index.stack_+n, -1);
  for(unsigned int i=0; i<seq.size; ++i) {
    index[seq[i].id()] = order.size;
    order.add(seq[i]);
  }
}

void Mistral::Lexicographic::initialise(Solver *s) {
  solver = s;
  index.initialise(s->variables.size);
  solver->add(this);
  last.initialise(0,s);
}

Mistral::Lexicographic::~Lexicographic() {}

Mistral::Variable Mistral::Lexicographic::select() {
  while(last<order.size && order[last].is_ground()) { 
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



std::ostream& operator<<(std::ostream& os, Mistral::MinDomain& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainOverDegree& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainOverWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxWeight& x) {
  return x.display(os);
}
