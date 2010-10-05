
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
#include <mistral_variable.hpp>



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

Mistral::VarOrdering::VarOrdering(Mistral::Solver *s) 
  : length(s->sequence.size), variables(s->sequence.list_) { solver = s; }

Mistral::VarOrdering::~VarOrdering() {}

Mistral::NoOrder::NoOrder(Solver *s) 
  : VarOrdering(s) {}

Mistral::NoOrder::~NoOrder() {}

Mistral::IntVar Mistral::NoOrder::select() {
  return variables[length-1];
}

// Mistral::IntVar Mistral::Search::backtrack() {
//   unsigned int i = sequence.size;
//   trail_.pop(sequence.size);
//   while(i < sequence.size) sequence[i++]->unassign();
//   return decision.pop();
// }

// void Mistral::Search::assign(IntVar x) { 
//   if(sequence.member(x)) sequence.erase(x); 
// }

// void Mistral::Search::initialise(Solver *s) {
//   sequence.initialise(s->variables, false);
//   trail_.initialise(0, s->variables.size);
//   decision.initialise(0, s->variables.size);
// }

// void Mistral::Search::init_search(Vector< IntVar >& seq, VarOrdering *h, RestartPolicy *p) {
//   for(unsigned int i=seq.size; i;) 
//     sequence.insert(seq[--i]);    
//   if(heuristic) delete heuristic;
//   if(policy) delete policy;
//   //heuristic = (h ? h : new NoOrder(this));
//   policy = p;
// }


