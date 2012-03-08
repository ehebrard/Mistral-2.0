
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

  The author can be contacted electronically at 
  ehebrard@cse.unsw.edu.au.
*/


#include <mistral_solver.hpp>
#include <mistral_backtrack.hpp>
#include <mistral_constraint.hpp>


//void _save_() { if(first_change) { solver->save(this); first_change=false; trail_.push(value); } }

//void Mistral::Environment::save(Mistral::Reversible *r) { saved_objs.add(r); }

Mistral::Constraint::Constraint(ConstraintImplementation* p) { 
  propagator = p; 
  propagator->ConstraintImplementation::initialise();
  data = propagator->type; 
}

void Mistral::Constraint::re_link_to(Mistral::Trigger* t) {
  propagator->on[index()] = t;
}

void Mistral::Constraint::set_rank(const int idx) {
  int rank = (data&CTYPE);
  propagator->index[rank] = idx;
}

void Mistral::Constraint::set_index(const int idx) {
  data &= ITYPE;
  data |= idx;
}

void Mistral::Constraint::set_id(const int idx) {
  propagator->id = idx;
}

int Mistral::Constraint::id() const {
  return propagator->id;
}

void Mistral::Constraint::consolidate() {
  propagator->consolidate();
}

void Mistral::Constraint::consolidate_var() {
  propagator->consolidate_var(index());
}

void Mistral::Constraint::initialise(Solver *s) {
  propagator->initialise_vars(s);
  propagator->initialise();
}

int Mistral::Constraint::priority() const {
  return (global() ? ((GlobalConstraint*)propagator)->priority : 2);
}

void Mistral::Constraint::post(Solver* solver) { propagator->initial_post(solver); }

void Mistral::Constraint::awaken() { 
  if(binary()) 
    ((BinaryConstraint*)propagator)->post(); 
  else if(ternary())
    ((TernaryConstraint*)propagator)->post(); 
  else
    ((GlobalConstraint*)propagator)->post(); 
}

void Mistral::Constraint::trigger() { 
  if(binary())  {
    //((Solver*)solver)->active_constraints.trigger((BinaryConstraint*)propagator);
    ((BinaryConstraint*)propagator)->trigger(); 
  } else if(ternary()) {
    //((Solver*)solver)->active_constraints.trigger((TernaryConstraint*)propagator);
    ((TernaryConstraint*)propagator)->trigger(); 
  } else {
    //((Solver*)solver)->active_constraints.trigger((GlobalConstraint*)propagator);
    ((GlobalConstraint*)propagator)->trigger(); 
  }
}

int Mistral::Constraint::num_active() const { 
  int n;
  if(binary()) 
    n = size_byte[((BinaryConstraint*)propagator)->active]; 
  else if(ternary())
    n = size_byte[((TernaryConstraint*)propagator)->active]; 
  else
    n = ((GlobalConstraint*)propagator)->active.size; 
  return n;
}


Mistral::ConstraintImplementation* Mistral::Constraint::freeze() {
  if(global()) ((GlobalConstraint*)propagator)->freeze();
  return (idempotent() ? propagator : NULL);
}

Mistral::ConstraintImplementation* Mistral::Constraint::defrost() {
  if(global()) ((GlobalConstraint*)propagator)->defrost();
  return NULL;
}

Mistral::Variable* Mistral::Constraint::get_scope() { 
  return propagator->_scope.stack_;
}

int Mistral::Constraint::arity() const {
  return propagator->_scope.size;
} 

int Mistral::Constraint::get_active(const int i) const { 
  int n=0;
  
  if(binary()) {
    if(i) n = 1; // the second element must be 1
    else n = (((BinaryConstraint*)propagator)->active >> 1); 
    // the first element is the min
  }
  else if(ternary()) {
    if(i==2) n = 2; // the third element must be 2
    else if(i) { // the second element can be either 1 or 2
      n = 1+(((BinaryConstraint*)propagator)->active & 1);
    } else { // the first element 
      int a = ((TernaryConstraint*)propagator)->active; 
      if(!(a&1)) n = (4-(a&2))/2;
    }
  }
  else
    n = ((GlobalConstraint*)propagator)->active[i]; 
  
  return n;
}

int Mistral::Constraint::check(int* sol) {
  return propagator->check(sol);
}

Mistral::PropagationOutcome Mistral::Constraint::propagate() {
  return propagator->propagate();
}

Mistral::PropagationOutcome Mistral::Constraint::propagate(const Event evt) {
  return propagator->propagate(index(), evt);
}

void Mistral::Constraint::restore() {
  unsigned int mytype = index();

  // std::cout << "Restore " << (*this) << " " ;

  if(binary()) {

    // print_bitset(((BinaryConstraint*)propagator)->active, 0, std::cout);
    // std::cout << " -> "; 
    
    //if(mytype&4) ((BinaryConstraint*)propagator)->active = 3;
    //if(mytype&2) ((BinaryConstraint*)propagator)->un_post();
    //if(mytype&1) ((BinaryConstraint*)propagator)->un_relax();
    //((BinaryConstraint*)propagator)->restore(index());
    ((BinaryConstraint*)propagator)->restore(data);

    // print_bitset(((BinaryConstraint*)propagator)->active, 0, std::cout);
    // std::cout << std::endl;

  } else if(ternary()) {
    //if(mytype&4) ((TernaryConstraint*)propagator)->re_activate(); 
    //if(mytype&2) ((TernaryConstraint*)propagator)->un_post(); 
    //if(mytype&1) ((TernaryConstraint*)propagator)->un_relax();
    ((TernaryConstraint*)propagator)->restore(data);
  } else {
    //if(mytype&4) ((GlobalConstraint*)propagator)->re_activate(); 
    if(mytype&2) ((GlobalConstraint*)propagator)->un_post(); 
    if(mytype&1) ((GlobalConstraint*)propagator)->un_relax();
    //((GlobalConstraint*)propagator)->restore(data);
  }
}

int Mistral::Constraint::get_backtrack_level() {
  if(global()) return ((GlobalConstraint*)propagator)->get_backtrack_level();
  else return propagator->solver->level-1;
}

Mistral::Decision Mistral::Constraint::get_decision() {
 if(global()) return ((GlobalConstraint*)propagator)->get_decision();
 else {
   Decision dec = ((Solver*)(propagator->solver))->decisions.back(0); 
   dec.invert();
   return dec;
 }
 //return propagator->solver->decisions.back(0);
}

std::ostream& Mistral::Constraint::display(std::ostream& os) const {
  propagator->display(os);
  // os 
  //   //<< ":" 
  //   << "[" << id() << "]";
  return os;
}

void Mistral::Environment::_restore_() {
  
  unsigned int previous_level;
  
  previous_level = trail_.pop();
  while( saved_ints.size > previous_level ) 
    saved_ints.pop()->restore();
  
  previous_level = trail_.pop();
  while( saved_lists.size > previous_level ) 
    saved_lists.pop()->restore();
  
  previous_level = trail_.pop();
  while( saved_bools.size > previous_level ) 
    *(saved_bools.pop()) = 3;
  
  previous_level = trail_.pop();
  previous_level = trail_.pop();
  
  --level;
  
}


std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::ReversibleBool& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::ReversibleBool* x) {
  return x->display(os);
}
