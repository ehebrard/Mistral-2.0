
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

#include <math.h>

#include <mistral_solver.hpp>
#include <mistral_variable.hpp>
#include <mistral_constraint.hpp>
 





std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Constraint& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Constraint* x) {
  return x->display(os);
}


std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::ConstraintImplementation& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::ConstraintImplementation* x) {
  return (x ? x->display(os) : os << "Null");
}


std::ostream& Mistral::operator<< (std::ostream& os,  Mistral::ConstraintTriggerArray& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os,  Mistral::ConstraintTriggerArray* x) {
  return x->display(os);
}



// void Mistral::Constraint::set_rank(const int idx) {
//   int rank = (data^CTYPE);
//   propagator->index[rank] = idx;
// }

// void Mistral::Constraint::set_id(const int idx) {
//   propagator->id = idx;
// }

// int Mistral::Constraint::id() {
//   return propagator->id;
// }


//   void Mistral::Constraint::post(solver*) { propagator->initial_post(solver); }

//   int Mistral::Constraint::check(int* sol) {
//     return propagator->check(sol);
//   }

//   Mistral::PropagationOutcome Mistral::Constraint::propagate() {
//     return propagator->propagate();
//   }

//   Mistral::PropagationOutcome Mistral::Constraint::propagate(const int changed_idx, const Event evt) {
//     return propagator->propagate(changed_idx, evt);
//   }

// void Mistral::Constraint::restore() {}


Mistral::ConstraintTriggerArray::ConstraintTriggerArray() {}

Mistral::ConstraintTriggerArray::ConstraintTriggerArray(const int size) {
  initialise(size);
}

void Mistral::ConstraintTriggerArray::initialise(const int size) {
  for(int i=0; i<3; ++i) on[i].initialise(size);
}

Mistral::ConstraintTriggerArray::~ConstraintTriggerArray() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete constraint trigger array" << std::endl;
#endif
}

std::ostream& Mistral::ConstraintTriggerArray::display(std::ostream& os) const {
  for(int i=2; i>=0; --i) {
    os << "["; 
    for(unsigned int j=0; j<on[i].size; ++j) {
      os << on[i][j] ; //<< ":"
	//<< on[i][j].index();
      if(j<on[i].size-1) os << ", ";
    }
    os << "]";
  }
  return os;
}


Mistral::ConstraintImplementation::ConstraintImplementation() {
  id = -1;
  //arity = 0;
  //trigger = NULL;
  self = NULL;
  index = NULL;
  enforce_nfc1 = true;
}

// Mistral::ConstraintImplementation::ConstraintImplementation(const int a) {
//   id = -1;
//    //arity = a;
// //   //num_triggers = 0;
// //   //self = new Constraint[arity];
// //   //index = new int[arity];
// //   //trigger = new Trigger*[arity];
// // }

Mistral::ConstraintImplementation::~ConstraintImplementation() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete constraint implementation" << std::endl;
#endif
  delete [] self;
  delete [] index;
  //delete [] trigger;
}

void Mistral::ConstraintImplementation::trigger_on(const int t, Variable x) {

  // display(std::cout);
  // std::cout << " triggers on " 
  // 	    << (t ? (t>1 ? "domain" : "range") : "value")  << " of " << x << std::endl;

  // //std::cout << "on[" << on.size << "] = ";

  
  //std::cout << x.domain_type << " " << CONST_VAR << std::endl;

  if(x.domain_type != CONST_VAR) {
    //Solver *solver = x.get_solver();

    // std::cout << (*this) << " triggers on "  << t << " "
    // 	    << (int*)(&(solver->constraint_graph[x.id()].on[t])) << std::endl;
    
    
    on.add(&(get_solver()->constraint_graph[x.id()].on[t]));
  } else {
    on.add(NULL);
  }
  _scope.add(x);
  //std::cout << on.back() << std::endl;

}

bool Mistral::ConstraintImplementation::is_triggered_on(const int i, const int t) {

  // std::cout << (*this) << " " << t << std::endl;

  // std::cout << (int*)(on[i]) << " >=? "
  // 	    << (int*)(&(get_solver()->constraint_graph[_scope[i].id()].on[t])) 
  // 	    << (on[i] >= &(get_solver()->constraint_graph[_scope[i].id()].on[t]) ? " yes" : " no") << std::endl;

  


  return( on[i] >= &(get_solver()->constraint_graph[_scope[i].id()].on[t]) );
}


// int Mistral::Constraint::get_backtrack_level() {
//   return solver->level-1;
// }

// Mistral::Decision Mistral::ConstraintImplementation::get_decision() { 
//   Decision dec = get_solver()->decisions.back(0); 
//   dec.invert();
//   return dec;
// }

void Mistral::GlobalConstraint::initialise_vars(Solver *s) {
  for(unsigned int i=0; i<scope.size; ++i) {
    scope[i].initialise(s, 1);

    //std::cout << "INITIALISE " << scope[i] << " => " << get_solver()->variables[scope[i].id()] << std::endl;

  }
}

void Mistral::ConstraintImplementation::initial_post(Solver *s) {

// #ifdef _DEBUG_RELAX
//   std::cout << "[" << std::setw(4) << id << "]: first post on: " ;
// #endif

//   solver = s;
//   // for each of its variables
//   self = new Constraint[on.size];
//   index = new int[on.size];
//   for(unsigned int i=0; i<on.size; ++i) {
//     //_scope[i].initialise(s, false);
//     self[i] = Constraint(this, i|type);
//     //post_on(i);
    
//     Constraint c = self[i];
//     c.data |= POSTED;

//     solver->save( c );

//     index[i] = -1;

//     //index[i] = on[i]->post(self[i]);
//     //if(_scope[i].domain_type != CONST_VAR) {
//     //if(!(_scope[i].is_constant())) {
//     if(!(_scope[i].is_ground())) {
      
//       //std::cout << "yes" << std::endl;

// #ifdef _DEBUG_RELAX
//       std::cout << _scope[i] << " " ;
// #endif

//       un_relax_from(i);
      
      
//     }  
//     else {
      
//       //std::cout << "no" << std::endl;

//       desactivate(i);
       
//     } 

//     // else {

//     //   //std::cout << "no" << std::endl;
//     // }
//   }

// #ifdef _DEBUG_RELAX
//   std::cout << std::endl;

// #endif

//   //mark_domain();


#ifdef _DEBUG_RELAX
  std::cout << "[" << std::setw(4) << id << "]: first post on: " ;
#endif

  solver = s;
  // for each of its variables
  self = new Constraint[on.size];
  index = new int[on.size];

  // First we go through the variables to check whether they are ground.
  // If so, we "desactivate" the corresponding var index.
  // Also, if there is only one active variable and the constraint enforces nfc1, we do not post it at all.
  int nb_actives = on.size;
  for(unsigned int i=0; i<on.size; ++i) {
    index[i] = -1;
    self[i] = Constraint(this, i|type);
    if(_scope[i].is_ground()) {
      --nb_actives;
      desactivate(i);
    }
  }

  // Now we post the constraint on the active variables (provided that there are at least 2)
  if(!enforce_nfc1 || nb_actives>1) {
    Constraint c;

    for(unsigned int i=0; i<on.size; ++i) {
      if(!(_scope[i].is_ground())) {
	
        c = self[i];
	c.data |= POSTED;

	solver->save( c );

#ifdef _DEBUG_RELAX
	std::cout << _scope[i] << " " ;
#endif

	un_relax_from(i);
      }  
    }
  }

#ifdef _DEBUG_RELAX
  std::cout << std::endl;
#endif

  //mark_domain();
}

int Mistral::ConstraintImplementation::get_trigger_type(const int i) {
  return (int)(on[i] - &(get_solver()->constraint_graph[_scope[i].id()].on[0]));
}
 
void Mistral::ConstraintImplementation::set_scope(const int i, Variable x) {
  // std::cout << "set scope[" << i << "] of ";
  // display(std::cout);
  // std::cout << " was " << _scope[i] << " now " << x << std::endl;

  // std::cout << "previous trigger list: " << on[i] << " / " << (int)(on[i] - &get_solver()->constraint_graph[_scope[i].id()].on[0]) << std::endl; 



  // std::cout << "trigger was: " << on[i] << " re-link to " ;

  // std::cout << "\njjj: " << get_solver()->constraint_graph[x.id()].on[0] << std::endl;
  // std::cout << "\njjj: " << get_solver()->constraint_graph[x.id()].on[1] << std::endl;
  // std::cout << "\njjj: " << get_solver()->constraint_graph[x.id()].on[2] << std::endl;

  // std::cout << get_trigger_type(i) << std::endl;

  on[i] = &(get_solver()->constraint_graph[x.id()].on[get_trigger_type(i)]);


  //std::cout << on[i] << std::endl;


  _scope.stack_[i] = x;

  // on[i] = 
  // on.add(&(get_solver()->constraint_graph[x.id()].on[t]));
}

int Mistral::Trigger::post(Constraint ct) {
  add(ct);
  return size-1;
}

void Mistral::Trigger::relax(const int idx) {
  //std::cout << "relax " << idx << " out of " << size << std::endl;
  if(idx != (int)(size-1)) stack_[size-1].set_rank(idx);
  remove(idx);
}




Mistral::PropagationOutcome Mistral::BinaryConstraint::propagate() {

  if(!support[0]) initialise_supports();

  PropagationOutcome wiped = CONSISTENT; 
  for(int i=0; IS_OK(wiped) && i<2; ++i) {
    wiped = BinaryConstraint::propagate(i, RANGE_EVENT);
  }
  return wiped;

}


Mistral::PropagationOutcome Mistral::BinaryConstraint::bound_propagate() {

  if(!support[0]) initialise_supports();

  PropagationOutcome wiped = CONSISTENT; 
  for(int i=0; IS_OK(wiped) && i<2; ++i) {
    wiped = BinaryConstraint::bound_propagate(i, RANGE_EVENT);
  }
  return wiped;

}


Mistral::PropagationOutcome Mistral::BinaryConstraint::bound_propagate(const int changed_idx, const Event evt) {
  PropagationOutcome wiped = CONSISTENT; 
  int vali, vnext;
  int revise_idx = 1-changed_idx;
  
  if(RANGE_CHANGED(evt)) {
    
#ifdef _DEBUG_GENPROPAG
    std::cout << " - (bc) revise " << scope[revise_idx] << " in " << scope[revise_idx].get_domain() 
	      << " w.r.t. " << scope[changed_idx] 
	      << std::endl;
#endif

    if(!support[0]) initialise_supports();
    
    vnext = scope[revise_idx].get_min();
    do {
      vali = vnext;
      vnext = scope[revise_idx].next(vali);
      
#ifdef _DEBUG_GENPROPAG
      std::cout << "  *find a support for " << scope[revise_idx] << " = " << vali << std::endl;
#endif
      
      if( !(scope[changed_idx].contain(support[revise_idx][vali][changed_idx])) &&
	  !find_bound_support(revise_idx, vali) ) {
	
#ifdef _DEBUG_GENPROPAG
	std::cout << "  => none found!" << std::endl;
#endif
	
	if(IS_FAIL(scope[revise_idx].remove(vali))) {
	  wiped = FAILURE(revise_idx);
	} 
      }
      else { 
#ifdef _DEBUG_GENPROPAG
     	std::cout << "  => ok" << std::endl;
#endif
	break;
      } 
    } while( vali < vnext );
    
    if( IS_OK(wiped) && vali < vnext ) {

      vnext = scope[revise_idx].get_max();
      do {
	vali = vnext;
	vnext = scope[revise_idx].prev(vali);
	
#ifdef _DEBUG_GENPROPAG
	std::cout << "  *find a support for " << scope[revise_idx] << " = " << vali << std::endl;
#endif
      
	if( !(scope[changed_idx].contain(support[revise_idx][vali][changed_idx])) &&
	    !find_bound_support(revise_idx, vali) ) {
	
#ifdef _DEBUG_GENPROPAG
	  std::cout << "  => none found!" << std::endl;
#endif
	  
	  if(IS_FAIL(scope[revise_idx].remove(vali))) {
	    wiped = FAILURE(revise_idx);
	  } 
	}
	else { 
#ifdef _DEBUG_GENPROPAG
	  std::cout << "  => ok" << std::endl;
#endif
	  break;
	} 
      } while( vali > vnext );
      
    }
    
#ifdef _DEBUG_GENPROPAG
    if(!(IS_OK(wiped))) {
      std::cout << "FAIL!!" << std::endl;
    }
#endif
  }

  return wiped;
}


Mistral::PropagationOutcome Mistral::BinaryConstraint::propagate(const int changed_idx, const Event evt) {
  PropagationOutcome wiped = CONSISTENT; 
  int vali, vnext;
  int revise_idx = 1-changed_idx;

#ifdef _DEBUG_GENPROPAG
    std::cout << " - revise " << scope[revise_idx] << " in " << scope[revise_idx].get_domain() 
	      << " w.r.t. " << scope[changed_idx] 
	      << std::endl;
#endif

  if(!support[0]) initialise_supports();

  vnext = scope[revise_idx].get_min();
  do {
    vali = vnext;
    vnext = scope[revise_idx].next(vali);

#ifdef _DEBUG_GENPROPAG
    std::cout << "  *find a support for " << scope[revise_idx] << " = " << vali << std::endl;
#endif

    if( !(scope[changed_idx].contain(support[revise_idx][vali][changed_idx])) &&
	!find_support(revise_idx, vali) ) {

#ifdef _DEBUG_GENPROPAG
	std::cout << "  => none found!" << std::endl;
#endif

      if(IS_FAIL(scope[revise_idx].remove(vali))) {
	wiped = FAILURE(revise_idx);
      } 
    } 
#ifdef _DEBUG_GENPROPAG
      else {
	std::cout << "  => ok" << std::endl;
      }
#endif

  } while( vali < vnext );

#ifdef _DEBUG_GENPROPAG
  if(!(IS_OK(wiped))) {
    std::cout << "FAIL!!" << std::endl;
  }
#endif

  return wiped;
}


// void Mistral::BinaryConstraint::use_residual_support() {
//   support = new bituple*[2];
//   for(int var=0; var<2; ++var) {
//     support[var] = new int*[]
//   }
// }

bool Mistral::BinaryConstraint::find_support(const int revise_idx, const int vli) 
{
  int changed_idx = 1-revise_idx, vali = -INFTY; //, *solution = support[revise_idx][vli];

  solution[revise_idx] = vli;
  solution[changed_idx] = scope[changed_idx].get_min();

  //std::cout << "   " << solution[0] << " " << solution[1] << std::endl;

  while( vali<solution[changed_idx] && check( solution ) ) {

    vali = solution[changed_idx];
    solution[changed_idx] = scope[changed_idx].next( vali );

    //if(vali<solution[changed_idx]) std::cout << "   " << solution[0] << " " << solution[1] << std::endl;

  } 
  return vali<solution[changed_idx];
}


// bool print_sol(int *sol, Mistral::Variable var, int val, Mistral::Variable against, 
// 	       int id, Mistral::BinaryConstraint* c) {
//   std::cout << "check [" << id << "] " << var << " in " << var.get_domain() << " = " << val 
// 	    << " against " << against << " in " << against.get_domain() << " | ";
//   c->display(std::cout);
//   std::cout << ": {" << sol[0] << ", " << sol[1] << "}\n";
//   return true;
// }


bool Mistral::BinaryConstraint::find_bound_support(const int revise_idx, const int vli) 
{
  int changed_idx = 1-revise_idx;
  int max_val = scope[changed_idx].get_max(); 

  solution[revise_idx] = vli;
  
  for( solution[changed_idx] = scope[changed_idx].get_min(); 
       //print_sol(solution, scope[revise_idx], vli, scope[changed_idx], id, this) &&
	 solution[changed_idx] <= max_val && check( solution );
       ++solution[changed_idx]);
  
  return solution[changed_idx]<=max_val;
}


Mistral::PropagationOutcome Mistral::TernaryConstraint::propagate() {
  if(!support[0]) initialise_supports();

  PropagationOutcome wiped = CONSISTENT; 

  for(int i=0; IS_OK(wiped) && i<3; ++i) {
    wiped = TernaryConstraint::propagate(i, RANGE_EVENT);
  }

  return wiped;
}

Mistral::PropagationOutcome Mistral::TernaryConstraint::bound_propagate() {
  if(!support[0]) initialise_supports();

  PropagationOutcome wiped = CONSISTENT; 

  for(int i=0; IS_OK(wiped) && i<3; ++i) {
    wiped = TernaryConstraint::bound_propagate(i, RANGE_EVENT);
  }

  return wiped;
}

Mistral::PropagationOutcome Mistral::TernaryConstraint::propagate(const int changed_idx, const Event evt) {
  PropagationOutcome wiped = CONSISTENT; 
  int vali, vnext;
  int revise_idx, other_idx;
  int x, y;

  if(!support[0]) initialise_supports();


#ifdef _DEBUG_GENPROPAG
  std::cout << " propagate " << event2str(evt) << " on " << scope[changed_idx] << std::endl;
#endif


  for(int i=1; i<3; ++i) {
    revise_idx = (changed_idx+i)%3;
    other_idx = (changed_idx+3-i)%3;


#ifdef _DEBUG_GENPROPAG
    std::cout << " -revise " << scope[revise_idx] << " in " << scope[revise_idx].get_domain() 
	      << " w.r.t. " << scope[changed_idx] << " & " << scope[other_idx]
	      << std::endl;
#endif

    vnext = scope[revise_idx].get_min();

    //    std::cout << vnext << std::endl;

    do {
      vali = vnext;
      vnext = scope[revise_idx].next(vali);

#ifdef _DEBUG_GENPROPAG
      std::cout << "  *find a support for " << scope[revise_idx] << " = " << vali << std::endl;
      std::cout.flush();
#endif
      x = support[revise_idx][vali][changed_idx];
      y = support[revise_idx][vali][other_idx];


      // std::cout << std::endl << x << " " << y << std::endl;


      
      if( (x == NOVAL 
	   || y == NOVAL
	   || !scope[changed_idx].contain(x)
	   || !scope[other_idx].contain(y)
	   ) &&
	  !find_support(revise_idx, vali) ) {

#ifdef _DEBUG_GENPROPAG
	std::cout << "  => none found!" << std::endl;
#endif

	if(IS_FAIL(scope[revise_idx].remove(vali))) {
	  wiped = FAILURE(revise_idx);
	} 
      } 
#ifdef _DEBUG_GENPROPAG
      else {
	std::cout << "  => ok" << std::endl;
      }
#endif

    } while( vali < vnext );
  }

#ifdef _DEBUG_GENPROPAG
  if(!(IS_OK(wiped))) {
    std::cout << "FAIL!!" << std::endl;
  }
#endif

  return wiped;
}



Mistral::PropagationOutcome Mistral::TernaryConstraint::bound_propagate(const int changed_idx, const Event evt) {
  PropagationOutcome wiped = CONSISTENT; 
  int vali, vnext;
  int revise_idx, other_idx;
  int x, y;

  if(!support[0]) initialise_supports();

#ifdef _DEBUG_GENPROPAG
  std::cout << " bound propagate " << event2str(evt) << " on " << scope[changed_idx] << std::endl;
#endif


  for(int i=1; i<3; ++i) {
    revise_idx = (changed_idx+i)%3;
    other_idx = (changed_idx+3-i)%3;


#ifdef _DEBUG_GENPROPAG
    std::cout << " - (bc) revise " << scope[revise_idx] << " in " << scope[revise_idx].get_domain() 
	      << " w.r.t. " << scope[changed_idx] << " & " << scope[other_idx]
	      << std::endl;
#endif

    vnext = scope[revise_idx].get_min();
    do {
      vali = vnext;
      vnext = scope[revise_idx].next(vali);

#ifdef _DEBUG_GENPROPAG
      std::cout << "  *find a support for " << scope[revise_idx] << " = " << vali << std::endl;
      std::cout.flush();
#endif
      x = support[revise_idx][vali][changed_idx];
      y = support[revise_idx][vali][other_idx];
      
      if( (x == NOVAL 
	   || y == NOVAL
	   || !scope[changed_idx].contain(x)
	   || !scope[other_idx].contain(y)
	   ) &&
	  !find_bound_support(revise_idx, vali) ) {

#ifdef _DEBUG_GENPROPAG
	std::cout << "  => none found!" << std::endl;
#endif

	if(IS_FAIL(scope[revise_idx].remove(vali))) {
	  wiped = FAILURE(revise_idx);
	} 
      } 
      else {
#ifdef _DEBUG_GENPROPAG
	std::cout << "  => ok" << std::endl;
#endif
	break;
      }
    } while( vali < vnext );

    if( IS_OK(wiped) && vali < vnext ) {

      vnext = scope[revise_idx].get_max();
      do {
	vali = vnext;
	vnext = scope[revise_idx].prev(vali);

#ifdef _DEBUG_GENPROPAG
	std::cout << "  *find a support for " << scope[revise_idx] << " = " << vali << std::endl;
	std::cout.flush();
#endif
	x = support[revise_idx][vali][changed_idx];
	y = support[revise_idx][vali][other_idx];
	
	if( (x == NOVAL 
	     || y == NOVAL
	   || !scope[changed_idx].contain(x)
	     || !scope[other_idx].contain(y)
	     ) &&
	    !find_bound_support(revise_idx, vali) ) {
	  
#ifdef _DEBUG_GENPROPAG
	  std::cout << "  => none found!" << std::endl;
#endif
	  
	  if(IS_FAIL(scope[revise_idx].remove(vali))) {
	    wiped = FAILURE(revise_idx);
	  } 
	} 
	else {
#ifdef _DEBUG_GENPROPAG
	  std::cout << "  => ok" << std::endl;
#endif
	  break;
	}
      } while( vali > vnext );
    }
  }

#ifdef _DEBUG_GENPROPAG
  if(!(IS_OK(wiped))) {
    std::cout << "FAIL!!" << std::endl;
  }
#endif

  return wiped;
}

bool Mistral::TernaryConstraint::find_support(const int revise_idx, const int vli) 
{
  bool no_support = true;

  int x_idx = (revise_idx+1)%3, y_idx = (revise_idx+2)%3;

  int valx = scope[x_idx].get_min(), valy; //, *solution = support[revise_idx][vli];
  
  solution[revise_idx] = vli;


  do {
    solution[x_idx] = valx;
    valy = scope[y_idx].get_min();

    do {
      solution[y_idx] = valy;

#ifdef _DEBUG_GENPROPAG
      std::cout << "    " << solution[0] << " " << solution[1] << " " << solution[2] << std::endl;
#endif

      no_support = check(solution); 

      if(no_support) 
	valy = scope[y_idx].next( valy );

    } while( solution[y_idx] < valy );
    
    if(no_support)
      valx = scope[x_idx].next( valx );

  } while( solution[x_idx] < valx );


  return !no_support;
}


bool Mistral::TernaryConstraint::find_bound_support(const int revise_idx, const int vli) 
{
  bool no_support = true;

  int x_idx = (revise_idx+1)%3, y_idx = (revise_idx+2)%3;


  int min_y_val = scope[y_idx].get_min();
  int max_y_val = scope[y_idx].get_max();
  int max_x_val = scope[x_idx].get_max();
  
  solution[revise_idx] = vli;


  for(solution[x_idx] = scope[x_idx].get_min();
      no_support && solution[x_idx] <= max_x_val;
      ++solution[x_idx]
      ) {
    for(solution[y_idx] = min_y_val;
	no_support && solution[y_idx] <= max_y_val;
	++solution[y_idx]
	) {
    
#ifdef _DEBUG_GENPROPAG
      std::cout << "    " << solution[0] << " " << solution[1] << " " << solution[2] << std::endl;
#endif

      no_support = check(solution); 
      
    }
  }

  return !no_support;
}


void Mistral::BinaryConstraint::trigger()  { get_solver()->active_constraints.trigger(this); }
void Mistral::TernaryConstraint::trigger() { get_solver()->active_constraints.trigger(this); }
void Mistral::GlobalConstraint::trigger()  { get_solver()->active_constraints.trigger(this); }


Mistral::GlobalConstraint::GlobalConstraint(Vector< Variable > scp) {
  for(unsigned int i=0; i<scp.size; ++i) scope.add(scp[i]);
}
Mistral::GlobalConstraint::GlobalConstraint(std::vector< Variable > scp) {
  for(std::vector< Variable >::iterator vi=scp.begin(); vi!=scp.end(); ++vi) scope.add(*vi);
}
Mistral::GlobalConstraint::GlobalConstraint(Variable* scp, const int n) {
  for(int i=0; i<n; ++i) scope.add(scp[i]);
}

int Mistral::GlobalConstraint::get_backtrack_level() {
  return solver->level-1;
}

Mistral::Decision Mistral::GlobalConstraint::get_decision() { 
  Decision dec = get_solver()->decisions.back(0); 
  dec.invert();
  return dec;
}


Mistral::PropagationOutcome Mistral::GlobalConstraint::propagate(const int changed_idx, const Event evt) {
  PropagationOutcome wiped = CONSISTENT; 
  unsigned int i;
  int vali, vnext;

  for( i=0; IS_OK(wiped) && i<scope.size; ++i ) {
    if( (int)i != changed_idx ) { 
      vnext = scope[i].get_min();
      do {
	vali = vnext;
	vnext = scope[i].next(vali);
	
#ifdef _DEBUG_GENPROPAG
	std::cout << "find a support for " << scope[i] << " = " << vali << std::endl;
#endif

	if( ( !first_support(i, vali) && !find_support(i, vali) ) ) {
	  if(IS_FAIL(scope[i].remove(vali))) {
	    wiped = FAILURE(i);
	  } 
	}
      } while( vali < vnext );
    }	
  }
  return wiped;
}

Mistral::PropagationOutcome Mistral::GlobalConstraint::propagate() {
  PropagationOutcome wiped = CONSISTENT; 
  unsigned int i, j;
  int vali, vnext;


  while(!changes.empty()) {
    j = changes.pop();
    for( i=0; IS_OK(wiped) && i<scope.size; ++i ) {
      if( i != j ) { 
	vnext = scope[i].get_min();
	do {
	  vali = vnext;
	  vnext = scope[i].next(vali);

#ifdef _DEBUG_GENPROPAG
	  // if(solver->parameters.verbosity) {
	  std::cout << "find a support for " << scope[i] << " = " << vali << std::endl;
	  // }
#endif

	  if( ( !first_support(i, vali) && !find_support(i, vali) ) ) {
	    if(IS_FAIL(scope[i].remove(vali))) {
	      wiped = FAILURE(i);
	    } else if(changes.list_ == events.list_) {
	      if(!changes.contain(i)) changes.add(i);
	    }
	  }
	} while( vali < vnext );
      }	
    }
  }
  return wiped;
}


Mistral::PropagationOutcome Mistral::GlobalConstraint::bound_propagate(const int changed_idx, const Event evt) {
  PropagationOutcome wiped = CONSISTENT; 
  unsigned int i;


  int vali, valmax;//, vnext;
  bool supported;
  for( i=0; IS_OK(wiped) && i<scope.size; ++i ) {
    if( (int)i != changed_idx ) { 
      supported = false;
      vali = scope[i].get_min();
      valmax = scope[i].get_max();
      while(!supported && IS_OK(wiped) && vali<=valmax) {
	
	//std::cout << "find support for " << scope[i] << "=" << vali << std::endl;
	
	if( ( !first_support(i, vali) && !find_bound_support(i, vali) ) ) {
	  if(IS_FAIL(scope[i].remove(vali))) {
	    wiped = FAILURE(i);
	  } else if(changes.list_ == events.list_) {
	    if(!changes.contain(i)) changes.add(i);
	  }
	} else supported = true;
	
	++vali;
      }
      
      if(supported && vali<=valmax) {
	supported = false;
	vali = valmax;
	while(!supported && IS_OK(wiped)) {
	  
	  //std::cout << "find support for " << scope[i] << "=" << vali << std::endl;
	  
	  if( ( !first_support(i, vali) && !find_bound_support(i, vali) ) ) {
	    if(IS_FAIL(scope[i].remove(vali))) {
		wiped = FAILURE(i);
	    }  else if(changes.list_ == events.list_) {
	      if(!changes.contain(i)) changes.add(i);
	    }
	  } else supported = true;
	  --vali;
	}
      }
    }	
  }
  return wiped;
}


Mistral::PropagationOutcome Mistral::GlobalConstraint::bound_propagate() {
  PropagationOutcome wiped = CONSISTENT;
  
  unsigned int i, j;

  int vali, valmax;//, vnext;
  bool supported;

  while(!changes.empty()) {
    j = changes.pop();

#ifdef _DEBUG_GENPROPAG
  std::cout << " propagate " << event2str(event_type[j]) << " on " << scope[j] << std::endl;
#endif

    for( i=0; IS_OK(wiped) && i<scope.size; ++i ) {
      //if( i != j ) { 

#ifdef _DEBUG_GENPROPAG
    std::cout << " - (bc) revise " << scope[i] << " in " << scope[i].get_domain() 
	      << " w.r.t. " << scope[j] 
	      << std::endl;
#endif

	supported = false;
	vali = scope[i].get_min();
	valmax = scope[i].get_max();
	while(!supported && IS_OK(wiped) && vali<=valmax) {

#ifdef _DEBUG_GENPROPAG
	  std::cout << "  *find a support for " << scope[i] << " = " << vali << std::endl;
	  std::cout.flush();
#endif
	  
	  if( ( !first_support(i, vali) && !find_bound_support(i, vali) ) ) {
	    
#ifdef _DEBUG_GENPROPAG
	    std::cout << "  => none found!" << std::endl;
#endif

	    if(IS_FAIL(scope[i].remove(vali))) {
	      wiped = FAILURE(i);
	    } else if(changes.list_ == events.list_) {
	      if(!changes.contain(i)) changes.add(i);
	    }
	  } else {

#ifdef _DEBUG_GENPROPAG
	    std::cout << "  => ok" << std::endl;
#endif

	    supported = true;
	  }
	  ++vali;
	}

	if(supported && vali<=valmax) {
	  supported = false;
	  vali = valmax;
	  while(!supported && IS_OK(wiped)) {

#ifdef _DEBUG_GENPROPAG
	  std::cout << "  *find a support for " << scope[i] << " = " << vali << std::endl;
	  std::cout.flush();
#endif

	    if( ( !first_support(i, vali) && !find_bound_support(i, vali) ) ) {

#ifdef _DEBUG_GENPROPAG
	  std::cout << "  *find a support for " << scope[i] << " = " << vali << std::endl;
	  std::cout.flush();
#endif

	      if(IS_FAIL(scope[i].remove(vali))) {
		wiped = FAILURE(i);
	      }  else if(changes.list_ == events.list_) {
		if(!changes.contain(i)) changes.add(i);
	      }
	    } else {

#ifdef _DEBUG_GENPROPAG
	      std::cout << "  => ok" << std::endl;
#endif
	      
	      supported = true;
	    }
	    --vali;
	  }
	}
	//}	
    }
  }
  
  return wiped;
}


bool Mistral::GlobalConstraint::first_support(const int vri, const int vli) 
{
  int j;
  if( supports && supports[vri][vli][0] != NOVAL ) {
    j=scope.size;
    while( j-- ) 
      if( vri != j )
	if (!scope[j].contain( supports[vri][vli][j] )) break;
    if( j < 0 ) 
      return true;
  } 
  j=scope.size;
  while( j-- ) 
    solution[j] = scope[j].get_min();
  solution[vri] = vli; 

  return false;
}

bool Mistral::GlobalConstraint::find_support(const int vri, const int vli) 
{
  int i=scope.size, vali;
  bool found=false;
  // sol is initialized: either to the value 
  // a variable is already assigned to
  // or to the first value in its domain
  while(i >= 0) {

#ifdef _DEBUG_GENPROPAG
    // if(solver->parameters.verbosity) {
      std::cout << "\t<" << solution[0] ;
      for(unsigned int k=1; k<scope.size; ++k) {
    	std::cout << " " << solution[k];
      }
      std::cout << "> ";
    // }
#endif

    // check this assignment
    if( !check( solution ) ) {

#ifdef _DEBUG_GENPROPAG
      // if(solver->parameters.verbosity) {
      std::cout << "OK!" << std::endl;
      // }
#endif

      found=true;
      if( supports ) {
	vali = scope.size;
	while( vali-- )
	  supports[vri][vli][vali] = solution[vali];
      }
      break;
    } 

#ifdef _DEBUG_GENPROPAG
    else 
      //if(solver->parameters.verbosity) {
      std::cout << "NO" << std::endl;
    //     }
#endif
 
    // try to assign more things
    // find the last var whose domain we have not exhausted
    --i;
    while( i >= 0 ) {
      if( i == vri || scope[i].is_ground() ) {
	--i;
	continue;
      }
      vali = scope[i].next( solution[i] );
      if(vali > solution[i]) {
	solution[i] = vali;
      	break;
      } else {
	solution[i] = scope[i].get_min();
      }
      --i;
    }
    if( i >= 0 )
      i = scope.size;
  } 
  return found;
}


bool Mistral::GlobalConstraint::find_bound_support(const int vri, const int vli) 
{
  int i=scope.size, vali;
  bool found=false;
  // sol is initialized: either to the value 
  // a variable is already assigned to
  // or to the first value in its domain
  while(i >= 0) {


    // std::cout << "\t<" << solution[0] ;
    // for(unsigned int k=1; k<scope.size; ++k) {
    //   std::cout << " " << solution[k];
    // }
    // std::cout << "> ";

    // check this assignment
    if( !check( solution ) ) {

      //std::cout << "OK!" << std::endl;

      found=true;
      if( supports ) {
	vali = scope.size;
	while( vali-- )
	  supports[vri][vli][vali] = solution[vali];
      }
      break;
    } //  else {
    //   std::cout << "NO" << std::endl;
    // }
 
    // try to assign more things
    // find the last var whose domain we have not exhausted
    --i;
    while( i >= 0 ) {
      if( i == vri || scope[i].is_ground() ) {
	--i;
	continue;
      }
      if(solution[i]>=scope[i].get_max())
	solution[i] = scope[i].get_min();
      else {
	++solution[i];
      	break;
      }
      --i;
    }
    if( i >= 0 )
      i = scope.size;
  } 
  return found;
}


Mistral::GlobalConstraint::~GlobalConstraint() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete global constraint" << std::endl;
#endif
  if(changes.list_ == events.list_)
    events.list_ = NULL;
  delete [] event_type;
  //delete [] scope;
}

void Mistral::GlobalConstraint::initialise() {
  // trail_.initialise(0,2*scope.size);
  // trail_.add(-1);

  event_type = new Event[on.size];
  // self = new int[scope.size];
  // index = new Constraint**[scope.size];
  // trigger = new int[scope.size];
  solution = new int[scope.size];
  changes.initialise(0, scope.size-1, false);
  supports = NULL;

  std::fill(event_type, event_type+on.size, NO_EVENT);
  // std::fill(self, self+scope.size, -1);
  // std::fill(index, index+scope.size, (Constraint**)NULL);
  // std::fill(trigger, trigger+scope.size, _DOMAIN_);
  std::fill(solution, solution+scope.size, 0);

  active.initialise(solver, 0, on.size-1, true);

  GlobalConstraint::set_idempotent();
}


std::ostream& Mistral::ConstraintImplementation::display(std::ostream& os) const {
  return os;
}

std::ostream& Mistral::GlobalConstraint::display(std::ostream& os) const {
  os << name() << "(" << scope[0];
  for(unsigned int i=1; i<scope.size; ++i)
    os << ", " << scope[i]/*.get_var()*/;
  os << ")";
  return os;
}

void Mistral::ConstraintNotEqual::initialise() {
  ConstraintImplementation::initialise();

  //std::cout << "TYPE AT INITIALISATION  " << type << std::endl;

  trigger_on(_VALUE_, scope[0]);
  trigger_on(_VALUE_, scope[1]);
  //set_idempotent(true);
  //postponed = false;
}

void Mistral::ConstraintNotEqual::mark_domain() {
  get_solver()->mark_non_convex(scope[0].id());
  get_solver()->mark_non_convex(scope[1].id());
}

Mistral::PropagationOutcome Mistral::ConstraintNotEqual::rewrite() {
  RewritingOutcome r_evt = NO_EVENT; 

  if(scope[0].is_ground() || scope[1].is_ground() || !(scope[0].intersect(scope[1]))) {

#ifdef _DEBUG_REWRITE
    std::cout << "    relax ";
    display(std::cout);
    std::cout << std::endl;
#endif

    r_evt = SUPPRESSED;
    relax();
  } else {

    if( scope[0].get_min() == 0 && 
	scope[1].get_min() == 0 &&
	scope[0].get_max() == 1 && 
	scope[1].get_max() == 1 &&
	scope[0].is_expression() &&
	scope[1].is_expression() ) {
      
      r_evt = SUPPRESSED;
      
#ifdef _DEBUG_REWRITE
      std::cout << "    relax ";
      display(std::cout);
      std::cout << std::endl;
#endif
      
      relax();
      
      Constraint con = Constraint( new PredicateNot(scope[0], scope[1]) );
      
#ifdef _DEBUG_REWRITE
      std::cout << "    post " << con << std::endl;
#endif

      get_solver()->add( con );
    }
  }  

  return r_evt;
}


Mistral::PropagationOutcome Mistral::ConstraintNotEqual::propagate() {

  // std::cout << "ACTIVE: " << active <<  " " ;
  // std::cout.flush();
  // std::cout << scope[0].get_domain() << " " << scope[1].get_domain() << std::endl;

  PropagationOutcome wiped_idx = CONSISTENT;
  if(!active) {
    if(scope[0].get_min() == scope[1].get_min()) wiped_idx = FAILURE(0);
  } else {
    
    // if(active != 3) {
    //   std::cout << "remove the value of " << scope[2-active] << ": " 
    // 		<< scope[2-active].get_domain() << " from " 
    // 		<< scope[active-1] << " in " << scope[active-1].get_domain() << std::endl;
    // }

    
    //print_active();
    
    // if(active != 3)
    //   std::cout << "remove " << scope[2-active].get_domain() << " from " 
    // 		<< scope[active-1].get_domain() << std::endl;


    if(active != 3 && (scope[active-1].remove(scope[2-active].get_min()) == FAIL_EVENT)) wiped_idx = FAILURE(active-1);
  }
  return wiped_idx;
}

Mistral::PropagationOutcome Mistral::ConstraintNotEqual::propagate(const int changed_idx, 
								   const Event evt) {

  //std::cout << "\n propagate changes on " << scope[changed_idx] << " to " << scope[1-changed_idx] << std::endl;

  //active^=(1 << changed_idx);
  int var = 1-changed_idx;
  PropagationOutcome wiped_idx = (scope[var].remove(scope[changed_idx].get_min()) == FAIL_EVENT ? FAILURE(var) : CONSISTENT);
  //if(IS_OK(wiped_idx)) relax_from(var);
  return wiped_idx;
}

std::ostream& Mistral::ConstraintNotEqual::display(std::ostream& os) const {
  os << scope[0]/*.get_var()*/ << " =/= " << scope[1]/*.get_var()*/;
  return os;
}

void Mistral::ConstraintEqual::initialise() {
  ConstraintImplementation::initialise();

  //Constraint::initialise();
  trigger_on(_DOMAIN_, scope[0]);
  trigger_on(_DOMAIN_, scope[1]);
  //set_idempotent(true);
}

Mistral::PropagationOutcome Mistral::ConstraintEqual::rewrite() {

// #ifdef _DEBUG_REWRITE
//       std::cout << "REWRITE EQUALITY " ;
//       display(std::cout);
//       std::cout << std::endl;
//       //std::cout << solver->variables << std::endl;
// #endif

  RewritingOutcome r_evt = SUPPRESSED; 

  //Mistral::PropagationOutcome wiped = propagate(); 


  /// CHANGE THAT:
  // relax all constraints around one of the variables and post clones instead


  // // std::cout << "propag: " << scope[0] << " in " << scope[0].get_domain() 
  // // 	    << " = " << scope[1] << " in " << scope[1].get_domain() << std::endl;


  // // std::cout << 11 << std::endl;

  // // std::cout << scope[0] << std::endl;

  // std::cout << domain2str(scope[0].domain_type) 
  // 	    << " "
  // 	    << (int*)(scope[0].variable)
  // 	    << std::endl;


  // // std::cout << scope[1] << std::endl;

  // std::cout << domain2str(scope[1].domain_type) 
  // 	    << " ";


  // std::cout << (int*)(scope[1].variable)
  // 	    << std::endl;

  // // 

  // std::cout << 11 << std::endl;

  // std::cout << (scope[1].is_ground()) << std::endl;

  // std::cout << 22 << std::endl;

  if( !(scope[1].is_ground()) && !(scope[0].is_ground()) ) {

    //std::cout << 33 << std::endl;

    if( scope[0].is_expression() && scope[1].is_expression() ) {

// #ifdef _DEBUG_REWRITE
//           std::cout << "ACTUALLY REWRITE " ;
//           display(std::cout);
//           std::cout << std::endl;
//           // //std::cout << solver->variables << std::endl;
// #endif 

      relax();

      int k[2], j;

      k[0] = scope[0].id();
      k[1] = scope[1].id();

      j = get_solver()->constraint_graph[k[0]].size() > 
   	get_solver()->constraint_graph[k[1]].size();
      // j is 1 iff degree(x0) > degree(x1), 
      // in that case we want to transfer x1's constraints on x0.
      // x[j] <- x[1-j]


      //std::cout << "k[j]=" << k[j] << std::endl;

      get_solver()->remove(scope[j]);
      //get_solver()->domain_types[k[j]] |= REMOVED_VAR;

#ifdef _DEBUG_REWRITE
      std::cout << "    relax constraints on " << scope[j]/*.get_var()*/ << " and post them on " << scope[1-j]/*.get_var()*/ << std::endl;
#endif

      Constraint con;
      for(Event trig = 0; trig<3; ++trig) {
	for(int i = get_solver()->constraint_graph[k[j]].on[trig].size; i--;) {
	  //for(int t=0; t<3; ++t) {
	  //for(int i=trigger[j]->on[t].size; --i;) {

	  //con = trigger[j]->on[t][i];


	  con = get_solver()->constraint_graph[k[j]].on[trig][i];

#ifdef _DEBUG_REWRITE
	  std::cout << "      relax " << con << std::endl;
#endif


	  con.relax();
	  con.set_scope(con.index(), scope[1-j]);
  	  get_solver()->add(con);


#ifdef _DEBUG_REWRITE
	  std::cout << "      post " << con  << std::endl;
#endif

	  // std::cout << "trig[" << k[1-j] << "][0]: " << get_solver()->constraint_graph[k[1-j]].on[0] << std::endl;
	  // std::cout << "trig[" << k[1-j] << "][1]: " << get_solver()->constraint_graph[k[1-j]].on[1] << std::endl;
	  // std::cout << "trig[" << k[1-j] << "][2]: " << get_solver()->constraint_graph[k[1-j]].on[2] << std::endl;

	  
	  
	}
      }

  //     ConstraintNode nd;
  //     nd = solver->constraint_graph[k[j]]->first(_VALUE_);
  //     while( solver->constraint_graph[k[j]]->next(nd) ) {	
  // 	nd.elt.constraint->relax();
  // 	nd.elt.constraint->scope[nd.elt.index] = scope[1-j];
  // 	solver->add(nd.elt.constraint);
  //     }
      
      //and now scope[j] points to scope[1-j]

      // SELF CHANGE
      //scope[j].expression->self = scope[1-j];
      
      scope[j].expression->id = scope[1-j].id();

// #ifdef _DEBUG_REWRITE
//       //std::cout << solver->variables << std::endl;
//       std::cout << "EQUALITY REWRITEN " ;
//       display(std::cout);
//       std::cout << std::endl;
// #endif

    } 
  } else {
    relax();
  }
  
  return r_evt;
}


Mistral::PropagationOutcome Mistral::ConstraintEqual::propagate() {
  PropagationOutcome wiped = CONSISTENT;
  if(active) {
    if(scope[1].set_domain(scope[0]) == FAIL_EVENT) wiped = FAILURE(1);
    if(IS_OK(wiped) && scope[0].set_domain(scope[1]) == FAIL_EVENT) wiped = FAILURE(0);
  } else if(scope[0].get_min() != scope[1].get_min()) {
    wiped = FAILURE(0);
  }

  return wiped;
}

Mistral::PropagationOutcome Mistral::ConstraintEqual::propagate(const int changed_idx, const Event evt) {
  return(scope[1-changed_idx].set_domain(scope[changed_idx]) == FAIL_EVENT ? FAILURE(1-changed_idx) : CONSISTENT);
}


std::ostream& Mistral::ConstraintEqual::display(std::ostream& os) const {
  os << (scope[0]/*.get_var()*/) << " == " << (scope[1]/*.get_var()*/);
  return os;
}

void Mistral::PredicateUpperBound::initialise() {
  ConstraintImplementation::initialise();
  trigger_on(_RANGE_, scope[0]);
  trigger_on(_VALUE_, scope[1]);
}

Mistral::Constraint Mistral::PredicateUpperBound::get_negation(const int var, Variable x) { 
  return Constraint( new PredicateLowerBound( scope[0], x, bound-1 ) );
}

Mistral::PropagationOutcome Mistral::PredicateUpperBound::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[1].is_ground() ) {
    if( scope[1].get_min() ) { // x[0] <= k 
      if( IS_FAIL(scope[0].set_max(bound)) )
	wiped = FAILURE(0);
    } else { // x[0] > k
      if( IS_FAIL(scope[0].set_min(bound+1)) )
	wiped = FAILURE(0);
    } 
  } else {
    if( scope[0].get_min() > bound ) {
      if( IS_FAIL(scope[1].set_domain(0)) ) return FAILURE(1);
    } else if( scope[0].get_max() <= bound ) {
      if( IS_FAIL(scope[1].remove(0)) ) return FAILURE(1);
    }
  }
  
  return wiped;
}


Mistral::PropagationOutcome Mistral::PredicateUpperBound::propagate(const int changed_idx, const Event evt) {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( changed_idx ) {
    if( LB_CHANGED(evt) ) { // x[0] <= k 
      if( IS_FAIL(scope[0].set_max(bound)) ) wiped = FAILURE(0);
    } else { // x[0] > k
      if( IS_FAIL(scope[0].set_min(bound+1)) ) wiped = FAILURE(0);
    } 
  } else {
    if( LB_CHANGED(evt) && scope[0].get_min() > bound ) {
      if( IS_FAIL(scope[1].set_domain(0)) ) return FAILURE(1);
    } else if( UB_CHANGED(evt) && scope[0].get_max() <= bound ) {
      if( IS_FAIL(scope[1].remove(0)) ) return FAILURE(1);
    }
  }
  
  return wiped;
}

std::ostream& Mistral::PredicateUpperBound::display(std::ostream& os) const {
  os << scope[1]/*.get_var()*/ << " <=> (" << scope[0]/*.get_var()*/ << " <= " << bound << ")";
  return os;
}


void Mistral::PredicateLowerBound::initialise() {
  ConstraintImplementation::initialise();
  trigger_on(_RANGE_, scope[0]);
  trigger_on(_VALUE_, scope[1]);
}

Mistral::PropagationOutcome Mistral::PredicateLowerBound::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[1].is_ground() ) {
    if( scope[1].get_min() ) { // x[0] >= k 
      if( IS_FAIL(scope[0].set_min(bound)) )
	wiped = FAILURE(0);
    } else { // x[0] < k
      if( IS_FAIL(scope[0].set_max(bound-1)) )
	wiped = FAILURE(0);
    } 
  } else {
    if( scope[0].get_max() < bound ) {
      if( IS_FAIL(scope[1].set_domain(0)) ) return FAILURE(1);
    } else if( scope[0].get_min() >= bound ) {
      if( IS_FAIL(scope[1].remove(0)) ) return FAILURE(1);
    }
  }
  
  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateLowerBound::propagate(const int changed_idx, const Event evt) {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( changed_idx ) {
    if( LB_CHANGED(evt) ) { // x[0] >= k 
      if( IS_FAIL(scope[0].set_min(bound)) )
	wiped = FAILURE(0);
    } else { // x[0] < k
      if( IS_FAIL(scope[0].set_max(bound-1)) )
	wiped = FAILURE(0);
    } 
  } else {
    if( UB_CHANGED(evt) && scope[0].get_max() < bound ) {
      if( IS_FAIL(scope[1].set_domain(0)) ) return FAILURE(1);
    } else if( LB_CHANGED(evt) && scope[0].get_min() >= bound ) {
      if( IS_FAIL(scope[1].remove(0)) ) return FAILURE(1);
    }
  }
  
  return wiped;
}

std::ostream& Mistral::PredicateLowerBound::display(std::ostream& os) const {
  os << scope[1]/*.get_var()*/ << " <=> (" << scope[0]/*.get_var()*/ << " >= " << bound << ")";
  return os;
}

void Mistral::PredicateLess::initialise() {
  ConstraintImplementation::initialise();
  trigger_on(_RANGE_, scope[0]);
  trigger_on(_RANGE_, scope[1]);
  trigger_on(_VALUE_, scope[2]);
}

//x0 + 3 > x1

Mistral::PropagationOutcome Mistral::PredicateLess::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[2].is_ground() ) {
    if( scope[2].get_min() ) { // x[0] + k <= x[1]
      if( IS_FAIL(scope[0].set_max(scope[1].get_max()-offset)) )
	wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].set_min(scope[0].get_min()+offset)) )
	wiped = FAILURE(1);
    } else if( scope[2].get_max() == 0 ) { 
      // x[0] + k > x[1]
      // x[0] + k - 1 >= x[1]
      // x[1] + (1-k) <= x[0]
      if( IS_FAIL(scope[0].set_min(scope[1].get_min()-offset+1)) )
	wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].set_max(scope[0].get_max()+offset-1)) )
	wiped = FAILURE(1);
    } 
  } else {
    if( scope[0].get_min() + offset > scope[1].get_max() ) {
      if( IS_FAIL(scope[2].set_domain(0)) ) return FAILURE(2);
    } else if( scope[0].get_max() + offset <= scope[1].get_min() ) {
      if( IS_FAIL(scope[2].remove(0)) ) return FAILURE(2);
    }
  }
  
  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateLess::propagate(const int changed_idx, const Event evt) {      
  Mistral::PropagationOutcome wiped = CONSISTENT;


  // std::cout << "HERE: "  << event2str(evt) << " on " << scope[changed_idx] << " in " << scope[changed_idx].get_domain() << std::endl;

  if( changed_idx == 2 ) {
    if( LB_CHANGED(evt) ) {
      // x[0] + offset <= x[1]
      if( IS_FAIL(scope[0].set_max(scope[1].get_max()-offset)) ) wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].set_min(scope[0].get_min()+offset)) ) wiped = FAILURE(1);
    } else {
      // x[0] + offset > x[1]
      if( IS_FAIL(scope[0].set_min(scope[1].get_min()-offset+1)) ) wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].set_max(scope[0].get_max()+offset-1)) ) wiped = FAILURE(1);
    }
  } else {
    if( scope[2].is_ground() ) {
      // it is already a precedence
      if( scope[2].get_min() ) {
	// x[0] + offset <= x[1]
	if(changed_idx) {
	  if(UB_CHANGED(evt) && IS_FAIL(scope[0].set_max(scope[1].get_max() - offset))) wiped = FAILURE(0);
	} else {
	  if(LB_CHANGED(evt) && IS_FAIL(scope[1].set_min(scope[0].get_min() + offset))) wiped = FAILURE(1);
	}
      } else {

	// std::cout << "APPLY " << scope[0] << " in " << scope[0].get_domain() << " + " << offset
	// 	  << " > "
	// 	  << scope[1] << " in " << scope[1].get_domain() << std::endl;

	// x[0] + offset > x[1]
	// x[1] + (1-k) <= x[0]	
	if(changed_idx) {
	  if(LB_CHANGED(evt) && IS_FAIL(scope[0].set_min(scope[1].get_min() - offset + 1))) wiped = FAILURE(0);
	} else {
	  if(UB_CHANGED(evt) && IS_FAIL(scope[1].set_max(scope[0].get_max() + offset - 1))) wiped = FAILURE(1);
	}
      }
    } else {
      // check if it must be a precedence
      if(changed_idx) {
	if(LB_CHANGED(evt)) {
	  // x[1] is too big to be less than or equal to  (x[0] + k - 1) ?
	  if(scope[1].get_min() > (scope[0].get_max() + offset -1) &&
	     IS_FAIL(scope[2].set_domain(1))) wiped = FAILURE(2); 
	}
	if(UB_CHANGED(evt)) {
	  // x[1] is too small to be greater than or equal to  (x[0] + k) ?
	  if(scope[1].get_max() < (scope[0].get_min() + offset) &&
	     IS_FAIL(scope[2].set_domain(0))) wiped = FAILURE(2); 
	}
      } else {
	if(LB_CHANGED(evt)) {
	  // x[0] is too big to be less than or equal to (x[1] - k) ?
	  if(scope[0].get_min() > (scope[1].get_max() - offset) &&
	     IS_FAIL(scope[2].set_domain(0))) wiped = FAILURE(2); 
	}  
      	if(UB_CHANGED(evt)) {
	  // x[0] is too small to be greater than or equal to (x[1] + 1 - k) ?
	  if(scope[0].get_max() < (scope[1].get_min() + 1 - offset) &&
	     IS_FAIL(scope[2].set_domain(1))) wiped = FAILURE(2); 
	}
      }
    }
  }

  return wiped;
}

std::ostream& Mistral::PredicateLess::display(std::ostream& os) const {
  os << scope[2]/*.get_var()*/ << " <=> (" << scope[0]/*.get_var()*/ ;
  if(offset==0) os << " <= ";
  else if(offset==1) os << " < ";
  else os << " + " << offset << " <= ";
  os << scope[1]/*.get_var()*/ << ")";
  return os;
}

void Mistral::ConstraintLess::initialise() {
  ConstraintImplementation::initialise();
  //Constraint::initialise();
  trigger_on(_RANGE_, scope[0]);
  trigger_on(_RANGE_, scope[1]);
  //set_idempotent(true);
}

Mistral::PropagationOutcome Mistral::ConstraintLess::propagate() {

//   std::cout << "changes: " << changes << std::endl;
//   std::cout << "events[0]: " << events[0] << " events[0]: " << events[1] << std::endl;

//   if(scope[0].id() == 6 && scope[1].id() == 9) {
//   std::cout << "propagate " << this << std::endl;
//   for(unsigned int i=0; i<scope.size; ++i)
//     std::cout << " " << scope[i].get_domain();
//   std::cout << std::endl;
//   }
  Mistral::PropagationOutcome wiped = CONSISTENT;
  if(scope[1].set_min(scope[0].get_min() + offset) == FAIL_EVENT) wiped = FAILURE(1);
  if(IS_OK(wiped) && scope[0].set_max(scope[1].get_max() - offset) == FAIL_EVENT) wiped = FAILURE(0);
  

//   if(scope[0].id() == 6 && scope[1].id() == 9) {
//     for(unsigned int i=0; i<scope.size; ++i)
//       std::cout << " " << scope[i].get_domain();
//     std::cout << std::endl;
//   }

  return wiped;
}

Mistral::PropagationOutcome Mistral::ConstraintLess::propagate(const int changed_idx, const Event evt) {

  Mistral::PropagationOutcome wiped = CONSISTENT;
  if(changed_idx==0) {
    if(LB_CHANGED(evt) && scope[1].set_min(scope[0].get_min() + offset) == FAIL_EVENT) wiped = FAILURE(1);
  } else {
    if(UB_CHANGED(evt) && scope[0].set_max(scope[1].get_max() - offset) == FAIL_EVENT) wiped = FAILURE(0);
  }

  return wiped;
}

std::ostream& Mistral::ConstraintLess::display(std::ostream& os) const {
  os << scope[0]/*.get_var()*/;
  if(offset < 0) os << " - " << (-offset+1) << " < ";
  else if(offset > 1) os << " + " << (offset-1) << " < ";
  else if(offset > 0) os << " < ";
  else os << " <= ";
  os << scope[1]/*.get_var()*/;
  return os;
}

Mistral::ConstraintDisjunctive::ConstraintDisjunctive(Variable x, Variable y, const int p0, const int p1) 
  : BinaryConstraint(x,y) { 
  processing_time[0] = p0; 
  processing_time[1] = p1;
  enforce_nfc1 = false;
}

Mistral::ConstraintDisjunctive::ConstraintDisjunctive(Vector< Variable >& scp, const int p0, const int p1) 
  : BinaryConstraint(scp) { 
  processing_time[0] = p0; 
  processing_time[1] = p1;
  enforce_nfc1 = false;
}

Mistral::ConstraintDisjunctive::ConstraintDisjunctive(std::vector< Variable >& scp, const int p0, const int p1) 
  : BinaryConstraint(scp) { 
  processing_time[0] = p0; 
  processing_time[1] = p1;
  enforce_nfc1 = false;
}

void Mistral::ConstraintDisjunctive::initialise() {
  
  precedence[0] = new ConstraintLess(scope[0], scope[1], processing_time[0]);
  precedence[1] = new ConstraintLess(scope[1], scope[0], processing_time[1]);

  ConstraintImplementation::initialise();
  trigger_on(_RANGE_, scope[0]);
  trigger_on(_RANGE_, scope[1]);
  //set_idempotent(true);
  //stress = 0;
}



void Mistral::ConstraintDisjunctive::decide(const int choice) {

#ifdef _DEBUG_RELAX
      std::cout << "[" << std::setw(4) << id << "](" << name() << "): force relax" << std::endl;
#endif

  relax();
  //active = 0;

  if(choice==1) {

#ifdef _DEBUG_DISJUNCTIVE
    std::cout << "c add precedence " << precedence[1] << std::endl;
#endif

    get_solver()->add(precedence[1]);
  } else {

#ifdef _DEBUG_DISJUNCTIVE
    std::cout << "c add precedence " << precedence[0] << std::endl;
#endif

    get_solver()->add(precedence[0]);
  }

}

Mistral::PropagationOutcome Mistral::ConstraintDisjunctive::propagate() {
  int hold = 3;

  // check if prec[1] is violated (x1 + p1 > x0).
  if(scope[1].get_min()+processing_time[1] > scope[0].get_max()) {
    hold &= 2;
  }

  // check if prec[1] is violated (x0 + p0 > x1).
  if(scope[0].get_min()+processing_time[0] > scope[1].get_max()) {
    hold &= 1;
  }

  if(!hold) return FAILURE(0);
  if(hold<3) {
    decide(hold);
  }
  return CONSISTENT;
}

Mistral::PropagationOutcome Mistral::ConstraintDisjunctive::propagate(const int changed_idx, 
								      const Event evt) {
  int hold = 3;


  // if(id == 11) {
  //   std::cout << "propagate " ;
  //   display(std::cout);
  //   std::cout << " because of " << event2str(evt) << " on " << scope[changed_idx] << std::endl;
  // }

// if(id == 11) 
//   std::cout << scope[1].get_min() << " + " << processing_time[1] << " ?> " << scope[0].get_max() << std::endl;

  // check if prec[1] is violated (x1 + p1 > x0).
  //if(changed_idx && LB_CHANGED(evt))
  if(scope[1].get_min()+processing_time[1] > scope[0].get_max()) {
    hold &= 2;
  }

// if(id == 11) 
//   std::cout << scope[0].get_min() << " + " << processing_time[0] << " ?> " << scope[1].get_max() << std::endl;

  // check if prec[1] is violated (x0 + p0 > x1).
  if(scope[0].get_min()+processing_time[0] > scope[1].get_max()) {
    hold &= 1;
  }

// if(id == 11) {
//   std::cout << "HOLD: " << hold << std::endl;
//  }

  if(!hold) return FAILURE(0);
  if(hold<3) {
    decide(hold);
  }
  return CONSISTENT;
}

void Mistral::ConstraintDisjunctive::consolidate() {
  for(unsigned int i=0; i<2; ++i) {
    scope[i] = scope[i]/*.get_var()*/;
    _scope[i] = _scope[i]/*.get_var()*/;
    precedence[i].consolidate();
    // precedence[0].get_scope()[i] = scope[i];
    // precedence[1].get_scope([1-i] = scope[i];
  }
}

std::ostream& Mistral::ConstraintDisjunctive::display(std::ostream& os) const {
  os << precedence[0] << " or " 
     << precedence[1] ;
  return os;
}



// Mistral::ConstraintTernaryDisjunctive::ConstraintTernaryDisjunctive(Vector< Variable >& scp, const int p0, const int p1) : Constraint(scp) { 
//   processing_time[0] = p0; 
//   processing_time[1] = p1;
// }

// Mistral::ConstraintTernaryDisjunctive::ConstraintTernaryDisjunctive(std::vector< Variable >& scp, const int p0, const int p1) : Constraint(scp) { 
//   processing_time[0] = p0; 
//   processing_time[1] = p1;
// }

// void Mistral::ConstraintTernaryDisjunctive::initialise() {
  
//   precedence[0] = new ConstraintLess(scope, processing_time[0]);
//   precedence[1] = new ConstraintLess(processing_time[1]);
//   precedence[1]->add(scope[1]);
//   precedence[1]->add(scope[0]);

//   Constraint::initialise();
//   trigger_on(_RANGE_, scope[0]);
//   trigger_on(_RANGE_, scope[1]);
//   trigger_on(_VALUE_, scope[2]);
//   set_idempotent(true);
//   stress = 0;
// }

// Mistral::PropagationOutcome Mistral::ConstraintTernaryDisjunctive::decide(const int choice) {
//   PropagationOutcome wiped = CONSISTENT;

//   if(choice==1) {
//     if(scope[2].set_domain(1) == FAIL_EVENT) wiped = FAILURE(2);
//     else {
//       relax();
//       solver->add(precedence[1]);
//     }
//   } else {
//     if(scope[2].set_domain(0) == FAIL_EVENT) wiped = FAILURE(2);
//     else {
//       relax();
//       solver->add(precedence[0]);
//     }
//   }

//   return wiped;
// }

// void Mistral::ConstraintTernaryDisjunctive::decide() {
//   relax();

//   if(scope[2].get_min()) {
//     solver->add(precedence[1]);
//   } else {
//     solver->add(precedence[0]);
//   }

// }

// Mistral::PropagationOutcome Mistral::ConstraintTernaryDisjunctive::propagate() {
//   PropagationOutcome wiped = CONSISTENT;
//   int hold = 3;

//   if(events.contain(2)) {
    
//     decide();
 
//   } else {
    
//     // check if prec[1] is violated (x1 + p1 > x0).
//     if(scope[1].get_min()+processing_time[1] > scope[0].get_max()) {
//       hold &= 2;
//     }
    
//     // check if prec[1] is violated (x0 + p0 > x1).
//     if(scope[0].get_min()+processing_time[0] > scope[1].get_max()) {
//       hold &= 1;
//     }
    
//     if(!hold) return FAILURE(0);
//     if(hold<3) {
//       wiped = decide(hold);
//     }

//   }

//   return wiped;
// }

// void Mistral::ConstraintTernaryDisjunctive::consolidate() {
//   for(unsigned int i=0; i<2; ++i) {
//     scope[i] = scope[i]/*.get_var()*/;
//     precedence[0]->scope[i] = scope[i];
//     precedence[1]->scope[1-i] = scope[i];
//   }
//   scope[2] = scope[2]/*.get_var()*/;
// }

// std::ostream& Mistral::ConstraintTernaryDisjunctive::display(std::ostream& os) const {
//   os << precedence[0] << " or " 
//      << precedence[1] ;
//   return os;
// }

Mistral::ConstraintReifiedDisjunctive::ConstraintReifiedDisjunctive(Variable x, Variable y, Variable z, const int p0, const int p1) 
  : TernaryConstraint(x,y,z) { 
  processing_time[0] = p0; 
  processing_time[1] = p1;
}

Mistral::ConstraintReifiedDisjunctive::ConstraintReifiedDisjunctive(Vector< Variable >& scp, const int p0, const int p1) 
  : TernaryConstraint(scp) { 
  processing_time[0] = p0; 
  processing_time[1] = p1;
}

Mistral::ConstraintReifiedDisjunctive::ConstraintReifiedDisjunctive(std::vector< Variable >& scp, const int p0, const int p1) 
  : TernaryConstraint(scp) { 
  processing_time[0] = p0; 
  processing_time[1] = p1;

}

void Mistral::ConstraintReifiedDisjunctive::initialise() {

  ConstraintImplementation::initialise();

  trigger_on(_RANGE_, scope[0]);
  trigger_on(_RANGE_, scope[1]);
  trigger_on(_VALUE_, scope[2]);
  //set_idempotent(true);
  //stress = 0;

  state = scope[2].get_var().bool_domain;

  if(scope[0].domain_type == RANGE_VAR) {
    min_t0_ptr = &(scope[0].range_domain->min);
    max_t0_ptr = &(scope[0].range_domain->max);
  } else if(scope[0].domain_type == BITSET_VAR) {
    min_t0_ptr = &(scope[0].bitset_domain->domain.min);
    max_t0_ptr = &(scope[0].bitset_domain->domain.max);
  } else {
    min_t0_ptr = &(scope[0].get_var().range_domain->min);
    max_t0_ptr = &(scope[0].get_var().range_domain->max);
    //std::cerr << "TODO " << scope[0].domain_type << std::endl;
    //exit(1);
  }

  if(scope[1].domain_type == RANGE_VAR) {
    min_t1_ptr = &(scope[1].range_domain->min);
    max_t1_ptr = &(scope[1].range_domain->max);
  } else if(scope[1].domain_type == BITSET_VAR) {
    min_t1_ptr = &(scope[1].bitset_domain->domain.min);
    max_t1_ptr = &(scope[1].bitset_domain->domain.max);
  } else {
    min_t1_ptr = &(scope[1].get_var().range_domain->min);
    max_t1_ptr = &(scope[1].get_var().range_domain->max);
    //std::cerr << "TODO " << scope[1].domain_type << std::endl;
    //exit(1);
  }

}

Mistral::PropagationOutcome Mistral::ConstraintReifiedDisjunctive::propagate() {
  return CONSISTENT;
}



Mistral::PropagationOutcome Mistral::ConstraintReifiedDisjunctive::propagate(const int changed_idx, const Event evt) {
  PropagationOutcome wiped = CONSISTENT;


#ifdef _DEBUG_RDISJUNCTIVE
  std::cout << std::endl << "propagate [" << id << "] " << this << " " << event2str(evt) 
	    << " on " << scope[changed_idx] << " in " << scope[changed_idx].get_domain() << std::endl;
  for(unsigned int i=0; i<3; ++i) {
    std::cout << scope[i] << " in " << scope[i].get_domain() << " ";
  }
  std::cout << std::endl;
#endif


  if(changed_idx == 2) {

#ifdef _DEBUG_RDISJUNCTIVE
    std::cout << "  -> now it's a precedence" << std::endl;
#endif
    
    // it becomes a precedence, we must enforce both rules no matter what and there cannot be further changes.
    
    if(LB_CHANGED(evt)) {
      if( scope[0].set_max( *max_t1_ptr-processing_time[0] ) == FAIL_EVENT) wiped = FAILURE(0);
      else if( scope[1].set_min( *min_t0_ptr+processing_time[0] ) == FAIL_EVENT) wiped = FAILURE(1);
    } else {
      if( scope[0].set_min( *min_t1_ptr+processing_time[1] ) == FAIL_EVENT) wiped = FAILURE(0);
      else if( scope[1].set_max( *max_t0_ptr-processing_time[1] ) == FAIL_EVENT) wiped = FAILURE(1);
    }
  } else if(*state!=3) {
    

#ifdef _DEBUG_RDISJUNCTIVE
    std::cout << "  -> it was already a precedence" << std::endl;
#endif

    // it was already a precedence, we enforce the changes

    if(*state==2) {
      if(// events.contain(1) 
	 changed_idx == 1
	 && UB_CHANGED(evt)) {
	if( scope[0].set_max( *max_t1_ptr-processing_time[0] ) == FAIL_EVENT) wiped = FAILURE(0);
      }
      if(IS_OK(wiped) && // events.contain(0)
	 changed_idx == 0
	 && LB_CHANGED(evt)) {
	if( scope[1].set_min( *min_t0_ptr+processing_time[0] ) == FAIL_EVENT) wiped = FAILURE(1);
      }
    } else {
      if(// events.contain(0)
	 changed_idx == 0
	 && UB_CHANGED(evt)) {
	if( scope[1].set_max( *max_t0_ptr-processing_time[1] ) == FAIL_EVENT) wiped = FAILURE(1);
      }
      if(IS_OK(wiped) && // events.contain(1)
	 changed_idx == 1
	 && LB_CHANGED(evt)) {
	if( scope[0].set_min( *min_t1_ptr+processing_time[1] ) == FAIL_EVENT) wiped = FAILURE(0);
      }
    } 
  } else {

#ifdef _DEBUG_RDISJUNCTIVE
    std::cout << "  -> could it be a precedence? (" 
	      << *min_t0_ptr << " + " << processing_time[0] 
	      << " > " <<  *max_t1_ptr << " || "
	      << *min_t1_ptr << " + " << processing_time[1] 
	      << " > " << *max_t0_ptr << ")" << std::endl;
#endif
    
    // the disjunction is not yet decided, we check if it should

    if( *min_t0_ptr + processing_time[0] > *max_t1_ptr ) {

      if( scope[2].set_domain(0) == FAIL_EVENT) wiped = FAILURE(2);
      // x[1]+p[1] <= x[0] because x[0]'s min is too high or x[1]'s max is too low
      else if( scope[0].set_min( *min_t1_ptr+processing_time[1] ) == FAIL_EVENT) wiped = FAILURE(0);
      else if( scope[1].set_max( *max_t0_ptr-processing_time[1] ) == FAIL_EVENT) wiped = FAILURE(1);

#ifdef _DEBUG_RDISJUNCTIVE
      std::cout << "  -> YES!: " ;
      display(std::cout);
      std::cout << std::endl;
#endif

    } else if( *min_t1_ptr + processing_time[1] > *max_t0_ptr ) {

      if( scope[2].set_domain(1) == FAIL_EVENT) wiped = FAILURE(2);
      else if( scope[0].set_max( *max_t1_ptr-processing_time[0] ) == FAIL_EVENT) wiped = FAILURE(0);
      else if( scope[1].set_min( *min_t0_ptr+processing_time[0] ) == FAIL_EVENT) wiped = FAILURE(1);
    
#ifdef _DEBUG_RDISJUNCTIVE
      std::cout << "  -> YES!: " ;
      display(std::cout);
      std::cout << std::endl;
#endif

    } 
#ifdef _DEBUG_RDISJUNCTIVE
    else 
      std::cout << "  -> NO!" << std::endl;
#endif

  }

#ifdef _DEBUG_RDISJUNCTIVE
  for(unsigned int i=0; i<3; ++i) {
    std::cout << scope[i] << " in " << scope[i].get_domain() << " ";
  }
  std::cout << std::endl;
  std::cout << "+-> " << (IS_OK(wiped) ? "OK" : "inconsistency!") << std::endl;
#endif
  
  //changes.clear();

  return wiped;
}

std::ostream& Mistral::ConstraintReifiedDisjunctive::display(std::ostream& os) const {
  if(scope[2].is_ground()) {
    if(!scope[2].get_max()) {
      os << scope[1] << " + " << processing_time[1] << " <= " << scope[0]/*.get_var()*/ ;
    } else {
      os << scope[0] << " + " << processing_time[0] << " <= " << scope[1]/*.get_var()*/ ;
    }
  } else {
    os << scope[0]/*.get_var()*/ << " + " << processing_time[0] << " <= " << scope[1]/*.get_var()*/ << " or " 
       << scope[1]/*.get_var()*/ << " + " << processing_time[1] << " <= " << scope[0]/*.get_var()*/ ;
  }
  return os;
}


void Mistral::PredicateEqual::initialise() {
    ConstraintImplementation::initialise();

  trigger_on(_DOMAIN_, scope[0]);
  trigger_on(_DOMAIN_, scope[1]);
  trigger_on(_VALUE_, scope[2]);

  //ConstraintImplementation::initialise();
  //set_idempotent(false);
}

void Mistral::PredicateEqual::mark_domain() {
  if(!spin) {
    get_solver()->mark_non_convex(scope[0].id());
    get_solver()->mark_non_convex(scope[1].id());
  }
}

Mistral::PropagationOutcome Mistral::PredicateEqual::rewrite() {

  // if(id==49) {
  //   std::cout << is_posted << " rewrite (" << scope[0] << " in " << scope[0].get_domain()
  // 	      << " == " << scope[1] << " in " << scope[1].get_domain()
  // 	      << ") <-> " << scope[2] << " in " << scope[2].get_domain() << std::endl;
  // }

  Mistral::PropagationOutcome wiped = propagate();

#ifdef _COMPLETE
  if( active.size == 2 ) {
    relax();

    if( scope[2].is_ground() ) {
      if( (spin + scope[2].get_min()) == 1 ) {
	solver->add(new ConstraintNotEqual(scope));
      } else {

	// std::cout << "   => change to an equal constraint" << std::endl;

	--scope.size;
	solver->add(new ConstraintEqual(scope));
      }
    } else {
      Vector< Variable > new_scope;

      int val = NOVAL, ground_idx = scope[1].is_ground();
      new_scope.add(scope[1-ground_idx]);
      new_scope.add(scope[2]);
      
      val = scope[ground_idx].get_min();
      solver->add(new PredicateConstantEqual(new_scope, val, spin));
    }
  }
#endif

  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateEqual::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[2].is_ground() ) {
    if( (spin + scope[2].get_min()) == 1 ) {
      if(scope[0].is_ground() && IS_FAIL(scope[1].remove(scope[0].get_min())))
	wiped = FAILURE(1);
      else {
	if(scope[1].is_ground() && IS_FAIL(scope[0].remove(scope[1].get_min())))
	  wiped = FAILURE(0);
      }
    } else {
      if( IS_FAIL(scope[0].set_domain(scope[1])) ) wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].set_domain(scope[0])) ) wiped = FAILURE(1);
    }
  } else {
    if( !(scope[0].intersect(scope[1])) ) {
      if( IS_FAIL(scope[2].remove(spin)) ) wiped = FAILURE(2);	    
    } else { 
      if( scope[0].is_ground() && scope[1].is_ground() ) {
	if( IS_FAIL(scope[2].set_domain(spin) )) wiped = FAILURE(2);
      }
    }
  }
  
  return wiped;
}


Mistral::PropagationOutcome Mistral::PredicateEqual::propagate(const int changed_idx, const Event evt) {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( changed_idx == 2 ) {
    if( LB_CHANGED(evt) == spin ) {
      // x[0] == x[1]
      if( IS_FAIL(scope[0].set_domain(scope[1])) ) wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].set_domain(scope[0])) ) wiped = FAILURE(1);
    } else {
      // x[0] != x[1]
      if(scope[0].is_ground() && IS_FAIL(scope[1].remove(scope[0].get_min()))) wiped = FAILURE(1);
      else if(scope[1].is_ground() && IS_FAIL(scope[0].remove(scope[1].get_min()))) wiped = FAILURE(0);
    }
  } else {
    if( scope[2].is_ground() ) {

      // std::cout << "HERE: " << spin << " " <<  scope[2].get_min() << " " << scope[changed_idx].get_domain() << std::endl;

      if( (spin != scope[2].get_min()) ) {
	if(scope[changed_idx].is_ground() && IS_FAIL(scope[1-changed_idx].remove(scope[changed_idx].get_min()))) 
	  wiped = FAILURE(1-changed_idx);
      } else {
	if( IS_FAIL(scope[1-changed_idx].set_domain(scope[changed_idx])) ) 
	  wiped = FAILURE(1-changed_idx);
      }
    } else {
      // check if (in)equality can be deduced
      if( !(scope[0].intersect(scope[1])) ) {
	if( IS_FAIL(scope[2].remove(spin)) ) wiped = FAILURE(2);	    
      } else { 
	if( scope[0].is_ground() && scope[1].is_ground() ) {
	  if( IS_FAIL(scope[2].set_domain(spin) )) wiped = FAILURE(2);
	}
      }
    }
  }
  
  return wiped;
}

std::ostream& Mistral::PredicateEqual::display(std::ostream& os) const {
  os << scope[2]/*.get_var()*/ << " <=> (";
  if(spin) os << scope[0]/*.get_var()*/ << " == " << scope[1]/*.get_var()*/;
  else os << scope[0]/*.get_var()*/ << " =/= " << scope[1]/*.get_var()*/;
  os << ")";
  return os;
}


void Mistral::PredicateConstantEqual::initialise() {
  ConstraintImplementation::initialise();
  trigger_on(_DOMAIN_, scope[0]);
  trigger_on(_VALUE_, scope[1]);
}

void Mistral::PredicateConstantEqual::mark_domain() {
  get_solver()->mark_non_convex(scope[0].id());
}

Mistral::PropagationOutcome Mistral::PredicateConstantEqual::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[1].is_ground() ) {
    if( (spin + scope[1].get_min()) == 1 ) {
      if( IS_FAIL(scope[0].remove(value)) ) wiped = FAILURE(0);
    } else {
      if( IS_FAIL(scope[0].set_domain(value)) ) wiped = FAILURE(0);
    }
  } else {
    if( !(scope[0].contain(value)) ) {
      if( IS_FAIL(scope[1].remove(spin)) ) wiped = FAILURE(1);	    
    } else if( scope[0].is_ground() ) {
      if( IS_FAIL(scope[1].set_domain(spin)) ) wiped = FAILURE(1);
    }
  }
  
  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateConstantEqual::propagate(const int changed_idx, const Event evt) {
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( changed_idx ) {
    if( (spin + LB_CHANGED(evt)) == 1 ) {
      if( IS_FAIL(scope[0].remove(value)) ) wiped = FAILURE(0);
    } else {
      if( IS_FAIL(scope[0].set_domain(value)) ) wiped = FAILURE(0);
    }
  } else {
    if( !(scope[0].contain(value)) ) {
      if( IS_FAIL(scope[1].remove(spin)) ) wiped = FAILURE(1);	    
    } else if( scope[0].is_ground() ) {
      if( IS_FAIL(scope[1].set_domain(spin)) ) wiped = FAILURE(1);
    }
  }

  return wiped;
}

std::ostream& Mistral::PredicateConstantEqual::display(std::ostream& os) const {
  os << scope[1]/*.get_var()*/ << " <=> (";
  if(spin) os << scope[0]/*.get_var()*/ << " == " << value;
  else os << scope[0]/*.get_var()*/ << " =/= " << value;
  os << ")";
  return os;
}


void Mistral::PredicateIntervalMember::initialise() {
  ConstraintImplementation::initialise();
  trigger_on(_DOMAIN_, scope[0]);
  trigger_on(_VALUE_, scope[1]);
}

void Mistral::PredicateIntervalMember::mark_domain() {
  get_solver()->mark_non_convex(scope[0].id());
}

Mistral::PropagationOutcome Mistral::PredicateIntervalMember::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[1].is_ground() ) {
    if( (spin + scope[1].get_min()) == 1 ) {
      if( IS_FAIL(scope[0].remove_interval(lower_bound, upper_bound)) ) wiped = FAILURE(0);
    } else {
      //if( IS_FAIL(scope[0].set_domain(lower_bound, upper_bound)) ) wiped = FAILURE(0);
      if( IS_FAIL(scope[0].set_min(lower_bound)) ) wiped = FAILURE(0);
      else if( IS_FAIL(scope[0].set_max(upper_bound)) ) wiped = FAILURE(0);
    }
  } else {

    //std::cout << scope[0].get_domain() << " <> [" << lower_bound << ".." << upper_bound << "]" << std::endl; 

    if( !(scope[0].intersect(lower_bound, upper_bound)) ) {
      if( IS_FAIL(scope[1].remove(spin)) ) wiped = FAILURE(1);	    
    } else if( scope[0].included(lower_bound, upper_bound) ) {
      if( IS_FAIL(scope[1].set_domain(spin)) ) wiped = FAILURE(1);
    }
  }
  
  //exit(1);

  return wiped;
}


Mistral::PropagationOutcome Mistral::PredicateIntervalMember::propagate(const int changed_idx, const Event evt) {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if(changed_idx == 1) {
    if( (spin + LB_CHANGED(evt)) == 1 ) {
      if( IS_FAIL(scope[0].remove_interval(lower_bound, upper_bound)) ) wiped = FAILURE(0);
    } else {
      if( IS_FAIL(scope[0].set_min(lower_bound)) ) wiped = FAILURE(0);
      else if( IS_FAIL(scope[0].set_max(upper_bound)) ) wiped = FAILURE(0);
    }
  } else {
    if( !(scope[0].intersect(lower_bound, upper_bound)) ) {
      if( IS_FAIL(scope[1].remove(spin)) ) wiped = FAILURE(1);	    
    } else if( scope[0].included(lower_bound, upper_bound) ) {
      if( IS_FAIL(scope[1].set_domain(spin)) ) wiped = FAILURE(1);
    }
  }

  return wiped;
}

std::ostream& Mistral::PredicateIntervalMember::display(std::ostream& os) const {
  os << scope[1]/*.get_var()*/ << " <=> (";
  if(spin) os << scope[0]/*.get_var()*/ << " in [" << lower_bound << ".." << upper_bound << "]";
  else os << scope[0]/*.get_var()*/ << " in ..]" << lower_bound << "," << upper_bound << "[..";
  os << ")";
  return os;
}


void Mistral::PredicateSetMember::initialise() {
  ConstraintImplementation::initialise();
  trigger_on(_DOMAIN_, scope[0]);
  trigger_on(_VALUE_, scope[1]);
  //set_idempotent(false);

  non_values.initialise(scope[0].get_min(), scope[0].get_max(), BitSet::full);
  non_values.setminus_with(values);
}

void Mistral::PredicateSetMember::mark_domain() {
  get_solver()->mark_non_convex(scope[0].id());
}

Mistral::PropagationOutcome Mistral::PredicateSetMember::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[1].is_ground() ) {
    if( (spin + scope[1].get_min()) == 1 ) {
      if( IS_FAIL(scope[0].set_domain(non_values)) ) wiped = FAILURE(0);
    } else {
      if( IS_FAIL(scope[0].set_domain(values)) ) wiped = FAILURE(0);
    }
  } else {
    if( !(scope[0].intersect(values)) ) {
      if( IS_FAIL(scope[1].remove(spin)) ) wiped = FAILURE(1);	    
    } else if( scope[0].included(values) ) {
      if( IS_FAIL(scope[1].set_domain(spin)) ) wiped = FAILURE(1);
    }
  }
  
  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateSetMember::propagate(const int changed_idx, const Event evt) {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( changed_idx == 1 ) {
    if( (spin + LB_CHANGED(evt)) == 1 ) {
      if( IS_FAIL(scope[0].set_domain(non_values)) ) wiped = FAILURE(0);
    } else {
      if( IS_FAIL(scope[0].set_domain(values)) ) wiped = FAILURE(0);
    }
  } else {
    if( !(scope[0].intersect(values)) ) {
      if( IS_FAIL(scope[1].remove(spin)) ) wiped = FAILURE(1);	    
    } else if( scope[0].included(values) ) {
      if( IS_FAIL(scope[1].set_domain(spin)) ) wiped = FAILURE(1);
    }
  }
  
  return wiped;
}

std::ostream& Mistral::PredicateSetMember::display(std::ostream& os) const {
  os << scope[1]/*.get_var()*/ << " <=> (" << scope[0]/*.get_var()*/ ;
  if(!spin) os << " not";
  os << " in " << values; 
  os << ")";
  return os;
}


void Mistral::PredicateOffset::initialise() {
  ConstraintImplementation::initialise();
  trigger_on(_RANGE_, scope[0]);
  trigger_on(_RANGE_, scope[1]);
}

Mistral::PropagationOutcome Mistral::PredicateOffset::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( IS_FAIL(scope[0].set_min( scope[1].get_min() - offset )) ) wiped = FAILURE(0); 
  else if( IS_FAIL(scope[0].set_max( scope[1].get_max() - offset )) ) wiped = FAILURE(0); 
  else if( IS_FAIL(scope[1].set_min( scope[0].get_min() + offset )) ) wiped = FAILURE(1); 
  else if( IS_FAIL(scope[1].set_max( scope[0].get_max() + offset )) ) wiped = FAILURE(1);
  
  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateOffset::propagate(const int changed_idx, const Event evt) {
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if(changed_idx == 1) {
    if( LB_CHANGED(evt) && IS_FAIL(scope[0].set_min( scope[1].get_min() - offset )) ) wiped = FAILURE(0); 
    else if( UB_CHANGED(evt) && IS_FAIL(scope[0].set_max( scope[1].get_max() - offset )) ) wiped = FAILURE(0);
  } else { 
    if( LB_CHANGED(evt) && IS_FAIL(scope[1].set_min( scope[0].get_min() + offset )) ) wiped = FAILURE(1); 
    else if( UB_CHANGED(evt) && IS_FAIL(scope[1].set_max( scope[0].get_max() + offset )) ) wiped = FAILURE(1);
  }

  return wiped;
}

std::ostream& Mistral::PredicateOffset::display(std::ostream& os) const {
  os << scope[1]/*.get_var()*/ << " == (" << scope[0]/*.get_var()*/ << " + " << offset << ")";
  return os;
}


void Mistral::PredicateFactor::initialise() {
  ConstraintImplementation::initialise();
  trigger_on(_RANGE_, scope[0]);
  trigger_on(_RANGE_, scope[1]);
}

Mistral::PropagationOutcome Mistral::PredicateFactor::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( IS_FAIL(scope[0].set_min( (factor>0 ? 
				 scope[1].get_min() : 
				 scope[1].get_max()) / factor )) ) {
    wiped = FAILURE(0); 

  } else if( IS_FAIL(scope[0].set_max( (factor>0 ? 
				      scope[1].get_max() :
					scope[1].get_min()) / factor )) ) {
    wiped = FAILURE(0); 

  } else if( IS_FAIL(scope[1].set_min( (factor>0 ? 
					scope[0].get_min() : 
					scope[0].get_max()) * factor )) ) {
    wiped = FAILURE(1); 

  } else if( IS_FAIL(scope[1].set_max( (factor>0 ? 
					scope[0].get_max() :
					scope[0].get_min()) * factor )) ) {
    wiped = FAILURE(1);

  }

  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateFactor::propagate(const int changed_idx, const Event evt) {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if(changed_idx == 1) {
    if(factor>0) {
      if( LB_CHANGED(evt) && IS_FAIL(scope[0].set_min( scope[1].get_min() / factor )) )
	wiped = FAILURE(0); 
      if( UB_CHANGED(evt) && IS_FAIL(scope[0].set_max( scope[1].get_max() / factor )) ) 
	wiped = FAILURE(0); 
    } else {
      if( UB_CHANGED(evt) && IS_FAIL(scope[0].set_min( scope[1].get_max() / factor )) )
	wiped = FAILURE(0); 
      if( LB_CHANGED(evt) && IS_FAIL(scope[0].set_max( scope[1].get_min() / factor )) )
	wiped = FAILURE(0); 
    }
  } else {
    if(factor>0) {
      if( LB_CHANGED(evt) && IS_FAIL(scope[1].set_min( scope[0].get_min() * factor )) )
	wiped = FAILURE(1); 
      if( UB_CHANGED(evt) && IS_FAIL(scope[1].set_max( scope[0].get_max() * factor )) ) 
	wiped = FAILURE(1); 
    } else {
      if( UB_CHANGED(evt) && IS_FAIL(scope[1].set_min( scope[0].get_max() * factor )) )
	wiped = FAILURE(1); 
      if( LB_CHANGED(evt) && IS_FAIL(scope[1].set_max( scope[0].get_min() * factor )) )
	wiped = FAILURE(1); 
    }
  }
  
  return wiped;
}

std::ostream& Mistral::PredicateFactor::display(std::ostream& os) const {
  os << scope[1]/*.get_var()*/ << " == (" << scope[0]/*.get_var()*/ << " * " << factor << ")";
  return os;
}


void Mistral::PredicateNot::initialise() {
  ConstraintImplementation::initialise();
  trigger_on(_VALUE_, scope[0]);
  trigger_on(_VALUE_, scope[1]);
}


Mistral::PropagationOutcome Mistral::PredicateNot::rewrite() {
  RewritingOutcome r_evt = NO_EVENT; 

  if(scope[0].is_ground() || scope[1].is_ground() ) {

#ifdef _DEBUG_REWRITE
    std::cout << "    relax ";
    display(std::cout);
    std::cout << std::endl;
#endif

    r_evt = SUPPRESSED;
    relax();
  } else {
    
    // find out if the constraints on x[0] can be negated

    bool hold = true;
    for(int var=0; hold && var<2; ++var) {

      int idx = scope[var].id();
      bool can_negate = true;
      Constraint con;
      for(Event trig = 0; can_negate && trig<3; ++trig) {
	for(int i = get_solver()->constraint_graph[idx].on[trig].size; can_negate && i--;) {
	  con = get_solver()->constraint_graph[idx].on[trig][i];
	  can_negate &= con.absorb_negation(con.index());
	}
      }
      
      if(can_negate) {
	
	hold = false;

	r_evt = SUPPRESSED;
	
#ifdef _DEBUG_REWRITE
	std::cout << "    relax ";
	display(std::cout);
	std::cout << " and remove " << scope[var] << std::endl;
#endif
	
	relax();
	get_solver()->remove(scope[var]);
	
	for(Event trig = 0; trig<3; ++trig) {
	  for(int i = get_solver()->constraint_graph[idx].on[trig].size; i--;) {
	    con = get_solver()->constraint_graph[idx].on[trig][i];
	    
#ifdef _DEBUG_REWRITE
	    std::cout << "      relax " << con << std::endl;
#endif
	    
	    con.relax();
	    con = Constraint(con.get_negation(con.index(), scope[1-var]));

#ifdef _DEBUG_REWRITE
	    std::cout << "      post " << con << std::endl;
#endif

	    get_solver()->add( con );
	  }
	}
      }
    }
  }  

  return r_evt;
}


Mistral::PropagationOutcome Mistral::PredicateNot::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[1].is_ground() ) {
    if( scope[1].get_min() ) { // x[0] == 0
      if( IS_FAIL(scope[0].set_domain(0)) ) wiped = FAILURE(0);
    } else if( scope[1].get_max() == 0 ) { //x[0] != 0
      if( IS_FAIL(scope[0].remove(0)) ) wiped = FAILURE(0);
    } 
  } else {
    if( !scope[0].contain(0) ) {
      if(IS_FAIL(scope[1].set_domain(0))) wiped = FAILURE(1);
    } else if( scope[0].is_ground() ) {
      if(IS_FAIL(scope[1].remove(0))) wiped = FAILURE(1);
    }
  }
  
  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateNot::propagate(const int changed_idx, const Event evt) {
  Mistral::PropagationOutcome wiped = CONSISTENT;

  //  std::cout << std::endl << scope[0].get_domain() << " " << scope[1].get_domain() << std::endl; 

  // if( changed_idx ) {
  //   if( LB_CHANGED(evt) ) { // x[0] == 0
  //     if( IS_FAIL(scope[0].set_domain(0)) ) wiped = FAILURE(0);
  //   } else { //x[0] != 0
  //     if( IS_FAIL(scope[0].remove(0)) ) wiped = FAILURE(0);
  //   } 
  // } else {
  //   if( LB_CHANGED(evt) ) {
  //     if(IS_FAIL(scope[1].set_domain(0))) wiped = FAILURE(1);
  //   } else {
  //     if(IS_FAIL(scope[1].remove(0))) wiped = FAILURE(1);
  //   }
  // }

  if( IS_FAIL( scope[1-changed_idx].set_domain(UB_CHANGED(evt)) ) ) wiped = FAILURE(1-changed_idx);


  //  std::cout << scope[0].get_domain() << " " << scope[1].get_domain() << std::endl << std::endl ; 
  
  return wiped;
}

std::ostream& Mistral::PredicateNot::display(std::ostream& os) const {
  os << scope[1]/*.get_var()*/ << " <=> not(" << scope[0]/*.get_var()*/ << ")";
  return os;
}


void Mistral::PredicateNeg::initialise() {
  ConstraintImplementation::initialise();
  trigger_on(_DOMAIN_, scope[0]);
  trigger_on(_DOMAIN_, scope[1]);
}

Mistral::PropagationOutcome Mistral::PredicateNeg::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  for(int changed_idx = 0; changed_idx < 2; ++changed_idx) {
    if(IS_FAIL(scope[1-changed_idx].set_max(-scope[changed_idx].get_min()))) wiped = FAILURE(1-changed_idx);
    if(IS_FAIL(scope[1-changed_idx].set_min(-scope[changed_idx].get_max()))) wiped = FAILURE(1-changed_idx);
    
    if( scope[changed_idx].get_size() < scope[1-changed_idx].get_size() ) {
      int vali, vnxt = scope[1-changed_idx].get_min();
      do {
	vali = vnxt;
	if(!scope[changed_idx].contain(-vali) && IS_FAIL(scope[1-changed_idx].remove(vali))) wiped = FAILURE(1-changed_idx);
	vnxt = scope[1-changed_idx].next(vali);
      } while(vali < vnxt);
    }
  }

  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateNeg::propagate(const int changed_idx, const Event evt) {
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( ASSIGNED(evt) ) {
    if(IS_FAIL(scope[1-changed_idx].set_domain(-scope[changed_idx].get_min()))) wiped = FAILURE(1-changed_idx);
  } else {
    if( LB_CHANGED(evt) ) {
      if(IS_FAIL(scope[1-changed_idx].set_max(-scope[changed_idx].get_min()))) wiped = FAILURE(1-changed_idx);
    }
    if( UB_CHANGED(evt) ) {
      if(IS_FAIL(scope[1-changed_idx].set_min(-scope[changed_idx].get_max()))) wiped = FAILURE(1-changed_idx);
    }
    if( scope[changed_idx].get_size() < scope[1-changed_idx].get_size() ) {
      int vali, vnxt = scope[1-changed_idx].get_min();
      do {
	vali = vnxt;
	if(!scope[changed_idx].contain(-vali) && IS_FAIL(scope[1-changed_idx].remove(vali))) wiped = FAILURE(1-changed_idx);
	vnxt = scope[1-changed_idx].next(vali);
      } while(vali < vnxt);
    }
  }
  
  return wiped;
}

std::ostream& Mistral::PredicateNeg::display(std::ostream& os) const {
  os << scope[1]/*.get_var()*/ << " = -(" << scope[0]/*.get_var()*/ << ")";
  return os;
}


void Mistral::PredicateAnd::initialise() {
  ConstraintImplementation::initialise();
  for(int i=0; i<3; ++i)
    trigger_on(_VALUE_, scope[i]);
}

Mistral::PropagationOutcome Mistral::PredicateAnd::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  // std::cout << std::endl << "PROPAGATE " << this << std::endl;
  // std::cout << scope[0].get_domain() << " & " 
  // 	    << scope[1].get_domain() << " <-> "
  // 	    << scope[2].get_domain() << std::endl; 


  if( scope[2].is_ground() ) {
    if( scope[2].get_min() ) {
      if( IS_FAIL(scope[0].remove(0)) ) 
	wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].remove(0)) ) 
	wiped = FAILURE(1);
    } else if( scope[2].get_max() == 0 ) {
      if( scope[1].get_min() ){
	if( IS_FAIL(scope[0].set_domain(0)) ) wiped = FAILURE(0);
      } else if( scope[0].get_min() ) {
	if( IS_FAIL(scope[1].set_domain(0)) ) wiped = FAILURE(1);
      }
    } 
  } else {
    if( scope[0].get_min() && scope[1].get_min() ) {
      if( IS_FAIL(scope[2].remove(0)) ) return FAILURE(2);
    } else if( !scope[0].get_max() || !scope[1].get_max() ) {
      if( IS_FAIL(scope[2].set_domain(0)) ) return FAILURE(2);
    }
  }

 // std::cout << scope[0].get_domain() << " & " 
 // 	    << scope[1].get_domain() << " <-> "
 // 	   << scope[2].get_domain() << std::endl << wiped << std::endl; 
  
  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateAnd::propagate(const int changed_idx, 
							     const Event evt) {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  // std::cout << std::endl << "PROPAGATE " << this << std::endl;
  // std::cout << scope[0].get_domain() << " & " 
  // 	    << scope[1].get_domain() << " <-> "
  // 	    << scope[2].get_domain() << std::endl; 

  if(changed_idx == 2) {
    if(LB_CHANGED(evt)) {
      if( IS_FAIL(scope[0].remove(0)) ) wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].remove(0)) ) wiped = FAILURE(1);   
    } else {
      if( scope[1].get_min() ){
	if( IS_FAIL(scope[0].set_domain(0)) ) wiped = FAILURE(0);
      } else if( scope[0].get_min() ) {
	if( IS_FAIL(scope[1].set_domain(0)) ) wiped = FAILURE(1);
      }
    }
  } else { // either z is not yet set, or it is a not(x and y) constraint
    if( scope[2].is_ground() ) {
      if(LB_CHANGED(evt)) {
	// it is an "not(AND)" and one of the variable was set to 1
	if(IS_FAIL(scope[1-changed_idx].set_domain(0))) wiped = FAILURE(1-changed_idx);
      } 
    } else {
      if(UB_CHANGED(evt)) {
	if(IS_FAIL(scope[2].set_domain(0))) wiped = FAILURE(2);
      } else if(scope[1-changed_idx].is_ground()) {
	if(IS_FAIL(scope[2].remove(0))) wiped = FAILURE(2);
      }
    }
  }

 // std::cout << scope[0].get_domain() << " & " 
 // 	    << scope[1].get_domain() << " <-> "
 // 	   << scope[2].get_domain() << std::endl << wiped << std::endl; 
  
  return wiped;
}

std::ostream& Mistral::PredicateAnd::display(std::ostream& os) const {
  os << scope[2]/*.get_var()*/ << " <=> (" << scope[0]/*.get_var()*/ << " and " << scope[1]/*.get_var()*/ << ")";
  return os;
}


void Mistral::PredicateOr::initialise() {
  ConstraintImplementation::initialise();
  for(int i=0; i<3; ++i)
    trigger_on(_VALUE_, scope[i]);
}

// Mistral::PropagationOutcome Mistral::PredicateOr::rewrite() {
//   Mistral::PropagationOutcome wiped = propagate();
//   if( scope[2].is_ground() && active.size == 2 ) {
//     relax();
//     if( scope[2].get_min() ) {
//       solver->add(new ConstraintOr(scope));
//     }
//   }
//   return wiped;
// }


Mistral::PropagationOutcome Mistral::PredicateOr::propagate(const int changed_idx, 
							    const Event evt) {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if(changed_idx == 2) {
    if(UB_CHANGED(evt)) {
      if( IS_FAIL(scope[0].set_domain(0)) ) wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].set_domain(0)) ) wiped = FAILURE(1);
    } else {
      if( scope[1].get_max() == 0 ){
	if( IS_FAIL(scope[0].remove(0)) ) wiped = FAILURE(0);
      } else if( scope[0].get_max() == 0 ) {
	if( IS_FAIL(scope[1].remove(0)) ) wiped = FAILURE(1);
      }
    } 
  } else { // either z is not yet set, or it is a not(x and y) constraint
    if( scope[2].is_ground() ) {
      if(UB_CHANGED(evt)) {
	// it is an "OR" constraint
	if(IS_FAIL(scope[1-changed_idx].remove(0))) wiped = FAILURE(1-changed_idx);
      }
    } else { 
      if(LB_CHANGED(evt)) {
	if( IS_FAIL(scope[2].remove(0)) ) wiped = FAILURE(2);
      } else if(scope[1-changed_idx].is_ground()) {
	if( IS_FAIL(scope[2].set_domain(0)) ) wiped = FAILURE(2);
      }
    }
  }

  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateOr::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[2].is_ground() ) {
    if( scope[2].get_max() == 0 ) {
      if( IS_FAIL(scope[0].set_domain(0)) ) 
	wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].set_domain(0)) ) 
	wiped = FAILURE(1);
    } else if( scope[2].get_min() ) {
      if( scope[1].get_max() == 0 ){
	if( IS_FAIL(scope[0].remove(0)) ) wiped = FAILURE(0);
      } else if( scope[0].get_max() == 0 ) {
	if( IS_FAIL(scope[1].remove(0)) ) wiped = FAILURE(1);
      }
    } 
  } else {
    if( scope[0].get_min() || scope[1].get_min() ) {
      if( IS_FAIL(scope[2].remove(0)) ) wiped = FAILURE(2);
    } else if( !scope[0].get_max() && !scope[1].get_max() ) {
      if( IS_FAIL(scope[2].set_domain(0)) ) wiped = FAILURE(2);
    }
  }

  return wiped;
}

std::ostream& Mistral::PredicateOr::display(std::ostream& os) const {
  os << scope[2]/*.get_var()*/ << " <=> (" << scope[0]/*.get_var()*/ << " or " << scope[1]/*.get_var()*/ << ")";
  return os;
}


void Mistral::ConstraintAnd::initialise() {
  ConstraintImplementation::initialise();
  
  trigger_on(_DOMAIN_, scope[0]);
  trigger_on(_DOMAIN_, scope[1]);
}

Mistral::PropagationOutcome Mistral::ConstraintAnd::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( IS_FAIL(scope[0].remove(0)) ) 
    wiped = FAILURE(0);
  else if( IS_FAIL(scope[1].remove(0)) ) 
    wiped = FAILURE(1);
  
  return wiped;
}

Mistral::PropagationOutcome Mistral::ConstraintAnd::propagate(const int changed_idx, 
							      const Event evt) {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( IS_FAIL(scope[0].remove(0)) ) 
    wiped = FAILURE(0);
  else if( IS_FAIL(scope[1].remove(0)) ) 
    wiped = FAILURE(1);
  
  return wiped;
}

std::ostream& Mistral::ConstraintAnd::display(std::ostream& os) const {
  os << "(" << scope[0]/*.get_var()*/ << " and " << scope[1]/*.get_var()*/ << ")";
  return os;
}

void Mistral::ConstraintOr::initialise() {
  ConstraintImplementation::initialise();
  trigger_on(_VALUE_, scope[0]);
  trigger_on(_VALUE_, scope[1]);
}

Mistral::PropagationOutcome Mistral::ConstraintOr::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[1].get_max() == 0 ){
    if( IS_FAIL(scope[0].remove(0)) ) wiped = FAILURE(0);
  } else if( scope[0].get_max() == 0 ) {
    if( IS_FAIL(scope[1].remove(0)) ) wiped = FAILURE(1);
  }
  
  return wiped;
}

Mistral::PropagationOutcome Mistral::ConstraintOr::propagate(const int changed_idx, 
							     const Event evt) {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if(UB_CHANGED(evt) && IS_FAIL(scope[1-changed_idx].remove(0))) wiped = FAILURE(1-changed_idx);
  
  return wiped;
}

std::ostream& Mistral::ConstraintOr::display(std::ostream& os) const {
  os << "(" << scope[0]/*.get_var()*/ << " or " << scope[1]/*.get_var()*/ << ")";
  return os;
}


void Mistral::ConstraintLex::initialise() {
  ConstraintImplementation::initialise();
  for(int i=0; i<4; ++i)
    trigger_on(_RANGE_, scope[i]);
  GlobalConstraint::initialise();
}




Mistral::PropagationOutcome Mistral::ConstraintLex::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;
  Event evt, aux;

  /*
    rules:
    
    1: IF b1=0 THEN x0 = x1     ACTIVATION: change on b1 or state&1 and change on x0 or x1
    2: IF x0 =/= x1 THEN b1=1   ACTIVATION: change on x0 or x1

    4: IF b0<b1 THEN x0<x1      ACTIVATION: change on b0 or b1 or state&4 and change on x0 or x1
    8: IF x0>=x1 THEN b0=b1     ACTIVATION: change on x0 or x1 or state&8 and change on b0 or b1

    16: IF b0=1 THEN b1=1       ACTIVATION: change on b0 or b1
    32: IF b1=0 THEN b0=0       ACTIVATION: change on b0 or b1
   */


  int state = (!scope[3].get_max() + 
	       4*(scope[2].get_max() < scope[3].get_min()) + 
	       8*(scope[0].get_min()>=scope[1].get_max()));


#ifdef _DEBUG_LEX
  std::cout << std::endl;
#endif


  while(IS_OK(wiped) && !changes.empty()) {
    evt = changes.pop();
    
#ifdef _DEBUG_LEX
    std::cout << " change on " << scope[evt] << " in " << scope[evt].get_domain() 
	      << " state: ( " ;
    if(state&1) std::cout << scope[3] << "=1 ";
    if(state&4) std::cout << scope[2] << "<" << scope[3] << " ";
    if(state&8) std::cout << scope[0] << ">=" << scope[1] << " ";
    std::cout << ")"<< std::endl;  
#endif
    
    if(evt < 2) {
      if(state&1) {

#ifdef _DEBUG_LEX
	std::cout << "Rule1: " << scope[0] << " = " << scope[1] << std::endl;
#endif

	// state&1 and change on x0 or x1 => R1
	FILTER2(aux,  1-evt, set_domain(scope[evt]) );
      }
      if(IS_OK(wiped) && state&4) {

#ifdef _DEBUG_LEX
	std::cout << "Rule3: " << scope[0] << " < " << scope[1] << std::endl;
#endif

	// state&4 and change on x0 or x1 => R3
	FILTER2(aux,  0, set_max(scope[1].get_max()-1) );
	FILTER2(aux,  1, set_min(scope[0].get_min()+1) );
      }

      // change on x0 or x1 => R2
      if(IS_OK(wiped) &&  !(state&1) && !scope[0].intersect(scope[1])) {

#ifdef _DEBUG_LEX
	std::cout << "Rule2: " << scope[3] << " = 1" << std::endl;
#endif

	FILTER2(aux,  3, set_domain(1) );
	//state |= 1;
	if(!scope[2].get_max()) state |= 4;
      }

      // change on x0 or x1 => R4
      if(IS_OK(wiped) &&  scope[0].get_min() >= scope[1].get_max()) {

#ifdef _DEBUG_LEX
	std::cout << "Rule4: " << scope[2] << " = " << scope[3] << std::endl;
#endif

	FILTER2(aux,  3, set_domain(scope[2]) );
	FILTER2(aux,  2, set_domain(scope[3]) );
	state |= 8;
	if(!scope[3].get_max()) state |= 1;
      }

    } else {
      
      if(evt == 3) {
	if(!scope[3].get_max()) {

#ifdef _DEBUG_LEX
	  std::cout << "Rule1: " << scope[0] << " = " << scope[1] << std::endl;
#endif
	  
	  state |= 1;
	  // change on b1 => R1
	  FILTER2(aux,  0, set_domain(scope[1]) );
	  FILTER2(aux,  1, set_domain(scope[0]) );


#ifdef _DEBUG_LEX
	  std::cout << "Rule6: " << scope[2] << " = 0 " << std::endl;
#endif

	  FILTER2(aux,  2, set_domain(0) );

	}
      } else {
	if( scope[2].get_min() ) {
	  
#ifdef _DEBUG_LEX
	  std::cout << "Rule5: " << scope[1] << " = 1 " << std::endl;
#endif

	  FILTER2(aux,  3, set_domain(1) );

	}
      }

      // state&8 and change on b0 or b1 => R4
      if(state&8) {

#ifdef _DEBUG_LEX
	std::cout << "Rule4: " << scope[2] << " = " << scope[3] << std::endl;
#endif

	FILTER2(aux,  5-evt, set_domain(scope[evt].get_min()) );

	if(!scope[3].get_max()) state |= 1;
      }

      // change on b0 or b1 => R3
      if(scope[2].get_max() < scope[3].get_min()) {

#ifdef _DEBUG_LEX
	std::cout << "Rule3: " << scope[0] << " < " << scope[1] << std::endl;
#endif


	FILTER2(aux,  0, set_max(scope[1].get_max()-1) );
	FILTER2(aux,  1, set_min(scope[0].get_min()+1) );
      }
    }    
  }

#ifdef _DEBUG_LEX
  if(!IS_OK(wiped))
    std::cout << "wipe out of x" << wiped << std::endl; 
#endif
  
  return wiped;
}

std::ostream& Mistral::ConstraintLex::display(std::ostream& os) const {
  os << "(" << scope[0]/*.get_var()*/ << " lex " << scope[1]/*.get_var()*/ << " - " << scope[2]/*.get_var()*/ << "." << scope[3]/*.get_var()*/ << ")";
  return os;
}


// void Mistral::ConstraintLexf::initialise() {
//   Constraint::initialise();
//   for(int i=0; i<3; ++i)
//     trigger_on(_RANGE_, i);
//   set_idempotent(true);
// }

// Mistral::PropagationOutcome Mistral::ConstraintLexf::propagate() {      
//   Mistral::PropagationOutcome wiped = CONSISTENT;
//   Event evt;

//   while(!changes.empty()) {
//     evt = changes.pop();

//     if(evt==2) {
//       if(scope[2].get_min()) {
// 	// => x0 < x1
// 	FILTER( 0, set_max(scope[1].get_max()-1) );
// 	FILTER( 1, set_min(scope[0].get_min()+1) );
//       } else {
// 	// => x0 = x1
// 	FILTER( 0, set_domain(scope[1]) );
// 	FILTER( 1, set_domain(scope[0]) );
//       }
//     } else {
//       if(scope[2].get_min()) {
// 	if(evt) FILTER( 0, set_max(scope[1].get_max()-1) );
// 	else FILTER( 1, set_min(scope[0].get_min()+1) );
//       } else if(!scope[2].get_max()) {
// 	FILTER( 1-evt, set_domain(scope[evt]) );
//       } else {
// 	// in anycase x0 <= x1
// 	FILTER( 0, set_max(scope[1].get_max()) );
// 	FILTER( 1, set_min(scope[0].get_min()) );
//       }
//     }
//   }  
  
//   return wiped;
// }

// std::ostream& Mistral::ConstraintLexf::display(std::ostream& os) const {
//   os << "(" << scope[0] << " lex " << scope[1] << ")";
//   return os;
// }


void Mistral::PredicateAdd::initialise() {
  ConstraintImplementation::initialise();

  trigger_on(_RANGE_, scope[0]);
  trigger_on(_RANGE_, scope[1]);
  trigger_on(_RANGE_, scope[2]);
}

Mistral::PropagationOutcome Mistral::PredicateAdd::rewrite() {
  // TODO: recode that
  
  Mistral::PropagationOutcome wiped = propagate();

  VarArray tmp;
  if(active == 3) {
    ConstraintImplementation *con;
    int i=0;
    for(; i<2; ++i)
      if(scope[i].is_ground()) {
	relax();
	tmp.add(scope[1-i]);
	tmp.add(scope[2]);
	if(scope[i].get_min() == 0) {
	  con = new ConstraintEqual(tmp);
	} else {
	  con = new PredicateOffset(tmp, scope[i].get_min());
	}
	get_solver()->add(Constraint(con, con->type));
      }
  }
  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateAdd::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  int ub_not_revised = 7;
  int lb_not_revised = 7;
  Event cevt;


#ifdef _DEBUG_ADD
  if(_DEBUG_ADD) {
    std::cout << "init prop "
	      << scope[0].get_domain() << " + " 
	      << scope[1].get_domain() << " = " 
	      << scope[2].get_domain() << std::endl;
  }
#endif

  while(ub_not_revised || lb_not_revised) {
    
    if(lb_not_revised&1) {
      // revise scope[0]'lb
      cevt = scope[0].set_min(scope[2].get_min() - scope[1].get_max());
      if(IS_FAIL(cevt)) {
	wiped = FAILURE(0);
	break;
      } else if(cevt != NO_EVENT) {
	// scope[0]'s lb changed and may trigger scope[2]'s lb or scope[1]'s ub
	lb_not_revised |= 4;
	ub_not_revised |= 2;
      }
      lb_not_revised ^= 1;
    }

    if(ub_not_revised&1) {
      // revise scope[0]'ub
      cevt = scope[0].set_max(scope[2].get_max() - scope[1].get_min());
      if(IS_FAIL(cevt)) {
	wiped = FAILURE(0);
	break;
      } else if(cevt != NO_EVENT) {
	// scope[0]'s ub changed and may trigger scope[2]'s ub or scope[1]'s lb
	ub_not_revised |= 4;
	lb_not_revised |= 2;
      }
      ub_not_revised ^= 1;
    }


    if(lb_not_revised&2) {
      // revise scope[1]'lb
      cevt = scope[1].set_min(scope[2].get_min() - scope[0].get_max());
      if(IS_FAIL(cevt)) {
	wiped = FAILURE(1);
	break;
      } else if(cevt != NO_EVENT) {
	// scope[1]'s lb changed and may trigger scope[2]'s lb or scope[0]'s ub
	lb_not_revised |= 4;
	ub_not_revised |= 1;
      }
      lb_not_revised ^= 2;
    }

    if(ub_not_revised&2) {
      // revise scope[1]'ub
      cevt = scope[1].set_max(scope[2].get_max() - scope[0].get_min());
      if(IS_FAIL(cevt)) {
	wiped = FAILURE(1);
	break;
      } else if(cevt != NO_EVENT) {
	// scope[1]'s ub changed and may trigger scope[2]'s ub or scope[0]'s lb
	ub_not_revised |= 4;
	lb_not_revised |= 1;
      }
      ub_not_revised ^= 2;
    }


    if(lb_not_revised&4) {
      // revise scope[2]'lb
      cevt = scope[2].set_min(scope[0].get_min() + scope[1].get_min());
      if(IS_FAIL(cevt)) {
	wiped = FAILURE(2);
	break;
      } else if(cevt != NO_EVENT) {
	// scope[2]'s lb changed and may trigger scope[0]'s lb or scope[1]'s lb
	lb_not_revised |= 2;
	lb_not_revised |= 1;
      }
      lb_not_revised ^= 4;
    }

    if(ub_not_revised&4) {
      // revise scope[2]'ub
      cevt = scope[2].set_max(scope[0].get_max() + scope[1].get_max());
      if(IS_FAIL(cevt)) {
	wiped = FAILURE(2);
	break;
      } else if(cevt != NO_EVENT) {
	// scope[2]'s ub changed and may trigger scope[0]'s ub or scope[1]'s ub
	ub_not_revised |= 2;
	ub_not_revised |= 1;
      }
      ub_not_revised ^= 4;
    }
  }

#ifdef _DEBUG_ADD
  if(_DEBUG_ADD) {
    std::cout << scope[0].get_domain() << " + " 
	      << scope[1].get_domain() << " = " 
	      << scope[2].get_domain() << (IS_OK(wiped) ? " ok" : " fail!") << std::endl << std::endl; 
  }
#endif

  return wiped;
}


Mistral::PropagationOutcome Mistral::PredicateAdd::propagate(const int changed_idx, 
							     const Event evt) {      
  Mistral::PropagationOutcome wiped = CONSISTENT;
  int cidx = changed_idx;
  int cevt = evt;

#ifdef _DEBUG_ADD
  if(_DEBUG_ADD) {
    std::cout << scope[0] << " + " 
	      << scope[1] << " = " 
	      << scope[2] << std::endl; 
    std::cout << scope[0].get_domain() << " + " 
	      << scope[1].get_domain() << " = " 
	      << scope[2].get_domain() << " " << event2str(evt) 
	      << " on " << scope[changed_idx].get_domain() << std::endl 
	      << on[0] << std::endl
	      << on[1] << std::endl
	      << on[2] << std::endl;
  }
#endif

  //_x_ + y = z
  
  //do {
  if(cidx == 0) {
    if(LB_CHANGED(evt)) {
      // max(y) = max(z) - min(x) 
      // update y's ub
      cevt = scope[1].set_max(scope[2].get_max() - scope[0].get_min());
      if(IS_FAIL(cevt)) wiped = FAILURE(1);
	
	if(IS_OK(wiped)) {
	  // min(z) = min(x) + min(y);
	  // update z's lb
	  cevt = scope[2].set_min(scope[0].get_min() + scope[1].get_min());
	  if(IS_FAIL(cevt)) wiped = FAILURE(2);
	}
      }
      if(UB_CHANGED(evt)) {
	// min(y) = min(z) - max(x) 
	// update y's lb
	cevt = scope[1].set_min(scope[2].get_min() - scope[0].get_max());
	if(IS_FAIL(cevt)) wiped = FAILURE(1);
	
	if(IS_OK(wiped)) {
	  // max(z) = max(x) + max(y)
	  // update z's ub
	  cevt = scope[2].set_max(scope[0].get_max() + scope[1].get_max());
	  if(IS_FAIL(cevt)) wiped = FAILURE(2);
	}
      }
    } else if(cidx == 1) {
      if(LB_CHANGED(evt)) {
	// max(x) = max(z) - min(y) 
	// update x's ub
	cevt = scope[0].set_max(scope[2].get_max() - scope[1].get_min());
	if(IS_FAIL(cevt)) wiped = FAILURE(0);
	
	if(IS_OK(wiped)) {
	  // min(z) = min(x) + min(y);
	  // update z's lb
	  cevt = scope[2].set_min(scope[0].get_min() + scope[1].get_min());
	  if(IS_FAIL(cevt)) wiped = FAILURE(2);
	}
      }
      if(UB_CHANGED(evt)) {
	// min(x) = min(z) - max(y) 
	// update x's lb
	cevt = scope[0].set_min(scope[2].get_min() - scope[1].get_max());
	if(IS_FAIL(cevt)) wiped = FAILURE(0);

	if(IS_OK(wiped)) {	
	  // max(z) = max(x) + max(y)
	  // update z's ub
	  cevt = scope[2].set_max(scope[0].get_max() + scope[1].get_max());
	  if(IS_FAIL(cevt)) wiped = FAILURE(2);
	}
      }
    } else {
      if(UB_CHANGED(evt)) {
	// max(x) = max(z) - min(y) 
	// update x's ub
	cevt = scope[0].set_max(scope[2].get_max() - scope[1].get_min());
	if(IS_FAIL(cevt)) wiped = FAILURE(0);

	if(IS_OK(wiped)) {	
	  // max(y) = max(z) - min(x);
	  // update y's ub
	  cevt = scope[1].set_max(scope[2].get_max() - scope[0].get_min());
	  if(IS_FAIL(cevt)) wiped = FAILURE(1);
	}
      }
      if(LB_CHANGED(evt)) {
	// min(x) = min(z) - max(y) 
	// update x's lb
	cevt = scope[0].set_min(scope[2].get_min() - scope[1].get_max());
	if(IS_FAIL(cevt)) wiped = FAILURE(0);

	if(IS_OK(wiped)) {	
	  // min(y) = min(z) - max(x) 
	  // update y's lb
	  cevt = scope[1].set_min(scope[2].get_min() - scope[0].get_max());
	  if(IS_FAIL(cevt)) wiped = FAILURE(1);
	}
      }
      //}
    
    if(IS_OK(wiped)) {
      update(cidx, cevt);
    }
  }

#ifdef _DEBUG_ADD
  if(_DEBUG_ADD) {
    std::cout << on[0] << std::endl
	      << on[1] << std::endl
	      << on[2] << std::endl
	      << scope[0].get_domain() << " + " 
	      << scope[1].get_domain() << " = " 
	      << scope[2].get_domain() << (IS_OK(wiped) ? " ok" : " fail!") << std::endl << std::endl; 
  }
#endif
  
  return wiped;
}

std::ostream& Mistral::PredicateAdd::display(std::ostream& os) const {
  os << scope[2]/*.get_var()*/ << " == (" << scope[0]/*.get_var()*/ << " + " << scope[1]/*.get_var()*/ << ")";
  return os;
}





// int __modulo_fct__(const int x, const int m) {
//   int mod = x%m;
//   if(mod && (mod<0) != (m<0))  mod += m;
//   return mod;
// }


// /*
//   [HERE WE ASSUME 0<k and 0<=a<=b]

// ================================
// |  min([c, d]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if k>b: a                   |
// |                              |
// |  if a<=k<=b: 0               |
// |                              |
// |  so we can assume that k < a |
// |                              |
// |  if (a%k + b - a) >= k: 0    |
// |  else: a%k                   |
// ================================

// [IF k>0 THEN:]

// ================================
// |  min([c, d]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if k>b: a                   |
// |                              |
// |  if a<=k<=b: 0               |
// |                              |
// |  so we can assume that k < a |
// |                              |
// |  if (a%k + b - a) >= k: 0    |
// |  else: a%k                   |
// ================================


// */
// int min_modulo(const int a, const int b, const int k) {
//   int value = a;
//   if(k<=b) {
//     if(k>=a) value = 0;
//     else {
//       int mod = __modulo_fct__(a,k);
//       if((mod + b - a) >= k) value = 0;
//       else value = mod;
//     }
//   }
//   return value;
// }
// /*
// ================================
// |  max([c, d]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if k>b: b                   |
// |                              |
// |  if a<=k<=b: k-1             |
// |                              |
// |  so we can assume that k < a |
// |                              |
// |  if (b%k - b + a) < 0: k-1   |
// |  else: b%k                   |
// ================================
// */
// int max_modulo(const int a, const int b, const int k) {
//   int value = b;
//   if(k<=b) {
//     if(k>=a) value = k-1;
//     else {
//       int mod = __modulo_fct__(b,k);
//       if((mod - b + a) < 0) value = k-1;
//       else value = mod;
//     }
//   }
//   return value;
// }
// /*
// ================================
// |  min([a, b]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if(a%k < c)                 |
// |    a + (c - a%k)             |
// |  if(a%k > d)                 |
// |    a + (d - a%k) + k         |
// |  otherwise: a                |
// ================================
// */
// int min_antimodulo(const int a, const int c, const int d, const int k) {
//   int value = a, mod = a%k;
//   if(mod < c) 
//     value = a + c - mod;
//   else if(mod > d) 
//     value = a + d - mod + k;
//   return value;
// }
// /*
// ================================
// |  max([a, b]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if(b%k < c)                 |
// |    b - (b%k - c) - k         |
// |  if(b%k > d)                 |
// |    b - (b%k - d)             |
// |  otherwise: b                |
// ================================
// */
// int max_antimodulo(const int b, const int c, const int d, const int k) {
//   int value = b, mod = b%k;
//   if(mod < c) 
//     value = b - mod + c - k;
//   else if(mod > d) 
//     value = b - mod + d;
//   return value;
// }

// void Mistral::PredicateModConstant::initialise() {
//   ConstraintImplementation::initialise();
//   trigger_on(_RANGE_, scope[0]);
//   trigger_on(_RANGE_, scope[1]);
// }

// Mistral::PropagationOutcome Mistral::PredicateModConstant::propagate() {      
//   Mistral::PropagationOutcome wiped = CONSISTENT;
  
//   if(IS_FAIL(scope[1].set_max(max_modulo(scope[0].get_min(), scope[0].get_max(), modulo)))) wiped = FAILURE(1);
//   if(IS_FAIL(scope[1].set_min(min_modulo(scope[0].get_min(), scope[0].get_max(), modulo)))) wiped = FAILURE(1);
//   if(IS_FAIL(scope[0].set_min(min_antimodulo(scope[0].get_min(), 
// 					     scope[1].get_min(), 
// 					     scope[1].get_max(), modulo)))) wiped = FAILURE(0);
//   if(IS_FAIL(scope[0].set_max(max_antimodulo(scope[0].get_max(), 
// 					     scope[1].get_min(), 
// 					     scope[1].get_max(), modulo)))) wiped = FAILURE(0);
//   return wiped;
// }

// Mistral::PropagationOutcome Mistral::PredicateModConstant::propagate(const int changed_idx, const Event evt) {
//   Mistral::PropagationOutcome wiped = CONSISTENT;


//   if(IS_FAIL(scope[1].set_max(max_modulo(scope[0].get_min(), scope[0].get_max(), modulo)))) wiped = FAILURE(1);
//   if(IS_FAIL(scope[1].set_min(min_modulo(scope[0].get_min(), scope[0].get_max(), modulo)))) wiped = FAILURE(1);
//   if(IS_FAIL(scope[0].set_min(min_antimodulo(scope[0].get_min(), 
// 					     scope[1].get_min(), 
// 					     scope[1].get_max(), modulo)))) wiped = FAILURE(0);
//   if(IS_FAIL(scope[0].set_max(max_antimodulo(scope[0].get_max(), 
// 					     scope[1].get_min(), 
// 					     scope[1].get_max(), modulo)))) wiped = FAILURE(0);

//   return wiped;
// }

// std::ostream& Mistral::PredicateModConstant::display(std::ostream& os) const {
//   os << scope[1]/*.get_var()*/ << " == (" << scope[0]/*.get_var()*/ << " % " << modulo << ")";
//   return os;
// }




// /*


//   [HERE WE ASSUME 0<=a<=b AND 0<=c<=d]

//   min([a,b] % [c,d]) =

//   if c>b: a

//   if [a,b] inter [c,d]: 0


//   so we can assume that d < a

//   if (b-a+1) >= c: 0 
     
// min(k \in [c,d])
//   if (a%k + b - a) >= k: 0
//   else: a%k

//  */


// void Mistral::PredicateMod::initialise() {
//   ConstraintImplementation::initialise();

//   trigger_on(_RANGE_, scope[0]);
//   trigger_on(_RANGE_, scope[1]);
//   trigger_on(_RANGE_, scope[2]);
// }

// Mistral::PropagationOutcome Mistral::PredicateMod::rewrite() {
//    Mistral::PropagationOutcome wiped = propagate();



//   // VarArray tmp;
//   // if(active == 3) {
//   //   ConstraintImplementation *con;
//   //   int i=0;
//   //   for(; i<2; ++i)
//   //     if(scope[i].is_ground()) {
//   // 	relax();
//   // 	tmp.add(scope[1-i]);
//   // 	tmp.add(scope[2]);
//   // 	if(scope[i].get_min() == 0) {
//   // 	  con = new ConstraintEqual(tmp);
//   // 	} else {
//   // 	  con = new PredicateOffset(tmp, scope[i].get_min());
//   // 	}
//   // 	get_solver()->add(Constraint(con, con->type));
//   //     }
//   // }
//   return wiped;
// }

// Mistral::PropagationOutcome Mistral::PredicateMod::propagate() { 
//   if(IS_FAIL(scope[1].remove(0))) return FAILURE(1);
//   return filter();
// }

// Mistral::PropagationOutcome Mistral::PredicateMod::filter() {      
//   Mistral::PropagationOutcome wiped = CONSISTENT;


// #ifdef _DEBUG_MOD
//   if(_DEBUG_MOD) {
//     std::cout << scope[0] << " % " 
// 	      << scope[1] << " = " 
// 	      << scope[2] << std::endl; 
//     std::cout << scope[0].get_domain() << " % " 
// 	      << scope[1].get_domain() << " = " 
// 	      << scope[2].get_domain() << " " << event2str(evt) 
// 	      << " on " << scope[changed_idx].get_domain() << std::endl 
// 	      << on[0] << std::endl
// 	      << on[1] << std::endl
// 	      << on[2] << std::endl;
//   }
// #endif


//   //x % y = z

//   // the initial bounds of z 
//   int min_z = scope[2].get_min();
//   int max_z = scope[2].get_max();
//   // the current bounds of z
//   int lb_z = +INFTY;
//   int ub_z = -INFTY;

  
//   // the initial bounds of x
//   int min_x = scope[0].get_min();
//   int max_x = scope[0].get_max();
//   // the current bounds of x
//   int lb_x = +INFTY;
//   int ub_x = -INFTY;


//   int modmin, modmax, antimodmin, antimodmax;
//   int modnext = scope[1].get_min();
//   int modulo = modnext-1;
//   for(;
//       IS_OK(wiped) && modulo<modnext; 
//       modnext = scope[1].next(modulo)) {

//     modulo = modnext;

//     //compute the min and max possible values for x%modulo, and update the current lb/ub of z accordingly
//     modmin = min_modulo(min_x, max_x, modulo);
//     if(modmin < lb_z) {
//       lb_z = modmin;
//     } 
//     modmax = max_modulo(min_x, max_x, modulo);
//     if(modmax > ub_z) {
//       ub_z = modmax;
//     } 

//     //compute the min and max possible values for x given that x%modulo = z, and update the current lb/ub of x accordingly
//     antimodmin = min_antimodulo(min_x, min_z, max_z, modulo);
//     if(antimodmin < lb_x) {
//       lb_x = antimodmin;
//     } 
//     antimodmax = max_antimodulo(min_x, min_z, max_z, modulo);
//     if(antimodmax > ub_x) {
//       ub_x = antimodmax;
//     } 

//     if(modmin>max_z || modmax<min_z || antimodmin>max_x || antimodmax<min_x) {
//       // forbidden value
//       if(IS_FAIL(scope[1].remove(modulo))) wiped = FAILURE(1);
//     }
//   }

//   if(IS_OK(wiped) && IS_FAIL(scope[2].set_max(max_z))) wiped = FAILURE(1);
//   if(IS_OK(wiped) && IS_FAIL(scope[2].set_min(min_z))) wiped = FAILURE(1);
//   if(IS_OK(wiped) && IS_FAIL(scope[0].set_min(min_x))) wiped = FAILURE(0);
//   if(IS_OK(wiped) && IS_FAIL(scope[0].set_max(max_x))) wiped = FAILURE(0);


// #ifdef _DEBUG_MOD
//   if(_DEBUG_MOD) {
//     std::cout << on[0] << std::endl
// 	      << on[1] << std::endl
// 	      << on[2] << std::endl
// 	      << scope[0].get_domain() << " % " 
// 	      << scope[1].get_domain() << " = " 
// 	      << scope[2].get_domain() << (IS_OK(wiped) ? " ok" : " fail!") << std::endl << std::endl; 
//   }
// #endif
  

//   return wiped;
// }

// Mistral::PropagationOutcome Mistral::PredicateMod::propagate(const int changed_idx, 
// 							     const Event evt) {    

//  return filter();
  
// //   Mistral::PropagationOutcome wiped = CONSISTENT;

// // #ifdef _DEBUG_MOD
// //   if(_DEBUG_MOD) {
// //     std::cout << scope[0] << " % " 
// // 	      << scope[1] << " = " 
// // 	      << scope[2] << std::endl; 
// //     std::cout << scope[0].get_domain() << " % " 
// // 	      << scope[1].get_domain() << " = " 
// // 	      << scope[2].get_domain() << " " << event2str(evt) 
// // 	      << " on " << scope[changed_idx].get_domain() << std::endl 
// // 	      << on[0] << std::endl
// // 	      << on[1] << std::endl
// // 	      << on[2] << std::endl;
// //   }
// // #endif



// //   //x % y = z

// //   // the initial bounds of z 
// //   int min_z = scope[2].get_min();
// //   int max_z = scope[2].get_max();
// //   // the current bounds of z
// //   int lb_z = +INFTY;
// //   int ub_z = -INFTY;

  
// //   // the initial bounds of x
// //   int min_x = scope[0].get_min();
// //   int max_x = scope[0].get_max();
// //   // the current bounds of x
// //   int lb_x = +INFTY;
// //   int ub_x = -INFTY;


// //   int modmin, modmax, antimodmin, antimodmax;
// //   int modnext = scope[1].get_min();
// //   int modulo = modnext-1;
// //   for(;
// //       IS_OK(wiped) && modulo<modnext; 
// //       modnext = scope[1].next(modulo)) {

// //     modulo = modnext;

// //     //compute the min and max possible values for x%modulo, and update the current lb/ub of z accordingly
// //     modmin = min_modulo(min_x, max_x, modulo);
// //     if(modmin < lb_z) {
// //       lb_z = modmin;
// //     } 
// //     modmax = max_modulo(min_x, max_x, modulo);
// //     if(modmax > ub_z) {
// //       ub_z = modmax;
// //     } 

// //     //compute the min and max possible values for x given that x%modulo = z, and update the current lb/ub of x accordingly
// //     antimodmin = min_antimodulo(min_x, min_z, max_z, modulo);
// //     if(antimodmin < lb_x) {
// //       lb_x = antimodmin;
// //     } 
// //     antimodmax = max_antimodulo(min_x, min_z, max_z, modulo);
// //     if(antimodmax > ub_x) {
// //       ub_x = antimodmax;
// //     } 

// //     if(modmin>max_z || modmax<min_z || antimodmin>max_x || antimodmax<min_x) {
// //       // forbidden value
// //       if(IS_FAIL(scope[1].remove(modulo))) wiped = FAILURE(1);
// //     }
// //   }

// //   if(IS_OK(wiped) && IS_FAIL(scope[2].set_max(max_z))) wiped = FAILURE(1);
// //   if(IS_OK(wiped) && IS_FAIL(scope[2].set_min(min_z))) wiped = FAILURE(1);
// //   if(IS_OK(wiped) && IS_FAIL(scope[0].set_min(min_x))) wiped = FAILURE(0);
// //   if(IS_OK(wiped) && IS_FAIL(scope[0].set_max(max_x))) wiped = FAILURE(0);



// //   // // int ub, lb, modulo, target, k, vali, vnxt;  
// //   // // if( changed_idx ) { // prune x0
// //   // //   if( scope[1].is_ground() ) { // modulo is known
// //   // //     modulo = scope[1].get_min();
// //   // //     if( modulo == 1 ) {// special case
// //   // // 	if(IS_FAIL(scope[0].set_domain(0))) wiped = FAILURE(0);
// //   // //     } else if( scope[2].is_ground() ) { // target is known
// //   // // 	target = scope[2].get_value();
// //   // // 	if(target >= modulo) wiped = FAILURE(0);
      
// //   // // 	if(IS_OK(wiped)) {
// //   // // 	  // positive/negative target
// //   // // 	  if( target < 0 ) {
// //   // // 	    if(IS_FAIL(scope[0].set_max(0))) wiped = FAILURE(0);
// //   // // 	  } else if( target > 0 ) {
// //   // // 	    if(IS_FAIL(scope[0].set_min(0))) wiped = FAILURE(0);
// //   // // 	    //consistent = scope[0]->setMin(0);
// //   // // 	  } else { 
// //   // // 	    if(IS_FAIL(scope[0].set_domain(0))) wiped = FAILURE(0);
// //   // // 	    //consistent = scope[0]->setDomain(0);
// //   // // 	  }

// //   // // 	  if(IS_OK(wiped)) {
// //   // // 	    // remove intervals [target+1+k*modulo..target+k*(modulo+1)-1]
// //   // // 	    k = (scope[0].get_max()-target-1)/modulo;
// //   // // 	    while(IS_OK(wiped)) {
// //   // // 	      lb = (target+1+k*modulo);
// //   // // 	      ub = (target+(k+1)*modulo-1);
// //   // // 	      if( ub < scope[0].get_min() ) break;
// //   // // 	      if(IS_FAIL(scope[0].remove_range( lb, ub ))) wiped = FAILURE(0);
// //   // // 	      --k;
// //   // // 	    }
// //   // // 	  }
// //   // // 	}      
// //   // //     } else {
// //   // // 	// prune x0 with respect to x2
	
// //   // // 	vnxt = scope[0].get_min();
// //   // // 	vali = vnxt-1;
// //   // // 	while(vali<vnxt) {
// //   // // 	  vali = vnxt;
// //   // // 	  k = (vali % modulo);
// //   // // 	  if(!scope[2].contain( k ) && IS_FAIL(scope[0].remove( vali ))) wiped = FAILURE(0);
// //   // // 	  vnxt = scope[0].next(vali);
// //   // // 	} 
	
// //   // //     }
// //   // //   } else {
// //   // //     // modulo is not known we want to prune x0
// //   // //     //[TODO!!]
// //   // //   }
// //   // // } else if( changed_idx != 2 ) {
// //   // //   // prune x2

// //   // //   if( scope[1].is_ground() ) { // modulo is known
// //   // //     if( scope[0].is_ground() ) {
// //   // // 	//consistent = scope[2]->setDomain( (scope[0].get_min() % scope[1].get_min()) );
// //   // // 	if(IS_FAIL(scope[2].set_domain( (scope[0].get_value() % scope[1].get_value()) ))) wiped = FAILURE(2);
// //   // //     } else {
	
// //   // // 	//[TODO: this is wrong (?!!)]
// //   // // 	modulo = scope[1].get_value();
// //   // // 	ub = scope[0].get_max();
// //   // // 	if( ub > 0 && modulo <= ub )
// //   // // 	  ub = modulo-1;
// //   // // 	else ub = 0;
// //   // // 	lb = scope[0].get_min();
// //   // // 	if( lb < 0 && 1-modulo > lb )
// //   // // 	  lb = 1-modulo;
// //   // // 	else lb = 0;
	
// //   // // 	if(IS_FAIL(scope[2].set_max(ub))) wiped = FAILURE(2);
// //   // // 	else if(IS_FAIL(scope[2].set_min(lb))) wiped = FAILURE(2);
// //   // // 	//consistent = scope[2]->setMax( ub ) && scope[2]->setMin( lb );
	
// //   // // 	vnxt = scope[2].get_min();
// //   // // 	vali = vnxt-1;
// //   // // 	while(IS_OK(wiped) && vali < vnxt) {
// //   // // 	  vali = vnxt;
	  
// //   // // 	  k = vali;
// //   // // 	  if( k > 0 ) {
// //   // // 	    lb = (scope[0].get_min()/modulo)*modulo;
// //   // // 	    k = std::max( k, lb+k );
// //   // // 	    while( !scope[0].contain(k) && k <= scope[0].get_max() )  {
// //   // // 	      k+=modulo;
// //   // // 	    }
// //   // // 	    if(k && k > scope[0].get_max()) {
// //   // // 	      if(IS_FAIL(scope[2].remove( vali ))) wiped = FAILURE(2);
// //   // // 	    } 
// //   // // 	  }
// //   // // 	  else {
// //   // // 	    ub = (scope[0].get_max()/modulo)*modulo;
// //   // // 	    if(k) k = std::min( k, ub+k );
// //   // // 	    else k = ub;
// //   // // 	    while( !scope[0].contain(k) && k >= scope[0].get_min() ) {
// //   // // 	      k-=modulo;
// //   // // 	    }
// //   // // 	    if(k < scope[0].get_min()) {
// //   // // 	      if(IS_FAIL(scope[2].remove( vali ))) wiped = FAILURE(2);
// //   // // 	    } 
	    
// //   // // 	    vnxt = scope[2].next(k);
// //   // // 	  }
// //   // // 	} 
// //   // //     }
// //   // //   } else {
// //   // //     // modulo is not known we want to prune x2
// //   // //   }
// //   // // } 
  


// // #ifdef _DEBUG_MOD
// //   if(_DEBUG_MOD) {
// //     std::cout << on[0] << std::endl
// // 	      << on[1] << std::endl
// // 	      << on[2] << std::endl
// // 	      << scope[0].get_domain() << " % " 
// // 	      << scope[1].get_domain() << " = " 
// // 	      << scope[2].get_domain() << (IS_OK(wiped) ? " ok" : " fail!") << std::endl << std::endl; 
// //   }
// // #endif
  
// //   return wiped;
// }

// std::ostream& Mistral::PredicateMod::display(std::ostream& os) const {
//   os << scope[2]/*.get_var()*/ << " == (" << scope[0]/*.get_var()*/ << " % " << scope[1]/*.get_var()*/ << ")";
//   return os;
// }


void Mistral::PredicateMul::initialise() {
  ConstraintImplementation::initialise();

  trigger_on(_RANGE_, scope[0]);
  trigger_on(_RANGE_, scope[1]);
  trigger_on(_RANGE_, scope[2]);

  GlobalConstraint::initialise();

  enforce_nfc1 = false;
}

Mistral::PropagationOutcome Mistral::PredicateMul::rewrite() {
  // TODO: recode that
  
Mistral::PropagationOutcome wiped = propagate();

  VarArray tmp;
  if(active.size == 2) {
    int i=0;
    for(; i<2; ++i)
      if(scope[i].is_ground()) {
	relax();
	tmp.add(scope[1-i]);
	tmp.add(scope[2]);
	if(scope[i].get_min() == 1) {
	  get_solver()->add(Constraint(new ConstraintEqual(tmp)));
	} else if(scope[i].get_min() != 0) {
	  get_solver()->add(Constraint(new PredicateFactor(tmp, scope[i].get_min())));
	}
      }
  }
  return wiped;
}

inline int xtimey( const int x, const int y, int& r )
{
  r = 0 ;

#ifdef _DEBUG_MUL
  std::cout << x << " * " << y << " = " << (x*y) << std::endl;
#endif

  return (x*y);
}

inline int xovery( const int x, const int y, int& r )
{

#ifdef _DEBUG_MUL
  std::cout << x << " / " << y << " = " ;
#endif

  if(y) {
    r = (x%y != 0);

#ifdef _DEBUG_MUL
    std::cout << (x/y) << std::endl;
#endif

    return x/y;
  }

#ifdef _DEBUG_MUL
    std::cout << (x<0 ? -INFTY : INFTY) << std::endl;
#endif

  return (x<0 ? -INFTY : INFTY);
}

inline int yoverx( const int x, const int y, int& r )
{

#ifdef _DEBUG_MUL
  std::cout << y << " / " << x << " = " ;
#endif

  if(x) {
    r = -(y%x != 0);

#ifdef _DEBUG_MUL
    std::cout << (y/x) << std::endl;
#endif

    return y/x;
  }

#ifdef _DEBUG_MUL
    std::cout << (y<0 ? -INFTY : INFTY) << std::endl;
#endif

  return (y<0 ? -INFTY : INFTY);
}



Mistral::PropagationOutcome Mistral::PredicateMul::revise_division(const int X, const int Y, const int Z) {
  // revise the domain of Z = X/Y (because Z*Y = X)
  Mistral::PropagationOutcome wiped = CONSISTENT;


#ifdef _DEBUG_MUL
  if(_DEBUG_MUL) {

  std::cout << "revise bounds of " << scope[Z].get_domain() << " = " 
	    << scope[X].get_domain() << "/" << scope[Y].get_domain() 
	    << " = [" << min_neg[Z] << ".." << max_neg[Z] << "|" 
	    << (zero[Z] ? "{0}" : "_") << "|" 
	    << min_pos[Z] << ".." << max_pos[Z] << "]" 
	    << std::endl; 
 
  }
#endif

  
  //int lb_pos=1, ub_pos=+INFTY, lb_neg=-INFTY, ub_neg=-1, lb_aux, ub_aux;
  
  // we start with all bounds at their previous values, and update them if necessary
  // (in which case we set up the pruning flag)
  int lb_pos=min_pos[Z], ub_pos=max_pos[Z], lb_neg=min_neg[Z], ub_neg=max_neg[Z], 
    nlb1, nlb2, nub1, nub2;
    //lb_aux, ub_aux;
  bool pruning_flag = false, pzero = false, ppos = max_pos[Z]<=0, pneg = min_neg[Z]>=0;
  
  // first rule: if X can be 0, and Y can be 0, then Z can be anything
  if(!zero[X]) { // if X cannot be 0, then neither can Y nor Z
    if(zero[Z]) {
      //zero[Z] = 0;
      pzero = true;
      pruning_flag = true;
    }
  } else { 
    if(!min_neg[X] && !max_pos[X] && !zero[Y]) {
      // if X must be 0 and Y cannot be 0, then Z must be 0.
      if(lb_neg<0 || ub_pos>0) {
	lb_neg = 0;
	ub_pos = 0;
	pruning_flag = true;
      }
    } 
  }

  if(lb_neg != ub_pos && (!zero[X] || !zero[Y])) { // if X and Y can both be 0, we cannot deduce anything
    if(IS_OK(wiped)) {
      if(max_pos[Z]>0) {
	// revise the positive part of Z's domain (if it has one)
	nlb1 = nlb2 = INFTY; //lb_neg;
	nub1 = nub2 = 0; //ub_neg;
	
	// it can either be the positive parts of X and Y:
	if(max_pos[X]>0 && max_pos[Y]>0) {
	  // compute the bounds
	  nlb1 = (int)(ceil((double)(min_pos[X])/(double)(max_pos[Y])));
	  nub1 = (int)(floor((double)(max_pos[X])/(double)(min_pos[Y])));
	}

	// or the negative parts of X and Y:
	if(min_neg[X]<0 && min_neg[Y]<0) {
	  // compute the bounds
	  nlb2 = (int)(ceil((double)(max_neg[X])/(double)(min_neg[Y])));
	  nub2 = (int)(floor((double)(min_neg[X])/(double)(max_neg[Y])));
	}
	if(nlb1>nlb2) nlb1 = nlb2;
	if(nub1<nub2) nub1 = nub2;

	if(lb_pos<nlb1) {
	  lb_pos = nlb1;
	  pruning_flag = true;
	}
	if(ub_pos>nub1) {
	  ub_pos = nub1;
	  pruning_flag = true;
	}

	if(lb_pos > max_pos[Z] || ub_pos < min_pos[Z]) ppos = true;

      } else if(pzero || !zero[Z]) // if(lb_pos || ub_pos)
	{
	  lb_pos = min_neg[Z];
	  ub_pos = max_neg[Z];
	} else {
	lb_pos = ub_pos = 0; 
      }

     if(min_neg[Z]<0) {
       // revise the negative part of Z's domain (if it has one)
       nlb1 = nlb2 = 0; //lb_pos;
       nub1 = nub2 = -INFTY; //ub_pos;
	
	// it can either be the negitive part of X and the positive part of Y:
	if(min_neg[X]<0 && max_pos[Y]>0) {
	  // compute the bounds

	  nlb1 = (int)(ceil((double)(min_neg[X])/(double)(min_pos[Y])));
	  nub1 = (int)(floor((double)(max_neg[X])/(double)(max_pos[Y])));
	}
	// or the negitive part of Y and the positive part of X:
	if(max_pos[X]>0 && min_neg[Y]<0) {
	  // compute the bounds
	  nlb2 = (int)(ceil((double)(max_pos[X])/(double)(max_neg[Y])));
	  nub2 = (int)(floor((double)(min_pos[X])/(double)(min_neg[Y])));
	}

	if(nlb1>nlb2) nlb1 = nlb2;
	if(nub1<nub2) nub1 = nub2;

	if(lb_neg<nlb1) {
	  lb_neg = nlb1;
	  pruning_flag = true;
	}
	if(ub_neg>nub1) {
	  ub_neg = nub1;
	  pruning_flag = true;
	}

	if(lb_neg > max_neg[Z] || ub_neg < min_neg[Z]) pneg = true;

     } else if(pzero || !zero[Z])// if(lb_neg || ub_neg)
       {
	 lb_neg = min_pos[Z];
	 ub_neg = max_pos[Z];
       } else {
       lb_neg = ub_neg = 0;
     }
    }
  }

  if(pneg && (pzero || !zero[Z]) && ppos) {
    wiped = FAILURE(Z);
  } else if(pruning_flag) {
    
#ifdef _DEBUG_MUL
  if(_DEBUG_MUL) {

    std::cout << "set bounds to " // << scope[Z].get_domain() << " = " 
	      << "[" << lb_neg << ".." << ub_neg << "|" 
	      << (!pzero&&zero[Z] ? "{0}" : "_") << "|" 
	      << lb_pos << ".." << ub_pos << "]" << std::endl; 
 
  }
#endif
    
    wiped = prune(lb_neg, ub_neg, lb_pos, ub_pos, pzero, Z);
  }

  return wiped;
}


Mistral::PropagationOutcome Mistral::PredicateMul::revise_multiplication(const int X, const int Y, const int Z) {
  // revise the domain of Z = X*Y 
  Mistral::PropagationOutcome wiped = CONSISTENT;


#ifdef _DEBUG_MUL
  if(_DEBUG_MUL) {

  std::cout << "revise bounds of " << scope[Z].get_domain() << " = " 
	    << scope[X].get_domain() << "*" << scope[Y].get_domain() 
	    << std::endl; 
 
  }
#endif

  int lb_pos=min_pos[Z], ub_pos=max_pos[Z], lb_neg=min_neg[Z], ub_neg=max_neg[Z], 
    nlb1, nlb2, nub1, nub2;
  
  bool pruning_flag = false, pzero = false;
  
  // if X = 0 or Y = 0, then Z = 0
  if( zero[Z] &&
      !zero[X] && !zero[Y]) {
    //zero[Z] = 0;
    pzero = true;
    pruning_flag = true;
  }
  else if((!min_neg[X] && !max_pos[X]) || (!min_neg[Y] && !max_pos[Y])) { 
    lb_neg = 0;
    ub_pos = 0;
    pruning_flag = true;
  }

  if(lb_neg != ub_pos) { 
    if(IS_OK(wiped)) {
      if(max_pos[Z]>0) {
	// revise the positive part of Z's domain (if it has one)
	nlb1 = nlb2 = INFTY; //lb_neg;
	nub1 = nub2 = 0; //ub_neg;

	// it can either be the positive parts of X and Y:
	if(max_pos[X]>0 && max_pos[Y]>0) {
	  // compute the bounds
	  nub1 = max_pos[X] * max_pos[Y];
	  nlb1 = min_pos[X] * min_pos[Y];
	}
	// or the negative parts of X and Y:
	if(min_neg[X]<0 && min_neg[Y]<0) {
	  // compute the bounds
	  nub2 = min_neg[X] * min_neg[Y];
	  nlb2 = max_neg[X] * max_neg[Y];
	}
	if(nlb1>nlb2) nlb1 = nlb2;
	if(nub1<nub2) nub1 = nub2;

	if(lb_pos<nlb1) {
	  lb_pos = nlb1;
	  pruning_flag = true;
	}
	if(ub_pos>nub1) {
	  ub_pos = nub1;
	  pruning_flag = true;
	}
      } else if(pzero || !zero[Z]) {
	lb_pos = min_neg[Z];
	ub_pos = max_neg[Z];
      } else {
	lb_pos = ub_pos = 0;
      }

     if(min_neg[Z]<0) {
       // revise the negative part of Z's domain (if it has one)
       nlb1 = nlb2 = 0; //lb_pos;
       nub1 = nub2 = -INFTY; //ub_pos;
	
	// it can either be the negitive part of X and the positive part of Y:
	if(min_neg[X]<0 && max_pos[Y]>0) {
	  // compute the bounds
	  nub1 = max_neg[X] * min_pos[Y];
	  nlb1 = min_neg[X] * max_pos[Y];
	}
	// or the negitive part of Y and the positive part of X:
	if(max_pos[X]>0 && min_neg[Y]<0) {
	  // compute the bounds
	  nub2 = max_neg[Y] * min_pos[X];
	  nlb2 = min_neg[Y] * max_pos[X];
	}

	if(nlb1>nlb2) nlb1 = nlb2;
	if(nub1<nub2) nub1 = nub2;
	
	if(lb_neg<nlb1) {
	  lb_neg = nlb1;
	  pruning_flag = true;
	}
	if(ub_neg>nub1) {
	  ub_neg = nub1;
	  pruning_flag = true;
	}
     }  else if(pzero || !zero[Z]) {
       lb_neg = min_pos[Z];
       ub_neg = max_pos[Z];
     } else {
       lb_neg = ub_neg = 0;
     }
    }
  }

  if(pruning_flag) {
#ifdef _DEBUG_MUL
  if(_DEBUG_MUL) {

  std::cout << "set bounds to " // << scope[Z].get_domain() << " = " 
	    << "[" << lb_neg << ".." << ub_neg << "|" 
	    << (zero[Z] ? "{0}" : "_") << "|" 
	    << lb_pos << ".." << ub_pos << "]" << std::endl; 
 
  }
#endif
  
    wiped = prune(lb_neg, ub_neg, lb_pos, ub_pos, pzero, Z);
  }
  
  return wiped;
}


Mistral::PropagationOutcome Mistral::PredicateMul::prune(const int lb_neg, 
							 const int ub_neg, 
							 const int lb_pos, 
							 const int ub_pos,
							 const bool pzero,
							 const int Z) {

  Event evt;
  Mistral::PropagationOutcome wiped = CONSISTENT;
  
  if(ub_pos < lb_neg) wiped = FAILURE(Z);
  else {
    if(lb_neg>min_neg[Z]) {
      evt = scope[Z].set_min( lb_neg );
      if( IS_FAIL(evt) ) wiped = FAILURE(Z);
      else {
	if(changes.contain(Z)) {
	  event_type[Z] |= evt;
	} else {
	  event_type[Z] = evt;
	  changes.add(Z);
	}
	min_neg[Z] = scope[Z].get_min();
      }
    }
    if(IS_OK(wiped) && ub_pos<max_pos[Z]) {
      evt = scope[Z].set_max( ub_pos );
      if( IS_FAIL(evt) ) wiped = FAILURE(Z);
      else {
	if(changes.contain(Z)) {
	  event_type[Z] |= evt;
	} else {
	  event_type[Z] = evt;
	  changes.add(Z);
	}
	max_pos[Z] = scope[Z].get_max();
      }
    }
    if(IS_OK(wiped) && (lb_pos>=min_neg[Z] || ub_neg<=max_pos[Z])) { 
      if(lb_pos-1>ub_neg && (pzero || (!zero[Z] && (min_pos[Z]<lb_pos || max_neg[Z]>ub_neg)))) {
	evt = scope[Z].remove_interval(ub_neg+1, lb_pos-1);
	if( IS_FAIL(evt) ) wiped = FAILURE(Z);
	else {
	  if(changes.contain(Z)) {
	    event_type[Z] |= evt;
	  } else {
	    event_type[Z] = evt;
	    changes.add(Z);
	  }
	  zero[Z] = 0;
	  min_pos[Z] = scope[Z].get_min_pos();
	  max_neg[Z] = scope[Z].get_max_neg();
	}
      } else {
	if(lb_pos>1 && min_pos[Z]<lb_pos) {
	  evt = scope[Z].remove_interval(1, lb_pos-1);
	  if( IS_FAIL(evt) ) wiped = FAILURE(Z);
	  else {
	    if(changes.contain(Z)) {
	      event_type[Z] |= evt;
	    } else {
	      event_type[Z] = evt;
	      changes.add(Z);
	    }
	    min_pos[Z] = scope[Z].get_min_pos();
	  }
	}
	if(ub_neg<-1 && max_neg[Z]>ub_neg) {
	  evt = scope[Z].remove_interval(ub_neg+1, -1);
	  if( IS_FAIL(evt) ) wiped = FAILURE(Z);
	  else {
	    if(changes.contain(Z)) {
	      event_type[Z] |= evt;
	    } else {
	      event_type[Z] = evt;
	      changes.add(Z);
	    }
	    max_neg[Z] = scope[Z].get_max_neg();
	  }
	}
      }
    }
  }

  if((min_neg[Z]>0 && min_neg[Z]<min_pos[Z]) ||
     (!min_neg[Z] && !zero[Z]))
    min_neg[Z] = min_pos[Z];

  if((max_pos[Z]<0 && max_pos[Z]>max_neg[Z]) ||
     (!max_pos[Z] && !zero[Z]))
    max_pos[Z] = max_neg[Z];


#ifdef _DEBUG_MUL
  if(_DEBUG_MUL) {

    if(IS_OK(wiped)) {
      std::cout << scope[0].get_domain() 
		<< " * " << scope[1].get_domain() 
		<< " = " << scope[2].get_domain() << std::endl;
    } else std::cout << "FAIL!" << std::endl ;
    
    std::cout
      << " now in [" << min_neg[Z] << ".." << max_neg[Z] << "|" 
      << (zero[Z] ? "{0}" : "_") << "|" 
      << min_pos[Z] << ".." << max_pos[Z] << "]" << std::endl;
    
    //std::cout << std::endl; 
    
  }
#endif


  return wiped;
}




Mistral::PropagationOutcome Mistral::PredicateMul::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;
  
#ifdef _DEBUG_MUL
  if(_DEBUG_MUL) {

  std::cout << std::endl << std::endl << "propagate " // << this 
	    << std::endl << scope[0].get_domain() 
 	    << " * " << scope[1].get_domain() 
     	    << " = " << scope[2].get_domain() << std::endl;
 
  }
#endif

  bool is_ground = true;
  int evt_idx, i;

  for(i=0; i<3; ++i) {
    max_pos[i] = scope[i].get_max();
    min_neg[i] = scope[i].get_min();

    is_ground &= (min_neg[i] == max_pos[i]);

    zero[i] = scope[i].contain(0);

    if(min_neg[i]<0)
      max_neg[i] = scope[i].get_max_neg();
    else
      max_neg[i] = min_neg[i];

    if(max_pos[i]>0)
      min_pos[i] = scope[i].get_min_pos();
    else
      min_pos[i] = max_pos[i];

  }

  
  // #ifdef _DEBUG_MUL
  //   if(_DEBUG_MUL) {
  
  //   std::cout << scope[0]/*.get_var()*/ << " in " << scope[0].get_domain() 
  // 	    << " * " << scope[1]/*.get_var()*/ << " in " << scope[1].get_domain() 
  //     	    << " = " << scope[2]/*.get_var()*/ << " in " << scope[2].get_domain() << std::endl;
  
  //   }
  // #endif
    
  
  if(!is_ground) {
    while(IS_OK(wiped) && !changes.empty()) {
      
      evt_idx = changes.pop();
      
#ifdef _DEBUG_MUL
      if(_DEBUG_MUL) {
	
	std::cout << std::endl << "react to " << scope[evt_idx]/*.get_var()*/ << " in " 
		  << scope[evt_idx].get_domain() 
		  << (LB_CHANGED(event_type[evt_idx]) ? " (change on LB) " : "")
		  << (UB_CHANGED(event_type[evt_idx]) ? " (change on UB) " : "")
		  << std::endl;
	
      }
#endif
      
      // x0 * x1 = x2 
      if(evt_idx == 2) {
	for(int i=0; IS_OK(wiped) && i<2; ++i)  // once for scope[0], once for scope[1]
	  // we update x0 and x1 :  xi = x2/x1-i
	  wiped = revise_division(2, 1-i, i);
      } else {
	// update x[2]
	wiped = revise_multiplication(evt_idx, 1-evt_idx, 2);
	// update x[1-evt_idx]
	if(IS_OK(wiped)) wiped = revise_division(2, evt_idx, 1-evt_idx);
      }
    }

  } else {
    
    if(min_neg[0] * min_neg[1] != min_neg[2]) wiped = 0;

  }

#ifdef _DEBUG_MUL
  if(_DEBUG_MUL) {

    std::cout << std::endl;
    if(IS_OK(wiped)) {
      
      std::cout << "  ===> " << scope[0]/*.get_var()*/ << " in " << scope[0].get_domain() 
		<< " * " << scope[1]/*.get_var()*/ << " in " << scope[1].get_domain() 
		<< " = " << scope[2]/*.get_var()*/ << " in " << scope[2].get_domain() << std::endl;

    } else {

      std::cout << "  ===> INCONSISTENT!!" << std::endl;

    }
 
  }
#endif
  
  return wiped;
}
  
std::ostream& Mistral::PredicateMul::display(std::ostream& os) const {
  os << scope[2]/*.get_var()*/ << " == (" << scope[0]/*.get_var()*/ << " * " << scope[1]/*.get_var()*/ << ")";
  return os;
}



// void Mistral::PredicateDiv::initialise() {
//   Constraint::initialise();
//   trigger_on(_RANGE_, scope[0]);
//   trigger_on(_RANGE_, scope[1]);
//   trigger_on(_RANGE_, scope[2]);
//   set_idempotent(true);
// }

// Mistral::PropagationOutcome Mistral::PredicateDiv::rewrite() {
//   Mistral::PropagationOutcome wiped = propagate();

//   VarArray tmp;
//   if(active.size == 2) {
//     int i=0;
//     for(; i<2; ++i)
//       if(scope[i].is_ground()) {
// 	relax();
// 	tmp.add(scope[1-i]);
// 	tmp.add(scope[2]);
// 	if(scope[i].get_min() == 1) {
// 	  solver->add(new ConstraintEqual(tmp));
// 	} else if(scope[i].get_min() != 0) {
// 	  solver->add(new PredicateFactor(tmp, scope[i].get_min()));
// 	}
//       }
//   }
//   return wiped;
// }




// Mistral::PropagationOutcome Mistral::PredicateDiv::revise_integer_division(const int X, const int Y, const int Z) {
//   // revise the domain of Z = X/Y (Z is the integer part of X/Y)
//   Mistral::PropagationOutcome wiped = CONSISTENT;


// #ifdef _DEBUG_DIV
//   std::cout << "revise bounds of " << scope[Z].get_domain() << " = " 
// 	    << scope[X].get_domain() << "/" << scope[Y].get_domain() 
// 	    << " = [" << min_neg[Z] << ".." << max_neg[Z] << "|" 
// 	    << (zero[Z] ? "{0}" : "_") << "|" 
// 	    << min_pos[Z] << ".." << max_pos[Z] << "]" 
// 	    << std::endl; 
// #endif

  
//   //int lb_pos=1, ub_pos=+INFTY, lb_neg=-INFTY, ub_neg=-1, lb_aux, ub_aux;
  
//   // we start with all bounds at their previous values, and update them if necessary
//   // (in which case we set up the pruning flag)
//   int lb_pos=min_pos[Z], ub_pos=max_pos[Z], lb_neg=min_neg[Z], ub_neg=max_neg[Z], 
//     nlb1, nlb2, nub1, nub2;
//     //lb_aux, ub_aux;
//   bool pruning_flag = false, pzero = false, ppos = max_pos[Z]<=0, pneg = min_neg[Z]>=0;
  
//   // // first rule: if X can be 0, and Y can be 0, then Z can be anything
//   // if(!zero[X]) { // if X cannot be 0, then neither can Y nor Z
//   //   if(zero[Z]) {
//   //     //zero[Z] = 0;
//   //     pzero = true;
//   //     pruning_flag = true;
//   //   }
//   // } else { 
//   //   if(!min_neg[X] && !max_pos[X] && !zero[Y]) {
//   //     // if X must be 0 and Y cannot be 0, then Z must be 0.
//   //     if(lb_neg<0 || ub_pos>0) {
//   // 	lb_neg = 0;
//   // 	ub_pos = 0;
//   // 	pruning_flag = true;
//   //     }
//   //   } 
//   // }

//   //	std::cout << pruning_flag << std::endl;

//   if(lb_neg != ub_pos && (!zero[X] || !zero[Y])) { // if X and Y can both be 0, we cannot deduce anything
//     if(IS_OK(wiped)) {
//       if(max_pos[Z]>0) {
// 	// revise the positive part of Z's domain (if it has one)
// 	nlb1 = nlb2 = INFTY; //lb_neg;
// 	nub1 = nub2 = 0; //ub_neg;
	
// 	// it can either be the positive parts of X and Y:
// 	if(max_pos[X]>0 && max_pos[Y]>0) {
// 	  // compute the bounds
// 	  //std::cout << "\t   lb = " << min_pos[X] << "/" << max_pos[Y] << std::endl;
// 	  //std::cout << "\t   ub = " << max_pos[X] << "/" << min_pos[Y] << std::endl;
// 	  nlb1 = min_pos[X]/max_pos[Y];
// 	  nub1 = max_pos[X]/min_pos[Y];
// 	}

// 	// or the negative parts of X and Y:
// 	if(min_neg[X]<0 && min_neg[Y]<0) {
// 	  // compute the bounds
// 	  nlb2 = max_neg[X]/min_neg[Y];
// 	  nub2 = min_neg[X]/max_neg[Y];
// 	}
// 	if(nlb1>nlb2) nlb1 = nlb2;
// 	if(nub1<nub2) nub1 = nub2;
	
// 	//std::cout << "positive bounds: [" << nlb1 << ".." << nub1 << "]" << std::endl;

// 	if(lb_pos<nlb1) {
// 	  lb_pos = nlb1;
// 	  pruning_flag = true;
// 	}
// 	if(ub_pos>nub1) {
// 	  ub_pos = nub1;
// 	  pruning_flag = true;
// 	}

// 	if(lb_pos > max_pos[Z] || ub_pos < min_pos[Z]) ppos = true;

//       } // else if(pzero || !zero[Z]) // if(lb_pos || ub_pos)
//       // 	{
//       // 	  lb_pos = min_neg[Z];
//       // 	  ub_pos = max_neg[Z];
//       // 	} else {
//       // 	lb_pos = ub_pos = 0; 
//       // }

//      if(min_neg[Z]<0) {
//        // revise the negative part of Z's domain (if it has one)
//        nlb1 = nlb2 = 0; //lb_pos;
//        nub1 = nub2 = -INFTY; //ub_pos;
	
// 	// it can either be the negitive part of X and the positive part of Y:
// 	if(min_neg[X]<0 && max_pos[Y]>0) {
// 	  // compute the bounds

// 	  nlb1 = min_neg[X]/min_pos[Y];
// 	  nub1 = max_neg[X]/max_pos[Y];
// 	}
// 	// or the negitive part of Y and the positive part of X:
// 	if(max_pos[X]>0 && min_neg[Y]<0) {
// 	  // compute the bounds
// 	  nlb2 = max_pos[X]/max_neg[Y];
// 	  nub2 = min_pos[X]/min_neg[Y];
// 	}

// 	if(nlb1>nlb2) nlb1 = nlb2;
// 	if(nub1<nub2) nub1 = nub2;

// 	//std::cout << "negative bounds: [" << nlb1 << ".." << nub1 << "]" << std::endl;
	
// 	if(lb_neg<nlb1) {
// 	  lb_neg = nlb1;
// 	  pruning_flag = true;
// 	}
// 	if(ub_neg>nub1) {
// 	  ub_neg = nub1;
// 	  pruning_flag = true;
// 	}

// 	if(lb_neg > max_neg[Z] || ub_neg < min_neg[Z]) pneg = true;

//      } // else if(pzero || !zero[Z])// if(lb_neg || ub_neg)
//      //   {
//      // 	 lb_neg = min_pos[Z];
//      // 	 ub_neg = max_pos[Z];
//      //   } else {
//      //   lb_neg = ub_neg = 0;
//      // }
//     }
//   }
  
//   //std::cout << pneg <<  (pzero || !zero[Z]) <<  ppos << std::endl;

//   if(pneg && (pzero || !zero[Z]) && ppos) {
//     wiped = FAILURE(Z);
//   } else if(pruning_flag) {
    
// #ifdef _DEBUG_DIV
//     std::cout << "set bounds to " // << scope[Z].get_domain() << " = " 
// 	      // << scope[X].get_domain() << "/" << scope[Y].get_domain() 
// 	      << "[" << lb_neg << ".." << ub_neg << "|" 
// 	      << (!pzero&&zero[Z] ? "{0}" : "_") << "|" 
// 	      << lb_pos << ".." << ub_pos << "]" << std::endl; 
// #endif
    
//     wiped = prune(lb_neg, ub_neg, lb_pos, ub_pos, pzero, Z);
//   }

//   return wiped;
// }


// Mistral::PropagationOutcome Mistral::PredicateDiv::revise_reverse_division(const int X, const int Y, const int Z) {
//   // revise the domain of Z = X/Y (Z is the integer part of X/Y)
//   Mistral::PropagationOutcome wiped = CONSISTENT;


// #ifdef _DEBUG_DIV
//   std::cout << "revise bounds of " << scope[Z].get_domain() << " = " 
// 	    << scope[X].get_domain() << "/" << scope[Y].get_domain() 
// 	    << " = [" << min_neg[Z] << ".." << max_neg[Z] << "|" 
// 	    << (zero[Z] ? "{0}" : "_") << "|" 
// 	    << min_pos[Z] << ".." << max_pos[Z] << "]" 
// 	    << std::endl; 
// #endif

  
//   //int lb_pos=1, ub_pos=+INFTY, lb_neg=-INFTY, ub_neg=-1, lb_aux, ub_aux;
  
//   // we start with all bounds at their previous values, and update them if necessary
//   // (in which case we set up the pruning flag)
//   int lb_pos=min_pos[Z], ub_pos=max_pos[Z], lb_neg=min_neg[Z], ub_neg=max_neg[Z], 
//     nlb1, nlb2, nub1, nub2, aux;
//     //lb_aux, ub_aux;
//   bool pruning_flag = false, pzero = false, ppos = max_pos[Z]<=0, pneg = min_neg[Z]>=0;
  
//   // // first rule: if X can be 0, and Y can be 0, then Z can be anything
//   // if(!zero[X]) { // if X cannot be 0, then neither can Y nor Z
//   //   if(zero[Z]) {
//   //     //zero[Z] = 0;
//   //     pzero = true;
//   //     pruning_flag = true;
//   //   }
//   // } else { 
//   //   if(!min_neg[X] && !max_pos[X] && !zero[Y]) {
//   //     // if X must be 0 and Y cannot be 0, then Z must be 0.
//   //     if(lb_neg<0 || ub_pos>0) {
//   // 	lb_neg = 0;
//   // 	ub_pos = 0;
//   // 	pruning_flag = true;
//   //     }
//   //   } 
//   // }

//   //	std::cout << pruning_flag << std::endl;

//   if(lb_neg != ub_pos && (!zero[X] || !zero[Y])) { // if X and Y can both be 0, we cannot deduce anything
//     if(IS_OK(wiped)) {
//       if(max_pos[Z]>0) {
// 	// revise the positive part of Z's domain (if it has one)
// 	nlb1 = nlb2 = INFTY; //lb_neg;
// 	nub1 = nub2 = 0; //ub_neg;
	
// 	// it can either be the positive parts of X and Y:
// 	if(max_pos[X]>0 && max_pos[Y]>0) {
// 	  // compute the bounds
// 	  //std::cout << "\t   lb = " << min_pos[X] << "/" << max_pos[Y] << std::endl;
// 	  //std::cout << "\t   ub = " << max_pos[X] << "/" << min_pos[Y] << std::endl;
// 	  //if(max_pos[Y] == 1) nlb1 = min_pos[X];
// 	  //else {
// 	  nlb1 = (min_pos[X]/(max_pos[Y]+1));
// 	  aux = ((min_pos[X]-1)/max_pos[Y]);
// 	  if(aux>nlb1) nlb1=aux;
// 	  //}
	  
// 	  if(min_pos[Y] == 1) {
// 	    nub1 = max_pos[X];
// 	  } else {
// 	    nub1 = (max_pos[X]/(min_pos[Y]-1));
// 	    aux = ((max_pos[X]+1)/min_pos[Y]);
// 	    if(aux<nub1) nub1=aux;
// 	  }
// 	}

// 	// or the negative parts of X and Y:
// 	if(min_neg[X]<0 && min_neg[Y]<0) {
// 	  // compute the bounds
// 	  nlb2 = max_neg[X]/(min_neg[Y]-1);
// 	  aux = (max_neg[X]+1)/min_neg[Y];
// 	  if(aux>nlb2) nlb2=aux;

// 	  if(max_neg[Y] == -1) {
// 	    nub2 = -min_neg[X];
// 	  } else {
// 	    nub2 = min_neg[X]/(max_neg[Y]+1);
// 	    aux = (min_neg[X]-1)/max_neg[Y];
// 	    if(aux<nub2) nub2=aux;
// 	  }
// 	}

// 	if(nlb1>nlb2) nlb1 = nlb2;
// 	if(nub1<nub2) nub1 = nub2;
	
// 	//std::cout << "positive bounds: [" << nlb1 << ".." << nub1 << "]" << std::endl;

// 	if(lb_pos<nlb1) {
// 	  lb_pos = nlb1;
// 	  pruning_flag = true;
// 	}
// 	if(ub_pos>nub1) {
// 	  ub_pos = nub1;
// 	  pruning_flag = true;
// 	}
	
// 	if(lb_pos > max_pos[Z] || ub_pos < min_pos[Z]) ppos = true;

//       } // else if(pzero || !zero[Z]) // if(lb_pos || ub_pos)
//       // 	{
//       // 	  lb_pos = min_neg[Z];
//       // 	  ub_pos = max_neg[Z];
//       // 	} else {
//       // 	lb_pos = ub_pos = 0; 
//       // }

//      if(min_neg[Z]<0) {
//        // revise the negative part of Z's domain (if it has one)
//        nlb1 = nlb2 = 0; //lb_pos;
//        nub1 = nub2 = -INFTY; //ub_pos;
	
// 	// it can either be the negitive part of X and the positive part of Y:
// 	if(min_neg[X]<0 && max_pos[Y]>0) {
// 	  // compute the bounds


// 	  nlb1 = min_neg[X]/min_pos[Y];
// 	  nub1 = max_neg[X]/max_pos[Y];
// 	}
// 	// or the negitive part of Y and the positive part of X:
// 	if(max_pos[X]>0 && min_neg[Y]<0) {
// 	  // compute the bounds
// 	  nlb2 = max_pos[X]/max_neg[Y];
// 	  nub2 = min_pos[X]/min_neg[Y];
// 	}

// 	if(nlb1>nlb2) nlb1 = nlb2;
// 	if(nub1<nub2) nub1 = nub2;

// 	//std::cout << "negative bounds: [" << nlb1 << ".." << nub1 << "]" << std::endl;
	
// 	if(lb_neg<nlb1) {
// 	  lb_neg = nlb1;
// 	  pruning_flag = true;
// 	}
// 	if(ub_neg>nub1) {
// 	  ub_neg = nub1;
// 	  pruning_flag = true;
// 	}

// 	if(lb_neg > max_neg[Z] || ub_neg < min_neg[Z]) pneg = true;

//      } // else if(pzero || !zero[Z])// if(lb_neg || ub_neg)
//      //   {
//      // 	 lb_neg = min_pos[Z];
//      // 	 ub_neg = max_pos[Z];
//      //   } else {
//      //   lb_neg = ub_neg = 0;
//      // }
//     }
//   }
  
//   //std::cout << pneg <<  (pzero || !zero[Z]) <<  ppos << std::endl;

//   if(pneg && (pzero || !zero[Z]) && ppos) {
//     wiped = FAILURE(Z);
//   } else if(pruning_flag) {
    
// #ifdef _DEBUG_DIV
//     std::cout << "set bounds to " // << scope[Z].get_domain() << " = " 
// 	      // << scope[X].get_domain() << "/" << scope[Y].get_domain() 
// 	      << "[" << lb_neg << ".." << ub_neg << "|" 
// 	      << (!pzero&&zero[Z] ? "{0}" : "_") << "|" 
// 	      << lb_pos << ".." << ub_pos << "]" << std::endl; 
// #endif
    
//     wiped = prune(lb_neg, ub_neg, lb_pos, ub_pos, pzero, Z);
//   }

//   return wiped;
// }


// Mistral::PropagationOutcome Mistral::PredicateDiv::revise_multiplication(const int X, const int Y, const int Z) {
//   // revise the domain of Z = X*Y 
//   Mistral::PropagationOutcome wiped = CONSISTENT;


// #ifdef _DEBUG_DIV
//   std::cout << "revise bounds of " << scope[Z].get_domain() << " = " 
// 	    << scope[X].get_domain() << "*" << scope[Y].get_domain() 
// 	    // << " = [" << min_neg[Z] << ".." << max_neg[Z] << "|" 
// 	    // << (zero[Z] ? "{0}" : "_") << "|" 
// 	    // << min_pos[Z] << ".." << max_pos[Z] << "]" 
// 	    << std::endl; 
// #endif

  
//   //int lb_pos=1, ub_pos=+INFTY, lb_neg=-INFTY, ub_neg=-1, lb_aux, ub_aux;
//   int lb_pos=min_pos[Z], ub_pos=max_pos[Z], lb_neg=min_neg[Z], ub_neg=max_neg[Z], 
//     nlb1, nlb2, nub1, nub2;
  
//   bool pruning_flag = false, pzero = false;
  

//   // std::cout << min_neg[X] << " " << max_pos[X] << " || " 
//   // 	    <<  min_neg[Y]  << " " << max_pos[Y] << std::endl;

//   // if X = 0 or Y = 0, then Z = 0
//   if( zero[Z] &&
//       !zero[X] && !zero[Y]) {
//     //zero[Z] = 0;
//     pzero = true;
//     pruning_flag = true;
//   }
//   else if((!min_neg[X] && !max_pos[X]) || (!min_neg[Y] && !max_pos[Y])) { 
//     lb_neg = 0;
//     ub_pos = 0;
//     pruning_flag = true;
//   }

//   // std::cout << lb_neg << "\\" << ub_pos << std::endl;

//   if(lb_neg != ub_pos) { 
//     if(IS_OK(wiped)) {
//       if(max_pos[Z]>0) {
// 	// revise the positive part of Z's domain (if it has one)
// 	//ub_pos = 0;
// 	//lb_pos = 1;
// 	//lb_aux = 1;
// 	nlb1 = nlb2 = INFTY; //lb_neg;
// 	nub1 = nub2 = 0; //ub_neg;

// 	// it can either be the positive parts of X and Y:
// 	if(max_pos[X]>0 && max_pos[Y]>0) {
// 	  // compute the bounds
// 	  // ub_pos = max_pos[X] * max_pos[Y];
// 	  // lb_pos = min_pos[X] * min_pos[Y];
// 	  nub1 = max_pos[X] * max_pos[Y];
// 	  nlb1 = min_pos[X] * min_pos[Y];
// 	}
// 	// or the negative parts of X and Y:
// 	if(min_neg[X]<0 && min_neg[Y]<0) {
// 	  // compute the bounds
// 	  // ub_aux = min_neg[X] * min_neg[Y];
// 	  // lb_aux = max_neg[X] * max_neg[Y];
// 	  nub2 = min_neg[X] * min_neg[Y];
// 	  nlb2 = max_neg[X] * max_neg[Y];
// 	}
// 	// if(lb_pos>lb_aux) lb_pos = lb_aux;
// 	// if(ub_pos<ub_aux) ub_pos = ub_aux;
// 	if(nlb1>nlb2) nlb1 = nlb2;
// 	if(nub1<nub2) nub1 = nub2;
	
// 	//std::cout << "positive bounds: [" << nlb1 << ".." << nub1 << "]" << std::endl;

// 	if(lb_pos<nlb1) {
// 	  lb_pos = nlb1;
// 	  pruning_flag = true;
// 	}
// 	if(ub_pos>nub1) {
// 	  ub_pos = nub1;
// 	  pruning_flag = true;
// 	}
//       } else if(pzero || !zero[Z]) {
// 	lb_pos = min_neg[Z];
// 	ub_pos = max_neg[Z];
//       } else {
// 	lb_pos = ub_pos = 0;
//       }

//      if(min_neg[Z]<0) {
//        // revise the negative part of Z's domain (if it has one)
//        //ub_neg = -1;
//        //lb_neg = 0;
//        //ub_aux = -1;
//        nlb1 = nlb2 = 0; //lb_pos;
//        nub1 = nub2 = -INFTY; //ub_pos;
	
// 	// it can either be the negitive part of X and the positive part of Y:
// 	if(min_neg[X]<0 && max_pos[Y]>0) {
// 	  // compute the bounds
// 	  // ub_neg = max_neg[X] * min_pos[Y];
// 	  // lb_neg = min_neg[X] * max_pos[Y];
// 	  nub1 = max_neg[X] * min_pos[Y];
// 	  nlb1 = min_neg[X] * max_pos[Y];

// 	  // std::cout << "\t  1 ub = " << max_neg[X] << "*" << min_pos[Y] << std::endl;
// 	  // std::cout << "\t  1 lb = " << min_neg[X] << "*" << max_pos[Y] << std::endl;

// 	}
// 	// or the negitive part of Y and the positive part of X:
// 	if(max_pos[X]>0 && min_neg[Y]<0) {
// 	  // compute the bounds
// 	  // ub_aux = max_neg[Y] * min_pos[X];
// 	  // lb_aux = min_neg[Y] * max_pos[X];
// 	  nub2 = max_neg[Y] * min_pos[X];
// 	  nlb2 = min_neg[Y] * max_pos[X];

// 	  // std::cout << "\t  2 ub = " << max_neg[Y] << "*" << min_pos[X] << std::endl;
// 	  // std::cout << "\t  2 lb = " << min_neg[Y] << "*" << max_pos[X] << std::endl;
// 	}
// 	// if(lb_neg>lb_aux) lb_neg = lb_aux;
// 	// if(ub_neg<ub_aux) ub_neg = ub_aux;

// 	if(nlb1>nlb2) nlb1 = nlb2;
// 	if(nub1<nub2) nub1 = nub2;


// 	  // std::cout << "change:" << std::endl;
// 	  // std::cout << "\t lbn = " << nlb1 << std::endl;
// 	  // std::cout << "\t ubn = " << nub1 << std::endl;


// 	//std::cout << "negative bounds: [" << nlb1 << ".." << nub1 << "]" << std::endl;
	
// 	if(lb_neg<nlb1) {
// 	  lb_neg = nlb1;
// 	  pruning_flag = true;
// 	}
// 	if(ub_neg>nub1) {
// 	  ub_neg = nub1;
// 	  pruning_flag = true;
// 	}
//      }  else if(pzero || !zero[Z]) {
//        lb_neg = min_pos[Z];
//        ub_neg = max_pos[Z];
//      } else {
//        lb_neg = ub_neg = 0;
//      }
//     }
//   }

//   if(pruning_flag) {
// #ifdef _DEBUG_DIV
//   std::cout << "set bounds to " // << scope[Z].get_domain() << " = " 
// 	    // << scope[X].get_domain() << "*" << scope[Y].get_domain() 
// 	    << "[" << lb_neg << ".." << ub_neg << "|" 
// 	    << (zero[Z] ? "{0}" : "_") << "|" 
// 	    << lb_pos << ".." << ub_pos << "]" << std::endl; 
// #endif
  
//     wiped = prune(lb_neg, ub_neg, lb_pos, ub_pos, pzero, Z);
//   }
  
//   return wiped;
// }


// Mistral::PropagationOutcome Mistral::PredicateDiv::prune(const int lb_neg, 
// 							 const int ub_neg, 
// 							 const int lb_pos, 
// 							 const int ub_pos,
// 							 const bool pzero,
// 							 const int Z) {


// // std::cout << changes << std::endl;
//   Event evt;
//   Mistral::PropagationOutcome wiped = CONSISTENT;
  
//   // int lb = lb_neg;
//   // if(lb>=0) lb = lb_pos;
//   // int ub = ub_pos;
//   // if(ub>=0) lb = lb_pos;

//   if(ub_pos < lb_neg) wiped = FAILURE(Z);
//   else {
//     if(lb_neg>min_neg[Z]) {

//       // std::cout  << lb_neg << ">" << min_neg[Z] 
//       // 		 << " -> update lb" << std::endl;

//       evt = scope[Z].set_min( lb_neg );
//       if( IS_FAIL(evt) ) wiped = FAILURE(Z);
//       else {
// 	if(changes.contain(Z)) {
// 	  event_type[Z] |= evt;
// 	} else {
// 	  event_type[Z] = evt;
// 	  changes.add(Z);
// 	}
// 	min_neg[Z] = scope[Z].get_min();
//       }
//     }
//     if(IS_OK(wiped) && ub_pos<max_pos[Z]) {
      
//       // std::cout  << ub_pos << ">" << max_pos[Z] 
//       // 		 << " -> update ub" << std::endl;

//       evt = scope[Z].set_max( ub_pos );
//       if( IS_FAIL(evt) ) wiped = FAILURE(Z);
//       else {
// 	if(changes.contain(Z)) {
// 	  event_type[Z] |= evt;
// 	} else {
// 	  event_type[Z] = evt;
// 	  changes.add(Z);
// 	}
// 	max_pos[Z] = scope[Z].get_max();
//       }
//     }
//     if(IS_OK(wiped) && (lb_pos>=min_neg[Z] || ub_neg<=max_pos[Z])) { 
      
//       // std::cout << lb_pos << ">=" << min_neg[Z] 
//       // 		<< " or " << ub_neg << "<=" << max_pos[Z]
//       // 		<< " -> may update inbounds" << std::endl;
      
//       if(lb_pos-1>ub_neg && (pzero || (!zero[Z] && (min_pos[Z]<lb_pos || max_neg[Z]>ub_neg)))) {
	  
// 	// std::cout  << lb_pos << ">" << min_pos[Z] 
// 	//  	   << " or " << ub_neg << "<" << max_neg[Z]
// 	// 	   << " -> update inbounds" << std::endl;
	
// 	evt = scope[Z].remove_interval(ub_neg+1, lb_pos-1);
// 	if( IS_FAIL(evt) ) wiped = FAILURE(Z);
// 	else {
// 	  if(changes.contain(Z)) {
// 	    event_type[Z] |= evt;
// 	  } else {
// 	    event_type[Z] = evt;
// 	    changes.add(Z);
// 	  }
// 	  zero[Z] = 0;
// 	  min_pos[Z] = scope[Z].get_min_pos();
// 	  max_neg[Z] = scope[Z].get_max_neg();
// 	}
//       } else {
// 	if(lb_pos>1 && min_pos[Z]<lb_pos) {


// 	  // std::cout << lb_pos << ">" << min_pos[Z] 
// 	  // 	    << " -> update negative ub" << std::endl;

// 	  evt = scope[Z].remove_interval(1, lb_pos-1);
// 	  if( IS_FAIL(evt) ) wiped = FAILURE(Z);
// 	  else {
// 	    //min_pos[Z] = lb_pos;
// 	    //max_neg[Z] = ub_neg;
// 	    if(changes.contain(Z)) {
// 	      event_type[Z] |= evt;
// 	    } else {
// 	      event_type[Z] = evt;
// 	      changes.add(Z);
// 	    }
// 	    min_pos[Z] = scope[Z].get_min_pos();
// 	  }
// 	}
// 	if(ub_neg<-1 && max_neg[Z]>ub_neg) {

// 	  // std::cout 
// 	  //   << ub_neg << "<" << max_neg[Z]
// 	  //   << " -> update positive lb" << std::endl;

// 	  evt = scope[Z].remove_interval(ub_neg+1, -1);
// 	  if( IS_FAIL(evt) ) wiped = FAILURE(Z);
// 	  else {
// 	    //min_pos[Z] = lb_pos;
// 	    //max_neg[Z] = ub_neg;
// 	    if(changes.contain(Z)) {
// 	      event_type[Z] |= evt;
// 	    } else {
// 	      event_type[Z] = evt;
// 	      changes.add(Z);
// 	    }
// 	    max_neg[Z] = scope[Z].get_max_neg();
// 	  }
// 	}
//       }
//     }
//   }

//   if((min_neg[Z]>0 && min_neg[Z]<min_pos[Z]) ||
//      (!min_neg[Z] && !zero[Z]))
//     min_neg[Z] = min_pos[Z];

//   if((max_pos[Z]<0 && max_pos[Z]>max_neg[Z]) ||
//      (!max_pos[Z] && !zero[Z]))
//     max_pos[Z] = max_neg[Z];




//   //std::cout << changes << std::endl;
// #ifdef _DEBUG_DIV
// 	if(IS_OK(wiped)) {
// 	  std::cout << scope[0].get_domain() 
// 		    << " * " << scope[1].get_domain() 
// 		    << " = " << scope[2].get_domain() << std::endl;
// 	} else std::cout << "FAIL!" << std::endl ;

// 	std::cout
// 	  << " now in [" << min_neg[Z] << ".." << max_neg[Z] << "|" 
// 	  << (zero[Z] ? "{0}" : "_") << "|" 
// 	  << min_pos[Z] << ".." << max_pos[Z] << "]" << std::endl;

//   std::cout << std::endl; 
// #endif


//   return wiped;
// }





// Mistral::PropagationOutcome Mistral::PredicateDiv::propagate() {      
//   Mistral::PropagationOutcome wiped = CONSISTENT;
  

// #ifdef _DEBUG_DIV
//   std::cout << std::endl << std::endl << "propagate " // << this 
// 	    << std::endl;
// #endif

//   for(int i=0; i<3; ++i) {
//     max_pos[i] = scope[i].get_max();
//     min_neg[i] = scope[i].get_min();

//     zero[i] = scope[i].contain(0);

//     if(min_neg[i]<0)
//       max_neg[i] = scope[i].get_max_neg();
//     else
//       max_neg[i] = min_neg[i];

//     if(max_pos[i]>0)
//       min_pos[i] = scope[i].get_min_pos();
//     else
//       min_pos[i] = max_pos[i];

//   }

 
  
// #ifdef _DEBUG_DIV
//   std::cout << scope[0]/*.get_var()*/ << " in " << scope[0].get_domain() 
// 	    << " / " << scope[1]/*.get_var()*/ << " in " << scope[1].get_domain() 
//     	    << " = " << scope[2]/*.get_var()*/ << " in " << scope[2].get_domain() << std::endl;
// #endif

//   //int i, j, lb, ub, evt_idx, rev_idx, aux_idx;
//   int evt_idx;
//   //Event evt;

//   if(zero[1] && scope[1].remove(0) == FAIL_EVENT) {
//     wiped = FAILURE(1);
//   }


//   while(IS_OK(wiped) && !changes.empty()) {

//     evt_idx = changes.pop();

// #ifdef _DEBUG_DIV
//     std::cout << "react to " << scope[evt_idx]/*.get_var()*/ << " in " 
// 	      << scope[evt_idx].get_domain() 
// 	      << (LB_CHANGED(event_type[evt_idx]) ? " (change on LB) " : "")
// 	      << (UB_CHANGED(event_type[evt_idx]) ? " (change on UB) " : "")
// 	      << std::endl;
// #endif

//     // x0 / x1 = x2 
//     if(evt_idx < 2) {
//       // we update x2 = x0/x1
//       wiped = revise_integer_division(0, 1, 2);
//       // revise x(1-evt)
//       if(IS_OK(wiped)) {
// 	if(evt_idx) // revise x0
// 	  wiped = revise_multiplication(1, 2, 0);
// 	else // revise x1 (x2 = x0/x1) -> x1 = x0/x2
// 	  wiped = revise_integer_division(0, 2, 1);
//       }
//     } else {
//       // update x0 = x1*x2 
//       wiped = revise_multiplication(1, 2, 0);
//       // update x1 = x0/x2
//       if(IS_OK(wiped)) wiped = revise_integer_division(0, 2, 1);
//     }
//   }
    
//   return wiped;
// }
  
// std::ostream& Mistral::PredicateDiv::display(std::ostream& os) const {
//   os << scope[2] << " == (" << scope[0] << " / " << scope[1] << ")";
//   return os;
// }




Mistral::ConstraintBoolSumEqual::ConstraintBoolSumEqual(Vector< Variable >& scp, const int t)
  : GlobalConstraint(scp) { 
  priority = 1;
  total = t; 
}

void Mistral::ConstraintBoolSumEqual::initialise() {
  ConstraintImplementation::initialise();
  for(unsigned int i=0; i<scope.size; ++i)
    trigger_on(_VALUE_, scope[i]);
  GlobalConstraint::initialise();
}

Mistral::ConstraintBoolSumEqual::~ConstraintBoolSumEqual() 
{
#ifdef _DEBUG_MEMORY
  std::cout << "c delete boolsumequal constraint" << std::endl;
#endif
}

Mistral::PropagationOutcome Mistral::ConstraintBoolSumEqual::propagate() 
{
  int _min_ = 0;
  int _max_ = 0;
  unsigned int i;

  for( i=0; i<scope.size; ++i ) {
    _min_ += scope[i].get_min();
    _max_ += scope[i].get_max();
  }
  if(_min_ > total || _max_ < total) return FAILURE(active[0]);
  else if(_min_ == total ) {
    for( i=0; i<scope.size; ++i ) 
      if( !scope[i].is_ground() ) scope[i].set_domain(0);
  } else if(_max_ == total ) {
    for( i=0; i<scope.size; ++i ) 
      if( !scope[i].is_ground() ) scope[i].set_domain(1);
  }
  return CONSISTENT;
}

int Mistral::ConstraintBoolSumEqual::check( const int* s ) const 
{
  int i=scope.size, t=0;
  while(--i) t+=s[i];
  return total != t; 
}

std::ostream& Mistral::ConstraintBoolSumEqual::display(std::ostream& os) const {
  os << "(" << scope[0]/*.get_var()*/ ;
  for(unsigned int i=1; i<scope.size; ++i) 
    os << " + " << scope[i]/*.get_var()*/;
  os << ") == " << total ;
  return os;
}


Mistral::ConstraintBoolSumInterval::ConstraintBoolSumInterval(Vector< Variable >& scp, const int l, const int u)
  : GlobalConstraint(scp) { 
  priority = 1;
  lb = l; 
  ub = u; 
}

void Mistral::ConstraintBoolSumInterval::initialise() {
  ConstraintImplementation::initialise();
  for(unsigned int i=0; i<scope.size; ++i)
    trigger_on(_VALUE_, scope[i]);
  GlobalConstraint::initialise();
}

Mistral::ConstraintBoolSumInterval::~ConstraintBoolSumInterval() 
{
#ifdef _DEBUG_MEMORY
  std::cout << "c delete boolsuminterval constraint" << std::endl;
#endif
}

Mistral::PropagationOutcome Mistral::ConstraintBoolSumInterval::propagate() 
{
  int _min_ = 0;
  int _max_ = 0;
  unsigned int i;

  for( i=0; i<scope.size; ++i ) {
    _min_ += scope[i].get_min();
    _max_ += scope[i].get_max();
  }
  if(_min_ > ub || _max_ < lb) return FAILURE(active[0]);
  else if(_min_ == ub ) {
    for( i=0; i<scope.size; ++i ) 
      if( !scope[i].is_ground() ) scope[i].set_domain(0);
  } else if(_max_ == lb ) {
    for( i=0; i<scope.size; ++i ) 
      if( !scope[i].is_ground() ) scope[i].set_domain(1);
  }
  return CONSISTENT;
}

int Mistral::ConstraintBoolSumInterval::check( const int* s ) const 
{
  int i=scope.size, t=0;
  while(--i) t+=s[i];
  return (t<lb || t>ub); 
}

std::ostream& Mistral::ConstraintBoolSumInterval::display(std::ostream& os) const {
  os << "(" << scope[0]/*.get_var()*/ ;
  for(unsigned int i=1; i<scope.size; ++i) 
    os << " + " << scope[i]/*.get_var()*/;
  os << ") in [" << lb << "," << ub << "]" ;
  return os;
}



Mistral::PredicateBoolSum::PredicateBoolSum(Vector< Variable >& scp, Variable tot)
  : GlobalConstraint(scp) { 
  scope.add(tot);
  priority = 1;
}

Mistral::PredicateBoolSum::PredicateBoolSum(std::vector< Variable >& scp, Variable tot)
  : GlobalConstraint(scp) {
  scope.add(tot);
  priority = 1;
}

Mistral::PredicateBoolSum::PredicateBoolSum(Vector< Variable >& scp)
  : GlobalConstraint(scp) { 
  priority = 1;
}

Mistral::PredicateBoolSum::PredicateBoolSum(std::vector< Variable >& scp)
  : GlobalConstraint(scp) { 
  priority = 1;
}

void Mistral::PredicateBoolSum::initialise() {
  ConstraintImplementation::initialise();
  for(unsigned int i=0; i<scope.size-1; ++i) {
    trigger_on(_VALUE_, scope[i]);
  }
  trigger_on(_RANGE_, scope[scope.size-1]);
  GlobalConstraint::initialise();
}

Mistral::PredicateBoolSum::~PredicateBoolSum() 
{
#ifdef _DEBUG_MEMORY
  std::cout << "c delete boolsum predicate" << std::endl;
#endif
}

Mistral::PropagationOutcome Mistral::PredicateBoolSum::propagate() 
{
  int _min_ = 0;
  int _max_ = 0;
  unsigned int // _min_=0, _max_=0,
    i, n=scope.size-1;
  int lb, ub;

  for( i=0; i<n; ++i ) {
    _min_ += scope[i].get_min();
    _max_ += scope[i].get_max();
  }
  if(_min_ > scope[n].get_max() || _max_ < scope[n].get_min()) return FAILURE(n);
  else {
    scope[n].set_min(_min_);
    scope[n].set_max(_max_);
    lb = scope[n].get_min();
    ub = scope[n].get_max();

    if(_min_ == ub) {
      for( i=0; i<n; ++i ) 
	if( !scope[i].is_ground() ) scope[i].set_domain(0);
    } else if(_max_ == lb) {
      for( i=0; i<n; ++i ) 
	if( !scope[i].is_ground() ) scope[i].set_domain(1);
    }
  }
  return CONSISTENT;
}

int Mistral::PredicateBoolSum::check( const int* s ) const 
{
  int i=scope.size-1, t=0;
  while(i--) t+=s[i];
  return (t != s[scope.size-1]); 
}

std::ostream& Mistral::PredicateBoolSum::display(std::ostream& os) const {
  os << "(" << scope[0]/*.get_var()*/ ;
  for(unsigned int i=1; i<scope.size-1; ++i) 
    os << " + " << scope[i]/*.get_var()*/;
  os << ") = " << scope[scope.size-1];
  return os;
}



Mistral::PredicateWeightedSum::PredicateWeightedSum(Vector< Variable >& scp, 
						    const int L, const int U)
  : GlobalConstraint(scp), lower_bound(L), upper_bound(U) { 
  priority = 1;
  for(unsigned int i=0; i<scope.size; ++i) {
    weight.add(1);
  }
}

Mistral::PredicateWeightedSum::PredicateWeightedSum(Vector< Variable >& scp, 
						    Vector< int >& wgt,
						    const int L, const int U)
  : GlobalConstraint(scp), lower_bound(L), upper_bound(U) { 
  priority = 1;
  for(unsigned int i=0; i<scope.size; ++i) {
    weight.add(wgt[i]);
  }
}

Mistral::PredicateWeightedSum::PredicateWeightedSum(std::vector< Variable >& scp, 
						    std::vector< int >& wgt,
						    const int L, const int U)
  : GlobalConstraint(scp), lower_bound(L), upper_bound(U) { 
  priority = 1;
  for(unsigned int i=0; i<scope.size; ++i) {
    weight.add(wgt[i]);
  }
}

void Mistral::PredicateWeightedSum::initialise() {
  ConstraintImplementation::initialise();
  //set_idempotent(true);
  //set_idempotent(false);

  wpos = 0;
  wneg = weight.size;

  // display(std::cout);
  // std::cout << std::endl;

  
  int aux_i;
  Variable aux_v;

  for(int i=0; i<wneg; ++i) {

    // std::cout << weight << std::endl;

    if(weight[i] == 1) { // swap i with wpos and increment wpos
      if(i>wpos) {
	weight[i] = weight[wpos];
	weight[wpos] = 1;
	
	aux_v = scope[i];
	scope[i] = scope[wpos];
	scope[wpos] = aux_v;

	--i;
      }
      ++wpos;
    } else if(weight[i] < 0) { // decrement wneg and swap i with wneg 
      --wneg;

      aux_i = weight[i];
      weight[i] = weight[wneg];
      weight[wneg] = aux_i;

      aux_v = scope[i];
      scope[i] = scope[wneg];
      scope[wneg] = aux_v;

      --i;
    }
  }

  for(unsigned int i=0; i<scope.size; ++i)
    trigger_on(_RANGE_, scope[i]);
  GlobalConstraint::initialise();



  lo_bound = new int[scope.size];
  up_bound = new int[scope.size];
  span = new int[scope.size];

  // //std::cout << (int*)env << std::endl;
  // std::cout  << "-- " << (int*)solver << std::endl;

  // exit(1);


  unknown_parity.initialise(solver, 0, scope.size-1, true);
  //parity.Reversible::initialise(scope[0].get_solver());
  parity.initialise(solver, ((lower_bound%2)!=0));


  //std::cout << "--parity-- " << parity << std::endl;


  for(int i=0; i<wpos; ++i) {
    if( scope[i].is_ground() ) {
      // the parity  only if the only one val is odd
      if( scope[i].get_min()%2 ) parity = 1-parity;
      unknown_parity.reversible_remove(i);
    }
  }

  for(unsigned int i=wpos; i<scope.size; ++i) {
    if( weight[i]%2 == 0 )
      unknown_parity.reversible_remove(i);
    else if( scope[i].is_ground() ) {
      unknown_parity.reversible_remove(i);
      if( scope[i].get_min()%2 ) parity = 1-parity;
    }
  }

   // display(std::cout);
   // std::cout << std::endl;

  // exit(1);

}

// void Mistral::PredicateWeightedSum::change_weight(const int i, const int w) {
//   if(weight[i]%2 =w%2))

//   weight[i] = w;
// }

Mistral::PredicateWeightedSum::~PredicateWeightedSum() 
{
#ifdef _DEBUG_MEMORY
  std::cout << "c delete weightedsum constraint" << std::endl;
#endif
  delete [] lo_bound;
  delete [] up_bound;
  delete [] span;
}




Mistral::PropagationOutcome Mistral::PredicateWeightedSum::rewrite() {

#ifdef _DEBUG_REWRITE
      std::cout << "REWRITE SUM " ;
      display(std::cout);
      std::cout << std::endl;
#endif

  RewritingOutcome r_evt = NO_EVENT; 


  //std::cout << scope.size << " " << wpos << " " << wneg << " " << lower_bound << " " << upper_bound << std::endl;

  // check if it can be rewritten as an ADD predicate
  if(scope.size == 3 && wpos == 1 && wneg == 1 && lower_bound == 0 && upper_bound == 0 ) {

    //std::cout <<  "RELAX" << std::endl;


    r_evt = SUPPRESSED;
    relax();
    get_solver()->add(Constraint(new PredicateAdd(scope[1], scope[2], scope[0])));
  }


  else if(scope.size == 2 && wpos == 1 && wneg == 1) {
    if(lower_bound == upper_bound) {
      r_evt = SUPPRESSED;
      relax();
      if(lower_bound == 0) {
	get_solver()->add(Constraint(new ConstraintEqual(scope[0], scope[1])));
      } else {
	get_solver()->add(Constraint(new PredicateOffset(scope[1], scope[0], lower_bound)));
      } 
    } else {
      if(upper_bound == INFTY) {
	r_evt = SUPPRESSED;
	relax();
	get_solver()->add(Constraint(new ConstraintLess(scope[1], scope[0], -lower_bound)));
      }
    }
  }
  
  return r_evt;
}



Mistral::PropagationOutcome Mistral::PredicateWeightedSum::propagate() 
{
  
  int i, j;
  // compute the max and th min
  int tmin, smin=0, tmax, smax=0// , maxspan=0
    , arity=scope.size;
  PropagationOutcome wiped = CONSISTENT;
  
#ifdef _DEBUG_WEIGHTEDSUM
  if(_DEBUG_WEIGHTEDSUM) {
    std::cout << std::endl << "propagate " << lower_bound << " <= " ;
    for(i=0; i<arity; ++i) {
      std::cout << " " << weight[i] << scope[i] << ":" << scope[i].get_domain();
    }
    std::cout << " <= " << upper_bound << std::endl << changes << std::endl;
  }
#endif
  
  for(i=0; i<wpos; ++i) {
    smax += (up_bound[i] = scope[i].get_max());
    smin += (lo_bound[i] = scope[i].get_min());
    span[i] = (up_bound[i]-lo_bound[i]);

#ifdef _DEBUG_WEIGHTEDSUM
    if(_DEBUG_WEIGHTEDSUM) {
      if(i)
	std::cout << " + [" << lo_bound[i] << "," << up_bound[i] << "] = [" << smin << "," << smax << "] ";
      else
	std::cout << "[" << smin << "," << smax << "] ";
    }
#endif
    
  }
  for(i=wpos; i<wneg; ++i) {
    smax += weight[i] * (up_bound[i] = scope[i].get_max());
    smin += weight[i] * (lo_bound[i] = scope[i].get_min());
    span[i] = weight[i] * (up_bound[i]-lo_bound[i]);
    
#ifdef _DEBUG_WEIGHTEDSUM
    if(_DEBUG_WEIGHTEDSUM) {
      if(i)
	std::cout << " + [" << lo_bound[i] << "," << up_bound[i] << "] = [" << smin << "," << smax << "] ";
      else
	std::cout << "[" << smin << "," << smax << "] ";
    }
#endif
    
  }
  for(i=wneg; i<arity; ++i) {
    smax += weight[i] * (lo_bound[i] = scope[i].get_min());
    smin += weight[i] * (up_bound[i] = scope[i].get_max());
    span[i] = weight[i] * (lo_bound[i]-up_bound[i]);
    
#ifdef _DEBUG_WEIGHTEDSUM
    if(_DEBUG_WEIGHTEDSUM) {
      if(i)
	std::cout << " + [" << lo_bound[i] << "," << up_bound[i] << "] = [" << smin << "," << smax << "] ";
      else
	std::cout << "[" << smin << "," << smax << "] ";
    }
#endif
  }

  
  while(IS_OK(wiped) && !events.empty()) {

#ifdef _DEBUG_WEIGHTEDSUM
    if(_DEBUG_WEIGHTEDSUM) {
      std::cout << "\nprocessing events: " << events << std::endl;
    }
#endif

    if(lower_bound == upper_bound) {
      j = events.size;
      while( j-- ) {
	i = events[j];

	//std::cout << i << ": " << (span[i]) << " " << (unknown_parity.contain(i)) << std::endl;


	if(span[i]==0 && unknown_parity.contain(i)) {
	  unknown_parity.reversible_remove(i);
	  if( lo_bound[i]%2 ) parity = 1-parity;
	}
      }
      
#ifdef _DEBUG_WEIGHTEDSUM
      if(_DEBUG_WEIGHTEDSUM) {
	display(std::cout);
	std::cout << std::endl << unknown_parity << ": " << (parity ? "odd" : "even") << std::endl;
      }
#endif   
      
      
      if(unknown_parity.size == 0) {

#ifdef _DEBUG_WEIGHTEDSUM
	if(_DEBUG_WEIGHTEDSUM) {
	  std::cout << "parity failure " << std::endl;
	}
#endif

	if(parity != 0) wiped = FAILURE(arity-1);
      } else if(unknown_parity.size == 1) { // it needs to be of parity "parity"
	i = unknown_parity[0];
	
#ifdef _DEBUG_WEIGHTEDSUM
	if(_DEBUG_WEIGHTEDSUM) {
	  std::cout << "parity pruning: " << (lo_bound[i]%2) << " " << parity << std::endl ;
	}
#endif
	
	while(IS_OK(wiped) && (lo_bound[i]%2==0) != (parity==0)) {
	  
#ifdef _DEBUG_WEIGHTEDSUM
	  if(_DEBUG_WEIGHTEDSUM) {
	    std::cout << scope[i] << " in " << scope[i].get_domain() << " => ";
	  }
#endif
	  
	  tmin = lo_bound[i];
	  if( IS_FAIL(scope[i].set_min(++lo_bound[i])) ) wiped = FAILURE(i);
	  else {
	    lo_bound[i] = scope[i].get_min();
	    if(i<wneg) smin += ((lo_bound[i] - tmin)*weight[i]);
	    else smax += ((lo_bound[i] - tmin)*weight[i]);
	  }

#ifdef _DEBUG_WEIGHTEDSUM
	  if(_DEBUG_WEIGHTEDSUM) {
	    std::cout << scope[i] << " in " << scope[i].get_domain() << std::endl;
	  }
#endif
	  
	}
	
	while(IS_OK(wiped) && (up_bound[i]%2==0) != (parity==0)) {

#ifdef _DEBUG_WEIGHTEDSUM
	  if(_DEBUG_WEIGHTEDSUM) {
	    std::cout << scope[i] << " in " << scope[i].get_domain() << " => ";
	  }
#endif
	  
	  tmin = up_bound[i];
	  if( IS_FAIL(scope[i].set_max(--up_bound[i])) ) wiped = FAILURE(i);
	  else {
	    up_bound[i] = scope[i].get_max();
	    if(i<wneg) smax -= ((tmin - up_bound[i])*weight[i]);
	    else smin -= ((tmin - up_bound[i])*weight[i]);
	  }

#ifdef _DEBUG_WEIGHTEDSUM
	  if(_DEBUG_WEIGHTEDSUM) {
	    std::cout << scope[i] << " in " << scope[i].get_domain() << std::endl;
	  }
#endif

	}
      }
    }

    if(IS_OK(wiped)) {
      
      events.clear();
      
#ifdef _DEBUG_WEIGHTEDSUM
      if(_DEBUG_WEIGHTEDSUM) {
	std::cout << " [" << smin << "," << smax << "]" << std::endl;
      }
#endif
      
      if( smax < lower_bound || smin > upper_bound ) wiped = FAILURE(arity-1);
      else {
	tmax = (smax - lower_bound);
	tmin = (upper_bound - smin);
	
	for(i=0; IS_OK(wiped) && i<wpos; ++i) {
	  
	  if( tmin < (up_bound[i] - lo_bound[i]) ) {

#ifdef _DEBUG_WEIGHTEDSUM
	  if(_DEBUG_WEIGHTEDSUM) {
	    std::cout << scope[i] << " in " << scope[i].get_domain() << " <= " << (lo_bound[i] + tmin) << std::endl;
	  }
#endif
	    
	    if(IS_FAIL(scope[i].set_max( lo_bound[i] + tmin ))) wiped = FAILURE(i);
	    else {
	      
	      events.add(i);
	      event_type[i] = UB_EVENT;
	    
	    }
	  }
	}
	
	for(i=wpos; IS_OK(wiped) && i<wneg; ++i) {
	  	  
	  if( tmin < (up_bound[i] - lo_bound[i]) * weight[i] ) {

#ifdef _DEBUG_WEIGHTEDSUM
	  if(_DEBUG_WEIGHTEDSUM) {
	    std::cout << scope[i] << " in " << scope[i].get_domain() << " <= " << (lo_bound[i] + tmin/weight[i]) << std::endl;
	  }
#endif
	  	    
	    if(IS_FAIL(scope[i].set_max( lo_bound[i] + tmin/weight[i] ))) wiped = FAILURE(i);
	    
	    else {
	      events.add(i);
	      event_type[i] = UB_EVENT;
	    }
	  }
	}
	
	for(i=wneg; IS_OK(wiped) && i<arity; ++i) {
	  
	  if( tmin < (lo_bound[i] - up_bound[i]) * weight[i] ) {

#ifdef _DEBUG_WEIGHTEDSUM
	  if(_DEBUG_WEIGHTEDSUM) {
	    std::cout << scope[i] << " in " << scope[i].get_domain() << " >= " << (up_bound[i] + tmin/weight[i]) << std::endl;
	  }
#endif
	    
	    if(IS_FAIL(scope[i].set_min( up_bound[i] + tmin/weight[i] ))) wiped = FAILURE(i);
	    
	    else {
	      events.add(i);
	      event_type[i] = LB_EVENT;
	    }
	  }
	}
	
	for(i=0; IS_OK(wiped) && i<wpos; ++i) {
	  
	  if( tmax < (up_bound[i] - lo_bound[i]) ) {
	    
#ifdef _DEBUG_WEIGHTEDSUM
	  if(_DEBUG_WEIGHTEDSUM) {
	    std::cout << scope[i] << " in " << scope[i].get_domain() << " >= " << (up_bound[i] - tmax) << std::endl;
	  }
#endif

	    if(IS_FAIL(scope[i].set_min( up_bound[i] - tmax ))) wiped = FAILURE(i);
	    
	    else {
	      
	      if(events.contain(i)) {
		event_type[i] |= LB_EVENT;
	      } else {
		events.add(i);
		event_type[i] = LB_EVENT;
	      }

	    }
	  }
	}
	
	for(i=wpos; IS_OK(wiped) && i<wneg; ++i) {
	  
	  if( tmax < (up_bound[i] - lo_bound[i]) * weight[i] ) {

#ifdef _DEBUG_WEIGHTEDSUM
	  if(_DEBUG_WEIGHTEDSUM) {
	    std::cout << scope[i] << " in " << scope[i].get_domain() << " >= " << (up_bound[i] - tmax/weight[i]) << std::endl;
	  }
#endif
	    
	    if(IS_FAIL(scope[i].set_min( up_bound[i] - tmax/weight[i] ))) wiped = FAILURE(i);
	    
	    else {
	      
	      if(events.contain(i)) {
		event_type[i] |= LB_EVENT;
	      } else {
		events.add(i);
		event_type[i] = LB_EVENT;
	      }
	      
	    }
	  }
	}


	for(i=wneg; IS_OK(wiped) && i<arity; ++i) {
	  
	  if( tmax < (lo_bound[i] - up_bound[i]) * weight[i] ) {  
	    
#ifdef _DEBUG_WEIGHTEDSUM
	  if(_DEBUG_WEIGHTEDSUM) {
	    std::cout << scope[i] << " in " << scope[i].get_domain() << " <= " << (lo_bound[i] - tmax/weight[i]) << std::endl;
	  }
#endif

	    if(IS_FAIL(scope[i].set_max( lo_bound[i] - tmax/weight[i] ))) wiped = FAILURE(i);
	    
	    else {
	      
	      if(events.contain(i)) {
		event_type[i] |= UB_EVENT;
	      } else {
		events.add(i);
		event_type[i] = UB_EVENT;
	      }

	    }
	  }
	}
      }
      
      /// update smin and smax
      for(unsigned int j=0; IS_OK(wiped) && j<events.size; ++j) {
	i = events[j];
	if(i<wpos) {
	  if(LB_CHANGED(event_type[i])){ 
	    smin -= lo_bound[i];
	    lo_bound[i] = scope[i].get_min();
	    span[i] = (up_bound[i]-lo_bound[i]);
	    smin += lo_bound[i];
	  } 
	  if(UB_CHANGED(event_type[i])){ 
	    smax -= up_bound[i];
	    up_bound[i] = scope[i].get_max();
	    span[i] = (up_bound[i]-lo_bound[i]);
	    smax += up_bound[i];
	  }
	} else if(i<wneg) {
	  if(LB_CHANGED(event_type[i])){ 
	    smin -= (lo_bound[i] * weight[i]);
	    lo_bound[i] = scope[i].get_min();
	    span[i] = weight[i] * (up_bound[i]-lo_bound[i]);
	    smin += (lo_bound[i] * weight[i]);
	  } 
	  if(UB_CHANGED(event_type[i])){ 
	    smax -= (up_bound[i] * weight[i]);
	    up_bound[i] = scope[i].get_max();
	    span[i] = weight[i] * (up_bound[i]-lo_bound[i]);
	    smax += (up_bound[i] * weight[i]);
	  }
	} else {
	  if(LB_CHANGED(event_type[i])){ 
	    smax -= (lo_bound[i] * weight[i]);
	    lo_bound[i] = scope[i].get_min();
	    span[i] = weight[i] * (lo_bound[i]-up_bound[i]);
	    smax += (lo_bound[i] * weight[i]);
	  } 
	  if(UB_CHANGED(event_type[i])){ 
	    smin -= (up_bound[i] * weight[i]);
	    up_bound[i] = scope[i].get_max();
	    span[i] = weight[i] * (lo_bound[i]-up_bound[i]);
	    smin += (up_bound[i] * weight[i]);
	  }
	}
      }
    }
  }

#ifdef _DEBUG_WEIGHTEDSUM
  if(_DEBUG_WEIGHTEDSUM) {
    std::cout << "result: ";
    for(i=0; i<arity; ++i) {
      std::cout << " " << weight[i] << scope[i] << ":" << scope[i].get_domain();
    }
    std::cout << std::endl;
  }
#endif

  return wiped;
}

int Mistral::PredicateWeightedSum::check( const int* s ) const 
{
  int i=scope.size, t=0;
  while(i--) {
    t+=(weight[i]*s[i]);
  }
  return (t < lower_bound || t > upper_bound); 
}

std::ostream& Mistral::PredicateWeightedSum::display(std::ostream& os) const {
  if(lower_bound > -INFTY) 
    os << lower_bound << " <= " ;
  os << weight[0] << "*" << scope[0]/*.get_var()*/ ;

  for(unsigned int i=1; i<scope.size; ++i) 
    os << " + " << weight[i] << "*" << scope[i]/*.get_var()*/;
  
  if(upper_bound < INFTY) 
    os << " <= " << upper_bound;
 

  return os;
}





Mistral::PredicateElement::PredicateElement(Vector< Variable >& scp, const int o)
  : GlobalConstraint(scp) {
  offset = o;
  priority = 1;
}

Mistral::PredicateElement::PredicateElement(std::vector< Variable >& scp, const int o)
  : GlobalConstraint(scp) { 
  offset = o;
  priority = 1;
}

void Mistral::PredicateElement::initialise() {
  ConstraintImplementation::initialise();

  int n = scope.size-1;
  aux_dom.initialise( std::min( 0, scope[n].get_min() ), 
		      std::max( n, scope[n].get_max() ), 
		      BitSet::empt );

  for(unsigned int i=0; i<scope.size; ++i)
    trigger_on(_DOMAIN_, scope[i]);
  //set_idempotent(true);

  GlobalConstraint::initialise();

  /////
  scope[n-1].set_min(0+offset);
  scope[n-1].set_max(n-2+offset);

}

Mistral::PredicateElement::~PredicateElement() 
{ 
#ifdef _DEBUG_MEMORY
  std::cout << "c delete element predicate" << std::endl;
#endif
}


Mistral::PropagationOutcome Mistral::PredicateElement::propagate() 
{

  PropagationOutcome wiped = CONSISTENT;
  int i, n = scope.size-2, evt, nxt; //, lb, ub, val;
  //Variable N = scope[n];
  //Variable V = scope[n+1];
  bool update_V = true;

#ifdef _DEBUG_ELEMENT 
  if(_DEBUG_ELEMENT) {
    std::cout << std::endl << std::endl << "X: " << scope[0].get_domain();
    for(i=1; i<n; ++i) {
      std::cout << " " << scope[i].get_domain();
    }
    std::cout << "[" << scope[n].get_domain() << "-" << offset << "] = " << scope[n+1].get_domain() << std::endl;
  } 
#endif

  while(IS_OK(wiped) && update_V) {
    update_V = false;
    while(IS_OK(wiped) && !changes.empty()) {
    
#ifdef _DEBUG_ELEMENT 
      if(_DEBUG_ELEMENT) {
	std::cout << changes << " " << events << std::endl;
      } 
#endif
      evt = changes.pop();
  
#ifdef _DEBUG_ELEMENT 
      if(_DEBUG_ELEMENT) {
	std::cout << "react to " << scope[evt] << " in " << scope[evt].get_domain() << std::endl;
      } 
#endif
      if(evt < n && scope[n].contain(evt+offset)) {
  
#ifdef _DEBUG_ELEMENT 
	if(_DEBUG_ELEMENT) {
	  std::cout << "  update " << scope[n] << " in " << scope[n].get_domain() << ": "
		    << scope[n+1] << " in " << scope[n+1].get_domain() << " inter "
		    << scope[evt] << " in " << scope[evt].get_domain() << "?" << std::endl;
	} 
#endif
	update_V = true;
	// scope[n] may change if scope[n+1] changes, or any X changes
	if( !scope[n+1].intersect(scope[evt]) ) {
	  
// 	  event_type[n] = scope[n].remove(evt+offset);
// 	  if( IS_FAIL(event_type[n]) ) {
// 	    wiped = FAILURE(n);
	    
// #ifdef _DEBUG_ELEMENT 
// 	    if(_DEBUG_ELEMENT) {
// 	      std::cout << "  => FAIL" << std::endl;
// 	    } 
// #endif	
// 	  } else if(!changes.contain(n)) {
// 	    changes.add(n);

// #ifdef _DEBUG_ELEMENT 
// 	    if(_DEBUG_ELEMENT) {
// 	      std::cout << "  => " << scope[n] << " in " << scope[n].get_domain() << std::endl;
// 	    } 
// #endif
// 	  }

	  FILTER1( n, remove(evt+offset) );


	}
      } else if(evt == n) {
	update_V = true;
	if(ASSIGNED(event_type[n])) {
	  // X may change if scope[n] becomes assigned 
	  i = scope[n].get_min()-offset;
	  
#ifdef _DEBUG_ELEMENT 
	  if(_DEBUG_ELEMENT) {
	    std::cout << scope[n+1] << " in " << scope[n+1].get_domain() << " == "
		      << scope[i] << " in " << scope[i].get_domain() << std::endl;
	  } 
#endif
	  
	  event_type[i] = scope[i].set_domain(scope[n+1]);
	  if( IS_FAIL(event_type[i]) ) { 
#ifdef _DEBUG_ELEMENT 
	    if(_DEBUG_ELEMENT) {
	      std::cout << "  => FAIL" << std::endl;
	    } 
#endif
	    wiped = FAILURE(i);
	  } else {
	    if( event_type[i] != NO_EVENT && !changes.contain(i) ) {
	      changes.add(i);
#ifdef _DEBUG_ELEMENT 
	      if(_DEBUG_ELEMENT) {
		std::cout << "  => " << scope[i] << " in " << scope[i].get_domain() << std::endl;
	      } 
#endif
	    }
// 	    event_type[n+1] = scope[n+1].set_domain(scope[i]);
// 	    if( IS_FAIL(event_type[n+1]) ) {
// #ifdef _DEBUG_ELEMENT 
// 	      if(_DEBUG_ELEMENT) {
// 		std::cout << "  => FAIL" << std::endl;
// 	      } 
// #endif
// 	      wiped = FAILURE(n+1);
// 	    } else if( event_type[n+1] != NO_EVENT && !changes.contain(n+1) ) {
// 	      changes.add(n+1);
// #ifdef _DEBUG_ELEMENT 
// 	      if(_DEBUG_ELEMENT) {
// 		std::cout << "  => " << scope[n+1] << " in " << scope[n+1].get_domain() << std::endl;
// 	      } 
// #endif
// 	    }

	    FILTER1( n+1 , set_domain(scope[i]) );

	  }
	}
      } else if(evt == n+1) {
#ifdef _DEBUG_ELEMENT 
	if(_DEBUG_ELEMENT) {
	  std::cout << "  update " << scope[n] << " in " << scope[n].get_domain() << ": "
		    << scope[n+1] << " in " << scope[n+1].get_domain() << " inter " << std::endl;
	} 
#endif
	
	if( scope[n].is_ground() ) {
#ifdef _DEBUG_ELEMENT 
	  if(_DEBUG_ELEMENT) {
	    std::cout << "  update " << scope[n+1] << " in " << scope[n+1].get_domain() << std::endl;
	  } 
#endif
	  
	  i = scope[n].get_min()-offset;
// 	  event_type[i] = scope[i].set_domain(scope[n+1]);
// 	  if( IS_FAIL(event_type[i]) ) { 
// #ifdef _DEBUG_ELEMENT 
// 	    if(_DEBUG_ELEMENT) {
// 	      std::cout << "  => FAIL" << std::endl;
// 	    } 
// #endif
// 	    wiped = FAILURE(i);
// 	  } else {
// 	    if( event_type[i] != NO_EVENT && !changes.contain(i) ) {
// 	      changes.add(i);
// #ifdef _DEBUG_ELEMENT 
// 	      if(_DEBUG_ELEMENT) {
// 		std::cout << "  => " << scope[i] << " in " << scope[i].get_domain() << std::endl;
// 	      } 
// #endif
// 	    }
// 	  }

	  FILTER1( i, set_domain(scope[n+1]) );

	}

	// nxt = N.get_min();
	// do {
	//   i = nxt;
	  
	//   std::cout << " " << i ;
	  
	//   nxt = N.next(i);
	// } while( i<nxt );
	
	// std::cout << std::endl;

	  
	nxt = scope[n].get_min();
	do {
	  i = nxt;
	  nxt = scope[n].next(i);

#ifdef _DEBUG_ELEMENT 
	  if(_DEBUG_ELEMENT) {
	    std::cout << "  " << scope[i-offset] << " in " << scope[i-offset].get_domain() << "?" ;
	  } 
#endif
	  if( !scope[n+1].intersect(scope[i-offset]) ) {
#ifdef _DEBUG_ELEMENT 
	    if(_DEBUG_ELEMENT) {
	      std::cout << " NO" << std::endl;
	    } 
#endif	 

	    //std::cout << n << " " << N << " " << N.get_domain() << " " << i << std::endl;

// 	    event_type[n] = scope[n].remove(i);
// 	    if( IS_FAIL(event_type[n]) ) {
// #ifdef _DEBUG_ELEMENT 
// 	      if(_DEBUG_ELEMENT) {
// 		std::cout << "  => FAIL" << std::endl;
// 	      } 
// #endif
// 	      wiped = FAILURE(n);
// 	    } else if( event_type[n] != NO_EVENT && !changes.contain(n) ) {
// 	      changes.add(n);
// #ifdef _DEBUG_ELEMENT 
// 	      if(_DEBUG_ELEMENT) {
// 		std::cout << "  => " << scope[n] << " in " << scope[n].get_domain() << std::endl;
// 	      } 
// #endif
// 	    }

	    FILTER1( n, remove(i) );
	    
#ifdef _DEBUG_ELEMENT 
	    if(_DEBUG_ELEMENT) {
	      if(event_type[n] == NO_EVENT) {
		std::cout << "  => NO EVENT" << std::endl;
	      }
	    } 
#endif	    
	  } 
#ifdef _DEBUG_ELEMENT 
	  
	  else if(_DEBUG_ELEMENT) {
	    std::cout << " YES" << std::endl;
	  }
#endif
	} while( i<nxt );
      }
    }
    
    if(IS_OK(wiped) && update_V) {
#ifdef _DEBUG_ELEMENT 
      if(_DEBUG_ELEMENT) {
	std::cout << "  update " << scope[n+1] << " in " << scope[n+1].get_domain() << std::endl;
      } 
#endif
      aux_dom.clear();
      //lb=INFTY; ub=-INFTY;
      
      // std::cout << "N " << n << ": " << scope[n] << " in " << scope[n].get_domain() << std::endl;
      // if(scope[n].domain_type == BITSET_VAR) {
      // 	scope[n].bitset_domain->debug_print();
      // }

      nxt = scope[n].get_min();
      do {
	i = nxt;
	// if(scope[i-offset].is_range()) {
	//   val = scope[i-offset].get_min();
	//   if(lb>val) lb = val;
	//   val = scope[i-offset].get_max();
	//   if(ub<val) ub = val;
	// } else {

	

	// std::cout << i << "-" << offset << " = " << (i-offset) << std::endl;

	// std::cout << scope << std::endl;


	scope[i-offset].union_to(aux_dom);
	//}
	nxt = scope[n].next(i);
#ifdef _DEBUG_ELEMENT 
	if(_DEBUG_ELEMENT) {
	  std::cout << " + " << scope[i-offset] << " in " << scope[i-offset].get_domain() << ": " << aux_dom << std::endl; 
	} 
#endif
      } while( i<nxt );
      // if(lb<=ub) aux_dom.fill(lb,ub);
      // std::cout << aux_dom << std::endl;
      
//       event_type[n+1] = scope[n+1].set_domain(aux_dom);
//       if( IS_FAIL(event_type[n+1]) ) { 
// #ifdef _DEBUG_ELEMENT 
// 	if(_DEBUG_ELEMENT) {
// 	  std::cout << "  => FAIL" << std::endl;
// 	} 
// #endif
// 	wiped = FAILURE(n+1);
//       } else if( event_type[n+1] != NO_EVENT && !changes.contain(n+1) ) {
// 	changes.add(n+1);
 #ifdef _DEBUG_ELEMENT 
 	if(_DEBUG_ELEMENT) {
 	  std::cout << "     " << scope[n+1] << " in " << scope[n+1].get_domain() << " = " << aux_dom << std::endl;
 	} 
 #endif
//       }

      FILTER1(n+1, set_domain(aux_dom) );

 #ifdef _DEBUG_ELEMENT 
 	if(_DEBUG_ELEMENT) {
 	  std::cout << "  => " << scope[n+1] << " in " << scope[n+1].get_domain() << std::endl;
 	} 
 #endif


    }
  }

#ifdef _DEBUG_ELEMENT 
  if(_DEBUG_ELEMENT) { 
    std::cout << "return " << wiped << std::endl;
    std::cout << "X: " << scope[0].get_domain();
    for(i=1; i<n; ++i) {
      std::cout << " " << scope[i].get_domain();
    }
    std::cout << "[" << scope[n].get_domain() << "] = " << scope[n+1].get_domain() << std::endl << std::endl ;
  } 
#endif 

  //exit(1);


  return wiped;
}

int Mistral::PredicateElement::check( const int* s ) const 
{
  return (s[s[scope.size-2]-offset] != s[scope.size-1]);
}

std::ostream& Mistral::PredicateElement::display(std::ostream& os) const {
  os << "(" << scope[0]/*.get_var()*/;
  for(unsigned int i=1; i<scope.size-2; ++i) {
    os << " " << scope[i]/*.get_var()*/;
  }
  os << ")[" << scope.back(2)/*.get_var()*/ << "] == " << scope.back()/*.get_var()*/;
  return os;
}


Mistral::ConstraintCliqueNotEqual::ConstraintCliqueNotEqual(Vector< Variable >& scp)
  : GlobalConstraint(scp) { priority = 2; }


void Mistral::ConstraintCliqueNotEqual::initialise() {
  ConstraintImplementation::initialise();
  for(unsigned int i=0; i<scope.size; ++i) {
    trigger_on(_VALUE_, scope[i]);
  }
  GlobalConstraint::initialise();
  //GlobalConstraint::set_idempotent(true);
}

void Mistral::ConstraintCliqueNotEqual::mark_domain() {
  for(unsigned int i=0; i<scope.size; ++i) {
    get_solver()->mark_non_convex(scope[i].id());
  }
}

Mistral::ConstraintCliqueNotEqual::~ConstraintCliqueNotEqual() 
{
#ifdef _DEBUG_MEMORY
  std::cout << "c delete cliquenotequal constraint" << std::endl;
#endif
}



Mistral::PropagationOutcome Mistral::ConstraintCliqueNotEqual::propagate() 
{

  unsigned int i, j, n=active.size, m;
  int active_var, value;
  Event evt;

#ifdef _DEBUG_CLIQUENOTEQUAL
  std::cout << "propagate " << this << std::endl;
  for(i=0; i<scope.size; ++i) 
    std::cout << scope[i].get_domain() << (changes.contain(i) ? "* " : " " );
  std::cout << std::endl;
#endif
  for(i=0; i<changes.size; ++i) {
    value = scope[changes[i]].get_min();

#ifdef _DEBUG_CLIQUENOTEQUAL
    std::cout << "react to " << scope[changes[i]] << "=" << value << std::endl; 
    std::cout << "  check pairs of recently assigned vars: " << std::endl;
#endif
    for(j=i+1; j<changes.size; ++j) {
#ifdef _DEBUG_CLIQUENOTEQUAL
      std::cout << "    " << scope[changes[j]] << "=" << scope[changes[j]].get_min() << " ";
#endif
      if(scope[changes[j]].get_min() == value) {
#ifdef _DEBUG_CLIQUENOTEQUAL
	std::cout << "FAIL!" << std::endl;
#endif
	return FAILURE(changes[j]);
      } 
#ifdef _DEBUG_CLIQUENOTEQUAL
      else std::cout << "OK!" << std::endl;
#endif
    }
    
#ifdef _DEBUG_CLIQUENOTEQUAL
    std::cout << "  remove this value from active vars: " << std::endl;
#endif
    // since the set of active variables might change while
    // processing this loop, we do it backward
    for(j=n; j; --j) {
      active_var = active[j-1];

#ifdef _DEBUG_CLIQUENOTEQUAL
      std::cout << "    " << scope[active_var] << "-" << value << " ";
#endif

      evt = scope[active_var].remove(value);
      if(evt == FAIL_EVENT) {
#ifdef _DEBUG_CLIQUENOTEQUAL
	std::cout << "FAIL!" << std::endl;
#endif
	//assigned.clear();
	return FAILURE(active_var);
      } else if(ASSIGNED(evt)) {
	//assigned.add(active_var);
	active.reversible_remove(active_var);
      }
#ifdef _DEBUG_CLIQUENOTEQUAL
      else std::cout << "OK!" << std::endl;
#endif
    }
  }

  /// The following is to ensure idempotency
  m = active.size;

  //std::cout << n << " " << m << std::endl;

  while(m < n) {
    for(i=m; i<n; ++i) {
      value = scope[active[i]].get_min();
#ifdef _DEBUG_CLIQUENOTEQUAL
      std::cout << "rreact to " << scope[active[i]] << "=" << value << std::endl; 
      std::cout << "  check pairs of recently assigned vars: " << std::endl;
#endif
      for(j=i+1; j<n; ++j) {
#ifdef _DEBUG_CLIQUENOTEQUAL
	std::cout << "    " << scope[changes[j]] << "=" << scope[changes[j]].get_min() << " ";
#endif
	if(scope[active[j]].get_min() == value) {
#ifdef _DEBUG_CLIQUENOTEQUAL
	  std::cout << "FAIL!" << std::endl;
#endif
	  return FAILURE(active[j]);
	}
#ifdef _DEBUG_CLIQUENOTEQUAL
	else std::cout << "OK!" << std::endl;
#endif
      }

#ifdef _DEBUG_CLIQUENOTEQUAL
      std::cout << "  remove this value from active vars: " << std::endl;
#endif

      for(j=m; j; --j) {
	active_var = active[j-1];
#ifdef _DEBUG_CLIQUENOTEQUAL
	std::cout << "    " << scope[active_var] << "-" << value << " ";
#endif
	if(scope[active_var].remove(value) == FAIL_EVENT) {
#ifdef _DEBUG_CLIQUENOTEQUAL
	std::cout << "FAIL!" << std::endl;
#endif
	  return FAILURE(active_var);
	}
#ifdef _DEBUG_CLIQUENOTEQUAL
	else std::cout << "OK!" << std::endl;
#endif
      }
    }
    n = m;
    m = active.size;
  }

  return CONSISTENT;
}

int Mistral::ConstraintCliqueNotEqual::check( const int* s ) const 
{
  int i=scope.size, j;
  while(--i) {
    j=i;
    while(j--)
      if( s[i] == s[j] ) return 1;
  }
  return 0; 
}

std::ostream& Mistral::ConstraintCliqueNotEqual::display(std::ostream& os) const {
  os << "=/=(" << scope[0]/*.get_var()*/ ;
  for(unsigned int i=1; i<scope.size; ++i) 
    os << " ," << scope[i]/*.get_var()*/;
  os << ")" ; //<< events << " " << event_type;
  return os;
}


/**********************************************
 * AllDiff Constraint 
 **********************************************/

const int INCONSISTENT = 0;
const int CHANGES      = 1;
const int NO_CHANGES   = 2;

Mistral::ConstraintAllDiff::ConstraintAllDiff(Vector< Variable >& scp)
  : GlobalConstraint(scp) { priority = 0; }

Mistral::ConstraintAllDiff::ConstraintAllDiff(std::vector< Variable >& scp)
  : GlobalConstraint(scp) { priority = 0; }


void Mistral::ConstraintAllDiff::initialise() {

  ConstraintImplementation::initialise();

  for(unsigned int i=0; i<scope.size; ++i) {
    trigger_on(_RANGE_, scope[i]);
  }

  GlobalConstraint::initialise();

  //GlobalConstraint::set_idempotent(true);
  unsigned int i;
  //level = &(scope[0].get_solver()->level);
  lastLevel = -2;
  nb = 0;

  iv        = new Interval[scope.size];
  minsorted = new Interval*[scope.size];
  maxsorted = new Interval*[scope.size];
  bounds    = new int[2*scope.size+2];
  std::fill(bounds, bounds+2*scope.size+2, 0);

  for( i=0; i<scope.size; ++i ) {
    minsorted[i] = maxsorted[i] = &iv[i];  
    iv[i].min = iv[i].max = NOVAL;
  }

  t = new int[2*scope.size+2];
  d = new int[2*scope.size+2];
  h = new int[2*scope.size+2];
}

void Mistral::ConstraintAllDiff::mark_domain() {
  for(unsigned int i=0; i<scope.size; ++i) {
    get_solver()->mark_non_convex(scope[i].id());
  }
}


Mistral::ConstraintAllDiff::~ConstraintAllDiff() 
{
#ifdef _DEBUG_MEMORY
  std::cout << "c delete alldiff constraint" << std::endl;
#endif
  delete [] bounds;
  delete [] maxsorted;
  delete [] minsorted;
  delete [] iv;
  delete [] h;
  delete [] d;
  delete [] t;
}

void sortmin( Mistral::Interval *v[], int n ) 
{
  int i, current;
  bool sorted;
  Mistral::Interval *t;

  current = n-1;
  sorted = false;
  while( !sorted ) {
    sorted = true;
    for( i = 0; i < current; i++ ) {
      if( v[i]->min > v[i+1]->min ) {
        t = v[i];
        v[i] = v[i+1];
        v[i+1] = t;
        sorted = false;
      }
    }
    current--;
  }
}

void sortmax( Mistral::Interval *v[], int n ) 
{
  int i, current;
  bool sorted;
  Mistral::Interval *t;

  current = 0;
  sorted = false;
  while( !sorted ) {
    sorted = true;
    for( i = n-1; i > current; i-- ) {
      if( v[i]->max < v[i-1]->max ) {
        t = v[i];
        v[i] = v[i-1];
        v[i-1] = t;
        sorted = false;
      }
    }
    current++;
  }
}

void Mistral::ConstraintAllDiff::sortit() 
{  
  int i,j,nb;
  int min,max,last;

  sortmin(minsorted, scope.size);
  sortmax(maxsorted, scope.size);

  min = minsorted[0]->min;
  max = maxsorted[0]->max + 1;
  bounds[0] = last = min-2;

  for (i=j=nb=0;;) { // merge minsorted[] and maxsorted[] into bounds[]
    if (i<(int)(scope.size) && min<=max) {	// make sure minsorted exhausted first
      if (min != last)
        bounds[++nb] = last = min;
      minsorted[i]->minrank = nb;
      if (++i < (int)(scope.size))
        min = minsorted[i]->min;
    } else {
      if (max != last)
	bounds[++nb] = last = max;
      maxsorted[j]->maxrank = nb;
      if (++j == (int)(scope.size)) break;
      max = maxsorted[j]->max + 1;
    }
  }
  ConstraintAllDiff::nb = nb;
  bounds[nb+1] = bounds[nb] + 2;
}


void pathset(int *t, int start, int end, int to) 
{
  int k, l;
  for (l=start; (k=l) != end; t[k]=to) 
    l = t[k];  
}

int pathmin(int *t, int i) 
{
  for (; t[i] < i; i=t[i]) ;
  return i;
}

int pathmax(int *t, int i) 
{
  for (; t[i] > i; i=t[i]) ;  
  return i;
}


int Mistral::ConstraintAllDiff::filterlower() 
{
  int i,j;
  int w,x,y,z,_changes_ = 0;

  for (i=1; i<=nb+1; i++)
    d[i] = bounds[i] - bounds[t[i]=h[i]=i-1];
  for (i=0; i<(int)(scope.size); i++) { // visit Intervals in increasing max order
    x = maxsorted[i]->minrank; y = maxsorted[i]->maxrank;
    j = t[z = pathmax(t, x+1)];
    if (--d[z] == 0)
      t[z = pathmax(t, t[z]=z+1)] = j;
    pathset(t, x+1, z, z); // path compression
    if (d[z] < bounds[z]-bounds[y]) return INCONSISTENT; // no solution
    if (h[x] > x) {
      maxsorted[i]->min = bounds[w = pathmax(h, h[x])];
      pathset(h, x, w, w); // path compression
      _changes_ = 1;
    }
    if (d[z] == bounds[z]-bounds[y]) {
      pathset(h, h[y], j-1, y); // mark hall Interval
      h[y] = j-1; //("hall Interval [%d,%d)\n",bounds[j],bounds[y]);
    }
  }
  if( _changes_ )
    return CHANGES;
  else
    return NO_CHANGES;
}


int Mistral::ConstraintAllDiff::filterupper()
{
  int i,j;
  int w,x,y,z,_changes_ = 0;

  for (i=0; i<=nb; i++)
    d[i] = bounds[t[i]=h[i]=i+1] - bounds[i];
  for (i=scope.size; --i>=0; ) { // visit Intervals in decreasing min order
    x = minsorted[i]->maxrank; y = minsorted[i]->minrank;
    j = t[z = pathmin(t, x-1)];
    if (--d[z] == 0)
      t[z = pathmin(t, t[z]=z-1)] = j;
    pathset(t, x-1, z, z);
    if (d[z] < bounds[y]-bounds[z]) return INCONSISTENT; // no solution
    if (h[x] < x) {
      minsorted[i]->max = bounds[w = pathmin(h, h[x])] - 1;
      pathset(h, x, w, w);
      _changes_ = 1;
    }
    if (d[z] == bounds[y]-bounds[z]) {
      pathset(h, h[y], j+1, y);
      h[y] = j+1;
    }
  }
  if( _changes_ )
    return CHANGES;
  else
    return NO_CHANGES;
}

Mistral::PropagationOutcome Mistral::ConstraintAllDiff::propagate() 
{
  unsigned int i, a, b;

  int status_lower, status_upper;
  int l, u;

   a = 0;
   b = scope.size;

  //if( lastLevel != ((solver->level) - 1) ) {
  if( lastLevel != ((solver->level) - 1) ) {
    // not incremental
    status_lower = CHANGES;
    status_upper = CHANGES;
    i = 0;
    while (i < scope.size) {
      iv[i].min = scope[i].get_min();
      iv[i].max = scope[i].get_max();
      i++;
    }
  }
  else {
    // incremental
    status_lower = NO_CHANGES;
    status_upper = NO_CHANGES;
    for( i = a; i < b; i++ ) {
      l = iv[i].min;
      u = iv[i].max;
      iv[i].min = scope[i].get_min();
      iv[i].max = scope[i].get_max();
      if( l != iv[i].min ) status_lower = CHANGES;
      if( u != iv[i].max ) status_upper = CHANGES;
    }
  }

//   //a = 0;
//   //b = scope.size;

//   //if( lastLevel != ((solver->level) - 1) ) {
//   if( lastLevel != ((*level) - 1) ) {
//     // not incremental
//     status_lower = CHANGES;
//     status_upper = CHANGES;
//     i = 0;
//     while (i < scope.size) {
//       iv[i].min = scope[i].get_min();
//       iv[i].max = scope[i].get_max();
//       i++;
//     }
//   }
//   else {
//     // incremental
//     status_lower = NO_CHANGES;
//     status_upper = NO_CHANGES;
//     //for( i = a; i < b; i++ ) {

//     //std::cout << "ccc " ;
//     for(unsigned int j=0; j<changes.size; ++j) {
//       i = changes[j];
// //       std::cout << i << " ";
// //     }
// //     std::cout << std::endl;

// //     std::cout << "rrr " ;
// //     for( i = a; i < b; i++ ) {

//       l = iv[i].min;
//       u = iv[i].max;
//       iv[i].min = scope[i].get_min();
//       iv[i].max = scope[i].get_max();

//       //      if( l != iv[i].min || u != iv[i].max ) std::cout << i << " ";

//       if( l != iv[i].min ) status_lower = CHANGES;
//       if( u != iv[i].max ) status_upper = CHANGES;
//     }
//     //    std::cout << std::endl << std::endl;
//   }

  //lastLevel = *level;//(solver->level);
  lastLevel = (solver->level);

  if( status_lower == NO_CHANGES && status_upper == NO_CHANGES ) 
    return CONSISTENT;

  sortit();



  status_lower = filterlower();
  if( status_lower != INCONSISTENT )
    status_upper = filterupper();  

  if( (status_lower == INCONSISTENT) || (status_upper == INCONSISTENT) ) 
    { return FAILURE(changes.back()); }
  else
    if( (status_lower == CHANGES) || (status_upper == CHANGES) ) {
      i = 0;
      while (i < scope.size) {
	if( scope[i].set_min( iv[i].min ) == FAIL_EVENT )  { return FAILURE(i); }
	if( scope[i].set_max( iv[i].max ) == FAIL_EVENT )  { return FAILURE(i); }
	i++;
      }
    }  
  return CONSISTENT;
}

int Mistral::ConstraintAllDiff::check( const int* s ) const 
{
  int i=scope.size, j;
  while(--i) {
    j=i;
    while(j--)
      if( s[i] == s[j] ) return 1;
  }
  return 0; 
}

std::ostream& Mistral::ConstraintAllDiff::display(std::ostream& os) const {
  os << "alldiff(" << scope[0]/*.get_var()*/ ;
  for(unsigned int i=1; i<scope.size; ++i) 
    os << " ," << scope[i]/*.get_var()*/;
  os << ")" ;
  return os;
}



Mistral::PredicateMin::PredicateMin(Vector< Variable >& scp) : GlobalConstraint(scp) { priority = 1; }

Mistral::PredicateMin::~PredicateMin() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete min predicate" << std::endl;
#endif
}

void Mistral::PredicateMin::initialise() {
  ConstraintImplementation::initialise();
  last_min = 0;
  for(unsigned int i=0; i<scope.size; ++i) {
    trigger_on(_RANGE_, scope[i]);
    if(scope[i].get_min() < scope[last_min].get_min())
      last_min = i;
  }
  GlobalConstraint::initialise();

  candidates.initialise(solver, 0, scope.size-2, true);
  int n = scope.size-1;
  for(int i=0; i<n; ++i) {
    if(scope[i].get_min() > scope[n].get_max()) {
      candidates.reversible_remove(i);
      //relax_from(i);
      if(!scope[i].is_ground()) un_post_from(i);
      //un_post_from(i);
    }
  }

  enforce_nfc1 = false;
  //set_idempotent(true);
}

int Mistral::PredicateMin::check( const int* sol ) const {
  int i=scope.size-1, the_min=sol[0], smin = sol[i];
  while( --i > 0 ) {
    if(sol[i] < the_min) the_min = sol[i];
    if(the_min < smin) break;
  }
  return the_min != smin;
}

Mistral::PropagationOutcome Mistral::PredicateMin::propagate() {
  PropagationOutcome wiped = CONSISTENT;
  unsigned int n = scope.size-1, evt;
  int i, val, aux, var; 


#ifdef _DEBUG_MIN
  if(_DEBUG_MIN) {
  std::cout << std::endl << active << std::endl << candidates << std::endl;
  }
#endif

  if(candidates.size != 1) {
    while(IS_OK(wiped) && !changes.empty()) {
      if(changes.contain(n)) { 

#ifdef _DEBUG_MIN
	if(_DEBUG_MIN) {
	  std::cout << " " << scope[n] << " has changed " << std::endl;
	}
#endif
	
	changes.remove(n);
	// the "self" variable has changed.
	
	// if the lower bound has been changed, we need to update the lb of every active variable
	if(LB_CHANGED(event_type[n])) {
	  val = scope[n].get_min();
	  for(i=candidates.size; IS_OK(wiped) && --i>=0; ) {
	    var = candidates[i];

	    
	    if(scope[var].get_min() < val) {
	      
#ifdef _DEBUG_MIN
	      if(_DEBUG_MIN) {
		std::cout << " => " << scope[var] << " in " << scope[var].get_domain() << " >= " << val << std::endl;
	      }
#endif
	      
	      event_type[var] |= scope[var].set_min(val);
	      if(FAILED(event_type[var])) wiped = FAILURE(var);
	      else if(!changes.contain(var)) {
		
#ifdef _DEBUG_MIN
		if(_DEBUG_MIN) {
		  std::cout << "    " << scope[var] << " has changed: " << scope[var].get_domain() << std::endl;
		}
#endif
		
		changes.add(var);
		
		if((scope[var].get_min() > scope[n].get_max()) && candidates.contain(var)) {
		  
		  candidates.reversible_remove(var);
		  
#ifdef _DEBUG_MIN
		  if(_DEBUG_MIN) {
		    std::cout << "    " << scope[var] << " can no longer be min" << std::endl;
		  }
#endif
		}
	      }	
	    }
	  }
	}
	
	// if the upper bound has been changed, some active variables may become unactive
	// The min(maximum) might also change because variables became unactive
	if(IS_OK(wiped) && UB_CHANGED(event_type[n])) {
	  aux = val = scope[n].get_max();
	
#ifdef _DEBUG_MIN
	  if(_DEBUG_MIN) {
	    std::cout << " relax the constraint from vars greater than " << val << std::endl;
	  }
#endif
	
	  for(i=candidates.size; IS_OK(wiped) && i;) {
	    --i;
	    var = candidates[i];

	    //std::cout << i << ": " << 

	    if(scope[var].get_min() > val) {
	    
#ifdef _DEBUG_MIN
	      if(_DEBUG_MIN) {
		std::cout << "    relax1 " << this << " from " << scope[var] << std::endl;
	      }
#endif 

	      //relax_from(var);
	      if(candidates.contain(var)) candidates.reversible_remove(var);
	      //std::cout << candidates << std::endl;

	    } else if(aux > scope[var].get_max()) aux = scope[var].get_max();
	  }

#ifdef _DEBUG_MIN
	  if(_DEBUG_MIN) {
	    std::cout << " => " << scope[n] << " in " << scope[n].get_domain() << " <= " << aux << std::endl;
	  }
#endif

	  if(scope[n].set_max(aux) == FAIL_EVENT) wiped = FAILURE(n);	
	}
      } 

      if(IS_OK(wiped)) {

	// store the min min of the changed vars in aux, and the min max in val
	aux = INFTY;
	val = INFTY;
	while(!changes.empty()) {
	  evt = changes.pop();
	
#ifdef _DEBUG_MIN
	  if(_DEBUG_MIN) {
	    std::cout << "-event on " << scope[evt] 
		      << (LB_CHANGED(event_type[evt]) ? " (lb) " : " ")
		      << (UB_CHANGED(event_type[evt]) ? "(ub) " : " ")
		      << std::endl;
	  }
#endif
	
	  if(LB_CHANGED(event_type[evt]) 
	     && scope[evt].get_min()<aux
	     ) {
	    aux = scope[evt].get_min();
	    if(aux > scope[n].get_max()) {

#ifdef _DEBUG_MIN
	      if(_DEBUG_MIN) {
		std::cout << "    relax2 " << this << " from " << scope[evt] << std::endl;
	      }
#endif 

	      //relax_from(evt);
	      if(candidates.contain(evt)) candidates.reversible_remove(evt);
	      //std::cout << candidates << std::endl;
	    }
	  }
	  if(UB_CHANGED(event_type[evt]) && scope[evt].get_max()<val)
	    val = scope[evt].get_max();
	  //if(ASSIGNED(event_type[evt]))
	  //candidates.reversible_remove(evt);
	}
      
#ifdef _DEBUG_MIN
	if(_DEBUG_MIN) {
	  if(val < INFTY) 
	    std::cout << " New min maximum: " << val << std::endl;
	  if(aux < INFTY) 
	    std::cout << " New min minimum: " << aux << std::endl;
	}
#endif
      
      
	if(val < INFTY) {
	
#ifdef _DEBUG_MIN
	  if(_DEBUG_MIN) {
	    std::cout << " => " << scope[n] << " in " << scope[n].get_domain() 
		      << " <= " << val << std::endl;
	  }
#endif
	
	  event_type[n] = scope[n].set_max(val);
	  if(event_type[n] == FAIL_EVENT) wiped = FAILURE(n);
	  else {
	    changes.add(n);
	  }
	}
      
	if(aux < INFTY) {
	  val = scope[n].get_min();
	
#ifdef _DEBUG_MIN
	  if(_DEBUG_MIN) {
	    std::cout << " previous min minimum witness: " << scope[last_min] 
		      << " in " << scope[last_min].get_domain() << std::endl;
	  }
#endif
      
	  if(aux > val && scope[last_min].get_min() >  val) {
	  
	    // std::cout << candidates << std::endl;
	    // for(i=0; i<scope.size; ++i) {
	    //   std::cout << " " << scope[i].get_domain();
	    // }
	    // std::cout << std::endl;

	    // look for a new witness for the min
	    i = candidates.size;
	    while( i > 0 ) {
	      --i;

	      // std::cout << i << std::endl;
	    
	      // std::cout << candidates[i] << std::endl;

	      var = scope[candidates[i]].get_min();
	      if(aux > var) aux = var;
	      if(var <= val) {
		last_min = candidates[i];
		break;
	      } 
	    }
	  
	    if(aux > val) {
#ifdef _DEBUG_MIN
	      if(_DEBUG_MIN) {
		std::cout << " => " << scope[n] << " in " << scope[n].get_domain() 
			  << " >= " << aux << std::endl;
	      }
#endif
	      event_type[n] = scope[n].set_min(aux);
	      if(event_type[n] == FAIL_EVENT) wiped = FAILURE(n);
	      else {
		if(!changes.contain(n)) changes.add(n);
	      }
	    }
	  } else {
	  
#ifdef _DEBUG_MIN
	    if(_DEBUG_MIN) {
	      std::cout << " no need to update " << scope[n] << "'s min" << std::endl;
	    }
#endif
	  
	  }
	}
      }
    }  
  }

  if(candidates.size == 1) {

#ifdef _DEBUG_MIN
    if(_DEBUG_MIN) {
      std::cout << scope[candidates[0]] << " is the last min var " << std::endl;
    }
#endif

    if(scope[candidates[0]].set_domain(scope[n]) == FAIL_EVENT) wiped = FAILURE(candidates[0]);
    if(scope[n].set_domain(scope[candidates[0]]) == FAIL_EVENT) wiped = FAILURE(n);
  }


  return wiped;
}

std::ostream& Mistral::PredicateMin::display(std::ostream& os) const {
  os << scope.back() << " == min(" << scope[0]/*.get_var()*/;
  for(unsigned int i=1; i<scope.size-1; ++i) {
    os << ", " << scope[i]/*.get_var()*/ ;
  }
  os << ")";
  return os;
}




Mistral::PredicateMax::PredicateMax(Vector< Variable >& scp) : GlobalConstraint(scp) { priority = 1; }

Mistral::PredicateMax::~PredicateMax() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete max predicate" << std::endl;
#endif
}

void Mistral::PredicateMax::initialise() {
  ConstraintImplementation::initialise();
  last_max = 0;
  for(unsigned int i=0; i<scope.size; ++i) {
    trigger_on(_RANGE_, scope[i]);
    if(scope[i].get_max() > scope[last_max].get_max())
      last_max = i;
  }

  GlobalConstraint::initialise();


//   //std::cout << id << std::endl;
// #ifdef _DEBUG_MAX
//   for(int i=0; i<scope.size-1; ++i) {
//       std::cout << scope[i] << " in " << scope[i].get_domain() << std::endl;
//   }
// #endif


  candidates.initialise(solver, 0, scope.size-2, true);
  int n = scope.size-1;
  for(int i=0; i<n; ++i) {
    if(scope[i].get_max() < scope[n].get_min()) {
      candidates.reversible_remove(i);
      //relax_from(i);
      if(!scope[i].is_ground()) un_post_from(i);
    }
  }

  enforce_nfc1 = false;
  //set_idempotent(true);
}

int Mistral::PredicateMax::check( const int* sol ) const {
  int i=scope.size-1, the_max=sol[0], smax = sol[i];
  while( --i > 0 ) {
    if(sol[i] > the_max) the_max = sol[i];
    if(the_max > smax) break;
  }
  return the_max != smax;
}

Mistral::PropagationOutcome Mistral::PredicateMax::propagate() {
  PropagationOutcome wiped = CONSISTENT;
  unsigned int n = scope.size-1, evt;
  int i, val, aux, var; 

#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
  for(i=0; i<changes.size; ++i) {
    std::cout << event2str(event_type[changes[i]]) << " " 
	      << LB_CHANGED(event_type[changes[i]]) << " " 
	      << UB_CHANGED(event_type[changes[i]]) << std::endl;
  }
  }
#endif

#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
  std::cout << "propagate max( " ;
  for(i=0; i<n; ++i) 
    std::cout << "[" << scope[i].get_min() << "," << scope[i].get_max() << "] ";
  std::cout << ") = [" << scope[n].get_min() << "," << scope[n].get_max() << "]" ;

  std::cout << std::endl << active << std::endl << candidates << " " << changes << std::endl;
  }
#endif

  if(candidates.size != 1) {
    while(IS_OK(wiped) && !changes.empty()) {
      if(changes.contain(n)) { 
	
#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
	std::cout << " " << scope[n] << " has changed " << std::endl;
  }
#endif
	
	changes.remove(n);
	// the "self" variable has changed.
	
	// if the upper bound has been changed, we need to update the ub of every active variable
	if(UB_CHANGED(event_type[n])) {
	  val = scope[n].get_max();
	  for(i=candidates.size; IS_OK(wiped) && --i>=0; ) {

	    var = candidates[i];

	    if(scope[var].get_max() > val) {
#ifdef _DEBUG_MAX
	      if(_DEBUG_MAX)
		{
		  std::cout << " => " << scope[var] << " in " << scope[var].get_domain() << " <= " << val << std::endl;
		}
#endif
	      event_type[var] |= scope[var].set_max(val);
	      //if(event_type[var] == FAIL_EVENT) wiped = FAILURE(var);
	      if(FAILED(event_type[var])) wiped = FAILURE(var);
	      else if(!changes.contain(var)) {
#ifdef _DEBUG_MAX
		if(_DEBUG_MAX)
		  {
		    std::cout << "    " << scope[var] << " has changed: " << scope[var].get_domain() << std::endl;
		  }
#endif
		changes.add(var);
		if((scope[var].get_max() < scope[n].get_min()) && candidates.contain(var)) {
		  candidates.reversible_remove(var);
#ifdef _DEBUG_MAX
		  if(_DEBUG_MAX)
		    {
		      std::cout << "    " << scope[var] << " can no longer be max" << std::endl;
		    }
#endif
		}
	      }
	    }
	  }
	}
	
	// if the lower bound has been changed, some active variables may become unactive
	// The max(maximum) might also change because variables became unactive
	if(IS_OK(wiped) && LB_CHANGED(event_type[n])) {
	  aux = val = scope[n].get_min();
	  
#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
	  std::cout << " relax the constraint from vars lower than " << val << std::endl;
  }
#endif
	  
	  for(i=candidates.size; IS_OK(wiped) && i;) {
	    --i;
	    var = candidates[i];
	    
	    //std::cout << i << ": " << 
	    
	    if(scope[var].get_max() < val) {
	      
#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
	      std::cout << "    relax1 " << this << " from " << scope[var] << std::endl;
  }
#endif 
	      
	      //relax_from(var);
	      if(candidates.contain(var)) candidates.reversible_remove(var);
	      //std::cout << candidates << std::endl;
	      
	    } else if(aux < scope[var].get_min()) aux = scope[var].get_min();
	  }
	  
#ifdef _DEBUG_MAX
	  if(_DEBUG_MAX)
	    {
	      std::cout << " => " << scope[n] << " in " << scope[n].get_domain() << " <= " << aux << std::endl;
	    }
#endif
	  
	  if(scope[n].set_min(aux) == FAIL_EVENT) wiped = FAILURE(n);	
	}
      } 
      
      if(IS_OK(wiped)) {
	
	// store the max max of the changed vars in aux, and the max max in val
	aux = -INFTY;
	val = -INFTY;
	while(!changes.empty()) {
	  evt = changes.pop();
	  
#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
	  std::cout << "-event on " << scope[evt] 
		    << (LB_CHANGED(event_type[evt]) ? " (lb) " : " ")
		    << (UB_CHANGED(event_type[evt]) ? "(ub) " : " ")
		    << std::endl;
  }
#endif
	  
	  if(UB_CHANGED(event_type[evt]) 
	     && scope[evt].get_max()>aux
	     ) {
	    aux = scope[evt].get_max();
	    if(aux < scope[n].get_min()) {
	      
#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
	      std::cout << "    relax2 " << this << " from " << scope[evt] << std::endl;
  }
#endif 

	      //relax_from(evt);
	      if(candidates.contain(evt)) candidates.reversible_remove(evt);
	      //std::cout << candidates << std::endl;
	    }
	  }
	  if(LB_CHANGED(event_type[evt]) && scope[evt].get_min()>val)
	    val = scope[evt].get_min();
	  //if(ASSIGNED(event_type[evt]))
	  //candidates.reversible_remove(evt);
	}
	
#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
	if(val > -INFTY) 
	  std::cout << " New min maximum: " << val << std::endl;
	if(aux > -INFTY) 
	  std::cout << " New max maximum: " << aux << std::endl;
  }
#endif
      
      
	if(val > -INFTY) {
	
#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
	  std::cout << " => " << scope[n] << " in " << scope[n].get_domain() 
		    << " >= " << val << std::endl;
  }
#endif
	
	  event_type[n] = scope[n].set_min(val);
	  if(event_type[n] == FAIL_EVENT) wiped = FAILURE(n);
	  else {
	    changes.add(n);
	  }
	}
	
	if(aux > -INFTY) {
	  val = scope[n].get_max();
	  
#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
	  std::cout << " previous max minimum witness: " << scope[last_max] 
		    << " in " << scope[last_max].get_domain() << std::endl;
  }
#endif
      
	  if(aux < val && scope[last_max].get_max() < val) {
	  
	    // std::cout << candidates << std::endl;
	    // for(i=0; i<scope.size; ++i) {
	    //   std::cout << " " << scope[i].get_domain();
	    // }
	    // std::cout << std::endl;
	    
	    // look for a new witness for the max
	    i = candidates.size;
	    while( i > 0 ) {
	      --i;
	      
	      // std::cout << i << std::endl;
	      
	      // std::cout << candidates[i] << std::endl;
	      
	      var = scope[candidates[i]].get_max();
	      if(aux < var) aux = var;
	      if(var >= val) {
		last_max = candidates[i];
		break;
	      } 
	    }
	  
	    if(aux < val) {
#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
	      std::cout << " => " << scope[n] << " in " << scope[n].get_domain() 
			<< " <= " << aux << std::endl;
  }
#endif
	      event_type[n] = scope[n].set_max(aux);
	      if(event_type[n] == FAIL_EVENT) wiped = FAILURE(n);
	      else {
		if(!changes.contain(n)) changes.add(n);
	      }
	    }
	  } else {
	    
#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
	    std::cout << " no need to update " << scope[n] << "'s max" << std::endl;
  }
#endif
	    
	  }
	}
      }
    }  
  }

#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
  std::cout << candidates << std::endl;
  }
#endif


  if(candidates.size == 1) {
    
#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
  {
    std::cout << scope[candidates[0]] << " is the last max var " << std::endl;
  }
#endif
    
    if(scope[candidates[0]].set_domain(scope[n]) == FAIL_EVENT) wiped = FAILURE(candidates[0]);
    if(scope[n].set_domain(scope[candidates[0]]) == FAIL_EVENT) wiped = FAILURE(n);
  } else if(candidates.empty()) {
    wiped = n;
  }

#ifdef _DEBUG_MAX
  if(_DEBUG_MAX)
    {
      std::cout << wiped << " max( ";
      for(int k=0; k<scope.size-1; ++k)
	std::cout // << scope[k] << " in " 
	  << scope[k].get_domain() << " ";
      std::cout << ") = " // << scope[scope.size-1] << " in "
		<< scope[scope.size-1].get_domain() << std::endl << std::endl;
    }
#endif
  
  return wiped;
}

std::ostream& Mistral::PredicateMax::display(std::ostream& os) const {
  os << scope.back()/*.get_var()*/ << " == max(" << scope[0]/*.get_var()*/;
  for(unsigned int i=1; i<scope.size-1; ++i) {
    os << ", " << scope[i]/*.get_var()*/ ;
  }
  os << ")";
  return os;
}


// // Mistral::ConstraintNogoodBase::ConstraintNogoodBase(Vector< Variable >& scp) 
// //   : Constraint(scp) { }

// // void Mistral::ConstraintNogoodBase::mark_domain() {
// //   for(unsigned int i=0; i<scope.size; ++i) {
// //     solver->mark_non_convex(scope[i].id());
// //   }
// // }

// // void Mistral::ConstraintNogoodBase::initialise() {
// //   Constraint::initialise();
// //   //int min_idx = INFTY;
// //   //int max_idx = -INFTY;
// //   for(unsigned int i=0; i<scope.size; ++i) {
// //     //if(scope[i].get_id() > max_idx) max_idx = scope[i].get_id();
// //     //if(scope[i].get_id() < min_idx) min_idx = scope[i].get_id();
// //     trigger_on(_VALUE_, i);
// //   }
// //   set_idempotent(true);

// //   //if(// min_idx < 
// //   //max_idx) {
// //   watch_structure.initialise(0,scope.size);
// //     //}

// //   for(unsigned int i=0; i<scope.size; ++i) {
// //     //add(scope[i]);
// //     watch_structure[i] = new Vector< Array< Literal >* > [scope[i].get_size()];
// //     watch_structure[i]-=scope[i].get_min();
// //   }
// // }

// // Mistral::ConstraintNogoodBase::~ConstraintNogoodBase() {
// // }

// // void Mistral::ConstraintNogoodBase::add(Variable x) {
// //   unsigned int idx = x.id();
// //   if(idx >= scope.size) {
// //     while(scope.capacity <= idx)
// //       scope.extendStack();
// //     scope[idx] = x;
// //     while(watch_structure.capacity <= idx)
// //       watch_structure.extendStack();
// //   }
// // }

// // void Mistral::ConstraintNogoodBase::add( Vector < Literal >& clause ) {
// //  if(clause.size > 1) {
// //    Array<Literal> *cl = Array<Literal>::Array_new(clause);
// //    nogood.add( cl );
// //    // std::cout << clause[0].get_atom() << " " << (1-clause[0].get_value()) << std::endl;
// //    // std::cout << watch_structure[clause[0].get_atom()][1-clause[0].get_value()] << std::endl;

// //    watch_structure[clause[0].get_atom()][1-clause[0].get_value()].add(cl);
// //    watch_structure[clause[1].get_atom()][1-clause[1].get_value()].add(cl);
// //  } else {
// //    clause[0].apply(scope[clause[0].get_atom()]);
// //  }
// // }

// // int Mistral::ConstraintNogoodBase::check( const int* sol ) const {
// //   unsigned int i, j;
// //   bool satisfied=true;

// //   for(i=0; i<nogood.size; ++i) {
// //     satisfied = true;
// //     Array < Literal >& clause = *(nogood[i]);
// //     for(j=0; satisfied && j<clause.size; ++j) {
// //       satisfied = clause[j].check(sol[j]);
// //     }
// //   }
  
// //   return !satisfied;
// // }

// // Mistral::PropagationOutcome Mistral::ConstraintNogoodBase::propagate() {
// //   Array< Literal > *conflict=NULL;

// //   int x, v, cw;
// //   Literal p;
// //   //unsigned int i;
// //   while( !changes.empty() ) {
// //     std::cout << changes << std::endl;
// //     //std::cout << scope << std::endl;
// //     x = changes.pop();
// //     //std::cout << scope[x] << " in " << scope[x].get_domain() << std::endl;
// //     v = scope[x].get_min();
// //     std::cout << x << "=" << v << ": " << watch_structure[x][v] << std::endl;

// //     p = Literal(x,EQ,v);


// //     cw = watch_structure[x][v].size;
// //     while(cw-- && !conflict) {
// //       conflict = update_watcher(cw, x, v);
// //     }
// //   }

// //   return CONSISTENT;
// // }


// // inline Array< Literal >* SatSolver::update_watcher(const int cw, const int x, const int v)
// // {
// //   Array< Literal > *cl = watch_structure[x][v][cw];
// //   Array< Literal >& clause = *cl;


// //   unsigned int j;
// //   Lit q, r;
// //   Atom v, w;

// // #ifdef _DEBUG_WATCH
// //   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// //   std::cout << " ";
// //   print_clause( std::cout, watch_structure[p][cw] );
// //   std::cout << std::endl;
// // #endif

// //   //ensure that p is the second watched lit
// //   if( clause[1] != p ) {
// //     q = clause[1];
// //     clause[0] = q;
// //     clause[1] = p;
// //   } else q = clause[0];
// //   v=state[UNSIGNED(q)];

// // #ifdef _DEBUG_WATCH
// //   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// //    std::cout << " second watched: ";
// //    printLiteral(std::cout, q);
// //    std::cout << std::endl;

// //    if( LEVEL(v) < assumptions.size && SIGN(v) == SIGN(q) ) {
// //      for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// //      std::cout << " satisfied by second watched" << std::endl;
// //    }
// // #endif

// //   //check if the other watched lit is assigned
// //   if( LEVEL(v) >= assumptions.size || SIGN(v) != SIGN(q) ) {
// //     for(j=2; j<clause.size; ++j) {
// //       // for each literal q of the clause,
// //       r = clause[j];
// //       w = state[UNSIGNED(r)];

// // #ifdef _DEBUG_WATCH
// //   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// //        std::cout << "\tcheck another lit: ";
// //        printLiteral(std::cout, r);
// //        std::cout << std::endl;
// // #endif

// //       if( LEVEL(w) >= assumptions.size ) { // this literal is not set
// // 	// then it is a good candidate to replace p
// // 	clause[1] = r;
// // 	clause[j] = p;
// // 	watch_structure[p].remove(cw);
// // 	watch_structure[r].add(cl);

// // #ifdef _DEBUG_WATCH
// //   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// // 	std::cout << " new watched" << std::endl;
// // #endif

// // 	break;	
// //       }
// //       // if it is set true, then the clause is satisfied
// //       else if( SIGN(w) == SIGN(r) ) {

// // #ifdef _DEBUG_WATCH
// //   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// // 	std::cout << " satisfied" << std::endl;
// // #endif

// // 	break;
// //       }
// //     }
      
// //     if( j == clause.size ) // no replacement could be found
// //       { 

// // #ifdef _DEBUG_WATCH
// //   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// // 	std::cout << " no more watched" << std::endl;
// // #endif

// // 	if( LEVEL(v) >= assumptions.size ) {
// // 	  // the last literal (other watched lit) is not set yet, we set it
// // 	  add_lit(q);
// // 	  reason[UNSIGNED(q)] = cl;

// // 	  //#ifdef _DEBUGNOGOOD
// // #ifdef _DEBUG_WATCH
// //   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// // 	  std::cout //<< "unit prune disjunct b" << y 
// // 	    << " because ";
// // 	  print_clause( std::cout, cl );
// // 	  std::cout << std::endl;
// // #endif

// // 	  //	  std::cout << " prop " << q << std::endl;
	  
// // 	} else 
// // 	  // it is set to false already, we fail
// // 	  if( // polarity[y] == -q
// // 	     SIGN(v) != SIGN(q)
// // 	      ) {

// // #ifdef _DEBUG_WATCH
// //   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// // 	    std::cout << " fail" << std::endl;
// // #endif

// // 	    return cl;
// // 	  }
// //       }
// //   }

// //   return NULL;
// // }




// // std::ostream& Mistral::ConstraintNogoodBase::display(std::ostream& os) const {
// //   os << " (";
// //   if(nogood.size>0) {
// //     os << nogood[0];
// //     for(unsigned int i=1; i<nogood.size; ++i)
// //       os << " " << nogood[i]  ;
// //   }
// //   os << ")";
// //   return os;
// // }



// // Mistral::ConstraintClauseBase::ConstraintClauseBase(Vector< Variable >& scp) 
// //   : Constraint(scp) { }

// // void Mistral::ConstraintClauseBase::mark_domain() {
// //   for(unsigned int i=0; i<scope.size; ++i) {
// //     solver->mark_non_convex(scope[i].id());
// //   }
// // }

// // void Mistral::ConstraintClauseBase::initialise() {
// //   Constraint::initialise();
// //   for(unsigned int i=0; i<scope.size; ++i) {
// //     trigger_on(_VALUE_, i);
// //   }
// //   set_idempotent(true);

// //   is_watched_by.initialise(0,2*scope.size);


// //   // for(unsigned int i=0; i<scope.size; ++i) {
// //   //   is_watched_by[i] = new Vector< Clause* >;
// //   //   is_watched_by[i]-=scope[i].get_min();
// //   // }
// // }

// // Mistral::ConstraintClauseBase::~ConstraintClauseBase() {
// // }

// // void Mistral::ConstraintClauseBase::add(Variable x) {
// //   unsigned int idx = x.id();
// //   if(idx == scope.size) {
// //     scope.add(x);
// //     while(is_watched_by.capacity <= 2*idx)
// //       is_watched_by.extendStack();
// //   } else if(idx > scope.size) {
// //     while(scope.capacity <= idx)
// //       scope.extendStack();
// //     scope[idx] = x;
// //     while(is_watched_by.capacity <= idx)
// //       is_watched_by.extendStack();
// //   }
// // }

// // void Mistral::ConstraintClauseBase::add( Vector < Literal >& clause ) {
// //  if(clause.size > 1) {
// //    Clause *cl = Clause::Array_new(clause);
// //    clause.add( cl );
// //    is_watched_by[clause[0]].add(cl);
// //    is_watched_by[clause[1]].add(cl);
// //  } else {
// //    scope[UNSIGNED(clause[0])].set_domain(SIGN(clause[0]));
// //  }
// // }

// // int Mistral::ConstraintClauseBase::check( const int* sol ) const {
// //   unsigned int i, j;
// //   bool satisfied=true;

// //   for(i=0; i<clause.size; ++i) {
// //     satisfied = true;
// //     Clause& clause = *(clause[i]);
// //     for(j=0; satisfied && j<clause.size; ++j) {
// //       satisfied = //clause[j].check(sol[j]);
// // 	(sol[j] == SIGN(clause[j]));
// //     }
// //   }
  
// //   return !satisfied;
// // }

// // Mistral::PropagationOutcome Mistral::ConstraintClauseBase::propagate() {
// //   Clause *conflict=NULL;


// //   int x, v, cw;
// //   Lit p;
// //   //unsigned int i;
// //   while( !changes.empty() ) {
// //     std::cout << changes << std::endl;
// //     //std::cout << scope << std::endl;
// //     x = changes.pop();
// //     //std::cout << scope[x] << " in " << scope[x].get_domain() << std::endl;
// //     v = scope[x].get_min();
// //     std::cout << x << "=" << v << ": " << is_watched_by[x][v] << std::endl;

// //     p = 2*x+v;

// //     cw = is_watched_by[p].size;
// //     while(cw-- && !conflict) {
// //       conflict = update_watcher(cw, p);
// //     }
// //   }

// //   return CONSISTENT;
// // }


// // inline Clause* SatSolver::update_watcher(const int cw, const Lit p)
// // {
// //   Clause *cl = is_watched_by[p][cw];
// //   Clause& clause = *cl;

// //   unsigned int j;
// //   Lit q, r;
// //   //Atom v, w;
// //   Variable v, w;

// //   //ensure that p is the second watched lit
// //   if( clause[1] != p ) {
// //     q = clause[1];
// //     clause[0] = q;
// //     clause[1] = p;
// //   } else q = clause[0];
// //   v=scope[UNSIGNED(q)].;


// //   //check if the other watched lit is assigned
// //   //if( LEVEL(v) >= assumptions.size || SIGN(v) != SIGN(q) ) {
// //   if( !v.is_ground() || v.get_min() != SIGN(q) ) {
// //     for(j=2; j<clause.size; ++j) {
// //       // for each literal q of the clause,
// //       r = clause[j];
// //       w = scope[UNSIGNED(r)];

// //       if( !w.is_ground() ) { // this literal is not set
// // 	// then it is a good candidate to replace p
// // 	clause[1] = r;
// // 	clause[j] = p;
// // 	is_watched_by[p].remove(cw);
// // 	is_watched_by[r].add(cl);

// // 	break;	
// //       }
// //       // if it is set true, then the clause is satisfied
// //       else if( w.get_min() == SIGN(r) ) {
// // 	break;
// //       }
// //     }
      
// //     if( j == clause.size ) // no replacement could be found
// //       { 
// // 	if( !v,is_ground() ) {
// // 	  // the last literal (other watched lit) is not set yet, we set it
// // 	  add_lit(q);
// // 	  //reason[UNSIGNED(q)] = cl;
// // 	} else 
// // 	  // it is set to false already, we fail
// // 	  if( v.get_min() != SIGN(q) ) {
// // 	    return cl;
// // 	  }
// //       }
// //   }

// //   return NULL;
// // }




// // std::ostream& Mistral::ConstraintClauseBase::display(std::ostream& os) const {
// //   os << " (";
// //   if(clauses.size>0) {
// //     os << clauses[0];
// //     for(unsigned int i=1; i<clauses.size; ++i)
// //       os << " " << clauses[i]  ;
// //   }
// //   os << ")";
// //   return os;
// // }




// std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Literal& x) {
//   return x.display(os);
// }

// std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Literal* x) {
//   return x->display(os);
// }


// std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::ConstraintTriggerArray& x) {
//   return x.display(os);
// }

// std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::ConstraintTriggerArray* x) {
//   return x->display(os);
// }
