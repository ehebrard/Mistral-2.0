
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


#include <mistral_variable.hpp>
#include <mistral_constraint.hpp>


Mistral::Constraint::Constraint() {
  priority = 2;
  id = -1;
  //weight = 1.0;
}

Mistral::Constraint::Constraint(Vector< Variable >& scp) {
  priority = 2;
  id = -1;
  //weight = 1.0;
  //initialise(scp);
  for(unsigned int i=0; i<scp.size; ++i) {
    //scope.add(scp[i].get_var());
    scope.add(scp[i]);
  }
}

Mistral::Constraint::~Constraint() {
  if(changes.list_ == events.list_)
    events.list_ = NULL;

  delete [] solution;
  delete [] self;
  delete [] trigger;
  delete [] event_type;
  //delete [] _scope;
  //scope = NULL;
  
}

// void Mistral::Constraint::initialise(Vector< Variable >& scp) {
//   for(unsigned int i=0; i<scp.size; ++i) {
//     scope.add(scp[i]);
//     //trigger.add(_domain_);
//   }
//   initialise();
// }

void Mistral::Constraint::initialise() {
  //arity = scp.size;
  //weight = w;

  //std::cout << "Initialise " << this << std::endl;

  trail_.initialise(0,2*scope.size);
  trail_.add(-1);

  //_scope = new IntVar[scope.size];
  event_type = new Event[scope.size];
  self = new int[scope.size];
  trigger = new int[scope.size];
  solution = new int[scope.size];
  changes.initialise(0, scope.size-1, false);
  supports = NULL;


  //std::memcpy(_scope, scp.stack_, scope.size*sizeof(IntVar));
  std::fill(event_type, event_type+scope.size, NO_EVENT);
  std::fill(self, self+scope.size, -1);
  std::fill(trigger, trigger+scope.size, _domain_);
  std::fill(solution, solution+scope.size, 0);

  active.initialise(0, scope.size-1, true);
  for(unsigned int i=0; i<scope.size; ++i) 
    if(scope[i].domain_type != BOOL_VAR && scope[i].is_ground()) active.erase(i);

  //scope = _scope;
}


void Mistral::Constraint::post(Solver *s) {
  solver = s;
  unsigned int i, j;
  // for each of its variables
  for(i=0; i<scope.size; ++i) {
    // add the constraint to the list of constraints on that variable
    j = scope[i].id();
    ConstraintTrigger ct(this, i);
    int trg = trigger[i];
    ConstraintList *lst = solver->constraint_graph[j];
    unsigned int elt = lst->reversible_add(ct, trg);
    self[i] = elt;
  }
}

void Mistral::Constraint::relax() {
  unsigned int i, j;
  int k;

  //#ifdef _DEBUG_AC
  std::cout << "relax " << this << " from ";
  //#endif


  for(i=0; i<active.size; ++i) {
    j = active[i];
    k = scope[j].id();

    //#ifdef _DEBUG_AC
    std::cout << scope[j] << "'s c-list " ;
    //#endif

    solver->constraint_graph[k]->reversible_erase(self[j], trigger[j]);
  }

  //#ifdef _DEBUG_AC
  std::cout << std::endl;
  //#endif

  std::cout << solver << std::endl;
  
}

Mistral::Constraint* Mistral::Constraint::notify_assignment(const int var, const int level) {
  //display(std::cout);
  //std::cout << std::endl;
  std::cout << "remove " << scope[var] << " from " ;
  display(std::cout);
  std::cout << active << " => " ;
  
  Constraint *r = NULL;
  if(trail_.back() != level) {
	trail_.add(active.size);
	trail_.add(level);
	r = this;
      }
  active.erase(var);
  
  std::cout << active << std::endl;
  
  return r;
}
// Mistral::Constraint* Mistral::Constraint::notify_assignment(const int var, const int level) {

//  //   std::cout << "notify var: rem " << var 
// //  	    << " from " << active << std::endl;


// #ifdef _DEBUG_AC
//    std::cout << this << " -= " << scope[var] << std::endl;
// #endif
//    Constraint *r = NULL;

//   assert(active.contain(var));

//   if(trail_.back() != level) {
//     trail_.add(active.size);
//     trail_.add(level);
//     r = this;
//   }
//   active.erase(var);

// #ifdef _DEBUG_AC
//   std::cout << "number of active variables left: " << active.size << std::endl;
// #endif

//       // BUGGY: the constraint might be relaxed from the 
//       //        last variable's list without getting a trigger
//       //        from it. The 'changes' list might thus not be complete
//       //        instead we relax it when freezing.

//   //if(active.size == 1) relax();
//   return r;
// }

// // void Mistral::Constraint::entail() {
// // //   unsigned int i, j;
// // //   for(i=0; i<active.size; ++i) {
// // //     j = active[i];
// // //     scope[j]->constraints.reversible_erase(self[j]);
// // //   }
// // }

// void Mistral::Constraint::restore() {
//   trail_.pop();
//   active.size = trail_.pop();
// }

// void Mistral::Constraint::trigger_on(const int t, const int x) {
//   trigger[x] = t;
// //   //self_list[x] = &(_scope[x]->constraints);
// //  ConstraintTrigger ct(this, x);
// //   self[x] = _scope[x]->constraints.create(ct, t);
// }

//void Mistral::Constraint::post() {
//   for(int i=0; i<scope.size; ++i)
//     _scope[i]->constraints.reversible_add(self[i], trigger[i]);
//}
 
//void Mistral::Constraint::relax() {
//   for(int i=0; i<scope.size; ++i)
//     _scope[i]->constraints.reversible_erase(self[i], trigger[i]);
//}

// void Mistral::Constraint::assign(const int var) {
// //   active.erase(var);
// //   if(active.size == 1) {
// //     int i = active.back();
// //     _scope[i]->constraints.erase(self[i], trigger[i]);
// //   }
// }

// void Mistral::Constraint::unassign(const int var) {
// //   if(active.size == 1) {
// //     int i = active.back();
// //     _scope[i]->constraints.add(self[i], trigger[i]);
// //   } 
// //   if(!(active.contain(var))) active.add(var);
// }

bool Mistral::Constraint::firstSupport(const int vri, const int vli) 
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

bool Mistral::Constraint::findSupport(const int vri, const int vli) 
{
  int i=scope.size, vali;
  bool found=false;
  // sol is initialized: either to the value 
  // a variable is already assigned to
  // or to the first value in its domain
  while(i >= 0) {
    // check this assignment
    if( !check( solution ) ) {
      found=true;
      if( supports ) {
	vali = scope.size;
	while( vali-- )
	  supports[vri][vli][vali] = solution[vali];
      }
      break;
    }
 
    // try to assign more things
    // find the last var whose domain we have not exhausted
    --i;
    while( i >= 0 ) {
      if( i == vri || scope[i].is_ground() ) {
	--i;
	continue;
      }
      solution[i] = scope[i].next( solution[i] );
      if(solution[i] != NOVAL)
      	break;
      else
	solution[i] = scope[i].get_min();
      --i;
    }
    if( i >= 0 )
      i = scope.size;
  } 
  return found;
}


// bool Mistral::Constraint::propagate()
// {
//   int i, consistent=1, evt = ( Constraint::RANGETRIGGER );
//   //int i, consistent=1, evt = ( Constraint::VALUETRIGGER );
//   for( i=0; consistent && i<scope.size; ++i )
//     consistent = propagate( i, evt );
//   return consistent;
// } 

// std::string Mistral::Constraint::getString() const {
//   std::string return_str = name()+"("+(scope[0]->getString());
//   for(int i=1; i<scope.size; ++i)
//     return_str += (", "+(scope[1]->getString()));
//   return return_str+")";
// }

std::ostream& Mistral::Constraint::display(std::ostream& os) const {
  os << name() << "(" << scope[0];
  for(unsigned int i=1; i<scope.size; ++i)
    os << ", " << scope[1];
  os << ")";
  return os;
}

// void Mistral::ConstraintNotEqual::post() {
//   for(unsigned int i=0; i<scope.size; ++i)
//     trigger_on(_value_, i);
// }

// Mistral::InlineConstraint* Mistral::InlineConstraint::Constraint_new(Vector< Variable >& X)
// {
//   void* mem = malloc(sizeof(Constraint) + 
// 		     sizeof(Variable)*X.size // + 
// // 		     sizeof(int)*5*(X.size)
// 		     );
//   return new (mem) InlineConstraint(X); 
// }



// void Mistral::InlineConstraint::initialise() {

//   Constraint::initialise();

// //   trail_.initialise(0,2*scope.size);
// //   trail_.add(-1);

// //   //_scope = new IntVar[scope.size];
// //   solution = new int[scope.size];
// //   self = new int[scope.size];
// //   trigger = new int[scope.size];
// //   //self = &data[scope.size*(sizeof(Variable)/sizeof(int))+scope.size];
// //   //trigger = &data[scope.size*(sizeof(Variable)/sizeof(int))+2*scope.size];
  
// //   event_type = &data[scope.size*(sizeof(Variable)/sizeof(int))];
// //   changes.list_ = &data[scope.size*(sizeof(Variable)/sizeof(int))+scope.size];
// //   changes.index_ = &((unsigned int *)data)[scope.size*(sizeof(Variable)/sizeof(int))+2*scope.size];
// //   changes.start_ = changes.index_;
// //   changes.size = 0;
// //   changes.capacity = scope.size;

// //   active.list_ = &data[scope.size*(sizeof(Variable)/sizeof(int))+3*scope.size];
// //   active.index_ = &((unsigned int *)data)[scope.size*(sizeof(Variable)/sizeof(int))+4*scope.size];
// //   active.start_ = active.index_;
// //   active.size = 0;
// //   active.capacity = scope.size;

// //   //std::memcpy(_scope, scp.stack_, scope.size*sizeof(IntVar));
// //   std::fill(event_type, event_type+scope.size, NO_EVENT);
// //   std::fill(self, self+scope.size, -1);
// //   std::fill(trigger, trigger+scope.size, _domain_);
// //   std::fill(solution, solution+scope.size, 0);

// //   for(unsigned int i=0; i<scope.size; ++i) {
// //     active.index_[i] = i;
// //     active.list_[i] = i;
// //     changes.index_[i] = i;
// //     changes.list_[i] = i;
// //   }

// //   for(unsigned int i=0; i<scope.size; ++i) {
// //     //active.index_[i] = i;
// //     if(scope[i].is_ground()) active.erase(i);
// //   }
// }
// Mistral::InlineConstraint::InlineConstraint(Vector< Variable >& scp) 
//   : Constraint() {

//   scope.stack_ = &((Variable*)data)[0];
//   //scope.size = scp.size;
//   scope.size = 0;
//   scope.capacity = scp.size;
//   //memcpy(scope.stack_, scp.stack_, sizeof(Variable));
//   for(unsigned int i=0; i<scp.size; ++i) {
//     scope.add(scp[i]);
//   }
// }
// Mistral::InlineConstraint::~InlineConstraint() {
// //   event_type = NULL;
// //   changes.list_ = NULL;
// //   active.list_ = NULL;
// //   changes.start_ = NULL;
// //   active.start_ = NULL;
// }


Mistral::PropagationOutcome Mistral::ConstraintNotEqual::propagate() {
  int var = 1-changes[0];
  if(active.size)
    return(scope[var].remove(scope[1-var].get_min()) == FAIL_EVENT ? FAILURE(var) : CONSISTENT);
  else
    return(scope[0].get_min() == scope[1].get_min() ? FAILURE(var) : CONSISTENT);
}

// std::string Mistral::ConstraintNotEqual::getString() const {
//   return (toString(scope[0])+" =/= "+toString(scope[1]));
// }

std::ostream& Mistral::ConstraintNotEqual::display(std::ostream& os) const {
  os << scope[0] << " =/= " << scope[1];
  return os;
}


Mistral::PropagationOutcome Mistral::ConstraintEqual::propagate() {

  unsigned int i;
  int var;

//   std::cout << "propagate " << (this) << " b/c" ;
//   for(i=0; i<this->changes.size; ++i) 
//     std::cout << " " << this->scope[this->changes[i]];
//   std::cout << std::endl;
//   for(i=0; i<this->scope.size; ++i) {
//     std::cout << this->scope[i] << ": " << (this->scope[i].get_domain()) << " ";
//   }
//   std::cout << std::endl;
  

  PropagationOutcome wiped = CONSISTENT;
  if(active.size) {
    for(i=0; IS_OK(wiped) && i<changes.size; ++i) {

      var = changes[i];

      //std::cout << scope[var] << " has changed" <<std::endl;

      if(scope[1-var].set_domain(scope[var]) == FAIL_EVENT)
	wiped = FAILURE(1-var);
    }
  } else if(scope[0].get_min() != scope[1].get_min()) {

    //std::cout << "all ground" <<std::endl;

    wiped = FAILURE(0);
  }


//   if(IS_OK(wiped)) {
//     for(unsigned int i=0; i<this->scope.size; ++i) {
//       std::cout << this->scope[i] << ": " << (this->scope[i].get_domain()) << " ";
//     }
//   } else {
//     std::cout << "inconsistent!!" << std::endl;
//   }  
//   std::cout << std::endl;


  return wiped;
}


std::ostream& Mistral::ConstraintEqual::display(std::ostream& os) const {
  os << scope[0] << " == " << scope[1];
  return os;
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

std::ostream& Mistral::PredicateUpperBound::display(std::ostream& os) const {
  os << scope[1] << " <=> (" << scope[0] << " <= " << bound << ")";
  return os;
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

std::ostream& Mistral::PredicateLowerBound::display(std::ostream& os) const {
  os << scope[1] << " <=> (" << scope[0] << " >= " << bound << ")";
  return os;
}


Mistral::PropagationOutcome Mistral::PredicateLess::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[2].is_ground() ) {
    if( scope[2].get_min() ) { // x[0] + k <= x[1]
      if( IS_FAIL(scope[0].set_max(scope[1].get_max()-offset)) )
	wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].set_min(scope[0].get_min()+offset)) )
	wiped = FAILURE(1);
    } else if( scope[2].get_max() == 0 ) { // x[1] + (1-k) <= x[0]
      if( IS_FAIL(scope[0].set_min(scope[1].get_min()-offset+1)) )
	wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].set_max(scope[0].get_max()-offset+1)) )
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

std::ostream& Mistral::PredicateLess::display(std::ostream& os) const {
  os << scope[2] << " <=> (" << scope[0] ;
  if(offset==0) os << " <= ";
  else if(offset==1) os << " < ";
  else os << " + " << offset << " <= ";
  os << scope[1] << ")";
  return os;
}
// void Mistral::ConstraintLess::post() {
//   for(unsigned int i=0; i<scope.size; ++i)
//     trigger_on(_range_, i);
// }

Mistral::PropagationOutcome Mistral::ConstraintLess::propagate() {
  if(changes.contain(0) && LB_CHANGED(event_type[0])) {
    if(scope[1].set_min(scope[0].get_min() + offset) == FAIL_EVENT) 
      return FAILURE(1);
  }
  if(changes.contain(1) && UB_CHANGED(event_type[1])) {
    if(scope[0].set_max(scope[1].get_max() - offset) == FAIL_EVENT) 
      return FAILURE(0);
  }
  return CONSISTENT;
}

// std::string Mistral::ConstraintLess::getString() const {
//   std::string return_str = toString(scope[0]);
//   if(offset < 0) return_str += (" - "+toString(-offset+1)+" < ");
//   else if(offset > 1) return_str += (" + "+toString(offset-1)+" < ");
//   else if(offset > 0) return_str += (" < ");
//   else return_str += (" <= ");
//   return_str += toString(scope[1]);
//   return return_str;
// }

std::ostream& Mistral::ConstraintLess::display(std::ostream& os) const {
  os << scope[0];
  if(offset < 0) os << " - " << (-offset+1) << " < ";
  else if(offset > 1) os << " + " << (offset-1) << " < ";
  else if(offset > 0) os << " < ";
  else os << " <= ";
  os << scope[1];
  return os;
}

// Mistral::PredicateEqual::post() {
//   for(unsigned int i=0; i<scope.size; ++i)
//     trigger_on(_domain_, i);
// }

Mistral::PropagationOutcome Mistral::PredicateEqual::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[2].is_ground() ) {
    if( (spin + scope[2].get_min()) != 1 ) {
      if(scope[0].is_ground() && !(scope[1].remove(scope[0].get_min())))
	wiped = FAILURE(1);
      else {
	if(scope[1].is_ground() && !(scope[0].remove(scope[1].get_min())))
	  wiped = FAILURE(0);
      }
    } else {
      if( !(scope[0].set_domain(scope[1])) ) wiped = FAILURE(0);
      else if( !(scope[1].set_domain(scope[0])) ) wiped = FAILURE(1);
    }
  } else {
    if( !(scope[0].intersect(scope[1])) ) {
      if(!(scope[2].remove(spin))) wiped = FAILURE(2);	    
    } else { 
      if( scope[0].is_ground() && scope[1].is_ground() ) {
	if(!(scope[2].set_domain(spin))) wiped = FAILURE(2);
      }
    }
  }
  
  return wiped;
}

std::ostream& Mistral::PredicateEqual::display(std::ostream& os) const {
  os << scope[2] << " <=> (";
  if(spin) os << scope[0] << " == " << scope[1];
  else os << scope[0] << " =/= " << scope[1];
  os << ")";
  return os;
}


Mistral::PropagationOutcome Mistral::PredicateConstantEqual::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[1].is_ground() ) {
    if( (spin + scope[1].get_min()) != 1 ) {
      if( !(scope[0].remove(value)) ) wiped = FAILURE(0);
    } else {
      if( !(scope[0].set_domain(value)) ) wiped = FAILURE(0);
    }
  } else {
    if( !(scope[0].contain(value)) ) {
      if(!(scope[1].remove(spin))) wiped = FAILURE(1);	    
    } else if( scope[0].is_ground() ) {
      if(!(scope[1].set_domain(spin))) wiped = FAILURE(1);
    }
  }
  
  return wiped;
}

std::ostream& Mistral::PredicateConstantEqual::display(std::ostream& os) const {
  os << scope[1] << " <=> (";
  if(spin) os << scope[0] << " == " << value;
  else os << scope[0] << " =/= " << value;
  os << ")";
  return os;
}


Mistral::PropagationOutcome Mistral::PredicateOffset::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( IS_FAIL(scope[0].set_min( scope[1].get_min() - offset )) ) wiped = FAILURE(0); 
  else if( IS_FAIL(scope[0].set_max( scope[1].get_max() - offset )) ) wiped = FAILURE(0); 
  else if( IS_FAIL(scope[1].set_min( scope[0].get_min() + offset )) ) wiped = FAILURE(1); 
  else if( IS_FAIL(scope[1].set_max( scope[0].get_max() + offset )) ) wiped = FAILURE(1);
  
  return wiped;
}

std::ostream& Mistral::PredicateOffset::display(std::ostream& os) const {
  os << scope[1] << " == (" << scope[0] << " + " << offset << ")";
  return os;
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

std::ostream& Mistral::PredicateNot::display(std::ostream& os) const {
  os << scope[1] << " <=> not(" << scope[0] << ")";
  return os;
}


Mistral::PropagationOutcome Mistral::PredicateAnd::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

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
  
  return wiped;
}

std::ostream& Mistral::PredicateAnd::display(std::ostream& os) const {
  os << scope[2] << " <=> (" << scope[0] << " and " << scope[1] << ")";
  return os;
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
  os << scope[2] << " <=> (" << scope[0] << " or " << scope[1] << ")";
  return os;
}


Mistral::PropagationOutcome Mistral::ConstraintAnd::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( IS_FAIL(scope[0].remove(0)) ) 
    wiped = FAILURE(0);
  else if( IS_FAIL(scope[1].remove(0)) ) 
    wiped = FAILURE(1);
  
  return wiped;
}

std::ostream& Mistral::ConstraintAnd::display(std::ostream& os) const {
  os << "(" << scope[0] << " and " << scope[1] << ")";
  return os;
}


// Mistral::ConstraintOr* Mistral::ConstraintOr::ConstraintOr_new(Vector< Variable >& X)
// {
 
//   void* mem = malloc(sizeof(Constraint) + 
// 		     sizeof(Variable)*X.size//  + 
// // 		     sizeof(int)*5*(X.size)
// 		     );
//   return new (mem) ConstraintOr(X); 
// }

// Mistral::ConstraintOr::ConstraintOr(Vector< Variable >& scp) 
//   : Constraint() {

//   scope.stack_ = &((Variable*)data)[0];
//   scope.size = 0;
//   scope.capacity = scp.size;
//   for(unsigned int i=0; i<scp.size; ++i) {
//     scope.add(scp[i]);
//   }
// }

Mistral::PropagationOutcome Mistral::ConstraintOr::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( scope[1].get_max() == 0 ){
    if( IS_FAIL(scope[0].remove(0)) ) wiped = FAILURE(0);
  } else if( scope[0].get_max() == 0 ) {
    if( IS_FAIL(scope[1].remove(0)) ) wiped = FAILURE(1);
  }
  
  return wiped;
}

std::ostream& Mistral::ConstraintOr::display(std::ostream& os) const {
  os << "(" << scope[0] << " or " << scope[1] << ")";
  return os;
}

// void Mistral::PredicateAdd::post() {
//   for(unsigned int i=0; i<scope.size; ++i)
//     trigger_on(_range_, i);
// }

Mistral::PropagationOutcome Mistral::PredicateAdd::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  int last_change = 0;
  Outcome evt;
  
  for(int i=0; i<3; i=(i+1)%3) {
    evt = NO_EVENT;
    if(i!=2) {
      evt |= scope[i].set_min(scope[2].get_min() - scope[1-i].get_max());
      if(IS_FAIL(evt)) { wiped = FAILURE(i); break; }
      evt |= scope[i].set_max(scope[2].get_max() - scope[1-i].get_min()); 
      if(IS_FAIL(evt)) { wiped = FAILURE(i); break; }
    } else {
      evt |= scope[2].set_min(scope[0].get_min() + scope[1].get_min()); 
      if(IS_FAIL(evt)) { wiped = FAILURE(2); break; }
      evt |= scope[2].set_max(scope[0].get_max() + scope[1].get_max()); 
      if(IS_FAIL(evt)) { wiped = FAILURE(2); break; }
    }
    if(evt != NO_EVENT) last_change = i;
    else if(last_change == (i+1)%3) break;
  }
  
  return wiped;
}

std::ostream& Mistral::PredicateAdd::display(std::ostream& os) const {
  os << scope[2] << " == (" << scope[0] << " + " << scope[1] << ")";
  return os;
}


Mistral::ConstraintBoolSumEqual::ConstraintBoolSumEqual(Vector< Variable >& scp, const int t)
  : Constraint(scp) { 
  priority = 2;
  total = t; 
}

void Mistral::ConstraintBoolSumEqual::initialise() {
  Constraint::initialise();
  for(unsigned int i=0; i<scope.size; ++i)
    trigger_on(_value_, i);
  //set_idempotent(false);
  set_idempotent(true);
}

Mistral::ConstraintBoolSumEqual::~ConstraintBoolSumEqual() 
{
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
  os << "(" << scope[0] ;
  for(unsigned int i=1; i<scope.size; ++i) 
    os << " + " << scope[i];
  os << ") == " << total ;
  return os;
}


Mistral::ConstraintCliqueNotEqual::ConstraintCliqueNotEqual(Vector< Variable >& scp)
  : Constraint(scp) {}


void Mistral::ConstraintCliqueNotEqual::initialise() {
  Constraint::initialise();
  for(unsigned int i=0; i<scope.size; ++i)
    trigger_on(_value_, i);
  //set_idempotent(false);
  set_idempotent(true);
}

Mistral::ConstraintCliqueNotEqual::~ConstraintCliqueNotEqual() 
{
}

Mistral::PropagationOutcome Mistral::ConstraintCliqueNotEqual::propagate() 
{
  
// //
// //    

//   if(scope[0].id() == 0) {
//     std::cout << std::endl;
//     std::cout << "propagate " << this << std::endl;
//     std::cout << scope[0] << " " << scope[0].get_domain() ;
//     for(int x=0; x<scope.size; ++x) {
//       std::cout << ", " << scope[x] << " " << scope[x].get_domain() ;
//     }
//     std::cout << std::endl;
//     std::cout << "changes: " << changes << std::endl;
//     std::cout << "active: " << active << std::endl;
//   }

  unsigned int i, j, n=active.size, m;
  int active_var, value;
  //Event evt;
  for(i=0; i<changes.size; ++i) {
    //ground_var = changes[i];
    value = scope[changes[i]].get_min();
    
    for(j=i+1; j<changes.size; ++j) {
      if(scope[changes[j]].get_min() == value) {
	return FAILURE(changes[j]);
      } 
    }

    // since the set of active variables my change while
    // processing this loop, we do it backward
    for(j=n; j; --j) {
    //for(j=0; j<active.size; ++j) {
      active_var = active[j-1];
      //evt = scope[active_var].remove(value);
      if(scope[active_var].remove(value) == FAIL_EVENT) 
	return FAILURE(active_var);
    }
  }


  /// The following is to ensure idempotency
  m = active.size;
//   if(scope[0].id() == 0) {
//   std::cout << m << " . " << n << std::endl; 
//   }
  while(m < n) {
    for(i=m; i<n; ++i) {
      value = scope[active[i]].get_min();
      for(j=i+1; j<n; ++j)
	if(scope[active[j]].get_min() == value)
	  return FAILURE(active[j]);
      for(j=m; j; --j) {
	active_var = active[j-1];
	if(scope[active_var].remove(value) == FAIL_EVENT) 
	  return FAILURE(active_var);
      }
    }
    n = m;
    m = active.size;
  }

//   if(scope[0].id() == 0) {
//   std::cout << scope[0] << " " << scope[0].get_domain() ;
//   for(int x=1; x<scope.size; ++x) {
//     std::cout << ", " << scope[x] << " " << scope[x].get_domain() ;
//   }
//   std::cout << std::endl;
//   }

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

// std::ostream& Mistral::ConstraintCliqueNotEqual::display(std::ostream& os) const {
//   os << "cliqueNE(" << scope[0] ;
//   for(unsigned int i=1; i<scope.size; ++i) 
//     os << " ," << scope[i];
//   os << ")" ;
//   return os;
// }



/**********************************************
 * AllDiff Constraint 
 **********************************************/

const int INCONSISTENT = 0;
const int CHANGES      = 1;
const int NO_CHANGES   = 2;

Mistral::ConstraintAllDiff::ConstraintAllDiff(Vector< Variable >& scp)
  : Constraint(scp) { priority = 0; }


void Mistral::ConstraintAllDiff::initialise() {
  Constraint::initialise();
  for(unsigned int i=0; i<scope.size; ++i)
    trigger_on(_range_, i);
  set_idempotent(true);
  //priority = 0;
  init();
}


// void Mistral::ConstraintAllDiff::post() {
//   for(unsigned int i=0; i<scope.size; ++i)
//     trigger_on(_range_, i);
// }

void Mistral::ConstraintAllDiff::init() 
{
  unsigned int i;
  level = &(scope[0].get_solver()->level);
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


Mistral::ConstraintAllDiff::~ConstraintAllDiff() 
{
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
  int w,x,y,z,changes = 0;

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
      changes = 1;
    }
    if (d[z] == bounds[z]-bounds[y]) {
      pathset(h, h[y], j-1, y); // mark hall Interval
      h[y] = j-1; //("hall Interval [%d,%d)\n",bounds[j],bounds[y]);
    }
  }
  if( changes )
    return CHANGES;
  else
    return NO_CHANGES;
}


int Mistral::ConstraintAllDiff::filterupper()
{
  int i,j;
  int w,x,y,z,changes = 0;

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
      changes = 1;
    }
    if (d[z] == bounds[y]-bounds[z]) {
      pathset(h, h[y], j+1, y);
      h[y] = j+1;
    }
  }
  if( changes )
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
  if( lastLevel != ((*level) - 1) ) {
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

  lastLevel = *level;//(solver->level);
  //lastLevel = (solver->level);

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

// std::string Mistral::ConstraintAllDiff::getString() const 
// {    
//   std::string return_str = ("alldiff("+toString(scope[0]));
//   for(int i=1; i<scope.size; ++i) 
//     return_str += (" ,"+toString(scope[i]));
//   return (return_str+")");
// }

std::ostream& Mistral::ConstraintAllDiff::display(std::ostream& os) const {
  os << "alldiff(" << scope[0] ;
  for(unsigned int i=1; i<scope.size; ++i) 
    os << " ," << scope[i];
  os << ")" ;
  return os;
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Constraint& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Constraint* x) {
  return x->display(os);
}
