
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


// bool Mistral::Literal::check(Variable x) {
//   unsigned int mtype = _data_&2;
//   if(mtype == 0)
//     return ((_data_&1) ^ (x->)
// }

bool Mistral::Literal::apply(Variable x) const {
  unsigned int dat = __data__&3;
  //std::cout << dat << " " << __value__ << std::endl;
  if(dat == 1) {
    return x.set_domain(__value__);
  } else if(dat == 0) {
    return x.remove(__value__);
  } else if(dat == 3) {
    return x.set_max(__value__);
  } else {
    return x.set_min(__value__);
  }
}

bool Mistral::Literal::check(const int x) const {
  unsigned int mtype = __data__&2;
  if(mtype == 0)
    // __data__&1 == 1 -> equality
    return ((__data__&1) == (x==__value__));
  else
    // __data__&1 == 1 -> upper bound
    return ((__data__&1) == (x<=__value__));
}

bool Mistral::Literal::is_true(Variable x) const {
  unsigned int mtype = __data__&2;
  if(mtype == 0)
    // __data__&1 == 1 -> equality
    return ((__data__&1) == (x.equal(__value__)));
  else
    // __data__&1 == 1 -> upper bound
    return ((__data__&1) == (x.get_max()<=__value__));
}

bool Mistral::Literal::is_false(Variable x) const {
  unsigned int mtype = __data__&2;
  if(mtype == 0)
    // __data__&1 == 1 -> equality
    return ((__data__&1) != (x.contain(__value__)));
  else
    // __data__&1 == 1 -> upper bound
    return ((__data__&1) != (x.get_min()>__value__));
}

std::ostream& Mistral::Literal::display(std::ostream& os) const {
  unsigned int dat = __data__&3;
  if(dat == 0) {
    if(__value__==0) 
      os << "x" << get_atom() ;
    else if(__value__==1)
      os << "~x" << get_atom() ;
    else
      os << "x" << get_atom() << "!=" << __value__;
  } else if(dat == 1) {
      if(__value__==0) 
      os << "~x" << get_atom() ;
    else if(__value__==1)
      os << "x" << get_atom() ;
    else
      os << "x" << get_atom() << "=" << __value__;
  } else if(dat == 2) {
    os << "x" << get_atom() << "<" << __value__+1;
  } else {
    os << "x" << get_atom() << ">" << __value__;
  }
  return os;
}

Mistral::Constraint::Constraint() {
  stress = 1;
  priority = 2;
  id = -1;
}

Mistral::Constraint::Constraint(Vector< Variable >& scp) {
  stress = 1;
  priority = 2;
  id = -1;
  for(unsigned int i=0; i<scp.size; ++i) {
    scope.add(scp[i]);
  }
}

void Mistral::Constraint::add(Variable X) {
  scope.add(X);
}

Mistral::Constraint::~Constraint() {
  if(changes.list_ == events.list_)
    events.list_ = NULL;

  delete [] solution;
  delete [] self;
  delete [] trigger;
  delete [] event_type;
}


void Mistral::Constraint::mark_domain() {
}

void Mistral::Constraint::initialise() {
  trail_.initialise(0,2*scope.size);
  trail_.add(-1);

  event_type = new Event[scope.size];
  self = new int[scope.size];
  trigger = new int[scope.size];
  solution = new int[scope.size];
  changes.initialise(0, scope.size-1, false);
  supports = NULL;

  std::fill(event_type, event_type+scope.size, NO_EVENT);
  std::fill(self, self+scope.size, -1);
  std::fill(trigger, trigger+scope.size, _DOMAIN_);
  std::fill(solution, solution+scope.size, 0);

  active.initialise(0, scope.size-1, true);

  
//   for(unsigned int i=0; i<scope.size; ++i) 
//     if(scope[i].domain_type != BOOL_VAR && scope[i].is_ground()) 
//       active.erase(i);
}

Mistral::PropagationOutcome Mistral::Constraint::rewrite() {
  return propagate();
}

void Mistral::Constraint::consolidate() {

//   std::cout << "consolidate " ;
//   display(std::cout);
//   std::cout << std::endl;

  for(unsigned int i=0; i<scope.size; ++i) {
    scope[i] = scope[i].get_var();
  }

  //display(std::cout);
  //std::cout << " consolidated" << std::endl;
}



void Mistral::Constraint::post(Solver *s) {
  solver = s;

  active.fill();

  unsigned int i, j;
  // for each of its variables
  for(i=0; i<scope.size; ++i) {
    if(!scope[i].is_ground()) {
      // add the constraint to the list of constraints on that variable
      j = scope[i].id();
      ConstraintTrigger ct(this, i);
      int trg = trigger[i];
      ConstraintList *lst = solver->constraint_graph[j];
      unsigned int elt = lst->reversible_add(ct, trg);
      self[i] = elt;
    } else {
      active.erase(i);
    }
//   for(unsigned int i=0; i<scope.size; ++i) 
//     if(scope[i].domain_type != BOOL_VAR && scope[i].is_ground()) 
//       active.erase(i);
  }


//   std::cout << "MARK DOMAIN " ;
//   display(std::cout);
//   std::cout << std::endl;
//   for(i=0; i<scope.size; ++i) {
//     std::cout << solver->domain_types[scope[i].id()] << " ";
//   }
//   std::cout << std::endl;

  mark_domain();

//   for(i=0; i<scope.size; ++i) {
//     std::cout << solver->domain_types[scope[i].id()] << " ";
//   }
//   std::cout << std::endl;

//   if(scope[0].id() == 5 || scope[1].id() == 5) {
//     std::cout << "X5CLIST AFTER POST ";
//     solver->print_clist(5);
//     std::cout << "C" << id << "ALIST AFTER POST "
// 	      << active << " " << trail_ << std::endl;
//   }
}

void Mistral::Constraint::relax() {
  if(!active.empty()) {
    unsigned int i, j;
    int k;
    
#ifdef _DEBUG_AC
    for(i=0; i<scope.size; ++i) {
      std::cout << " " << scope[i] << " in " << scope[i].get_domain();
    }
    std::cout << std::endl << "relax " << this << " from";
#endif
      
    for(i=0; i<active.size; ++i) {
      j = active[i];
      k = scope[j].id();
      
#ifdef _DEBUG_AC
      std::cout << " " << scope[j] << "'s c-list ("; 
      
      bool is_in = false;
      ConstraintNode nd;
      nd = solver->constraint_graph[k]->first(_VALUE_);
      while( !is_in && solver->constraint_graph[k]->next(nd) ) {
	std::cout << " " << nd.elt.constraint;
	is_in = (nd.elt.constraint->id == id);
      }
      
      if(!is_in) {
	std::cout << ") WAS NOT IN!" << std::endl;
	exit(1);
      }
#endif

      //std::cout << "INDEX " << k << " " << j << std::endl;
      //std::cout << "ERASE FROM " << scope[j] << std::endl;
      solver->constraint_graph[k]->reversible_erase(self[j], trigger[j]);
//       std::cout << "TEST" << std::endl;
//       if(solver->constraint_graph[k]->contain(self[j])) {
// 	exit(1);
//       }
    }
    
#ifdef _DEBUG_AC
    std::cout << " )" << std::endl;
#endif

    if(trail_.back() != solver->level) {
      solver->saved_cons.add(this);
      trail_.add(active.size);
      trail_.add(solver->level);
    }
    active.clear();
  }//  else {
    
//     std::cout << "do not relax because there is no active variable" 
// 	      << std::endl << scope[0] << " in " << scope[0].get_domain()
// 	      << " " << scope[1] << " in " << scope[1].get_domain() << std::endl;
//   }

//   if(scope[0].id() == 5 || scope[1].id() == 5) {
//     std::cout << "X5CLIST AFTER RELAX ";
//     solver->print_clist(5);
//     std::cout << "C" << id << "ALIST AFTER RELAX "
// 	      << active << " " << trail_ << std::endl;
//     //<< solver->constraint_graph[5] << std::endl;
//   }

//   for(unsigned int i=0; i<scope.size; ++i) {
//     if(!scope[i].is_ground()) {
//       int j = scope[i].id();
//       std::cout << "INDEX " << j << " " << i << std::endl;
//       if(solver->constraint_graph[j]->contain(self[i])) {
// 	std::cout << std::endl << this << " (" << self[i] << ") is in " << scope[i] << std::endl;
// 	std::cout << solver->constraint_graph[j]->data.stack_[self[i]].prev << std::endl;
// 	std::cout << solver->constraint_graph[j]->data.stack_[self[i]].next << std::endl;

// 	// std::cout << solver->constraint_graph[j]->data.stack_[self[i]].prev.next
// // 		  << " " << solver->constraint_graph[j]->data.stack_[self[i]].next.prev 
// //		  << std::endl;

// 	ConstraintNode nd;
// 	std::cout << std::endl;
// 	nd = solver->constraint_graph[j]->first(_VALUE_);
// 	while( solver->constraint_graph[j]->next(nd) ) {
// 	  std::cout << " " << nd.elt.constraint;
// 	}
// 	std::cout << std::endl;
// 	nd = solver->constraint_graph[j]->first(_VALUE_);
// 	while( solver->constraint_graph[j]->next(nd) ) {
// 	  std::cout << " " << nd.prev;
// 	}
// 	std::cout << std::endl;
// 	exit(1);
//       }
//     }
//   }

}

// Mistral::Constraint* Mistral::Constraint::notify_assignment(const int var, const int level) {
//   Constraint *r = NULL;
//   if(trail_.back() != level) {
//     trail_.add(active.size);
//     trail_.add(level);
//     r = this;
//   }
//   active.erase(var);
//   return r;
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

std::ostream& Mistral::Constraint::display(std::ostream& os) const {
  os << name() << "(" << scope[0];
  for(unsigned int i=1; i<scope.size; ++i)
    os << ", " << scope[i];
  os << ")";
  return os;
}

void Mistral::ConstraintNotEqual::initialise() {
  Constraint::initialise();
  trigger_on(_VALUE_, 0);
  trigger_on(_VALUE_, 1);
  set_idempotent(true);
}

void Mistral::ConstraintNotEqual::mark_domain() {
  solver->mark_non_convex(scope[0].id());
  solver->mark_non_convex(scope[1].id());
}

Mistral::PropagationOutcome Mistral::ConstraintNotEqual::propagate() {
  int var = 1-changes[0];
  if(active.size)
    return(scope[var].remove(scope[1-var].get_min()) == FAIL_EVENT ? FAILURE(var) : CONSISTENT);
  else
    return(scope[0].get_min() == scope[1].get_min() ? FAILURE(var) : CONSISTENT);
}

std::ostream& Mistral::ConstraintNotEqual::display(std::ostream& os) const {
  os << scope[0] << " =/= " << scope[1];
  return os;
}

void Mistral::ConstraintEqual::initialise() {
  Constraint::initialise();
  trigger_on(_DOMAIN_, 0);
  trigger_on(_DOMAIN_, 1);
  set_idempotent(true);
}

Mistral::PropagationOutcome Mistral::ConstraintEqual::rewrite() {
  Mistral::PropagationOutcome wiped = propagate();


 

  if( active.size == 2 ) {
    if( scope[0].is_expression() && scope[1].is_expression() ) {


//  std::cout << "REWRITE EQUALITY " ;
//   display(std::cout);
//   std::cout << std::endl;
//       std::cout << solver->variables << std::endl;


      relax();

      int k[2], j;

      k[0] = scope[0].id();
      k[1] = scope[1].id();

      j = solver->constraint_graph[k[0]]->degree > 
	solver->constraint_graph[k[1]]->degree;
      // j is 1 iff degree(x0) > degree(x1), 
      // in that case we want to transfer x1's constraints on x0.
      // x[j] <- x[1-j]
      
      ConstraintNode nd;
      nd = solver->constraint_graph[k[j]]->first(_VALUE_);
      while( solver->constraint_graph[k[j]]->next(nd) ) {	
	nd.elt.constraint->relax();
	nd.elt.constraint->scope[nd.elt.index] = scope[1-j];
	solver->add(nd.elt.constraint);
      }
      
      //and now scope[j] points to scope[1-j]
      scope[j].expression->self = scope[1-j];
      scope[j].expression->id = scope[1-j].id();

//       std::cout << solver->variables << std::endl;
//       std::cout << "EQUALITY REWRITEN " ;
//       display(std::cout);
//       std::cout << std::endl;

    }
  }

  return wiped;
}


Mistral::PropagationOutcome Mistral::ConstraintEqual::propagate() {

  unsigned int i;
  int var;

  PropagationOutcome wiped = CONSISTENT;
  if(active.size) {
    for(i=0; IS_OK(wiped) && i<changes.size; ++i) {
      var = changes[i];
      if(scope[1-var].set_domain(scope[var]) == FAIL_EVENT)
	wiped = FAILURE(1-var);
    }
  } else if(scope[0].get_min() != scope[1].get_min()) {
    wiped = FAILURE(0);
  }

  return wiped;
}


std::ostream& Mistral::ConstraintEqual::display(std::ostream& os) const {
  os << scope[0] << " == " << scope[1];
  return os;
}

void Mistral::PredicateUpperBound::initialise() {
  Constraint::initialise();
  trigger_on(_RANGE_, 0);
  trigger_on(_VALUE_, 1);
  set_idempotent(false);
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


void Mistral::PredicateLowerBound::initialise() {
  Constraint::initialise();
  trigger_on(_RANGE_, 0);
  trigger_on(_VALUE_, 1);
  set_idempotent(false);
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

void Mistral::PredicateLess::initialise() {
  Constraint::initialise();
  trigger_on(_RANGE_, 0);
  trigger_on(_RANGE_, 1);
  trigger_on(_VALUE_, 2);
  set_idempotent(false);
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

void Mistral::ConstraintLess::initialise() {
  Constraint::initialise();
  trigger_on(_RANGE_, 0);
  trigger_on(_RANGE_, 1);
  set_idempotent(true);
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

  if(changes.contain(0) && LB_CHANGED(event_type[0])) {
    if(scope[1].set_min(scope[0].get_min() + offset) == FAIL_EVENT) 
      return FAILURE(1);
  }
  if(changes.contain(1) && UB_CHANGED(event_type[1])) {
    if(scope[0].set_max(scope[1].get_max() - offset) == FAIL_EVENT) 
      return FAILURE(0);
  }

//   if(scope[0].id() == 6 && scope[1].id() == 9) {
//     for(unsigned int i=0; i<scope.size; ++i)
//       std::cout << " " << scope[i].get_domain();
//     std::cout << std::endl;
//   }

  return CONSISTENT;
}

std::ostream& Mistral::ConstraintLess::display(std::ostream& os) const {
  os << scope[0];
  if(offset < 0) os << " - " << (-offset+1) << " < ";
  else if(offset > 1) os << " + " << (offset-1) << " < ";
  else if(offset > 0) os << " < ";
  else os << " <= ";
  os << scope[1];
  return os;
}

Mistral::ConstraintDisjunctive::ConstraintDisjunctive(Vector< Variable >& scp, const int p0, const int p1) : Constraint(scp) { 
  processing_time[0] = p0; 
  processing_time[1] = p1;
  
  precedence[0] = new ConstraintLess(scope, processing_time[0]);
  precedence[1] = new ConstraintLess(processing_time[1]);
  precedence[1]->add(scope[1]);
  precedence[1]->add(scope[0]);
}

void Mistral::ConstraintDisjunctive::initialise() {
  Constraint::initialise();
  trigger_on(_RANGE_, 0);
  trigger_on(_RANGE_, 1);
  set_idempotent(true);
  stress = 0;
}

void Mistral::ConstraintDisjunctive::decide(const int choice) {

  //std::cout << "decide " << this << " => " ;


  relax();
  

  if(choice==1) {
    solver->add(precedence[1]);
    //std::cout << precedence[1] ;
  } else {
    solver->add(precedence[0]);
    //std::cout << precedence[0] ;
  }

  //std::cout << std::endl;

//   std::cout << solver->active_constraints << std::endl;

}

Mistral::PropagationOutcome Mistral::ConstraintDisjunctive::propagate() {
  //(x0 + p0 <= x1 || x1 + p1 <= x0).
  int hold = 3;

 //  if(scope[0].id() == 6 && scope[1].id() == 9) {
//     std::cout << "DISJUNCTIVE " << this << std::endl;
    
//  for(unsigned int i=0; i<scope.size; ++i)
//    std::cout << " " << scope[i] << " in " << scope[i].get_domain();
//   std::cout << std::endl;

//     std::cout << scope[1].get_min() << " + " << processing_time[1] 
// 	      << " ?> " << scope[0].get_max() ;
//   }

  // check is prec[1] is violated (x1 + p1 > x0).
  if(scope[1].get_min()+processing_time[1] > scope[0].get_max()) {
    hold &= 2;

//     if(scope[0].id() == 6 && scope[1].id() == 9) {
//       std::cout << " YES" << std::endl;
//     }
  } // else   if(scope[0].id() == 6 && scope[1].id() == 9) {    std::cout << " NO" << std::endl;
//   }

//   if(scope[0].id() == 6 && scope[1].id() == 9) {
//     std::cout << scope[0].get_min() << " + " << processing_time[0] 
// 	      << " ?> " << scope[1].get_max() ;
//   }

  // check is prec[1] is violated (x0 + p0 > x1).
  if(scope[0].get_min()+processing_time[0] > scope[1].get_max()) {
    hold &= 1;
//     if(scope[0].id() == 6 && scope[1].id() == 9) {
//       std::cout << " YES" << std::endl;
//     }
  } //  else  if(scope[0].id() == 6 && scope[1].id() == 9) {     std::cout << " NO" << std::endl;
//   }

  if(!hold) return FAILURE(0);
  if(hold<3) {
    decide(hold);
  }
  return CONSISTENT;
}

void Mistral::ConstraintDisjunctive::consolidate() {
  for(unsigned int i=0; i<scope.size; ++i) {
    scope[i] = scope[i].get_var();
    precedence[0]->scope[i] = scope[i];
    precedence[1]->scope[1-i] = scope[i];
  }
}

std::ostream& Mistral::ConstraintDisjunctive::display(std::ostream& os) const {
  os << precedence[0] << " or " 
     << precedence[1] ;
  return os;
}

void Mistral::PredicateEqual::initialise() {
  Constraint::initialise();
  trigger_on(_DOMAIN_, 0);
  trigger_on(_DOMAIN_, 1);
  trigger_on(_VALUE_, 2);
  set_idempotent(false);
}

void Mistral::PredicateEqual::mark_domain() {
  if(!spin) {
    solver->mark_non_convex(scope[0].id());
    solver->mark_non_convex(scope[1].id());
  }
}

Mistral::PropagationOutcome Mistral::PredicateEqual::rewrite() {
  Mistral::PropagationOutcome wiped = propagate();
  if( active.size == 2 ) {
    relax();

    if( scope[2].is_ground() ) {
      if( (spin + scope[2].get_min()) != 1 ) {
	solver->add(new ConstraintNotEqual(scope));
      } else {
	solver->add(new ConstraintNotEqual(scope));
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

  return wiped;
}

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


void Mistral::PredicateConstantEqual::initialise() {
  Constraint::initialise();
  trigger_on(_DOMAIN_, 0);
  trigger_on(_VALUE_, 1);
  set_idempotent(false);
}

void Mistral::PredicateConstantEqual::mark_domain() {
  if(!spin && 
     value > scope[0].get_min() && 
     value < scope[0].get_max()) {
    solver->mark_non_convex(scope[0].id());
  }
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


void Mistral::PredicateOffset::initialise() {
  Constraint::initialise();
  for(int i=0; i<2; ++i)
    trigger_on(_RANGE_, i);
  set_idempotent(false);
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

void Mistral::PredicateNot::initialise() {
  Constraint::initialise();
  for(int i=0; i<2; ++i)
    trigger_on(_VALUE_, i);
  set_idempotent(false);
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


void Mistral::PredicateAnd::initialise() {
  Constraint::initialise();
  for(int i=0; i<3; ++i)
    trigger_on(_VALUE_, i);
  set_idempotent(false);
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


void Mistral::PredicateOr::initialise() {
  Constraint::initialise();
  for(int i=0; i<3; ++i)
    trigger_on(_VALUE_, i);
  set_idempotent(false);
}

Mistral::PropagationOutcome Mistral::PredicateOr::rewrite() {
  Mistral::PropagationOutcome wiped = propagate();
  if( scope[2].is_ground() && active.size == 2 ) {
    relax();
    if( scope[2].get_min() ) {
      solver->add(new ConstraintOr(scope));
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
  os << scope[2] << " <=> (" << scope[0] << " or " << scope[1] << ")";
  return os;
}


void Mistral::ConstraintAnd::initialise() {
  Constraint::initialise();
  for(int i=0; i<2; ++i)
    trigger_on(_RANGE_, i);
  set_idempotent(false);
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

void Mistral::ConstraintOr::initialise() {
  Constraint::initialise();
  for(int i=0; i<2; ++i)
    trigger_on(_RANGE_, i);
  set_idempotent(true);
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

std::ostream& Mistral::ConstraintOr::display(std::ostream& os) const {
  os << "(" << scope[0] << " or " << scope[1] << ")";
  return os;
}

void Mistral::PredicateAdd::initialise() {
  Constraint::initialise();
  trigger_on(_RANGE_, 0);
  trigger_on(_RANGE_, 1);
  trigger_on(_RANGE_, 2);
  set_idempotent(true);
}

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
  priority = 1;
  total = t; 
}

void Mistral::ConstraintBoolSumEqual::initialise() {
  Constraint::initialise();
  for(unsigned int i=0; i<scope.size; ++i)
    trigger_on(_VALUE_, i);
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
  for(unsigned int i=0; i<scope.size; ++i) {
    trigger_on(_VALUE_, i);
  }
  set_idempotent(true);
}

void Mistral::ConstraintCliqueNotEqual::mark_domain() {
  for(unsigned int i=0; i<scope.size; ++i) {
    solver->mark_non_convex(scope[i].id());
  }
}

Mistral::ConstraintCliqueNotEqual::~ConstraintCliqueNotEqual() 
{
}

Mistral::PropagationOutcome Mistral::ConstraintCliqueNotEqual::propagate() 
{
  unsigned int i, j, n=active.size, m;
  int active_var, value;
  for(i=0; i<changes.size; ++i) {
    value = scope[changes[i]].get_min();
    for(j=i+1; j<changes.size; ++j) {
      if(scope[changes[j]].get_min() == value) {
	return FAILURE(changes[j]);
      } 
    }

    // since the set of active variables my change while
    // processing this loop, we do it backward
    for(j=n; j; --j) {
      active_var = active[j-1];
      if(scope[active_var].remove(value) == FAIL_EVENT) 
	return FAILURE(active_var);
    }
  }

  /// The following is to ensure idempotency
  m = active.size;
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
  for(unsigned int i=0; i<scope.size; ++i) {
    trigger_on(_RANGE_, i);
  }
  set_idempotent(true);
  init();
}

void Mistral::ConstraintAllDiff::mark_domain() {
  for(unsigned int i=0; i<scope.size; ++i) {
    solver->mark_non_convex(scope[i].id());
  }
}


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

std::ostream& Mistral::ConstraintAllDiff::display(std::ostream& os) const {
  os << "alldiff(" << scope[0] ;
  for(unsigned int i=1; i<scope.size; ++i) 
    os << " ," << scope[i];
  os << ")" ;
  return os;
}


// Mistral::ConstraintNogoodBase::ConstraintNogoodBase(Vector< Variable >& scp) 
//   : Constraint(scp) { }

// void Mistral::ConstraintNogoodBase::mark_domain() {
//   for(unsigned int i=0; i<scope.size; ++i) {
//     solver->mark_non_convex(scope[i].id());
//   }
// }

// void Mistral::ConstraintNogoodBase::initialise() {
//   Constraint::initialise();
//   //int min_idx = INFTY;
//   //int max_idx = -INFTY;
//   for(unsigned int i=0; i<scope.size; ++i) {
//     //if(scope[i].get_id() > max_idx) max_idx = scope[i].get_id();
//     //if(scope[i].get_id() < min_idx) min_idx = scope[i].get_id();
//     trigger_on(_VALUE_, i);
//   }
//   set_idempotent(true);

//   //if(// min_idx < 
//   //max_idx) {
//   watch_structure.initialise(0,scope.size);
//     //}

//   for(unsigned int i=0; i<scope.size; ++i) {
//     //add(scope[i]);
//     watch_structure[i] = new Vector< Array< Literal >* > [scope[i].get_size()];
//     watch_structure[i]-=scope[i].get_min();
//   }
// }

// Mistral::ConstraintNogoodBase::~ConstraintNogoodBase() {
// }

// void Mistral::ConstraintNogoodBase::add(Variable x) {
//   unsigned int idx = x.id();
//   if(idx >= scope.size) {
//     while(scope.capacity <= idx)
//       scope.extendStack();
//     scope[idx] = x;
//     while(watch_structure.capacity <= idx)
//       watch_structure.extendStack();
//   }
// }

// void Mistral::ConstraintNogoodBase::add( Vector < Literal >& clause ) {
//  if(clause.size > 1) {
//    Array<Literal> *cl = Array<Literal>::Array_new(clause);
//    nogood.add( cl );
//    // std::cout << clause[0].get_atom() << " " << (1-clause[0].get_value()) << std::endl;
//    // std::cout << watch_structure[clause[0].get_atom()][1-clause[0].get_value()] << std::endl;

//    watch_structure[clause[0].get_atom()][1-clause[0].get_value()].add(cl);
//    watch_structure[clause[1].get_atom()][1-clause[1].get_value()].add(cl);
//  } else {
//    clause[0].apply(scope[clause[0].get_atom()]);
//  }
// }

// int Mistral::ConstraintNogoodBase::check( const int* sol ) const {
//   unsigned int i, j;
//   bool satisfied=true;

//   for(i=0; i<nogood.size; ++i) {
//     satisfied = true;
//     Array < Literal >& clause = *(nogood[i]);
//     for(j=0; satisfied && j<clause.size; ++j) {
//       satisfied = clause[j].check(sol[j]);
//     }
//   }
  
//   return !satisfied;
// }

// Mistral::PropagationOutcome Mistral::ConstraintNogoodBase::propagate() {
//   Array< Literal > *conflict=NULL;

//   int x, v, cw;
//   Literal p;
//   //unsigned int i;
//   while( !changes.empty() ) {
//     std::cout << changes << std::endl;
//     //std::cout << scope << std::endl;
//     x = changes.pop();
//     //std::cout << scope[x] << " in " << scope[x].get_domain() << std::endl;
//     v = scope[x].get_min();
//     std::cout << x << "=" << v << ": " << watch_structure[x][v] << std::endl;

//     p = Literal(x,EQ,v);


//     cw = watch_structure[x][v].size;
//     while(cw-- && !conflict) {
//       conflict = update_watcher(cw, x, v);
//     }
//   }

//   return CONSISTENT;
// }


// inline Array< Literal >* SatSolver::update_watcher(const int cw, const int x, const int v)
// {
//   Array< Literal > *cl = watch_structure[x][v][cw];
//   Array< Literal >& clause = *cl;


//   unsigned int j;
//   Lit q, r;
//   Atom v, w;

// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
//   std::cout << " ";
//   print_clause( std::cout, watch_structure[p][cw] );
//   std::cout << std::endl;
// #endif

//   //ensure that p is the second watched lit
//   if( clause[1] != p ) {
//     q = clause[1];
//     clause[0] = q;
//     clause[1] = p;
//   } else q = clause[0];
//   v=state[UNSIGNED(q)];

// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
//    std::cout << " second watched: ";
//    printLiteral(std::cout, q);
//    std::cout << std::endl;

//    if( LEVEL(v) < assumptions.size && SIGN(v) == SIGN(q) ) {
//      for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
//      std::cout << " satisfied by second watched" << std::endl;
//    }
// #endif

//   //check if the other watched lit is assigned
//   if( LEVEL(v) >= assumptions.size || SIGN(v) != SIGN(q) ) {
//     for(j=2; j<clause.size; ++j) {
//       // for each literal q of the clause,
//       r = clause[j];
//       w = state[UNSIGNED(r)];

// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
//        std::cout << "\tcheck another lit: ";
//        printLiteral(std::cout, r);
//        std::cout << std::endl;
// #endif

//       if( LEVEL(w) >= assumptions.size ) { // this literal is not set
// 	// then it is a good candidate to replace p
// 	clause[1] = r;
// 	clause[j] = p;
// 	watch_structure[p].remove(cw);
// 	watch_structure[r].add(cl);

// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// 	std::cout << " new watched" << std::endl;
// #endif

// 	break;	
//       }
//       // if it is set true, then the clause is satisfied
//       else if( SIGN(w) == SIGN(r) ) {

// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// 	std::cout << " satisfied" << std::endl;
// #endif

// 	break;
//       }
//     }
      
//     if( j == clause.size ) // no replacement could be found
//       { 

// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// 	std::cout << " no more watched" << std::endl;
// #endif

// 	if( LEVEL(v) >= assumptions.size ) {
// 	  // the last literal (other watched lit) is not set yet, we set it
// 	  add_lit(q);
// 	  reason[UNSIGNED(q)] = cl;

// 	  //#ifdef _DEBUGNOGOOD
// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// 	  std::cout //<< "unit prune disjunct b" << y 
// 	    << " because ";
// 	  print_clause( std::cout, cl );
// 	  std::cout << std::endl;
// #endif

// 	  //	  std::cout << " prop " << q << std::endl;
	  
// 	} else 
// 	  // it is set to false already, we fail
// 	  if( // polarity[y] == -q
// 	     SIGN(v) != SIGN(q)
// 	      ) {

// #ifdef _DEBUG_WATCH
//   for(unsigned int i=0; i<decisions.size; ++i) std::cout << " " ;
// 	    std::cout << " fail" << std::endl;
// #endif

// 	    return cl;
// 	  }
//       }
//   }

//   return NULL;
// }




// std::ostream& Mistral::ConstraintNogoodBase::display(std::ostream& os) const {
//   os << " (";
//   if(nogood.size>0) {
//     os << nogood[0];
//     for(unsigned int i=1; i<nogood.size; ++i)
//       os << " " << nogood[i]  ;
//   }
//   os << ")";
//   return os;
// }



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


  // for(unsigned int i=0; i<scope.size; ++i) {
  //   is_watched_by[i] = new Vector< Clause* >;
  //   is_watched_by[i]-=scope[i].get_min();
  // }
}

Mistral::ConstraintClauseBase::~ConstraintClauseBase() {
}

void Mistral::ConstraintClauseBase::add(Variable x) {
  unsigned int idx = x.id();
  if(idx == scope.size) {
    scope.add(x);
    while(is_watched_by.capacity <= 2*idx)
      is_watched_by.extendStack();
  } else if(idx > scope.size) {
    while(scope.capacity <= idx)
      scope.extendStack();
    scope[idx] = x;
    while(is_watched_by.capacity <= idx)
      is_watched_by.extendStack();
  }
}

void Mistral::ConstraintClauseBase::add( Vector < Literal >& clause ) {
 if(clause.size > 1) {
   Clause *cl = Clause::Array_new(clause);
   clause.add( cl );
   is_watched_by[clause[0]].add(cl);
   is_watched_by[clause[1]].add(cl);
 } else {
   scope[UNSIGNED(clause[0])].set_domain(SIGN(clause[0]));
 }
}

int Mistral::ConstraintClauseBase::check( const int* sol ) const {
  unsigned int i, j;
  bool satisfied=true;

  for(i=0; i<clause.size; ++i) {
    satisfied = true;
    Clause& clause = *(clause[i]);
    for(j=0; satisfied && j<clause.size; ++j) {
      satisfied = //clause[j].check(sol[j]);
	(sol[j] == SIGN(clause[j]));
    }
  }
  
  return !satisfied;
}

Mistral::PropagationOutcome Mistral::ConstraintClauseBase::propagate() {
  Clause *conflict=NULL;


  int x, v, cw;
  Lit p;
  //unsigned int i;
  while( !changes.empty() ) {
    std::cout << changes << std::endl;
    //std::cout << scope << std::endl;
    x = changes.pop();
    //std::cout << scope[x] << " in " << scope[x].get_domain() << std::endl;
    v = scope[x].get_min();
    std::cout << x << "=" << v << ": " << is_watched_by[x][v] << std::endl;

    p = 2*x+v;

    cw = is_watched_by[p].size;
    while(cw-- && !conflict) {
      conflict = update_watcher(cw, p);
    }
  }

  return CONSISTENT;
}


inline Clause* SatSolver::update_watcher(const int cw, const Lit p)
{
  Clause *cl = is_watched_by[p][cw];
  Clause& clause = *cl;

  unsigned int j;
  Lit q, r;
  //Atom v, w;
  Variable v, w;

  //ensure that p is the second watched lit
  if( clause[1] != p ) {
    q = clause[1];
    clause[0] = q;
    clause[1] = p;
  } else q = clause[0];
  v=scope[UNSIGNED(q)].;


  //check if the other watched lit is assigned
  //if( LEVEL(v) >= assumptions.size || SIGN(v) != SIGN(q) ) {
  if( !v.is_ground() || v.get_min() != SIGN(q) ) {
    for(j=2; j<clause.size; ++j) {
      // for each literal q of the clause,
      r = clause[j];
      w = scope[UNSIGNED(r)];

      if( !w.is_ground() ) { // this literal is not set
	// then it is a good candidate to replace p
	clause[1] = r;
	clause[j] = p;
	is_watched_by[p].remove(cw);
	is_watched_by[r].add(cl);

	break;	
      }
      // if it is set true, then the clause is satisfied
      else if( w.get_min() == SIGN(r) ) {
	break;
      }
    }
      
    if( j == clause.size ) // no replacement could be found
      { 
	if( !v,is_ground() ) {
	  // the last literal (other watched lit) is not set yet, we set it
	  add_lit(q);
	  //reason[UNSIGNED(q)] = cl;
	} else 
	  // it is set to false already, we fail
	  if( v.get_min() != SIGN(q) ) {
	    return cl;
	  }
      }
  }

  return NULL;
}




std::ostream& Mistral::ConstraintClauseBase::display(std::ostream& os) const {
  os << " (";
  if(clauses.size>0) {
    os << clauses[0];
    for(unsigned int i=1; i<clauses.size; ++i)
      os << " " << clauses[i]  ;
  }
  os << ")";
  return os;
}


std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Constraint& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Constraint* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Literal& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Literal* x) {
  return x->display(os);
}
