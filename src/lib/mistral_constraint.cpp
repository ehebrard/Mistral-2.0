
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
  is_posted = false;
  stress = 1;
  priority = 2;
  id = -1;
}

Mistral::Constraint::Constraint(Vector< Variable >& scp) {
  is_posted = false;
  stress = 1;
  priority = 2;
  id = -1;
  for(unsigned int i=0; i<scp.size; ++i) {
    scope.add(scp[i]);
  }
}

Mistral::Constraint::Constraint(std::vector< Variable >& scp) {
  is_posted = false;
  stress = 1;
  priority = 2;
  id = -1;
  for(unsigned int i=0; i<scp.size(); ++i) {
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
//       active.remove(i);
}

void Mistral::Constraint::restore() {

  //std::cout << "restore " << this << " " << is_posted << std::endl;

  // if(active.size && active.size <= stress) {
  //   solver->notify_post(this);
  // }

  if(!is_posted) {
    is_posted = true;
    solver->notify_post(this);
  }

  trail_.pop();
  active.size = trail_.pop();

  //std::cout << " -> " << active << std::endl;

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

int Mistral::Constraint::get_backtrack_level() {
  return solver->level-1;
}

Mistral::Decision Mistral::Constraint::get_decision() { 
  Decision dec = solver->decisions.back(0); 
  dec.invert();
  return dec;
}

void Mistral::Constraint::post(Solver *s) {
  solver = s;
  is_posted = true;

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
      active.remove(i);
    }
//   for(unsigned int i=0; i<scope.size; ++i) 
//     if(scope[i].domain_type != BOOL_VAR && scope[i].is_ground()) 
//       active.remove(i);
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

  //std::cout << "try to relax " << this ;


  if(is_posted && !active.empty()) {

    //std::cout << " ok" << std::endl;

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
      //std::cout << "REMOVE FROM " << scope[j] << std::endl;
      solver->constraint_graph[k]->reversible_remove(self[j], trigger[j]);
//       std::cout << "TEST" << std::endl;
//       if(solver->constraint_graph[k]->contain(self[j])) {
// 	exit(1);
//       }
    }
    
#ifdef _DEBUG_AC
    std::cout << " )" << std::endl;
#endif

    // need to save the constraint so that it is restored
    if(trail_.back() != solver->level) {
      solver->saved_cons.add(this);
      trail_.add(active.size);
      trail_.add(solver->level);
    }

    //std::cout << "relax " << this << std::endl;

    solver->notify_relax(this);

    is_posted = false;
    

    //active.clear();

  }//   else {
  //   std::cout << " not posted!" << std::endl;
  // }   
 
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
//   // std::cout << "remove " << scope[var] << " from " << this //<< std::endl;
//   // 	    << " " << active << " -> "; 
  
//   Constraint *r = NULL;
//   if(trail_.back() != level) {
//     trail_.add(active.size);
//     trail_.add(level);
// 	r = this;
//   }
  
//   active.remove(var);
  
//   //std::cout << active << std::endl;
  
//   return r;
// }

// Mistral::Constraint* Mistral::Constraint::notify_assignment(const int var, const int level) {
//   Constraint *r = NULL;
//   if(trail_.back() != level) {
//     trail_.add(active.size);
//     trail_.add(level);
//     r = this;
//   }
//   active.remove(var);
//   return r;
// }

Mistral::PropagationOutcome Mistral::Constraint::bound_propagate() {
  // PropagationOutcome wiped = CONSISTENT;
  
  // unsigned int i, j;

  // // for( j=0; j<scope.size; ++j) {
  // //   std::cout << scope[j] << " in " << scope[j].get_domain() << " ";
  // // }
  // // std::cout << std::endl;

  // int vali, valmax;//, vnext;
  // bool supported;
  // for( j=0; j<changes.size; ++j) {

  //   //std::cout << "change on " << scope[changes[j]] << std::endl;

  //   for( i=0; IS_OK(wiped) && i<scope.size; ++i ) {

  //     //std::cout << scope[i] << "?" << std::endl; 

  //     if( (int)i != changes[j] ) { 

  // 	supported = false;
  // 	vali = scope[i].get_min();
  // 	valmax = scope[i].get_max();

  // 	while(!supported && IS_OK(wiped) && vali<=valmax) {
  // 	  //std::cout << "find support for " << scope[i] << "=" << vali << std::endl;

  // 	  if( ( !first_support(i, vali) && !find_bound_support(i, vali) ) ) {

  // 	    //std::cout << "\tNo support found, remove!" << std::endl;

  // 	    if(IS_FAIL(scope[i].remove(vali))) {
  // 	      //std::cout << "FAIL!!" << std::endl;
  // 	      wiped = FAILURE(i);
  // 	    } 
  // 	  } else supported = true;

  // 	  ++vali;
  // 	}

  // 	if(supported && vali<=valmax) {
  // 	  supported = false;
  // 	  vali = valmax;

  // 	  while(!supported && IS_OK(wiped)) {
  // 	    //std::cout << "find support for " << scope[i] << "=" << vali << std::endl;
	    
  // 	    if( ( !first_support(i, vali) && !find_bound_support(i, vali) ) ) {
	      
  // 	      //std::cout << "\tNo support found, remove!" << std::endl;
	      
  // 	      if(IS_FAIL(scope[i].remove(vali))) {
  // 		//std::cout << "FAIL!!" << std::endl;
  // 		wiped = FAILURE(i);
  // 	      } 
  // 	    } else supported = true;
  // 	    --vali;
  // 	  }
  // 	}
  //     }	
  //   }
  // }
  
  // return wiped;


  PropagationOutcome wiped = CONSISTENT;
  
  unsigned int i, j;

  int vali, valmax;//, vnext;
  bool supported;
  while(!changes.empty()) {
    j = changes.pop();
    for( i=0; IS_OK(wiped) && i<scope.size; ++i ) {
      if( i != j ) { 
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
  }
  
  return wiped;
}

Mistral::PropagationOutcome Mistral::Constraint::propagate() {
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

	  // if(solver->parameters.verbosity) {
	  //   std::cout << "find a support for " << scope[i] << " = " << vali << std::endl;
	  // }

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

bool Mistral::Constraint::first_support(const int vri, const int vli) 
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

bool Mistral::Constraint::find_support(const int vri, const int vli) 
{
  int i=scope.size, vali;
  bool found=false;
  // sol is initialized: either to the value 
  // a variable is already assigned to
  // or to the first value in its domain
  while(i >= 0) {

    // if(solver->parameters.verbosity) {
    //   std::cout << "\t<" << solution[0] ;
    //   for(unsigned int k=1; k<scope.size; ++k) {
    // 	std::cout << " " << solution[k];
    //   }
    //   std::cout << "> ";
    // }

    // check this assignment
    if( !check( solution ) ) {

      // if(solver->parameters.verbosity) {
      // 	std::cout << "OK!" << std::endl;
      // }

      found=true;
      if( supports ) {
	vali = scope.size;
	while( vali-- )
	  supports[vri][vli][vali] = solution[vali];
      }
      break;
    } 
// else if(solver->parameters.verbosity) {
//       std::cout << "NO" << std::endl;
//     }
 
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


bool Mistral::Constraint::find_bound_support(const int vri, const int vli) 
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


      // std::cout << "REWRITE EQUALITY " ;
      // display(std::cout);
      // std::cout << std::endl;
      // //std::cout << solver->variables << std::endl;


  Mistral::PropagationOutcome wiped = propagate(); 

  // std::cout << "propag: " << scope[0] << " in " << scope[0].get_domain() 
  // 	    << " = " << scope[1] << " in " << scope[1].get_domain() << std::endl;


  if( active.size == 2 ) {
    if( scope[0].is_expression() && scope[1].is_expression() ) {


      // std::cout << "ACTUALLY REWRITE " ;
      // display(std::cout);
      // std::cout << std::endl;
      // // //std::cout << solver->variables << std::endl;


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

      // //std::cout << solver->variables << std::endl;
      // std::cout << "EQUALITY REWRITEN " ;
      // display(std::cout);
      // std::cout << std::endl;

    } 
  } else {
    relax();
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
  
  // precedence[0] = new ConstraintLess(scope, processing_time[0]);
  // precedence[1] = new ConstraintLess(processing_time[1]);
  // precedence[1]->add(scope[1]);
  // precedence[1]->add(scope[0]);
}

Mistral::ConstraintDisjunctive::ConstraintDisjunctive(std::vector< Variable >& scp, const int p0, const int p1) : Constraint(scp) { 
  processing_time[0] = p0; 
  processing_time[1] = p1;
  
  // precedence[0] = new ConstraintLess(scope, processing_time[0]);
  // precedence[1] = new ConstraintLess(processing_time[1]);
  // precedence[1]->add(scope[1]);
  // precedence[1]->add(scope[0]);
}

void Mistral::ConstraintDisjunctive::initialise() {

  //processing_time[0] = p0; 
  //processing_time[1] = p1;
  
  precedence[0] = new ConstraintLess(scope, processing_time[0]);
  precedence[1] = new ConstraintLess(processing_time[1]);
  precedence[1]->add(scope[1]);
  precedence[1]->add(scope[0]);

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
    //std::cout << "add " << precedence[1] << std::endl;
  } else {
    solver->add(precedence[0]);
    //std::cout << "add " << precedence[0] << std::endl;
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
      if( (spin + scope[2].get_min()) == 1 ) {
	solver->add(new ConstraintNotEqual(scope));
      } else {
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

  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateEqual::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;


  //std::cout << spin << std::endl;

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


void Mistral::PredicateFactor::initialise() {
  Constraint::initialise();
  for(int i=0; i<2; ++i)
    trigger_on(_RANGE_, i);
  set_idempotent(false);
}

Mistral::PropagationOutcome Mistral::PredicateFactor::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( IS_FAIL(scope[0].set_min( (factor>0 ? 
				 scope[1].get_min() : 
				 scope[1].get_max()) / factor )) ) 
    wiped = FAILURE(0); 

  else if( IS_FAIL(scope[0].set_max( (factor>0 ? 
				      scope[1].get_max() :
				      scope[1].get_min()) / factor )) ) 
    wiped = FAILURE(0); 

  else if( IS_FAIL(scope[1].set_min( (factor>0 ? 
				      scope[1].get_min() : 
				      scope[1].get_max()) * factor )) ) 
    wiped = FAILURE(1); 

  else if( IS_FAIL(scope[1].set_max( (factor>0 ? 
				      scope[1].get_max() :
				      scope[1].get_min()) * factor )) ) 
    wiped = FAILURE(1);
  
  return wiped;
}

std::ostream& Mistral::PredicateFactor::display(std::ostream& os) const {
  os << scope[1] << " == (" << scope[0] << " * " << factor << ")";
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
    trigger_on(_DOMAIN_, i);
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
    trigger_on(_VALUE_, i);
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

Mistral::PropagationOutcome Mistral::PredicateAdd::rewrite() {
  Mistral::PropagationOutcome wiped = propagate();

  VarArray tmp;
  if(active.size == 2) {
    int i=0;
    for(; i<2; ++i)
      if(scope[i].is_ground()) {
	relax();
	tmp.add(scope[1-i]);
	tmp.add(scope[2]);
	if(scope[i].get_min() == 0) {
	  solver->add(new ConstraintEqual(tmp));
	} else {
	  solver->add(new PredicateOffset(tmp, scope[i].get_min()));
	}
      }
  }
  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateAdd::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;

  // // x[0] + x[1] = x[2]

   // std::cout << std::endl
   // 	    << scope[0] << " in " << scope[0].get_domain() << " + "
   // 	    << scope[1] << " in " << scope[1].get_domain() << " = "
   // 	    << scope[2] << " in " << scope[2].get_domain() << std::endl;

  // int i, j;
  // while( IS_OK(wiped) && !events.empty() ) {
  //   while( IS_OK(wiped) && i-- ) {
  //     j = events.pop();

  //     std::cout << "react to event on " << scope[j] << " in " << scope[j].get_domain() << std::endl;

  //     if(j==2) { 
  // 	if(LB_CHANGED(event_type[2])) { // the lower bound of x[2] has changed
  // 	  if( IS_FAIL(scope[0].set_min(scope[2].get_min() - scope[1].get_max())) ) wiped = FAILURE(0);
  // 	  else if( IS_FAIL(scope[1].set_min(scope[2].get_min() - scope[0].get_max())) ) wiped = FAILURE(1);
  // 	} 
  // 	if(UB_CHANGED(event_type[2])) { // the upper bound of x[2] has changed
  // 	  if( IS_FAIL(scope[0].set_max(scope[2].get_max() - scope[1].get_min())) ) wiped = FAILURE(0);
  // 	  else if( IS_FAIL(scope[1].set_max(scope[2].get_max() - scope[0].get_min())) ) wiped = FAILURE(1);
  // 	}
  //     } else {
  // 	if(LB_CHANGED(event_type[j])) { // the lower bound of x[j] has changed
  // 	  if( IS_FAIL(scope[1-j].set_max(scope[2].get_max() - scope[j].get_min())) ) wiped = FAILURE(0);
  // 	  else if( IS_FAIL(scope[2].set_min(scope[0].get_min() + scope[1].get_min())) ) wiped = FAILURE(1);
  // 	}
  // 	if(UB_CHANGED(event_type[j])) { // the upper bound of x[j] has changed
  // 	  if( IS_FAIL(scope[1-j].set_min(scope[2].get_min() - scope[j].get_max())) ) wiped = FAILURE(0);
  // 	  else if( IS_FAIL(scope[2].set_max(scope[0].get_max() + scope[1].get_max())) ) wiped = FAILURE(1);
  // 	}
  //     }

  //     std::cout << scope[0] << " in " << scope[0].get_domain() << " + "
  // 		<< scope[1] << " in " << scope[1].get_domain() << " = "
  // 		<< scope[2] << " in " << scope[2].get_domain() << std::endl;

  //   }
  // }

  // std::cout << scope[0] << " in " << scope[0].get_domain() << " + "
  // 	    << scope[1] << " in " << scope[1].get_domain() << " = "
  // 	    << scope[2] << " in " << scope[2].get_domain() << std::endl << std::endl;

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
  os << scope[2] << " === (" << scope[0] << " + " << scope[1] << ")";
  return os;
}


void Mistral::PredicateMul::initialise() {
  Constraint::initialise();
  trigger_on(_RANGE_, 0);
  trigger_on(_RANGE_, 1);
  trigger_on(_RANGE_, 2);
  set_idempotent(true);
}

Mistral::PropagationOutcome Mistral::PredicateMul::rewrite() {
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
	  solver->add(new ConstraintEqual(tmp));
	} else if(scope[i].get_min() != 0) {
	  solver->add(new PredicateFactor(tmp, scope[i].get_min()));
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


Mistral::PropagationOutcome Mistral::PredicateMul::pruneZeros(const int changedIdx)
{

#ifdef _DEBUG_MUL
  std::cout << "check the zeros" << std::endl;
#endif

//   bool sprint = ( scope[0]->id == 5740 && scope[1]->id == 5741 );
    
//   if( sprint )
//     {
//       std::cout << "check the zeros" << std::endl;
//       print( std::cout );
//       std::cout << " " << changedIdx << std::endl;
//     }
  Mistral::PropagationOutcome wiped = CONSISTENT;

  if( !scope[changedIdx].contain(0) ) { // scope[changedIdx] cannot be 0
    if( changedIdx == 2 ) { // the product is not 0

 //      if( sprint )
// 	std::cout << "scope[0]->remove(0) && scope[1]->remove(0)" << std::endl;

      if( IS_FAIL(scope[0].remove(0)) ) wiped = FAILURE(0);
      else if( IS_FAIL(scope[1].remove(0)) ) wiped = FAILURE(1);
    } else { // one factor is not zero    
      if( !scope[1-changedIdx].contain(0) && IS_FAIL(scope[2].remove(0)) ) wiped = FAILURE(2);
      if( scope[2].is_ground() && scope[2].get_min() == 0 && 
	  IS_FAIL(scope[1-changedIdx].set_domain(0)) ) wiped = FAILURE(1-changedIdx);
    }
  } else if( scope[changedIdx].is_ground() ) { // scope[changedIdx] must be 0
    if( changedIdx == 2 ) { // the product is 0
      if( !scope[0].contain(0) && IS_FAIL(scope[1].set_domain(0)) ) wiped = FAILURE(1);
      else if( !scope[1].contain(0) && IS_FAIL(scope[0].set_domain(0)) ) wiped = FAILURE(0);
    } else { // one factor is 0
      if(IS_FAIL(scope[2].set_domain(0))) wiped = FAILURE(2);
    }
  }

#ifdef _DEBUG_MUL
  std::cout << scope[0].get_domain() 
	    << " * " << scope[1].get_domain() 
    	    << " = " << scope[2].get_domain() << std::endl;
#endif

  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateMul::pruneUnary(const int last) 
{
  Mistral::PropagationOutcome wiped = CONSISTENT;

#ifdef _DEBUG_MUL
  std::cout << "prune unary (" << scope[last] << " in " 
	    << scope[last].get_domain() << ")" << std::endl;
#endif

  if( last < 2 ) {
    int v = scope[1-last].get_min();
    int w = scope[2].get_min();
    if(v) {
      if(w%v || IS_FAIL(scope[last].set_domain( w/v ))) wiped = FAILURE(last);
    } else if(w) wiped = FAILURE(last);
  } else
    if( IS_FAIL(scope[last].set_domain( scope[0].get_min() * scope[1].get_min() )) ) 
      wiped = FAILURE(last);
  
#ifdef _DEBUG_MUL
  std::cout << scope[0].get_domain() 
	    << " * " << scope[1].get_domain() 
    	    << " = " << scope[2].get_domain() << std::endl;
#endif

  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateMul::pruneBinary(const int otherIdx, const int reviseIdx, const int v)
{
  Mistral::PropagationOutcome wiped = CONSISTENT;
//   std::cout << "revise " ;
//   scope[reviseIdx]->print(std::cout);
//   std::cout << " with respect to ";
//   scope[otherIdx]->print(std::cout);
//   std::cout << std::endl;

  // bool yox = false, xoy = false;

#ifdef _DEBUG_MUL
  std::cout << "prune binary (" << scope[reviseIdx] << " in " 
	    << scope[reviseIdx].get_domain() << " <- "
	    << scope[otherIdx] << " in " 
	    << scope[otherIdx].get_domain()
	    << ")" << std::endl;
#endif
  
  int w1, w2, lb, ub, r;
  int (*oper)(const int, const int, int&);
  if( otherIdx+reviseIdx != 1 ) {
     w1 = scope[otherIdx].get_min();
     w2 = scope[otherIdx].get_max();
    if( reviseIdx == 2 ) {
      oper = xtimey;
    } else {
      //yox = true;
      oper = yoverx;
    }
  } else {
    w1 = scope[otherIdx].get_max_neg();
    w2 = scope[otherIdx].get_min_pos();
    
    //xoy = true;
    oper = xovery;
  }



  // int w1 = scope[otherIdx].get_min(), w2 = scope[otherIdx].get_max(), lb, ub, r;
  if(v<0) {
    lb = w1;
    w1 = w2;
    w2 = lb;
  }


 //   if( (xoy && !w1) || (yox && !v) ) {

//     print( std::cout );
//     std::cout << std::endl;

#ifdef _DEBUG_MUL
  std::cout << "revise " << scope[reviseIdx] << " in "
	    << scope[reviseIdx].get_domain() << " with respect to "
	    << scope[otherIdx] << " in "
    	    << scope[otherIdx].get_domain() 
	    << " (" << v << ")" << std::endl;
#endif
//   }

#ifdef _DEBUG_MUL
  std::cout << "lb = " ;
  if(oper == yoverx) {
    std::cout << w1 << "/" << v << std::endl;
  } else if(oper == xovery) {
    std::cout << v << "/" << w1 << std::endl;
  } else {
    std::cout << v << "*" << w1 << std::endl;
  }
#endif

  lb = oper(v, w1, r);
//   if(r) {
//     std::cout << "**" << std::endl;
//     ++lb;
//   }
  switch( r ) {
  case -1: {
#ifdef _DEBUG_MUL
    std::cout << "lb is not exact" << std::endl;
    std::cout << w1 << "/" << v << " =/= " << lb << std::endl;
#endif    
    if(v>0) {
#ifdef _DEBUG_MUL
      std::cout << v << "*" << lb << " < " << w1 << std::endl;
#endif
      if(v*lb < w1) ++lb;
    } else {
#ifdef _DEBUG_MUL
      std::cout << v << "*" << lb << " > " << w1 << std::endl;
#endif
      if(v*lb > w1) ++lb;
    }
#ifdef _DEBUG_MUL
    std::cout << " => " << lb << std::endl;
#endif 
  } break; 
  case 1: {
#ifdef _DEBUG_MUL
    std::cout << "lb is not exact" << std::endl;
    std::cout << v << "/" << w1 << " =/= " << lb << std::endl;    
#endif   
    if(v>0) {
#ifdef _DEBUG_MUL
      std::cout << v << "/" << lb << " > " << w1 << "?" << std::endl;
#endif
      if((double)v/(double)lb < (double)w1) ++lb;
      //if(w1*lb > v) ++lb;
    } else {
#ifdef _DEBUG_MUL
      std::cout << v << "/" << lb << " > " << w1 << "?" << std::endl;
#endif
      if((double)v/(double)lb > (double)w1) ++lb;
      //if(w1*lb < v) ++lb;
    }
#ifdef _DEBUG_MUL
    std::cout << " => " << lb << std::endl;
#endif  
  } break;
#ifdef _DEBUG_MUL
  default : {
    std::cout << "lb is exact" << std::endl;
  }
#endif
  }

#ifdef _DEBUG_MUL
  std::cout << "ub = " ;
  if(oper == yoverx) {
    std::cout << w2 << "/" << v << std::endl;
  } else if(oper == xovery) {
    std::cout << v << "/" << w2 << std::endl;
  } else {
    std::cout << v << "*" << w2 << std::endl;
  }
#endif

  ub = oper(v, w2, r);
//   if(r) {
//     std::cout << "gg" << std::endl;
//     --ub;
//   }
  switch( r ) {
  case -1: {
#ifdef _DEBUG_MUL
    std::cout << "ub is not exact" << std::endl;
    std::cout << w2 << "/" << v << " =/= " << ub << std::endl;
#endif   
    if(v>0) {
#ifdef _DEBUG_MUL
      std::cout << v << "*" << ub << " > " << w2 << "?" << std::endl;
#endif      
      if(v*ub > w2) --ub;
    } else {
#ifdef _DEBUG_MUL
      std::cout << v << "*" << ub << " < " << w2 << "?" << std::endl;
#endif
      if(v*ub < w2) --ub;
    }
#ifdef _DEBUG_MUL
    std::cout << " => " << ub << std::endl;
#endif
  } break; 
  case 1: {
#ifdef _DEBUG_MUL
    std::cout << "ub is not exact" << std::endl;
    std::cout << v << "/" << w2 << " =/= " << ub << std::endl;    
#endif   
    if(v>0) {
#ifdef _DEBUG_MUL
      std::cout << v << "/" << ub << " > " << w2 << "?" << std::endl;
#endif     
      if(v/ub > w2) --ub;
    } else {
#ifdef _DEBUG_MUL
      std::cout << v << "/" << ub << " < " << w2 << "?" << std::endl;
#endif     
      if(v/ub < w2) --ub;
    }
#ifdef _DEBUG_MUL
    std::cout << " => " << ub << std::endl;
#endif
  } break;
#ifdef _DEBUG_MUL
  default : {
    std::cout << "ub is exact" << std::endl;
  }
#endif
  }

#ifdef _DEBUG_MUL
  std::cout << "set " << scope[reviseIdx] << " in " 
	    << scope[reviseIdx].get_domain() << " to [" 
	    << lb << ".." << ub << "]" << std::endl;
#endif
  
     event_type[reviseIdx] = scope[reviseIdx].set_min( lb );
      if(event_type[reviseIdx] != NO_EVENT) {
	if( IS_FAIL(event_type[reviseIdx]) ) wiped = FAILURE(reviseIdx);
	else if(!changes.contain(reviseIdx)) {
	  changes.add(reviseIdx);
	}
      }
  if(IS_OK(wiped)) {
	event_type[reviseIdx] = scope[reviseIdx].set_max( ub );
	if(event_type[reviseIdx] != NO_EVENT) {
	  if( IS_FAIL(event_type[reviseIdx]) ) wiped = FAILURE(reviseIdx);
	  else if(!changes.contain(reviseIdx)) {
	    changes.add(reviseIdx);
	  }
	}
      }

  // if( IS_FAIL(scope[reviseIdx].set_min(lb)) ) wiped = FAILURE(reviseIdx);
  // else if( IS_FAIL(scope[reviseIdx].set_max(ub)) ) wiped = FAILURE(reviseIdx);
  // else if( xdom.table && !scope[otherIdx]->is_range() && !scope[reviseIdx]->is_range() )
  //   {
  //     //std::cout << "there are holes, we achieve AC" << std::endl; 
  //     xdom.clear();
  //     DomainIterator *valit = scope[otherIdx]->begin();
  //     do xdom.insert( oper(v, *valit, r) );
  //     while( valit->next() ); 
  //     if( !scope[reviseIdx].set_domain( xdom ) ) return false;
  //   }


#ifdef _DEBUG_MUL
  std::cout << scope[0].get_domain() 
	    << " * " << scope[1].get_domain() 
    	    << " = " << scope[2].get_domain() << std::endl;
#endif

  std::cout << wiped << " " << IS_OK(wiped) << std::endl;

  return wiped;
}

Mistral::PropagationOutcome Mistral::PredicateMul::pruneTernary(const int reviseIdx)
{
  Mistral::PropagationOutcome wiped = CONSISTENT;  

  
  if( reviseIdx == 2 || !scope[2].contain(0) || !scope[1-reviseIdx].contain(0) )
    {

#ifdef _DEBUG_MUL
  std::cout << "prune ternary (" << scope[reviseIdx] << " in " 
	    << scope[reviseIdx].get_domain() << ") " << reviseIdx << std::endl;
#endif

      int (*oper)(const int, const int, int&);
      int x, y, bound[6], zero[3], lb, ub, v[4], i, r;
      
      for(i=0; i<3; ++i)
	{     
	  //if(reviseIdx==2 || i==reviseIdx || i==2) {
	  if(i!=1-reviseIdx) {
	    std::cout << "ext" << std::endl;
 
	    bound[2*i] = scope[i].get_min(); 
	    bound[2*i+1] = scope[i].get_max(); 
	  } else {

	    std::cout << "abs" << std::endl;
	    bound[2*i] = scope[i].get_min_pos(); 
	    bound[2*i+1] = scope[i].get_max_neg(); 
	  }

	  // zero[i] = 0;
	  // if(!bound[2*i]) { ++bound[2*i]; ++zero[i]; }
	  // if(!bound[2*i+1]) { --bound[2*i+1]; zero[i]+=2; }

	  // 	  if(scope[0]->id == 1079) {
#ifdef _DEBUG_MUL
	  std::cout << "{" << bound[2*i] << "," << bound[2*i+1] << "} " << zero[i] << std::endl;
#endif
	  // 	  }
	}

      if(reviseIdx == 2) {
	oper = xtimey;
	x = 0;
	y = 1;
      } else {
	oper = xovery;
	x = 2;
	y = 1-reviseIdx;
      }

      for(i=0; i<4; ++i) {
	v[i] = oper(bound[2*x+i/2], bound[2*y+i%2], r);

	// 	  if(scope[0]->id == 1079) {
#ifdef _DEBUG_MUL
	std::cout << bound[2*x+i/2] << (reviseIdx == 2 ? " * " : " / ") << bound[2*y+i%2] << " = " << v[i] << std::endl;
#endif
	// 	  }
      }

      ub = 0;
      lb = 0;

      for(i=1; i<4; ++i) 
	if( v[i] > v[ub] ) ub = i;
	else if( v[i] < v[lb] ) lb = i;

      lb = v[lb];
      ub = v[ub];
  
      // 	  if(scope[0]->id == 1079) {
#ifdef _DEBUG_MUL
      std::cout << "[" << lb << "," << ub << "]" << std::endl;
#endif
      // 	  }

      // if(zero[reviseIdx]) {
      // 	if(lb > 0) lb = 0;
      // 	if(ub < 0) ub = 0;
      // }
      // // 	    // WARNING CHANGE, MAY BE BUGGY!!
      // //       if( (zero[reviseIdx] & 1) && lb > 0 ) lb = 0;
      // //       if( (zero[reviseIdx] & 2) && ub < 0 ) ub = 0;

      // 	  if(scope[0]->id == 1079) {
#ifdef _DEBUG_MUL
      std::cout << "[" << lb << "," << ub << "]" << std::endl;
#endif
      // 	  }

      event_type[reviseIdx] = scope[reviseIdx].set_min( lb );
      if(event_type[reviseIdx] != NO_EVENT) {
	if( IS_FAIL(event_type[reviseIdx]) ) wiped = FAILURE(reviseIdx);
	else if(!changes.contain(reviseIdx)) {
	  changes.add(reviseIdx);
	}
      }
      
      if(IS_OK(wiped)) {
	event_type[reviseIdx] = scope[reviseIdx].set_max( ub );
	if(event_type[reviseIdx] != NO_EVENT) {
	  if( IS_FAIL(event_type[reviseIdx]) ) wiped = FAILURE(reviseIdx);
	  else if(!changes.contain(reviseIdx)) {
	    changes.add(reviseIdx);
	  }
	}
      }
      //if( IS_FAIL(scope[reviseIdx].set_max( ub )) ) wiped = FAILURE(reviseIdx);
 
#ifdef _DEBUG_MUL
  std::cout << scope[0].get_domain() 
	    << " * " << scope[1].get_domain() 
    	    << " = " << scope[2].get_domain() << std::endl;
#endif

   }

  return wiped;
}


// Mistral::PropagationOutcome Mistral::pmul1::propagate() {   


// }

void xdivy( const double x, const double y, int& lb, int& ub) {
  double dval = x/y;
  int cval = ceil(dval);
  int fval = floor(dval);
  if(cval < lb) lb = cval;
  if(fval > ub) ub = fval;
}

void xmuly( const int x, const int y, int& lb, int& ub) {
  int val = x*y;
  if(val < lb) lb = val;
  if(val > ub) ub = val;
}

void Mistral::PredicateMul::compute_division(Variable X, Variable Y, int& lb, int& ub) {
  lb = +INFTY;
  ub = -INFTY;

  if(X.contain(0)) {
    if(Y.contain(0)) {
      // case 1: both x and y contain 0, the interval is then [-INFTY, +INFTY]
      lb = -INFTY;
      ub = +INFTY;
    } else {
      // case 5: y does not contain 0
      // then x/y = min/max of all 4 combinations
      xdivy((double)(X.get_min()), (double)(Y.get_min()), lb, ub);
      xdivy((double)(X.get_min()), (double)(Y.get_max()), lb, ub);
      xdivy((double)(X.get_min()), (double)(Y.get_min_pos()), lb, ub);
      xdivy((double)(X.get_min()), (double)(Y.get_max_neg()), lb, ub);
      xdivy((double)(X.get_max()), (double)(Y.get_min()), lb, ub);
      xdivy((double)(X.get_max()), (double)(Y.get_max()), lb, ub);
      xdivy((double)(X.get_max()), (double)(Y.get_min_pos()), lb, ub);
      xdivy((double)(X.get_max()), (double)(Y.get_max_neg()), lb, ub);
    }
  } else {
    if(!Y.contain(0))  {
      // case 5: y does not contain 0
      // then x/y = min/max of all 4 combinations
      xdivy((double)(X.get_min()), (double)(Y.get_min()), lb, ub);
      xdivy((double)(X.get_min()), (double)(Y.get_max()), lb, ub);
      xdivy((double)(X.get_min()), (double)(Y.get_min_pos()), lb, ub);
      xdivy((double)(X.get_min()), (double)(Y.get_max_neg()), lb, ub);
      xdivy((double)(X.get_max()), (double)(Y.get_min()), lb, ub);
      xdivy((double)(X.get_max()), (double)(Y.get_max()), lb, ub);
      xdivy((double)(X.get_max()), (double)(Y.get_min_pos()), lb, ub);
      xdivy((double)(X.get_max()), (double)(Y.get_max_neg()), lb, ub);
    } else if(Y.get_size()>1) { // case 2: x does not contain 0 and y = 0, the interval is void
      if(Y.get_min() == 0) {
	// case 4: x does not contain 0 and either min(y) = 0 or max(y) = 0
	// then x/y = x/z where d(z) = d(y)-{0}
	xdivy((double)(X.get_min()), (double)(Y.get_min_pos()), lb, ub);
	xdivy((double)(X.get_min()), (double)(Y.get_max()), lb, ub);
	xdivy((double)(X.get_max()), (double)(Y.get_min_pos()), lb, ub);
	xdivy((double)(X.get_max()), (double)(Y.get_max()), lb, ub);	
      } else if(Y.get_max() == 0) {
	// case 4: x does not contain 0 and either min(y) = 0 or max(y) = 0
	// then x/y = x/z where d(z) = d(y)-{0}
	xdivy((double)(X.get_min()), (double)(Y.get_min()), lb, ub);
	xdivy((double)(X.get_min()), (double)(Y.get_max_neg()), lb, ub);
	xdivy((double)(X.get_max()), (double)(Y.get_min()), lb, ub);
	xdivy((double)(X.get_max()), (double)(Y.get_max_neg()), lb, ub);
      } else {
	// case 3: x does not contain 0 and min(y) < 0 < max(y), the interval is
	// [min(-max(x),-min(x),max(x),min(x)), max(-max(x),-min(x),max(x),min(x))]
  	xdivy((double)(X.get_min()), (double)(Y.get_min_pos()), lb, ub);
	xdivy((double)(X.get_min()), (double)(Y.get_max_neg()), lb, ub);
	xdivy((double)(X.get_max()), (double)(Y.get_min_pos()), lb, ub);
	xdivy((double)(X.get_max()), (double)(Y.get_max_neg()), lb, ub);
      }
    }
  }
}



Mistral::PropagationOutcome Mistral::PredicateMul::revise_division(const int X, const int Y, const int Z) {
  // revise the domain of Z = X/Y (because Z*Y = X)
  Mistral::PropagationOutcome wiped = CONSISTENT;


#ifdef _DEBUG_MUL
  std::cout << "revise bounds of " << scope[Z].get_domain() << " = " 
	    << scope[X].get_domain() << "/" << scope[Y].get_domain() 
	    // << " = [" << min_neg[Z] << ".." << max_neg[Z] << "|" 
	    // << (zero[Z] ? "{0}" : "_") << "|" 
	    // << min_pos[Z] << ".." << max_pos[Z] << "]" 
	    << std::endl; 
#endif

  
  //int lb_pos=1, ub_pos=+INFTY, lb_neg=-INFTY, ub_neg=-1, lb_aux, ub_aux;
  
  // we start with all bounds at their previous values, and update them if necessary
  // (in which case we set up the pruning flag)
  int lb_pos=min_pos[Z], ub_pos=max_pos[Z], lb_neg=min_neg[Z], ub_neg=max_neg[Z], 
    nlb1, nlb2, nub1, nub2;
    //lb_aux, ub_aux;
  bool pruning_flag = false, pzero = false;
  
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

  //	std::cout << pruning_flag << std::endl;

  if(lb_neg != ub_pos && (!zero[X] || !zero[Y])) { // if X and Y can both be 0, we cannot deduce anything
    if(IS_OK(wiped)) {
      if(max_pos[Z]>0) {
	// revise the positive part of Z's domain (if it has one)
	nlb1 = nlb2 = INFTY; //lb_neg;
	nub1 = nub2 = 0; //ub_neg;
	
	// it can either be the positive parts of X and Y:
	if(max_pos[X]>0 && max_pos[Y]>0) {
	  // compute the bounds
	  //std::cout << "\t   lb = " << min_pos[X] << "/" << max_pos[Y] << std::endl;
	  //std::cout << "\t   ub = " << max_pos[X] << "/" << min_pos[Y] << std::endl;
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
	
	//std::cout << "positive bounds: [" << nlb1 << ".." << nub1 << "]" << std::endl;

	if(lb_pos<nlb1) {
	  lb_pos = nlb1;
	  pruning_flag = true;
	}
	if(ub_pos>nub1) {
	  ub_pos = nub1;
	  pruning_flag = true;
	}
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

	//std::cout << "negative bounds: [" << nlb1 << ".." << nub1 << "]" << std::endl;

	if(lb_neg<nlb1) {
	  lb_neg = nlb1;
	  pruning_flag = true;
	}
	if(ub_neg>nub1) {
	  ub_neg = nub1;
	  pruning_flag = true;
	}
     } else if(pzero || !zero[Z])// if(lb_neg || ub_neg)
       {
       lb_neg = min_pos[Z];
       ub_neg = max_pos[Z];
       } else {
       lb_neg = ub_neg = 0;
     }
    }
  }

  if(pruning_flag) {
    
#ifdef _DEBUG_MUL
    std::cout << "set bounds to " // << scope[Z].get_domain() << " = " 
	      // << scope[X].get_domain() << "/" << scope[Y].get_domain() 
	      << "[" << lb_neg << ".." << ub_neg << "|" 
	      << (!pzero&&zero[Z] ? "{0}" : "_") << "|" 
	      << lb_pos << ".." << ub_pos << "]" << std::endl; 
#endif
    
    wiped = prune(lb_neg, ub_neg, lb_pos, ub_pos, pzero, Z);
  }

  return wiped;
}


Mistral::PropagationOutcome Mistral::PredicateMul::revise_multiplication(const int X, const int Y, const int Z) {
  // revise the domain of Z = X*Y 
  Mistral::PropagationOutcome wiped = CONSISTENT;


#ifdef _DEBUG_MUL
  std::cout << "revise bounds of " << scope[Z].get_domain() << " = " 
	    << scope[X].get_domain() << "*" << scope[Y].get_domain() 
	    // << " = [" << min_neg[Z] << ".." << max_neg[Z] << "|" 
	    // << (zero[Z] ? "{0}" : "_") << "|" 
	    // << min_pos[Z] << ".." << max_pos[Z] << "]" 
	    << std::endl; 
#endif

  
  //int lb_pos=1, ub_pos=+INFTY, lb_neg=-INFTY, ub_neg=-1, lb_aux, ub_aux;
  int lb_pos=min_pos[Z], ub_pos=max_pos[Z], lb_neg=min_neg[Z], ub_neg=max_neg[Z], 
    nlb1, nlb2, nub1, nub2;
  
  bool pruning_flag = false, pzero = false;
  

  // std::cout << min_neg[X] << " " << max_pos[X] << " || " 
  // 	    <<  min_neg[Y]  << " " << max_pos[Y] << std::endl;

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

  // std::cout << lb_neg << "\\" << ub_pos << std::endl;

  if(lb_neg != ub_pos) { 
    if(IS_OK(wiped)) {
      if(max_pos[Z]>0) {
	// revise the positive part of Z's domain (if it has one)
	//ub_pos = 0;
	//lb_pos = 1;
	//lb_aux = 1;
	nlb1 = nlb2 = INFTY; //lb_neg;
	nub1 = nub2 = 0; //ub_neg;

	// it can either be the positive parts of X and Y:
	if(max_pos[X]>0 && max_pos[Y]>0) {
	  // compute the bounds
	  // ub_pos = max_pos[X] * max_pos[Y];
	  // lb_pos = min_pos[X] * min_pos[Y];
	  nub1 = max_pos[X] * max_pos[Y];
	  nlb1 = min_pos[X] * min_pos[Y];
	}
	// or the negative parts of X and Y:
	if(min_neg[X]<0 && min_neg[Y]<0) {
	  // compute the bounds
	  // ub_aux = min_neg[X] * min_neg[Y];
	  // lb_aux = max_neg[X] * max_neg[Y];
	  nub2 = min_neg[X] * min_neg[Y];
	  nlb2 = max_neg[X] * max_neg[Y];
	}
	// if(lb_pos>lb_aux) lb_pos = lb_aux;
	// if(ub_pos<ub_aux) ub_pos = ub_aux;
	if(nlb1>nlb2) nlb1 = nlb2;
	if(nub1<nub2) nub1 = nub2;
	
	//std::cout << "positive bounds: [" << nlb1 << ".." << nub1 << "]" << std::endl;

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
       //ub_neg = -1;
       //lb_neg = 0;
       //ub_aux = -1;
       nlb1 = nlb2 = 0; //lb_pos;
       nub1 = nub2 = -INFTY; //ub_pos;
	
	// it can either be the negitive part of X and the positive part of Y:
	if(min_neg[X]<0 && max_pos[Y]>0) {
	  // compute the bounds
	  // ub_neg = max_neg[X] * min_pos[Y];
	  // lb_neg = min_neg[X] * max_pos[Y];
	  nub1 = max_neg[X] * min_pos[Y];
	  nlb1 = min_neg[X] * max_pos[Y];

	  // std::cout << "\t  1 ub = " << max_neg[X] << "*" << min_pos[Y] << std::endl;
	  // std::cout << "\t  1 lb = " << min_neg[X] << "*" << max_pos[Y] << std::endl;

	}
	// or the negitive part of Y and the positive part of X:
	if(max_pos[X]>0 && min_neg[Y]<0) {
	  // compute the bounds
	  // ub_aux = max_neg[Y] * min_pos[X];
	  // lb_aux = min_neg[Y] * max_pos[X];
	  nub2 = max_neg[Y] * min_pos[X];
	  nlb2 = min_neg[Y] * max_pos[X];

	  // std::cout << "\t  2 ub = " << max_neg[Y] << "*" << min_pos[X] << std::endl;
	  // std::cout << "\t  2 lb = " << min_neg[Y] << "*" << max_pos[X] << std::endl;
	}
	// if(lb_neg>lb_aux) lb_neg = lb_aux;
	// if(ub_neg<ub_aux) ub_neg = ub_aux;

	if(nlb1>nlb2) nlb1 = nlb2;
	if(nub1<nub2) nub1 = nub2;


	  // std::cout << "change:" << std::endl;
	  // std::cout << "\t lbn = " << nlb1 << std::endl;
	  // std::cout << "\t ubn = " << nub1 << std::endl;


	//std::cout << "negative bounds: [" << nlb1 << ".." << nub1 << "]" << std::endl;
	
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
  std::cout << "set bounds to " // << scope[Z].get_domain() << " = " 
	    // << scope[X].get_domain() << "*" << scope[Y].get_domain() 
	    << "[" << lb_neg << ".." << ub_neg << "|" 
	    << (zero[Z] ? "{0}" : "_") << "|" 
	    << lb_pos << ".." << ub_pos << "]" << std::endl; 
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

// #ifdef _DEBUG_MUL
//   std::cout
// 	      << " was in [" << min_neg[Z] << ".." << max_neg[Z] << "|" 
// 	      << (zero[Z] ? "{0}" : "_") << "|" 
// 	      << min_pos[Z] << ".." << max_pos[Z] << "]" << std::endl; 

//   std::cout
// 	      << " set to [" << lb_neg << ".." << ub_neg << "|" 
// 	      << (zero[Z]&&!pzero ? "{0}" : "_") << "|" 
// 	      << lb_pos << ".." << ub_pos << "]" << std::endl; 

// 	  std::cout << scope[0].get_domain() 
// 		    << " * " << scope[1].get_domain() 
// 		    << " = " << scope[2].get_domain() << std::endl;
// #endif 

// std::cout << changes << std::endl;
  Event evt;
  Mistral::PropagationOutcome wiped = CONSISTENT;
  
  // int lb = lb_neg;
  // if(lb>=0) lb = lb_pos;
  // int ub = ub_pos;
  // if(ub>=0) lb = lb_pos;

  if(ub_pos < lb_neg) wiped = FAILURE(Z);
  else {
    if(lb_neg>min_neg[Z]) {

      // std::cout // << lb_neg << ">" << min_neg[Z] 
      // 		<< " -> update lb" << std::endl;

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

      // std::cout // << ub_pos << ">" << max_pos[Z] 
      // 		<< " -> update ub" << std::endl;

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
      
      // std::cout << lb_pos << ">=" << min_neg[Z] 
      // 		<< " or " << ub_neg << "<=" << max_pos[Z]
      // 		<< " -> may update inbounds" << std::endl;
      
      if(lb_pos-1>ub_neg && (pzero || (!zero[Z] && (min_pos[Z]<lb_pos || max_neg[Z]>ub_neg)))) {
	  
	// std::cout // << lb_pos << ">" << min_pos[Z] 
	// 	  // << " or " << ub_neg << "<" << max_neg[Z]
	// 	  << " -> update inbounds" << std::endl;

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


	  // std::cout //<< lb_pos << ">" << min_pos[Z] 
	  // 	  << " -> update negative ub" << std::endl;

	  evt = scope[Z].remove_interval(1, lb_pos-1);
	  if( IS_FAIL(evt) ) wiped = FAILURE(Z);
	  else {
	    //min_pos[Z] = lb_pos;
	    //max_neg[Z] = ub_neg;
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

	  // std::cout 
	  //   //<< ub_neg << "<" << max_neg[Z]
	  // 	  << " -> update positive lb" << std::endl;

	  evt = scope[Z].remove_interval(ub_neg+1, -1);
	  if( IS_FAIL(evt) ) wiped = FAILURE(Z);
	  else {
	    //min_pos[Z] = lb_pos;
	    //max_neg[Z] = ub_neg;
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




  //std::cout << changes << std::endl;
#ifdef _DEBUG_MUL
	if(IS_OK(wiped)) {
	  std::cout << scope[0].get_domain() 
		    << " * " << scope[1].get_domain() 
		    << " = " << scope[2].get_domain() << std::endl;
	} else std::cout << "FAIL!" << std::endl ;

	std::cout
	  << " now in [" << min_neg[Z] << ".." << max_neg[Z] << "|" 
	  << (zero[Z] ? "{0}" : "_") << "|" 
	  << min_pos[Z] << ".." << max_pos[Z] << "]" << std::endl;

  std::cout << std::endl; 
#endif


  return wiped;
}



void Mistral::PredicateMul::compute_multiplication(Variable X, Variable Y, int& lb, int& ub) {
  lb = +INFTY;
  ub = -INFTY;

  xmuly((X.get_min()), (Y.get_min()), lb, ub);
  xmuly((X.get_min()), (Y.get_max()), lb, ub);
  xmuly((X.get_max()), (Y.get_min()), lb, ub);
  xmuly((X.get_max()), (Y.get_max()), lb, ub);
}


// void Mistral::PredicateMul::refine_lbound(const int x, const int y, const int z, 
// 					  int (*oper)(const int, const int, 
// 						      int&, int&)) {
  
//   int lbound = scope[x].get_min();
//   int rest;
//   int lb, ub;

//   do {
//     oper(lbound, scope[y].get_min())
//   }

// }


void Mistral::PredicateMul::refine_bounds(const int evt_idx) {
  double dval, db1, db2, da1, da2, daux;
  int ival, ib1, ib2, ia1, ia2, iaux;
    // 
    if(evt_idx == 2) {
      // check the bounds of scope[2] 
      for(int i=0; i<2; ++i) {
	da1 = (double)(scope[i].get_min());	
	da2 = (double)(scope[i].get_max());
	if(LB_CHANGED(event_type[evt_idx])) {
#ifdef _DEBUG_MUL
	  std::cout << " check lower bound of " << scope[2] << " in " << scope[2].get_domain() 
		    << " with respect to " << scope[1-i] << " in " << scope[1-i].get_domain() << std::endl;
#endif
	  // check that [min(2)/max(i), min(2)/min(i)] intersects [min(1-i), max(1-i)]
	  dval = ival = scope[2].get_min();
	  do {
	    db1 = (da1 != 0 ? dval/da1 : -INFTY);
	    db2 = (da2 != 0 ? dval/da2 : +INFTY);
	    if(db1>db2) {
	      daux = db1;
	      db1 = db2;
	      db2 = daux;
	    }
	    ib1 = ceil(db1);
	    ib2 = floor(db2);
#ifdef _DEBUG_MUL
	    std::cout << "   : " << dval << " -> [" << ib1 << ".." << ib2 << "] inter "
		      << scope[1-i].get_domain() << std::endl;
#endif
	    ++dval;
	  } while( ib2<ib1 || ib2<scope[1-i].get_min() || ib1>scope[1-i].get_max() );
	  ival = dval-1;
	  if(ival > scope[2].get_min()) {
	    exit(1);
	  }
	}

	if(UB_CHANGED(event_type[evt_idx])) {
#ifdef _DEBUG_MUL
	  std::cout << " check upper bound of " << scope[2] << " in " << scope[2].get_domain() 
		    << " with respect to " << scope[1-i] << " in " << scope[1-i].get_domain() << std::endl;
#endif
	  // check that [max(2)/max(i), max(2)/min(i)] intersects [min(1-i), max(1-i)]
	  dval = ival = scope[2].get_max();
	  do {
	    db1 = (da1 != 0 ? dval/da1 : -INFTY);
	    db2 = (da2 != 0 ? dval/da2 : +INFTY);
	    if(db1>db2) {
	      daux = db1;
	      db1 = db2;
	      db2 = daux;
	    }
	    ib1 = ceil(db1);
	    ib2 = floor(db2);
#ifdef _DEBUG_MUL
	    std::cout << "   : " << dval << " -> [" << ib1 << ".." << ib2 << "] inter "
		      << scope[1-i].get_domain() << std::endl;
#endif
	    --dval;
	  } while( ib2<ib1 || ib2<scope[1-i].get_min() || ib1>scope[1-i].get_max() );
	  ival = dval+1;
	  if(ival < scope[2].get_max()) {
	    exit(1);
	  }
	}
      }
    } else {
      // check the bounds of scope[evt_idx] with respect to scope[2]
      ia1 = scope[1-evt_idx].get_min();
      ia2 = scope[1-evt_idx].get_max();
      
      if(LB_CHANGED(event_type[evt_idx])) {
#ifdef _DEBUG_MUL
	std::cout << " check lower bound of " << scope[evt_idx] << " in " << scope[evt_idx].get_domain() 
		  << " with respect to " << scope[2] << " in " << scope[2].get_domain() << std::endl;
#endif
	
	ival = scope[evt_idx].get_min();
	
	do {
	  ib1 = ival*ia1;
	  ib2 = ival*ia2;
	  if(ib1>ib2) {
	    iaux = ib1;
	    ib1 = ib2;
	    ib2 = iaux;
	  }
#ifdef _DEBUG_MUL
	  std::cout << "   : " << ival << " -> [" << ib1 << ".." << ib2 << "] inter "
		    << scope[2].get_domain() << std::endl;
#endif
	  ++ival;
	} while( ib2<scope[2].get_min() || ib1>scope[2].get_max() );
	
 	--ival;
	if(ival > scope[evt_idx].get_min()) {
	  exit(1);
	}
      }


      if(UB_CHANGED(event_type[evt_idx])) {
#ifdef _DEBUG_MUL
	std::cout << " check upper bound of " << scope[evt_idx] << " in " << scope[evt_idx].get_domain() 
		  << " with respect to " << scope[2] << " in " << scope[2].get_domain() << std::endl;
#endif
	
	ival = scope[evt_idx].get_max();
	
	do {
	  ib1 = ival*ia1;
	  ib2 = ival*ia2;
	  if(ib1>ib2) {
	    iaux = ib1;
	    ib1 = ib2;
	    ib2 = iaux;
	  }
#ifdef _DEBUG_MUL
	  std::cout << "   : " << ival << " -> [" << ib1 << ".." << ib2 << "] inter "
		    << scope[2].get_domain() << std::endl;
#endif
	  --ival;
	} while( ib2<scope[2].get_min() || ib1>scope[2].get_max() );
	
 	++ival;
	if(ival < scope[evt_idx].get_max()) {
	  exit(1);
	}
      }



      // check the bounds of scope[evt_idx] with respect to scope[1-evt_idx]
      da1 = (double)(scope[2].get_min());
      da2 = (double)(scope[2].get_max());
      
      if(LB_CHANGED(event_type[evt_idx])) {
#ifdef _DEBUG_MUL
	std::cout << " check lower bound of " << scope[evt_idx] << " in " << scope[evt_idx].get_domain() 
		  << " with respect to " << scope[1-evt_idx] << " in " << scope[1-evt_idx].get_domain() << std::endl;
#endif
	
	dval = ival = scope[evt_idx].get_min();
	do {
	  db1 = (dval != 0 ? da1/dval : -INFTY);
	  db2 = (dval != 0 ? da2/dval : +INFTY);
	  if(db1>db2) {
	    daux = db1;
	    db1 = db2;
	    db2 = daux;
	  }
	  ib1 = ceil(db1);
	  ib2 = floor(db2);
#ifdef _DEBUG_MUL
	  std::cout << "   : " << dval << " -> [" << ib1 << ".." << ib2 << "] inter "
		    << scope[1-evt_idx].get_domain() << std::endl;
#endif
	  ++dval;
	} while( ib2<ib1 || ib2<scope[1-evt_idx].get_min() || ib1>scope[1-evt_idx].get_max() );
	
	ival = dval-1;
	if(ival > scope[evt_idx].get_min()) {
	  exit(1);
	}
      }


      if(UB_CHANGED(event_type[evt_idx])) {
#ifdef _DEBUG_MUL
	std::cout << " check upper bound of " << scope[evt_idx] << " in " << scope[evt_idx].get_domain() 
		  << " with respect to " << scope[1-evt_idx] << " in " << scope[1-evt_idx].get_domain() << std::endl;
#endif
	
	dval = ival = scope[evt_idx].get_max();
	do {
	  db1 = (dval != 0 ? da1/dval : -INFTY);
	  db2 = (dval != 0 ? da2/dval : +INFTY);
	  if(db1>db2) {
	    daux = db1;
	    db1 = db2;
	    db2 = daux;
	  }
	  ib1 = ceil(db1);
	  ib2 = floor(db2);
#ifdef _DEBUG_MUL
	  std::cout << "   : " << dval << " -> [" << ib1 << ".." << ib2 << "] inter "
		    << scope[1-evt_idx].get_domain() << std::endl;
#endif
	  --dval;
	} while( ib2<ib1 || ib2<scope[1-evt_idx].get_min() || ib1>scope[1-evt_idx].get_max() );
	
	ival = dval+1;
	if(ival < scope[evt_idx].get_max()) {
	  exit(1);
	}
      }
      
    }
}


Mistral::PropagationOutcome Mistral::PredicateMul::propagate() {      
  Mistral::PropagationOutcome wiped = CONSISTENT;
  
// #ifdef _DEBUG_MUL
//   std::cout << std::endl << std::endl << scope[0].get_domain() 
// 	    << " * " << scope[1].get_domain() 
//     	    << " = " << scope[2].get_domain() << std::endl;
// #endif

//   int i, j, k, v, changedIdx;
  
//   while(!changes.empty()) {

//     changedIdx = changes.pop();

// #ifdef _DEBUG_MUL
//            std::cout << "react to " << scope[changedIdx] << " in " 
// 		     << scope[changedIdx].get_domain() << std::endl;
// #endif

//     i = (changedIdx+1)%3;
//     j = (changedIdx+2)%3;
//     k=3;
    
//     do wiped = pruneZeros(--k);
//     while( IS_OK(wiped) && k );
    
    
//     if( IS_OK(wiped) ) {
//       /// first, we check the particular cases:
//       /*************************************************
//        * Case 1: two variables are ground 
//        * Case 2: the calling variable is ground
//        * Case 3: one variable is ground
//        *************************************************/
      
//       if( scope[changedIdx].is_ground() )
// 	{
// 	  if( scope[i].is_ground() ) {
// 	    //std::cout << "FC" << std::endl;
// 	    wiped = pruneUnary(j);
// 	  } else if( scope[j].is_ground() ) {
// 	    //std::cout << "FC" << std::endl;
// 	    wiped = pruneUnary(i);
// 	  } else {
// 	    v = scope[changedIdx].get_min();
// 	    //std::cout << "AC-1" << std::endl;
// 	    wiped = pruneBinary(i, j, v);
// 	    if(IS_OK(wiped)) wiped = pruneBinary(j, i, v);
// 	  }
// 	}
//       else {
// 	switch( scope[i].is_ground() + 2*scope[j].is_ground() ) {
// 	case 0: {
// 	  wiped = pruneTernary(i);
// 	  if(IS_OK(wiped)) wiped = pruneTernary(j);
// 	} break;
// 	case 1: {
// 	  v = scope[i].get_min();
// 	  wiped = pruneBinary(changedIdx, j, v);
// 	} break;
// 	case 2: {
// 	  v = scope[j].get_min();
// 	  wiped = pruneBinary(changedIdx, i, v);
// 	} break;
// 	  //default: {
// 	  //}
// 	}
//       }
//     }
//   }






// #ifdef _DEBUG_MUL
//   std::cout << std::endl << std::endl << scope[0] << " in " << scope[0].get_domain() 
// 	    << " * " << scope[1] << " in " << scope[1].get_domain() 
//     	    << " = " << scope[2] << " in " << scope[2].get_domain() << std::endl;
// #endif

//   int i, j, k, lb, ub, evt_idx, rev_idx, aux_idx;

//   int (*oper)(const int, const int, int&);
  
//   int lbs[3], ubs[3], min_ps[3], max_ns[3], zero[3], sign[3], vals[4], rest[4], *x[2], *y[2];
  
//   for(i=0; i<3; ++i) {
//     lbs[i] = scope[i].get_min();
//     ubs[i] = scope[i].get_max();
//     min_ps[i] = scope[i].get_min_pos(); //(scope[i].contain(0) ? 0 : scope[i].get_min_pos());
//     max_ns[i] = scope[i].get_max_neg(); //(scope[i].contain(0) ? 0 : scope[i].get_max_neg());
//     zero[i] = scope[i].contain(0);
//     sign[i] = (lbs[i]<0 + 2*(ubs[i]>0));
//   }

  


//   while(IS_OK(wiped) && !changes.empty()) {

//     evt_idx = changes.pop();

// #ifdef _DEBUG_MUL
//     std::cout << "react to " << scope[evt_idx] << " in " 
// 	      << scope[evt_idx].get_domain() << std::endl;
// #endif
    
//     if(evt_idx == 2) {
//       for(i=0; IS_OK(wiped) && i<2; ++i) { // once for scope[0], once for scope[1]
// 	// x0 * x1 = x2
       
// 	// we update x0 and x1 :  xi = x2/x1-i
// 	rev_idx = i;
// 	aux_idx = 1-i;
// 	oper = xovery; // evt / aux


// #ifdef _DEBUG_MUL
// 	std::cout << "check bounds of " << scope[rev_idx] << " = " 
// 		  << scope[evt_idx] << "/" << scope[aux_idx] ;
// #endif


	
	

// 	// the lower value is obtained using x[2]'s min absolute value 
// 	// divided by x[i]'s max absolute value 
// 	x[0] = min_ps;
// 	x[1] = max_ns;
// 	y[0] = ubs;
// 	y[1] = lbs;

	
// 	for(j=0; j<2; ++j)
// 	  for(k=0; k<2; ++k) 
// 	    vals[2*j+k] = oper(x[j][evt_idx],y[k][aux_idx],rest[2*j+k]);

// 	lb = vals[0];
// 	for(j=1; j<4; ++j)
// 	  if(lb>vals[j]) lb = vals[j];


// 	// the higher value is obtained using x[2]'s max absolute value 
// 	// divided by x[i]'s min absolute value 
// 	x[0] = ubs;
// 	x[1] = lbs;
// 	y[0] = min_ps;
// 	y[1] = max_ns;

// 	for(j=0; j<2; ++j)
// 	  for(k=0; k<2; ++k) 
// 	    vals[2*j+k] = oper(x[j][evt_idx],y[k][aux_idx],rest[2*j+k]);

// 	ub = vals[0];
// 	for(j=1; j<4; ++j)
// 	  if(ub<vals[j]) ub = vals[j];

// #ifdef _DEBUG_MUL
// 	std::cout << " = [" << lb << ".." << ub << "]" << std::endl; 
// #endif

// 	event_type[rev_idx] = scope[rev_idx].set_min( lb );
// 	if(event_type[rev_idx] != NO_EVENT) {
// 	  if( IS_FAIL(event_type[rev_idx]) ) wiped = FAILURE(rev_idx);
// 	  else if(!changes.contain(rev_idx)) {
// 	    changes.add(rev_idx);
// 	  }
// 	}
	
// 	if(IS_OK(wiped)) {
// 	  event_type[rev_idx] = scope[rev_idx].set_max( ub );
// 	  if(event_type[rev_idx] != NO_EVENT) {
// 	    if( IS_FAIL(event_type[rev_idx]) ) wiped = FAILURE(rev_idx);
// 	    else if(!changes.contain(rev_idx)) {
// 	      changes.add(rev_idx);
// 	    }
// 	  }
// 	}


// #ifdef _DEBUG_MUL
// 	if(IS_OK(wiped)) {
// 	  std::cout << scope[0].get_domain() 
// 		    << " * " << scope[1].get_domain() 
// 		    << " = " << scope[2].get_domain() << std::endl;
// 	} else std::cout << "FAIL!" << std::endl;
// #endif

//       }
//     } else {
//       // update x[2]
      
//       // we update x0 and x1 :  xi = x2/x1-i
//       rev_idx = 2;
//       aux_idx = 1-evt_idx;
//       oper = xtimey;  
      
      
// #ifdef _DEBUG_MUL
//       std::cout << "check bounds of " << scope[rev_idx] << " = " 
// 		<< scope[evt_idx] << "*" << scope[aux_idx] ;
// #endif
      
//       // the lower value is obtained using x[2]'s min absolute value 
//       // divided by x[i]'s max absolute value 
//       x[0] = lbs;
//       x[1] = ubs;
//       y[0] = ubs;
//       y[1] = lbs;
      
//       for(j=0; j<2; ++j)
// 	for(k=0; k<2; ++k) 
// 	  vals[2*j+k] = oper(x[j][evt_idx],y[k][aux_idx],rest[2*j+k]);
      
//       lb = vals[0];
//       ub = vals[0];
//       for(j=1; j<4; ++j) {
// 	if(lb>vals[j]) lb = vals[j];
// 	if(ub<vals[j]) ub = vals[j];
//       }
      
// #ifdef _DEBUG_MUL
//       std::cout << " = [" << lb << ".." << ub << "]" << std::endl; 
// #endif

      
//       event_type[rev_idx] = scope[rev_idx].set_min( lb );
//       if(event_type[rev_idx] != NO_EVENT) {
// 	if( IS_FAIL(event_type[rev_idx]) ) wiped = FAILURE(rev_idx);
// 	else if(!changes.contain(rev_idx)) {
// 	  changes.add(rev_idx);
// 	}
//       }
      
//       if(IS_OK(wiped)) {
// 	event_type[rev_idx] = scope[rev_idx].set_max( ub );
// 	if(event_type[rev_idx] != NO_EVENT) {
// 	  if( IS_FAIL(event_type[rev_idx]) ) wiped = FAILURE(rev_idx);
// 	  else if(!changes.contain(rev_idx)) {
// 	    changes.add(rev_idx);
// 	  }
// 	}
//       }

// #ifdef _DEBUG_MUL
//       if(IS_OK(wiped)) {
// 	std::cout << scope[0].get_domain() 
// 		  << " * " << scope[1].get_domain() 
// 		  << " = " << scope[2].get_domain() << std::endl;
//       } else std::cout << "FAIL!" << std::endl;
// #endif

//       if(IS_OK(wiped)) {
	
// 	// we update x[1-evt_idx]
// 	rev_idx = 1-evt_idx;
// 	aux_idx = 2;
// 	oper = yoverx; // aux / evt

// #ifdef _DEBUG_MUL
// 	std::cout << "check bounds of " << scope[rev_idx] << " = " 
// 		  << scope[aux_idx] << "/" << scope[evt_idx] ;
// #endif
	  
// 	// the lower value is obtained using x[2]'s min absolute value 
// 	// divided by x[i]'s max absolute value 
// 	y[0] = min_ps;
// 	y[1] = max_ns;
// 	x[0] = ubs;
// 	x[1] = lbs;
	
// 	for(j=0; j<2; ++j)
// 	  for(k=0; k<2; ++k) 
// 	    vals[2*j+k] = oper(x[j][evt_idx],y[k][aux_idx],rest[2*j+k]);
	
// 	lb = vals[0];
// 	for(j=1; j<4; ++j)
// 	  if(lb>vals[j]) lb = vals[j];
	
	
// 	// the higher value is obtained using x[2]'s max absolute value 
// 	// divided by x[i]'s min absolute value 
// 	y[0] = ubs;
// 	y[1] = lbs;
// 	x[0] = min_ps;
// 	x[1] = max_ns;
	
// 	for(j=0; j<2; ++j)
// 	  for(k=0; k<2; ++k) 
// 	    vals[2*j+k] = oper(x[j][evt_idx],y[k][aux_idx],rest[2*j+k]);
	
// 	ub = vals[0];
// 	for(j=1; j<4; ++j)
// 	  if(ub<vals[j]) ub = vals[j];
	
// #ifdef _DEBUG_MUL
// 	std::cout << " = [" << lb << ".." << ub << "]" << std::endl; 
// #endif
	
	
// 	event_type[rev_idx] = scope[rev_idx].set_min( lb );
// 	if(event_type[rev_idx] != NO_EVENT) {
// 	  if( IS_FAIL(event_type[rev_idx]) ) wiped = FAILURE(rev_idx);
// 	  else if(!changes.contain(rev_idx)) {
// 	    changes.add(rev_idx);
// 	  }
// 	}
	
// 	if(IS_OK(wiped)) {
// 	  event_type[rev_idx] = scope[rev_idx].set_max( ub );
// 	  if(event_type[rev_idx] != NO_EVENT) {
// 	    if( IS_FAIL(event_type[rev_idx]) ) wiped = FAILURE(rev_idx);
// 	    else if(!changes.contain(rev_idx)) {
// 	      changes.add(rev_idx);
// 	    }
// 	  }
// 	}
	
// #ifdef _DEBUG_MUL
// 	if(IS_OK(wiped)) {
// 	  std::cout << scope[0].get_domain() 
// 		    << " * " << scope[1].get_domain() 
// 		    << " = " << scope[2].get_domain() << std::endl;
// 	} else std::cout << "FAIL!" << std::endl;
// #endif

//       }
//     }
//   }


#ifdef _DEBUG_MUL
  std::cout << std::endl << std::endl << "propagate " // << this 
	    << std::endl;
#endif

  for(int i=0; i<3; ++i) {
    max_pos[i] = scope[i].get_max();
    min_neg[i] = scope[i].get_min();

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

  
#ifdef _DEBUG_MUL
  std::cout << scope[0].get_var() << " in " << scope[0].get_domain() 
	    << " * " << scope[1].get_var() << " in " << scope[1].get_domain() 
    	    << " = " << scope[2].get_var() << " in " << scope[2].get_domain() << std::endl;
#endif

  //int i, j, lb, ub, evt_idx, rev_idx, aux_idx;
  int evt_idx;
  //Event evt;

  while(IS_OK(wiped) && !changes.empty()) {

    evt_idx = changes.pop();

#ifdef _DEBUG_MUL
    std::cout << "react to " << scope[evt_idx].get_var() << " in " 
	      << scope[evt_idx].get_domain() 
	      << (LB_CHANGED(event_type[evt_idx]) ? " (change on LB) " : "")
	      << (UB_CHANGED(event_type[evt_idx]) ? " (change on UB) " : "")
	      << std::endl;
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
  
    




// #ifdef _DEBUG_MUL
//   std::cout << std::endl << std::endl << scope[0] << " in " << scope[0].get_domain() 
// 	    << " * " << scope[1] << " in " << scope[1].get_domain() 
//     	    << " = " << scope[2] << " in " << scope[2].get_domain() << std::endl;
// #endif

//   int i, j, lb, ub, evt_idx, rev_idx, aux_idx;
//   Event evt;

//   while(IS_OK(wiped) && !changes.empty()) {

//     evt_idx = changes.pop();

// #ifdef _DEBUG_MUL
//     std::cout << "react to " << scope[evt_idx] << " in " 
// 	      << scope[evt_idx].get_domain() 
// 	      << (LB_CHANGED(event_type[evt_idx]) ? " (change on LB) " : "")
// 	      << (UB_CHANGED(event_type[evt_idx]) ? " (change on UB) " : "")
// 	      << std::endl;
// #endif

//     //refine_bounds(evt_idx);

    
//     if(evt_idx == 2) {
//       for(i=0; IS_OK(wiped) && i<2; ++i) { // once for scope[0], once for scope[1]
// 	// x0 * x1 = x2
       
// 	// we update x0 and x1 :  xi = x2/x1-i
// 	rev_idx = i;
// 	aux_idx = 1-i;
// 	compute_division(scope[evt_idx], scope[aux_idx], lb, ub);

// #ifdef _DEBUG_MUL
// 	std::cout << "check bounds of " << scope[rev_idx] << " = " 
// 		  << scope[evt_idx] << "/" << scope[aux_idx] 
// 		  << " = [" << lb << ".." << ub << "]" << std::endl; 
// #endif

// 	evt = scope[rev_idx].set_min( lb );
// 	if( IS_FAIL(evt) ) wiped = FAILURE(rev_idx);
// 	else if(evt != NO_EVENT) {
// 	  if(changes.contain(rev_idx)) {
// 	    event_type[rev_idx] |= evt;
// 	  } else {
// 	    event_type[rev_idx] = evt;
// 	    changes.add(rev_idx);
// 	  }
// 	}
	
// 	if(IS_OK(wiped)) {
// 	  evt = scope[rev_idx].set_max( ub );
// 	  if( IS_FAIL(evt) ) wiped = FAILURE(rev_idx);
// 	  else if(evt != NO_EVENT) {
// 	    if(changes.contain(rev_idx)) {
// 	      event_type[rev_idx] |= evt;
// 	    } else {
// 	      event_type[rev_idx] = evt;
// 	      changes.add(rev_idx);
// 	    }
// 	  }
// 	}

// #ifdef _DEBUG_MUL
// 	if(IS_OK(wiped)) {
// 	  std::cout << scope[0].get_domain() 
// 		    << " * " << scope[1].get_domain() 
// 		    << " = " << scope[2].get_domain() << std::endl;
// 	} else std::cout << "FAIL!" << std::endl;
// #endif

//       }
//     } else {
//       // update x[2]
      
//       // we update x0 and x1 :  xi = x2/x1-i
//       rev_idx = 2;
//       aux_idx = 1-evt_idx;
//       compute_multiplication(scope[evt_idx], scope[aux_idx], lb, ub);
      
      
// #ifdef _DEBUG_MUL
//       std::cout << "check bounds of " << scope[rev_idx] << " = " 
// 		<< scope[evt_idx] << "*" << scope[aux_idx] 
// 		<< " = [" << lb << ".." << ub << "]" << std::endl; 
// #endif
      
// 	evt = scope[rev_idx].set_min( lb );
// 	if( IS_FAIL(evt) ) wiped = FAILURE(rev_idx);
// 	else if(evt != NO_EVENT) {
// 	  if(changes.contain(rev_idx)) {
// 	    event_type[rev_idx] |= evt;
// 	  } else {
// 	    event_type[rev_idx] = evt;
// 	    changes.add(rev_idx);
// 	  }
// 	}
	
// 	if(IS_OK(wiped)) {
// 	  evt = scope[rev_idx].set_max( ub );
// 	  if( IS_FAIL(evt) ) wiped = FAILURE(rev_idx);
// 	  else if(evt != NO_EVENT) {
// 	    if(changes.contain(rev_idx)) {
// 	      event_type[rev_idx] |= evt;
// 	    } else {
// 	      event_type[rev_idx] = evt;
// 	      changes.add(rev_idx);
// 	    }
// 	  }
// 	}

// #ifdef _DEBUG_MUL
//       if(IS_OK(wiped)) {
// 	std::cout << scope[0].get_domain() 
// 		  << " * " << scope[1].get_domain() 
// 		  << " = " << scope[2].get_domain() << std::endl;
//       } else std::cout << "FAIL!" << std::endl;
// #endif

//       if(IS_OK(wiped)) {
	
// 	// we update x[1-evt_idx]
// 	rev_idx = 1-evt_idx;
// 	aux_idx = 2;
// 	compute_division(scope[aux_idx], scope[evt_idx], lb, ub);

// #ifdef _DEBUG_MUL
// 	std::cout << "check bounds of " << scope[rev_idx] << " = " 
// 		  << scope[aux_idx] << "/" << scope[evt_idx] 
// 		  << " = [" << lb << ".." << ub << "]" << std::endl; 
// #endif

// 	evt = scope[rev_idx].set_min( lb );
// 	if( IS_FAIL(evt) ) wiped = FAILURE(rev_idx);
// 	else if(evt != NO_EVENT) {
// 	  if(changes.contain(rev_idx)) {
// 	    event_type[rev_idx] |= evt;
// 	  } else {
// 	    event_type[rev_idx] = evt;
// 	    changes.add(rev_idx);
// 	  }
// 	}
	
// 	if(IS_OK(wiped)) {
// 	  evt = scope[rev_idx].set_max( ub );
// 	  if( IS_FAIL(evt) ) wiped = FAILURE(rev_idx);
// 	  else if(evt != NO_EVENT) {
// 	    if(changes.contain(rev_idx)) {
// 	      event_type[rev_idx] |= evt;
// 	    } else {
// 	      event_type[rev_idx] = evt;
// 	      changes.add(rev_idx);
// 	    }
// 	  }
// 	}	
	
// #ifdef _DEBUG_MUL
// 	if(IS_OK(wiped)) {
// 	  std::cout << scope[0].get_domain() 
// 		    << " * " << scope[1].get_domain() 
// 		    << " = " << scope[2].get_domain() << std::endl;
// 	} else std::cout << "FAIL!" << std::endl;
// #endif
	
//       }
//     }
//   }
  
  
  return wiped;
}
  
std::ostream& Mistral::PredicateMul::display(std::ostream& os) const {
  os << scope[2] << " == (" << scope[0] << " * " << scope[1] << ")";
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


Mistral::ConstraintBoolSumInterval::ConstraintBoolSumInterval(Vector< Variable >& scp, const int l, const int u)
  : Constraint(scp) { 
  priority = 1;
  lb = l; 
  ub = u; 
}

void Mistral::ConstraintBoolSumInterval::initialise() {
  Constraint::initialise();
  for(unsigned int i=0; i<scope.size; ++i)
    trigger_on(_VALUE_, i);
  set_idempotent(true);
}

Mistral::ConstraintBoolSumInterval::~ConstraintBoolSumInterval() 
{
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
  os << "(" << scope[0] ;
  for(unsigned int i=1; i<scope.size; ++i) 
    os << " + " << scope[i];
  os << ") in [" << lb << "," << ub << "]" ;
  return os;
}



Mistral::PredicateBoolSum::PredicateBoolSum(Vector< Variable >& scp, Variable tot)
  : Constraint(scp) { 
  add(tot);
  priority = 1;
}

Mistral::PredicateBoolSum::PredicateBoolSum(std::vector< Variable >& scp, Variable tot)
  : Constraint(scp) {
  add(tot);
  priority = 1;
}

Mistral::PredicateBoolSum::PredicateBoolSum(Vector< Variable >& scp)
  : Constraint(scp) { 
  priority = 1;
}

Mistral::PredicateBoolSum::PredicateBoolSum(std::vector< Variable >& scp)
  : Constraint(scp) { 
  priority = 1;
}

void Mistral::PredicateBoolSum::initialise() {
  Constraint::initialise();
  for(unsigned int i=0; i<scope.size-1; ++i)
    trigger_on(_VALUE_, i);
  trigger_on(_RANGE_, scope.size-1);
  set_idempotent(true);
}

Mistral::PredicateBoolSum::~PredicateBoolSum() 
{
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
  os << "(" << scope[0] ;
  for(unsigned int i=1; i<scope.size-1; ++i) 
    os << " + " << scope[i];
  os << ") = " << scope[scope.size-1];
  return os;
}



// Mistral::PredicateWeightedSum::PredicateWeightedSum(Vector< Variable >& scp, 
// 						    Vector< int >& wgt, Variable tot)
//   : Constraint(scp), lower_bound(0), upper_bound(0) { 
//   add(tot);
//   priority = 1;
//   for(unsigned int i=0; i<wgt.size; ++i) {
//     weight.add(wgt[i]);
//   }
//   weight.add(-1);
// }

// Mistral::PredicateWeightedSum::PredicateWeightedSum(std::vector< Variable >& scp, 
// 						    std::vector< int >& wgt, Variable tot)
//   : Constraint(scp), lower_bound(0), upper_bound(0) {
//   add(tot);
//   priority = 1;
//   for(unsigned int i=0; i<wgt.size; ++i) {
//     weight.add(wgt[i]);
//   }
//   weight.add(-1);
// }


Mistral::PredicateWeightedSum::PredicateWeightedSum(Vector< Variable >& scp, 
						    const int L, const int U)
  : Constraint(scp), lower_bound(L), upper_bound(U) { 
  priority = 1;
  for(unsigned int i=0; i<scope.size; ++i) {
    weight.add(1);
  }
}

Mistral::PredicateWeightedSum::PredicateWeightedSum(Vector< Variable >& scp, 
						    Vector< int >& wgt,
						    const int L, const int U)
  : Constraint(scp), lower_bound(L), upper_bound(U) { 
  priority = 1;
  for(unsigned int i=0; i<scope.size; ++i) {
    weight.add(wgt[i]);
  }
}

Mistral::PredicateWeightedSum::PredicateWeightedSum(std::vector< Variable >& scp, 
						    std::vector< int >& wgt,
						    const int L, const int U)
  : Constraint(scp), lower_bound(L), upper_bound(U) { 
  priority = 1;
  for(unsigned int i=0; i<scope.size; ++i) {
    weight.add(wgt[i]);
  }
}

void Mistral::PredicateWeightedSum::initialise() {
  Constraint::initialise();
  for(unsigned int i=0; i<scope.size; ++i)
    trigger_on(_RANGE_, i);
  set_idempotent(true);
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

  lo_bound = new int[scope.size];
  up_bound = new int[scope.size];
  span = new int[scope.size];

  // //std::cout << (int*)env << std::endl;
  // std::cout  << "-- " << (int*)solver << std::endl;

  // exit(1);

  unknown_parity.initialise(0,scope.size-1,true,scope[0].get_solver());
  parity.Reversible::initialise(scope[0].get_solver());
  parity.initialise(lower_bound%2);

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
  delete [] lo_bound;
  delete [] up_bound;
  delete [] span;
}

Mistral::PropagationOutcome Mistral::PredicateWeightedSum::propagate() 
{


  int i, j;
  // compute the max and th min
  int tmin, smin=0, tmax, smax=0// , maxspan=0
    , arity=scope.size;
  //int unknown_parity = -1;
  //int parity = 0;
  PropagationOutcome wiped = CONSISTENT;

  // //if(lower_bound == 72) {
  //   std::cout << std::endl << "propagate ";
  //   for(i=0; i<arity; ++i) {
  //     std::cout << " " << weight[i] << scope[i] << ":" << scope[i].get_domain();
  //   }
  //   std::cout << std::endl << changes << std::endl;
  //   //}

  //unknown_parity.fill();

  for(i=0; i<wpos; ++i) {
    smax += (up_bound[i] = scope[i].get_max());
    smin += (lo_bound[i] = scope[i].get_min());
    span[i] = (up_bound[i]-lo_bound[i]);
    //if( span[i] > maxspan ) maxspan = span[i];

    // if( unknown_parity>-2 ) { // otherwise there are already 2 vars with unknown parity, hence we can't conclude
    //   // parity reasonning: coef is 1, so the parity is known only if the variable is ground
    // if( span[i] == 0 ) {
    //   // the parity changes only if the only one val is odd
    //   if( lo_bound[i]%2 ) parity = 1-parity;
    //   unknown_parity.remove(i);
    // }//  else {
    //   if( unknown_parity==-1 )
    //  	  unknown_parity = i;
    //   else
    // 	unknown_parity = -2;
    // }
    // }

     // if(i)
     //   std::cout << " + [" << lo_bound[i] << "," << up_bound[i] << "] = [" << smin << "," << smax << "] ";
     // else
     //   std::cout << "[" << smin << "," << smax << "] ";
  }
  for(i=wpos; i<wneg; ++i) {
    smax += weight[i] * (up_bound[i] = scope[i].get_max());
    smin += weight[i] * (lo_bound[i] = scope[i].get_min());
    span[i] = weight[i] * (up_bound[i]-lo_bound[i]);
    //if( span[i] > maxspan ) maxspan = span[i];

    // if( unknown_parity>-2 ) { // otherwise there are already 2 vars with unknown parity, hence we can't conclude
    //   // parity reasonning: coef is >1, so the parity is known only if the variable is ground, or if the coef is even
    //if( span[i] == 0 ) {
    // 	// the parity changes only if the only one val is odd
    // 	if( (weight[i]*lo_bound[i])%2 ) parity = 1-parity;
    //   } else if( weight[i]%2 ) {
    // 	if( unknown_parity==-1 )
    // 	  unknown_parity = i;
    // 	else
    // 	  unknown_parity = -2;
    //   }
    // }

    // std::cout << " + [" << lo_bound[i] << "," << up_bound[i] << "] = [" << smin << "," << smax << "] ";
    // std::cout << "[" << smin << "," << smax << "] ";
  }
  for(i=wneg; i<arity; ++i) {
    smax += weight[i] * (lo_bound[i] = scope[i].get_min());
    smin += weight[i] * (up_bound[i] = scope[i].get_max());
    span[i] = weight[i] * (lo_bound[i]-up_bound[i]);
    //if( span[i] > maxspan ) maxspan = span[i];

    // if( unknown_parity>-2 ) { // otherwise there are already 2 vars with unknown parity, hence we can't conclude
    //   // parity reasonning: coef is >1, so the parity is known only if the variable is ground, or if the coef is even
    //   if( span[i] == 0 ) {
    // 	// the parity changes only if the only one val is odd
    // 	if( (weight[i]*lo_bound[i])%2 ) parity = 1-parity;
    //   } else if( weight[i]%2 ) {
    // 	if( unknown_parity==-1 )
    // 	  unknown_parity = i;
    // 	else
    // 	  unknown_parity = -2;
    //   }
    // }

    // std::cout << " + [" << lo_bound[i] << "," << up_bound[i] << "] = [" << smin << "," << smax << "] ";
    // std::cout << "[" << smin << "," << smax << "] ";
  }


  //std::cout << std::endl << events << std::endl << changes << std::endl << IS_OK(wiped) << std::endl;;
  
  while(IS_OK(wiped) && !events.empty()) {

    
    if(lower_bound == upper_bound) {
      j = events.size;
      while( j-- ) {
	// for(j=0; j<events.size; ++j) {
	// 	--j;

	//std::cout << j <<  " > " << events.size << std::endl;
	//std::cout.flush();

	i = events[j];

	//std::cout << j <<  " ? " << events.size << std::endl;
	if(span[i]==0 && unknown_parity.contain(i)) {
	  unknown_parity.reversible_remove(i);
	  if( lo_bound[i]%2 ) parity = 1-parity;
	}

	//std::cout << j <<  " * " << events.size << std::endl;

	//if(j==0) break;
      }

     
      // display(std::cout);
      // std::cout << std::endl << unknown_parity << ": " << (parity ? "odd" : "even") << std::endl;
      
      
      
      if(unknown_parity.size == 0) {
	//std::cout << "parity failure " << std::endl;
	
	if(parity != 0) wiped = FAILURE(arity-1);
      } else if(unknown_parity.size == 1) { // it needs to be of parity "parity"
	i = unknown_parity[0];
	
	//std::cout << "parity pruning: " << (lo_bound[i]%2) << " " << parity << std::endl ;
	
	while(IS_OK(wiped) && (lo_bound[i]%2==0) != (parity==0)) {
	  
	  //std::cout << scope[i] << " in " << scope[i].get_domain() << " => ";
	  
	  tmin = lo_bound[i];
	  if( IS_FAIL(scope[i].set_min(++lo_bound[i])) ) wiped = FAILURE(i);
	  else {
	    lo_bound[i] = scope[i].get_min();
	    if(i<wneg) smin += ((lo_bound[i] - tmin)*weight[i]);
	    else smax += ((lo_bound[i] - tmin)*weight[i]);
	  }
	  //std::cout << scope[i] << " in " << scope[i].get_domain() << std::endl;
	  
	}

	// std::cout << "parity pruning: " << scope[i].get_max() << std::endl ;
	// std::cout << "parity pruning: " << up_bound[i] << std::endl ;
	// std::cout << "parity pruning: " << (up_bound[i]%2) << " " << parity << std::endl ;
	// std::cout << "parity pruning: " << (up_bound[i]%2==0) << " " << (parity==0) << std::endl ;

	while(IS_OK(wiped) && (up_bound[i]%2==0) != (parity==0)) {

	  //std::cout << scope[i] << " in " << scope[i].get_domain() << " => ";

	  tmin = up_bound[i];
	  if( IS_FAIL(scope[i].set_max(--up_bound[i])) ) wiped = FAILURE(i);
	  else {
	    up_bound[i] = scope[i].get_max();
	    if(i<wneg) smax -= ((tmin - up_bound[i])*weight[i]);
	    else smin -= ((tmin - up_bound[i])*weight[i]);
	  }

	  //std::cout << scope[i] << " in " << scope[i].get_domain() << std::endl;
	}
      }
    }

    // //if(lower_bound == 72) {
    // std::cout << "events: " << events << std::endl;
    // //}

    if(IS_OK(wiped)) {

      events.clear();

      // //if(lower_bound == 72) {
      // 	std::cout << " [" << smin << "," << smax << "]" << std::endl;
      // 	//}

      if( smax < lower_bound || smin > upper_bound ) wiped = FAILURE(arity-1);
      else {
	tmax = (smax - lower_bound);
	tmin = (upper_bound - smin);

	//std::cout << "tmin=" << tmin << "  tmax=" << tmax << std::endl;
      
	/// prune with respect to the lower bound?
	//if( tmin < maxspan ) {
	
	for(i=0; IS_OK(wiped) && i<wpos; ++i) {
	  
	  // std::cout << tmin << " <? " << (up_bound[i] - lo_bound[i]) << std::endl;
	  
	  if( tmin < (up_bound[i] - lo_bound[i]) ) {
	    
	    // //if(lower_bound == 72) {
	    //   std::cout << "prune " << scope[i] << " events before " << events << std::endl;
	    //   //}
	    if(IS_FAIL(scope[i].set_max( lo_bound[i] + tmin ))) wiped = FAILURE(i);
	    else {
	      
	      events.add(i);
	      event_type[i] = UB_EVENT;
	      
	      // //if(lower_bound == 72) {
	      // 	std::cout << "events after " << events  << " (add " << i << ")" << std::endl;
	      // 	//}
	    }
	  }
	}
	
	for(i=wpos; IS_OK(wiped) && i<wneg; ++i) {
	  
	  //std::cout << tmin << " <? " << (up_bound[i] - lo_bound[i]) * weight[i] << std::endl;
	  
	  if( tmin < (up_bound[i] - lo_bound[i]) * weight[i] ) {
	    
	    // //if(lower_bound == 72) {
	    //   std::cout << "prune " << scope[i] << " events before " << events << std::endl;
	    //   //}
	    
	    if(IS_FAIL(scope[i].set_max( lo_bound[i] + tmin/weight[i] ))) wiped = FAILURE(i);
	    
	    else {
	      events.add(i);
	      event_type[i] = UB_EVENT;
	      
	      // //if(lower_bound == 72) {
	      // 	std::cout << "events after " << events  << " (add " << i << ")" << std::endl;
	      // 	//}
	      
	    }
	  }
	}
	
	for(i=wneg; IS_OK(wiped) && i<arity; ++i) {
	  
	  //std::cout << tmin << " <? " << (lo_bound[i] - up_bound[i]) * weight[i] << std::endl;
	  
	  if( tmin < (lo_bound[i] - up_bound[i]) * weight[i] ) {
	    
	    // //if(lower_bound == 72) {
	    //   std::cout << "prune " << scope[i] << " events before " << events << std::endl;
	    //   //}
	    
	    if(IS_FAIL(scope[i].set_min( up_bound[i] + tmin/weight[i] ))) wiped = FAILURE(i);
	    
	    else {
	      events.add(i);
	      event_type[i] = LB_EVENT;
	      
	      // //if(lower_bound == 72) {
	      // 	std::cout << "events after " << events  << " (add " << i << ")" << std::endl;
	      // 	//}
	      
	    }
	  }
	}
	//}
	
	/// prune with respect to the upwer bound?
	//if( tmax < maxspan ) {
	
	for(i=0; IS_OK(wiped) && i<wpos; ++i) {
	  
	  // std::cout << tmax << " <? " << (up_bound[i] - lo_bound[i]) << std::endl;
	  
	  if( tmax < (up_bound[i] - lo_bound[i]) ) {
	    
	    // //if(lower_bound == 72) {
	    //   std::cout << "prune " << scope[i] << " events before " << events << std::endl;
	    //   //}
	    
	    if(IS_FAIL(scope[i].set_min( up_bound[i] - tmax ))) wiped = FAILURE(i);
	    
	    else {
	      
	      if(events.contain(i)) {
		event_type[i] |= LB_EVENT;
	      } else {
		events.add(i);
		event_type[i] = LB_EVENT;
	      }


	      // // if(lower_bound == 72) {
	      // 	std::cout << "events after " << events  << " (add " << i << ")" << std::endl;
	      // 	// }
	    }
	  }
	}
	
	for(i=wpos; IS_OK(wiped) && i<wneg; ++i) {
	  
	  //std::cout << tmax << " <? " << (up_bound[i] - lo_bound[i]) * weight[i] << std::endl;
	  
	  if( tmax < (up_bound[i] - lo_bound[i]) * weight[i] ) {
	    
	    // //if(lower_bound == 72) {
	    //   std::cout << "prune " << scope[i] << " events before " << events << std::endl;
	    //   //}
	    
	    if(IS_FAIL(scope[i].set_min( up_bound[i] - tmax/weight[i] ))) wiped = FAILURE(i);
	    
	    else {
	      
	      if(events.contain(i)) {
		event_type[i] |= LB_EVENT;
	      } else {
		events.add(i);
		event_type[i] = LB_EVENT;
	      }
	      
	      // //if(lower_bound == 72) {
	      // 	std::cout << "events after " << events  << " (add " << i << ")" << std::endl;
	      // 	// }
	    }
	  }
	}


	for(i=wneg; IS_OK(wiped) && i<arity; ++i) {
	  
	  //std::cout << tmax << " <? " << ((lo_bound[i] - up_bound[i]) * weight[i]) << std::endl;
	  
	  if( tmax < (lo_bound[i] - up_bound[i]) * weight[i] ) {  
	    
	    // //if(lower_bound == 72) {
	    //   std::cout << "prune " << scope[i] << " events before " << events << std::endl;
	    //   //}
	    
	    if(IS_FAIL(scope[i].set_max( lo_bound[i] - tmax/weight[i] ))) wiped = FAILURE(i);
	    
	    
	    
	    else {
	      
	      if(events.contain(i)) {
		event_type[i] |= UB_EVENT;
	      } else {
		events.add(i);
		event_type[i] = UB_EVENT;
	      }

	      
	      // //if(lower_bound == 72) {
	      // 	std::cout << "events after " << events  << " (add " << i << ")" << std::endl;
	      // 	//}
	    }
	  }
	}
      }
      //}
      
      /// update smin and smax
      for(unsigned int j=0; IS_OK(wiped) && j<events.size; ++j) {
	i = events[j];
	
	// // if(lower_bound == 72) {
	// std::cout << " even on " << scope[i] << " in " <<  scope[i].get_domain() 
	// //     //   // 		<< " " << event_type[i] << " " << LB_EVENT << "/" << UB_EVENT 
	// //     //   // 		<< " " << (event_type[i]&LB_EVENT) 
	// //     //   // 		<< (event_type[i]&UB_EVENT)  
	// 	  << std::endl;
	// // }
	
	if(i<wpos) {
	  if(LB_CHANGED(event_type[i])){ 
	    smin -= lo_bound[i];
	    lo_bound[i] = scope[i].get_min();
	    smin += lo_bound[i];
	  } 
	  if(UB_CHANGED(event_type[i])){ 
	    smax -= up_bound[i];
	    up_bound[i] = scope[i].get_max();
	    smax += up_bound[i];
	  }
	} else if(i<wneg) {
	  if(LB_CHANGED(event_type[i])){ 
	    smin -= (lo_bound[i] * weight[i]);
	    lo_bound[i] = scope[i].get_min();
	    smin += (lo_bound[i] * weight[i]);
	  } 
	  if(UB_CHANGED(event_type[i])){ 
	    smax -= (up_bound[i] * weight[i]);
	    up_bound[i] = scope[i].get_max();
	    smax += (up_bound[i] * weight[i]);
	  }
	} else {
	  if(LB_CHANGED(event_type[i])){ 
	    smax -= (lo_bound[i] * weight[i]);
	    lo_bound[i] = scope[i].get_min();
	    smax += (lo_bound[i] * weight[i]);
	  } 
	  if(UB_CHANGED(event_type[i])){ 
	    smin -= (up_bound[i] * weight[i]);
	    up_bound[i] = scope[i].get_max();
	    smin += (up_bound[i] * weight[i]);
	  }
	}
      }
      //std::cout << std::endl;
    }
  }
  
  // //if(lower_bound == 72) {
  //   std::cout << "result: ";
  //   for(i=0; i<arity; ++i) {
  //     std::cout << " " << weight[i] << scope[i] << ":" << scope[i].get_domain();
  //   }
  //   std::cout << std::endl;
  //   // }

  return wiped;
}

int Mistral::PredicateWeightedSum::check( const int* s ) const 
{
  int i=scope.size, t=0;
  while(i--) t+=(weight[i]*s[i]);
  return (t < lower_bound || t > upper_bound); 
}

std::ostream& Mistral::PredicateWeightedSum::display(std::ostream& os) const {
  if(lower_bound > -INFTY) 
    os << lower_bound << " <= " ;
  os << weight[0] << "*" << scope[0] ;

  for(unsigned int i=1; i<scope.size; ++i) 
    os << " + " << weight[i] << "*" << scope[i];
  
  if(upper_bound < INFTY) 
    os << " <= " << upper_bound;
 

  return os;
}





Mistral::PredicateElement::PredicateElement(Vector< Variable >& scp, const int o)
  : Constraint(scp) {
  offset = o;
  priority = 1;
}

Mistral::PredicateElement::PredicateElement(std::vector< Variable >& scp, const int o)
  : Constraint(scp) { 
  offset = o;
  priority = 1;
}

void Mistral::PredicateElement::initialise() {
  Constraint::initialise();

  int n = scope.size-1;
  aux_dom.initialise( std::min( 0, scope[n].get_min() ), 
		      std::max( n, scope[n].get_max() ), 
		      BitSet::empt );

  for(unsigned int i=0; i<scope.size; ++i)
    trigger_on(_DOMAIN_, i);
  set_idempotent(true);

  /////
  scope[n-1].set_min(0+offset);
  scope[n-1].set_max(n-2+offset);

}

Mistral::PredicateElement::~PredicateElement() 
{ 
}

//#define _DEBUG_ELEMENT true

Mistral::PropagationOutcome Mistral::PredicateElement::propagate() 
{
  PropagationOutcome wiped = CONSISTENT;
  int i, n = scope.size-2, evt, nxt; //, lb, ub, val;
  Variable N = scope[n];
  Variable V = scope[n+1];
  bool update_V = true;

#ifdef _DEBUG_ELEMENT
  std::cout << std::endl << std::endl << "X: " << scope[0].get_domain();
  for(i=1; i<n; ++i) {
    std::cout << " " << scope[i].get_domain();
  }
  std::cout << "[" << scope[n].get_domain() << "-" << offset << "] = " << scope[n+1].get_domain() << std::endl;
#endif

  while(IS_OK(wiped) && update_V) {
    update_V = false;
    while(IS_OK(wiped) && !changes.empty()) {

#ifdef _DEBUG_ELEMENT
           std::cout << changes << std::endl;
#endif
      evt = changes.pop();

#ifdef _DEBUG_ELEMENT
           std::cout << "react to " << scope[evt] << " in " << scope[evt].get_domain() << std::endl;
#endif
      if(evt < n && N.contain(evt)) {

#ifdef _DEBUG_ELEMENT
	std::cout << "  update " << N << " in " << N.get_domain() << ": "
		  << V << " in " << V.get_domain() << " inter "
		  << scope[evt] << " in " << scope[evt].get_domain() << "?" << std::endl;
#endif
	update_V = true;
		// N may change if V changes, or any X changes
	if( !V.intersect(scope[evt]) ) {

	  event_type[n] = N.remove(evt+offset);
	  if( IS_FAIL(event_type[n]) ) {
	    wiped = FAILURE(n);

#ifdef _DEBUG_ELEMENT
	    	    std::cout << "  => FAIL" << std::endl;
#endif	
	  } else if(!changes.contain(n)) {
	    changes.add(n);

#ifdef _DEBUG_ELEMENT
	    	    std::cout << "  => " << N << " in " << N.get_domain() << std::endl;
#endif
	  }
	}
      } else if(evt == n) {
	update_V = true;
	if(ASSIGNED(event_type[n])) {
	  // X may change if N becomes assigned 
	  i = N.get_min()-offset;
	  
#ifdef _DEBUG_ELEMENT
	  std::cout << V << " in " << V.get_domain() << " == "
		    << scope[i] << " in " << scope[i].get_domain() << std::endl;
#endif
	  
	  event_type[i] = scope[i].set_domain(V);
	  if( IS_FAIL(event_type[i]) ) { 
#ifdef _DEBUG_ELEMENT
	    std::cout << "  => FAIL" << std::endl;
#endif
	    wiped = FAILURE(i);
	  } else {
	    if( event_type[i] != NO_EVENT && !changes.contain(i) ) {
	      changes.add(i);
#ifdef _DEBUG_ELEMENT
	      std::cout << "  => " << scope[i] << " in " << scope[i].get_domain() << std::endl;
#endif
	    }
	    event_type[n+1] = V.set_domain(scope[i]);
	    if( IS_FAIL(event_type[n+1]) ) {
#ifdef _DEBUG_ELEMENT
	      std::cout << "  => FAIL" << std::endl;
#endif
	      wiped = FAILURE(n+1);
	  } else if( event_type[n+1] != NO_EVENT && !changes.contain(n+1) ) {
	    changes.add(n+1);
#ifdef _DEBUG_ELEMENT
	    std::cout << "  => " << V << " in " << V.get_domain() << std::endl;
#endif
	    }
	  }
	}
      } else if(evt == n+1) {
#ifdef _DEBUG_ELEMENT
	std::cout << "  update " << N << " in " << N.get_domain() << ": "
		  << V << " in " << V.get_domain() << " inter " << std::endl;
#endif
	
	if( N.is_ground() ) {
#ifdef _DEBUG_ELEMENT
	  std::cout << "  update " << V << " in " << V.get_domain() << std::endl;
#endif
	  
	  i = N.get_min()-offset;
	  event_type[i] = scope[i].set_domain(V);
	  if( IS_FAIL(event_type[i]) ) { 
#ifdef _DEBUG_ELEMENT
	    std::cout << "  => FAIL" << std::endl;
#endif
	    wiped = FAILURE(i);
	  } else {
	    if( event_type[i] != NO_EVENT && !changes.contain(i) ) {
	      changes.add(i);
#ifdef _DEBUG_ELEMENT
	      std::cout << "  => " << scope[i] << " in " << scope[i].get_domain() << std::endl;
#endif
	    }
	  }
	}

	// nxt = N.get_min();
	// do {
	//   i = nxt;
	  
	//   std::cout << " " << i ;
	  
	//   nxt = N.next(i);
	// } while( i<nxt );
	
	// std::cout << std::endl;

	  
	nxt = N.get_min();
	do {
	  i = nxt;
	  nxt = N.next(i);

#ifdef _DEBUG_ELEMENT
	  std::cout << "  " << scope[i-offset] << " in " << scope[i-offset].get_domain() << "?" ;
#endif
	  if( !V.intersect(scope[i-offset]) ) {
#ifdef _DEBUG_ELEMENT
	    std::cout << " NO" << std::endl;
#endif	 
	    event_type[n] = N.remove(i);
	    if( IS_FAIL(event_type[n]) ) {
#ifdef _DEBUG_ELEMENT
	      std::cout << "  => FAIL" << std::endl;
#endif
	      wiped = FAILURE(n);
	    } else if( event_type[n] != NO_EVENT && !changes.contain(n) ) {
	      changes.add(n);
#ifdef _DEBUG_ELEMENT
	      std::cout << "  => " << N << " in " << N.get_domain() << std::endl;
#endif
	    } 
#ifdef _DEBUG_ELEMENT
	    if(event_type[n] == NO_EVENT) {
	      std::cout << "  => NO EVENT" << std::endl;
	    }
#endif	    
	  } 
#ifdef _DEBUG_ELEMENT
	  else {
	    std::cout << " YES" << std::endl;
	  }
#endif
	} while( i<nxt );
      }
    }
    
    if(IS_OK(wiped) && update_V) {
#ifdef _DEBUG_ELEMENT
      std::cout << "  update " << V << " in " << V.get_domain() << std::endl;
#endif
      aux_dom.clear();
      //lb=INFTY; ub=-INFTY;
      
      nxt = N.get_min();
      do {
	i = nxt;
	// if(scope[i-offset].is_range()) {
	//   val = scope[i-offset].get_min();
	//   if(lb>val) lb = val;
	//   val = scope[i-offset].get_max();
	//   if(ub<val) ub = val;
	// } else {
	  scope[i-offset].union_to(aux_dom);
	  //}
	nxt = N.next(i);
#ifdef _DEBUG_ELEMENT
	std::cout << " + " << scope[i-offset] << " in " << scope[i-offset].get_domain() << ": " << aux_dom << std::endl; 
#endif
      } while( i<nxt );
      // if(lb<=ub) aux_dom.fill(lb,ub);
      // std::cout << aux_dom << std::endl;

      event_type[n+1] = V.set_domain(aux_dom);

      if( IS_FAIL(event_type[n+1]) ) { 
#ifdef _DEBUG_ELEMENT
	std::cout << "  => FAIL" << std::endl;
#endif
	wiped = FAILURE(n+1);
      } else if( event_type[n+1] != NO_EVENT && !changes.contain(n+1) ) {
	changes.add(n+1);
#ifdef _DEBUG_ELEMENT
	std::cout << "  => " << V << " in " << V.get_domain() << std::endl;
#endif
      }
    }
  }

#ifdef _DEBUG_ELEMENT
  std::cout << "return " << wiped << std::endl;
  std::cout << "X: " << scope[0].get_domain();
  for(i=1; i<n; ++i) {
    std::cout << " " << scope[i].get_domain();
  }
  std::cout << "[" << scope[n].get_domain() << "] = " << scope[n+1].get_domain() << std::endl << std::endl ;
#endif

  return wiped;
}

int Mistral::PredicateElement::check( const int* s ) const 
{
  return (s[s[scope.size-2]-offset] != s[scope.size-1]);
}

std::ostream& Mistral::PredicateElement::display(std::ostream& os) const {
  os << "(" << scope[0];
  for(unsigned int i=1; i<scope.size-2; ++i) {
    os << " " << scope[i];
  }
  os << ")[" << scope.back(2) << "] == " << scope.back();
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

Mistral::ConstraintAllDiff::ConstraintAllDiff(std::vector< Variable >& scp)
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



// Mistral::ConstraintClauseBase::ConstraintClauseBase(Vector< Variable >& scp) 
//   : Constraint(scp) { }

// void Mistral::ConstraintClauseBase::mark_domain() {
//   for(unsigned int i=0; i<scope.size; ++i) {
//     solver->mark_non_convex(scope[i].id());
//   }
// }

// void Mistral::ConstraintClauseBase::initialise() {
//   Constraint::initialise();
//   for(unsigned int i=0; i<scope.size; ++i) {
//     trigger_on(_VALUE_, i);
//   }
//   set_idempotent(true);

//   is_watched_by.initialise(0,2*scope.size);


//   // for(unsigned int i=0; i<scope.size; ++i) {
//   //   is_watched_by[i] = new Vector< Clause* >;
//   //   is_watched_by[i]-=scope[i].get_min();
//   // }
// }

// Mistral::ConstraintClauseBase::~ConstraintClauseBase() {
// }

// void Mistral::ConstraintClauseBase::add(Variable x) {
//   unsigned int idx = x.id();
//   if(idx == scope.size) {
//     scope.add(x);
//     while(is_watched_by.capacity <= 2*idx)
//       is_watched_by.extendStack();
//   } else if(idx > scope.size) {
//     while(scope.capacity <= idx)
//       scope.extendStack();
//     scope[idx] = x;
//     while(is_watched_by.capacity <= idx)
//       is_watched_by.extendStack();
//   }
// }

// void Mistral::ConstraintClauseBase::add( Vector < Literal >& clause ) {
//  if(clause.size > 1) {
//    Clause *cl = Clause::Array_new(clause);
//    clause.add( cl );
//    is_watched_by[clause[0]].add(cl);
//    is_watched_by[clause[1]].add(cl);
//  } else {
//    scope[UNSIGNED(clause[0])].set_domain(SIGN(clause[0]));
//  }
// }

// int Mistral::ConstraintClauseBase::check( const int* sol ) const {
//   unsigned int i, j;
//   bool satisfied=true;

//   for(i=0; i<clause.size; ++i) {
//     satisfied = true;
//     Clause& clause = *(clause[i]);
//     for(j=0; satisfied && j<clause.size; ++j) {
//       satisfied = //clause[j].check(sol[j]);
// 	(sol[j] == SIGN(clause[j]));
//     }
//   }
  
//   return !satisfied;
// }

// Mistral::PropagationOutcome Mistral::ConstraintClauseBase::propagate() {
//   Clause *conflict=NULL;


//   int x, v, cw;
//   Lit p;
//   //unsigned int i;
//   while( !changes.empty() ) {
//     std::cout << changes << std::endl;
//     //std::cout << scope << std::endl;
//     x = changes.pop();
//     //std::cout << scope[x] << " in " << scope[x].get_domain() << std::endl;
//     v = scope[x].get_min();
//     std::cout << x << "=" << v << ": " << is_watched_by[x][v] << std::endl;

//     p = 2*x+v;

//     cw = is_watched_by[p].size;
//     while(cw-- && !conflict) {
//       conflict = update_watcher(cw, p);
//     }
//   }

//   return CONSISTENT;
// }


// inline Clause* SatSolver::update_watcher(const int cw, const Lit p)
// {
//   Clause *cl = is_watched_by[p][cw];
//   Clause& clause = *cl;

//   unsigned int j;
//   Lit q, r;
//   //Atom v, w;
//   Variable v, w;

//   //ensure that p is the second watched lit
//   if( clause[1] != p ) {
//     q = clause[1];
//     clause[0] = q;
//     clause[1] = p;
//   } else q = clause[0];
//   v=scope[UNSIGNED(q)].;


//   //check if the other watched lit is assigned
//   //if( LEVEL(v) >= assumptions.size || SIGN(v) != SIGN(q) ) {
//   if( !v.is_ground() || v.get_min() != SIGN(q) ) {
//     for(j=2; j<clause.size; ++j) {
//       // for each literal q of the clause,
//       r = clause[j];
//       w = scope[UNSIGNED(r)];

//       if( !w.is_ground() ) { // this literal is not set
// 	// then it is a good candidate to replace p
// 	clause[1] = r;
// 	clause[j] = p;
// 	is_watched_by[p].remove(cw);
// 	is_watched_by[r].add(cl);

// 	break;	
//       }
//       // if it is set true, then the clause is satisfied
//       else if( w.get_min() == SIGN(r) ) {
// 	break;
//       }
//     }
      
//     if( j == clause.size ) // no replacement could be found
//       { 
// 	if( !v,is_ground() ) {
// 	  // the last literal (other watched lit) is not set yet, we set it
// 	  add_lit(q);
// 	  //reason[UNSIGNED(q)] = cl;
// 	} else 
// 	  // it is set to false already, we fail
// 	  if( v.get_min() != SIGN(q) ) {
// 	    return cl;
// 	  }
//       }
//   }

//   return NULL;
// }




// std::ostream& Mistral::ConstraintClauseBase::display(std::ostream& os) const {
//   os << " (";
//   if(clauses.size>0) {
//     os << clauses[0];
//     for(unsigned int i=1; i<clauses.size; ++i)
//       os << " " << clauses[i]  ;
//   }
//   os << ")";
//   return os;
// }



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
