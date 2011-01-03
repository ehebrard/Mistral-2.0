
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
    if(scope[i].is_ground()) active.erase(i);

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

#ifdef _DEBUG_AC
  std::cout << "relax " << this << " from ";
#endif


  for(i=0; i<active.size; ++i) {
    j = active[i];
    k = scope[j].id();

#ifdef _DEBUG_AC
    std::cout << scope[j] << "'s c-list " ;
#endif

    solver->constraint_graph[k]->reversible_erase(self[j], trigger[j]);
  }

#ifdef _DEBUG_AC
  std::cout << std::endl;
#endif
  
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


// void Mistral::ConstraintLess::post() {
//   for(unsigned int i=0; i<scope.size; ++i)
//     trigger_on(_range_, i);
// }

Mistral::PropagationOutcome Mistral::ConstraintLess::propagate() {
  if(changes.contain(0) && LB_CHANGED(event_type[0])) {
    if(scope[1].setMin(scope[0].get_min() + offset) == FAIL_EVENT) 
      return FAILURE(1);
  }
  if(changes.contain(1) && UB_CHANGED(event_type[1])) {
    if(scope[0].setMax(scope[1].get_max() - offset) == FAIL_EVENT) 
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
      if( !(scope[0].setDomain(scope[1])) ) wiped = FAILURE(0);
      else if( !(scope[1].setDomain(scope[0])) ) wiped = FAILURE(1);
    } else {
      if(scope[0].is_ground() && !(scope[1].remove(scope[0].get_min())))
	wiped = FAILURE(1);
      else {
	if(scope[1].is_ground() && !(scope[0].remove(scope[1].get_min())))
	  wiped = FAILURE(0);
      }
    }
  } else {
    if( !(scope[0].intersect(scope[1])) ) {
      if(!(scope[2].remove(spin))) wiped = FAILURE(2);	    
    } else { 
      if( scope[0].is_ground() && scope[1].is_ground() ) {
	if(!(scope[2].setDomain( spin ))) wiped = FAILURE(2);
      }
    }
  }
  
  return wiped;
}

// std::string Mistral::PredicateEqual::getString() const {
//   std::string return_str = (toString(scope[2])+" <=> (");
//   if(spin) return_str += (toString(scope[0])+" == "+toString(scope[1]));
//   else return_str += (toString(scope[0])+" =/= "+toString(scope[1]));
//   return_str += ")";
//   return return_str;
// }

std::ostream& Mistral::PredicateEqual::display(std::ostream& os) const {
  os << scope[2] << " <=> (";
  if(spin) os << scope[0] << " == " << scope[1];
  else os << scope[0] << " =/= " << scope[1];
  os << ")";
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
    //   std::cout << "prune " << scope[i] << " in " << scope[i].get_domain() << std::endl;

    if(i!=2) {
      evt |= scope[i].setMin(scope[2].get_min() - scope[1-i].get_max());
      if(IS_FAIL(evt)) { wiped = FAILURE(i); break; }
      evt |= scope[i].setMax(scope[2].get_max() - scope[1-i].get_min()); 
      if(IS_FAIL(evt)) { wiped = FAILURE(i); break; }
    } else {
      evt |= scope[2].setMin(scope[0].get_min() + scope[1].get_min()); 
      if(IS_FAIL(evt)) { wiped = FAILURE(2); break; }
      evt |= scope[2].setMax(scope[0].get_max() + scope[1].get_max()); 
      if(IS_FAIL(evt)) { wiped = FAILURE(2); break; }
    }

//     std::cout << "=> " << scope[i] << " in " << scope[i].get_domain() 
// 	      << " (" << evt << ")" << std::endl;

    if(evt != NO_EVENT) last_change = i;
    else if(last_change == (i+1)%3) break;
  }


//   Mistral::PropagationOutcome wiped = CONSISTENT;
//   int i, n_changes=2, min[3], max[3], bmin[3], bmax[3], smin, smax;
//   std::cout << std::endl;
//   for(i=0; i<3; ++i) {
//     bmin[i] = min[i] = scope[i].get_min();
//     bmax[i] = max[i] = scope[i].get_max();
//     std::cout << min[i] << "--" << max[i] << std::endl;
//   }

//   while(wiped == CONSISTENT && n_changes >= 2) {
//     n_changes = 0;
//     for(i=0; i<2; ++i) {
//       min[i] = min[2]-max[1-i];
//       max[i] = max[2]-min[1-i];
//       if(min[i] > bmax[i] || bmin[i] > max[i]) { wiped = FAILURE(i); break; }
//     }
//     min[2] = min[0]+min[1];
//     max[2] = max[0]+max[1];
//     if(min[2] > bmax[2] || bmin[2] > max[2]) { wiped = FAILURE(2); break; }

//      std::cout << "[" << bmin[0] << "," << bmax[0] << "] + "
// 	      << "[" << bmin[1] << "," << bmax[1] << "] = "
// 	      << "[" << bmin[2] << "," << bmax[2] << "]" << std::endl;

//     std::cout << "[" << min[0] << "," << max[0] << "] + "
// 	      << "[" << min[1] << "," << max[1] << "] = "
// 	      << "[" << min[2] << "," << max[2] << "]" << std::endl;

//     for(i=0; i<3; ++i) {      
//       smin = bmin[i] < min[i];
//       smax = bmax[i] > max[i];
//       if(smin) {
// 	scope[i].setMin(min[i]);
// 	bmin[i] = min[i];
//       } else min[i] = bmin[i];
//       if(smax) {
// 	scope[i].setMax(max[i]);
// 	bmax[i] = max[i];
//       } else max[i] = bmax[i];
//       n_changes += (smin || smax);
//     }
//     std::cout << n_changes << std::endl;
//   }

//   while(wiped == CONSISTENT && n_changes < 2)
//   // update scope[0] and scope[1]
//   for(i=0; i<2; ++i) {
//     //if(changes.contain(2) || changes.contain(i)) {
//     //if(is_lower_bound(evt_type[2]) || is_upper_bound(evt_type[i]))
//     //new_bound = scope[2].get_min() - scope[i].get_max();
//     //if(scope[1-i].setMin() == FAIL_EVENT) 
//     if(scope[1-i].setMin(scope[2].get_min() - scope[i].get_max()) == FAIL_EVENT) 
//       { wiped = FAILURE(1-i); break; }
//     //if(is_upper_bound(evt_type[2]) || is_lower_bound(evt_type[i]))
//     if(scope[1-i].setMax(scope[2].get_max() - scope[i].get_min()) == FAIL_EVENT) 
//       { wiped = FAILURE(1-i); break; }
//     //}
//   }
  
//   if(scope[2].setMin(scope[0].get_min() + scope[1].get_min()) == FAIL_EVENT) 
//     { wiped = FAILURE(2); }
//   else if(scope[2].setMax(scope[0].get_max() + scope[1].get_max()) == FAIL_EVENT) 
//     { wiped = FAILURE(2); }
  
  return wiped;
}

// std::string Mistral::PredicateAdd::getString() const {
//   return (toString(scope[2])+" = ("+toString(scope[0])+" + "+toString(scope[1])+")");
// }

std::ostream& Mistral::PredicateAdd::display(std::ostream& os) const {
  os << scope[2] << " = (" << scope[0] << " + " << scope[1] << ")";
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

std::ostream& Mistral::ConstraintCliqueNotEqual::display(std::ostream& os) const {
  os << "cliqueNE(" << scope[0] ;
  for(unsigned int i=1; i<scope.size; ++i) 
    os << " ," << scope[i];
  os << ")" ;
  return os;
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
  for(unsigned int i=0; i<scope.size; ++i)
    trigger_on(_range_, i);
  set_idempotent(true);
  priority = 0;
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
	if( scope[i].setMin( iv[i].min ) == FAIL_EVENT )  { return FAILURE(i); }
	if( scope[i].setMax( iv[i].max ) == FAIL_EVENT )  { return FAILURE(i); }
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
