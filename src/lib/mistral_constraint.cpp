
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


#include <mistral_constraint.hpp>
#include <mistral_variable.hpp>


Mistral::Constraint::Constraint() {
  priority = 2;
  id = -1;
  weight = 1.0;
}

Mistral::Constraint::Constraint(Vector< IntVar >& scp, const float w) {
  priority = 2;
  id = -1;
  weight = 1.0;
  initialise(scp, w);
}

Mistral::Constraint::~Constraint() {
  if(changes.list_ == events.list_)
    events.list_ = NULL;

  delete [] solution;
  delete [] self;
  delete [] trigger;
  delete [] event_type;
  delete [] _scope;
  scope = NULL;
  
}

void Mistral::Constraint::initialise(Vector< IntVar >& scp, const float w) {
  arity = scp.size;
  weight = w;

  _scope = new IntVar[arity];
  event_type = new Event[arity];
  self = new int[arity];
  trigger = new int[arity];
  solution = new int[arity];
  changes.initialise(0, arity-1, false);
  supports = NULL;

  std::memcpy(_scope, scp.stack_, arity*sizeof(IntVar));
  std::fill(event_type, event_type+arity, NO_EVENT);
  std::fill(self, self+arity, -1);
  std::fill(trigger, trigger+arity, _domain_);
  std::fill(solution, solution+arity, 0);

  active.initialise(0, arity-1, true);
  for(int i=0; i<arity; ++i) 
    if(_scope[i]->isGround()) active.erase(i);

  scope = _scope;
}

void Mistral::Constraint::notifyAssignment(const int var) {
  assert(active.member(var));

  active.erase(var);
  if(active.size == 1) entail();
}

void Mistral::Constraint::entail() {
  unsigned int i, j;
  for(i=0; i<active.size; ++i) {
    j = active[i];
    _scope[j]->constraints.reversible_erase(self[j]);
  }
}

void Mistral::Constraint::restore() {
  active.size = trail_act.pop();
}

void Mistral::Constraint::triggerOn(const int t, const int x) {
  trigger[x] = t;
  //self_list[x] = &(_scope[x]->constraints);
  ConstraintTrigger ct(this, x);
  self[x] = _scope[x]->constraints.create(ct, t);
}

void Mistral::Constraint::post() {
  for(int i=0; i<arity; ++i)
    _scope[i]->constraints.reversible_add(self[i], trigger[i]);
}
 
void Mistral::Constraint::relax() {
  for(int i=0; i<arity; ++i)
    _scope[i]->constraints.reversible_erase(self[i], trigger[i]);
}

void Mistral::Constraint::assign(const int var) {
  active.erase(var);
  if(active.size == 1) {
    int i = active.back();
    _scope[i]->constraints.erase(self[i], trigger[i]);
  }
}

void Mistral::Constraint::unassign(const int var) {
  if(active.size == 1) {
    int i = active.back();
    _scope[i]->constraints.insert(self[i], trigger[i]);
  } 
  if(!(active.member(var))) active.insert(var);
}

bool Mistral::Constraint::firstSupport(const int vri, const int vli) 
{
  int j;
  if( supports && supports[vri][vli][0] != NOVAL ) {
    j=arity;
    while( j-- ) 
      if( vri != j )
	if (!scope[j]->contain( supports[vri][vli][j] )) break;
    if( j < 0 ) 
      return true;
  } 
  j=arity;
  while( j-- )
    solution[j] = scope[j]->domain.min;
  solution[vri] = vli; 

  return false;
}

bool Mistral::Constraint::findSupport(const int vri, const int vli) 
{
  int i=arity, vali;
  bool found=false;
  // sol is initialized: either to the value 
  // a variable is already assigned to
  // or to the first value in its domain
  while(i >= 0) {
    // check this assignment
    if( !check( solution ) ) {
      found=true;
      if( supports ) {
	vali = arity;
	while( vali-- )
	  supports[vri][vli][vali] = solution[vali];
      }
      break;
    }
 
    // try to assign more things
    // find the last var whose domain we have not exhausted
    --i;
    while( i >= 0 ) {
      if( i == vri || scope[i]->isGround() ) {
	--i;
	continue;
      }
      solution[i] = scope[i]->next( solution[i] );
      if(solution[i] != NOVAL)
      	break;
      else
	solution[i] = scope[i]->domain.min;
      --i;
    }
    if( i >= 0 )
      i = arity;
  } 
  return found;
}


// bool Mistral::Constraint::propagate()
// {
//   int i, consistent=1, evt = ( Constraint::RANGETRIGGER );
//   //int i, consistent=1, evt = ( Constraint::VALUETRIGGER );
//   for( i=0; consistent && i<arity; ++i )
//     consistent = propagate( i, evt );
//   return consistent;
// } 

// std::string Mistral::Constraint::getString() const {
//   std::string return_str = name()+"("+(scope[0]->getString());
//   for(int i=1; i<arity; ++i)
//     return_str += (", "+(scope[1]->getString()));
//   return return_str+")";
// }

std::ostream& Mistral::Constraint::display(std::ostream& os) const {
  os << name() << "(" << scope[0];
  for(int i=1; i<arity; ++i)
    os << ", " << scope[1];
  os << ")";
  return os;
}

Mistral::IntVar Mistral::ConstraintNotEqual::propagate() {
  int var = 1-changes[0];
  if(active.size)
    return(scope[var]->remove(scope[1-var]->domain.min) == FAIL_EVENT ? scope[var] : NULL);
  else
    return(scope[0]->domain.min == scope[1]->domain.min ? scope[var] : NULL);
}

// std::string Mistral::ConstraintNotEqual::getString() const {
//   return (toString(scope[0])+" =/= "+toString(scope[1]));
// }

std::ostream& Mistral::ConstraintNotEqual::display(std::ostream& os) const {
  os << scope[0] << " =/= " << scope[1];
  return os;
}

Mistral::IntVar Mistral::ConstraintLess::propagate() {
  if(changes.member(0) && (trigger[0] & LB_EVENT)) {
    if(scope[1]->setMin(scope[0]->domain.min + offset) == FAIL_EVENT) return scope[1];
  }
  if(changes.member(1) && (trigger[1] & UB_EVENT)) {
    if(scope[0]->setMax(scope[1]->domain.max - offset) == FAIL_EVENT) return scope[0];
  }
  return NULL;
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

Mistral::IntVar Mistral::PredicateEqual::propagate() {      
  Mistral::IntVar not_consistent = NULL;

  if( scope[2]->isGround() ) {
    if( (spin + scope[2]->domain.min) != 1 ) {
      if( !(scope[0]->setDomain(scope[1])) ) not_consistent = scope[0];
      else if( !(scope[1]->setDomain(scope[0])) ) not_consistent = scope[1];
    } else {
      if(scope[0]->isGround() && !(scope[1]->remove(scope[0]->domain.min)))
	not_consistent = scope[1];
      else {
	if(scope[1]->isGround() && !(scope[0]->remove(scope[1]->domain.min)))
	  not_consistent = scope[0];
      }
    }
  } else {
    if( !(scope[0]->intersect(scope[1])) ) {
      if(!(scope[2]->remove(spin))) not_consistent = scope[2];	    
    } else { 
      if( scope[0]->isGround() && scope[1]->isGround() ) {
	if(!(scope[2]->setDomain( spin ))) not_consistent = scope[2];
      }
    }
  }
  
  return not_consistent;
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

Mistral::IntVar Mistral::PredicateAdd::propagate() {      
  Mistral::IntVar wiped_out = NULL;
  int i;

  // update scope[0] and scope[1]
  for(i=0; i<2; ++i) {
    //if(changes.member(2) || changes.member(i)) {
    //if(is_lower_bound(evt_type[2]) || is_upper_bound(evt_type[i]))
    if(scope[1-i]->setMin(scope[2]->domain.min - scope[i]->domain.max) == FAIL_EVENT) 
      { wiped_out = scope[1-i]; break; }
    //if(is_upper_bound(evt_type[2]) || is_lower_bound(evt_type[i]))
    if(scope[1-i]->setMax(scope[2]->domain.max - scope[i]->domain.min) == FAIL_EVENT) 
      { wiped_out = scope[1-i]; break; }
    //}
  }
  
  if(scope[2]->setMin(scope[0]->domain.min + scope[1]->domain.min) == FAIL_EVENT) 
    { wiped_out = scope[2]; }
  else if(scope[2]->setMax(scope[0]->domain.max + scope[1]->domain.max) == FAIL_EVENT) 
    { wiped_out = scope[2]; }
  
  return wiped_out;
}

// std::string Mistral::PredicateAdd::getString() const {
//   return (toString(scope[2])+" = ("+toString(scope[0])+" + "+toString(scope[1])+")");
// }

std::ostream& Mistral::PredicateAdd::display(std::ostream& os) const {
  os << scope[2] << " = (" << scope[0] << " + " << scope[1] << ")";
  return os;
}

/**********************************************
 * AllDiff Constraint 
 **********************************************/

const int INCONSISTENT = 0;
const int CHANGES      = 1;
const int NO_CHANGES   = 2;

Mistral::ConstraintAllDiff::ConstraintAllDiff(Vector< IntVar >& scp)
  : Constraint(scp) {
  for(int i=0; i<arity; ++i)
    triggerOn(_range_, i);
  set_idempotent(true);
  priority = 0;
  level = &(scp[0]->solver->level);
  init();
}

void Mistral::ConstraintAllDiff::init() 
{
  int i;
  lastLevel = -1;
  nb = 0;

  iv        = new Interval[arity];
  minsorted = new Interval*[arity];
  maxsorted = new Interval*[arity];
  bounds    = new int[2*arity+2];
  std::fill(bounds, bounds+2*arity+2, 0);

  for( i=0; i<arity; ++i ) {
    minsorted[i] = maxsorted[i] = &iv[i];  
    iv[i].min = iv[i].max = NOVAL;
  }

  t = new int[2*arity+2];
  d = new int[2*arity+2];
  h = new int[2*arity+2];
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
  int i,j,nb,min,max,last;

  sortmin(minsorted, arity);
  sortmax(maxsorted, arity);

  min = minsorted[0]->min;
  max = maxsorted[0]->max + 1;
  bounds[0] = last = min-2;

  for (i=j=nb=0;;) { // merge minsorted[] and maxsorted[] into bounds[]
    if (i<arity && min<=max) {	// make sure minsorted exhausted first
      if (min != last)
        bounds[++nb] = last = min;
      minsorted[i]->minrank = nb;
      if (++i < arity)
        min = minsorted[i]->min;
    } else {
      if (max != last)
	bounds[++nb] = last = max;
      maxsorted[j]->maxrank = nb;
      if (++j == arity) break;
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
  int i,j,w,x,y,z;
  int changes = 0;

  for (i=1; i<=nb+1; i++)
    d[i] = bounds[i] - bounds[t[i]=h[i]=i-1];
  for (i=0; i<arity; i++) { // visit Intervals in increasing max order
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
  int i,j,w,x,y,z;
  int changes = 0;

  for (i=0; i<=nb; i++)
    d[i] = bounds[t[i]=h[i]=i+1] - bounds[i];
  for (i=arity; --i>=0; ) { // visit Intervals in decreasing min order
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

Mistral::IntVar Mistral::ConstraintAllDiff::propagate() 
{
  int i, a, b;

  int status_lower, status_upper;
  int l, u;

  a = 0;
  b = arity;

  //if( lastLevel != ((solver->level) - 1) ) {
  if( lastLevel != ((*level) - 1) ) {
    // not incremental
    status_lower = CHANGES;
    status_upper = CHANGES;
    i = 0;
    while (i < arity) {
      iv[i].min = scope[i]->domain.min;
      iv[i].max = scope[i]->domain.max;
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
      iv[i].min = scope[i]->domain.min;
      iv[i].max = scope[i]->domain.max;
      if( l != iv[i].min ) status_lower = CHANGES;
      if( u != iv[i].max ) status_upper = CHANGES;
    }
  }

  lastLevel = *level;//(solver->level);
  //lastLevel = (solver->level);

  if( status_lower == NO_CHANGES && status_upper == NO_CHANGES ) 
    return NULL;

  sortit();

  status_lower = filterlower();
  if( status_lower != INCONSISTENT )
    status_upper = filterupper();  

  if( (status_lower == INCONSISTENT) || (status_upper == INCONSISTENT) ) 
    { return scope[changes.back()]; }
  else
    if( (status_lower == CHANGES) || (status_upper == CHANGES) ) {
      i = 0;
      while (i < arity) {
	if( !scope[i]->setMin( iv[i].min ) )  { return scope[i]; }
	if( !scope[i]->setMax( iv[i].max ) )  { return scope[i]; }
	i++;
      }
    }  
  return NULL;
}

int Mistral::ConstraintAllDiff::check( const int* s ) const 
{
  int i=arity, j;
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
//   for(int i=1; i<arity; ++i) 
//     return_str += (" ,"+toString(scope[i]));
//   return (return_str+")");
// }

std::ostream& Mistral::ConstraintAllDiff::display(std::ostream& os) const {
  os << "alldiff(" << scope[0] ;
  for(int i=1; i<arity; ++i) 
    os << " ," << scope[i];
  os << ")" ;
  return os;
}
