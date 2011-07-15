
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
#include <mistral_solver.hpp>



Mistral::Variable::Variable() {
  domain_type = DYN_VAR;
  variable = NULL;
}

Mistral::Variable::Variable(const int value) {
  domain_type = CONST_VAR;
  constant_value = value;
}

Mistral::Variable::Variable(VariableImplementation* impl, const int type) {
  domain_type = type;
  variable = impl;
}

Mistral::Variable::Variable(Expression* exp) {
  domain_type = EXPRESSION;
  expression = exp;
}

int BOOL_DOM = 3;

Mistral::Variable::Variable(Vector< int >& values, const int type) {

  //std::cout << values << std::endl;

  if(values.back() - values.front() + 1 == (int)(values.size)) {
    initialise_domain(values.front(), values.back(), type);
  } else {
    initialise_domain(values, type);

    //std::cout << "variable constructor cannot yet take a vector of values as argument!" << std::endl;
    //exit(1);
  }
}

Mistral::Variable::Variable(const int lo, const int up, const int type) {
  initialise_domain(lo, up, type);
}

Mistral::Variable::Variable(Variable X, bool h) {
  if(X.domain_type == RANGE_VAR) {
    int lb = X.get_initial_min();
    int ub = X.get_initial_max();
    initialise_domain(lb, ub, BITSET_VAR);
    variable->id = X.variable->id;
    variable->solver = X.variable->solver;
    ((VariableRange*)X.variable)->set_history((VariableBitmap*)variable);
  }
}

void Mistral::Variable::initialise_domain(const int lo, const int up, const int type) {
  if(lo == up) {
    domain_type = CONST_VAR;
    constant_value = lo;
  } else if(type == EXPRESSION) {
    domain_type = EXPRESSION;
    expression = new Expression(lo, up);
  } else if((type & BOOL_VAR) && lo==0 && up==1) {
    bool_domain = &BOOL_DOM;
    variable = new VariableImplementation();
  } else if(type & RANGE_VAR) {
    domain_type = RANGE_VAR;
    range_domain = new VariableRange(lo, up);
  } else {
    domain_type = BITSET_VAR;

    int nwords = 1+(up >> BitSet::EXP)-(lo >> BitSet::EXP);
#ifdef _BIT64
    // int nwords = ((up>=0 ? ((up+1)/64 + ((up+1)%64)!=0) : 0) + 
    // 		  (lo<0 ? (-lo/64 + ((lo%64)!=0)) : 0)) ;//((up-lo) / 64);
    if(nwords == 1) bitset_domain = new VariableWord<unsigned long long int, 1>(lo, up);
    else if(nwords == 2) bitset_domain = new VariableWord<unsigned long long int, 2>(lo, up);
    else if(nwords == 3) bitset_domain = new VariableWord<unsigned long long int, 3>(lo, up);
#else
    // int nwords = ((up>=0 ? ((up+1)/32 + ((up+1)%32)!=0) : 0) + 
    // 		  (lo<0 ? (-lo/32 + ((lo%32)!=0)) : 0)) ;//((up-lo) / 32);

    
    //std::cout <<lo << "-" << up << " " <<  nwords << std::endl;

    //int nwords = ((up-lo) / 32);
    if(nwords == 1) bitset_domain = new VariableWord<unsigned int, 1>(lo, up);
    else if(nwords == 2) bitset_domain = new VariableWord<unsigned int, 2>(lo, up);
    else if(nwords == 3) bitset_domain = new VariableWord<unsigned int, 3>(lo, up);
    else if(nwords == 4) bitset_domain = new VariableWord<unsigned int, 4>(lo, up);
    else if(nwords == 5) bitset_domain = new VariableWord<unsigned int, 5>(lo, up);
#endif
    else bitset_domain = new VariableBitmap(lo, up);
  }
}


void Mistral::Variable::initialise_domain(Vector< int >& values, const int type) {
  if(type == EXPRESSION) {
    domain_type = EXPRESSION;
    expression = new Expression(values);
  } else {
    domain_type = BITSET_VAR;
    
    int nwords = 1+(values[0] >> BitSet::EXP)-(values.back() >> BitSet::EXP);
#ifdef _BIT64
    if(nwords == 1) bitset_domain = new VariableWord<unsigned long long int, 1>(values);
    else if(nwords == 2) bitset_domain = new VariableWord<unsigned long long int, 2>(values);
    else if(nwords == 3) bitset_domain = new VariableWord<unsigned long long int, 3>(values);
#else
    if(nwords == 1) bitset_domain = new VariableWord<unsigned int, 1>(values);
    else if(nwords == 2) bitset_domain = new VariableWord<unsigned int, 2>(values);
    else if(nwords == 3) bitset_domain = new VariableWord<unsigned int, 3>(values);
    else if(nwords == 4) bitset_domain = new VariableWord<unsigned int, 4>(values);
    else if(nwords == 5) bitset_domain = new VariableWord<unsigned int, 5>(values);
#endif
    else bitset_domain = new VariableBitmap(values);
  }
}



Mistral::Variable Mistral::Variable::get_var() {
  if(domain_type == EXPRESSION) {
    return expression->self.get_var();
  } else if(domain_type == CONST_VAR || !variable->solver) {
    return *this;
  }
  return variable->solver->variables[variable->id];
}

std::ostream& Mistral::Variable::display(std::ostream& os) const {
  if(domain_type == EXPRESSION) {
    expression->display(os);
    //}
  } else if(domain_type == CONST_VAR) {
    os << constant_value;
  } else {
    int id = variable->id;
    if       (domain_type ==  BITSET_VAR) {
      os << "x" << id // << (VariableBitmap*)displayDomain(os)
	;
    } // else if(domain_type ==    LIST_VAR) {
      //       os << "y" << id;
      //     } 
    else if(domain_type ==   RANGE_VAR) {
      os << "r" << id;
    }
      //     } else if(domain_type == VIRTUAL_VAR) {
      //       return ((VariableVirtual *)implementation)->display(os);
      //     }
    else  {
      os << "b" << id;
    }
  }
  return os;
  
  //if(*((int*)domain_type)  ==   CONST_VAR) os << constant_value;
  //else os = implementation->display(os);
  //return os;
}

void Mistral::Variable::initialise(Solver *s, const bool top) {
  if(domain_type == EXPRESSION) {
    if(!expression->is_initialised()) {
      for(unsigned int i=0; i<expression->children.size; ++i) {
	expression->children[i].initialise(s, false);
      }

      if(top && !expression->children.empty()) {
	expression->extract_constraint(s);
      } else {
	expression->extract_variable(s);
	expression->id = expression->self.id();
	expression->solver = s;	
	expression->extract_predicate(s);//);
      }
      s->expression_store.add(expression);
    }
  } else {
    if(domain_type != CONST_VAR && variable->solver != s) {
      s->declare(this->get_var());
      s->sequence.declare(*this);
    } 
  }
}



// Mistral::Event Mistral::Variable::setValue( const int val ) 
//   {
// //     std::cout << "SET VALUE" << std::endl;
// //     std::cout << (int*)domain_type << std::endl;
// //     std::cout << *bool_domain << std::endl;

//     //int nstat = (*bool_domain & val);
//     //std::cout << (*bool_domain) << "&" << val << ": " << nstat << std::endl;
//     int dom = *bool_domain;

//     if( val == dom ) return NO_EVENT;
//     if( val < 1 || val > 2 || !(val&dom) ) return FAIL_EVENT;

//     *bool_domain = val;

//     variable->trigger_value_event_and_save(this);
//     return VALUE_EVENT;
//   }

Mistral::Event Mistral::Variable::setValue( const int val ) 
{
  // val should be 1 or 2

  int dom = *bool_domain;
  //int nstat = val&3;

  //std::cout << "1 set " << id() << " dom (" << dom << ") to " << val << std::endl;
  
  if( val == dom ) return NO_EVENT;
  else if( dom<3 // || val<1 || val>2
	   ) return FAIL_EVENT;
  
  *bool_domain = val;
  //std::cout << "1 VALUE EVENT" << std::endl;
  
  //variable->trigger_value_event_and_save(this);
  variable->trigger_value_event_and_save();


  return VALUE_EVENT;
}


Mistral::Event Mistral::Variable::setState( const int vals ) 
{
  // vals should be 1, 2 or 3

  int dom = *bool_domain;
  int ndom = vals&dom;
  //int nstat = val&3;

  //std::cout << "2 set " << id() << " dom (" << dom << ") to " << vals << std::endl;
  
  if( ndom == dom ) return NO_EVENT;
  else if( !ndom ) return FAIL_EVENT;
  
  *bool_domain = ndom;
  //std::cout << "2 VALUE EVENT" << std::endl;
  
  //variable->trigger_value_event_and_save(this);
  variable->trigger_value_event_and_save();

  return VALUE_EVENT;
}

int Mistral::Variable::get_solution_value() const {
  return(domain_type == CONST_VAR ? constant_value :
	 variable->get_solution_value());
}

int Mistral::Variable::get_solution_min() const {
  return(domain_type == CONST_VAR ? constant_value :
	 variable->get_solution_min());
}

int Mistral::Variable::get_solution_max() const {
  return(domain_type == CONST_VAR ? constant_value :
	 variable->get_solution_max());
}

  int Mistral::Variable::get_value() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->get_value();
    else if(domain_type ==    LIST_VAR) return list_domain->get_value();
    else if(domain_type ==   RANGE_VAR) return range_domain->get_value();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_value();
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.get_value();
    else  return (*bool_domain-1);
  }

unsigned int Mistral::Variable::get_size() const {
  if     (domain_type ==  BITSET_VAR) return bitset_domain->get_size();
  else if(domain_type ==    LIST_VAR) return list_domain->get_size();
  else if(domain_type ==   RANGE_VAR) return range_domain->get_size();
  //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_size();
  else if(domain_type ==   CONST_VAR) return 1;
  else if(domain_type ==   EXPRESSION) return expression->self.get_size();
  else  return ((*bool_domain+1)/2);
}

/// Returns the degree (number of constraints)
unsigned int Mistral::Variable::get_degree() const {
  if(domain_type ==   CONST_VAR) return 0;
  else if(domain_type ==   EXPRESSION) return expression->self.get_degree();
  else return variable->solver->constraint_graph[variable->id]->degree;
}


std::string Mistral::Variable::get_domain() const {
  std::ostringstream buf;
  if     (domain_type ==  BITSET_VAR) buf << bitset_domain->domain;
  else if(domain_type ==    LIST_VAR) buf << list_domain->domain;
  else if(domain_type ==   RANGE_VAR) {
    //Mistral::VariableRange* r = (Mistral::VariableRange*)implementation;
    if(range_domain->get_min() == range_domain->get_max())
      buf << "[" << range_domain->get_min() << "]";
    // else if(range_domain->get_min() == range_domain->get_max()-1)
    //   buf << "{" << range_domain->get_min() << "," <<  range_domain->get_max() << "}";
    else 
      buf << "[" << range_domain->get_min() << ".." <<  range_domain->get_max() << "]";
  }
  //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_domain();
  else if(domain_type ==   CONST_VAR)  {
    buf << "{" << constant_value << "}";
  }
  else if(domain_type ==   BOOL_VAR)  {
    buf << "{0,1}";
  } 
  else if(domain_type ==   EXPRESSION) return expression->self.get_domain();
  else {
    if(*bool_domain == 3) buf << "{0,1}";
    else if(*bool_domain == 2) buf << "{1}";
    else buf << "{0}";
  }
  return buf.str();
}

  int Mistral::Variable::get_min() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->get_min();
    else if(domain_type ==    LIST_VAR) return list_domain->get_min();
    else if(domain_type ==   RANGE_VAR) return range_domain->get_min();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_min();
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.get_min();
    else  return (!(*bool_domain & 1));
  }

  int Mistral::Variable::get_max() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->get_max();
    else if(domain_type ==    LIST_VAR) return list_domain->get_max();
    else if(domain_type ==   RANGE_VAR) return range_domain->get_max();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_max();
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.get_max();
    else  return (*bool_domain >> 1);
  }

  int Mistral::Variable::get_initial_min() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->get_initial_min();
    else if(domain_type ==    LIST_VAR) return list_domain->get_initial_min();
    else if(domain_type ==   RANGE_VAR) return range_domain->get_initial_min();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_initial_min();
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.get_initial_min();
    else  return 0;
  }

  int Mistral::Variable::get_initial_max() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->get_initial_max();
    else if(domain_type ==    LIST_VAR) return list_domain->get_initial_max();
    else if(domain_type ==   RANGE_VAR) return range_domain->get_initial_max();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_initial_max();
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.get_initial_max();
    else  return 1;
  }

  int Mistral::Variable::get_min_pos() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->get_min_pos();
    else if(domain_type ==    LIST_VAR) return list_domain->get_min_pos();
    else if(domain_type ==   RANGE_VAR) return range_domain->get_min_pos();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_min_pos();
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.get_min_pos();
    else  return (*bool_domain >> 1); //(!(*bool_domain & 1));
  }

  int Mistral::Variable::get_max_neg() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->get_max_neg();
    else if(domain_type ==    LIST_VAR) return list_domain->get_max_neg();
    else if(domain_type ==   RANGE_VAR) return range_domain->get_max_neg();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_max_neg();
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.get_max_neg();
    else  return (!(*bool_domain & 1));
  }

  int Mistral::Variable::next(const int v) const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->next(v);
    else if(domain_type ==    LIST_VAR) return list_domain->next(v);
    else if(domain_type ==   RANGE_VAR) return range_domain->next(v);
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->next(v);
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.next(v);
    else  return (*bool_domain >> 1);
  }

  int Mistral::Variable::prev(const int v) const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->prev(v);
    else if(domain_type ==    LIST_VAR) return list_domain->prev(v);
    else if(domain_type ==   RANGE_VAR) return range_domain->prev(v);
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->prev(v);
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.prev(v);
    else  return (!(*bool_domain & 1));
  }

  bool Mistral::Variable::is_range() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->is_range();
    else if(domain_type ==    LIST_VAR) return list_domain->is_range();
    else if(domain_type ==   RANGE_VAR) return range_domain->is_range();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->is_range();
    else if(domain_type ==   CONST_VAR) return true;
    else if(domain_type ==   EXPRESSION) return expression->self.is_range();
    else return true;
  }

// bool Mistral::Variable::is_ground(const Expression *x) const {
//   return x->self.is_ground();
// }
  bool Mistral::Variable::is_ground() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->is_ground();
    else if(domain_type ==    LIST_VAR) return list_domain->is_ground();
    else if(domain_type ==   RANGE_VAR) return range_domain->is_ground();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->is_ground();
    else if(domain_type ==   CONST_VAR) return true;
    else if(domain_type ==   EXPRESSION) return expression->self.is_ground();
    else  return (*bool_domain != 3);
  }

  bool Mistral::Variable::equal(const int v) const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->equal(v);
    else if(domain_type ==    LIST_VAR) return list_domain->equal(v);
    else if(domain_type ==   RANGE_VAR) return range_domain->equal(v);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->equal(v);
    else if(domain_type ==   CONST_VAR) return (constant_value == v);
    else if(domain_type ==   EXPRESSION) return expression->self.equal(v);
    else  return (*bool_domain-1 == v);
  }

  bool Mistral::Variable::contain(const int v) const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->contain(v);
    else if(domain_type ==    LIST_VAR) return list_domain->contain(v);
    else if(domain_type ==   RANGE_VAR) return range_domain->contain(v);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->contain(v);
    else if(domain_type ==   CONST_VAR) return (constant_value == v);
    else if(domain_type ==   EXPRESSION) return expression->self.contain(v);
    else  return (!(v >> 1) && (*bool_domain & (v+1)));
  }

  bool Mistral::Variable::intersect(const int lo, const int up) const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->intersect(lo, up);
    else if(domain_type ==    LIST_VAR) return list_domain->intersect(lo, up);
    else if(domain_type ==   RANGE_VAR) return range_domain->intersect(lo, up);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->intersect(lo, up);
    else if(domain_type ==   CONST_VAR) return (constant_value >= lo && constant_value <= up);
    else if(domain_type ==   EXPRESSION) return expression->self.intersect(lo, up);
    else  return (((lo<=0) | (2*(up>0))) & *bool_domain);
  }

  bool Mistral::Variable::included(const int lo, const int up) const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->included(lo, up);
    else if(domain_type ==    LIST_VAR) return list_domain->included(lo, up);
    else if(domain_type ==   RANGE_VAR) return range_domain->included(lo, up);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->included(lo, up);
    else if(domain_type ==   CONST_VAR) return (constant_value >= lo && constant_value <= up);
    else if(domain_type ==   EXPRESSION) return expression->self.included(lo, up);
    else  {
      int state = *bool_domain;
      return ( up >= (state >> 1) && (lo <= !(state & 1)) );
    }
  }

  bool Mistral::Variable::includes(const int lo, const int up) const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->includes(lo, up);
    else if(domain_type ==    LIST_VAR) return list_domain->includes(lo, up);
    else if(domain_type ==   RANGE_VAR) return range_domain->includes(lo, up);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->includes(lo, up);
    else if(domain_type ==   CONST_VAR) return (constant_value == lo && constant_value == up);
    else if(domain_type ==   EXPRESSION) return expression->self.includes(lo, up);
    else  {
      int state = *bool_domain;
      return ( up <= (state >> 1) && (lo >= !(state & 1)) );
    }
  }

  bool Mistral::Variable::intersect(const BitSet& s) const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->intersect(s);
    else if(domain_type ==    LIST_VAR) return list_domain->intersect(s);
    else if(domain_type ==   RANGE_VAR) return range_domain->intersect(s);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->intersect(s);
    else if(domain_type ==   CONST_VAR) return (s.contain(constant_value));
    else if(domain_type ==   EXPRESSION) return expression->self.intersect(s);
    else  return s.intersect(*bool_domain);
  }

  bool Mistral::Variable::included(const BitSet& s) const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->included(s);
    else if(domain_type ==    LIST_VAR) return list_domain->included(s);
    else if(domain_type ==   RANGE_VAR) return range_domain->included(s);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->included(s);
    else if(domain_type ==   CONST_VAR) return (s.contain(constant_value));
    else if(domain_type ==   EXPRESSION) return expression->self.included(s);
    else  return s.includes(*bool_domain);
  }

  bool Mistral::Variable::includes(const BitSet& s) const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->includes(s);
    else if(domain_type ==    LIST_VAR) return list_domain->includes(s);
    else if(domain_type ==   RANGE_VAR) return range_domain->includes(s);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->includes(s);
    else if(domain_type ==   CONST_VAR) return (s.size() == 1 && s.contain(constant_value));
    else if(domain_type ==   EXPRESSION) return expression->self.includes(s);
    else  return s.included(*bool_domain);
  }

  bool Mistral::Variable::intersect(const Mistral::Variable& x) const {
    if(is_ground()) return x.contain(get_min());
    else if(is_range()) return x.intersect(get_min(), get_max());
    else if(domain_type ==  BITSET_VAR) return x.intersect(bitset_domain->domain.values);
    std::cout << "TODO!" << std::endl;
    return true;
  }

  bool Mistral::Variable::included(const Mistral::Variable& x) const {
    if(is_ground()) return x.contain(get_min());
    else if(x.is_ground()) return equal(x.get_min());
    else if(is_range()) return x.includes(get_min(), get_max());
    else if(x.is_range()) return included(x.get_min(), x.get_max()); 
    else if(domain_type ==  BITSET_VAR) return x.includes(bitset_domain->domain.values);
    std::cout << "TODO!" << std::endl;
    return true;
  }

  bool Mistral::Variable::includes(const Mistral::Variable& x) const {
    if(is_ground()) return x.equal(get_min());
    else if(x.is_ground()) return contain(x.get_min());
    else if(is_range()) return x.included(get_min(), get_max());
    else if(x.is_range()) return includes(x.get_min(), x.get_max()); 
    else if(domain_type ==  BITSET_VAR) return x.included(bitset_domain->domain.values);
    std::cout << "TODO!" << std::endl;
    return true;
  }

  void Mistral::Variable::intersect_to( BitSet& s ) const {
    if     (domain_type ==  BITSET_VAR) bitset_domain->intersect_to(s);
    else if(domain_type ==    LIST_VAR) list_domain->intersect_to(s);
    else if(domain_type ==   RANGE_VAR) range_domain->intersect_to(s);
    else if(domain_type == VIRTUAL_VAR) virtual_domain->intersect_to(s);
    else if(domain_type ==   CONST_VAR) {
      if(s.contain(constant_value)) {
	s.clear();
	s.add(constant_value);
      } else s.clear();
    }
    else if(domain_type ==   EXPRESSION) return expression->self.intersect_to(s);
    else  return s.intersect_with(*bool_domain);
  }

  void Mistral::Variable::union_to( BitSet& s ) const {
    if     (domain_type ==  BITSET_VAR) bitset_domain->union_to(s);
    else if(domain_type ==    LIST_VAR) list_domain->union_to(s);
    else if(domain_type ==   RANGE_VAR) range_domain->union_to(s);
    else if(domain_type == VIRTUAL_VAR) virtual_domain->union_to(s);
    else if(domain_type ==   CONST_VAR) s.add(constant_value);
    else if(domain_type ==   EXPRESSION) expression->self.union_to(s);
    else s.union_with(*bool_domain);
  }

  Mistral::Event Mistral::Variable::remove(const int v) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->remove(v);
    else if(domain_type ==    LIST_VAR) return list_domain->remove(v);
    else if(domain_type ==   RANGE_VAR) return range_domain->remove(v);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->remove(v);
    else if(domain_type ==   CONST_VAR) return (constant_value == v ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.remove(v);
    else {

      return (v<0||v>1 ? NO_EVENT : setValue(2-v));
      
      //return ((v<0 || v>1) ? NO_EVENT : setValue(2-v));

      //int dom = *bool_domain;


      //std::cout << dom << " " << v << " " << (dom&(v+1)) << 
      //return ((v < (!(dom&1)) || v > (dom>>1)) ? NO_EVENT : setValue(2-v));

      //return (((dom&(v+1)) == (v+1)) ? setValue(2-v) : NO_EVENT);
    }
  }

  Mistral::Event Mistral::Variable::set_domain(const int v) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->set_domain(v);
    else if(domain_type ==    LIST_VAR) return list_domain->set_domain(v);
    else if(domain_type ==   RANGE_VAR) return range_domain->set_domain(v);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->set_domain(v);
    else if(domain_type ==   CONST_VAR) return (constant_value != v ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.set_domain(v);
    else {

      return (v<0||v>1 ? FAIL_EVENT : setValue(1+v));
      //int dom = *bool_domain;
      //return ((dom==1+v) ? NO_EVENT : setValue(1+v));
    }
  }

  Mistral::Event Mistral::Variable::set_min(const int lo) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->set_min(lo);
    else if(domain_type ==    LIST_VAR) return list_domain->set_min(lo);
    else if(domain_type ==   RANGE_VAR) return range_domain->set_min(lo);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->set_min(lo);
    else if(domain_type ==   CONST_VAR) return (constant_value < lo ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.set_min(lo);
    else {

      // 1 [1, inf[
      // 2 [2, inf[
      // 3 [2, inf[
      return(lo==1 ? setValue(2) : (lo>1 ? FAIL_EVENT : NO_EVENT));

//       int dom = *bool_domain;
//       return (lo<=(!(dom&1)) ? NO_EVENT : (lo>(dom>>1) ? FAIL_EVENT : setValue(2)));
//       //return (lo<1 ? NO_EVENT : (lo>1 ? FAIL_EVENT : setValue(2)));
    }
  }

  Mistral::Event Mistral::Variable::set_max(const int up) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->set_max(up);
    else if(domain_type ==    LIST_VAR) return list_domain->set_max(up);
    else if(domain_type ==   RANGE_VAR) return range_domain->set_max(up);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->set_max(up);
    else if(domain_type ==   CONST_VAR) return (constant_value > up ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.set_max(up);
    else {

      return(up==0 ? setValue(1) : (up<0 ? FAIL_EVENT : NO_EVENT));

//       int dom = *bool_domain;
//       return (up>=(dom>>1) ? NO_EVENT : (up<(dom&1) ? FAIL_EVENT : setValue(1)));
      //return (up>0 ? NO_EVENT : (up<0 ? FAIL_EVENT : setValue(1)));
    }
  }

  Mistral::Event Mistral::Variable::set_domain(const BitSet& s) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->set_domain(s);
    else if(domain_type ==    LIST_VAR) return list_domain->set_domain(s);
    else if(domain_type ==   RANGE_VAR) return range_domain->set_domain(s);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->set_domain(s);
    else if(domain_type ==   CONST_VAR) return (s.contain(constant_value) ? NO_EVENT : FAIL_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.set_domain(s);
    else {
      //return ((s.pos_words<1 || s.neg_words>0) ? FAIL_EVENT : ((s.table[0]&3)==3 ? NO_EVENT : setValue(s.table[0])));
      return ((s.pos_words<1 || s.neg_words>0) ? FAIL_EVENT : setState(s.table[0]&*bool_domain));
    }
  }

  Mistral::Event Mistral::Variable::set_domain(Mistral::Variable& x) {


    //if(x.is_ground()) return set_domain(x.get_min());
    if(x.is_ground()) {
      return  set_domain(x.get_min());
    }
    else if(x.is_range()) {
      Event evt = (set_min(x.get_min()) | set_max(x.get_max()));
      return evt;
    }
    else if(x.domain_type ==  BITSET_VAR) return set_domain(x.bitset_domain->domain.values);
    else if(x.domain_type ==  EXPRESSION) return set_domain(x.expression->self);
    else {
      std::cout << "TODO! (set_domain(var))" << std::endl;
      exit(1);
    }
    return NO_EVENT;
  }

  Mistral::Event Mistral::Variable::removeSet(const BitSet& s) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->removeSet(s);
    else if(domain_type ==    LIST_VAR) return list_domain->removeSet(s);
    else if(domain_type ==   RANGE_VAR) return range_domain->removeSet(s);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->removeSet(s);
    else if(domain_type ==   CONST_VAR) return (s.contain(constant_value) ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.removeSet(s);
    else {

      return ((s.pos_words<1 || s.neg_words>0 || (s.table[0]^3)==3) ? NO_EVENT : setValue(s.table[0]^3));
    }
  }

  Mistral::Event Mistral::Variable::remove_interval(const int lo, const int up) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->remove_interval(lo, up);
    else if(domain_type ==    LIST_VAR) return list_domain->remove_interval(lo, up);
    else if(domain_type ==   RANGE_VAR) return range_domain->remove_interval(lo, up);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->remove_interval(lo, up);
    else if(domain_type ==   CONST_VAR) return ((constant_value < lo || constant_value > up) ? NO_EVENT : FAIL_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.remove_interval(lo, up);
    else {

      return (lo==1 ? setValue(1) : (up==0 ? setValue(2) : ((lo>1 || up<0) ? NO_EVENT : FAIL_EVENT)));
    }
  }

  Mistral::Event Mistral::Variable::restore() {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->restore();
    //else if(domain_type ==    LIST_VAR) return list_domain->restore();
    else if(domain_type ==   RANGE_VAR) return range_domain->restore();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->restore();
    else if(domain_type ==   CONST_VAR) return NO_EVENT;
    else if(domain_type ==   EXPRESSION) {
      //std::cout << "RESTORE EXPRESSION" << std::endl;
      //exit(1);
      return expression->self.restore();
    }
    else {
      *bool_domain = 3;
      return NO_EVENT;
    } 
  }


// bool Mistral::VariableImplementation::is_new(Solver *s) {
//   return (solver != s);
// }

// void Mistral::VariableImplementation::initialise(Solver *s) {
//   solver = s;
//   id = s->declare(*this);
// }

int Mistral::VariableImplementation::get_solution_value() const { 
  return solver->last_solution_lb[id] ;
}  
int Mistral::VariableImplementation::get_solution_min() const { 
  return solver->last_solution_lb[id] ; 
} 
int Mistral::VariableImplementation::get_solution_max() const { 
  return solver->last_solution_ub[id] ; 
}  


//void Mistral::VariableImplementation::trigger_value_event_and_save(Mistral::Variable *x) {
void Mistral::VariableImplementation::trigger_value_event_and_save() {
  solver->trigger_event(id, VALUE_EVENT);
  //solver->save(*x);
  //solver->save(this, BOOL_VAR);
  solver->save(id);
}

Mistral::BitsetDomain::BitsetDomain(const int lb, const int ub) {
  min = lb;
  max = ub;
  size = ub-lb+1;
  //values.initialise();
}

void Mistral::BitsetDomain::initialise(const int lb, const int ub, const bool vals) {
  min = lb;
  max = ub;
  size = ub-lb+1;
  if(vals) values.initialise(lb, ub, BitSet::full);
}

Mistral::BitsetDomain::BitsetDomain(Vector< int >& vals) {
  initialise(vals);
}

void Mistral::BitsetDomain::initialise(Vector< int >& vals) {
  min = vals.front();
  max = vals.back();
  size = vals.size;
  values.initialise(vals);
}

bool Mistral::Decision::propagateRelation() {
  //return !((Constraint*)(var.implementation))->propagate();
  return !var.constraint->propagate();
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Variable& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Variable* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Decision& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Decision* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::BitsetDomain& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::BitsetDomain* x) {
  return x->display(os);
}

Mistral::Expression::Expression(const int lo, const int up) 
  : VariableImplementation() {
  id=-1; 
  Variable x(lo, up, DYN_VAR);
  self = x;
}

Mistral::Expression::Expression(Vector< int >& values) 
  : VariableImplementation() {
  id=-1; 
  Variable x(values, DYN_VAR);
  self = x;
}

Mistral::Expression::Expression(Vector< Variable >& args) 
  : VariableImplementation() {
  id=-1; 
  for(unsigned int i=0; i<args.size; ++i)
    children.add(args[i]);
}
Mistral::Expression::Expression(Vector< Variable >& args, const int lo, const int up) 
  : VariableImplementation() {
  id=-1; 
  Variable x(lo, up, DYN_VAR);
  self = x;
  for(unsigned int i=0; i<args.size; ++i)
    children.add(args[i]);
}
Mistral::Expression::Expression(std::vector< Variable >& args) 
  : VariableImplementation() {
  id=-1; 
  for(unsigned int i=0; i<args.size(); ++i)
    children.add(args[i]);
}
Mistral::Expression::Expression(std::vector< Variable >& args, const int lo, const int up) 
  : VariableImplementation() {
  id=-1; 
  Variable x(lo, up, DYN_VAR);
  self = x;
  for(unsigned int i=0; i<args.size(); ++i)
    children.add(args[i]);
}
Mistral::Expression::Expression(Variable X, Variable Y) 
  : VariableImplementation() {
     children.add(X);
     children.add(Y);
  }
Mistral::Expression::Expression(Variable X) 
  : VariableImplementation() {
     children.add(X);
  }
Mistral::Expression::~Expression() {

  // std::cout << "delete exp subvar: " << self << std::endl;

  // int domain_type = self.domain_type;
  // if     (domain_type ==  BITSET_VAR) delete self.bitset_domain;
  // else if(domain_type ==    LIST_VAR) delete self.list_domain;
  // else if(domain_type ==   RANGE_VAR) delete self.range_domain;
  // else if(domain_type == VIRTUAL_VAR) delete self.virtual_domain;
  // else if(domain_type ==  EXPRESSION) delete self.expression;
  // else if(domain_type !=   CONST_VAR) delete self.variable;
}

void Mistral::Expression::extract_variable(Solver *s) {
  self.initialise(s, false);
  self = self.get_var();
}

// void Mistral::Expression::extract_variable(Solver *s) {
//   self.initialise(s, false);
//   self = self.get_var();
// }

std::ostream& Mistral::Expression::display(std::ostream& os) const {
  if(is_initialised())
    os << self << ":";
  os << get_name() << "(" ;
  if(children.empty()) os << self;
  else {
    os << children[0];
    for(unsigned int i=1; i<children.size-is_initialised(); ++i) {
      os << ", " << children[i];
    }
  }
  os << ")";
  return os;
}
//   Mistral::BinaryExpression::BinaryExpression(Variable X, Variable Y) 
//     : Expression() {
//     children.add(X);
//     children.add(Y);
//   }
//   Mistral::BinaryExpression::~BinaryExpression() {}


Mistral::AddExpression::AddExpression(Variable X, Variable Y) 
  : Expression(X,Y) {
}
Mistral::AddExpression::~AddExpression() {}

void Mistral::AddExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Add predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::AddExpression::extract_variable(Solver *s) {
//void Mistral::AddExpression::reify(Solver *s, Variable X) {
    //int lb = children[0].get_var().get_min()+children[1].get_var().get_min();
    //int ub = children[0].get_var().get_max()+children[1].get_var().get_max();

    int lb = children[0].get_min()+children[1].get_min();
    int ub = children[0].get_max()+children[1].get_max();

    Variable aux(lb, ub, DYN_VAR);
    self = aux;

    self.initialise(s, false);
    self = self.get_var();
    children.add(self);
  }

const char* Mistral::AddExpression::get_name() const {
  return "add";
}

  void Mistral::AddExpression::extract_predicate(Solver *s) {
    s->add(new PredicateAdd(children));
  }


Mistral::OffsetExpression::OffsetExpression(Variable X, const int ofs) 
  : Expression(X) { 
  offset=ofs; 
}
Mistral::OffsetExpression::~OffsetExpression() {}
  
void Mistral::OffsetExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Offset predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::OffsetExpression::extract_variable(Solver *s) {
//void Mistral::OffsetExpression::reify(Solver *s, Variable X) {
  int lb = children[0].get_min()+offset;
  int ub = children[0].get_max()+offset;

    Variable aux(lb, ub, DYN_VAR);
    self = aux;

  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
}

const char* Mistral::OffsetExpression::get_name() const {
  return "offset";
}

void Mistral::OffsetExpression::extract_predicate(Solver *s) {
  s->add(new PredicateOffset(children, offset));
}

Mistral::Variable Mistral::Variable::operator+(Variable x) {
  Variable exp(new AddExpression(*this,x));
  return exp;
}

Mistral::Variable Mistral::Variable::operator+(int k) {
  Variable exp(new OffsetExpression(*this,k));
  return exp;
}


Mistral::MulExpression::MulExpression(Variable X, Variable Y) 
  : Expression(X,Y) {
}
Mistral::MulExpression::~MulExpression() {}

void Mistral::MulExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Mul predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::MulExpression::extract_variable(Solver *s) {

  int z[4], lb = INFTY, ub = -INFTY, i = 4;
  z[0] = children[0].get_min() * children[1].get_min();
  z[1] = children[0].get_min() * children[1].get_max();
  z[2] = children[0].get_max() * children[1].get_min();
  z[3] = children[0].get_max() * children[1].get_max();  
  while( i-- ) {
    if( z[i] > ub ) ub = z[i];
    if( z[i] < lb ) lb = z[i];
  }

  Variable aux(lb, ub, DYN_VAR);
  self = aux;

  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
}

const char* Mistral::MulExpression::get_name() const {
  return "mul";
}

void Mistral::MulExpression::extract_predicate(Solver *s) {
  s->add(new PredicateMul(children));
}


Mistral::FactorExpression::FactorExpression(Variable X, const int fct) 
  : Expression(X) { 
  factor=fct; 
}
Mistral::FactorExpression::~FactorExpression() {}
  
void Mistral::FactorExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Factor predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::FactorExpression::extract_variable(Solver *s) {

  int lb = (factor<0 ? children[0].get_max() : children[0].get_min())*factor;
  int ub = (factor<0 ? children[0].get_min() : children[0].get_max())*factor;

  Variable aux(lb, ub, DYN_VAR);
  self = aux;
  
  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
}

const char* Mistral::FactorExpression::get_name() const {
  return "factor";
}

void Mistral::FactorExpression::extract_predicate(Solver *s) {
  s->add(new PredicateFactor(children, factor));
}

Mistral::Variable Mistral::Variable::operator*(Variable x) {
  Variable exp(new MulExpression(*this,x));
  return exp;
}

Mistral::Variable Mistral::Variable::operator*(int k) {
  Variable exp(new FactorExpression(*this,k));
  return exp;
}


  Mistral::SubExpression::SubExpression(Variable X, Variable Y) 
  : Expression(X,Y) {
}
  Mistral::SubExpression::~SubExpression() {}
  
void Mistral::SubExpression::extract_constraint(Solver *s) {
      std::cerr << "Error: Sub predicate can't be used as a constraint" << std::endl;
  exit(0);
  }

void Mistral::SubExpression::extract_variable(Solver *s) {
  //void Mistral::SubExpression::reify(Solver *s, Variable X) {
  //int lb = children[0].get_var().get_min()-children[1].get_var().get_max();
  //int ub = children[0].get_var().get_max()-children[1].get_var().get_min();
  
  int lb = children[0].get_min()-children[1].get_max();
  int ub = children[0].get_max()-children[1].get_min();
  
    Variable aux(lb, ub, DYN_VAR);
    self = aux;

  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
 
}

  void Mistral::SubExpression::extract_predicate(Solver *s) {
    Constraint *sub = new PredicateAdd();
    for(int i=3; i;) sub->scope.add(children[--i]// .get_var()
				    );
    //sub->initialise();
    //return sub;
    s->add(sub);
  }

const char* Mistral::SubExpression::get_name() const {
  return "sub";
}

Mistral::Variable Mistral::Variable::operator-(Variable x) {
  Variable exp(new SubExpression(*this,x));
  return exp;
}


Mistral::NotExpression::NotExpression(Variable X) 
  : Expression(X) { 
}
Mistral::NotExpression::~NotExpression() {}
  
void Mistral::NotExpression::extract_constraint(Solver *s) {
  children[0].remove(0);
  // std::cerr << "Error: Not predicate can't be used as a constraint" << std::endl;
  // exit(0);
}

void Mistral::NotExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  self = aux;

  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
}

const char* Mistral::NotExpression::get_name() const {
  return "not";
}

void Mistral::NotExpression::extract_predicate(Solver *s) {
  s->add(new PredicateNot(children));
}

Mistral::Variable Mistral::Variable::operator!() {
  Variable exp(new NotExpression(*this));
  return exp;
}

Mistral::AndExpression::AndExpression(Variable X, Variable Y) 
  : Expression(X,Y) {}
Mistral::AndExpression::~AndExpression() {}
  
void Mistral::AndExpression::extract_constraint(Solver *s) {
  //s->add(new ConstraintAnd(children));
  //std::cout << "not implemented" << std::endl;
  //exit(1);
  children[0].remove(0);
  children[1].remove(0);
}

void Mistral::AndExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  self = aux;

  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
}

void Mistral::AndExpression::extract_predicate(Solver *s) {
  s->add(new PredicateAnd(children));
}

const char* Mistral::AndExpression::get_name() const {
  return "and";
}

Mistral::Variable Mistral::Variable::operator&&(Variable x) {
  Variable exp(new AndExpression(*this,x));
  return exp;
}



Mistral::OrExpression::OrExpression(Variable X, Variable Y) 
  : Expression(X,Y) {}
Mistral::OrExpression::~OrExpression() {}
  
void Mistral::OrExpression::extract_constraint(Solver *s) {
  s->add(new ConstraintOr(children));
  //s->add(ConstraintOr::ConstraintOr_new(children));
  //std::cout << "not implemented" << std::endl;
  //exit(1);
}

void Mistral::OrExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  self = aux;
  
  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
}

void Mistral::OrExpression::extract_predicate(Solver *s) {
  s->add(new PredicateOr(children));
}

const char* Mistral::OrExpression::get_name() const {
  return "or";
}

Mistral::Variable Mistral::Variable::operator||(Variable x) {
  Variable exp(new OrExpression(*this,x));
  return exp;
}




//   Mistral::NeqExpression::NeqExpression(Variable X, Variable Y) 
//   : BinaryExpression(X,Y) {}
//   Mistral::NeqExpression::~NeqExpression() {}
  
// void Mistral::NeqExpression::extract_constraint(Solver *s) {
// //     Constraint *neq = new ConstraintNotEqual(children);
// //     neq->initialise();
// //     return neq;
//   s->add(new ConstraintNotEqual(children));
//   }

//   void Mistral::NeqExpression::extract_variable(Solver *s) {
//     Variable aux(0, 1, BOOL_VAR);
//     self = aux;

//     children.add(self);
//     self.initialise(s, false);
//   }

//   void Mistral::NeqExpression::extract_predicate(Solver *s) {
// //     Constraint *neq = new PredicateEqual(children, false);
// //     neq->initialise();
// //     return neq;
//     s->add(new PredicateEqual(children, false));
//   }

// const char* Mistral::NeqExpression::get_name() const {
//   return "neq";
// }



Mistral::EqualExpression::EqualExpression(Variable X, Variable Y, const int sp) 
  : Expression(X,Y) { spin=sp; value=NOVAL; }
Mistral::EqualExpression::EqualExpression(Variable X, const int y, const int sp) 
  : Expression() { children.add(X); value=y; spin=sp; }
Mistral::EqualExpression::~EqualExpression() {}

void Mistral::EqualExpression::extract_constraint(Solver *s) {
  if(spin) {
    if(children.size==2) s->add(new ConstraintEqual(children));
    else children[0].set_domain(value);
  } else {
    if(children.size==2) s->add(new ConstraintNotEqual(children));
    else children[0].remove(value);
  }
}

void Mistral::EqualExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  self = aux;
  
  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
}

void Mistral::EqualExpression::extract_predicate(Solver *s) {

  // if(children.size == 3)
  //   std::cout << children[0] << " in " << children[0].get_domain() << std::endl
  // 	      << children[1] << " in " << children[1].get_domain() << std::endl
  // 	      << children[2] << " in " << children[2].get_domain() << std::endl;
  // else if(children.size == 2)
  //   std::cout << children[0] << " in " << children[0].get_domain() << std::endl
  // 	      << children[1] << " in " << children[1].get_domain() << std::endl;
  // else 
  //   std::cout << "??" << std::endl;

  if(children.size==3) s->add(new PredicateEqual(children, spin));
  else s->add(new PredicateConstantEqual(children, value, spin));
}

const char* Mistral::EqualExpression::get_name() const {
  return "equal";
}

Mistral::Variable Mistral::Variable::operator==(Variable x) {
  Variable exp(new EqualExpression(*this,x,1));
  return exp;
}

Mistral::Variable Mistral::Variable::operator!=(Variable x) {
  Variable exp(new EqualExpression(*this,x,0));
  return exp;
}

Mistral::Variable Mistral::Variable::operator==(const int x) {
  Variable exp(new EqualExpression(*this,x,1));
  return exp;
}

Mistral::Variable Mistral::Variable::operator!=(const int x) {
  Variable exp(new EqualExpression(*this,x,0));
  return exp;
}


// Mistral::Variable Mistral::Variable::operator==(Variable x) {
//   Variable exp(new EqualExpression(*this,x));
//   return exp;
// }


Mistral::PrecedenceExpression::PrecedenceExpression(Variable X, Variable Y, 
						    const int of, const int sp) 
  : Expression(X,Y) { spin = sp; offset = of; }
Mistral::PrecedenceExpression::PrecedenceExpression(Variable X,  
						    const int of, const int sp) 
  : Expression(X) { spin = sp; offset = of; }
Mistral::PrecedenceExpression::~PrecedenceExpression() {}
  
void Mistral::PrecedenceExpression::extract_constraint(Solver *s) {
  if(children.size==2) {
    s->add(new ConstraintLess(children, offset));
  }
  else {
    if(spin) {
      children[0].set_max(offset);
    } else {     
      children[0].set_min(offset);
    }
  }
}

void Mistral::PrecedenceExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  self = aux;
  
  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
}

void Mistral::PrecedenceExpression::extract_predicate(Solver *s) {
  if(children.size==3) s->add(new PredicateLess(children, offset));
  else if(spin) {
    s->add(new PredicateUpperBound(children, offset));
  } else {
    s->add(new PredicateLowerBound(children, offset));
  }
}

const char* Mistral::PrecedenceExpression::get_name() const {
  return "prec";
}

Mistral::Variable Mistral::Precedence(Variable X, const int d, Variable Y) 
{
  Variable exp(new PrecedenceExpression(X,Y,d));
  return exp;
}


Mistral::Variable Mistral::Variable::operator<(Variable x) {
  Variable exp(new PrecedenceExpression(*this,x,1));
  return exp;
}

Mistral::Variable Mistral::Variable::operator>(Variable x) {
  Variable exp(new PrecedenceExpression(x,*this,1));
  return exp;
}

Mistral::Variable Mistral::Variable::operator<=(Variable x) {
  Variable exp(new PrecedenceExpression(*this,x));
  return exp;
}

Mistral::Variable Mistral::Variable::operator>=(Variable x) {
  Variable exp(new PrecedenceExpression(x,*this));
  return exp;
}

Mistral::Variable Mistral::Variable::operator<(const int k) {
  Variable exp(new PrecedenceExpression(*this,k-1,1));
  return exp;
}

Mistral::Variable Mistral::Variable::operator>(const int k) {
  Variable exp(new PrecedenceExpression(*this,k+1,0));
  return exp;
}

Mistral::Variable Mistral::Variable::operator<=(const int k) {
  Variable exp(new PrecedenceExpression(*this,k,1));
  return exp;
}

Mistral::Variable Mistral::Variable::operator>=(const int k) {
  Variable exp(new PrecedenceExpression(*this,k,0));
  return exp;
}


Mistral::DisjunctiveExpression::DisjunctiveExpression(Variable X, Variable Y, 
						      const int px, const int py) 
  : Expression(X,Y) { processing_time[0] = px; processing_time[1] = py; }
Mistral::DisjunctiveExpression::~DisjunctiveExpression() {}
  
void Mistral::DisjunctiveExpression::extract_constraint(Solver *s) {
  s->add(new ConstraintDisjunctive(children, processing_time[0], processing_time[1]));
}

void Mistral::DisjunctiveExpression::extract_variable(Solver *s) {
  std::cerr << "Error: Disjunctive constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

void Mistral::DisjunctiveExpression::extract_predicate(Solver *s) {
  std::cerr << "Error: Disjunctive constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::DisjunctiveExpression::get_name() const {
  return "dsijunct";
}

Mistral::Variable Mistral::Disjunctive(Variable X, Variable Y, 
				       const int px, const int py) 
{
  Variable exp(new DisjunctiveExpression(X,Y,px,py));
  return exp;
}


Mistral::AllDiffExpression::AllDiffExpression(Vector< Variable >& args, const int ct) 
  : Expression(args) { consistency_level = ct; }

Mistral::AllDiffExpression::~AllDiffExpression() {}

void Mistral::AllDiffExpression::extract_constraint(Solver *s) { 
//   Constraint *con = new ConstraintAllDiff(children); 
//   con->initialise();
//   return con;
  if(consistency_level == BOUND_CONSISTENCY)
    s->add(new ConstraintAllDiff(children)); 
  s->add(new ConstraintCliqueNotEqual(children)); 
//   Vector< Variable > pair;
//   for(unsigned int i=0; i<children.size-1; ++i)
//     for(unsigned int j=i+1; j<children.size; ++j) {
//       pair.clear();
//       pair.add(children[i]);
//       pair.add(children[j]);
//       s->add(new ConstraintNotEqual(pair));
//     }
}

void Mistral::AllDiffExpression::extract_variable(Solver *s) {
      std::cerr << "Error: AllDiff constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

void Mistral::AllDiffExpression::extract_predicate(Solver *s) { 
      std::cerr << "Error: AllDiff constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::AllDiffExpression::get_name() const {
  return "alldiff";
}

Mistral::Variable Mistral::AllDiff(Vector< Variable >& args, const int ct) {
  Variable exp(new AllDiffExpression(args,ct));
  return exp;
}


Mistral::BoolSumExpression::BoolSumExpression(const int l, const int u) 
  : Expression() {
  lb = l;
  ub = u;
}

Mistral::BoolSumExpression::BoolSumExpression(Vector< Variable >& args, const int l, const int u) 
  : Expression(args) {
  lb = l;
  ub = u;
}

Mistral::BoolSumExpression::~BoolSumExpression() {}

void Mistral::BoolSumExpression::extract_constraint(Solver *s) { 
  s->add(new ConstraintBoolSumInterval(children,lb,ub)); 
}

void Mistral::BoolSumExpression::extract_variable(Solver *s) {
  Variable aux(lb, ub, DYN_VAR);
  self = aux;
  
  self.initialise(s, false);
  self = self.get_var();
  //children.add(self);
}

void Mistral::BoolSumExpression::extract_predicate(Solver *s) { 
  s->add(new PredicateBoolSum(children, self)); 
}

const char* Mistral::BoolSumExpression::get_name() const {
  return "bool_sum";
}

Mistral::Variable Mistral::BoolSum(Vector< Variable >& args, const int l, const int u) {
  Variable exp(new BoolSumExpression(args,l,u));
  return exp;
}

Mistral::Variable Mistral::BoolSum(Vector< Variable >& args) {
  Variable exp(new BoolSumExpression(args,0,args.size));
  return exp;
}


Mistral::LinearExpression::LinearExpression(Vector< Variable >& args, 
					    Vector< int >& wgts, 
					    const int l, const int u) 
  : Expression(args) {
  lower_bound = l;
  upper_bound = u;
  for(unsigned int i=0; i<wgts.size; ++i) {
    weight.add(wgts[i]);
  }
}

Mistral::LinearExpression::LinearExpression(std::vector< Variable >& args, 
					    std::vector< int >& wgts, 
					    const int l, const int u) 
  : Expression(args) {
  lower_bound = l;
  upper_bound = u;
  for(unsigned int i=0; i<wgts.size(); ++i) {
    weight.add(wgts[i]);
  }
}

// Mistral::LinearExpression::LinearExpression(Variable X, const int coef) 
//   : Expression(X) { 
//   weight.add(coef);
// }
Mistral::LinearExpression::~LinearExpression() {}
  
// -1 -1 -1: - x - y =  z: 
// -1 -1  1: - x - y = -z: z = y + x
// -1  1 -1: - x + y =  z: y = x + z
// -1  1  1: - x + y = -z: x = y + z
//  1 -1 -1: + x - y =  z: x = y + z
//  1 -1  1: + x - y = -z: y = x + z
//  1  1 -1: + x + y =  z: z = x + y
//  1  1  1: + x + y = -z: 

void Mistral::LinearExpression::extract_constraint(Solver *s) {
  // check if we can use an 'Add' or 'Sub' predicate
  int post_add = false;
  if(children.size == 3 && 
     std::abs(weight[0]) == 1 &&
     std::abs(weight[1]) == 1 &&
     std::abs(weight[2]) == 1
     ) {
    int i=0;
    for(; i<3; ++i)
      if(weight[i] != weight[(i+1)%3] && weight[(i+1)%3] == weight[(i+2)%3]) {
	Variable x = children[i];
	children[i] = children[2];
	children[2] = x;
	weight[0] = weight[1] = 1;
	weight[2] = -1;
	break;
      }
    if(i<3) {
      post_add = true;
      s->add( new PredicateAdd(children) );
    }  
  }
  
  if(!post_add)
    s->add(new PredicateWeightedSum(children, weight, lower_bound, upper_bound));
}


void Mistral::LinearExpression::initialise_bounds() {
  int tlb=0;
  int tub=0;

  int lb = children[0].get_min()*weight[0];
  int ub = children[0].get_max()*weight[0];
  
  if(lb < ub) {
    tlb += lb;
    tub += ub;
  } else {
    tlb += ub;
    tub += lb;
  }

  for(unsigned int i=1; i<children.size; ++i) {
    lb = children[i].get_min()*weight[i];
    ub = children[i].get_max()*weight[i];
  
    if(lb < ub) {
      tlb += lb;
      tub += ub;
    } else {
      tlb += ub;
      tub += lb;
    }
  }

  if(tlb > lower_bound) lower_bound = tlb;
  if(tub < upper_bound) upper_bound = tub;
}

void Mistral::LinearExpression::extract_variable(Solver *s) {
  initialise_bounds();

  Variable aux(lower_bound, upper_bound, DYN_VAR);
  self = aux;

  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
  weight.add(-1);
}

const char* Mistral::LinearExpression::get_name() const {
  return "linear_sum";
}

void Mistral::LinearExpression::extract_predicate(Solver *s) {
  s->add(new PredicateWeightedSum(children, weight, 0, 0));
}

Mistral::Variable Mistral::Sum(Vector< Variable >& args, Vector< int >& wgts, Variable T) {
  LinearExpression *lexpr = new LinearExpression(args,wgts,0,0);
  lexpr->children.add(T);
  lexpr->weight.add(-1);
  Variable exp(lexpr);
  return exp;
}
Mistral::Variable Mistral::Sum(std::vector< Variable >& args, std::vector< int >& wgts, Variable T) {
  LinearExpression *lexpr = new LinearExpression(args,wgts,0,0);
  lexpr->children.add(T);
  lexpr->weight.add(-1);
  Variable exp(lexpr);
  return exp;
}
Mistral::Variable Mistral::Sum(Vector< Variable >& args, Vector< int >& wgts, const int l, const int u) {
  Variable exp( new LinearExpression(args, wgts, l, u) );
  return exp;
}
Mistral::Variable Mistral::Sum(std::vector< Variable >& args, std::vector< int >& wgts, const int l, const int u) {
  Variable exp( new LinearExpression(args, wgts, l, u) );
  return exp;
}


Mistral::ElementExpression::ElementExpression(Vector< Variable >& args, 
					      Variable X, int ofs) 
  : Expression(), offset(ofs) {
  for(unsigned int i=0; i<args.size; ++i) children.add(args[i]);
  children.add(X);
}

Mistral::ElementExpression::~ElementExpression() {}


void Mistral::ElementExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Element predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::ElementExpression::initialise_domain() {

  lower_bound = INFTY;
  upper_bound = -INFTY;

  int arity = children.size-1;
  //Vector< int > values;

  BitSet domain;

  int i, nxt = children[arity].get_min();
  do {
    i = nxt-offset;
    
    if(i>=0 && i<arity) {
      if(children[i].get_min() < lower_bound) lower_bound = children[i].get_min();
      if(children[i].get_max() > upper_bound) upper_bound = children[i].get_max();
    }

    nxt = children[arity].next(i+offset);
  } while(i+offset<nxt);
  

  domain.initialise(lower_bound, upper_bound, BitSet::empt);
  
  nxt = children[arity].get_min();
  do {
    i = nxt-offset;
    
    if(i>=0 && i<arity) {
      children[i].union_to(domain);
    }

    nxt = children[arity].next(i+offset);
  } while(i+offset<nxt);


  nxt = lower_bound;
  do {
    i = nxt;
    values.add(i);
    nxt = domain.next(i);
  } while(i<nxt);
  
}

void Mistral::ElementExpression::extract_variable(Solver *s) {
  initialise_domain();

  Variable aux(values, DYN_VAR);
  self = aux;

  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
}

const char* Mistral::ElementExpression::get_name() const {
  return "element";
}

void Mistral::ElementExpression::extract_predicate(Solver *s) {
  s->add(new PredicateElement(children, offset));
}


Mistral::Variable Mistral::VarArray::operator[](Variable X) {
  Variable exp( new ElementExpression(*this, X, 0) );
  return exp;
}

Mistral::Variable Mistral::VarArray::operator[](const int X) {
  return stack_[X];
}

Mistral::Variable Mistral::Element(Vector<Variable> X, Variable selector, int offset) {
  Variable exp( new ElementExpression(X, selector, offset) );
  return exp;
}

Mistral::Variable Mistral::Element(VarArray X, Variable selector, int offset) {
  Variable exp( new ElementExpression(X, selector, offset) );
  return exp;
}


Mistral::SetExpression::SetExpression(const int lelt, const int uelt, 
				      const int clb, const int cub) 
  : BoolSumExpression(clb, cub) {
  elements.initialise(0, uelt-lelt+1);
  for(int elt=lelt; elt<=uelt; ++elt) {
    elements.add(elt);
    Variable x(0, 1, BOOL_VAR);
    children.add(x);
  }
}

Mistral::SetExpression::~SetExpression() {}

const char* Mistral::SetExpression::get_name() const {
  return "set";
}

int Mistral::SetExpression::get_element_index(const int vali)  {
  int target, lb = 0, ub = elements.size-1, idx = -1;
  while(lb <= ub) {

    //std::cout << vali << " [" << lb << "," << ub << "] -> " ; //<< std::endl;

    target = (lb+ub)/2;

    // std::cout << target << ": " << elements[target] << " " << idx << " -- " <<
    // 	      elements << std::endl;

    if(elements[target] == vali) { idx=target; break; }
    else if(elements[target] > vali) ub = target-1;
    else lb = target+1;
  }
  //  std::cout << " ==> " << idx << std::endl; 
  return idx;
}

Mistral::Variable Mistral::SetExpression::get_elt_var(const int vali) {
  // int target = 0, lb = -1, ub = elements.size;
  // while(lb < ub) {
  //   target = lb+ub/2;
  //   if(elements[target] == vali) break;
  //   else if(elements[target] > vali) ub = target;
  //   else lb = target;
  // }
  // return children[target];
  int target = get_element_index(vali);
  if(target>=0) return children[target];
  Variable x;
  return x;
}

std::ostream& Mistral::SetExpression::display(std::ostream& os) const {
  os << "{" ;
  if(children[0].get_min()) os << elements[0];
  for(unsigned int i=1; i<elements.size; ++i) {
    if(children[i].get_min()) os << ", " << elements[i];
  }
  os << "} <= S" << id << " <= {"; 
  if(children[0].get_max()) os << elements[0];
  for(unsigned int i=1; i<elements.size; ++i) {
    if(children[i].get_max()) os << ", " << elements[i];
  }
  os << "}";
  return os;
}

int Mistral::SetExpression::get_solution_value() const { 
  int i=elements.size;
  int t = 0;
  int m = 1;
  while(--i>=0) {
    if(solver->last_solution_lb[children[i].id()]) {
      t += m*elements[i];
      m *= 10;
    }
  }
  return t;
} 

Mistral::Variable Mistral::SetVariable(const int lelt, const int uelt, const int clb, const int cub) {
  int ub = cub;
  if(cub > uelt-lelt+1) ub = (uelt-lelt+1);
  SetExpression *se = new SetExpression(lelt, uelt, clb, ub);
  //std::cout << "create " << se << std::endl;
  Variable exp(se);
  return exp;
}

Mistral::Variable Mistral::Card(Variable S) { return S; }

Mistral::SubsetExpression::SubsetExpression(Variable X, Variable Y) 
  : Expression(X,Y) { 
}

Mistral::SubsetExpression::~SubsetExpression() {}

void Mistral::SubsetExpression::extract_constraint(Solver *s) { 
  unsigned int i = 0, j = 0;
  SetExpression *x = (SetExpression*)(children[0].expression);
  SetExpression *y = (SetExpression*)(children[1].expression);

  s->add(children[0] <= children[1]);

  while(i<x->elements.size && j<y->elements.size) {
    if(x->elements[i] == y->elements[j]) {
      s->add(x->children[i] <= y->children[j]);
      ++i;
      ++j;
    } else if(x->elements[i] <= y->elements[j]) {
      s->add(x->children[i] == 0);
      ++i;
    } else {
      ++j;
    }
  }
  while(i<x->elements.size) {
    s->add(x->children[i] == 0);
    ++i;
  }
}

void Mistral::SubsetExpression::extract_variable(Solver *s) {
  std::cerr << "Error: Subset constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

void Mistral::SubsetExpression::extract_predicate(Solver *s) { 
  std::cerr << "Error: Subset constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::SubsetExpression::get_name() const {
  return "subset";
}

Mistral::Variable Mistral::Subset(Variable X, Variable Y) {
  Variable exp(new SubsetExpression(X, Y));
  return exp;
}


Mistral::MemberExpression::MemberExpression(Variable X, Variable Y) 
  : Expression(X,Y) { 
}

Mistral::MemberExpression::~MemberExpression() {}

void Mistral::MemberExpression::extract_constraint(Solver *s) { 

  //unsigned int i = 0, j = 0;
  SetExpression *y = (SetExpression*)(children[1].expression);

  s->add(children[1] >= 1); // at least one element in the set

  int vali, vnxt = children[0].get_min(), idx;
  do {

    vali = vnxt;

    idx = y->get_element_index(vali);

    if(idx>=0) {
      s->add((children[0] == vali) <= y->children[idx]);
    } else {
      s->add(children[0] != vali);
    }
    vnxt = children[0].next(vali);


  } while(vali < vnxt);

}

void Mistral::MemberExpression::extract_variable(Solver *s) {
  std::cerr << "Error: Member constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

void Mistral::MemberExpression::extract_predicate(Solver *s) { 
  std::cerr << "Error: Member constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::MemberExpression::get_name() const {
  return "member";
}

Mistral::Variable Mistral::Member(Variable X, Variable Y) {
  Variable exp(new MemberExpression(X, Y));
  return exp;
}


Mistral::Goal::Goal(method t) : type(t) {
  lower_bound = 0;
  upper_bound = 0;
}

Mistral::Goal::Goal(method t, Variable X) : type(t) {
  objective = X;

  // std::cout << "OBJECTIVE=" << objective << " in " << objective.get_domain() << std::endl;

  lower_bound = objective.get_min()-1;
  upper_bound = objective.get_max()+1;
}

Mistral::Goal::~Goal() {}

bool Mistral::Goal::enforce() {
  if(type == MINIMIZATION) {
    return IS_FAIL(objective.set_max(upper_bound-1));
  } else if(type == MAXIMIZATION) {
    return IS_FAIL(objective.set_min(lower_bound+1));
  }
  return false;
}
    

Mistral::Outcome Mistral::Goal::notify_exhausted() {
  if(type == SATISFACTION)
    return UNSAT;
  return OPT;
}

Mistral::Outcome Mistral::Goal::notify_solution(Solver *solver) {
  if(type == MINIMIZATION) {

    //std::cout << "minimization algorithm: new solution (ub=" ;
    upper_bound = objective.get_min();
    //std::cout << upper_bound << ")" << std::endl;
    
    //std::cout << objective << " in " << objective.get_domain() << std::endl;
    if(!solver->level) lower_bound = upper_bound;
    if(upper_bound == lower_bound) return OPT;


    //std::cout << solver->level << " " << solver->decisions.size << std::endl;
    solver->branch_right();
    return UNKNOWN; //(upper_bound == lower_bound ? OPT : UNKNOWN);
  }
  else if(type == MAXIMIZATION) {

    lower_bound = objective.get_max();

    if(!solver->level) upper_bound = lower_bound;
    if(upper_bound == lower_bound) return OPT;


    //std::cout << solver->level << " " << solver->decisions.size << std::endl;
    solver->branch_right();
    return UNKNOWN; //(upper_bound == lower_bound ? OPT : UNKNOWN);

    // solver->branch_right();
    // return (upper_bound == lower_bound ? OPT : UNKNOWN);
  }
  else if(type == ALLSOLUTIONS) {
    //solver->store_solution();

    if(!solver->level) return OPT;
    solver->branch_right();
    return UNKNOWN;
  }

  //std::cout << "satisfaction algorithm: new solution" << std::endl;
  return SAT;
}




