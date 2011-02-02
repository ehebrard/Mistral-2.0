
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

}

Mistral::Variable::Variable(const int lo, const int up, const int type) {
    if(lo == up) {
      domain_type = CONST_VAR;
      constant_value = lo;
    } else if(type == EXPRESSION) {
      domain_type = EXPRESSION;
      expression = new Expression(lo, up);
    } else if((type & BOOL_VAR) && lo==0 && up==1) {
      //domain_type = BOOL_VAR;
      bool_domain = &BOOL_DOM;
      variable = new VariableImplementation();
    } else if(type & RANGE_VAR) {
      domain_type = RANGE_VAR;
      range_domain = new VariableRange(lo, up);
    } else {
      domain_type = BITSET_VAR;
#ifdef _BIT64
      int nwords = ((up-lo) / 64);
      if(nwords < 1) bitset_domain = new VariableWord<unsigned long long int, 1>(lo, up);
      else if(nwords < 2) bitset_domain = new VariableWord<unsigned long long int, 2>(lo, up);
      else if(nwords < 3) bitset_domain = new VariableWord<unsigned long long int, 3>(lo, up);
#else
      int nwords = ((up-lo) / 32);
      if(nwords < 1) bitset_domain = new VariableWord<unsigned int, 1>(lo, up);
      else if(nwords < 2) bitset_domain = new VariableWord<unsigned int, 2>(lo, up);
      else if(nwords < 3) bitset_domain = new VariableWord<unsigned int, 3>(lo, up);
      else if(nwords < 4) bitset_domain = new VariableWord<unsigned int, 4>(lo, up);
      else if(nwords < 5) bitset_domain = new VariableWord<unsigned int, 5>(lo, up);
#endif
      else bitset_domain = new VariableBitmap(lo, up);
    }
}



Mistral::Variable Mistral::Variable::get_var() {
 //  display(std::cout);
//   std::cout << ".get_var()" << std::endl;

  if(domain_type == EXPRESSION) {
    //   std::cout << "expression" << std::endl;
    return expression->self.get_var();
  } else if(domain_type == CONST_VAR || !variable->solver) {
    //    std::cout << "self" << std::endl;
    return *this;
  }
  //  std::cout << "solver's" << std::endl;
  return variable->solver->variables[variable->id];
    //*this;
}

std::ostream& Mistral::Variable::display(std::ostream& os) const {
  //os << domain_type << "|";
  if(domain_type == EXPRESSION) {
    //Expression *exp = (Expression*)implementation;
    if(expression->is_initialised())
      os << expression->self << ":";
      //os << expression->self;
    //else {
    os << expression->get_name() << "(" ;
    if(expression->children.empty()) os << expression->self;
    else {
      os << expression->children[0];
      for(unsigned int i=1; i<expression->children.size-expression->is_initialised(); ++i) {
	os << ", " << expression->children[i];
      }
    }
    os << ")";
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
    
//     std::cout << std::endl << "Add a new expression: " << this 
//     	      << (top ? " at top-level" : " nested")
// 	      << std::endl;
    
    //Expression *exp = (Expression*)implementation;
    if(!expression->is_initialised()) {

      for(unsigned int i=0; i<expression->children.size; ++i) {
	//std::cout << this << " add child" << std::endl; 
	expression->children[i].initialise(s, false);
      }

      if(top && !expression->children.empty()) {
	expression->extract_constraint(s);
	//std::cout << this << " extract constraint: " << s->constraints.back() << std::endl; 
      } else {
	
	//std::cout << this << " extract variable " << std::endl;
	expression->extract_variable(s);
	expression->id = expression->self.id();
	expression->solver = s;	
	expression->extract_predicate(s);//);
	//std::cout << "now: " << this << std::endl;
      }
    }

    //else std::cout << this << " is known " << std::endl;
    
  } else {

//      std::cout << std::endl << "declare a new variable: " << this 
//    		<< " in " << get_domain() << std::endl;
     

    if(domain_type != CONST_VAR && variable->solver != s) {

 
      //variable->id = 
      s->declare(this->get_var());
      //variable->solver = s;
//       if(domain_type == BOOL_VAR) {
// 	s->booleans.add(this);
//       }
//       if(domain_type == BOOL_VAR) { bool_domain = 
// // 	  //(int)(new int[1]);
// // 	  (int)(s->getNextBooleanSlot());
// // 	//std::cout << "NEW DOMAIN: " << *((int*)domain_type) << std::endl;
// 	  }

      s->sequence.declare(*this);

//       std::cout << "declared a new variable: " << this 
//    		<< " in " << get_domain() << std::endl;
    } 

    //else std::cout << this << " is known " << std::endl;
  }

  //  std::cout << "end initialise (" << id() << ")" << std::endl << std::endl ;
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


Mistral::Event Mistral::Variable::setDomain( const int vals ) 
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


std::string Mistral::Variable::get_domain() const {
  std::ostringstream buf;
  if     (domain_type ==  BITSET_VAR) buf << bitset_domain->domain;
  else if(domain_type ==    LIST_VAR) buf << list_domain->domain;
  else if(domain_type ==   RANGE_VAR) {
    //Mistral::VariableRange* r = (Mistral::VariableRange*)implementation;
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

  int Mistral::Variable::get_minCapacity() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->get_minCapacity();
    else if(domain_type ==    LIST_VAR) return list_domain->get_minCapacity();
    else if(domain_type ==   RANGE_VAR) return range_domain->get_minCapacity();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_minCapacity();
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.get_minCapacity();
    else  return 0;
  }

  int Mistral::Variable::get_maxCapacity() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->get_maxCapacity();
    else if(domain_type ==    LIST_VAR) return list_domain->get_maxCapacity();
    else if(domain_type ==   RANGE_VAR) return range_domain->get_maxCapacity();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_maxCapacity();
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.get_maxCapacity();
    else  return 1;
  }

  int Mistral::Variable::get_minPosAbs() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->get_minPosAbs();
    else if(domain_type ==    LIST_VAR) return list_domain->get_minPosAbs();
    else if(domain_type ==   RANGE_VAR) return range_domain->get_minPosAbs();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_minPosAbs();
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.get_minPosAbs();
    else  return (!(*bool_domain & 1));
  }

  int Mistral::Variable::get_minNegAbs() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->get_minNegAbs();
    else if(domain_type ==    LIST_VAR) return list_domain->get_minNegAbs();
    else if(domain_type ==   RANGE_VAR) return range_domain->get_minNegAbs();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_minNegAbs();
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->self.get_minNegAbs();
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

  void Mistral::Variable::intersectTo( BitSet& s ) const {
    if     (domain_type ==  BITSET_VAR) bitset_domain->intersectTo(s);
    else if(domain_type ==    LIST_VAR) list_domain->intersectTo(s);
    else if(domain_type ==   RANGE_VAR) range_domain->intersectTo(s);
    else if(domain_type == VIRTUAL_VAR) virtual_domain->intersectTo(s);
    else if(domain_type ==   CONST_VAR) {
      if(s.contain(constant_value)) {
	s.clear();
	s.add(constant_value);
      } else s.clear();
    }
    else if(domain_type ==   EXPRESSION) return expression->self.intersectTo(s);
    else  return s.intersectWith(*bool_domain);
  }

  void Mistral::Variable::unionTo( BitSet& s ) const {
    if     (domain_type ==  BITSET_VAR) bitset_domain->unionTo(s);
    else if(domain_type ==    LIST_VAR) list_domain->unionTo(s);
    else if(domain_type ==   RANGE_VAR) range_domain->unionTo(s);
    else if(domain_type == VIRTUAL_VAR) virtual_domain->unionTo(s);
    else if(domain_type ==   CONST_VAR) s.add(constant_value);
    else if(domain_type ==   EXPRESSION) expression->self.unionTo(s);
    else s.unionWith(*bool_domain);
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
      return ((s.pos_words<1 || s.neg_words>0) ? FAIL_EVENT : setDomain(s.table[0]&*bool_domain));
    }
  }

  Mistral::Event Mistral::Variable::set_domain(Mistral::Variable& x) {
    if(x.is_ground()) return set_domain(x.get_min());
    else if(x.is_range()) return (set_min(x.get_min()) | set_max(x.get_max()));
    else if(x.domain_type ==  BITSET_VAR) return set_domain(bitset_domain->domain.values);
    std::cout << "TODO!" << std::endl;
    return NO_EVENT;
  }

  Mistral::Event Mistral::Variable::removeSet(const BitSet& s) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->removeSet(s);
    else if(domain_type ==    LIST_VAR) return list_domain->removeSet(s);
    else if(domain_type ==   RANGE_VAR) return range_domain->removeSet(s);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->removeSet(s);
    else if(domain_type ==   CONST_VAR) return (s.contain(constant_value) ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.removeSet(s);
    else  return ((s.pos_words<1 || s.neg_words>0 || (s.table[0]^3)==3) ? NO_EVENT : setValue(s.table[0]^3));
  }

  Mistral::Event Mistral::Variable::removeRange(const int lo, const int up) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->removeRange(lo, up);
    else if(domain_type ==    LIST_VAR) return list_domain->removeRange(lo, up);
    else if(domain_type ==   RANGE_VAR) return range_domain->removeRange(lo, up);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->removeRange(lo, up);
    else if(domain_type ==   CONST_VAR) return ((constant_value < lo || constant_value > up) ? NO_EVENT : FAIL_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.removeRange(lo, up);
    else  return (lo==1 ? setValue(1) : (up==0 ? setValue(2) : ((lo>1 || up<0) ? NO_EVENT : FAIL_EVENT)));
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
    int Mistral::VariableImplementation::get_solution_lb() const { 
      return solver->last_solution_lb[id] ; 
    } 
    int Mistral::VariableImplementation::get_solution_ub() const { 
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

void Mistral::Expression::extract_variable(Solver *s) {
  self.initialise(s, false);
  self = self.get_var();
}

Mistral::Expression::Expression(Vector< Variable >& args) 
  : VariableImplementation() {
    id=-1; 
    for(unsigned int i=0; i<args.size; ++i)
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
Mistral::Expression::~Expression() {}

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
    //int lb = children[0].get_var().get_min()+children[1].get_var().get_min();
    //int ub = children[0].get_var().get_max()+children[1].get_var().get_max();

    int lb = children[0].get_min()+children[1].get_min();
    int ub = children[0].get_max()+children[1].get_max();

    Variable aux(lb, ub, DYN_VAR
		 );
    self = aux;

    self.initialise(s, false);
    self = self.get_var();
    children.add(self);
  }

const char* Mistral::AddExpression::get_name() {
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
  int lb = children[0].get_min()+offset;
  int ub = children[0].get_max()+offset;

  Variable aux(lb, ub, DYN_VAR);
  self = aux;

  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
}

const char* Mistral::OffsetExpression::get_name() {
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


  Mistral::SubExpression::SubExpression(Variable X, Variable Y) 
  : Expression(X,Y) {
}
  Mistral::SubExpression::~SubExpression() {}
  
void Mistral::SubExpression::extract_constraint(Solver *s) {
      std::cerr << "Error: Sub predicate can't be used as a constraint" << std::endl;
  exit(0);
  }

  void Mistral::SubExpression::extract_variable(Solver *s) {
    //int lb = children[0].get_var().get_min()-children[1].get_var().get_max();
    //int ub = children[0].get_var().get_max()-children[1].get_var().get_min();
    
    int lb = children[0].get_min()-children[1].get_max();
    int ub = children[0].get_max()-children[1].get_min();
    Variable aux(lb, ub, DYN_VAR// , (RANGE_VAR | BOOL_VAR)
		 );
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

const char* Mistral::SubExpression::get_name() {
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
  std::cerr << "Error: Not predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::NotExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  self = aux;

  self.initialise(s, false);
  self = self.get_var();
  children.add(self);
}

const char* Mistral::NotExpression::get_name() {
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
  std::cout << "not implemented" << std::endl;
  exit(1);
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

const char* Mistral::AndExpression::get_name() {
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

const char* Mistral::OrExpression::get_name() {
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

// const char* Mistral::NeqExpression::get_name() {
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
  if(children.size==3) s->add(new PredicateEqual(children, spin));
  else s->add(new PredicateConstantEqual(children, value, spin));
}

const char* Mistral::EqualExpression::get_name() {
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

const char* Mistral::PrecedenceExpression::get_name() {
  return "prec";
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

const char* Mistral::DisjunctiveExpression::get_name() {
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

const char* Mistral::AllDiffExpression::get_name() {
  return "alldiff";
}

Mistral::Variable Mistral::AllDiff(Vector< Variable >& args, const int ct) {
  Variable exp(new AllDiffExpression(args,ct));
  return exp;
}



Mistral::BoolSumExpression::BoolSumExpression(Vector< Variable >& args, const int t) 
  : Expression(args) { total = t; }

Mistral::BoolSumExpression::~BoolSumExpression() {}

void Mistral::BoolSumExpression::extract_constraint(Solver *s) { 
  s->add(new ConstraintBoolSumEqual(children,total)); 
}

void Mistral::BoolSumExpression::extract_variable(Solver *s) {
      std::cerr << "Error: BoolSum constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

void Mistral::BoolSumExpression::extract_predicate(Solver *s) { 
      std::cerr << "Error: BoolSum constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::BoolSumExpression::get_name() {
  return "bool_sum";
}

Mistral::Variable Mistral::BoolSum(Vector< Variable >& args, const int t) {
  Variable exp(new BoolSumExpression(args,t));
  return exp;
}







// class MistralSum : public Predicate {

// private : 
  
//   inline int check_add() {
//     is_add == (is_add || 
// 	       (coefficients.size == 5 &&
// 		coefficients[4] == 0 && 
// 		coefficients[3] == 0));
//     if(is_add) {
//       if(coefficients[0]*coefficients[1]*coefficients[2] == -1 &&
// 	 coefficients[0]+coefficients[1]+coefficients[2] == 1) {
// 	int neg_idx = 0;
// 	while(coefficients[neg_idx] != -1) ++neg_idx;
//       }
//     } 
//     return -1;
//   }

// public :

//   // coefficients[0..children.size-1] stand for the coeffs of the linear sum
//   // coefficients[children.size] stands for the lower bound of the sum
//   // coefficients[children.size+1] stands for the upper bound of the sum
//   Vector< int > coefficients;
//   bool is_add;

//   MistralSum(Vector< Variable > args, Vector< int > coef) {
//     for(int i=0; i<args.size; ++i) {
//       children.add(args[i]);
//       coefficients.add(coef[i]);
//     }
// //     if(coef.size>args.size) {
// //       coefficients.add(coef[args.size]);
// //       coefficients.add(coef[args.size+1]);
// //     } else {
// //       coefficients.add(0);
// //       coefficients.add(0);
// //     }
//     is_add = false;
//   }
//   MistralSum(Variable X, Variable Y, Variable Z) {
//     children.add(X);
//     children.add(Y);
//     children.add(Z);
// //     coefficients.add(1);
// //     coefficients.add(1);
// //     coefficients.add(-1);
// //     coefficients.add(0);
// //     coefficients.add(0);
//     is_add = true;
//   }
//   ~MistralSum() {}
  
//   virtual Constraint* extract_constraint() { return NULL; }
//   virtual Constraint* extract_predicate() { 
//     // check if it is a simple 'Add'
//     int neg_idx = check_add();
//     if(neg_idx>=0) {
//       Vector< Variable > scope;
//       scope.add(children[(neg_idx+1)%3]);
//       scope.add(children[(neg_idx+2)%3]);
//       scope.add(children[neg_idx]);
//       return PredicateAdd(scope);
//     } else {
//       std::cout << "Not implemented" << std::endl;
//       exit(1);
//     } 


//   };
//   virtual void extract_variable(Variable& x) {
//     x.domain_type = 
//   }

// };


//  /// Remove all values that do not appear in the current domain of the Variable "x"
// Mistral::Event Mistral::VariableRange::setDomain(Variable x) {
//     //return setDomain(x.next(min-1), x.prev(max+1));
//     return NO_EVENT;
//   }


  
//   /// Whether the domain has a nonempty intersection with the Variable x
//   bool Mistral::VariableRange::intersect(Variable x) const { return x.intersect(min,max); }
//   /// Whether the domain is included in the Variable x 
//   bool Mistral::VariableRange::included(Variable x) const { return x.includes(min,max); }
//   /// Whether the domain is included in the Variable x 
//   bool Mistral::VariableRange::includes(Variable x) const { return x.included(min,max); }
 



