
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

// Mistral::Variable::Variable() {
//   domain_type = DYN_VAR;
//   implementation = NULL;
// }

// Mistral::Variable::Variable(const int value) {
//   domain_type = CONST_VAR;
//   implementation = (VariableImplementation*)value;
// }

// Mistral::Variable::Variable(VariableImplementation* impl, const int type) {
//   domain_type = type;
//   implementation = impl;
// }

// Mistral::Variable::Variable(Expression* impl) {
//   domain_type = EXPRESSION;
//   implementation = (VariableImplementation*)impl;

  

// //   std::cout << "create a new expression: "
// // 	    << impl << " " << domain_type << std::endl;
// //   //this << std::endl;

// }

// Mistral::Variable::Variable(const int lo, const int up, const int type) {

//   //if(type == CONST_VAR || (type == DYN_VAR && lo == up)) {
//   if((type & CONST_VAR) && lo == up) {

//     domain_type = CONST_VAR;
//     implementation = (VariableImplementation*)lo;
//     //} else if(type == BOOL_VAR || (type == DYN_VAR && lo==0 && up==1)) {
//   }//  else if((type & BOOL_VAR) && lo==0 && up==1) {

// //     //domain_type = (int)(s->getNextBooleanSlot());
// //     domain_type = BOOL_VAR;
// //     implementation = new VariableImplementation();
// //     //} else if(type == RANGE_VAR) {
// //   }
//   else if(type & RANGE_VAR) {

//     domain_type = RANGE_VAR;
//     implementation = new VariableRange(lo, up);
//   //       //} else if(type == LIST_VAR) {
//   //     } else if(type & LIST_VAR) {
//   //       domain_type = LIST_VAR;
//   //       implementation = new VariableList(lo, up);
//   } 
//   else {
//     domain_type = BITSET_VAR;
    
// #ifdef _BIT64
//     int nwords = ((up-lo) / 64);
//     if(nwords < 1) implementation = new VariableWord<unsigned long long int, 1>(lo, up);
//     else if(nwords < 2) implementation = new VariableWord<unsigned long long int, 2>(lo, up);
//     else if(nwords < 3) implementation = new VariableWord<unsigned long long int, 3>(lo, up);
// #else
//     int nwords = ((up-lo) / 32);
//     if(nwords < 1) implementation = new VariableWord<unsigned int, 1>(lo, up);
//     else if(nwords < 2) implementation = new VariableWord<unsigned int, 2>(lo, up);
//     else if(nwords < 3) implementation = new VariableWord<unsigned int, 3>(lo, up);
//     else if(nwords < 4) implementation = new VariableWord<unsigned int, 4>(lo, up);
//     else if(nwords < 5) implementation = new VariableWord<unsigned int, 5>(lo, up);
// #endif
//     else implementation = new VariableBitmap(lo, up);
//   }
// }


// Mistral::Variable Mistral::Variable::get_var() {
//   if(domain_type == EXPRESSION)
//     return ((Expression*)implementation)->self;
//   return *this;
// }

// std::ostream& Mistral::Variable::display(std::ostream& os) const {
//   if(domain_type == EXPRESSION) {
//     Expression *exp = (Expression*)implementation;
//     if(exp->is_initialised())
//       os << exp->self << ":";
//     os << exp->get_name()
//        << "(" << exp->children[0];
//     for(unsigned int i=1; i<exp->children.size-exp->is_initialised(); ++i) {
//       os << ", " << exp->children[i];
//     }
//     os << ")";
//   } else if(domain_type == CONST_VAR) {
//     os << implementation;
//   } else {
//     int id = implementation->id;
//     if       (domain_type ==  BITSET_VAR) {
//       os << "x" << id // << (VariableBitmap*)displayDomain(os)
// 	;
//     } // else if(domain_type ==    LIST_VAR) {
//       //       os << "y" << id;
//       //     } 
//     else if(domain_type ==   RANGE_VAR) {
//       os << "r" << id;
//     }
//       //     } else if(domain_type == VIRTUAL_VAR) {
//       //       return ((VariableVirtual *)implementation)->display(os);
//       //     }
//     else  {
//       os << "b" << id;
//     }
//   }
//   return os;
  
//   //if(*((int*)domain_type)  ==   CONST_VAR) os << (int)implementation;
//   //else os = implementation->display(os);
//   //return os;
// }

// void Mistral::Variable::initialise(Solver *s, const bool top) {
//   if(domain_type == EXPRESSION) {
    
// //     std::cout << "Add a new expression: " << this 
// //    	      << (top ? " at top-level" : " nested")
// //    	      << std::endl;
    
//     Expression *exp = (Expression*)implementation;
//     if(!exp->is_initialised()) {
//       for(unsigned int i=0; i<exp->children.size; ++i) {
// 	//std::cout << this << " add child" << std::endl; 
// 	exp->children[i].initialise(s, false);
//       }
//       if(top) {
// 	exp->extract_constraint(s);
// 	//std::cout << this << " extract constraint: " << s->constraints.back() << std::endl; 
//       } else {
	
// 	//std::cout << this << " extract variable " << std::endl;
// 	exp->extract_variable(s);
// 	exp->id = exp->self.id();
// 	exp->solver = s;
// 	exp->extract_predicate(s);//);
// 	//std::cout << "now: " << this << std::endl;
//       }
//     }

//     //else std::cout << this << " is known " << std::endl;
    
//   } else {
//     if(domain_type != CONST_VAR && implementation->solver != s) {

// //       std::cout << "declare a new variable: " << this 
// //   		<< " in " << get_domain() << std::endl;

//       implementation->solver = s;
//       implementation->id = s->declare(*this);

// //       if(domain_type == BOOL_VAR) {domain_type = 
// // 	  //(int)(new int[1]);
// // 	  (int)(s->getNextBooleanSlot());
// // 	//std::cout << "NEW DOMAIN: " << *((int*)domain_type) << std::endl;
// //       }

//       s->sequence.declare(*this);

// //       std::cout << "declare a new variable: " << this 
// //   		<< " in " << get_domain() << std::endl;
//     } 

//     //else std::cout << this << " is known " << std::endl;
//   }

//   //std::cout << "end initialise (" << id() << ")" << std::endl;
// }



// Mistral::Event Mistral::Variable::setValue( const int val ) 
//   {
// //     std::cout << "SET VALUE" << std::endl;
// //     std::cout << (int*)domain_type << std::endl;
// //     std::cout << *((int*)domain_type) << std::endl;

//     int nstat = (*((int*)domain_type) & val);
//     if( !nstat ) return FAIL_EVENT;

//     *((int*)domain_type) = nstat;

//     implementation->trigger_value_event_and_save(this);
//     return VALUE_EVENT;
//   }

// int Mistral::Variable::get_solution_value() const {
//   return(domain_type == CONST_VAR ? (int)implementation :
// 	 implementation->get_solution_value());
// }

//   int Mistral::Variable::get_value() const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->get_value();
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->get_value();
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->get_value();
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->get_value();
//     else if(domain_type ==   CONST_VAR) return (int)implementation;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.get_value();
//     else  return (*((int*)domain_type)-1);
//   }

//   unsigned int Mistral::Variable::get_size() const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->get_size();
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->get_size();
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->get_size();
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->get_size();
//     else if(domain_type ==   CONST_VAR) return 1;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.get_size();
//     else  return ((*((int*)domain_type)+1)/2);
//   }

// // Mistral::BitsetDomain Mistral::Variable::get_domain() const {
// //     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->domain;
// //     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->domain;
// //     else if(domain_type ==   RANGE_VAR) {

// //       Mistral::VariableRange* r = (Mistral::VariableRange*)implementation;
// //       BitsetDomain d(r->get_min(), r->get_max());

// //       std::cout << "HERE " << (r->get_min()) << " " << r->get_max() << std::endl;
// //       std::cout << d << std::endl;
      
// //       return d;
// //     }
// //     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->get_domain();
// //     else if(domain_type ==   CONST_VAR)  {
// //       BitsetDomain d((int)implementation, (int)implementation);
// //       return d;
// //     }
// //     else {
// //       int dom = *((int*)domain_type);
// //       BitsetDomain d((!(dom&1)),(dom/2));
// //       return d;
// //     }
// //   }

// std::string Mistral::Variable::get_domain() const {
//   std::ostringstream buf;
//   if     (domain_type ==  BITSET_VAR) buf << ((Mistral::VariableBitmap  *)implementation)->domain;
//   else if(domain_type ==    LIST_VAR) buf << ((Mistral::VariableList    *)implementation)->domain;
//   else if(domain_type ==   RANGE_VAR) {
//     Mistral::VariableRange* r = (Mistral::VariableRange*)implementation;
//     buf << "[" << r->get_min() << ".." <<  r->get_max() << "]";
//   }
//   //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->get_domain();
//   else if(domain_type ==   CONST_VAR)  {
//     buf << "{" << (int)implementation << "}";
//   }
//   else if(domain_type ==   BOOL_VAR)  {
//     buf << "{0,1}";
//   } 
//   else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.get_domain();
//   else {
//     int dom = *((int*)domain_type);
//     if(dom == 3) buf << "{0,1}";
//     else if(dom == 2) buf << "{1}";
//     else buf << "{0}";
//   }
//   return buf.str();
// }

//   int Mistral::Variable::get_min() const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->get_min();
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->get_min();
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->get_min();
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->get_min();
//     else if(domain_type ==   CONST_VAR) return (int)implementation;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.get_min();
//     else  return (!(*((int*)domain_type) & 1));
//   }

//   int Mistral::Variable::get_max() const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->get_max();
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->get_max();
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->get_max();
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->get_max();
//     else if(domain_type ==   CONST_VAR) return (int)implementation;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.get_max();
//     else  return (*((int*)domain_type) >> 1);
//   }

//   int Mistral::Variable::get_minCapacity() const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->get_minCapacity();
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->get_minCapacity();
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->get_minCapacity();
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->get_minCapacity();
//     else if(domain_type ==   CONST_VAR) return (int)implementation;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.get_minCapacity();
//     else  return 0;
//   }

//   int Mistral::Variable::get_maxCapacity() const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->get_maxCapacity();
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->get_maxCapacity();
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->get_maxCapacity();
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->get_maxCapacity();
//     else if(domain_type ==   CONST_VAR) return (int)implementation;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.get_maxCapacity();
//     else  return 1;
//   }

//   int Mistral::Variable::get_minPosAbs() const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->get_minPosAbs();
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->get_minPosAbs();
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->get_minPosAbs();
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->get_minPosAbs();
//     else if(domain_type ==   CONST_VAR) return (int)implementation;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.get_minPosAbs();
//     else  return (!(*((int*)domain_type) & 1));
//   }

//   int Mistral::Variable::get_minNegAbs() const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->get_minNegAbs();
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->get_minNegAbs();
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->get_minNegAbs();
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->get_minNegAbs();
//     else if(domain_type ==   CONST_VAR) return (int)implementation;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.get_minNegAbs();
//     else  return (!(*((int*)domain_type) & 1));
//   }

//   int Mistral::Variable::next(const int v) const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->next(v);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->next(v);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->next(v);
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->next(v);
//     else if(domain_type ==   CONST_VAR) return (int)implementation;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.next(v);
//     else  return (*((int*)domain_type) >> 1);
//   }

//   int Mistral::Variable::prev(const int v) const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->prev(v);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->prev(v);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->prev(v);
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->prev(v);
//     else if(domain_type ==   CONST_VAR) return (int)implementation;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.prev(v);
//     else  return (!(*((int*)domain_type) & 1));
//   }

//   bool Mistral::Variable::is_range() const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->is_range();
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->is_range();
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->is_range();
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->is_range();
//     else if(domain_type ==   CONST_VAR) return true;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.is_range();
//     else return true;
//   }

//   bool Mistral::Variable::is_ground() const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->is_ground();
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->is_ground();
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->is_ground();
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->is_ground();
//     else if(domain_type ==   CONST_VAR) return true;
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.is_ground();
//     else  return (*((int*)domain_type) != 3);
//   }

//   bool Mistral::Variable::equal(const int v) {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->equal(v);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->equal(v);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->equal(v);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->equal(v);
//     else if(domain_type ==   CONST_VAR) return ((int)implementation == v);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.equal(v);
//     else  return (*((int*)domain_type)-1 == v);
//   }

//   bool Mistral::Variable::contain(const int v) {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->contain(v);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->contain(v);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->contain(v);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->contain(v);
//     else if(domain_type ==   CONST_VAR) return ((int)implementation == v);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.contain(v);
//     else  return (!(v >> 1) && (*((int*)domain_type) & (v+1)));
//   }

//   bool Mistral::Variable::intersect(const int lo, const int up) const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->intersect(lo, up);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->intersect(lo, up);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->intersect(lo, up);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->intersect(lo, up);
//     else if(domain_type ==   CONST_VAR) return ((int)implementation >= lo && (int)implementation <= up);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.intersect(lo, up);
//     else  return ((lo<=0 | 2*(up>0)) & *((int*)domain_type));
//   }

//   bool Mistral::Variable::included(const int lo, const int up) const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->included(lo, up);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->included(lo, up);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->included(lo, up);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->included(lo, up);
//     else if(domain_type ==   CONST_VAR) return ((int)implementation >= lo && (int)implementation <= up);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.included(lo, up);
//     else  {
//       int state = *((int*)domain_type);
//       return ( up >= (state >> 1) && (lo <= !(state & 1)) );
//     }
//   }

//   bool Mistral::Variable::includes(const int lo, const int up) const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->includes(lo, up);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->includes(lo, up);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->includes(lo, up);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->includes(lo, up);
//     else if(domain_type ==   CONST_VAR) return ((int)implementation == lo && (int)implementation == up);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.includes(lo, up);
//     else  {
//       int state = *((int*)domain_type);
//       return ( up <= (state >> 1) && (lo >= !(state & 1)) );
//     }
//   }

//   bool Mistral::Variable::intersect(const BitSet& s) const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->intersect(s);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->intersect(s);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->intersect(s);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->intersect(s);
//     else if(domain_type ==   CONST_VAR) return (s.contain((int)implementation));
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.intersect(s);
//     else  return s.intersect(*((int*)domain_type));
//   }

//   bool Mistral::Variable::included(const BitSet& s) const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->included(s);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->included(s);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->included(s);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->included(s);
//     else if(domain_type ==   CONST_VAR) return (s.contain((int)implementation));
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.included(s);
//     else  return s.includes(*((int*)domain_type));
//   }

//   bool Mistral::Variable::includes(const BitSet& s) const {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->includes(s);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->includes(s);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->includes(s);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->includes(s);
//     else if(domain_type ==   CONST_VAR) return (s.size() == 1 && s.contain((int)implementation));
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.includes(s);
//     else  return s.included(*((int*)domain_type));
//   }

//   bool Mistral::Variable::intersect(const Mistral::Variable& x) const {
//     //     if     (domain_type ==  BITSET_VAR) return x->intersect(((Mistral::VariableBitmap  *)implementation)->domain);
//     //     else if(domain_type ==    LIST_VAR) return x->intersect(((Mistral::VariableList    *)implementation)->domain); //((Mistral::VariableList    *)implementation)->intersect(x);
//     //     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->intersect(x);
//     //     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->intersect(x);
//     //     else if(domain_type ==   CONST_VAR) return x.contain((int)implementation);
//     //     else  {
//     //       int state = *((int*)domain_type);
//     //       if(state == 3) return x.intersect(0,1); 
//     //       else if(state == 2) return x.contain(1); 
//     //       else return x.contain(0);
//     //     }
//     return false;
//   }

//   bool Mistral::Variable::included(const Mistral::Variable& x) const {
//     //     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->included(x);
//     //     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->included(x);
//     //     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->included(x);
//     //     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->included(x);
//     //     else if(domain_type ==   CONST_VAR) return x.contain((int)implementation);
//     //     else  {
//     //       int state = *((int*)domain_type);
//     //       if(state == 3) return x.includes(0,1); 
//     //       else if(state == 2) return x.contain(1); 
//     //       else return x.contain(0);
//     //     }
//     return false;
//   }

//   bool Mistral::Variable::includes(const Mistral::Variable& x) const {
//     //     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->includes(x);
//     //     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->includes(x);
//     //     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->includes(x);
//     //     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->includes(x);
//     //     else if(domain_type ==   CONST_VAR) return (x.size() == 1 && x.contain((int)implementation));
//     //     else  {
//     //       int state = *((int*)domain_type);
//     //       if(state == 3) return x.included(0,1); 
//     //       else if(state == 2) return x.equal(1); 
//     //       else return x.equal(0);
//     //     }
//     return false;
//   }

//   void Mistral::Variable::intersectTo( BitSet& s ) const {
//     if     (domain_type ==  BITSET_VAR) ((Mistral::VariableBitmap  *)implementation)->intersectTo(s);
//     else if(domain_type ==    LIST_VAR) ((Mistral::VariableList    *)implementation)->intersectTo(s);
//     else if(domain_type ==   RANGE_VAR) ((Mistral::VariableRange   *)implementation)->intersectTo(s);
//     else if(domain_type == VIRTUAL_VAR) ((Mistral::VariableVirtual *)implementation)->intersectTo(s);
//     else if(domain_type ==   CONST_VAR) {
//       if(s.contain((int)implementation)) {
// 	s.clear();
// 	s.add((int)implementation);
//       } else s.clear();
//     }
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.intersectTo(s);
//     else  return s.intersectWith(*((int*)domain_type));
//   }

//   void Mistral::Variable::unionTo( BitSet& s ) const {
//     if     (domain_type ==  BITSET_VAR) ((Mistral::VariableBitmap  *)implementation)->unionTo(s);
//     else if(domain_type ==    LIST_VAR) ((Mistral::VariableList    *)implementation)->unionTo(s);
//     else if(domain_type ==   RANGE_VAR) ((Mistral::VariableRange   *)implementation)->unionTo(s);
//     else if(domain_type == VIRTUAL_VAR) ((Mistral::VariableVirtual *)implementation)->unionTo(s);
//     else if(domain_type ==   CONST_VAR) s.add((int)implementation);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.unionTo(s);
//     else  return s.unionWith(*((int*)domain_type));
//   }

//   Mistral::Event Mistral::Variable::remove(const int v) {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->remove(v);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->remove(v);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->remove(v);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->remove(v);
//     else if(domain_type ==   CONST_VAR) return ((int)implementation == v ? FAIL_EVENT : NO_EVENT);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.remove(v);
//     else  return ((v>1 || v<0) ? NO_EVENT : setValue(2-v));
//   }

//   Mistral::Event Mistral::Variable::setDomain(const int v) {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->setDomain(v);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->setDomain(v);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->setDomain(v);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->setDomain(v);
//     else if(domain_type ==   CONST_VAR) return ((int)implementation != v ? FAIL_EVENT : NO_EVENT);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.setDomain(v);
//     else  return ((v>1 || v<0) ? NO_EVENT : setValue(1+v));
//   }

//   Mistral::Event Mistral::Variable::setMin(const int lo) {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->setMin(lo);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->setMin(lo);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->setMin(lo);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->setMin(lo);
//     else if(domain_type ==   CONST_VAR) return ((int)implementation < lo ? FAIL_EVENT : NO_EVENT);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.setMin(lo);
//     else  return (lo<1 ? NO_EVENT : (lo>1 ? FAIL_EVENT : setValue(2)));
//   }

//   Mistral::Event Mistral::Variable::setMax(const int up) {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->setMax(up);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->setMax(up);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->setMax(up);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->setMax(up);
//     else if(domain_type ==   CONST_VAR) return ((int)implementation > up ? FAIL_EVENT : NO_EVENT);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.setMax(up);
//     else  return (up>0 ? NO_EVENT : (up<0 ? FAIL_EVENT : setValue(1)));
//   }

//   Mistral::Event Mistral::Variable::setDomain(const BitSet& s) {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->setDomain(s);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->setDomain(s);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->setDomain(s);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->setDomain(s);
//     else if(domain_type ==   CONST_VAR) return (s.contain((int)implementation) ? NO_EVENT : FAIL_EVENT);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.setDomain(s);
//     else  return ((s.pos_words<1 || s.neg_words>0) ? FAIL_EVENT : ((s.table[0]&3)==3 ? NO_EVENT : setValue(s.table[0])));
//   }

//   Mistral::Event Mistral::Variable::setDomain(Mistral::Variable& x) {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->setDomain(x);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->setDomain(x);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->setDomain(x);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->setDomain(x);
//     else if(domain_type ==   CONST_VAR) return (x.contain((int)implementation) ? NO_EVENT : FAIL_EVENT);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.setDomain(x);
//     else  {
//       if(!x.contain(0)) {
// 	if(!x.contain(1)) return FAIL_EVENT;
// 	else return setValue(2);
//       } else {
// 	if(!x.contain(1)) return setValue(1);
// 	else return NO_EVENT;
//       }
//     }
//   }

//   Mistral::Event Mistral::Variable::removeSet(const BitSet& s) {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->removeSet(s);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->removeSet(s);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->removeSet(s);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->removeSet(s);
//     else if(domain_type ==   CONST_VAR) return (s.contain((int)implementation) ? FAIL_EVENT : NO_EVENT);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.removeSet(s);
//     else  return ((s.pos_words<1 || s.neg_words>0 || (s.table[0]^3)==3) ? NO_EVENT : setValue(s.table[0]^3));
//   }

//   Mistral::Event Mistral::Variable::removeRange(const int lo, const int up) {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->removeRange(lo, up);
//     else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->removeRange(lo, up);
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->removeRange(lo, up);
//     else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->removeRange(lo, up);
//     else if(domain_type ==   CONST_VAR) return (((int)implementation < lo || (int)implementation > up) ? NO_EVENT : FAIL_EVENT);
//     else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.removeRange(lo, up);
//     else  return (lo==1 ? setValue(1) : (up==0 ? setValue(2) : ((lo>1 || up<0) ? NO_EVENT : FAIL_EVENT)));
//   }

//   Mistral::Event Mistral::Variable::restore() {
//     if     (domain_type ==  BITSET_VAR) return ((Mistral::VariableBitmap  *)implementation)->restore();
//     //else if(domain_type ==    LIST_VAR) return ((Mistral::VariableList    *)implementation)->restore();
//     else if(domain_type ==   RANGE_VAR) return ((Mistral::VariableRange   *)implementation)->restore();
//     //else if(domain_type == VIRTUAL_VAR) return ((Mistral::VariableVirtual *)implementation)->restore();
//     else if(domain_type ==   CONST_VAR) return NO_EVENT;
//     //else if(domain_type ==   EXPRESSION) return ((Expression*)implementation)->self.restore();
//     else {
//       *((int*)domain_type) = 3;
//       return NO_EVENT;
//     } 
//   }


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

  

//   std::cout << "create a new expression: "
// 	    << impl << " " << domain_type << std::endl;
//   //this << std::endl;

}

Mistral::Variable::Variable(const int lo, const int up, const int type) {

  //if(type == CONST_VAR || (type == DYN_VAR && lo == up)) {
  if((type & CONST_VAR) && lo == up) {

    domain_type = CONST_VAR;
    constant_value = lo;
    //} else if(type == BOOL_VAR || (type == DYN_VAR && lo==0 && up==1)) {
  }//  else if((type & BOOL_VAR) && lo==0 && up==1) {

//     //domain_type = (int)(s->getNextBooleanSlot());
//     domain_type = BOOL_VAR;
//     implementation = new VariableImplementation();
//     //} else if(type == RANGE_VAR) {
//   }
  else if(type & RANGE_VAR) {

    domain_type = RANGE_VAR;
    range_domain = new VariableRange(lo, up);
  //       //} else if(type == LIST_VAR) {
  //     } else if(type & LIST_VAR) {
  //       domain_type = LIST_VAR;
  //       implementation = new VariableList(lo, up);
  } 
  else {
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
  if(domain_type == EXPRESSION)
    return expression->self;
  return *this;
}

std::ostream& Mistral::Variable::display(std::ostream& os) const {
  if(domain_type == EXPRESSION) {
    //Expression *exp = (Expression*)implementation;
    if(expression->is_initialised())
      os << expression->self << ":";
    os << expression->get_name()
       << "(" << expression->children[0];
    for(unsigned int i=1; i<expression->children.size-expression->is_initialised(); ++i) {
      os << ", " << expression->children[i];
    }
    os << ")";
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
    
//     std::cout << "Add a new expression: " << this 
//    	      << (top ? " at top-level" : " nested")
//    	      << std::endl;
    
    //Expression *exp = (Expression*)implementation;
    if(!expression->is_initialised()) {
      for(unsigned int i=0; i<expression->children.size; ++i) {
	//std::cout << this << " add child" << std::endl; 
	expression->children[i].initialise(s, false);
      }
      if(top) {
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
    if(domain_type != CONST_VAR && variable->solver != s) {

//       std::cout << "declare a new variable: " << this 
//   		<< " in " << get_domain() << std::endl;

      variable->solver = s;
      variable->id = s->declare(*this);

//       if(domain_type == BOOL_VAR) {domain_type = 
// 	  //(int)(new int[1]);
// 	  (int)(s->getNextBooleanSlot());
// 	//std::cout << "NEW DOMAIN: " << *((int*)domain_type) << std::endl;
//       }

      s->sequence.declare(*this);

//       std::cout << "declare a new variable: " << this 
//   		<< " in " << get_domain() << std::endl;
    } 

    //else std::cout << this << " is known " << std::endl;
  }

  //std::cout << "end initialise (" << id() << ")" << std::endl;
}



Mistral::Event Mistral::Variable::setValue( const int val ) 
  {
//     std::cout << "SET VALUE" << std::endl;
//     std::cout << (int*)domain_type << std::endl;
//     std::cout << *bool_domain << std::endl;

    int nstat = (*bool_domain & val);
    if( !nstat ) return FAIL_EVENT;

    *bool_domain = nstat;

    variable->trigger_value_event_and_save(this);
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

// Mistral::BitsetDomain Mistral::Variable::get_domain() const {
//     if     (domain_type ==  BITSET_VAR) return bitset_domain->domain;
//     else if(domain_type ==    LIST_VAR) return list_domain->domain;
//     else if(domain_type ==   RANGE_VAR) {

//       Mistral::VariableRange* r = (Mistral::VariableRange*)implementation;
//       BitsetDomain d(r->get_min(), r->get_max());

//       std::cout << "HERE " << (r->get_min()) << " " << r->get_max() << std::endl;
//       std::cout << d << std::endl;
      
//       return d;
//     }
//     //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_domain();
//     else if(domain_type ==   CONST_VAR)  {
//       BitsetDomain d(constant_value, constant_value);
//       return d;
//     }
//     else {
//       int dom = *bool_domain;
//       BitsetDomain d((!(dom&1)),(dom/2));
//       return d;
//     }
//   }

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

  bool Mistral::Variable::is_ground() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->is_ground();
    else if(domain_type ==    LIST_VAR) return list_domain->is_ground();
    else if(domain_type ==   RANGE_VAR) return range_domain->is_ground();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->is_ground();
    else if(domain_type ==   CONST_VAR) return true;
    else if(domain_type ==   EXPRESSION) return expression->self.is_ground();
    else  return (*bool_domain != 3);
  }

  bool Mistral::Variable::equal(const int v) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->equal(v);
    else if(domain_type ==    LIST_VAR) return list_domain->equal(v);
    else if(domain_type ==   RANGE_VAR) return range_domain->equal(v);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->equal(v);
    else if(domain_type ==   CONST_VAR) return (constant_value == v);
    else if(domain_type ==   EXPRESSION) return expression->self.equal(v);
    else  return (*bool_domain-1 == v);
  }

  bool Mistral::Variable::contain(const int v) {
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
    else  return ((lo<=0 | 2*(up>0)) & *bool_domain);
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
    //     if     (domain_type ==  BITSET_VAR) return x->intersect(bitset_domain->domain);
    //     else if(domain_type ==    LIST_VAR) return x->intersect(list_domain->domain); //list_domain->intersect(x);
    //     else if(domain_type ==   RANGE_VAR) return range_domain->intersect(x);
    //     else if(domain_type == VIRTUAL_VAR) return virtual_domain->intersect(x);
    //     else if(domain_type ==   CONST_VAR) return x.contain(constant_value);
    //     else  {
    //       int state = *bool_domain;
    //       if(state == 3) return x.intersect(0,1); 
    //       else if(state == 2) return x.contain(1); 
    //       else return x.contain(0);
    //     }
    return false;
  }

  bool Mistral::Variable::included(const Mistral::Variable& x) const {
    //     if     (domain_type ==  BITSET_VAR) return bitset_domain->included(x);
    //     else if(domain_type ==    LIST_VAR) return list_domain->included(x);
    //     else if(domain_type ==   RANGE_VAR) return range_domain->included(x);
    //     else if(domain_type == VIRTUAL_VAR) return virtual_domain->included(x);
    //     else if(domain_type ==   CONST_VAR) return x.contain(constant_value);
    //     else  {
    //       int state = *bool_domain;
    //       if(state == 3) return x.includes(0,1); 
    //       else if(state == 2) return x.contain(1); 
    //       else return x.contain(0);
    //     }
    return false;
  }

  bool Mistral::Variable::includes(const Mistral::Variable& x) const {
    //     if     (domain_type ==  BITSET_VAR) return bitset_domain->includes(x);
    //     else if(domain_type ==    LIST_VAR) return list_domain->includes(x);
    //     else if(domain_type ==   RANGE_VAR) return range_domain->includes(x);
    //     else if(domain_type == VIRTUAL_VAR) return virtual_domain->includes(x);
    //     else if(domain_type ==   CONST_VAR) return (x.size() == 1 && x.contain(constant_value));
    //     else  {
    //       int state = *bool_domain;
    //       if(state == 3) return x.included(0,1); 
    //       else if(state == 2) return x.equal(1); 
    //       else return x.equal(0);
    //     }
    return false;
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
    else  return ((v>1 || v<0) ? NO_EVENT : setValue(2-v));
  }

  Mistral::Event Mistral::Variable::setDomain(const int v) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->setDomain(v);
    else if(domain_type ==    LIST_VAR) return list_domain->setDomain(v);
    else if(domain_type ==   RANGE_VAR) return range_domain->setDomain(v);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->setDomain(v);
    else if(domain_type ==   CONST_VAR) return (constant_value != v ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.setDomain(v);
    else  return ((v>1 || v<0) ? NO_EVENT : setValue(1+v));
  }

  Mistral::Event Mistral::Variable::setMin(const int lo) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->setMin(lo);
    else if(domain_type ==    LIST_VAR) return list_domain->setMin(lo);
    else if(domain_type ==   RANGE_VAR) return range_domain->setMin(lo);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->setMin(lo);
    else if(domain_type ==   CONST_VAR) return (constant_value < lo ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.setMin(lo);
    else  return (lo<1 ? NO_EVENT : (lo>1 ? FAIL_EVENT : setValue(2)));
  }

  Mistral::Event Mistral::Variable::setMax(const int up) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->setMax(up);
    else if(domain_type ==    LIST_VAR) return list_domain->setMax(up);
    else if(domain_type ==   RANGE_VAR) return range_domain->setMax(up);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->setMax(up);
    else if(domain_type ==   CONST_VAR) return (constant_value > up ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.setMax(up);
    else  return (up>0 ? NO_EVENT : (up<0 ? FAIL_EVENT : setValue(1)));
  }

  Mistral::Event Mistral::Variable::setDomain(const BitSet& s) {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->setDomain(s);
    else if(domain_type ==    LIST_VAR) return list_domain->setDomain(s);
    else if(domain_type ==   RANGE_VAR) return range_domain->setDomain(s);
    else if(domain_type == VIRTUAL_VAR) return virtual_domain->setDomain(s);
    else if(domain_type ==   CONST_VAR) return (s.contain(constant_value) ? NO_EVENT : FAIL_EVENT);
    else if(domain_type ==   EXPRESSION) return expression->self.setDomain(s);
    else  return ((s.pos_words<1 || s.neg_words>0) ? FAIL_EVENT : ((s.table[0]&3)==3 ? NO_EVENT : setValue(s.table[0])));
  }

  Mistral::Event Mistral::Variable::setDomain(Mistral::Variable& x) {
//     if     (domain_type ==  BITSET_VAR) return bitset_domain->setDomain(x);
//     else if(domain_type ==    LIST_VAR) return list_domain->setDomain(x);
//     else if(domain_type ==   RANGE_VAR) return range_domain->setDomain(x);
//     else if(domain_type == VIRTUAL_VAR) return virtual_domain->setDomain(x);
//     else if(domain_type ==   CONST_VAR) return (x.contain(constant_value) ? NO_EVENT : FAIL_EVENT);
//     else if(domain_type ==   EXPRESSION) return expression->self.setDomain(x);
//     else  {
//       if(!x.contain(0)) {
// 	if(!x.contain(1)) return FAIL_EVENT;
// 	else return setValue(2);
//       } else {
// 	if(!x.contain(1)) return setValue(1);
// 	else return NO_EVENT;
//       }
//     }
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
    //else if(domain_type ==   EXPRESSION) return expression->self.restore();
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


void Mistral::VariableImplementation::trigger_value_event_and_save(Mistral::Variable *x) {
  solver->trigger_event(id, VALUE_EVENT);
  solver->save(*x);
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


  Mistral::Expression::Expression(Vector< Variable >& args) {
    id=-1; 
    for(unsigned int i=0; i<args.size; ++i)
      children.add(args[i]);
  }
Mistral::Expression::~Expression() {}

  Mistral::BinaryExpression::BinaryExpression(Variable X, Variable Y) 
    : Expression() {
    children.add(X);
    children.add(Y);
  }
  Mistral::BinaryExpression::~BinaryExpression() {}


  Mistral::AddExpression::AddExpression(Variable X, Variable Y) 
  : BinaryExpression(X,Y) {}
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

    Variable aux(lb, ub// , (RANGE_VAR | BOOL_VAR)
		 );
    self = aux;

    children.add(self);
    self.initialise(s, false);
  }

const char* Mistral::AddExpression::get_name() {
  return "add";
}

  void Mistral::AddExpression::extract_predicate(Solver *s) {
    s->add(new PredicateAdd(children));
  }

Mistral::Variable Mistral::Variable::operator+(Variable x) {
  Variable exp(new AddExpression(*this,x));
  return exp;
}


  Mistral::SubExpression::SubExpression(Variable X, Variable Y) 
  : BinaryExpression(X,Y) {}
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
    Variable aux(lb, ub// , (RANGE_VAR | BOOL_VAR)
		 );
    self = aux;

    children.add(self);
    self.initialise(s, false);
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




  Mistral::NeqExpression::NeqExpression(Variable X, Variable Y) 
  : BinaryExpression(X,Y) {}
  Mistral::NeqExpression::~NeqExpression() {}
  
void Mistral::NeqExpression::extract_constraint(Solver *s) {
//     Constraint *neq = new ConstraintNotEqual(children);
//     neq->initialise();
//     return neq;
  s->add(new ConstraintNotEqual(children));
  }

  void Mistral::NeqExpression::extract_variable(Solver *s) {
    Variable aux(0, 1, BOOL_VAR);
    self = aux;

    children.add(self);
    self.initialise(s, false);
  }

  void Mistral::NeqExpression::extract_predicate(Solver *s) {
//     Constraint *neq = new PredicateEqual(children, false);
//     neq->initialise();
//     return neq;
    s->add(new PredicateEqual(children, false));
  }

const char* Mistral::NeqExpression::get_name() {
  return "neq";
}

Mistral::Variable Mistral::Variable::operator!=(Variable x) {
  Variable exp(new NeqExpression(*this,x));
  return exp;
}




Mistral::PrecedenceExpression::PrecedenceExpression(Variable X, Variable Y, const int of) 
  : BinaryExpression(X,Y) { offset = of; }
Mistral::PrecedenceExpression::~PrecedenceExpression() {offset = 0;}
  
void Mistral::PrecedenceExpression::extract_constraint(Solver *s) {
  s->add(new ConstraintLess(children, offset));
  }

  void Mistral::PrecedenceExpression::extract_variable(Solver *s) {
    Variable aux(0, 1, BOOL_VAR);
    self = aux;

    children.add(self);
    self.initialise(s, false);
  }

  void Mistral::PrecedenceExpression::extract_predicate(Solver *s) {
//     Constraint *prec = new PredicateLess(children, offset);
//     prec->initialise();
//     return prec;
//    return NULL;
      std::cerr << "Error: Precedence constraint can't yet be used as a predicate" << std::endl;
  exit(0);
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
      std::cerr << "Error: Precedence constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

void Mistral::AllDiffExpression::extract_predicate(Solver *s) { 
      std::cerr << "Error: Precedence constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::AllDiffExpression::get_name() {
  return "alldiff";
}

Mistral::Variable Mistral::AllDiff(Vector< Variable >& args, const int ct) {
  Variable exp(new AllDiffExpression(args,ct));
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
 



