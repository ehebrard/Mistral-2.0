
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


#include <mistral_constraint.hpp>
#include <mistral_variable.hpp>
#include <mistral_solver.hpp>


#define _DEBUG_BUILD true



#define PROFILING_HEAD \
double __t_tag__ = 0;		 \
  if(domain_type != CONST_VAR) { \
  __t_tag__ = get_run_time(); \
  } \

#define PROFILING_FOOT(method) \
if(domain_type != CONST_VAR) {	\
  __t_tag__ = get_run_time() - __t_tag__; \
  int __idx__ = VARTYPE[(domain_type < DYN_VAR ? domain_type : BOOL_VAR)]; \
  variable->solver->statistics.prof_time[method][__idx__] += __t_tag__; \
  ++variable->solver->statistics.prof_num[method][__idx__]; \
  } \



/**
Constructors for variable
*/

// 
Mistral::Variable::Variable() {
  domain_type = DYN_VAR;
  variable = NULL;
}

// Constant variable
Mistral::Variable::Variable(const int value) {
  domain_type = NULL;
  domain_type = CONST_VAR;
  variable = NULL;
  constant_value = value;
}

Mistral::Variable::Variable(VariableImplementation* impl, const int type) {
  domain_type = NULL;
  domain_type = type;
  variable = impl;
}

Mistral::Variable::Variable(Expression* exp) {
  domain_type = NULL;
  domain_type = EXPRESSION;
  expression = exp;
}

int BOOL_DOM = 3;

Mistral::Variable::Variable(Vector< int >& values, const int type) {
  if(values.back() - values.front() + 1 == (int)(values.size)) {
    initialise_domain(values.front(), values.back(), type);
  } else {
    initialise_domain(values, type);
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
    if(nwords == 1) bitset_domain = new VariableWord<unsigned long long int, 1>(lo, up);
    else if(nwords == 2) bitset_domain = new VariableWord<unsigned long long int, 2>(lo, up);
    else if(nwords == 3) bitset_domain = new VariableWord<unsigned long long int, 3>(lo, up);
#else
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
  
  // SELF CHANGE
  // if(domain_type == EXPRESSION) {
  //   //Variable x = expression->self;
  //   //while(x.is_expression)
  //   return expression->self.get_var();
  // } else 
 
  if(domain_type == CONST_VAR || !variable->solver) {
    return *this;
  } else if(variable->id == -2) return expression->_self.get_var();
  return ((Solver*)(variable->solver))->variables[variable->id];
}

const Mistral::Variable Mistral::Variable::get_var() const {
  // SELF CHANGE
  // if(domain_type == EXPRESSION) {
  //   return expression->self.get_var();
  // } else 
  
  if(domain_type == CONST_VAR || !variable->solver) {
    return *this;
  }
  return ((Solver*)(variable->solver))->variables[variable->id];
}

Mistral::Solver* Mistral::Variable::get_solver() {
  //std::cout << (int*)(variable) << std::endl;
  return (Solver*)(variable->solver);
}



std::ostream& Mistral::Variable::display(std::ostream& os) const {
  if(domain_type == EXPRESSION) {
    expression->display(os);
    //}
  } else if(domain_type == CONST_VAR) {
    os << constant_value;
  } else {
    int id = variable->id;
    
    if (domain_type ==  BITSET_VAR) {
      os << "x" ;
    } // else if(domain_type ==    LIST_VAR) {
      //       os << "y" << id;
      //     } 
    else if(domain_type ==   RANGE_VAR) {
      os << "r" ;
    }
    else  {
      os << "b" ;
    }

    if(variable->is_initialised()) {
      os << id ;
    } else {
      os << "_";
    }

  }
  return os;
  
  //if(*((int*)domain_type)  ==   CONST_VAR) os << constant_value;
  //else os = implementation->display(os);
  //return os;
}

void Mistral::Variable::initialise(Solver *s, const int level) {


  // std::cout << "beg initialise " ;
  // display(std::cout);
  // std::cout << std::endl;

  if(domain_type == EXPRESSION) {
    if(!expression->is_initialised()) {

      // std::cout << "initialise children of " << std::endl;
      // //Variable X(expression);

      // std::cout << *this << std::endl;

      // for(unsigned int i=0; i<expression->children.size; ++i) {
      //  std::cout << "  " << expression->children[i] << std::endl;
      // }
     

      for(unsigned int i=0; i<expression->children.size; ++i) {
	expression->children[i].initialise(s, level+1);
      }

      // std::cout << "**initialise " << std::endl;
      // //Variable X(expression);

      // std::cout << *this << std::endl;

      if(level == 0 && !expression->children.empty()) {

	//std::cout << "-> constraint!" << std::endl;

	expression->extract_constraint(s);
      } else {

	//std::cout << "-> predicate! (extract var)" << std::endl;

	expression->extract_variable(s);



	// SELF CHANGE
	Variable X = expression->_self;
	//expression->id = (X.domain_type == CONST_VAR ? s->variables.size-1 : X.id());
	expression->id = (X.domain_type == CONST_VAR ? -2 : X.id());

	//expression->id = X.id();
	expression->solver = s;	


	//std::cout << " extracted " << X << " -> predicate " << std::endl;

	expression->extract_predicate(s);//);

	//std::cout << "done" << std::endl;

      }
      s->expression_store.add(expression);
    }
    //s->add_var(*this);
    //s->declared_variables[id()] = *this;
  } else {
    if(domain_type != CONST_VAR && variable->solver != s) {
      s->declare(this->get_var());
      s->sequence.declare(*this);
    } 
    //s->add_var(*this);
  }



  // std::cout << "end initialise " ;
  // display(std::cout);
  // std::cout << std::endl;

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
  Event evt = VALUE_C;
  int dom = *bool_domain;
  //int nstat = val&3;

  //std::cout << "1 set " << id() << " dom (" << dom << ") to " << val << std::endl;
  
  if( val == dom ) return NO_EVENT;
  else if( dom<3 // || val<1 || val>2
	   ) return FAIL_EVENT;
  
  *bool_domain = val;
  //std::cout << "1 VALUE EVENT" << std::endl;
  
  //variable->trigger_value_event_and_save(this);
  //variable->trigger_value_event_and_save();

  evt |= (val==1 ? UB_EVENT : LB_EVENT);
  //solver->trigger_event(id, evt);
  //solver->save(*x);
  //solver->save(this, BOOL_VAR);
  //solver->save(id);
  variable->trigger_event_and_save(evt);

  return evt;
}


Mistral::Event Mistral::Variable::setState( const int vals ) 
{
  // vals should be 1, 2 or 3
  Event evt = VALUE_C;
  int dom = *bool_domain;
  int ndom = vals&dom;
  //int nstat = val&3;

  //std::cout << "2 set " << id() << " dom (" << dom << ") to " << vals << std::endl;
  
  if( ndom == dom ) return NO_EVENT;
  else if( !ndom ) return FAIL_EVENT;
  
  *bool_domain = ndom;
  //std::cout << "2 VALUE EVENT" << std::endl;
  

  evt |= (ndom==1 ? UB_EVENT : LB_EVENT);
  //solver->trigger_event(id, evt);
  //solver->save(*x);
  //solver->save(this, BOOL_VAR);
  //solver->save(id);
  //variable->trigger_value_event_and_save(this);
  variable->trigger_event_and_save(evt);

  return evt;
}

int Mistral::Variable::get_solution_int_value() const {
  int value = 0;
  if(is_initialised()) {
    // return(domain_type == CONST_VAR ? constant_value :
    // 	   variable->get_solution_int_value());
    value = variable->get_solution_int_value();
  } else {
    value = get_min();
  }
  return value;
}

std::string Mistral::Variable::get_solution_str_value() const {
  std::ostringstream ret_str;

  if(is_initialised()) {
    // if(domain_type == CONST_VAR)
    //   ret_str << constant_value ;
    // else
    ret_str << variable->get_solution_str_value();
  } else {
    ret_str << get_min();
  }

  return ret_str.str();
}

int Mistral::Variable::get_solution_min() const {
  int value = 0;
  if(is_initialised()) {
    value = variable->get_solution_min();
  } else {
    value = get_min();
  }
  return value;
  // return(domain_type == CONST_VAR ? constant_value :
  // 	 variable->get_solution_min());
}

int Mistral::Variable::get_solution_max() const {
  int value = 0;
  if(is_initialised()) {
    value = variable->get_solution_max();
  } else {
    value = get_max();
  }
  return value;
  // return(domain_type == CONST_VAR ? constant_value :
  // 	 variable->get_solution_max());
}

  int Mistral::Variable::get_value() const {
    if     (domain_type ==  BITSET_VAR) return bitset_domain->get_value();
    else if(domain_type ==    LIST_VAR) return list_domain->get_value();
    else if(domain_type ==   RANGE_VAR) return range_domain->get_value();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_value();
    else if(domain_type ==   CONST_VAR) return constant_value;
    else if(domain_type ==   EXPRESSION) return expression->get_self().get_value();
    else  return (*bool_domain-1);
  }

unsigned int Mistral::Variable::get_size() const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif
      
      unsigned int r_size = 0;

  if     (domain_type ==  BITSET_VAR) r_size = bitset_domain->get_size();
  else if(domain_type ==    LIST_VAR) r_size = list_domain->get_size();
  else if(domain_type ==   RANGE_VAR) r_size = range_domain->get_size();
  //else if(domain_type == VIRTUAL_VAR) r_size = virtual_domain->get_size();
  else if(domain_type ==   CONST_VAR) r_size = 1;
  else if(domain_type ==   EXPRESSION) r_size = expression->get_self().get_size();
  else  r_size = ((*bool_domain+1)/2);

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_get_size_)
#endif

    return r_size;
}

/// Returns the degree (number of constraints)
unsigned int Mistral::Variable::get_degree() const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      unsigned int r_degree = NO_EVENT;

  if(domain_type ==   CONST_VAR) r_degree = 0;
  else if(domain_type ==   EXPRESSION) r_degree = expression->get_self().get_degree();
  else r_degree = ((Solver*)(variable->solver))->constraint_graph[variable->id].size();


#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_get_degree_)
#endif

    return r_degree;
}


std::string Mistral::Variable::get_domain() const {
  std::ostringstream buf;
  if     (domain_type ==  BITSET_VAR) buf << bitset_domain->domain;
  else if(domain_type ==    LIST_VAR) buf << list_domain->domain;
  else if(domain_type ==   RANGE_VAR) {
    //Mistral::VariableRange* r = (Mistral::VariableRange*)implementation;
    if(range_domain->get_min() == range_domain->get_max())
      buf << range_domain->get_min() ;
    else if(range_domain->get_min() == range_domain->get_max()-1)
      buf << "[" << range_domain->get_min() << "," <<  range_domain->get_max() << "]";
    else 
      buf << "[" << range_domain->get_min() << ".." <<  range_domain->get_max() << "]";
  }
  //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_domain();
  else if(domain_type ==   CONST_VAR)  {
    buf << constant_value ;
  }
  else if(domain_type ==   BOOL_VAR)  {
    buf << "[0,1]";
  } 
  else if(domain_type ==   EXPRESSION) {
    buf << "*" << expression->get_self().get_domain();
  } else {
    if(*bool_domain == 3) buf << "[0,1]";
    else if(*bool_domain == 2) buf << "1";
    else buf << "0";
  }
  return buf.str();
}


std::string Mistral::Variable::get_history() const {
  std::ostringstream buf;
  if     (domain_type ==  BITSET_VAR) buf << bitset_domain->get_history();
  else if(domain_type ==   RANGE_VAR) buf << range_domain->get_history();
  else if(domain_type >      DYN_VAR && *bool_domain < 3) buf << "[0,1]";
  return buf.str();
}

  int Mistral::Variable::get_min() const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

    // if     (domain_type ==  BITSET_VAR) {
    //   //std::cout << "bitset " << bitset_domain->get_min() << std::endl;
    //   return bitset_domain->get_min();
    // } else if(domain_type ==    LIST_VAR)  {
    //   //std::cout << "list" << std::endl;
    //   return list_domain->get_min();
    // } else if(domain_type ==   RANGE_VAR)  {
    //   //std::cout << "range " << range_domain->get_min() << std::endl;
    //   return range_domain->get_min();
    // //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_min();
    // } else if(domain_type ==   CONST_VAR)  {
    //   //std::cout << "constant" << std::endl;
    //   return constant_value;
    // } else if(domain_type ==   EXPRESSION)  {
    //   //std::cout << "expression" << expression->get_self().get_min() << std::endl;
    //   return expression->get_self().get_min();
    // } else  return (!(*bool_domain & 1));

      
    //   std::cout << "get min of "  ;
    // std::cout.flush();
    // display(std::cout);
    // std::cout << std::endl;
    
    int of_the_living_dead = 0;
    if     (domain_type ==  BITSET_VAR) {
      //std::cout << "bitset " << bitset_domain->get_min() << std::endl;
      of_the_living_dead =  bitset_domain->get_min();
    } else if(domain_type ==    LIST_VAR)  {
      //std::cout << "list" << std::endl;
      of_the_living_dead =  list_domain->get_min();
    } else if(domain_type ==   RANGE_VAR)  {
      //std::cout << "range " << range_domain->get_min() << std::endl;
      of_the_living_dead =  range_domain->get_min();
      //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_min();
    } else if(domain_type ==   CONST_VAR)  {
      //std::cout << "constant" << std::endl;
      of_the_living_dead =  constant_value;
    } else if(domain_type ==   EXPRESSION)  {
      //std::cout << "expression" << expression->get_self().get_min() << std::endl;
      of_the_living_dead =  expression->get_self().get_min();
    } else  of_the_living_dead = !(*bool_domain & 1);


#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_get_min_)
#endif

    return of_the_living_dead;
  }

  int Mistral::Variable::get_max() const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

    int of_the_mummy = 0;

    if     (domain_type ==  BITSET_VAR) of_the_mummy = bitset_domain->get_max();
    else if(domain_type ==    LIST_VAR) of_the_mummy = list_domain->get_max();
    else if(domain_type ==   RANGE_VAR) of_the_mummy = range_domain->get_max();
    //else if(domain_type == VIRTUAL_VAR) of_the_mummy = virtual_domain->get_max();
    else if(domain_type ==   CONST_VAR) of_the_mummy = constant_value;
    else if(domain_type ==   EXPRESSION) of_the_mummy = expression->get_self().get_max();
    else  of_the_mummy = (*bool_domain >> 1);

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_get_max_)
#endif

    return of_the_mummy;
  }

  int Mistral::Variable::get_initial_min() const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      int r_min = 0;

    if     (domain_type ==  BITSET_VAR) r_min = bitset_domain->get_initial_min();
    else if(domain_type ==    LIST_VAR) r_min = list_domain->get_initial_min();
    else if(domain_type ==   RANGE_VAR) r_min = range_domain->get_initial_min();
    //else if(domain_type == VIRTUAL_VAR) r_min = virtual_domain->get_initial_min();
    else if(domain_type ==   CONST_VAR) r_min = constant_value;
    else if(domain_type ==   EXPRESSION) r_min = expression->get_self().get_initial_min();
    //else  r_min = 0;

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_get_initial_min_)
#endif

    return r_min;
  }

  int Mistral::Variable::get_initial_max() const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      int r_max = 1;

    if     (domain_type ==  BITSET_VAR) r_max = bitset_domain->get_initial_max();
    else if(domain_type ==    LIST_VAR) r_max = list_domain->get_initial_max();
    else if(domain_type ==   RANGE_VAR) r_max = range_domain->get_initial_max();
    //else if(domain_type == VIRTUAL_VAR) r_max = virtual_domain->get_initial_max();
    else if(domain_type ==   CONST_VAR) r_max = constant_value;
    else if(domain_type ==   EXPRESSION) r_max = expression->get_self().get_initial_max();
    //else  r_max = 1;

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_get_initial_max_)
#endif

    return r_max;
  }

int Mistral::Variable::get_min_pos() const {
  
#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    
    int r_min = 0;

  if     (domain_type ==  BITSET_VAR) r_min = bitset_domain->get_min_pos();
  else if(domain_type ==    LIST_VAR) r_min = list_domain->get_min_pos();
  else if(domain_type ==   RANGE_VAR) r_min = range_domain->get_min_pos();
  //else if(domain_type == VIRTUAL_VAR) r_min = virtual_domain->get_min_pos();
  else if(domain_type ==   CONST_VAR) r_min = constant_value;
  else if(domain_type ==   EXPRESSION) r_min = expression->get_self().get_min_pos();
  else  r_min = (*bool_domain >> 1); //(!(*bool_domain & 1));

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_get_min_pos_)
#endif

    return r_min;

}

  int Mistral::Variable::get_max_neg() const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      int r_max = 0;

    if     (domain_type ==  BITSET_VAR) r_max = bitset_domain->get_max_neg();
    else if(domain_type ==    LIST_VAR) r_max = list_domain->get_max_neg();
    else if(domain_type ==   RANGE_VAR) r_max = range_domain->get_max_neg();
    //else if(domain_type == VIRTUAL_VAR) r_max = virtual_domain->get_max_neg();
    else if(domain_type ==   CONST_VAR) r_max = constant_value;
    else if(domain_type ==   EXPRESSION) r_max = expression->get_self().get_max_neg();
    else  r_max = (!(*bool_domain & 1));

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_get_max_neg_)
#endif

    return r_max;
  }

  int Mistral::Variable::next(const int v) const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      int r_val = v;

    if     (domain_type ==  BITSET_VAR) r_val = bitset_domain->next(v);
    else if(domain_type ==    LIST_VAR) r_val = list_domain->next(v);
    else if(domain_type ==   RANGE_VAR) r_val = range_domain->next(v);
    //else if(domain_type == VIRTUAL_VAR) r_val = virtual_domain->next(v);
    else if(domain_type ==   CONST_VAR) r_val = constant_value;
    else if(domain_type ==   EXPRESSION) r_val = expression->get_self().next(v);
    else  r_val = (*bool_domain >> 1);

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_next_)
#endif

    return r_val;
  }

  int Mistral::Variable::prev(const int v) const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      int r_val = v;

    if     (domain_type ==  BITSET_VAR) r_val = bitset_domain->prev(v);
    else if(domain_type ==    LIST_VAR) r_val = list_domain->prev(v);
    else if(domain_type ==   RANGE_VAR) r_val = range_domain->prev(v);
    //else if(domain_type == VIRTUAL_VAR) r_val = virtual_domain->prev(v);
    else if(domain_type ==   CONST_VAR) r_val = constant_value;
    else if(domain_type ==   EXPRESSION) r_val = expression->get_self().prev(v);
    else  r_val = (!(*bool_domain & 1));

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_prev_)
#endif

    return r_val;
  }

  bool Mistral::Variable::is_range() const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      bool answer = true;

    if     (domain_type ==  BITSET_VAR) answer = bitset_domain->is_range();
    else if(domain_type ==    LIST_VAR) answer = list_domain->is_range();
    //else if(domain_type ==   RANGE_VAR) answer = range_domain->is_range();
    //else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->is_range();
    //else if(domain_type ==   CONST_VAR) answer = true;
    else if(domain_type ==   EXPRESSION) answer = expression->get_self().is_range();
    //else answer = true;

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_is_range_)
#endif

    return answer;
  }

// // bool Mistral::Variable::is_ground(const Expression *x) const {
// //   return x->self.is_ground();
// // }
//   bool Mistral::Variable::is_constant() const {

// #ifdef _PROFILING_PRIMITIVE
//     PROFILING_HEAD
// #endif

//       bool answer = true;

//     // std::cout << (int*)(variable)
//     // 	      << std::endl
//     // 	      << domain_type << " == " 
//     // 	      << BITSET_VAR << "? "
//     // 	      << LIST_VAR << "? "
//     // 	      << RANGE_VAR << "? "
//     // 	      << CONST_VAR << "? " << std::endl;


//     //std::cout << domain2str(domain_type) << std::endl;

//     if     (domain_type ==  BITSET_VAR) answer = bitset_domain->is_ground();
//     else if(domain_type ==    LIST_VAR) answer = list_domain->is_ground();
//     else if(domain_type ==   RANGE_VAR) answer = range_domain->is_ground();
//     //else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->is_ground();
//     else if(domain_type ==   CONST_VAR) answer = true;
//     else if(domain_type ==   EXPRESSION) answer = expression->get_self().is_ground();
//     else  answer = (*bool_domain != 3);

// #ifdef _PROFILING_PRIMITIVE
//     PROFILING_FOOT(_m_is_ground_)
// #endif

//     return answer;
//   }


  bool Mistral::Variable::is_ground() const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      bool answer = true;

    // std::cout << (int*)(variable)
    // 	      << std::endl
    // 	      << domain_type << " == " 
    // 	      << BITSET_VAR << "? "
    // 	      << LIST_VAR << "? "
    // 	      << RANGE_VAR << "? "
    // 	      << CONST_VAR << "? " << std::endl;


    //std::cout << domain2str(domain_type) << std::endl;

    if     (domain_type ==  BITSET_VAR) answer = bitset_domain->is_ground();
    else if(domain_type ==    LIST_VAR) answer = list_domain->is_ground();
    else if(domain_type ==   RANGE_VAR) answer = range_domain->is_ground();
    //else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->is_ground();
    else if(domain_type ==   CONST_VAR) answer = true;
    else if(domain_type ==   EXPRESSION) answer = expression->get_self().is_ground();
    else  answer = (*bool_domain != 3);

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_is_ground_)
#endif

    return answer;
  }

  bool Mistral::Variable::equal(const int v) const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      bool answer = true;

    if     (domain_type ==  BITSET_VAR) answer = bitset_domain->equal(v);
    else if(domain_type ==    LIST_VAR) answer = list_domain->equal(v);
    else if(domain_type ==   RANGE_VAR) answer = range_domain->equal(v);
    else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->equal(v);
    else if(domain_type ==   CONST_VAR) answer = (constant_value == v);
    else if(domain_type ==   EXPRESSION) answer = expression->get_self().equal(v);
    else  answer = (*bool_domain-1 == v);

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_equal_)
#endif

    return answer;
  }

  bool Mistral::Variable::contain(const int v) const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      bool answer = true;

    if     (domain_type ==  BITSET_VAR) answer = bitset_domain->contain(v);
    else if(domain_type ==    LIST_VAR) answer = list_domain->contain(v);
    else if(domain_type ==   RANGE_VAR) answer = range_domain->contain(v);
    else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->contain(v);
    else if(domain_type ==   CONST_VAR) answer = (constant_value == v);
    else if(domain_type ==   EXPRESSION) answer = expression->get_self().contain(v);
    else  answer = (!(v >> 1) && (*bool_domain & (v+1)));

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_contain_)
#endif

    return answer;
  }

  bool Mistral::Variable::intersect(const int lo, const int up) const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      bool answer = true;

    if     (domain_type ==  BITSET_VAR) answer = bitset_domain->intersect(lo, up);
    else if(domain_type ==    LIST_VAR) answer = list_domain->intersect(lo, up);
    else if(domain_type ==   RANGE_VAR) answer = range_domain->intersect(lo, up);
    else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->intersect(lo, up);
    else if(domain_type ==   CONST_VAR) answer = (constant_value >= lo && constant_value <= up);
    else if(domain_type ==   EXPRESSION) answer = expression->get_self().intersect(lo, up);
    else  answer = (((lo<=0) | (2*(up>0))) & *bool_domain);

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_intersect_range_)
#endif

    return answer;
  }

  bool Mistral::Variable::included(const int lo, const int up) const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      bool answer = true;

    if     (domain_type ==  BITSET_VAR) answer = bitset_domain->included(lo, up);
    else if(domain_type ==    LIST_VAR) answer = list_domain->included(lo, up);
    else if(domain_type ==   RANGE_VAR) answer = range_domain->included(lo, up);
    else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->included(lo, up);
    else if(domain_type ==   CONST_VAR) answer = (constant_value >= lo && constant_value <= up);
    else if(domain_type ==   EXPRESSION) answer = expression->get_self().included(lo, up);
    else  {
      int state = *bool_domain;
      answer = ( up >= (state >> 1) && (lo <= !(state & 1)) );
    }

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_included_range_)
#endif

    return answer;
  }

  bool Mistral::Variable::includes(const int lo, const int up) const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      bool answer = true;

    if     (domain_type ==  BITSET_VAR) answer = bitset_domain->includes(lo, up);
    else if(domain_type ==    LIST_VAR) answer = list_domain->includes(lo, up);
    else if(domain_type ==   RANGE_VAR) answer = range_domain->includes(lo, up);
    else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->includes(lo, up);
    else if(domain_type ==   CONST_VAR) answer = (constant_value == lo && constant_value == up);
    else if(domain_type ==   EXPRESSION) answer = expression->get_self().includes(lo, up);
    else  {
      int state = *bool_domain;
      answer = ( up <= (state >> 1) && (lo >= !(state & 1)) );
    }

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_includes_range_)
#endif

    return answer;
  }

  bool Mistral::Variable::intersect(const BitSet& s) const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

      bool answer = true;

    if     (domain_type ==  BITSET_VAR) answer = bitset_domain->intersect(s);
    else if(domain_type ==    LIST_VAR) answer = list_domain->intersect(s);
    else if(domain_type ==   RANGE_VAR) answer = range_domain->intersect(s);
    else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->intersect(s);
    else if(domain_type ==   CONST_VAR) answer = (s.contain(constant_value));
    else if(domain_type ==   EXPRESSION) answer = expression->get_self().intersect(s);
    else  answer = s.intersect(*bool_domain);

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_intersect_set_)
#endif

    return answer;
  }

  bool Mistral::Variable::included(const BitSet& s) const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif
      bool answer = true;

    if     (domain_type ==  BITSET_VAR) answer = bitset_domain->included(s);
    else if(domain_type ==    LIST_VAR) answer = list_domain->included(s);
    else if(domain_type ==   RANGE_VAR) answer = range_domain->included(s);
    else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->included(s);
    else if(domain_type ==   CONST_VAR) answer = (s.contain(constant_value));
    else if(domain_type ==   EXPRESSION) answer = expression->get_self().included(s);
    else  answer = s.includes(*bool_domain);

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_included_set_)
#endif

    return answer;
  }

  bool Mistral::Variable::includes(const BitSet& s) const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif
      bool answer = true;

    if     (domain_type ==  BITSET_VAR) answer = bitset_domain->includes(s);
    else if(domain_type ==    LIST_VAR) answer = list_domain->includes(s);
    else if(domain_type ==   RANGE_VAR) answer = range_domain->includes(s);
    else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->includes(s);
    else if(domain_type ==   CONST_VAR) answer = (s.size() == 1 && s.contain(constant_value));
    else if(domain_type ==   EXPRESSION) answer = expression->get_self().includes(s);
    else  answer = s.included(*bool_domain);

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_includes_set_)
#endif

    return answer;
  }

  bool Mistral::Variable::intersect(const Mistral::Variable& x) const {

      bool answer = true;

    if(is_ground()) answer = x.contain(get_min());
    else if(x.is_ground()) answer = contain(x.get_min());
    else if(is_range()) answer = x.intersect(get_min(), get_max());
    else if(x.is_range()) answer = intersect(x.get_min(), x.get_max());
    else if(domain_type ==  BITSET_VAR)
      answer = x.intersect(bitset_domain->domain.values);
    // else if(domain_type ==    LIST_VAR) {
    //   std::cout << "TODO! (intersect - list)" << std::endl;
    //   exit(1);
    // }
    // // else if(domain_type == VIRTUAL_VAR) {

    // // }
    else 
      answer = //x.intersect(expression->get_self()); //
	expression->get_self().intersect(x);

   return answer;
  }

  bool Mistral::Variable::included(const Mistral::Variable& x) const {

      bool answer = true;

    if(is_ground()) answer = x.contain(get_min());
    else if(x.is_ground()) answer = equal(x.get_min());
    else if(is_range()) answer = x.includes(get_min(), get_max());
    else if(x.is_range()) answer = included(x.get_min(), x.get_max()); 
    else if(domain_type ==  BITSET_VAR)
      answer = x.includes(bitset_domain->domain.values);
    //std::cout << "TODO! (included)" << std::endl;
    //answer = true;
    else 
      answer = //x.includes(expression->get_self());//.included(x);
	expression->get_self().included(x);

   return answer;
  }

  bool Mistral::Variable::includes(const Mistral::Variable& x) const {

      bool answer = true;

    if(is_ground()) answer = x.equal(get_min());
    else if(x.is_ground()) answer = contain(x.get_min());
    else if(is_range()) answer = x.included(get_min(), get_max());
    else if(x.is_range()) answer = includes(x.get_min(), x.get_max()); 
    else if(domain_type ==  BITSET_VAR)
      answer = x.included(bitset_domain->domain.values);
    // std::cout << "TODO!" << std::endl;
    // answer = true;
    else 
      answer = 
	expression->get_self().includes(x);

   return answer;
  }

  void Mistral::Variable::intersect_to( BitSet& s ) const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

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
    else if(domain_type ==   EXPRESSION) expression->get_self().intersect_to(s);
    else  s.intersect_with(*bool_domain);

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_intersect_to_)
#endif

  }

  void Mistral::Variable::union_to( BitSet& s ) const {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif

    if     (domain_type ==  BITSET_VAR) bitset_domain->union_to(s);
    else if(domain_type ==    LIST_VAR) list_domain->union_to(s);
    else if(domain_type ==   RANGE_VAR) range_domain->union_to(s);
    else if(domain_type == VIRTUAL_VAR) {
      
      // std::cout << domain_type << std::endl;

      // std::cout << "make the union of " << *this // << " in " << get_domain()
      // 		<< " into " << s << std::endl;

      virtual_domain->union_to(s);
    } else if(domain_type ==   CONST_VAR) s.add(constant_value);
    else if(domain_type ==   EXPRESSION) expression->get_self().union_to(s);
    else s.union_with(*bool_domain);

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_union_to_)
#endif

  }

  Mistral::Event Mistral::Variable::remove(const int v) {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif
      Event evt = NO_EVENT;

    if     (domain_type ==  BITSET_VAR) evt = bitset_domain->remove(v);
    else if(domain_type ==    LIST_VAR) evt = list_domain->remove(v);
    else if(domain_type ==   RANGE_VAR) evt = range_domain->remove(v);
    else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->remove(v);
    else if(domain_type ==   CONST_VAR) evt = (constant_value == v ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) evt = expression->get_self().remove(v);
    else {
      evt = (v<0||v>1 ? NO_EVENT : setValue(2-v));
    }

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_remove_)
#endif

    return evt;
  }

  Mistral::Event Mistral::Variable::set_domain(const int v) {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif
      Event evt = NO_EVENT;

    if     (domain_type ==  BITSET_VAR) evt = bitset_domain->set_domain(v);
    else if(domain_type ==    LIST_VAR) evt = list_domain->set_domain(v);
    else if(domain_type ==   RANGE_VAR) evt = range_domain->set_domain(v);
    else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->set_domain(v);
    else if(domain_type ==   CONST_VAR) evt = (constant_value != v ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) evt = expression->get_self().set_domain(v);
    else {

      evt = (v<0||v>1 ? FAIL_EVENT : setValue(1+v));
      //int dom = *bool_domain;
      //evt = ((dom==1+v) ? NO_EVENT : setValue(1+v));
    }

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_set_domain_value_)
#endif

    return evt;
  }

  Mistral::Event Mistral::Variable::set_min(const int lo) {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif
      Event evt = NO_EVENT;

    if     (domain_type ==  BITSET_VAR) evt = bitset_domain->set_min(lo);
    else if(domain_type ==    LIST_VAR) evt = list_domain->set_min(lo);
    else if(domain_type ==   RANGE_VAR) evt = range_domain->set_min(lo);
    else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->set_min(lo);
    else if(domain_type ==   CONST_VAR) evt = (constant_value < lo ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) evt = expression->get_self().set_min(lo);
    else {

      // 1 [1, inf[
      // 2 [2, inf[
      // 3 [2, inf[
      evt = (lo==1 ? setValue(2) : (lo>1 ? FAIL_EVENT : NO_EVENT));

//       int dom = *bool_domain;
//       evt = (lo<=(!(dom&1)) ? NO_EVENT : (lo>(dom>>1) ? FAIL_EVENT : setValue(2)));
//       //evt = (lo<1 ? NO_EVENT : (lo>1 ? FAIL_EVENT : setValue(2)));
    }

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_set_min_)
#endif

    return evt;
  }

  Mistral::Event Mistral::Variable::set_max(const int up) {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif
      Event evt = NO_EVENT;

    if     (domain_type ==  BITSET_VAR) evt = bitset_domain->set_max(up);
    else if(domain_type ==    LIST_VAR) evt = list_domain->set_max(up);
    else if(domain_type ==   RANGE_VAR) evt = range_domain->set_max(up);
    else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->set_max(up);
    else if(domain_type ==   CONST_VAR) evt = (constant_value > up ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) evt = expression->get_self().set_max(up);
    else {

      evt =(up==0 ? setValue(1) : (up<0 ? FAIL_EVENT : NO_EVENT));

//       int dom = *bool_domain;
//       evt = (up>=(dom>>1) ? NO_EVENT : (up<(dom&1) ? FAIL_EVENT : setValue(1)));
      //evt = (up>0 ? NO_EVENT : (up<0 ? FAIL_EVENT : setValue(1)));
    }

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_set_max_)
#endif

    return evt;
  }

  Mistral::Event Mistral::Variable::set_domain(const BitSet& s) {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif
      Event evt = NO_EVENT;

    if     (domain_type ==  BITSET_VAR) evt = bitset_domain->set_domain(s);
    else if(domain_type ==    LIST_VAR) evt = list_domain->set_domain(s);
    else if(domain_type ==   RANGE_VAR) evt = range_domain->set_domain(s);
    else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->set_domain(s);
    else if(domain_type ==   CONST_VAR) evt = (s.contain(constant_value) ? NO_EVENT : FAIL_EVENT);
    else if(domain_type ==   EXPRESSION) evt = expression->get_self().set_domain(s);
    else {
      //evt = ((s.pos_words<1 || s.neg_words>0) ? FAIL_EVENT : ((s.table[0]&3)==3 ? NO_EVENT : setValue(s.table[0])));
      evt = ((s.pos_words<1 || s.neg_words>0) ? FAIL_EVENT : setState(s.table[0]&*bool_domain));
    }

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_set_domain_set_)
#endif

    return evt;
  }

  Mistral::Event Mistral::Variable::set_domain(Mistral::Variable& x) {

    Event evt = NO_EVENT;
    
    //if(x.is_ground()) evt = set_domain(x.get_min());
    if(x.is_ground()) {
      evt =  set_domain(x.get_min());
    }
    else if(x.is_range()) {
      Event evt = (set_min(x.get_min()) | set_max(x.get_max()));
      evt = evt;
    }
    else if(x.domain_type ==  BITSET_VAR) evt = set_domain(x.bitset_domain->domain.values);
    else if(x.domain_type ==  EXPRESSION) {
      Variable y = x.expression->get_self();
      evt = set_domain(y);
    } else {
      std::cout << "TODO! (set_domain(var))" << std::endl;
      exit(1);
    }

    return evt;

    //evt = NO_EVENT;
  }

  Mistral::Event Mistral::Variable::removeSet(const BitSet& s) {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif
      Event evt = NO_EVENT;

    if     (domain_type ==  BITSET_VAR) evt = bitset_domain->removeSet(s);
    else if(domain_type ==    LIST_VAR) evt = list_domain->removeSet(s);
    else if(domain_type ==   RANGE_VAR) evt = range_domain->removeSet(s);
    else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->removeSet(s);
    else if(domain_type ==   CONST_VAR) evt = (s.contain(constant_value) ? FAIL_EVENT : NO_EVENT);
    else if(domain_type ==   EXPRESSION) evt = expression->get_self().removeSet(s);
    else {

      evt = ((s.pos_words<1 || s.neg_words>0 || (s.table[0]^3)==3) ? NO_EVENT : setValue(s.table[0]^3));
    }

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_remove_set_)
#endif

    return evt;
  }

  Mistral::Event Mistral::Variable::remove_interval(const int lo, const int up) {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif
      Event evt = NO_EVENT;

    if     (domain_type ==  BITSET_VAR) evt = bitset_domain->remove_interval(lo, up);
    else if(domain_type ==    LIST_VAR) evt = list_domain->remove_interval(lo, up);
    else if(domain_type ==   RANGE_VAR) evt = range_domain->remove_interval(lo, up);
    else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->remove_interval(lo, up);
    else if(domain_type ==   CONST_VAR) evt = ((constant_value < lo || constant_value > up) ? NO_EVENT : FAIL_EVENT);
    else if(domain_type ==   EXPRESSION) evt = expression->get_self().remove_interval(lo, up);
    else {

      evt = (lo==1 ? setValue(1) : (up==0 ? setValue(2) : ((lo>1 || up<0) ? NO_EVENT : FAIL_EVENT)));
    }

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_remove_interval_)
#endif

    return evt;
  }




  Mistral::Event Mistral::Variable::restore() {

#ifdef _PROFILING_PRIMITIVE
    PROFILING_HEAD
#endif
      Event evt = NO_EVENT;

    // if(id() == 17) {
    //   for(int i=0; i<((VariableImplementation*)(variable))->solver->level; ++i)
    // 	std::cout << " ";
    //   std::cout << "restore " ;
    //   display(std::cout);
    //   std::cout << " in " << get_domain() << " " 
    // 		<< ((VariableBitmap*)(variable))->trail_ << std::endl;
    // }

    if     (domain_type ==  BITSET_VAR) evt = bitset_domain->restore();
    //else if(domain_type ==    LIST_VAR) evt = list_domain->restore();
    else if(domain_type ==   RANGE_VAR) evt = range_domain->restore();
    //else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->restore();
    else if(domain_type ==   CONST_VAR) evt = NO_EVENT;
    else if(domain_type ==   EXPRESSION) {
      //std::cout << "RESTORE EXPRESSION" << std::endl;
      //exit(1);
      evt = expression->get_self().restore();
    }
    else {
      *bool_domain = 3;
      //evt = NO_EVENT;
    } 

#ifdef _PROFILING_PRIMITIVE
    PROFILING_FOOT(_m_restore_)
#endif

    return evt;
  }


Mistral::Event Mistral::VariableRange::remove(const int v) {
  Event removal = DOMAIN_EVENT;


  //std::cout << "remove " << v << " from " << this << " in [" << min << ".." << max << "]" << std::endl;

  
  // first check if we can abort early
  if(min>v || max<v) {
    return NO_EVENT;
  }
  if(min!=v && max!=v) {
    ((Solver*)solver)->make_non_convex(id);
    removal = ((Solver*)solver)->variables[id].remove(v);
    return removal;
    //return NO_EVENT;
  }
  if(min==max) return FAIL_EVENT;
  
  save();
  
  if(min==v) {
    ++min;
    removal |= LB_EVENT;
  } else {
    --max;
    removal |= UB_EVENT;
  }
  
  if(min == max) removal |= VALUE_C; 
  solver->trigger_event(id, removal);

  //std::cout << removal << std::endl;



  return removal; 
}


   /// Remove all values that do not appear in the set "s"
Mistral::Event Mistral::VariableRange::set_domain(const BitSet& s) {
      Event setdomain = NO_EVENT;

      // std::cout << "include " << s.includes(min, max) << std::endl;
      // std::cout << "intersect " << s.intersect(min, max) << std::endl;
      // std::cout << "[" << s.next(min-1) << ".." << s.prev(max+1) << "]" << std::endl;
      

      if(s.includes(min, max)) return NO_EVENT;
      if(!s.intersect(min, max)) return FAIL_EVENT;
      int lb = s.next(min-1);
      int ub = s.prev(max+1);
      if(s.includes(lb, ub)) {
	if(lb>min) {
	  setdomain |= set_min(lb);
	} 
	if(ub<max) {
	  setdomain |= set_max(ub);
	}
      } else {

	// std::cout << "the intersection is not convex" << std::endl;

	((Solver*)solver)->make_non_convex(id);


	// std::cout << solver->variables[id] << " in " 
	// 	  << solver->variables[id].get_domain() << std::endl;

	return ((Solver*)solver)->variables[id].set_domain(s);
      }

      //return set_domain(s.next(min-1), s.prev(max+1));
      return setdomain;
    }


// bool Mistral::VariableImplementation::is_new(Solver *s) {
//   return (solver != s);
// }

// void Mistral::VariableImplementation::initialise(Solver *s) {
//   solver = s;
//   id = s->declare(*this);
// }

int Mistral::VariableImplementation::get_solution_int_value() const { 
  return ((Solver*)solver)->last_solution_lb[id] ;
}  
std::string Mistral::VariableImplementation::get_solution_str_value() const { 
  std::ostringstream ret_str;
  ret_str <<  ((Solver*)solver)->last_solution_lb[id] ;
  return ret_str.str();
}  
int Mistral::VariableImplementation::get_solution_min() const { 
  return ((Solver*)solver)->last_solution_lb[id] ; 
} 
int Mistral::VariableImplementation::get_solution_max() const { 
  return ((Solver*)solver)->last_solution_ub[id] ; 
}  


//void Mistral::VariableImplementation::trigger_value_event_and_save(Mistral::Variable *x) {
void Mistral::VariableImplementation::trigger_value_event_and_save() {
  solver->trigger_event(id, VALUE_EVENT);
  //solver->save(*x);
  //solver->save(this, BOOL_VAR);
  solver->save(id);
}

void Mistral::VariableImplementation::trigger_event_and_save(const Event evt) {
  solver->trigger_event(id, evt);
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

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Goal& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Goal* x) {
  return x->display(os);
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
  _self = x;
}

Mistral::Expression::Expression(Vector< int >& values) 
  : VariableImplementation() {
  id=-1; 
  Variable x(values, DYN_VAR);
  _self = x;
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
  _self = x;
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
  _self = x;
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
#ifdef _DEBUG_MEMORY
  std::cout << "c delete expression" << std::endl;
#endif
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
  // SELF CHANGE

  _self.initialise(s, 1);
  _self = _self.get_var();

  //solver = s;
  //Variable self(lb, ub, );
}

// void Mistral::Expression::extract_variable(Solver *s) {
//   self.initialise(s, 1);
//   self = self.get_var();
// }

Mistral::Variable Mistral::Expression::get_self() { return (id >= 0 ? ((Solver*)solver)->variables[id] : _self); }

std::ostream& Mistral::Expression::display(std::ostream& os) const {
  if(is_initialised())
    os << "e" << id << ":";
  os << get_name() << "(" ;
  if(children.empty()) os << (id>=0 ? ((Solver*)solver)->variables[id] : _self);
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
Mistral::AddExpression::~AddExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete add expression" << std::endl;
#endif
}

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
    _self = aux;

    _self.initialise(s, 1);
    _self = _self.get_var();
    children.add(_self);
  }

const char* Mistral::AddExpression::get_name() const {
  return "add";
}

void Mistral::AddExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateAdd(children)));
}


Mistral::OffsetExpression::OffsetExpression(Variable X, const int ofs) 
  : Expression(X) { 
  offset=ofs; 
}
Mistral::OffsetExpression::~OffsetExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete offset expression" << std::endl;
#endif
}
  
void Mistral::OffsetExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Offset predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::OffsetExpression::extract_variable(Solver *s) {
//void Mistral::OffsetExpression::reify(Solver *s, Variable X) {
  int lb = children[0].get_min()+offset;
  int ub = children[0].get_max()+offset;

    Variable aux(lb, ub, DYN_VAR);
    _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::OffsetExpression::get_name() const {
  return "offset";
}

void Mistral::OffsetExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateOffset(children, offset)));
}

Mistral::Variable Mistral::Variable::operator+(Variable x) {
  Variable exp(new AddExpression(*this,x));
  return exp;
}

Mistral::Variable Mistral::Variable::operator+(int k) {
  Variable exp(new OffsetExpression(*this,k));
  return exp;
}

Mistral::Variable Mistral::Variable::operator-(int k) {
  Variable exp(new OffsetExpression(*this,-k));
  return exp;
}

Mistral::MulExpression::MulExpression(Variable X, Variable Y) 
  : Expression(X,Y) {
}
Mistral::MulExpression::~MulExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete mul expression" << std::endl;
#endif
}

void Mistral::MulExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Mul predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::MulExpression::extract_variable(Solver *s) {

  int z[4], lb = INFTY, ub = -INFTY, i = 4, j, k = 0, b[4];


  // std::cout << "extract variable for " << children[0].get_domain() << " x "
  //  	    << children[1].get_domain() << std::endl;



  while(i-->2) { // std::cout << i-2 << std::endl;
    b[i] = children[i-2].get_min(); }
  do { // std::cout << i << std::endl;
    b[i] = children[i].get_max(); } while(i--);

  //0.1
  //0.3
  

  //std::cout << b[0] << " " << b[1] << " " << b[2] << " " << b[3] << std::endl; 

  for(i=0; i<4; i+=2) 
    for(j=1; j<4; j+=2) {
      z[k] = b[i]*b[j];
      if(((b[i]>0 && b[j]>0) || (b[i]<0 && b[j]<0)) && z[k]<0) // int overflow
	z[k] = INFTY;
      else if(((b[i]>0 && b[j]<0) || (b[i]<0 && b[j]>0)) && z[k]>0) // int overflow
	z[k] = -INFTY;


      //std::cout << i << "." << j << " " << b[i] << " * " << b[j] << std::endl;

      ++k;
    }
  
  i = 4;
  while( i-- ) {
    if( z[i] > ub ) ub = z[i];
    if( z[i] < lb ) lb = z[i];
  }

  // std::cout << children[0].get_domain() << " x "
  //  	    << children[1].get_domain() << " = ["
  //  	    << lb << ".." << ub << "]" << std::endl;


  // exit(1);

  Variable aux(lb, ub, DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::MulExpression::get_name() const {
  return "mul";
}

void Mistral::MulExpression::extract_predicate(Solver *s) {

  //std::cout << "extract predicate " << children[0] << " * " << children[1] << " = " << children[2] << std::endl;


  s->add(Constraint(new PredicateMul(children)));
}


// Mistral::DivExpression::DivExpression(Variable X, Variable Y) 
//   : Expression(X,Y) {
// }
// Mistral::DivExpression::~DivExpression() {
//#ifdef _DEBUG_MEMORY
//  std::cout << "c delete div expression" << std::endl;
//#endif
//}

// void Mistral::DivExpression::extract_constraint(Solver *s) {
//   std::cerr << "Error: Div predicate can't be used as a constraint" << std::endl;
//   exit(0);
// }

// void Mistral::DivExpression::extract_variable(Solver *s) {
//   // z = x/y

//   int max_pos_x = children[0].get_max();
//   int min_neg_x = children[0].get_min();

//   int min_pos_x = children[0].get_min_pos();
//   int max_neg_x = children[0].get_max_neg();

//   int max_pos_y = children[1].get_max();
//   int min_neg_y = children[1].get_min();

//   int min_pos_y = children[1].get_min_pos();
//   int max_neg_y = children[1].get_max_neg();

//   // positive part of z's domain
//   int nlb1, nlb2, nub1, nub2;

//   nlb1 = nlb2 = INFTY; //lb_neg;
//   nub1 = nub2 = 0; //ub_neg;
	
//   // it can either be the positive parts of X and Y:
//   if(max_pos_x>0 && max_pos_y>0) {
//     // compute the bounds
//     //std::cout << "\t   lb = " << min_pos_x << "/" << max_pos_y << std::endl;
//     //std::cout << "\t   ub = " << max_pos_x << "/" << min_pos_y << std::endl;
//     // nlb1 = (int)(ceil((double)(min_pos_x)/(double)(max_pos_y)));
//     // nub1 = (int)(floor((double)(max_pos_x)/(double)(min_pos_y)));
//     nlb1 = min_pos_x/max_pos_y;
//     nub1 = max_pos_x/min_pos_y;
//   }

//   // or the negative parts of X and Y:
//   if(min_neg_x<0 && min_neg_y<0) {
//     // compute the bounds
//     // nlb2 = (int)(ceil((double)(max_neg_x)/(double)(min_neg_y)));
//     // nub2 = (int)(floor((double)(min_neg_x)/(double)(max_neg_y)));
//     nlb2 = max_neg_x/min_neg_y;
//     nub2 = min_neg_x/max_neg_y;
//   }
//   if(nlb1>nlb2) nlb1 = nlb2;
//   if(nub1<nub2) nub1 = nub2;
  
//   int lb_pos = nlb1;
//   int ub_pos = nub1;

  
//   nlb1 = nlb2 = 0; //lb_pos;
//   nub1 = nub2 = -INFTY; //ub_pos;
  
//   // it can either be the negative part of X and the positive part of Y:
//   if(min_neg_x<0 && max_pos_y>0) {
//     // compute the bounds  
//     // nlb1 = (int)(ceil((double)(min_neg_x)/(double)(min_pos_y)));
//     // nub1 = (int)(floor((double)(max_neg_x)/(double)(max_pos_y)));
//     nlb1 = min_neg_x/min_pos_y;
//     nub1 = max_neg_x/max_pos_y;
//   }
//   // or the negative part of Y and the positive part of X:
//   if(max_pos_x>0 && min_neg_y<0) {
//     // compute the bounds
//     // nlb2 = (int)(ceil((double)(max_pos_x)/(double)(max_neg_y)));
//     // nub2 = (int)(floor((double)(min_pos_x)/(double)(min_neg_y)));
//     nlb2 = max_pos_x/max_neg_y;
//     nub2 = min_pos_x/min_neg_y;
//   }
  
//   if(nlb1>nlb2) nlb1 = nlb2;
//   if(nub1<nub2) nub1 = nub2;

//   int lb_neg = nlb1;
//   int ub_neg = nub1;

//   //std::cout << "[" << lb_neg <<".." << ub_neg << "] u [" << lb_pos << ".." << ub_pos << "]" << std::endl;

//   if(lb_neg>ub_neg)
//     lb_neg = lb_pos;
//   else if(lb_pos>ub_pos)
//     ub_pos = ub_neg;

//   if(children[0].contain(0)) {
//     if(lb_neg > 0) lb_neg = 0;
//     if(ub_pos < 0) ub_pos = 0;
//   }


//   Variable aux(lb_neg, ub_pos, DYN_VAR);
//   _self = aux;

//   _self.initialise(s, 1);
//   _self = _self.get_var();
//   children.add(_self);
// }

// const char* Mistral::DivExpression::get_name() const {
//   return "div";
// }

// void Mistral::DivExpression::extract_predicate(Solver *s) {
//   VarArray scope;
//   scope.add(children[1]);
//   scope.add(children[2]);
//   scope.add(children[0]);

//   s->add(new PredicateMul(children));
//   children[1].remove(0);
// }


Mistral::FactorExpression::FactorExpression(Variable X, const int fct) 
  : Expression(X) { 
  factor=fct; 
}
Mistral::FactorExpression::~FactorExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete factor expression" << std::endl;
#endif
}
  
void Mistral::FactorExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Factor predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::FactorExpression::extract_variable(Solver *s) {

  int lb = (factor<0 ? children[0].get_max() : children[0].get_min())*factor;
  int ub = (factor<0 ? children[0].get_min() : children[0].get_max())*factor;

  Variable aux(lb, ub, DYN_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::FactorExpression::get_name() const {
  return "factor";
}

void Mistral::FactorExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateFactor(children, factor)));
}

Mistral::Variable Mistral::Variable::operator*(Variable x) {

  // std::cout << "ADD " ;
  // display(std::cout);
  // std::cout << " * " << x << std::endl;

  Variable exp(new MulExpression(*this,x));
  return exp;
}

Mistral::Variable Mistral::Variable::operator*(int k) {
  Variable exp(new FactorExpression(*this,k));
  return exp;
}


// Mistral::Variable Mistral::Variable::operator/(Variable x) {
//   Variable exp(new DivExpression(*this,x));
//   return exp;
// }

// Mistral::Variable Mistral::Variable::operator/(int k) {
//   Variable exp(new QuotientExpression(*this,k));
//   return exp;
// }


  Mistral::SubExpression::SubExpression(Variable X, Variable Y) 
  : Expression(X,Y) {
}
  Mistral::SubExpression::~SubExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete sub expression" << std::endl;
#endif
}
  
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
    _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
 
}

  void Mistral::SubExpression::extract_predicate(Solver *s) {
    VarArray tmp;
    for(int i=3; i;) tmp.add(children[--i]);

    Constraint sub(new PredicateAdd(tmp));
    
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
Mistral::NotExpression::~NotExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete not expression" << std::endl;
#endif
}
  
void Mistral::NotExpression::extract_constraint(Solver *s) {
  children[0].remove(0);
  // std::cerr << "Error: Not predicate can't be used as a constraint" << std::endl;
  // exit(0);
}

void Mistral::NotExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::NotExpression::get_name() const {
  return "not";
}

void Mistral::NotExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateNot(children)));
}

Mistral::Variable Mistral::Variable::operator!() {
  Variable exp(new NotExpression(*this));
  return exp;
}


Mistral::NegExpression::NegExpression(Variable X) 
  : Expression(X) { 
}
Mistral::NegExpression::~NegExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete neg expression" << std::endl;
#endif
}
  
void Mistral::NegExpression::extract_constraint(Solver *s) {
  //children[0].remove(0);
  std::cerr << "Error: Neg predicate can't be used as a constraint" << std::endl;
  // exit(0);
}

void Mistral::NegExpression::extract_variable(Solver *s) {
  Variable aux(-children[0].get_max(), -children[0].get_min(), DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::NegExpression::get_name() const {
  return "neg";
}

void Mistral::NegExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateNeg(children)));
}

Mistral::Variable Mistral::Variable::operator-() {
  Variable exp(new NegExpression(*this));
  return exp;
}

Mistral::AndExpression::AndExpression(Variable X, Variable Y) 
  : Expression(X,Y) {}
Mistral::AndExpression::~AndExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete and expression" << std::endl;
#endif
}
  
void Mistral::AndExpression::extract_constraint(Solver *s) {
  //s->add(new ConstraintAnd(children));
  //std::cout << "not implemented" << std::endl;
  //exit(1);
  children[0].remove(0);
  children[1].remove(0);
}

void Mistral::AndExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::AndExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateAnd(children)));
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
Mistral::OrExpression::~OrExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete or expression" << std::endl;
#endif
}
  
void Mistral::OrExpression::extract_constraint(Solver *s) {
  s->add(Constraint(new ConstraintOr(children)));
  //s->add(ConstraintOr::ConstraintOr_new(children));
  //std::cout << "not implemented" << std::endl;
  //exit(1);
}

void Mistral::OrExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::OrExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateOr(children)));
}

const char* Mistral::OrExpression::get_name() const {
  return "or";
}

Mistral::Variable Mistral::Variable::operator||(Variable x) {
  Variable exp(new OrExpression(*this,x));
  return exp;
}




// //   Mistral::NeqExpression::NeqExpression(Variable X, Variable Y) 
// //   : BinaryExpression(X,Y) {}
// //   Mistral::NeqExpression::~NeqExpression() {}
  
// // void Mistral::NeqExpression::extract_constraint(Solver *s) {
// // //     Constraint *neq = new ConstraintNotEqual(children);
// // //     neq->initialise();
// // //     return neq;
// //   s->add(new ConstraintNotEqual(children));
// //   }

// //   void Mistral::NeqExpression::extract_variable(Solver *s) {
// //     Variable aux(0, 1, BOOL_VAR);
// //     _self = aux;

// //     children.add(_self);
// //     _self.initialise(s, 1);
// //   }

// //   void Mistral::NeqExpression::extract_predicate(Solver *s) {
// // //     Constraint *neq = new PredicateEqual(children, false);
// // //     neq->initialise();
// // //     return neq;
// //     s->add(new PredicateEqual(children, false));
// //   }

// // const char* Mistral::NeqExpression::get_name() const {
// //   return "neq";
// // }



Mistral::EqualExpression::EqualExpression(Variable X, Variable Y, const int sp) 
  : Expression(X,Y) { spin=sp; value=NOVAL; }
Mistral::EqualExpression::EqualExpression(Variable X, const int y, const int sp) 
  : Expression() { children.add(X); value=y; spin=sp; }
Mistral::EqualExpression::~EqualExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete equal expression" << std::endl;
#endif
}

void Mistral::EqualExpression::extract_constraint(Solver *s) {

  //std::cout << "extract constraint from expression equal" << std::endl;
  
  if(spin) {
    
    // std::cout << "EQUAL " 
    // 	      << (children[0].domain_type == EXPRESSION) << "/" << (children[0].expression->is_set()) 
    // 	      << " " 
    // 	      << (children[1].domain_type == EXPRESSION && children[1].expression->is_set())
    // 	      <<std::endl;
    
    // exit(1);
    //#ifdef _COMPLETE
    if(children[0].domain_type == EXPRESSION && children[0].expression->is_set()
       && children.size>1 &&
       children[1].domain_type == EXPRESSION && children[1].expression->is_set()
       ) {
      
      unsigned int i = 0, j = 0;
      SetExpression *x = (SetExpression*)(children[0].expression);
      SetExpression *y = (SetExpression*)(children[1].expression);
      Variable bx;
      Variable by;
      
      s->add(new ConstraintEqual(children));
      
      Vector<Variable> scp;

      //std::cout << x->elts_ub << " == " << y->elts_ub << std::endl;
      
      
      while(i<x->elts_ub.size && j<y->elts_ub.size) {
	
	// std::cout << i << " - " << j << std::endl;
	
	// std::cout << x->get_index_var(i) << std::endl;

	// std::cout << y->get_index_var(j) << std::endl;

	scp.clear();
	bx = x->get_index_var(i);
	by = y->get_index_var(j);
	if(x->elts_ub[i] == y->elts_ub[j]) {
	  if(bx.is_ground()) {
	    by.set_domain(bx.get_min());
	  } else if(by.is_ground()) {
	    bx.set_domain(by.get_min());
	  } else {
	    scp.add(by);
	    scp.add(bx);
	    s->add(new ConstraintEqual(scp));
	  }
	  ++i;
	  ++j;
	} else if(x->elts_ub[i] <= y->elts_ub[j]) {
	  // x->elts_ub[i] can't be in Y
	  bx.set_domain(0);
	  ++i;
	} else {
	  // y->elts_ub[j] can't be in X
	  by.set_domain(0);
	  ++j;
	}
      }
    
      // remaining values of X
      while(i<x->elts_ub.size) {
	bx = x->get_index_var(i);
	bx.set_domain(0);
	++i;
      }
      
      // remaining values of Y
      while(j<y->elts_ub.size) {
	by = y->get_index_var(j);
	by.set_domain(0);
	++j;
      }
    }
    else

      // #endif

      if(children.size==2)
    	  s->add(Constraint(new ConstraintEqual(children), (BINARY|IDEMPOTENT)));
      else
    	  {
	    if(children[0].set_domain(value) == FAIL_EVENT)
	      s->fail();
    	  }

    //#ifdef _COMPLETE
 
  } else {
    if(children[0].domain_type == EXPRESSION && children[0].expression->is_set()
       && children.size>1 &&
       children[1].domain_type == EXPRESSION && children[1].expression->is_set()
       ) {
      unsigned int i = 0, j = 0;
      SetExpression *x = (SetExpression*)(children[0].expression);
      SetExpression *y = (SetExpression*)(children[1].expression);
      Variable bx;
      Variable by;

      Vector<Variable> scp;
      Vector<Variable> disjunction;

      while(i<x->elts_ub.size && j<y->elts_ub.size) {
	Variable b(0,1);

	scp.clear();
	bx = x->get_index_var(i);
	by = y->get_index_var(j);
	if(x->elts_ub[i] == y->elts_ub[j]) {
	  if(bx.is_ground()) {
	    scp.add(by);
	    scp.add(b);
	    s->add(Constraint(new PredicateConstantEqual(scp,bx.get_min())));
	    disjunction.add(b);
	  } else if(by.is_ground()) {
	    scp.add(bx);
	    scp.add(b);
	    s->add(Constraint(new PredicateConstantEqual(scp,by.get_min())));
	    disjunction.add(b);
	  } else {
	    scp.add(by);
	    scp.add(bx);
	    scp.add(b);
	    s->add(Constraint(new PredicateConstantEqual(scp)));
	    disjunction.add(b);
	  }
	  ++i;
	  ++j;
	} else if(x->elts_ub[i] <= y->elts_ub[j]) {
	  // x->elts_ub[i] can't be in Y
	  scp.add(bx);
	  scp.add(b);
	  s->add(Constraint(new PredicateConstantEqual(scp,0)));
	  disjunction.add(b);
	  ++i;
	} else {
	  // y->elts_ub[j] can't be in X
	  scp.add(by);
	  scp.add(b);
	  s->add(Constraint(new PredicateConstantEqual(scp,0)));
	  disjunction.add(b);
	  ++j;
	}
      }
    
      // remaining values of X
      while(i<x->elts_ub.size) {
	Variable b(0,1);
	bx = x->get_index_var(i);
	scp.add(bx);
	scp.add(b);
	s->add(Constraint(new PredicateConstantEqual(scp,0)));
	disjunction.add(b);
	++i;
      }
      
      // remaining values of Y
      while(j<y->elts_ub.size) {
	Variable b(0,1);
	by = y->get_index_var(j);
	scp.add(by);
	scp.add(b);
	s->add(Constraint(new PredicateConstantEqual(scp,0)));
	disjunction.add(b);
	++j;
      }

      s->add(Constraint(new ConstraintBoolSumInterval(disjunction,0,disjunction.size-1)));
    }
    else 

      //#endif
      
      if(children.size==2) {
	
	//std::cout << "not equal " << children << std::endl;

	s->add(Constraint(new ConstraintNotEqual(children)));

	//exit(1);

	
      } else if(children[0].remove(value) == FAIL_EVENT) s->fail();

  }
}


// void Mistral::EqualExpression::extract_constraint(Solver *s) {
//   if(spin) {
    
//     if(children.size==2) s->add(Constraint(new ConstraintEqual(children)// , (BINARY|IDEMPOTENT)
// 					   ));
//     else children[0].set_domain(value);
    
//   } else {
    
//     if(children.size==2) {
//       s->add(Constraint(new ConstraintNotEqual(children)// , (BINARY|IDEMPOTENT)
// 			));
//     } else children[0].remove(value);
//   }
// }

void Mistral::EqualExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::EqualExpression::extract_predicate(Solver *s) {
  
  if(children.size==3) {

    //#ifdef _COMPLETE

    if(children[0].domain_type == EXPRESSION && children[0].expression->is_set()) {
      
      s->add(Constraint(new PredicateEqual(children, 1)));

      unsigned int i = 0, j = 0;
      SetExpression *x = (SetExpression*)(children[0].expression);
      SetExpression *y = (SetExpression*)(children[1].expression);
      Variable bx;
      Variable by;

      Vector<Variable> aux;
      Vector<Variable> scp;

      while(i<x->elts_ub.size && j<y->elts_ub.size) {
	scp.clear();
	Variable b(0,1);
	b.initialise(s, 1);
	aux.add(b);
	scp.add(children[2]);
	scp.add(aux.back());
	s->add(new ConstraintLess(scp));
	scp.clear();
	if(x->elts_ub[i] == y->elts_ub[j]) {
	  bx = x->get_index_var(i);
	  by = y->get_index_var(j);
	  if(bx.is_ground()) {
	    scp.add(by);
	    scp.add(aux.back());
	    s->add(new PredicateConstantEqual(scp, bx.get_min(), 1));
	  } else if(by.is_ground()) {
	    scp.add(bx);
	    scp.add(aux.back());
	    s->add(new PredicateConstantEqual(scp, by.get_min(), 1));
	  } else {
	    scp.add(bx);
	    scp.add(by);
	    scp.add(aux.back());
	    s->add(new PredicateEqual(scp, 1));
	  }
	  ++i;
	  ++j;
	} else if(x->elts_ub[i] <= y->elts_ub[j]) {
	  // x->elts_ub[i] can't be in Y
	  scp.add(by);
	  scp.add(aux.back());
	  s->add(new PredicateConstantEqual(scp, 0, 1));
	  ++i;
	} else {
	  // y->elts_ub[j] can't be in X
	  scp.add(bx);
	  scp.add(aux.back());
	  s->add(new PredicateConstantEqual(scp, 0, 1));
	  ++j;
	}
      }
    
      // remaining values of X
      while(i<x->elts_ub.size) {
	
	scp.add(by);
	scp.add(aux.back());
	s->add(new PredicateConstantEqual(scp, 0, 1));
	++i;
      }
      
      // remaining values of Y
      while(j<y->elts_ub.size) {
	
	scp.add(bx);
	scp.add(aux.back());
	s->add(new PredicateConstantEqual(scp, 0, 1));
	++j;
      }

      //s->add( Sum() )

    } else {

      s->add(Constraint(new PredicateEqual(children, spin)));

    }

  } else s->add(Constraint(new PredicateConstantEqual(children, value, spin)));
}

const char* Mistral::EqualExpression::get_name() const {
  return "equal";
}

Mistral::EqualSetExpression::EqualSetExpression(Variable X, Variable Y, const int sp) 
  : Expression(X,Y) { spin=sp; value=NOVAL; }
Mistral::EqualSetExpression::EqualSetExpression(Variable X, const int y, const int sp) 
  : Expression() { children.add(X); value=y; spin=sp; }
Mistral::EqualSetExpression::~EqualSetExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete set expression" << std::endl;
#endif
}

void Mistral::EqualSetExpression::extract_constraint(Solver *s) {
  if(spin) {
    if(children.size==2) s->add(new ConstraintEqual(children));
    else children[0].set_domain(value);
  } else {
    if(children.size==2) s->add(new ConstraintNotEqual(children));
    else children[0].remove(value);
  }
}

void Mistral::EqualSetExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::EqualSetExpression::extract_predicate(Solver *s) {

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

const char* Mistral::EqualSetExpression::get_name() const {
  return "equal";
}

Mistral::Variable Mistral::Variable::operator==(Variable x) {
  // Variable exp;
  // if(domain_type == EXPRESSION && is_set()) {
  //   exp = Variable(new EqualSetExpression(*this,x,1));
  // } else {
  //   exp = Variable(new EqualExpression(*this,x,1));
  // }
  // return exp;
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
Mistral::PrecedenceExpression::~PrecedenceExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete precedence expression" << std::endl;
#endif
}
  
void Mistral::PrecedenceExpression::extract_constraint(Solver *s) {
  if(children.size==2) {
    s->add(Constraint(new ConstraintLess(children, offset)// , (BINARY|IDEMPOTENT)
		      ) );
  }
  else {
    if(spin) {
      
      //std::cout << "HERE" << std::endl;

      if(children[0].set_max(offset) == FAIL_EVENT) 
	{ s->fail(); }
    } else {     
      if(children[0].set_min(offset) == FAIL_EVENT)
	{ s->fail(); }
    }
  }
}

void Mistral::PrecedenceExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::PrecedenceExpression::extract_predicate(Solver *s) {
  if(children.size==3) {
    s->add(Constraint(new PredicateLess(children, offset)));
  } else if(spin) {
    s->add(Constraint(new PredicateUpperBound(children, offset)));
  } else {
    s->add(Constraint(new PredicateLowerBound(children, offset)));
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
Mistral::DisjunctiveExpression::~DisjunctiveExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete disjunctive expression" << std::endl;
#endif
}
  
void Mistral::DisjunctiveExpression::extract_constraint(Solver *s) {
  s->add(Constraint(new ConstraintDisjunctive(children, processing_time[0], processing_time[1])));
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
  return "disjunct";
}

Mistral::Variable Mistral::Disjunctive(Variable X, Variable Y, 
				       const int px, const int py) 
{
  Variable exp(new DisjunctiveExpression(X,Y,px,py));
  return exp;
}



Mistral::ReifiedDisjunctiveExpression::ReifiedDisjunctiveExpression(Variable X, Variable Y, 
								    const int px, const int py) 
  : Expression(X,Y) { processing_time[0] = px; processing_time[1] = py; }

Mistral::ReifiedDisjunctiveExpression::~ReifiedDisjunctiveExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete reified disjunctive expression" << std::endl;
#endif
}
  
void Mistral::ReifiedDisjunctiveExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: ReifiedDisjunctive constraint can't yet be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::ReifiedDisjunctiveExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::ReifiedDisjunctiveExpression::extract_predicate(Solver *s) {
  //s->add(new ConstraintTernaryDisjunctive(children, processing_time[0], processing_time[1]));
  s->add(Constraint(new ConstraintReifiedDisjunctive(children, processing_time[0], processing_time[1])));
}

const char* Mistral::ReifiedDisjunctiveExpression::get_name() const {
  return "r-disjunct";
}

Mistral::Variable Mistral::ReifiedDisjunctive(Variable X, Variable Y, 
					      const int px, const int py) 
{
  Variable exp(new ReifiedDisjunctiveExpression(X,Y,px,py));
  return exp;
}



Mistral::FreeExpression::FreeExpression(Variable X) 
  : Expression(X) { };

Mistral::FreeExpression::~FreeExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete free expression" << std::endl;
#endif
}
  
const char* Mistral::FreeExpression::get_name() const {
  return "free";
}

Mistral::Variable Mistral::Free(Variable X) 
{
  Variable exp(new FreeExpression(X));
  return exp;
}


Mistral::AllDiffExpression::AllDiffExpression(Vector< Variable >& args, const int ct) 
  : Expression(args) { consistency_level = ct; }

Mistral::AllDiffExpression::~AllDiffExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete alldiff expression" << std::endl;
#endif
}

void Mistral::AllDiffExpression::extract_constraint(Solver *s) { 
//   Constraint *con = new ConstraintAllDiff(children); 
//   con->initialise();
//   return con;
  if(consistency_level == BOUND_CONSISTENCY)
    s->add(Constraint(new ConstraintAllDiff(children))); 
  s->add(Constraint(new ConstraintCliqueNotEqual(children))); 
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


Mistral::LexExpression::LexExpression(Vector< Variable >& r1, Vector< Variable >& r2, const int st_)
  : Expression() { 
  int row_size = r1.size;
  for(int i=0; i<row_size; ++i)
    children.add(r1[i]);
  for(int i=0; i<row_size; ++i)
    children.add(r2[i]);
  for(int i=0; i<=row_size; ++i)
    children.add( Variable(0,1) );
  strict = st_; 
}

Mistral::LexExpression::~LexExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete lex expression" << std::endl;
#endif
}

void Mistral::LexExpression::extract_constraint(Solver *s) { 
  int arity = (children.size-1)/3;
  VarArray scp;
  // scp.add(children[0]);
  // scp.add(children[arity]);
  // scp.add(children[2*arity]);

  // s->add( new ConstraintLexf(scp) );

  for(int i=0; i<arity; ++i) {
    scp.clear();
    
    scp.add(children[i]);
    scp.add(children[i+arity]);
    scp.add(children[i+2*arity]);
    scp.add(children[i+2*arity+1]);

    s->add( new ConstraintLex(scp) );
  }

  if(children[2*arity].set_domain(0) == FAIL_EVENT)
    { s->fail(); }
  if(strict) 
    if(children[3*arity].set_domain(1) == FAIL_EVENT)
      { s->fail(); }
}

void Mistral::LexExpression::extract_variable(Solver *s) {
      std::cerr << "Error: Lex constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

void Mistral::LexExpression::extract_predicate(Solver *s) { 
      std::cerr << "Error: Lex constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::LexExpression::get_name() const {
  return "lex";
}

Mistral::Variable Mistral::LexLeq(VarArray& r1, VarArray& r2) {
  Variable exp(new LexExpression(r1, r2, false));
  return exp;
}

Mistral::Variable Mistral::LexLess(VarArray& r1, VarArray& r2) {
  Variable exp(new LexExpression(r1, r2, true));
  return exp;
}

Mistral::Variable Mistral::VarArray::operator<(VarArray& X) {
  return LexLess(*this, X);
}

Mistral::Variable Mistral::VarArray::operator>(VarArray& X) {
  return LexLess(X, *this);
}

Mistral::Variable Mistral::VarArray::operator<=(VarArray& X) {
  return LexLeq(*this, X);
}

Mistral::Variable Mistral::VarArray::operator>=(VarArray& X) {
  return LexLeq(X, *this);
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

Mistral::BoolSumExpression::~BoolSumExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete boolsum expression" << std::endl;
#endif
}

void Mistral::BoolSumExpression::extract_constraint(Solver *s) { 
  s->add(new ConstraintBoolSumInterval(children,lb,ub)); 
}

void Mistral::BoolSumExpression::extract_variable(Solver *s) {
  Variable aux(lb, ub, DYN_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  //children.add(_self);
}

void Mistral::BoolSumExpression::extract_predicate(Solver *s) { 
  //s->add(new PredicateBoolSum(children, self)); 
  s->add(new PredicateBoolSum(children, s->variables[id])); 
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
Mistral::LinearExpression::~LinearExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete linear expression" << std::endl;
#endif
}
  
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
  if(lower_bound == 0 &&
     upper_bound == 0 && 
     children.size == 3 && 
     abs(weight[0]) == 1 &&
     abs(weight[1]) == 1 &&
     abs(weight[2]) == 1
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
      s->add(Constraint(new PredicateAdd(children)));
    }  
  }
  
  if(!post_add)
    s->add(Constraint(new PredicateWeightedSum(children, weight, lower_bound, upper_bound)));
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
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
  weight.add(-1);
}

const char* Mistral::LinearExpression::get_name() const {
  return "linear_sum";
}

void Mistral::LinearExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateWeightedSum(children, weight, 0, 0)));
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


Mistral::MinExpression::MinExpression(Vector< Variable >& args) 
  : Expression(args) {}

Mistral::MinExpression::~MinExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete min expression" << std::endl;
#endif
}

void Mistral::MinExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Min predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::MinExpression::extract_variable(Solver *s) {
  int lower_bound = children[0].get_min();
  int upper_bound = children[0].get_max();

  int arity = children.size;
  for(int i=1; i<arity; ++i) {
    if(children[i].get_min() < lower_bound) lower_bound = children[i].get_min();
    if(children[i].get_max() < upper_bound) upper_bound = children[i].get_max();
  }

  Variable aux(lower_bound, upper_bound, DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::MinExpression::get_name() const {
  return "min";
}

void Mistral::MinExpression::extract_predicate(Solver *s) {
  s->add(new PredicateMin(children));
}

Mistral::Variable Mistral::Min(Vector<Variable>& X) {
  Variable exp( new MinExpression(X) );
  return exp;
}

Mistral::Variable Mistral::Min(VarArray& X) {
  Variable exp( new MinExpression(X) );
  return exp;
}

Mistral::Variable Mistral::Min(Variable X, Variable Y) {
  MinExpression *mexp = new MinExpression();
  mexp->add(X);
  mexp->add(Y);
  Variable exp( mexp );
  return exp;
}



Mistral::Variable Mistral::Max(Vector<Variable>& X) {
  Variable exp( new MaxExpression(X) );
  return exp;
}

Mistral::Variable Mistral::Max(VarArray& X) {
  Variable exp( new MaxExpression(X) );
  return exp;
}

Mistral::Variable Mistral::Max(Variable X, Variable Y) {
  MaxExpression *mexp = new MaxExpression();
  mexp->add(X);
  mexp->add(Y);
  Variable exp( mexp );
  return exp;
}


Mistral::MaxExpression::MaxExpression(Vector< Variable >& args) 
  : Expression(args) {}

Mistral::MaxExpression::~MaxExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete max expression" << std::endl;
#endif
}

void Mistral::MaxExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Max predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::MaxExpression::extract_variable(Solver *s) {
  int lower_bound = children[0].get_min();
  int upper_bound = children[0].get_max();

  int arity = children.size;
  for(int i=1; i<arity; ++i) {
    if(children[i].get_min() > lower_bound) lower_bound = children[i].get_min();
    if(children[i].get_max() > upper_bound) upper_bound = children[i].get_max();
  }

  Variable aux(lower_bound, upper_bound, DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::MaxExpression::get_name() const {
  return "max";
}

void Mistral::MaxExpression::extract_predicate(Solver *s) {
  s->add(new PredicateMax(children));
}


Mistral::ElementExpression::ElementExpression(Vector< Variable >& args, 
					      Variable X, int ofs) 
  : Expression(), offset(ofs) {
  for(unsigned int i=0; i<args.size; ++i) children.add(args[i]);
  children.add(X);
}

Mistral::ElementExpression::~ElementExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete element expression" << std::endl;
#endif
}


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
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::ElementExpression::get_name() const {
  return "element";
}

void Mistral::ElementExpression::extract_predicate(Solver *s) {
  int arity = children.size-2;
  if(children[arity].set_min(offset) == FAIL_EVENT)
    { s->fail(); }
  if(children[arity].set_max(arity-1+offset) == FAIL_EVENT)
    { s->fail(); }
  s->add(Constraint(new PredicateElement(children, offset)));
}


Mistral::ElementSetExpression::ElementSetExpression(Vector< Variable >& args, 
						    Variable X, int ofs) 
  : SetExpression(), offset(ofs) {
  for(unsigned int i=0; i<args.size; ++i) children.add(args[i]);
  children.add(X);
  num_args = children.size;
  initialise_domain();
  initialise_elements();
}

Mistral::ElementSetExpression::~ElementSetExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete element set expression" << std::endl;
#endif
}


void Mistral::ElementSetExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: ElementSet predicate can't be used as a constraint" << std::endl;
  exit(0);
}

// Mistral::Variable Mistral::ElementSetExpression::get_index_var(const int idx) {
//   ((SetExpression*)(_self.variable))->children[i];
// }

void Mistral::ElementSetExpression::initialise_domain() {

  // lb and ub stand for the cardinality's bounds

  // val_ub and val_lb stand for the values' bounds
  int val_lb = INFTY;
  int val_ub = -INFTY;
  lb = INFTY;
  ub = -INFTY;

  int arity = num_args-1;


  // first, we go through the children to get the cardinality's and values' bounds

  SetExpression *S;
  int i, // j,
    nxt = children[arity].get_min(), bound;
  unsigned int j;
  do {
    i = nxt-offset;
    
    if(i>=0 && i<arity) {
      S = (SetExpression*)(children[i].variable);

      bound = S->elts_ub.front();
      if(bound<val_lb) val_lb = bound;

      bound = S->elts_ub.back();
      if(bound>val_ub) val_ub = bound;

      bound = S->lb;
      if(bound<lb) lb = bound;

      bound = S->ub;
      if(bound>ub) ub = bound;
    }

    nxt = children[arity].next(i+offset);
  } while(i+offset<nxt);
  


  // then we go once more over the children to get the lower and upper bounds of the set var
  BitSet lb_domain(0,arity-1,BitSet::full);
  BitSet ub_domain(0,arity-1,BitSet::empt);
  BitSet aux(0,arity-1,BitSet::empt);
  nxt = children[arity].get_min();

  do {
    i = nxt-offset;
    
    if(i>=0 && i<arity) {
      S = (SetExpression*)(children[i].variable);
      for(j=0; j<S->elts_ub.size; ++j) {
	ub_domain.add(S->elts_ub[j]);
      }

      for(j=0; j<S->elts_lb.size; ++j) {
	aux.add(S->elts_lb[j]);
      }

      lb_domain.intersect_with(aux);
      aux.clear();
    }

    nxt = children[arity].next(i+offset);
  } while(i+offset<nxt);


  //std::cout << lb_domain << " << S << " << ub_domain << std::endl;

  i = ub_domain.min();
  //j = ub_domain.max();

  do {
    elts_ub.add(i);
    nxt = i;
    i = ub_domain.next(nxt);
  } while(nxt<i);

  
 //  for(int elt=0; elt<arity; ++elt) { // for each element 
 //    nxt = children[arity].get_min();
 //  do {
 //    i = nxt-offset;
    
 //    if(i>=0 && i<arity) {
 //      if(children[i].get_min()==0) {
 // 	ub_domain.add()
 //      }
 // < lower_bound) lower_bound = children[i].get_min();
 //      if(children[i].get_max() > upper_bound) upper_bound = children[i].get_max();
 //    }

 //    nxt = children[arity].next(i+offset);
 //  } while(i+offset<nxt);
  

 //  domain.initialise(lower_bound, upper_bound, BitSet::empt);
  
 //  nxt = children[arity].get_min();
 //  do {
 //    i = nxt-offset;
    
 //    if(i>=0 && i<arity) {
 //      children[i].union_to(domain);
 //    }

 //    nxt = children[arity].next(i+offset);
 //  } while(i+offset<nxt);


 //  nxt = lower_bound;
 //  do {
 //    i = nxt;
 //    values.add(i);
 //    nxt = domain.next(i);
 //  } while(i<nxt);
  
}

// void Mistral::ElementSetExpression::extract_variable(Solver *s) {
//   initialise_domain();

//   Variable aux(new SetExpression(elts_lb, elts_ub, lb, ub));
//   _self = aux;

//   _self.initialise(s, 1);
//   //_self = _self.get_var();
//   //children.add(_self);
// }

const char* Mistral::ElementSetExpression::get_name() const {
  return "set_element";
}

void Mistral::ElementSetExpression::extract_predicate(Solver *s) {
 int arity = num_args-1;
 int i;
 unsigned int j, k;

 // std::cout << "HERE " << std::endl;
 // std::cout << children[arity] << " in " << children[arity].get_domain() << std::endl;

 //int nxt = children[arity].get_min();
 int lb_index = children[arity].get_min();
 int ub_index = children[arity].get_max();
 SetExpression *S;
 
 //VarArray scp[ub_index-lb_index+1];
 Vector<Variable> scp[elts_ub.size];
 
 for(i=lb_index; i<=ub_index; ++i) {
   if(i-offset>=0 && i-offset<arity) {
     if(children[arity].contain(i)) {

       //std::cout << children[i-offset] << std::endl;

       S = (SetExpression*)(children[i-offset].variable);

	   // std::cout << S->children << std::endl;
	   // std::cout << S->elts_ub << std::endl;
	   // std::cout << S->num_args << std::endl;


       k=0;
       for(j=0; j<elts_ub.size; ++j) {

	 //   std::cout << k << std::endl;
	 
	 // int hh = elts_ub[j];
	 // int jj = S->elts_ub[k];

	 // if(hh == jj) {

	 if(k<S->elts_ub.size && elts_ub[j] == S->elts_ub[k]) {
	   scp[j].add(S->get_index_var(k++));
	 } else {
	   Variable x(0);
	   scp[j].add(x);
	 }
       }
     } else {
       for(j=0; j<elts_ub.size; ++j) {
	 Variable x(0);
	 scp[j].add(x);
       }
     }
   }
 }

 
 for(j=0; j<elts_ub.size; ++j) {
   //std::cout << scp[j] << std::endl;

   scp[j].add(children[arity]);
   scp[j].add(get_index_var(j));
   
   s->add(new PredicateElement(scp[j], offset));
   
   
   //std::cout << std::endl << std::endl << s << std::endl;
   
   
   // s->add(Element(scp[j], children[arity].get_var(), offset) == 
   // 	   ((SetExpression*)(_self.variable))->children[j].get_var());
 }
 
 scp[0].clear();
 for(i=0; i<=arity; ++i) {
   scp[0].add(children[i]);
 }

 //scp[0].add(self);
 scp[0].add(s->variables[id]);
 s->add(new PredicateElement(scp[0], offset));

  // do {
  //   i = nxt-offset;

  //   if(i>=0 && i<arity) {
  //     S = (SetExpression*)(children[i].variable);

  //     k=0;
  //     for(unsigned int j=0; j<elts_ub.size; ++j) {
  // 	if(elts_ub[j] == S->elts_ub[k]) {
	  
  // 	}
  //     }

  //   nxt = children[arity].next(i+offset);
  // } while(i+offset<nxt);

  //s->add(new PredicateElementSet(children, offset));
}

Mistral::Variable Mistral::VarArray::operator[](Variable X) {
  Variable exp( new ElementExpression(*this, X, 0) );
  return exp;
}

Mistral::Variable Mistral::VarArray::operator[](const int X) {
  return stack_[X];
}

// // Mistral::Variable& Mistral::VarArray::operator[](const int X) {
// //   return &(stack_[X]);
// // }

void Mistral::VarArray::set(const int X, Variable x) {
  stack_[X] = x;
}

Mistral::Variable Mistral::Element(Vector<Variable>& X, Variable selector, int offset) {
  Variable exp( new ElementExpression(X, selector, offset) );
  return exp;
}

Mistral::Variable Mistral::Element(VarArray& X, Variable selector, int offset) {
  Variable exp( new ElementExpression(X, selector, offset) );
  return exp;
}


Mistral::Variable Mistral::ElementSet(Vector<Variable>& X, Variable selector, int offset) {
  Variable exp( new ElementSetExpression(X, selector, offset) );
  return exp;
}

Mistral::Variable Mistral::ElementSet(VarArray& X, Variable selector, int offset) {
  Variable exp( new ElementSetExpression(X, selector, offset) );
  return exp;
}



Mistral::IntersectionExpression::IntersectionExpression(Variable X, Variable Y) 
  : SetExpression() {
  children.add(X);
  children.add(Y);
  num_args = 2;
  initialise_domain();
  initialise_elements();
}

Mistral::IntersectionExpression::~IntersectionExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete intersection expression" << std::endl;
#endif
}


void Mistral::IntersectionExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Intersection predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::IntersectionExpression::initialise_domain() { 
  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);

  unsigned int i=0;
  unsigned int j=0;

  int x, y;
  
  while(i<X->elts_ub.size && j<Y->elts_ub.size) {

    x = X->elts_ub[i];
    y = Y->elts_ub[j];
    if(x == y) {
      elts_ub.add(x);
      ++i;++j;
    } else if(x > y) { 
      ++j;
    } else {
      ++i;
    }
  }

  lb = 0;
  ub = elts_ub.size;
}

const char* Mistral::IntersectionExpression::get_name() const {
  return "set_intersection";
}

void Mistral::IntersectionExpression::extract_predicate(Solver *s) {
  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);

  unsigned int i=0;
  unsigned int j=0;
  unsigned int k=0;

  int x, y;

  VarArray scp;
  
  while(i<X->elts_ub.size && j<Y->elts_ub.size) {

    x = X->elts_ub[i];
    y = Y->elts_ub[j];
    if(x == y) {
      scp.add(X->get_index_var(i));
      scp.add(Y->get_index_var(j));
      scp.add(get_index_var(k++));

      s->add( Constraint(new PredicateAnd(scp)) );

      scp.clear();
      ++i;++j;
    } else if(x > y) { 
      ++j;
    } else {
      ++i;
    }
  }  

  //scp.add(self);
  scp.add(s->variables[id]);
  scp.add(X);
  s->add( Constraint(new ConstraintLess(scp)) );
 
  scp.clear();
 
  //scp.add(self);
  scp.add(s->variables[id]);
  scp.add(Y);
  s->add( Constraint(new ConstraintLess(scp)) );
 
}

Mistral::Variable Mistral::Intersection(Variable X, Variable Y) {
  Variable exp(new IntersectionExpression(X, Y));
  return exp;
}


Mistral::SetExpression::SetExpression(const int lelt, const int uelt, 
				      const int clb, const int cub) 
  : BoolSumExpression(clb, cub) {
  num_args = 0;
  elts_ub.initialise(0, uelt-lelt+1);
  for(int elt=lelt; elt<=uelt; ++elt) {
    elts_ub.add(elt);
    // Variable x(0, 1, BOOL_VAR);
    // children.add(x);
  }
  initialise_elements();
}


// Mistral::SetExpression::SetExpression(const BitSet& lb, const BitSet& ub, const int clb, const int cub);

Mistral::SetExpression::SetExpression(const Vector<int>& lb, const Vector<int>& ub, const int clb, const int cub)
  : BoolSumExpression(clb, cub) {
  // elts_ub.initialise(0, ub.back()-ub.front()+1);
  // unsigned int i=0, j=0;
  // while(i<ub.size) {
  //   elts_ub.add(ub[i]);
  //   if(lb.size>j && lb[j] == ub[i]) {
  //     elts_lb.add(lb[j++]);
  //     Variable c(1, 1, CONST_VAR);
  //     children.add(c);
  //   } else {
  //     Variable x(0, 1, BOOL_VAR);
  //     children.add(x);
  //   }
  //   ++i;
  // }
  num_args = 0;
  elts_ub.initialise(0, ub.back()-ub.front()+1);
  elts_ub.initialise(0, lb.back()-lb.front()+1);
  unsigned int i=0;
  while(i<ub.size) {
    elts_ub.add(ub[i++]);
  }

  i=0;
  while(i<lb.size) {
    elts_lb.add(lb[i++]);
  }
  initialise_elements();
}


void Mistral::SetExpression::initialise_elements() {
  unsigned int i=0, j=0;
  while(i<elts_ub.size) {
    if(elts_lb.size>j && elts_lb[j] == elts_ub[i]) {
      ++j;
      Variable c(1);
      children.add(c);
    } else {
      Variable x(0, 1);
      children.add(x);
    }
    ++i;
  }
}


Mistral::SetExpression::~SetExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete set expression" << std::endl;
#endif
}

const char* Mistral::SetExpression::get_name() const {
  return "set";
}

int Mistral::SetExpression::get_element_index(const int vali)  {
  int target, xlb = 0, xub = elts_ub.size-1, idx = -1;
  while(xlb <= xub) {

    //std::cout << vali << " [" << xlb << "," << xub << "] -> " ; //<< std::endl;

    target = (xlb+xub)/2;

    // std::cout << target << ": " << elts_ub[target] << " " << idx << " -- " <<
    // 	      elts_ub << std::endl;

    if(elts_ub[target] == vali) { idx=target; break; }
    else if(elts_ub[target] > vali) xub = target-1;
    else xlb = target+1;
  }
  //  std::cout << " ==> " << idx << std::endl; 
  return idx;
}

// Mistral::Variable Mistral::SetExpression::get_index_var(const int idx) {
//   return children[idx];
// }


Mistral::Variable Mistral::SetExpression::get_elt_var(const int vali) {
  // int target = 0, lb = -1, ub = elts_ub.size;
  // while(lb < ub) {
  //   target = lb+ub/2;
  //   if(elts_ub[target] == vali) break;
  //   else if(elts_ub[target] > vali) ub = target;
  //   else lb = target;
  // }
  // return children[target];
  int target = get_element_index(vali);
  if(target>=0) return children[target];
  Variable x;
  return x;
}

std::ostream& Mistral::SetExpression::display(std::ostream& os) const {
  // os << std::endl << elts_lb << std::endl
  //    << elts_ub << std::endl;
  // for(unsigned int i=0; i<elts_ub.size; ++i) {
  //   os << " " << children[i].get_domain();
  // }
  // os << std::endl;

  
  //std::cout << children << std::endl;


  // bool first=true;
  // os << "{" ;
  // //if(children[0].get_min()) os << elts_ub[0];
  // for(unsigned int i=0; i<elts_ub.size; ++i) {
  //   if(children[num_args+i].get_min()) {
  //     if(!first) os << ", " ;
  //     else first = false; 
  //     os << elts_ub[i];
  //   }
  // }
  // first=true;
  // os << "} <= S" << id << " <= {"; 
  //if(children[0].get_max()) os << elts_ub[0];
  // for(unsigned int i=0; i<elts_ub.size; ++i) {
  //   if(children[num_args+i].get_max()) {
  //     if(!first) os << ", " ;
  //     else first = false;
  //     os << elts_ub[i];
  //   }
  // }
  // os << "}";

  os << "S" << id;
  return os;
}

int Mistral::SetExpression::get_solution_int_value() const { 
  int i=elts_ub.size;
  int t = 0;
  int m = 1;
  while(--i>=0) {
    if((children[i].domain_type == CONST_VAR && children[i].constant_value == 1) ||
       ((Solver*)solver)->last_solution_lb[children[i].id()]) {
      t += m*elts_ub[i];
      m *= 10;
    }
  }
  return t;
} 

std::string Mistral::SetExpression::get_solution_str_value() const { 
  std::ostringstream  ret_str;

  ret_str << "{";
  bool first = true;

  for(unsigned int i=0; i<elts_ub.size; ++i) {
    //std::cout << i << ": " << elts_ub[i] << " " << children[i].id() << " " << solver->last_solution_lb[children[i].id()] << std::endl;

    if((children[i].domain_type == CONST_VAR && children[i].constant_value == 1) ||
       ((Solver*)solver)->last_solution_lb[children[i].id()]) {
      if(!first) {
	ret_str << ", " ;
      } else {
	first = false;
      }
      ret_str << elts_ub[i];
    }
  }

  ret_str << "}";

  return ret_str.str();
} 

Mistral::Variable Mistral::SetVariable(const int lelt, const int uelt, const int clb, const int cub) {
  int ub = cub;
  if(cub > uelt-lelt+1) ub = (uelt-lelt+1);
  SetExpression *se = new SetExpression(lelt, uelt, clb, ub);
  //std::cout << "create " << se << std::endl;
  Variable exp(se);
  return exp;
}

Mistral::Variable Mistral::SetVariable(const Vector<int>& slb, const Vector<int>& sub, 
				       const int clb, const int cub) {
  SetExpression *se = new SetExpression(slb, sub, clb, cub);
  Variable exp(se);
  return exp;
}

Mistral::Variable Mistral::Card(Variable S) { return S; }

Mistral::SubsetExpression::SubsetExpression(Variable X, Variable Y) 
  : Expression(X,Y) { 
}

Mistral::SubsetExpression::~SubsetExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete subset expression" << std::endl;
#endif
}

void Mistral::SubsetExpression::extract_constraint(Solver *s) { 
  unsigned int i = 0, j = 0;
  SetExpression *x = (SetExpression*)(children[0].expression);
  SetExpression *y = (SetExpression*)(children[1].expression);


  
  // if there are k elements of Y that can't be in X, then card(X) + k <= card(Y)
  Vector<Variable> extras;
  int c_offset=0;
  
  while(i<x->elts_ub.size && j<y->elts_ub.size) {

     // std::cout << std::endl << x->elts_ub[i] << " - " <<  y->elts_ub[j] << std::endl;


    // //std::cout << std::endl;
    // for(int nx=0; nx<x->children.size; ++nx) {
    //   std::cout << "X: " << x->children[nx] << " in " << x->children[nx].get_domain() << ", ";
    // }

    // std::cout << std::endl;
    // for(int ny=0; ny<y->children.size; ++ny) {
    //   std::cout  << "Y: " << y->children[ny] << " in " << y->children[ny].get_domain() << ", ";
    // }


    if(x->elts_ub[i] == y->elts_ub[j]) {

       // std::cout << " -> equal, post prec constraint " << x->get_index_var(i) << " in " << x->get_index_var(i).get_domain() << " <= " << y->get_index_var(j) << " in " << y->get_index_var(j).get_domain() << std::endl;

      s->add(x->get_index_var(i) <= y->get_index_var(j));

      //      std::cout << s << std::endl << std::endl;

      ++i;
      ++j;
    } else if(x->elts_ub[i] <= y->elts_ub[j]) {

       // std::cout << " -> " << x->elts_ub[i] << " can't be in Y, hence in X " << x->get_index_var(i) << " in "
       // 		<< x->get_index_var(i).get_domain() << std::endl;

      s->add(x->get_index_var(i) == 0);
      //x->get_index_var(i).set_domain(0);


      //      std::cout << s << std::endl << std::endl;

      ++i;
    } else {

      //    std::cout << " -> " << y->elts_ub[j] << " is not in X, don't care" << std::endl;

      if(y->get_index_var(j).get_min()>0) {
	//std::cout << " (but must be in Y)" << std::endl;
	++c_offset;
      } else {
	//std::cout << " (may be in Y)" << std::endl;
	extras.add(y->get_index_var(j));
      }

      ++j;
    }
  }

  // remaining values of X
  while(i<x->elts_ub.size) {

    // std::cout << x->elts_ub[i] << " - " <<  y->elts_ub[j] << std::endl;
    // std::cout << " -> " << x->elts_ub[i] << " can't be in Y, hence in X* " << x->get_index_var(i) << std::endl;

    s->add(x->get_index_var(i) == 0);

    //std::cout << s << std::endl << std::endl;

    ++i;
  }

  // remaining values of Y
  while(j<y->elts_ub.size) {

    // std::cout << x->elts_ub[i] << " - " <<  y->elts_ub[j] << std::endl;
    // std::cout << " -> " << y->elts_ub[j] << " is not in X, don't care" << std::endl;

    if(y->get_index_var(j).get_min()>0) {
      //std::cout << " (but must be in Y)" << std::endl;
      ++c_offset;
    } else {
      //std::cout << " (may be in Y)" << std::endl;
      extras.add(y->get_index_var(j));
    }
    //std::cout << s << std::endl << std::endl;

    ++j;
  }


  // std::cout << extras << std::endl;

  if(extras.size>0) {
    //extras.add(children[0].get_var());
    s->add(Precedence((BoolSum(extras) + children[0]), c_offset, children[1]));
    //s->add(children[0].get_var() + Sum() <= children[1].get_var());
  } else {
    s->add(Precedence(children[0], c_offset, children[1]));
    //s->add(children[0].get_var() <= children[1].get_var());
  }

  //s->add(children[0].get_var() <= children[1].get_var());
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
  lb = +INFTY;
  ub = -INFTY;
  size = 0;
}

Mistral::MemberExpression::MemberExpression(Variable X, const int lo, const int up) 
  : Expression(X) { 
  lb = lo;
  ub = up;
  size = ub-lb+1;
}

Mistral::MemberExpression::MemberExpression(Variable X, const BitSet& s) 
  : Expression(X) {
  lb = s.min();
  ub = s.max();
  size = s.size();
}

Mistral::MemberExpression::MemberExpression(Variable X, const Vector<int>& s) 
  : Expression(X) {
  lb = s.front();
  ub = s.back();
  size = s.size; 
}

Mistral::MemberExpression::~MemberExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete member expression" << std::endl;
#endif
}

void Mistral::MemberExpression::extract_constraint(Solver *s) { 

  //unsigned int i = 0, j = 0;
  if(children.size == 2) {
    // set variable
    
    VarArray scp;
    SetExpression *y = (SetExpression*)(children[1].expression);
    s->add(children[1] >= 1); // at least one element in the set
    int vali, vnxt = children[0].get_min(), idx;
    do {
      scp.clear();
      vali = vnxt;
      idx = y->get_element_index(vali);
      if(idx>=0) {
	Variable x(0,1);
	scp.add(children[0]);
	scp.add(x);
	s->add(Constraint(new PredicateConstantEqual(scp,vali)));
	scp.clear();
	scp.add(x);
	scp.add(y->children[idx]);
	s->add(Constraint(new ConstraintLess(scp)));
	//s->add((children[0] == vali) <= y->children[idx]);
      } else {
	children[0].remove(vali);

	//s->add(children[0] != vali);
      }
      vnxt = children[0].next(vali);
    } while(vali < vnxt);
  } else if(size == (ub-lb+1)) {
    // interval 
    
    if(children[0].set_min(lb) == FAIL_EVENT) { s->fail(); }
    if(children[0].set_max(ub) == FAIL_EVENT) { s->fail(); }
  } else {
    // set

    if(children[0].set_domain(values) == FAIL_EVENT) { s->fail(); }
  }   

}

void Mistral::MemberExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::MemberExpression::extract_predicate(Solver *s) { 

  //unsigned int i = 0, j = 0;
  if(children.size == 3) {
    // set variable
    
    VarArray scp;
    VarArray conjunction;
    SetExpression *y = (SetExpression*)(children[1].expression);

    Variable empty(0,1);
    scp.add(children[1]);
    scp.add(empty);
    s->add(Constraint(new PredicateLowerBound(scp, 1)));

    scp.clear();
    scp.add(children[2]);
    scp.add(empty);
    s->add(Constraint(new ConstraintLess(scp)));

    //s->add(children[1] >= 1); // at least one element in the set
    int vali, vnxt = children[0].get_min(), idx;
    do {
      scp.clear();
      vali = vnxt;
      idx = y->get_element_index(vali);
      if(idx>=0) {
	Variable x(0,1);
	scp.add(children[0]);
	scp.add(x);
	s->add(Constraint(new PredicateConstantEqual(scp,vali,1)));
	
    	scp.clear();

	Variable b(0,1);
	scp.add(x);
	scp.add(y->children[idx]);
	scp.add(b);
	s->add(Constraint(new PredicateLess(scp)));

	conjunction.add(b);

	//s->add((children[0] == vali) <= y->children[idx]);
      } else {

	Variable b(0,1);
	scp.add(children[0]);
	scp.add(b);
	s->add(Constraint(new PredicateConstantEqual(scp,vali,0)));

	conjunction.add(b);

	//children[0].remove(vali);

	//s->add(children[0] != vali);
      }
      vnxt = children[0].next(vali);
    } while(vali < vnxt);

    Variable N(0, conjunction.size);
    conjunction.add(N);
    s->add(new PredicateBoolSum(conjunction));

    scp.clear();
    scp.add(N);
    scp.add(children[2]);
    s->add(Constraint(new PredicateConstantEqual(scp,conjunction.size-1,1)));

  } else if(size == (ub-lb+1)) {
    // interval 
    
    s->add(Constraint(new PredicateIntervalMember(children,lb,ub)));
    //children[0]->set_min(lb);
    //children[0]->set_max(ub);
  } else {
    // set
    
    s->add(Constraint(new PredicateSetMember(children,values)));
    //children[0]->set_domain(values);
  }
}

const char* Mistral::MemberExpression::get_name() const {
  return "member";
}

Mistral::Variable Mistral::Member(Variable X, Variable Y) {
  Variable exp(new MemberExpression(X, Y));
  return exp;
}

Mistral::Variable Mistral::Member(Variable X, const int lo, const int up) {
  Variable exp(new MemberExpression(X, lo, up));
  return exp;
}

Mistral::Variable Mistral::Member(Variable X, const BitSet& s) {
  Variable exp(new MemberExpression(X, s));
  return exp;
}

Mistral::Variable Mistral::Member(Variable X, const Vector<int>& s) {
  Variable exp(new MemberExpression(X, s));
  return exp;
}


Mistral::Goal::Goal(method t) : type(t) {
  lower_bound = 0;
  upper_bound = 0;
}

Mistral::Goal::Goal(method t, Variable X) : type(t) {
  objective = X;

  // std::cout << "OBJECTIVE=" << objective << " in " << objective.get_domain() << std::endl;

  lower_bound = objective.get_min()-1; //(type == MAXIMIZATION);
  upper_bound = objective.get_max()+1; //(type == MINIMIZATION);
}

Mistral::Goal::~Goal() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete goal" << std::endl;
#endif
}

bool Mistral::Goal::enforce() {
  if(type == MINIMIZATION) {
    return IS_FAIL(objective.set_max(upper_bound-1));
  } else if(type == MAXIMIZATION) {
    return IS_FAIL(objective.set_min(lower_bound+1));
  }
  return false;
}


int Mistral::Goal::value() const {
  return(type == MINIMIZATION ? upper_bound : lower_bound);
  // if(type == MINIMIZATION) {
  //   return upper_bound;
  // } else if(type == MAXIMIZATION) {
  //   return lower_bound;
  // } else return -1;
  // return false;
}

// int Mistral::Goal::best() const {
//   return(type = MINIMIZATION ? upper_bound : lower_bound);
//   // if(type == MINIMIZATION) {
//   //   return upper_bound;
//   // } else if(type == MAXIMIZATION) {
//   //   return lower_bound;
//   // } else return -1;
//   // return false;
// }


bool Mistral::Goal::is_optimization() const {
  return (type == MINIMIZATION || type == MAXIMIZATION);
}

bool Mistral::Goal::improving(const int val) const {
  return( type == SATISFACTION ? false :
	  (type == MINIMIZATION ? 
	   (val < upper_bound) : 
	   (val > lower_bound)) );
}
    

Mistral::Outcome Mistral::Goal::notify_exhausted() {
  if(type == SATISFACTION)
    return UNSAT;
  return OPT;
}

std::ostream& Mistral::Goal::display(std::ostream& os) const {
  if(type == MINIMIZATION) {
    os << "minimize " << objective ;
  } else   if(type == MAXIMIZATION) {
    os << "maximize " << objective ;
  } else if(type == ALLSOLUTIONS) {
    os << "find all solutions" ;
  } else {
    os << "find any solution" ;
  }
  return os;
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




