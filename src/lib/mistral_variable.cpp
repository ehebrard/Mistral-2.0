
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

void Mistral::FiniteDomain::initialise(const int lb, const int ub, const bool vals) {
  min = lb;
  max = ub;
  size = ub-lb+1;
  if(vals) values.initialise(lb, ub, BitSet::full);
}

// std::string Mistral::toString(const Mistral::VariableBitset& x) {
//   return x.getString();
// }

// std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::VariableBitset& x) {
//   return os << x.getString();
// }


// std::string Mistral::toString(const Mistral::FiniteDomain& x) {
//   return x.getString();
// }

//void Mistral::Variable_assert_constraints_();


std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::FiniteDomain& x) {
  //return os << x.getString();
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::FiniteDomain* x) {
  //return os << x.getString();
  return x->display(os);
}

void Mistral::FiniteDomain::_assert_() {
  if( values.table && (values.min() != min) ) {
    std::cerr << (*this) << " Error: inconsistent domain struct (min)" << std::endl;
    exit(0);
  }
  if( values.table && (values.max() != max) ) {
    std::cerr << (*this) << " Error: inconsistent domain struct (max)" << std::endl;
    exit(0);
  }
  if( (values.table || (size != max-min+1)) && ((unsigned int)size != values.size()) ) {
    std::cerr << (*this) << " Error: inconsistent domain struct (size)" << std::endl;
    exit(0);
  }
}

//std::ostream& operator<< (std::ostream& os, const VariableInt* x);


