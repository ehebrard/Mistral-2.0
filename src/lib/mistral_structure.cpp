
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


#include <mistral_structure.hpp>



std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::IntStack& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Queue& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::MultiSet& x) {
  return x.display(os);
}

// std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::VariableQueue& x) {
//   return x.display(os);
// }

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::IntStack* x) {
  return (x ? x->display(os) : os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Queue* x) {
  return (x ? x->display(os) : os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::MultiSet* x) {
  return (x ? x->display(os) : os);
}

// std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::VariableQueue* x) {
//   return x->display(os);
// }
