
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




// int __modulo_fct__(const int x, const int m) {
//   int mod = x%m;
//   if(mod && (mod<0) != (m<0))  mod += m;
//   return mod;
// }


// /*
//   [HERE WE ASSUME 0<k and 0<=a<=b]

// ================================
// |  min([c, d]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if k>b: a                   |
// |                              |
// |  if a<=k<=b: 0               |
// |                              |
// |  so we can assume that k < a |
// |                              |
// |  if (a%k + b - a) >= k: 0    |
// |  else: a%k                   |
// ================================

// [IF k>0 THEN:]

// ================================
// |  min([c, d]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if k>b: a                   |
// |                              |
// |  if a<=k<=b: 0               |
// |                              |
// |  so we can assume that k < a |
// |                              |
// |  if (a%k + b - a) >= k: 0    |
// |  else: a%k                   |
// ================================


// */
// int min_modulo(const int a, const int b, const int k) {
//   int value = a;
//   if(k<=b) {
//     if(k>=a) value = 0;
//     else {
//       int mod = __modulo_fct__(a,k);
//       if((mod + b - a) >= k) value = 0;
//       else value = mod;
//     }
//   }
//   return value;
// }
// /*
// ================================
// |  max([c, d]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if k>b: b                   |
// |                              |
// |  if a<=k<=b: k-1             |
// |                              |
// |  so we can assume that k < a |
// |                              |
// |  if (b%k - b + a) < 0: k-1   |
// |  else: b%k                   |
// ================================
// */
// int max_modulo(const int a, const int b, const int k) {
//   int value = b;
//   if(k<=b) {
//     if(k>=a) value = k-1;
//     else {
//       int mod = __modulo_fct__(b,k);
//       if((mod - b + a) < 0) value = k-1;
//       else value = mod;
//     }
//   }
//   return value;
// }
// /*
// ================================
// |  min([a, b]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if(a%k < c)                 |
// |    a + (c - a%k)             |
// |  if(a%k > d)                 |
// |    a + (d - a%k) + k         |
// |  otherwise: a                |
// ================================
// */
// int min_antimodulo(const int a, const int c, const int d, const int k) {
//   int value = a, mod = __modulo_fct__(a,k);
//   if(mod < c) 
//     value = a + c - mod;
//   else if(mod > d) 
//     value = a + d - mod + k;
//   return value;
// }
// /*
// ================================
// |  max([a, b]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if(b%k < c)                 |
// |    b - (b%k - c) - k         |
// |  if(b%k > d)                 |
// |    b - (b%k - d)             |
// |  otherwise: b                |
// ================================
// */
// int max_antimodulo(const int b, const int c, const int d, const int k) {
//   int value = b, mod = __modulo_fct__(b,k);
//   if(mod < c) 
//     value = b - mod + c - k;
//   else if(mod > d) 
//     value = b - mod + d;
//   return value;
// }


  
// Mistral::Interval::Interval(const int _min, const int _max) {min=_min; max=_max;}
// Mistral::Interval::~Interval() {}

// Mistral::Interval Mistral::Interval::get_union(Mistral::Interval dom) {
//   Interval I((min<=dom.min_ ? min : dom.min), (max>=dom.max ? max : dom.max));
//   return I;
// }
  
// bool Mistral::Interval::contain(const int x) {
//   return (min <= x && x <= max);
// }

// bool Mistral::Interval::empty() { return min>max; }



// Interval operator*(interval dom) {
  
// }


//   Interval operator*(interval dom);

//   Interval anti_mul(Interval dom);
//   Interval anti_mul(Interval dom);

//   Interval operator/(Interval dom);
//   Interval operator/(Interval dom);

//   Interval operator%(const int mod);
//   Interval operator%(const int mod);

  
// Mistral::NegativeHalfDomain::NegativeHalfDomain(const int _min, const int _max) : Interval(_min, _max) {}
// Mistral::NegativeHalfDomain::~NegativeHalfDomain() {}

// // the interval I such that this*arg = I
// Mistral::Interval Mistral::NegativeHalfDomain::operator*(Mistral::PositiveHalfDomain arg) {
//   return Interval(max*arg.max, min*arg.min);
// }
// Mistral::Interval Mistral::NegativeHalfDomain::operator*(Mistral::NegativeHalfDomain arg) {
//   return Interval(min*arg.min, max*arg.max);
// }

// // the interval I such that this*I = arg
// Mistral::Interval Mistral::NegativeHalfDomain::anti_mul(Mistral::PositiveHalfDomain arg) {
//   return Interval((int)(ceil((double)max/(double)(arg.min))), (int)(floor((double)min/(double)(arg.max))));
// }
// Mistral::Interval Mistral::NegativeHalfDomain::anti_mul(Mistral::NegativeHalfDomain arg) {
//   return Interval((int)(ceil((double)min/(double)(arg.max))), (int)(floor((double)max/(double)(arg.min))));
// }

// // the interval I such that this/arg = I
// Mistral::Interval Mistral::NegativeHalfDomain::operator/(Mistral::PositiveHalfDomain arg) {
//   return Interval(max/arg.min, min/arg.max);
// }
// Mistral::Interval Mistral::NegativeHalfDomain::operator/(Mistral::NegativeHalfDomain arg) {
//   return Interval(min/arg.max, max/arg.min);
// }

// // the interval I such that I/arg = this
// //Mistral::Interval Mistral::NegativeHalfDomain::get_dividand(Mistral::PositiveHalfDomain arg);
// //Mistral::Interval Mistral::NegativeHalfDomain::get_dividand(Mistral::NegativeHalfDomain arg);
// ///----> this->operator*(arg) <----///


// // the interval I such that arg/I = this
// //Mistral::Interval Mistral::NegativeHalfDomain::get_divisor(Mistral::PositiveHalfDomain arg);
// //Mistral::Interval Mistral::NegativeHalfDomain::get_divisor(Mistral::NegativeHalfDomain arg);
// ///----> arg->operator/(this) <----///

// // the interval I such that this%mod = I
// Mistral::Interval Mistral::NegativeHalfDomain::operator%(const int mod);
// Mistral::Interval Mistral::NegativeHalfDomain::operator%(const int mod);

// // the interval I such that I%mod = this
// Mistral::Interval Mistral::NegativeHalfDomain::anti_modulo(const int mod);
// Mistral::Interval Mistral::NegativeHalfDomain::anti_modulo(const int mod);



// Mistral::PositiveHalfDomain::PositiveHalfDomain(const int _min, const int _max) : Interval(_min, _max) {}
// Mistral::PositiveHalfDomain::~PositiveHalfDomain() {}

// Mistral::Interval Mistral::PositiveHalfDomain::operator*(Mistral::PositiveHalfDomain dom);
// Mistral::Interval Mistral::PositiveHalfDomain::operator*(Mistral::NegativeHalfDomain dom);

// Mistral::Interval Mistral::PositiveHalfDomain::anti_mul(Mistral::PositiveHalfDomain dom);
// Mistral::Interval Mistral::PositiveHalfDomain::anti_mul(Mistral::NegativeHalfDomain dom);

// Mistral::Interval Mistral::PositiveHalfDomain::operator/(Mistral::PositiveHalfDomain dom);
// Mistral::Interval Mistral::PositiveHalfDomain::operator/(Mistral::NegativeHalfDomain dom);

// Mistral::Interval Mistral::PositiveHalfDomain::anti_div(Mistral::PositiveHalfDomain dom);
// Mistral::Interval Mistral::PositiveHalfDomain::anti_div(Mistral::NegativeHalfDomain dom);

// Mistral::Interval Mistral::PositiveHalfDomain::operator%(const int mod) {
//   int I_min;
  
//   if(mod>0) {
//     I_min = min; // this is the case where mod>max, then there is an equality
//     if(mod<=max) { // otherwise
//       if(mod>=min) I_min = 0; // if mod >= min and mod mod <= max, then mod%mod = 0
//       else {
// 	int remainder = __modulo_fct__(min,mod);
// 	if((remainder + max - min) >= mod) I_min = 0;
// 	else I_min = remainder;
//       }
//     }
//   } else {
//     I_min = mod+1;
//     if(-mod<=min) {

//     }
//   }
//   return value;
// }

// Mistral::Interval Mistral::PositiveHalfDomain::operator%(const int mod);

//   // the return value is the interval I such that I%mod = this
// Mistral::Interval Mistral::PositiveHalfDomain::anti_modulo(const int mod);
// Mistral::Interval Mistral::PositiveHalfDomain::anti_modulo(const int mod);



// Mistral::BiInterval::BiInterval(const int n_min, const int n_max, const int p_min, const int p_max, const bool z) {
//   positive.min = p_min;
//   positive.max = p_max;
//   negative.min = n_min;
//   negative.max = n_max;
//   zero = z;
// }
// Mistral::BiInterval::BiInterval(Interval _neg, Interval _pos, const bool z) {
//   positive.min = _pos.min;
//   positive.max = _pos.max;
//   if(positive.min <  1) positive.min =  1;

//   negative.min = _neg.min;
//   negative.max = _neg.max;
//   if(negative.max > -1) negative.max = -1;

//   zero = (z || _pos.contain(0) || _neg.contain(0));
// }
// Mistral::BiInterval::BiInterval(const int min, const int max) {
//   if(min>=0) {
//     positive.min = min+(min==0);
//     positive.max = max;

//     negative.min = +INFTY;
//     negative.max = -INFTY;
    
//     zero = (min==0);
//   } else if(max<=0) {
//     positive.min = +INFTY
//     positive.max = -INFTY;

//     negative.min = min;
//     negative.max = max-(max==0);
    
//     zero = (max==0);
//   } else { // min<0 & max>0
//     positive.min = 1;
//     positive.max = max;

//     negative.min = min;
//     negative.max = -1;

//     zero = true;
//   }
// }
// Mistral::BiInterval::~BiInterval() {}
    
// Mistral::BiInterval Mistral::BiInterval::operator*(Mistral::BiInterval arg) {
//   Interval pospos = positive * arg.positive;
//   Interval negneg = negative * arg.negative;

//   Interval posneg = positive * arg.negative;
//   Interval negpos = negative * arg.positive;

//   BiInterval I(posneg.get_union(negpos), pospos.get_union(negneg), zero || arg.zero);

//   return I;
// }
// // the interval I such that this*I = arg     (I = this/arg)
// Mistral::BiInterval Mistral::BiInterval::anti_mul(Mistral::BiInterval arg) {
//   BiInterval I();
//   if(zero && arg.zero) {
//     // if this and arg can be 0, then I is not constrained
//     I = BiInterval(-INFTY, INFTY);
//   } else if(arg == 0) {
//     I = 0;
//   } else {
//     Interval pospos = positive.anti_mul(arg.positive);
//     Interval negneg = negative.anti_mul(arg.negative);

//     Interval posneg = positive.anti_mul(arg.negative);
//     Interval negpos = negative.anti_mul(arg.positive);

//     // if arg cannot be 0, then I cannot be 0
//     I = BiInterval(posneg.get_union(negpos), pospos.get_union(negneg), arg.zero);
//   }
//   return I;
// }
// // the interval I such that this/arg = I
// Mistral::BiInterval Mistral::BiInterval::operator/(Mistral::BiInterval arg) {
//   BiInterval I();
//   if(*this == 0) {
//     I = 0;
//   } else {
//     Interval pospos = positive/(arg.positive);
//     Interval negneg = negative/(arg.negative);
//     //Interval neg = posneg.get_union(negpos);

//     Interval posneg = positive/(arg.negative);
//     Interval negpos = negative/(arg.positive);
//     //Interval pos = pospos.get_union(negneg);
    
//     // if this can be 0, then I can be 0
//     I = BiInterval(posneg.get_union(negpos), pospos.get_union(negneg), zero);
//   }
//   return I;
// }
// Mistral::BiInterval Mistral::BiInterval::operator%(const int mod) {

// }

// Mistral::BiInterval Mistral::BiInterval::anti_modulo(const int mod);

