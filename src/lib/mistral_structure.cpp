
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

#include <mistral_structure.hpp>
#include <mistral_global.hpp>



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




int Mistral::__modulo_fct__(const int x, const int m) {
  int mod = x%m;
  if(mod && (mod<0) != (m<0))  mod += m;
  return mod;
}


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

Mistral::Interval::Interval() {min=+INFTY; max=-INFTY;}
Mistral::Interval::Interval(const BiInterval b) {
  if(b.positive.empty()) max = (b.zero ? 0 : b.negative.max);
  else max = b.positive.max;

  if(b.negative.empty()) min = (b.zero ? 0 : b.positive.min);
  else min = b.negative.min;
}
Mistral::Interval::Interval(const int _min, const int _max) {min=_min; max=_max;}
Mistral::Interval::~Interval() {}

Mistral::Interval Mistral::Interval::get_union(Mistral::Interval dom) {

  //std::cout << "[" << min << "," << max << "] U [" << dom.min << "," << dom.max << "] = [";

  Interval I((min<=dom.min ? min : dom.min), (max>=dom.max ? max : dom.max));

  //std::cout << I.min << "," << I.max << "]\n";

  return I;
}
  
bool Mistral::Interval::contain(const int x) const {
  return (min <= x && x <= max);
}

bool Mistral::Interval::empty() const { return min>max; }

// the interval I such that this%mod = I
Mistral::Interval Mistral::Interval::operator%(const int mod) {
  // the value of the modulo increase with higher values, unless 
  // unless they are distant enough to go back to the next cycle
  
  int I_min = __modulo_fct__(min,mod);
  int I_max = __modulo_fct__(max,mod);

  if(mod>0) { // positive modulo, hence I \in [0, INFTY]
    if(min + mod - I_min <= max) I_min = 0; // can we reach the next 0?
    if(max - I_max - 1 >= min) I_max = mod-1; // can we reach the next mod-1?
  } else { // negtive modulo, hence I \in [-INFTY, 0]
    if(min - I_min + 1 <= max) I_min = 1+mod;
    if(max + mod - I_max >= min) I_max = 0;
  }

  return Interval(I_min, I_max);
}


Mistral::Interval Mistral::Interval::operator*(const Mistral::Interval arg) {
  BiInterval bi(*this);

  // std::cout << "[" << min << "," << max << "] => [" 
  // 	    << bi.negative.min << "," << bi.negative.max << "|"
  // 	    << (bi.zero ? "0|" : ".|")
  // 	    << bi.positive.min << "," << bi.positive.max << "]\n";
    

  BiInterval bj(arg);

  // std::cout << "[" << arg.min << "," << arg.max << "] => [" 
  // 	    << bj.negative.min << "," << bj.negative.max << "|"
  // 	    << (bj.zero ? "0|" : ".|")
  // 	    << bj.positive.min << "," << bj.positive.max << "]\n";



  BiInterval bk = bi*bj;
  Interval I(bk);


  // std::cout << "[" << I.min << "," << I.max << "] <= [" 
  // 	    << bk.negative.min << "," << bk.negative.max << "|"
  // 	    << (bk.zero ? "0|" : ".|")
  // 	    << bk.positive.min << "," << bk.positive.max << "]\n";



  return I;
}

Mistral::Interval Mistral::Interval::anti_mul(const Mistral::Interval arg) {


  BiInterval bi(*this);


  // std::cout << "[" << min << "," << max << "] => [" 
  // 	    << bi.negative.min << "," << bi.negative.max << "|"
  // 	    << (bi.zero ? "0|" : ".|")
  // 	    << bi.positive.min << "," << bi.positive.max << "]\n";
    

  BiInterval bj(arg);


  // std::cout << "[" << arg.min << "," << arg.max << "] => [" 
  // 	    << bj.negative.min << "," << bj.negative.max << "|"
  // 	    << (bj.zero ? "0|" : ".|")
  // 	    << bj.positive.min << "," << bj.positive.max << "]\n";


  BiInterval bk = bi.anti_mul(bj);
  Interval I(bk);

  // std::cout << "[" << I.min << "," << I.max << "] <= [" 
  // 	    << bk.negative.min << "," << bk.negative.max << "|"
  // 	    << (bk.zero ? "0|" : ".|")
  // 	    << bk.positive.min << "," << bk.positive.max << "]\n";


  return I;
}

Mistral::Interval Mistral::Interval::operator/(const Mistral::Interval arg) {
  BiInterval bi(*this);
  BiInterval bj(arg);
  BiInterval bk = bi/bj;
  return Interval(bk);
}

// Mistral::Interval Mistral::Interval::operator%(const int mod) {
//   BiInterval bi(*this);
//   BiInterval bk = bi%mod;
//   return Interval(bk);
// }


Mistral::NegativeHalfDomain::NegativeHalfDomain() : Interval(INFTY, -INFTY) {}  
Mistral::NegativeHalfDomain::NegativeHalfDomain(const int _min, const int _max) : Interval(_min, _max) {}
Mistral::NegativeHalfDomain::~NegativeHalfDomain() {}

// the interval I such that this*arg = I
Mistral::Interval Mistral::NegativeHalfDomain::operator*(Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(min*arg.max, max*arg.min);
}
Mistral::Interval Mistral::NegativeHalfDomain::operator*(Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(max*arg.max, min*arg.min);
}

// the interval I such that this*I = arg
Mistral::Interval Mistral::NegativeHalfDomain::anti_mul(Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval((int)(ceil((double)min/(double)(arg.min))), (int)(floor((double)max/(double)(arg.max))));
}
Mistral::Interval Mistral::NegativeHalfDomain::anti_mul(Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval((int)(ceil((double)max/(double)(arg.min))), (int)(floor((double)min/(double)(arg.max))));
}

// the interval I such that this/arg = I
Mistral::Interval Mistral::NegativeHalfDomain::operator/(Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();

  //std::cout << min << "/" << arg.min << " , " << max << "/" << arg.max << std::endl;

  return Interval(min/arg.min, max/arg.max);
}
Mistral::Interval Mistral::NegativeHalfDomain::operator/(Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(max/arg.min, min/arg.max);
}

// the interval I such that I/arg = this
//Mistral::Interval Mistral::NegativeHalfDomain::get_dividand(Mistral::PositiveHalfDomain arg);
//Mistral::Interval Mistral::NegativeHalfDomain::get_dividand(Mistral::NegativeHalfDomain arg);
///----> this->operator*(arg) <----///


// the interval I such that arg/I = this
//Mistral::Interval Mistral::NegativeHalfDomain::get_divisor(Mistral::PositiveHalfDomain arg);
//Mistral::Interval Mistral::NegativeHalfDomain::get_divisor(Mistral::NegativeHalfDomain arg);
///----> arg->operator/(this) <----///

// // the interval I such that this%mod = I
// Mistral::Interval Mistral::NegativeHalfDomain::operator%(const int mod) {
//   // the value of the modulo increase with higher values, unless 
//   // unless they are distant enough to go back to the next cycle
  
//   int I_min = __modulo_fct__(min,mod);
//   int I_max = __modulo_fct__(max,mod);

//   if(mod>0) { // positive modulo, hence I \in [0, INFTY]
//     if(min + mod - I_min <= max) I_min = 0; // can we reach the next 0?
//     if(max - I_max - 1 >= min) I_max = mod-1; // can we reach the next mod-1?
//   } else { // negtive modulo, hence I \in [-INFTY, 0]
//     if(min - I_min + 1 <= max) I_min = 1+mod;
//     if(max + mod - I_max >= min) I_max = 0
//   }

//   return Interval(I_min, I_max);
// }

// the interval I such that I%mod = this
//Mistral::Interval Mistral::NegativeHalfDomain::anti_modulo(const int mod);


Mistral::PositiveHalfDomain::PositiveHalfDomain() : Interval(INFTY, -INFTY) {}
Mistral::PositiveHalfDomain::PositiveHalfDomain(const int _min, const int _max) : Interval(_min, _max) {}
Mistral::PositiveHalfDomain::~PositiveHalfDomain() {}

Mistral::Interval Mistral::PositiveHalfDomain::operator*(Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(min*arg.min, max*arg.max);
}
Mistral::Interval Mistral::PositiveHalfDomain::operator*(Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(max*arg.min, min*arg.max);
}

Mistral::Interval Mistral::PositiveHalfDomain::anti_mul(Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();

  // std::cout << "[" << (double)min/(double)(arg.max) << "," << (double)max/(double)(arg.min) << "] => ["
  // 	    << (int)(ceil((double)min/(double)(arg.max))) << "," << (int)(floor((double)max/(double)(arg.min))) << "]\n";

  return Interval((int)(ceil((double)min/(double)(arg.max))), (int)(floor((double)max/(double)(arg.min))));
}
Mistral::Interval Mistral::PositiveHalfDomain::anti_mul(Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();

  return Interval((int)(ceil((double)max/(double)(arg.max))), (int)(floor((double)min/(double)(arg.min))));
}

Mistral::Interval Mistral::PositiveHalfDomain::operator/(Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(min/arg.max, max/arg.min);
}

Mistral::Interval Mistral::PositiveHalfDomain::operator/(Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(max/arg.max, min/arg.min);
}

// Mistral::Interval Mistral::PositiveHalfDomain::anti_div(Mistral::PositiveHalfDomain arg);
// Mistral::Interval Mistral::PositiveHalfDomain::anti_div(Mistral::NegativeHalfDomain arg);

// Mistral::Interval Mistral::PositiveHalfDomain::operator%(const int mod) {
//   // the value of the modulo increase with higher values, unless 
//   // unless they are distant enough to go back to the next cycle
  
//   int I_min = __modulo_fct__(min,mod);
//   int I_max = __modulo_fct__(max,mod);

//   if(mod>0) { // positive modulo, hence I \in [0, INFTY]
//     if(min + mod - I_min <= max) I_min = 0; // can we reach the next 0?
//     if(max - I_max - 1 >= min) I_max = mod-1; // can we reach the next mod-1?
//   } else { // negtive modulo, hence I \in [-INFTY, 0]
//     if(min - I_min + 1 <= max) I_min = 1+mod;
//     if(max + mod - I_max >= min) I_max = 0
//   }

//   return Interval(I_min, I_max);
// }

// // // the return value is the interval I such that I%mod = this
// // // 
// // Mistral::Interval Mistral::PositiveHalfDomain::anti_modulo(const int mod) {
// //   int mod_min =  __modulo_fct__(min,mod);

// // }



Mistral::BiInterval::BiInterval(const int n_min, const int n_max, const int p_min, const int p_max, const bool z) {
  positive.min = p_min;
  positive.max = p_max;
  negative.min = n_min;
  negative.max = n_max;
  zero = z;
}
Mistral::BiInterval::BiInterval(const Interval _neg, const Interval _pos, const bool z) {
  positive.min = _pos.min;
  positive.max = _pos.max;
  if(positive.min <  1) positive.min =  1;

  negative.min = _neg.min;
  negative.max = _neg.max;
  if(negative.max > -1) negative.max = -1;

  zero = (z || _pos.contain(0) || _neg.contain(0));
}
Mistral::BiInterval::BiInterval(const Interval I) {
  initialise(I.min, I.max);  
}

Mistral::BiInterval::BiInterval(const int min, const int max) {
  initialise(min, max);  
}

void Mistral::BiInterval::initialise(const int min, const int max) {

  if(min>=0) {
    positive.min = min+(min==0);
    positive.max = max;

    negative.min = +INFTY;
    negative.max = -INFTY;
    
    zero = (min==0);
  } else if(max<=0) {
    positive.min = +INFTY;
    positive.max = -INFTY;

    negative.min = min;
    negative.max = max-(max==0);
    
    zero = (max==0);
  } else { // min<0 & max>0
    positive.min = 1;
    positive.max = max;

    negative.min = min;
    negative.max = -1;

    zero = true;
  }
}
Mistral::BiInterval::BiInterval() {
  positive.min = +INFTY;
  positive.max = -INFTY;
  
  positive.min = +INFTY;
  positive.max = -INFTY;
  
  zero = false;
}
Mistral::BiInterval::~BiInterval() {}
    
Mistral::BiInterval Mistral::BiInterval::operator*(Mistral::BiInterval arg) {
  Interval pospos = positive * arg.positive;

  //  std::cout << "pospos: [" << pospos.min << "," << pospos.max << "]\n";

  Interval negneg = negative * arg.negative;

  //  std::cout << "negneg: [" << negneg.min << "," << negneg.max << "]\n";


  Interval posneg = positive * arg.negative;

  //  std::cout << "posneg: [" << posneg.min << "," << posneg.max << "]\n";

  Interval negpos = negative * arg.positive;

  //  std::cout << "negpos: [" << negpos.min << "," << negpos.max << "]\n";

  BiInterval I(posneg.get_union(negpos), pospos.get_union(negneg), zero || arg.zero);

  return I;
}
// the interval I such that this*I = arg     (I = this/arg)
Mistral::BiInterval Mistral::BiInterval::anti_mul(Mistral::BiInterval arg) {
  BiInterval I;
  if(zero && arg.zero) {
    // if this and arg can be 0, then I is not constrained
    I = BiInterval(-INFTY, INFTY);
  } else if(arg == 0) {
    I = 0;
  } else {
    Interval pospos = positive.anti_mul(arg.positive);

    //std::cout << "pospos: [" << pospos.min << "," << pospos.max << "]\n";

    Interval negneg = negative.anti_mul(arg.negative);

    //std::cout << "negneg: [" << negneg.min << "," << negneg.max << "]\n";


    Interval posneg = positive.anti_mul(arg.negative);

    //std::cout << "posneg: [" << posneg.min << "," << posneg.max << "]\n";

    Interval negpos = negative.anti_mul(arg.positive);

    //std::cout << "negpos: [" << negpos.min << "," << negpos.max << "]\n";


    // if arg cannot be 0, then I cannot be 0
    I = BiInterval(posneg.get_union(negpos), pospos.get_union(negneg), zero);
  }
  return I;
}
// the interval I such that this/arg = I
Mistral::BiInterval Mistral::BiInterval::operator/(Mistral::BiInterval arg) {
  BiInterval I;
  if(!(arg == 0)){
    if(*this == 0) {
      I = 0;
    } else { 
    
      Interval pospos = positive/(arg.positive);
      
      //  std::cout << "pospos: [" << pospos.min << "," << pospos.max << "]\n";
      
      Interval negneg = negative/(arg.negative);
      
      //  std::cout << "negneg: [" << negneg.min << "," << negneg.max << "]\n";
      
      
      Interval posneg = positive/(arg.negative);
      
      //  std::cout << "posneg: [" << posneg.min << "," << posneg.max << "]\n";
      
      Interval negpos = negative/(arg.positive);
      
      //  std::cout << "negpos: [" << negpos.min << "," << negpos.max << "]\n";
      
      
      // if this can be 0, then I can be 0
      I = BiInterval(posneg.get_union(negpos), pospos.get_union(negneg), zero);
    }
  }
  return I;
}

bool Mistral::BiInterval::operator==(const int x) {
  if(x) {
    if(x<0) return x == negative.min && x == negative.max;
    else return x == positive.min && x == positive.max;
  } 
  return zero && positive.empty() && negative.empty();
}

void Mistral::BiInterval::operator=(const int x) {
  if(x) {
    if(x<0) {
      negative.min = x;
      negative.max = x;
      positive.min = +INFTY;
      positive.max = -INFTY;
      zero = false;
    } else {
      positive.min = x;
      positive.max = x;
      negative.min = +INFTY;
      negative.max = -INFTY;
      zero = false;
    }
  } else {
    positive.min = +INFTY;
    positive.max = -INFTY;
    negative.min = +INFTY;
    negative.max = -INFTY;
    zero = true;
  }
}

// Mistral::BiInterval Mistral::BiInterval::operator%(const int mod) {

// }

// Mistral::BiInterval Mistral::BiInterval::anti_modulo(const int mod);

