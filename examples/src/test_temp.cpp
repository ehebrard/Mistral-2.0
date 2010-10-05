

#include <iostream>
#include <sys/times.h>

using namespace std;

double getRunTime() {
  double df = 0;
  struct tms usage;
  static int clock_ticks = sysconf(_SC_CLK_TCK);
  times(&usage);
  df=((double)usage.tms_utime+(double)usage.tms_stime)/clock_ticks;
  return df;
}




class thing {

public:
  thing() {}
  virtual ~thing() {}
  
  int transform(int x) {
    return -(x+1);
  }

  int do_something() {
    int total = 2;
    for(int i=0; i<1000000007; ++i) {
      total = this->transform(total);
    }
    return total;
  }

};



class virtual_thing {

public:
  virtual_thing() {}
  virtual ~virtual_thing() {}
  
  virtual int transform(int x) = 0;

  int do_something() {
    int total = 2;
    for(int i=0; i<1000000007; ++i) {
      total = this->transform(total);
    }
    return total;
  }

};


class virtual_thing_child : public virtual_thing {

public:
  virtual_thing_child() {}
  virtual ~virtual_thing_child() {}

  int transform(int x) {
    return -(x+1);
  }

};



class case_thing {

public:

  int _t_;

  case_thing(int t) {_t_ = t;}
  virtual ~case_thing() {}
  
  int transform(int x) {
    switch(_t_) {
    case 0: return -(x+1);
    case 1: return -(x+3);
    case 2: return -(x+2);
    case 3: return -(x+1);
    case 4: return -(x+3);
    case 5: return -(x+2);
    case 6: return -(x+1);
    }
    return 0;
  }

  int do_something() {
    int total = 2;
    for(int i=0; i<1000000007; ++i) {
      total = this->transform(total);
    }
    return total;
  }

};


class if_thing {

public:

  int _t_;

  if_thing(int t) {_t_ = t;}
  virtual ~if_thing() {}
  
  int transform(int x) {
    if(_t_==0) return -(x+1);
    else if(_t_==1) return -(x+3);
    else if(_t_==2) return -(x+2);
    else if(_t_==3) return -(x+2);
    else if(_t_==4) return -(x+3);
    else if(_t_==5) return -(x+2);
    else if(_t_==6) return -(x+1);
    return 0;
  }

  int do_something() {
    int total = 2;
    for(int i=0; i<1000000007; ++i) {
      total = this->transform(total);
    }
    return total;
  }

};


class impl_thing {
public:
  int transform(int x) { std::cout << "should display that" << std::endl; return 0; }
};

class ithing0 : public impl_thing {
public:
  int transform(int x) { return -(x+1); } 
};

class ithing1 : public impl_thing {
public:
  int transform(int x) { return -(x+3); } 
};

class ithing2 : public impl_thing {
public:
  int transform(int x) { return -(x+2); } 
};

class ithing3 : public impl_thing {
public:
  int transform(int x) { return -(x+2); } 
};

class ithing4 : public impl_thing {
public:
  int transform(int x) { return -(x+3); } 
};

class ithing5 : public impl_thing {
public:
  int transform(int x) { return -(x+2); } 
};

class ithing6 : public impl_thing {
public:
  int transform(int x) { return -(x+1); } 
};



class cast_thing {

public:

  int _t_;
  impl_thing* impl;

  cast_thing(int t, impl_thing* i) {_t_ = t; impl=i;}
  virtual ~cast_thing() {}
  
  int transform(int x) {
    if(_t_==0) return ((ithing0*)impl)->transform(x);
    else if(_t_==1) return ((ithing1*)impl)->transform(x);
    else if(_t_==2) return ((ithing2*)impl)->transform(x);
    else if(_t_==3) return ((ithing3*)impl)->transform(x);
    else if(_t_==4) return ((ithing4*)impl)->transform(x);
    else if(_t_==5) return ((ithing5*)impl)->transform(x);
    else if(_t_==6) return ((ithing6*)impl)->transform(x);
    return 0;
  }

  int do_something() {
    int total = 2;
    for(int i=0; i<1000000007; ++i) {
      total = transform(total);
    }
    return total;
  }

};






int main(int argc, char *argv[])
{  
  
  thing *truc_normal = new thing();
  case_thing *truc_a_cas0 = new case_thing(0);
  case_thing *truc_a_cas6 = new case_thing(6);
  if_thing *truc_a_if0 = new if_thing(0);
  if_thing *truc_a_if3 = new if_thing(3);
  if_thing *truc_a_if6 = new if_thing(6);
//   template_thing *truc_a_template0 = new template_thing_child<0>();
//   template_thing *truc_a_template6 = new template_thing_child<6>();
  virtual_thing *truc_virtuel = new virtual_thing_child();

  ithing0 *it0 = new ithing0();
  cast_thing truc_cast0(0,it0);

  ithing3 *it3 = new ithing3();
  cast_thing truc_cast3(3,it3);

  ithing6 *it6 = new ithing6();
  cast_thing truc_cast6(6,it6);


  double time[11];

  int result;

  time[0] = getRunTime();
  //result = truc_normal->do_something();
  result = 2;
  for(int i=0; i<1000000007; ++i) {
    result = truc_normal->transform(result);
  }

  time[1] = getRunTime();
  cout << "normal (stat.) \t" << result << " \t" ;
  cout << (time[1] - time[0]) << endl;

  time[1] = getRunTime();
  result = truc_virtuel->do_something();
  time[2] = getRunTime();
  cout << "virtual \t" << result << " \t" ;
  cout << (time[2] - time[1]) << endl;

  time[2] = getRunTime();
  result = truc_a_cas0->do_something();
  time[3] = getRunTime();
  cout << "case (first) \t" << result << " \t" ;
  cout << (time[3] - time[2]) << endl;

  time[3] = getRunTime();
  result = truc_a_cas6->do_something();
  time[4] = getRunTime();
  cout << "case (last) \t" << result << " \t" ;
  cout << (time[4] - time[3]) << endl;

  time[4] = getRunTime();
  result = truc_a_if0->do_something();
  time[5] = getRunTime();
  cout << "if (first) \t" << result << " \t" ;
  cout << (time[5] - time[4]) << endl;

  time[5] = getRunTime();
  result = truc_a_if3->do_something();
  time[6] = getRunTime();
  cout << "if (mid) \t" << result << " \t" ;
  cout << (time[6] - time[5]) << endl;

  time[6] = getRunTime();
  result = truc_a_if6->do_something();
  time[7] = getRunTime();
  cout << "if (last) \t" << result << " \t" ;
  cout << (time[7] - time[6]) << endl;


  time[7] = getRunTime();
  result = truc_cast0.do_something();
  time[8] = getRunTime();
  cout << "cast (first) \t" << result << " \t" ;
  cout << (time[8] - time[7]) << endl;

  time[8] = getRunTime();
  result = truc_cast3.do_something();
  time[9] = getRunTime();
  cout << "cast (mid) \t" << result << " \t" ;
  cout << (time[9] - time[8]) << endl;

  time[9] = getRunTime();
  result = truc_cast6.do_something();
  time[10] = getRunTime();
  cout << "cast (last) \t" << result << " \t" ;
  cout << (time[10] - time[9]) << endl;

//   time[6] = getRunTime();
//   result = truc_a_template0->do_something();
//   time[7] = getRunTime();
//   cout << "temp. (first) \t" << result << " \t" ;
//   cout << (time[7] - time[6]) << endl;

//   time[7] = getRunTime();
//   result = truc_a_template6->do_something();
//   time[8] = getRunTime();
//   cout << "temp. (last) \t" << result << " \t" ;
//   cout << (time[8] - time[7]) << endl;

}


