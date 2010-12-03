
#include <mistral_variable.hpp>
#include <climits>
#include <iomanip>


using namespace std;
using namespace Mistral;


class ptr_type {
public:
  int id;
  unsigned int stack_id;
};

class wrap_type {
public:
  int type;
  ptr_type *truc;

  wrap_type(ptr_type* x) {truc = x; type=randint(12);}
  wrap_type() {truc = NULL; type=0;}

  void set_stack_id(const int idx) {truc->stack_id = idx;}
  unsigned int get_stack_id() const {return truc->stack_id;}
};


std::ostream& operator<< (std::ostream& os, const ptr_type* x) {
  os << x->id << "." << x->stack_id ;
  return os;
}

std::ostream& operator<< (std::ostream& os, const wrap_type& x) {
  os << x.truc->id << "." << x.truc->stack_id ;
  return os;
}


int main(int argc, char *argv[])
{  
  

  /*
  IntegerVar x0(17,20);
  IntegerVar x1(0,10);
  IntegerVar x2(-50,20);
  IntegerVar x3(58,107);
  IntegerVar x4(4,8);
  IntegerVar x5(8,10);
  IntegerVar x6(10,17);


  Solver s;
  s.add(&x0);
  s.add(&x1);
  s.add(&x2);
  s.add(&x3);
  s.add(&x4);
  s.add(&x5);
  s.add(&x6);



  Vector< IntVar > vars;
  
  vars.add( &x0 );
  vars.add( &x1 );
  vars.add( &x2 );
  vars.add( &x3 );
  vars.add( &x4 );
  vars.add( &x5 );
  vars.add( &x6 );

  //  cout << vars << endl;


  Stack< IntVar > seq( vars, true );

  cout << seq << endl;

  seq.remove(&x3);

  cout << seq << endl;

  seq.remove(&x1);

  cout << seq << endl;

  seq.remove(&x5);

  cout << seq << endl;
  
  seq.size = 7;

  cout << seq << endl;


  Vector< IntVar > first;
  first.add(&x0);
  first.add(&x1);
  first.add(&x2);

  //  cout << first << endl;

  Vector< IntVar > second;
  second.add(&x3);
  second.add(&x4);

  //  cout << second << endl;

  Vector< IntVar > third;
  third.add(&x5);
  third.add(&x6);

  //  cout << third << endl;


  Vector< Vector< IntVar > > strategy(3);
  
  strategy.add(first);
  strategy.add(second);
  strategy.add(third);
  
  cout << strategy << endl;
  */


  int N = 100;
  int M = 50000000;

  ptr_type **x = new ptr_type*[N];
  for(int i=0; i<N; ++i)
    {
      x[i] = new ptr_type();
      x[i]->id = i;
      x[i]->stack_id = 0;
    }

  Vector< ptr_type* > vars;
  for(int i=0; i<N; ++i)
    {
      vars.add(x[i]);
    }

  Stack< ptr_type* > seq( vars, true );

  //cout << seq << endl;



  wrap_type *y = new wrap_type[N];
  for(int i=0; i<N; ++i)
    {
      wrap_type w(x[i]);
      y[i] = w;
    }

  Vector< wrap_type > wars;
  for(int i=0; i<N; ++i)
    {
      wars.add(y[i]);
    }

  IdStack< wrap_type > weq( wars, true );





  int ri;
  double tb, ta;




  usrand(12);
  tb = getRunTime();  
  for(int i=0; i<M; ++i) {
    ri = randint(N);
    if(weq.member(y[ri])) {
      //std::cout << "remove " << y[ri] << endl;
      weq.remove(y[ri]);
      //std::cout << weq << endl << endl;
    } else {
      //std::cout << "add " << y[ri] << endl;
      weq.add(y[ri]);
      //std::cout << weq << endl << endl;
    }
  }
  ta = getRunTime();
  std::cout << (ta-tb) << " " << weq << endl << endl;





  usrand(12);
  tb = getRunTime();
  for(int i=0; i<M; ++i) {
    ri = randint(N);
    if(seq.member(x[ri])) {
      //std::cout << "remove " << x[ri] << endl;
      seq.remove(x[ri]);
      //std::cout << seq << endl << endl;
    } else {
      //std::cout << "add " << x[ri] << endl;
      seq.add(x[ri]);
    }
  }
  ta = getRunTime();
  std::cout << (ta-tb) << " " << seq << endl << endl;





}


