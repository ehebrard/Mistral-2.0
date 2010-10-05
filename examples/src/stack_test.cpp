
#include <mistral_variable.hpp>
#include <climits>
#include <iomanip>


using namespace std;
using namespace Mistral;




int main(int argc, char *argv[])
{  
  
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

  seq.erase(&x3);

  cout << seq << endl;

  seq.erase(&x1);

  cout << seq << endl;

  seq.erase(&x5);

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
}


