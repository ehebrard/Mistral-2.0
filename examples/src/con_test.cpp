
#include <mistral_constraint.hpp>
#include <mistral_variable.hpp>
#include <climits>
#include <iomanip>


using namespace std;
using namespace Mistral;




int main(int argc, char *argv[])
{  
  
  IntegerVar x0(0,10);
  IntegerVar x1(0,9);
  IntegerVar x2(0,8);
  IntegerVar x3(0,1);

//   s.add(&x0);
//   s.add(&x1);
//   s.add(&x2);
//   s.add(&x3);


  Vector< IntVar > scp(3);
  scp.add(&x0);
  scp.add(&x1);

  ConstraintNotEqual c1(scp);

  scp.clear();
  scp.add(&x1);
  scp.add(&x2);

  ConstraintNotEqual c2(scp);

  scp.clear();
  scp.add(&x0);
  scp.add(&x2);

  ConstraintLess c3(scp, 2);

  scp.clear();
  scp.add(&x0);
  scp.add(&x1);
  scp.add(&x3);

  PredicateEqual c4(scp);

  
  
//   cout << x0 << endl
//        << x1 << endl
//        << x2 << endl
//        << x3 << endl
//        << &c1 << endl
//        << &c2 << endl
//        << &c3 << endl
//        << &c4 << endl;

  Solver s;
  
  s.add(&c1);
  s.add(&c2);
  s.add(&c3);
  s.add(&c4);

  s.initialise();

  cout << s.variables << endl;

  
  cout << "INIT:" << endl;

  s.debug_print();

  cout << "\nSET MAX(x_0) TO 5" << endl;

  x0.setMax(5);

  s.debug_print();

  cout << "\nSET MIN(x_0) TO 5" << endl;

  x0.setMin(5);

  s.debug_print();

//   cout << "\nREMOVE 5 FROM (x_1)" << endl;

//   x1.remove(5);

//   s.debug_print();

  cout << "\nPROPAGATE" << endl;

  s.propagate();

  s.debug_print();

  
  cout << "\nBACKTRACK" << endl;

  s.backtrack();

  s.debug_print();

  //ACQueue< 3 > bidule;



}


