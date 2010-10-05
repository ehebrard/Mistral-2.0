
#include <mistral_variable.hpp>
#include <climits>
#include <iomanip>


using namespace std;
using namespace Mistral;





int main(int argc, char *argv[])
{  
  
  IntegerVar x0(0,10);
  IntegerVar x1(0,10);
  IntegerVar x2(0,1);
  IntegerVar x3(0,1);
  IntegerVar x4(4,8);
  IntegerVar x5(6,10);
  IntegerVar x6(8,12);


  Solver s;

  s.add(&x0);
  s.add(&x1);
  s.add(&x2);
  s.add(&x3);
  s.add(&x4);
  s.add(&x5);
  s.add(&x6);

  
  Vector< IntVar > scope(3);

  scope.add(&x0);
  scope.add(&x1);
  s.add(new ConstraintNotEqual(scope));

  scope.clear();
  scope.add(&x4);
  scope.add(&x5);
  scope.add(&x2);
  s.add(new PredicateEqual(scope));

  scope.clear();
  scope.add(&x0);
  scope.add(&x1);
  scope.add(&x2);
  s.add(new PredicateEqual(scope));

  scope.clear();
  scope.add(&x5);
  scope.add(&x6);
  scope.add(&x3);
  s.add(new PredicateEqual(scope));


  scope.clear();
  scope.add(&x5);
  scope.add(&x6);
  s.add(new ConstraintLess(scope));



  Vector< IntVar > first;
  first.add(&x0);
  first.add(&x1);
  first.add(&x2);
  first.add(&x4);

  Vector< IntVar > second;
  second.add(&x6);
  second.add(&x5);
  second.add(&x3);


  cout << s << endl;

  s.initialise();
  
  cout << "search on " << first << endl;
  s.depth_first_search(first);

  cout << s << endl;
  
  cout << "search on " << second << endl;
  s.depth_first_search(second);

  cout << s << endl;

}


