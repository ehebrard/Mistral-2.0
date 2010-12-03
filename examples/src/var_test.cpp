
#include <mistral_variable.hpp>
#include <climits>
#include <iomanip>


using namespace std;
using namespace Mistral;




int main(int argc, char *argv[])
{  
  
  IntegerVar x0(0,10);
  IntegerVar x1(-50,20);
  IntegerVar x2(58,107);

  cout << "start (level 0)" << endl;

  Solver s;

  s.add(&x0);
  s.add(&x1);
  s.add(&x2);

  s.initialise();


  cout << x0 << " in " << x0.domain << endl
       << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl << endl;

  x0._assert_();
  x1._assert_();
  x2._assert_();

  cout << "up one level (level 1)" << endl;

  s.make_node();

  cout << "change x0, x1 and x2" << endl;


  x0.remove(5);
  x1.setMax(50);
  x2.removeRange(-10, 75);
    

  cout << x0 << " in " << x0.domain << endl
       << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl << endl;

  x0._assert_();
  x1._assert_();
  x2._assert_();

  cout << "change x1 and x2" << endl;

  x1.remove(10);
  x2.removeRange(91, 101);
    

  cout << x0 << " in " << x0.domain << endl
       << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl << endl;

  x0._assert_();
  x1._assert_();
  x2._assert_();

  cout << "up one level (level 2)" << endl;

  s.make_node();


  cout << "change x0 and x2" << endl;

  BitSet e(76,107,BitSet::empt);
  e.insert(77);
  e.insert(87);
  e.insert(97);
  e.insert(107);  


  cout << e << endl;

  x0.setMin(3);
  x2.removeSet(e);
  x2.setMin(80);


  cout << x0 << " in " << x0.domain << endl
       << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl << endl;

  x0._assert_();
  x1._assert_();
  x2._assert_();

  cout << "up one level (level 3)" << endl;

  s.make_node();
    

  cout << "change x0 and x1" << endl;

  x1.setDomain(11);
  x0.setMax(2);
    

  cout << x0 << " in " << x0.domain << endl
       << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl << endl;

  x0._assert_();
  x1._assert_();
  x2._assert_();

  cout << "backtrack to level 2" << endl;

  s.backtrack();
    

  cout << x0 << " in " << x0.domain << endl
       << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl << endl;

  x0._assert_();
  x1._assert_();
  x2._assert_();


  cout << "backtrack to level 1" << endl;

  s.backtrack();
    

  cout << x0 << " in " << x0.domain << endl
       << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl << endl;

  x0._assert_();
  x1._assert_();
  x2._assert_();


  cout << "backtrack to level 0" << endl;

  s.backtrack();

  cout << x0 << " in " << x0.domain << endl
       << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl << endl;

  x0._assert_();
  x1._assert_();
  x2._assert_();








  Variable y0(0,10);
  Variable y1(-50,20);
  Variable y2(58,107);


//   cout << "start (level 0)" << endl;

//   Solver s;

//   s.add(&x0);
//   s.add(&x1);
//   s.add(&x2);

//   s.initialise();


//   cout << x0 << " in " << x0.domain << endl
//        << x1 << " in " << x1.domain << endl
//        << x2 << " in " << x2.domain << endl << endl;

//   x0._assert_();
//   x1._assert_();
//   x2._assert_();

//   cout << "up one level (level 1)" << endl;

//   s.make_node();

//   cout << "change x0, x1 and x2" << endl;


//   x0.remove(5);
//   x1.setMax(50);
//   x2.removeRange(-10, 75);
    

//   cout << x0 << " in " << x0.domain << endl
//        << x1 << " in " << x1.domain << endl
//        << x2 << " in " << x2.domain << endl << endl;

//   x0._assert_();
//   x1._assert_();
//   x2._assert_();

//   cout << "change x1 and x2" << endl;

//   x1.remove(10);
//   x2.removeRange(91, 101);
    

//   cout << x0 << " in " << x0.domain << endl
//        << x1 << " in " << x1.domain << endl
//        << x2 << " in " << x2.domain << endl << endl;

//   x0._assert_();
//   x1._assert_();
//   x2._assert_();

//   cout << "up one level (level 2)" << endl;

//   s.make_node();


//   cout << "change x0 and x2" << endl;

//   BitSet e(76,107,BitSet::empt);
//   e.insert(77);
//   e.insert(87);
//   e.insert(97);
//   e.insert(107);  


//   cout << e << endl;

//   x0.setMin(3);
//   x2.removeSet(e);
//   x2.setMin(80);


//   cout << x0 << " in " << x0.domain << endl
//        << x1 << " in " << x1.domain << endl
//        << x2 << " in " << x2.domain << endl << endl;

//   x0._assert_();
//   x1._assert_();
//   x2._assert_();

//   cout << "up one level (level 3)" << endl;

//   s.make_node();
    

//   cout << "change x0 and x1" << endl;

//   x1.setDomain(11);
//   x0.setMax(2);
    

//   cout << x0 << " in " << x0.domain << endl
//        << x1 << " in " << x1.domain << endl
//        << x2 << " in " << x2.domain << endl << endl;

//   x0._assert_();
//   x1._assert_();
//   x2._assert_();

//   cout << "backtrack to level 2" << endl;

//   s.backtrack();
    

//   cout << x0 << " in " << x0.domain << endl
//        << x1 << " in " << x1.domain << endl
//        << x2 << " in " << x2.domain << endl << endl;

//   x0._assert_();
//   x1._assert_();
//   x2._assert_();


//   cout << "backtrack to level 1" << endl;

//   s.backtrack();
    

//   cout << x0 << " in " << x0.domain << endl
//        << x1 << " in " << x1.domain << endl
//        << x2 << " in " << x2.domain << endl << endl;

//   x0._assert_();
//   x1._assert_();
//   x2._assert_();


//   cout << "backtrack to level 0" << endl;

//   s.backtrack();

//   cout << x0 << " in " << x0.domain << endl
//        << x1 << " in " << x1.domain << endl
//        << x2 << " in " << x2.domain << endl << endl;

//   x0._assert_();
//   x1._assert_();
//   x2._assert_();


}


