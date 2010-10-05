
#include <mistral_variable.hpp>
#include <climits>
#include <iomanip>


using namespace std;
using namespace Mistral;




int main(int argc, char *argv[])
{  
  
  IntegerVar x1(0,10);
  IntegerVar x2(-50,20);
  IntegerVar x3(58,107);

  cout << "start (level 0)" << endl;

  Solver s;

  s.add(&x1);
  s.add(&x2);
  s.add(&x3);

  s.initialise();

  //x1.debug_print();
  //x2.debug_print();
  //x3.debug_print();

  cout << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl
       << x3 << " in " << x3.domain << endl << endl;

  x1._assert_();
  x2._assert_();
  x3._assert_();

  cout << "change x0, x1 and x2" << endl;


  x1.remove(5);
  x2.setMax(50);
  x3.removeRange(-10, 75);
    
  //x1.debug_print();
  //x2.debug_print();
  //x3.debug_print();

  cout << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl
       << x3 << " in " << x3.domain << endl << endl;

  x1._assert_();
  x2._assert_();
  x3._assert_();

  cout << "change x1 and x2" << endl;

  x2.remove(10);
  x3.removeRange(91, 101);
    
  //x1.debug_print();
  //x2.debug_print();
  //x3.debug_print();

  cout << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl
       << x3 << " in " << x3.domain << endl << endl;

  x1._assert_();
  x2._assert_();
  x3._assert_();

  cout << "up one level (level 1)" << endl;

  s.make_node();

  //s.debug_print();

  //cout << (s.trail_size) << endl;

  //x1.debug_print();
  //x2.debug_print();
  //x3.debug_print();

//   cout << x1 << endl
//        << x2 << endl
//        << x3 << endl << endl;

//   x1._assert_();
//   x2._assert_();
//   x3._assert_();

  cout << "change x0 and x2" << endl;

  BitSet e(76,107,BitSet::empt);
  e.insert(77);
  e.insert(87);
  e.insert(97);
  e.insert(107);  

  x1.setMin(3);
  x3.removeSet(e);
  x3.setMin(80);

  //x1.debug_print();
  //x2.debug_print();
  //x3.debug_print();

  cout << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl
       << x3 << " in " << x3.domain << endl << endl;

  x1._assert_();
  x2._assert_();
  x3._assert_();

  cout << "up one level (level 2)" << endl;

  s.make_node();
    
  //s.debug_print();

  //cout << (s.trail_size) << endl;

  //x1.debug_print();
  //x2.debug_print();
  //x3.debug_print();

//   cout << x1 << endl
//        << x2 << endl
//        << x3 << endl << endl;

//   x1._assert_();
//   x2._assert_();
//   x3._assert_();

  cout << "change x0 and x1" << endl;

  //x2.setMin(11);
  //x2.setMax(12);
  x2.setDomain(11);
  x1.setMax(2);
    
  //x1.debug_print();
  //x2.debug_print();
  //x3.debug_print();

  cout << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl
       << x3 << " in " << x3.domain << endl << endl;

  x1._assert_();
  x2._assert_();
  x3._assert_();

  cout << "backtrack to level 1" << endl;

  s.backtrack();
    
  //x1.debug_print();
  //x2.debug_print();
  //x3.debug_print();

  cout << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl
       << x3 << " in " << x3.domain << endl << endl;

  x1._assert_();
  x2._assert_();
  x3._assert_();


  cout << "backtrack to level 0" << endl;

  s.backtrack();
    
  //x1.debug_print();
  //x2.debug_print();
  //x3.debug_print();

  cout << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl
       << x3 << " in " << x3.domain << endl << endl;

  x1._assert_();
  x2._assert_();
  x3._assert_();

  // cout << "backtrack!!" << endl;

//   s.backtrack();
    
//   //x1.debug_print();
//   //x2.debug_print();
//   //x3.debug_print();

//   cout << x1 << endl
//        << x2 << endl
//        << x3 << endl << endl;

//   x1._assert_();
//   x2._assert_();
//   x3._assert_();


  cout << "backtrack to level -1" << endl;

  s.backtrack();
    
  //x1.debug_print();
  //x2.debug_print();
  //x3.debug_print();

  cout << x1 << " in " << x1.domain << endl
       << x2 << " in " << x2.domain << endl
       << x3 << " in " << x3.domain << endl << endl;

  x1._assert_();
  x2._assert_();
  x3._assert_();

  // cout << "backtrack!!" << endl;

//   s.backtrack();
    
//   //x1.debug_print();
//   //x2.debug_print();
//   //x3.debug_print();

//   cout << x1 << endl
//        << x2 << endl
//        << x3 << endl << endl;

//   x1._assert_();
//   x2._assert_();
//   x3._assert_();

  
  

}


