
#include <mistral_structure.hpp>
#include <climits>
#include <iomanip>


using namespace std;
using namespace Mistral;



// class C {

// public:

//   void foo() { cout << "hi, I'm C" << endl; }

// };

// class A : public C {

// public:
//   A() {}
//   void foo() { cout << "hello, I'am A" << endl; }

// };


// class B : public C {

// public:
//   B() {}
//   void foo() { cout << "hello, I'am B" << endl; }

// };


// template < class T >
// class executor {

// public:

//   executor(T* y) {
//     x = y;
//   }

//   T* x;
//   void execute_foo() {
//     (*x).foo();
//   }

// };







int main(int argc, char *argv[])
{  

//   executor *exec_a;
//   executor *exec_b;
//   executor *exec_c;

//   A* a = new A(exec_a);
//   B* b = new B(exec_b);
//   C* c = new A(exec_c);


//   exec_a->execute_foo();
//   exec_b->execute_foo();
//   exec_c->execute_foo();

//   //execute_foo(a);

//   //execute_foo(b);

//   //execute_foo(*x);

//   exit(1);

  iBitset64 *is1 = iBitset64::make_new(-50, 50, iBitset64::full);
  is1->printBits( cout );
  cout << endl;
  //is1->print(cout);
  cout << (*is1) << endl
       << (is1->size()) 
       << " " << (is1->min())
       << " " << (is1->max())
       << endl;

  cout << "\n remove interval [-24,14], [26,49]" << endl;
  is1->removeInterval(-24,14);
  is1->removeInterval(26,49);

  is1->printBits( cout );
  cout << endl;
  //is1->print( cout );
  cout << (*is1) << endl
       << (is1->size()) 
       << " " << (is1->min())
       << " " << (is1->max())
       << endl;

  cout << "\n insert -1, 1, 5" << endl;
  is1->insert(-1);
  is1->insert(1);
  is1->insert(5);

  is1->printBits( cout );
  cout << endl;
  //is1->print( cout );
  cout << (*is1) << endl
       << (is1->size()) 
       << " " << (is1->min())
       << " " << (is1->max())
       << endl;

  cout << "\n insert interval [-10,-5]" << endl;
  is1->addInterval(-10,-5);

  is1->printBits( cout );
  cout << endl;
  //is1->print( cout );
  cout << (*is1) << endl
       << (is1->size()) 
       << " " << (is1->min())
       << " " << (is1->max())
       << endl;

  cout << "\n set min -35" << endl;
  is1->setMin(-35);

  is1->printBits( cout );
  cout << endl;
  //is1->print( cout );
  cout << (*is1) << endl
       << (is1->size()) 
       << " " << (is1->min())
       << " " << (is1->max())
       << endl;


  cout << "\n set max 40" << endl;
  is1->setMax(40);

  is1->printBits( cout );
  cout << endl;
  //is1->print( cout );
  cout << (*is1) << endl
       << (is1->size()) 
       << " " << (is1->min())
       << " " << (is1->max())
       << endl;

  cout << "\n set max 20" << endl;
  is1->setMax(20);

  is1->printBits( cout );
  cout << endl;
  //is1->print( cout );
  cout << (*is1) << endl
       << (is1->size()) 
       << " " << (is1->min())
       << " " << (is1->max())
       << endl;


  cout << "\n increment by 10" << endl;
  is1->increment(10);

  is1->printBits( cout );
  cout << endl;
  //is1->print( cout );
  cout << (*is1) << endl
       << (is1->size()) 
       << " " << (is1->min())
       << " " << (is1->max())
       << endl;

  cout << "\n decrement by 5" << endl;
  is1->decrement(5);

  is1->printBits( cout );
  cout << endl;
  //is1->print( cout );
  cout << (*is1) << endl
       << (is1->size()) 
       << " " << (is1->min())
       << " " << (is1->max())
       << endl;

  iBitset64 *is2 = iBitset64::make_new(-50, 50, iBitset64::empt); 

  cout << "\n negate" << endl;
  is1->negate(*is2);
  is2->printBits( cout );
  cout << endl;
  //is2->print( cout );
  cout << (*is2) << endl
       << (is2->size()) 
       << " " << (is2->min())
       << " " << (is2->max())
       << endl;



  BitSet s1(-50, 50, BitSet::full);

  s1.printBits( cout );
  cout << endl;
  //s1.print( cout );
  cout << s1 << endl
       << (s1.size()) 
       << " " << (s1.min())
       << " " << (s1.max())
       << endl;

  cout << "\n remove interval [-24,14], [26,49]" << endl;
  s1.removeInterval(-24,14);
  s1.removeInterval(26,49);

  s1.printBits( cout );
  cout << endl;
  //s1.print( cout );
  cout << s1 << endl
       << (s1.size()) 
       << " " << (s1.min())
       << " " << (s1.max())
       << endl;

  cout << "\n insert -1, 1, 5" << endl;
  s1.insert(-1);
  s1.insert(1);
  s1.insert(5);

  s1.printBits( cout );
  cout << endl;
  //s1.print( cout );
  cout << s1 << endl
       << (s1.size()) 
       << " " << (s1.min())
       << " " << (s1.max())
       << endl;

  cout << "\n insert interval [-10,-5]" << endl;
  s1.addInterval(-10,-5);

  s1.printBits( cout );
  cout << endl;
  //s1.print( cout );
  cout << s1 << endl
       << (s1.size()) 
       << " " << (s1.min())
       << " " << (s1.max())
       << endl;

  cout << "\n set min -35" << endl;
  s1.setMin(-35);

  s1.printBits( cout );
  cout << endl;
  //s1.print( cout );
  cout << s1 << endl
       << (s1.size()) 
       << " " << (s1.min())
       << " " << (s1.max())
       << endl;


  cout << "\n set max 40" << endl;
  s1.setMax(40);

  s1.printBits( cout );
  cout << endl;
  //s1.print( cout );
  cout << s1 << endl
       << (s1.size()) 
       << " " << (s1.min())
       << " " << (s1.max())
       << endl;

  cout << "\n set max 20" << endl;
  s1.setMax(20);

  s1.printBits( cout );
  cout << endl;
  //s1.print( cout );
  cout << s1 << endl
       << (s1.size()) 
       << " " << (s1.min())
       << " " << (s1.max())
       << endl;


  cout << "\n increment by 10" << endl;
  s1.increment(10);

  s1.printBits( cout );
  cout << endl;
  //s1.print( cout );
  cout << s1 << endl
       << (s1.size()) 
       << " " << (s1.min())
       << " " << (s1.max())
       << endl;

  cout << "\n decrement by 5" << endl;
  s1.decrement(5);

  s1.printBits( cout );
  cout << endl;
  //s1.print( cout );
  cout << s1 << endl
       << (s1.size()) 
       << " " << (s1.min())
       << " " << (s1.max())
       << endl;


  BitSet s2(-50, 50, BitSet::empt); 

  cout << "\n negate" << endl;
  s1.negate(s2);
  s2.printBits( cout );
  cout << endl;
  //s2.print( cout );
  cout << s2 << endl
       << (s2.size()) 
       << " " << (s2.min())
       << " " << (s2.max())
       << endl;

}


