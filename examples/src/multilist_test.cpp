
#include <mistral_backtrack.hpp>
#include <climits>
#include <iomanip>


using namespace std;
using namespace Mistral;


int main(int argc, char *argv[])
{  

  Solver s;
  s.initialise();

  ReversibleMultiList<int,2> alist;

  alist.solver = &s;

  //s.initialise();

  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "add 3 to 1" << endl;
  unsigned int elt3 = alist.create(3);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "add 12 to 1" << endl;
  unsigned int elt12 = alist.create(12);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "add 5 to 1" << endl;
  unsigned int elt5 = alist.create(5);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "add 7 to 2" << endl;
  unsigned int elt7 = alist.create(7,1);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "add 8 to 2" << endl;
  unsigned int elt8 = alist.create(8,1);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "add 4 to 1" << endl;
  unsigned int elt4 = alist.create(4);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "add 9 to 2" << endl;
  unsigned int elt9 = alist.create(9,1);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "erase 12 from 1" << endl;
  alist.reversible_erase(elt12, 0);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "erase 3 from 1" << endl;
  alist.reversible_erase(elt3, 0);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "erase 8 from 2" << endl;
  alist.reversible_erase(elt8, 1);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "erase 4 from 1" << endl;
  alist.reversible_erase(elt4, 0);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;


  cout << "restore" << endl;
  //alist.restore();
  s.backtrack();
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;


  cout << "up 1 level" << endl << endl;
  //alist.first_change = true;
  s.make_node();


  cout << "erase 12 from 1" << endl;
  alist.reversible_erase(elt12, 0);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "erase 3 from 1" << endl;
  alist.reversible_erase(elt3, 0);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "up 1 level" << endl << endl;
  //alist.first_change = true;
  s.make_node();

  cout << "erase 8 from 2" << endl;
  alist.reversible_erase(elt8, 1);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "erase 9 from 2" << endl;
  alist.reversible_erase(elt9, 1);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "restore" << endl;
  //alist.restore();
  s.backtrack();
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;


  cout << "up 1 level" << endl << endl;
  //alist.first_change = true;
  s.make_node();

  cout << "add elt88 to 1" << endl;
  unsigned int elt88 = alist.reversible_add(88,0);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "erase 5 from 1" << endl;
  alist.reversible_erase(elt5, 0);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "erase 7 from 2" << endl;
  alist.reversible_erase(elt7, 1);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "add elt99 to 2" << endl;
  unsigned int elt99 = alist.reversible_add(99,1);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "up 1 level" << endl << endl;
  //alist.first_change = true;
  s.make_node();

  cout << "erase 88 from 1" << endl;
  alist.reversible_erase(elt88, 0);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "erase 4 from 1" << endl;
  alist.reversible_erase(elt4, 0);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "erase 8 from 2" << endl;
  alist.reversible_erase(elt8, 1);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "up 1 level" << endl << endl;
  //alist.first_change = true;
  s.make_node();

  cout << "erase 99 from 2" << endl;
  alist.reversible_erase(elt99, 1);
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "restore" << endl;
  //alist.restore();
  s.backtrack();
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "restore" << endl;
  //alist.restore();
  s.backtrack();
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "restore" << endl;
  //alist.restore();
  s.backtrack();
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;

  cout << "restore" << endl;
  //alist.restore();
  s.backtrack();
  alist.reversible_debug_print(cout);
  cout << endl;
  alist.flat_print(cout);
  cout << endl;
  alist.print(cout);
  cout << endl;
  cout << endl;



//   alist.insert(elt3,1);
//   alist.reversible_debug_print(cout);
//   cout << endl;
//   alist.flat_print(cout);
//   cout << endl;
//   alist.print(cout);
//   cout << endl;
//   cout << endl;


  

  //int* x_2 = int_list.erase(&all_int[2])



  //int* x_3 = int_list.erase(&all_int[3])


//   cout << "\n erase 21, 22, 49" << endl;
//   int_list.erase(21);
//   int_list.erase(22);
//   int_list.erase(49);

  
//   int_list.print(cout);
//   cout << endl;


//   cout << "\n undo" << endl;
//   int_list.fill();
  
//   int_list.print(cout);
//   cout << endl;


//   cout << "\n clear" << endl;
//   int_list.clear();
  
//   int_list.print(cout);
//   cout << endl;


//   cout << "\n insert 21, 22, 49, 2, 3" << endl;
//   int_list.insert(21);
//   int_list.insert(22);
//   int_list.insert(49);
//   int_list.insert(2);
//   int_list.insert(3);

//   int_list.print(cout);
//   cout << endl;


//   cout << "\n extend" << endl;
//   int_list.extend();

//   int_list.print(cout);
//   cout << endl;

//   cout << "\n insert 58, 90" << endl;
//   int_list.insert(58);
//   int_list.insert(90);

//   int_list.print(cout);
//   cout << endl;

//   for(unsigned int i=0; i<int_list.size; ++i)
//     cout << " " << int_list[i] ;
//   cout << endl << endl;

//   delete [] all_int;

}


