
#include <mistral_backtrack.hpp>
#include <climits>
#include <iomanip>


using namespace std;
using namespace Mistral;



int main(int argc, char *argv[])
{  

  int lb[5] = {1,2,0,1,3};
  int ub[5] = {3,3,3,4,5};

  MultiSet m(-1,3,lb,ub);

  cout << m << endl;

  cout << "add 1" << endl;
  m.insert(1);
  cout << m << endl;
  cout << endl;

  cout << "add 2" << endl;
  m.insert(2);
  cout << m << endl;
  cout << endl;

  cout << "add 2" << endl;
  m.insert(2);
  cout << m << endl;
  cout << endl;

  cout << "remove 3" << endl;
  m.erase(3);
  cout << m << endl;
  cout << endl;

  cout << "remove 1" << endl;
  m.erase(1);
  cout << m << endl;
  cout << endl;

  cout << "remove -1" << endl;
  m.erase(-1);
  cout << m << endl;
  cout << endl;

}


