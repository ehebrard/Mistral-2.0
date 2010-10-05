
#include <mistral_variable.hpp>
#include <climits>
#include <iomanip>


using namespace std;
using namespace Mistral;





int main(int argc, char *argv[])
{  
  int N = (argc > 1 ? atoi(argv[1]) : 8);

  IntVar* pigeons = new IntVar[N];
  
  for(int i=0; i<N; ++i)
    //pigeons[i] = new IntegerVar(1,N-1);
    pigeons[i] = new VariableWord< unsigned int, 1 >(1,N-1);

  Solver s;

  Vector< IntVar > scope(2);
  for(int i=1; i<N; ++i) 
    for(int j=0; j<i; ++j) {
      scope.add(pigeons[i]);
      scope.add(pigeons[j]);
      s.add(new ConstraintNotEqual(scope));
      scope.clear();
    }

  s.initialise();
  s.depth_first_search();
  cout << s.statistics << endl;
}


