
#include <mistral_variable.hpp>
#include <climits>
#include <iomanip>


using namespace std;
using namespace Mistral;





int main(int argc, char *argv[])
{  
  int N = (argc > 1 ? atoi(argv[1]) : 8);


  Vector< VariableInt< unsigned int > > pigeons;
  

  for(int i=0; i<N/2; ++i) {
    VariableWord< unsigned int, 1 > pigeon(1,N-1);
    pigeons.add(pigeon);
  }

  for(int i=N/2; i<N; ++i) {
    VariableInt< unsigned int > pigeon(1,N-1);
    pigeons.add(pigeon);
  }

  cout << "P: " << pigeons << endl;

  exit(1);


  Solver s;

  Vector< IntVar > scope(2);
  for(int i=1; i<N; ++i) 
    for(int j=0; j<i; ++j) {
      scope.add(pigeons.stack_+i);
      scope.add(pigeons.stack_+j);
      s.add(new ConstraintNotEqual(scope));
      scope.clear();
    }

  s.initialise();
  //s.depth_first_search(NULL, new Geometric());
  s.depth_first_search();
  std::cout << s.statistics << std::endl;

}


