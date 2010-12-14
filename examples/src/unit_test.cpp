
#include <vector>
#include <iomanip>

#include <mistral_search.hpp>
#include <mistral_variable.hpp>
#include <mistral_constraint.hpp>


using namespace std;
using namespace Mistral;


// create X random unit tests of each type, 
// with verbosity V and thoroughness T.
// Then run them


#define LOW     1
#define MEDIUM  2
#define HIGH    3
#define EXTREME 4


class UnitTest {

protected:

  void checkDomainIntegrity(Variable X) {
    // get all values, min, max and size through iteration
    int nxt = X.get_min();
    int v=nxt;
  
    int xmin = INFTY;
    int xmax = -INFTY;
    Vector<int> values;
    int max_iteration = 10000;

    do {
      if(--max_iteration < 0) {
	cout << "Error while iterating (infinite loop?)!" << endl
	     << X << " in " << X.get_domain() << " / " << values << endl;
	exit(1);
      }
    
      for(int k=v+1; k<nxt; ++k) {
	if(X.contain(k)) {
	  cout << "Error while iterating (some values were missed)!" << endl
	       << X << " in " << X.get_domain() << " / " << values << endl;
	  exit(1);
	}
      }

      v=nxt;
      if(xmin > v) xmin = v;
      if(xmax < v) xmax = v;
      values.add(v);

      if(!X.contain(v)) {
	cout << "Error while iterating (got non-member values)!" << endl
	     << X << " in " << X.get_domain() << " / " << values << endl;
	exit(1);
      }

      nxt = X.next(v);    
    } while(nxt != v);

    if(xmin != X.get_min()) {
      cout << "lower bound not properly maintained!" << endl
	   << X << " in " << X.get_domain() << " / " << values << endl;
      exit(1);
    }

    if(xmax != X.get_max()) {
      cout << "upper bound not properly maintained!" << endl
	   << X << " in " << X.get_domain() << " / " << values << endl;
      exit(1);
    }

    if(values.size != X.get_size()) {
      cout << "domain size not properly maintained!" << endl
	   << X << " in " << X.get_domain() << " / " << values << endl;
      exit(1);
    }
  }

public:
  
  UnitTest() {
    Verbosity=LOW;
    Quality=MEDIUM;
    Quantity=MEDIUM;
  }
  
  UnitTest(const int vb, const int ql, const int qt) 
  {
    Verbosity = vb;
    Quality = ql;
    Quantity = qt;
  }


  virtual ~UnitTest() {}

  int Verbosity;
  int Quality;
  int Quantity;

  virtual void run() = 0;

};


class RandomDomainRandomRemove : public UnitTest {

public:
  
  RandomDomainRandomRemove(const int ql=MEDIUM, const int qt=MEDIUM) : UnitTest(LOW, ql, qt) {}
  ~RandomDomainRandomRemove() {}

  virtual void run() {
    /// create vars with random domains and remove a random set of values.
    Solver s;

    Vector<int> values;
    Vector<int> rand_values;

    int N = (1<<(1<<Quantity));

    int max_size = (1<<(1<<Quality)); // 2, 4, 16, 256, 65536
    int max_min = (1<<(1<<Quality));
  
    int step = (Quantity ? (1<<Quality)/(Quantity) : INFTY);

    if(Verbosity) cout << "Run random removal checks [2.." 
		       << max_size << "] x [2.." << max_min 
		       << "] (" << step << ") ";// << endl;
    if(Verbosity>LOW) cout << endl;

    for(int dom_size = 2; dom_size <= max_size ; dom_size += step) {
      if(Verbosity > LOW) {
	cout << " do " << (N*(max_min-2)/step) 
	     << " checks with domain size = " << dom_size << endl;
      }
      
      for(int dom_min = 2; dom_min <= max_min ; dom_min += step) {
      
	if(Verbosity > LOW) {
	  cout << "  do " << N << " checks with domains in [-" << (dom_min/2) 
	       << ".." << (dom_size-(dom_min/2)) << "] to [" << (dom_min/2) 
	       << ".." << (dom_size+(dom_min/2)) << "]" << endl;
	}

	for(int i=0; i<N; ++i) {
	  int lb = randint(dom_min) - (dom_min/2);
	  int ub = randint(dom_size) + lb;
	
	  int nvalues = dom_size;
	  for(int j=0; j<nvalues; ++j) 
	    rand_values.add(randint(2*dom_size)-(ub-lb+1));
	
	  Variable X(lb,ub);
	  s.add(X);
	
	  if(Verbosity > MEDIUM) {
	    cout << "    " << X << " in " << X.get_domain() << ": " << endl;
	    cout << "    remove " << rand_values << endl; 
	  }

	  for(int j=0; j<nvalues; ++j) { 

	    if(Verbosity > HIGH) {
	      cout << "      remove " << rand_values[j] << " from " 
		   << X << " in " << X.get_domain() << " => "; 
	    }
	    if( X.remove(rand_values[j]) != FAIL_EVENT) {
	      if(Verbosity > HIGH) {
		cout << "    " << X << " in " << X.get_domain() << endl; 
	      }
	  
	      checkDomainIntegrity(X);
	  
	      if(X.contain(rand_values[j])) {
		cout << "Error on [remove " << rand_values[j]
		     << " from " << X << " in " << X.get_domain() << endl;
		exit(1);
	      }
	    } else {
	      if(!X.equal(rand_values[j]))  {
		cout << "Error on [remove " << rand_values[j]
		     << " from " << X << " in " << X.get_domain() << endl;
		exit(1);
	      }
	      
	      if(Verbosity > HIGH) {
		cout << "    Wipe out!" << endl;
	      }
	    }
	  }
	
	  if(Verbosity > MEDIUM) {
	    cout << " OK (" << X << " in " << X.get_domain() << ")" << endl;
	  }

	  rand_values.clear();
	}      
      }
    }
  }
};


class RandomDomainRandomSetDomainAndRestore : public UnitTest {

public:
  
  RandomDomainRandomSetDomainAndRestore(const int ql=MEDIUM, const int qt=MEDIUM) : UnitTest(LOW, ql, qt) {}
  ~RandomDomainRandomSetDomainAndRestore() {}

  virtual void run() {
    Solver s;
    int N = (1<<(1<<Quantity));


    if(Verbosity) cout << "Run " << N << " random setDomain checks ";// << endl;
    if(Verbosity>LOW) cout << endl;

    for(int i=0; i<N; ++i) {
      int dom_size = randint(1<<(1<<Quality))+1;

      if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

      for(int j=0; j<dom_size; ++j) {

	int dom_min = randint(dom_size) - (dom_size/2);
	int lb = dom_min;
	int ub = dom_min+dom_size;
	Variable X(lb,ub);
	s.add(X);

	if(Verbosity > MEDIUM) cout << "    " << X << " in " << X.get_domain() << endl;

	for(int k=0; k<50; ++k) {

	  if(Verbosity > HIGH) cout << "      <- " << lb-k-1 ;

	  if(X.setDomain(lb-k-1) != FAIL_EVENT) {
	    cout << "Error on [setDomain] " << lb-k-1
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }
	  if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	    cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
		 << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	    exit(1);
	  }   

	  if(Verbosity > HIGH) cout << " (" << X << " in " 
				    << X.get_domain() << ")" << endl;  

	  checkDomainIntegrity(X);
	}

	for(int k=0; k<50; ++k) {

	  if(Verbosity > HIGH) cout << "      <- " << lb+k+1 ;

	  if(X.setDomain(ub+k+1) != FAIL_EVENT) {
	    cout << "Error on [setDomain] " << ub+k+1
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }
	  if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	    cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
		 << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	    exit(1);
	  }

	  if(Verbosity > HIGH) cout << " (" << X << " in " 
				    << X.get_domain() << ")" << endl;       

	  checkDomainIntegrity(X);
	}

	for(int k=0; k<dom_size; ++k) {

	  if(Verbosity > HIGH) cout << "      <- " << lb+k ;

	  X.setDomain(lb+k);
	  if(!X.equal(lb+k)) {
	    cout << "Error on [setDomain] " << lb+k
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }

	  if(Verbosity > HIGH) cout << " (" << X << " in " 
				    << X.get_domain() << ")" << endl;
	  
	  checkDomainIntegrity(X);

	  X.restore();
	  
	  if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	    cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	    exit(1);
	  }
	  
	  checkDomainIntegrity(X);
	}
      }
    }
  }


};


class RandomDomainRandomSetMinAndRestore : public UnitTest {

public:
  
  RandomDomainRandomSetMinAndRestore(const int ql=MEDIUM, const int qt=MEDIUM) : UnitTest(LOW, ql, qt) {}
  ~RandomDomainRandomSetMinAndRestore() {}

  virtual void run() {
    Solver s;
    int N = (1<<(1<<Quantity));


    if(Verbosity) cout << "Run " << N << " random setMin checks " ;// << endl;
    if(Verbosity>LOW) cout << endl;

    for(int i=0; i<N; ++i) {
      int dom_size = randint(1<<(1<<Quality))+1;

      if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

      for(int j=0; j<dom_size; ++j) {

	int dom_min = randint(dom_size) - (dom_size/2);
	int lb = dom_min;
	int ub = dom_min+dom_size;
	Variable X(lb,ub);
	s.add(X);

	if(Verbosity > MEDIUM) cout << "    " << X << " in " << X.get_domain() << endl;

	for(int k=0; k<50; ++k) {

	  if(Verbosity > HIGH) cout << "      <- " << lb-k-1 ;

	  if(X.setDomain(lb-k-1) != FAIL_EVENT) {
	    cout << "Error on [setDomain] " << lb-k-1
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }
	  if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	    cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
		 << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	    exit(1);
	  }   
	  
	  if(Verbosity > HIGH) cout << " (" << X << " in " 
				    << X.get_domain() << ")" << endl;  
	  
	  checkDomainIntegrity(X);
	}
	
	for(int k=0; k<50; ++k) {
	  
	  if(Verbosity > HIGH) cout << "      <- " << lb+k+1 ;
	  
	  if(X.setDomain(ub+k+1) != FAIL_EVENT) {
	    cout << "Error on [setDomain] " << ub+k+1
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }
	  if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	    cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
		 << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	    exit(1);
	  }

	  if(Verbosity > HIGH) cout << " (" << X << " in " 
				    << X.get_domain() << ")" << endl;       

	  checkDomainIntegrity(X);
	}

	for(int k=0; k<dom_size; ++k) {

	  if(Verbosity > HIGH) cout << "      <- " << lb+k ;

	  X.setDomain(lb+k);
	  if(!X.equal(lb+k)) {
	    cout << "Error on [setDomain] " << lb+k
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }

	  if(Verbosity > HIGH) cout << " (" << X << " in " 
				    << X.get_domain() << ")" << endl;
	  
	  checkDomainIntegrity(X);

	  X.restore();
	  
	  if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	    cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	    exit(1);
	  }
	  
	  checkDomainIntegrity(X);
	}
      }
    }
  }

};


class RandomDomainRandomSetMaxAndRestore : public UnitTest {

public:
  
  RandomDomainRandomSetMaxAndRestore(const int ql=MEDIUM, const int qt=MEDIUM) : UnitTest(LOW, ql, qt) {}
  ~RandomDomainRandomSetMaxAndRestore() {}

  virtual void run() {
    Solver s;
    int N = (1<<(1<<Quantity));


    if(Verbosity) cout << "Run " << N << " random setMax checks " ; // << endl;
    if(Verbosity>LOW) cout << endl;

    for(int i=0; i<N; ++i) {
      int dom_size = randint(1<<(1<<Quality))+1;

      if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

      for(int j=0; j<dom_size; ++j) {

	int dom_min = randint(dom_size) - (dom_size/2);
	int lb = dom_min;
	int ub = dom_min+dom_size;
	Variable X(lb,ub);
	s.add(X);

	if(Verbosity > MEDIUM) cout << "    " << X << " in " << X.get_domain() << endl;

	for(int k=0; k<50; ++k) {

	  if(Verbosity > HIGH) cout << "      <= " << lb-k-1 ;

	  if(X.setMax(lb-k-1) != FAIL_EVENT) {
	    cout << "Error on [setMax] " << lb-k-1
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }
	  if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	    cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
		 << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	    exit(1);
	  }   

	  if(Verbosity > HIGH) cout << " (" << X << " in " 
				    << X.get_domain() << ")" << endl;

	  checkDomainIntegrity(X);  
	}

	for(int k=0; k<50; ++k) {

	  if(Verbosity > HIGH) cout << "      <= " << ub+k ;

	  if(X.setMax(ub+k) != NO_EVENT) {
	    cout << "Error on [setMax] " << ub+k
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }
	  if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	    cout << "Error on [NO EVENT] " << X << " in " << X.get_domain() 
		 << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	    exit(1);
	  }

	  if(Verbosity > HIGH) cout << " (" << X << " in " 
				    << X.get_domain() << ")" << endl;

	  checkDomainIntegrity(X);       
	}

	for(int k=0; k<dom_size; ++k) {

	  if(Verbosity > HIGH) cout << "      <= " << lb+k ;

	  X.setMax(lb+k);
	  if(X.get_max() != lb+k) {
	    cout << "Error on [setMax] " << lb+k
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }

	  if(Verbosity > HIGH) cout << " (" << X << " in " 
				    << X.get_domain() << ")" << endl;
	  
	  checkDomainIntegrity(X);

	  X.restore();
	  
	  if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	    cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	    exit(1);
	  }
	  
	  checkDomainIntegrity(X);
	}
      }
    }
  }

};



class RandomDomainRandomSetDomainBitsetAndRestore : public UnitTest {

public:
  
  RandomDomainRandomSetDomainBitsetAndRestore(const int ql=MEDIUM, const int qt=MEDIUM) : UnitTest(LOW, ql, qt) {}
  ~RandomDomainRandomSetDomainBitsetAndRestore() {}

  virtual void run() {
    Solver s;
    int N = (1<<(1<<Quantity));


    if(Verbosity) cout << "Run " << N << " random setDomain (bitset) checks " ; // << endl;
    if(Verbosity>LOW) cout << endl;

    for(int i=0; i<N; ++i) {
      int dom_size = randint(1<<(1<<Quality))+1;

      if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

      for(int j=0; j<dom_size; ++j) {

	int dom_min = randint(dom_size) - (dom_size/2);
	int lb = dom_min;
	int ub = dom_min+dom_size;

	Variable X(lb,ub);
	s.add(X);

	if(Verbosity > LOW) cout << "    " << X << " in " << X.get_domain() << endl;

	BitSet vals(lb-50, ub+50, BitSet::empt);

	
	for(int k=0; k<dom_size; ++k) {

	  int set_size = randint(dom_size+10);
	  
	  for(int l=0; l<set_size; ++l) {
	    vals.add(randint(ub-lb+101) + (lb-50));
	  }
	  
	  if(Verbosity > MEDIUM) cout << "    Set: " << vals << ": ";

	  Event evt = X.setDomain(vals);

	  if(Verbosity > MEDIUM) cout << X << " in " << X.get_domain() << endl;

	  if(evt != FAIL_EVENT) {
	    int z = 0;
	    int nxt = X.get_min();
	    do {
	      z = nxt;
	      if(!vals.contain(z)) {
		cout << "Error on [setDomain] " 
		     << " on " << X << " in " << X.get_domain() << endl;
		exit(1);
	      }
	      nxt = X.next(z);
	    } while( nxt != z );
	    
	  
	    vals.flip();
	    z = 0;
	    nxt = vals.min();
	    do {
	      z = nxt;
	      if(X.contain(z)) {
		cout << "Error on [setDomain] " 
		     << " on " << X << " in " << X.get_domain() << endl;
		exit(1);
	      }
	      nxt = vals.next(z);
	    } while( nxt != z );
	    
	    checkDomainIntegrity(X);
	  
	    if(evt != NO_EVENT) X.restore();

	    checkDomainIntegrity(X);

	  }
	  
	  if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	    cout << X.get_min() << " " << lb << endl;
	    cout << X.get_max() << " " << ub << endl;
	    cout << X.get_size() << " " << (ub-lb+1) << endl;
	    cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	    exit(1);
	  }
	  
	  checkDomainIntegrity(X);
	}
      }
    }
  }

};



class RandomDomainRandomRemoveRangeAndRestore : public UnitTest {

public:
  
  RandomDomainRandomRemoveRangeAndRestore(const int ql=MEDIUM, const int qt=HIGH) : UnitTest(LOW, ql, qt) {}
  ~RandomDomainRandomRemoveRangeAndRestore() {}

  virtual void run() {
    Solver s;
    int N = (1<<(1<<Quantity));


    if(Verbosity) cout << "Run " << N << " random setDomain + removeRange checks " ;// << endl;
    if(Verbosity>LOW) cout << endl;

    for(int i=0; i<N; ++i) {
      int dom_size = randint(1<<(1<<Quality))+1;

      if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

      for(int j=0; j<dom_size; ++j) {

	int dom_min = randint(dom_size) - (dom_size/2);
	int lb = dom_min;
	int ub = dom_min+dom_size;

	Variable X(lb,ub);
	s.add(X);

	if(Verbosity > MEDIUM) cout << "    " << X << " in " << X.get_domain() << endl;

	BitSet vals(lb-50, ub+50, BitSet::empt);
	
	int set_size = randint(dom_size+10);
	  
	for(int l=0; l<set_size; ++l) {
	  vals.add(randint(ub-lb+101) + (lb-50));
	}
	  
	if(Verbosity > MEDIUM) cout << "    Set: " << vals << ": " ;
	s.save();
	
	Event evt1 = X.setDomain(vals);

	if(Verbosity > MEDIUM) cout << X << " in " << X.get_domain() << endl;

	if(evt1 != FAIL_EVENT) {	
	  int z = 0;
	  int nxt = X.get_min();
	  do {
	    z = nxt;
	    if(!vals.contain(z)) {
	      cout << "Error on [setDomain] " 
		   << " on " << X << " in " << X.get_domain() << endl;
	      exit(1);
	    }
	    nxt = X.next(z);
	  } while( nxt != z );
	  
	  
	  vals.flip();
	  z = 0;
	  nxt = vals.min();
	  do {
	    z = nxt;
	    if(X.contain(z)) {
	      cout << "Error on [setDomain] " 
		   << " on " << X << " in " << X.get_domain() << endl;
	      exit(1);
	    }
	    nxt = vals.next(z);
	  } while( nxt != z );
	  
	  checkDomainIntegrity(X);
	  
	  for(int l=lb; l<ub; ++l) {
	    for(int m=l+1; m<=ub; ++m) {
	      s.save();
	      if(Verbosity > HIGH) cout << "        remove [" << l << ".." 
					<< m << "] from " << X << " in " 
					<< X.get_domain() << ": ";
	      Event evt2 = X.removeRange(l,m);
	      if(Verbosity > HIGH) cout << X << " in " << X.get_domain() << endl;
	      if(evt2 != NO_EVENT && evt2 != FAIL_EVENT) {
		for(int v=l; v<=m; ++v) {
		  if(X.contain(v)) {
		    cout << "Error on [removeRange] " 
			 << " on " << X << " in " << X.get_domain() << endl;
		    exit(1);
		  }
		}
		checkDomainIntegrity(X);
	      }

	      s.restore();
	    }
	  }
	}
      
	s.restore();
	
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << X.get_min() << " " << lb << endl;
	  cout << X.get_max() << " " << ub << endl;
	  cout << X.get_size() << " " << (ub-lb+1) << endl;
	  cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	  exit(1);
	}
	
	checkDomainIntegrity(X);
      }
    }
  }


};

template< int NUM_HEAD >
class RandomCListRandomEraseAndRestore : public UnitTest {

public:

  int NUM_ELT;
  
  RandomCListRandomEraseAndRestore(const int ql=HIGH, const int qt=MEDIUM) 
    : UnitTest(LOW, ql, qt) { NUM_ELT = 100*(1 << (1 << Quality)); }
  ~RandomCListRandomEraseAndRestore() {}

  class NaiveImplementation : public Reversible {

  public:
    Vector<int> elements;
    Vector<int> whichlist;
    Vector<int> whichleveladded;
    Vector<int> whichlevelremoved;


    NaiveImplementation(Solver *s) 
      : Reversible(s) {}

    void create(int elt, int wl=0) {
      elements.add(elt);
      whichlist.add(wl);
      whichleveladded.add(env->level);
      whichlevelremoved.add(INFTY);
    }

    void reversible_erase(int elt, int wl) {
      //unsigned int i=0;
      //while(i<elements.size && elements[i] != elt) ++i;
      whichlevelremoved[elt] = env->level;
    }

    void reversible_add(int elt, int wl=0) {
      create(elt, wl);
    }

    void save() {}
    void restore() {}

    void print(int num_lists) {
      for(int i=0; i<num_lists; ++i) {
	cout << "[";
	for(unsigned int j=0; j<elements.size; ++j)
	  if(whichlist[j] >= i)
	    cout << " " << elements[j];
	cout << " ] ";
      }
    }

  };


  void checkEquality(ReversibleMultiList<int,NUM_HEAD>& alist, 
		     NaiveImplementation& blist) {
    BitSet inA(0, NUM_ELT, BitSet::empt);
    BitSet inB(0, NUM_ELT, BitSet::empt);

    for(int i=0; i<NUM_HEAD; ++i) {
      
      //alist.MultiList<int,NUM_HEAD>::debug_print(cout);
      //alist.reversible_debug_print(cout);
      //cout << endl;

      Node<int> nd = alist.first(i);
      while(alist.next(nd)) {
	inA.add((int)nd);
      }
      for(unsigned int j=0; j<blist.elements.size; ++j) {
// 	cout << blist.elements[j] << ": "
// 	     << blist.whichlist[j] << ">=" << i << " & "
// 	     << blist.whichlevelremoved[j] << ">" 
// 	     << blist.solver->level << " & "
// 	     << blist.whichleveladded[j] << "<=" 
// 	     << blist.solver->level << endl; 
	if(blist.whichlist[j] >= i && 
	   blist.whichlevelremoved[j] > blist.env->level &&
	   blist.whichleveladded[j] <= blist.env->level) 
	  inB.add(blist.elements[j]);
      }

      if(inA != inB) {
	cout << "Discrepancy between the lists!" << endl
	     << inA << endl
	     << inB << endl;
	exit(1);
      }//  else {
// 	cout 
// 	     << inA << endl
// 	     << inB << endl;
//       }
      
      inA.clear();
      inB.clear();
    }

    //    cout << endl;
  }

  virtual void run();
  

};


template< int NUM_HEAD >
void RandomCListRandomEraseAndRestore< NUM_HEAD >::run() {

    Solver s;
    s.initialise();

    int N = (1<<(1<<Quantity));
    if(Verbosity) cout << "Run " << N << " random ConstraintList checks " ;
    if(Verbosity>LOW) cout << endl;
    
    for(int iteration=0; iteration<N; ++iteration) {
      
      ReversibleMultiList<int,NUM_HEAD> alist;
      alist.env = &s;
      
      NaiveImplementation blist(&s);
      
      Vector<unsigned int> wl;
      
      Vector<unsigned int> elts;
      
      Vector<unsigned int> idx1;
      Vector<unsigned int> idx2;
      for(int i=0; i<NUM_ELT; ++i)
	elts.add(i+1);
      for(int i=0; i<NUM_ELT; ++i) {
	int j=randint(NUM_ELT-i)+i;
	elts[i] = elts[j];
	elts[j] = i+1;
      }
      
      
      int i_elts = (1 << (1 << Quality));
      int n_elts = i_elts;
      int head;
      for(int i=0; i<n_elts; ++i) {
	head = randint(NUM_HEAD);
	wl.add(head);

	if( i >= (int)(elts.size) ) {
	  cout << i << " " << NUM_ELT << endl;
	  exit(1);
	}

	unsigned int the_elt = alist.create( elts[i], head );
	idx1.add( the_elt );
	idx2.add( i );
	blist.create( elts[i], head );
      }
      
      checkEquality(alist, blist);
      
      
      int branch_length = randint(2*i_elts);
      
      BitSet isIn(0, NUM_ELT, BitSet::full);
      isIn.setMax(i_elts-1);
      

      if(Verbosity>LOW) cout << "Do a run on a branch of length " << branch_length << endl;

      for(int i=0; i<branch_length; ++i) {
	s.save();
	
	if(!randint(2))
	  for(int k=0; k<2; ++k) {
	    if(isIn.empty()) break;
	    
	    int j = randint(n_elts);
	    while(!isIn.contain(j))
	      j = randint(n_elts);
	    
	    alist.reversible_erase(idx1[j], wl[j]);
	    blist.reversible_erase(idx2[j], wl[j]);
	    
	    if(Verbosity>MEDIUM) {
	      for(int l=0; l<s.level; ++l) cout << " ";
	      cout << s.level << " remove " << (elts[idx2[j]]) << endl;
	      for(int l=0; l<s.level; ++l) cout << " ";
	      cout << alist << endl;
	    }

	  isIn.erase(j);
	  
	}

      if(!randint(3) && n_elts<NUM_ELT-1) {
	// add
	head = randint(NUM_HEAD);
	wl.add(head);
	idx1.add( alist.reversible_add( elts[n_elts], head ) );
	idx2.add( n_elts );
	blist.create( elts[n_elts], head );

	if(Verbosity>MEDIUM) {
	  for(int l=0; l<s.level; ++l) cout << " ";
	  cout << s.level << " add " << (elts[idx2[n_elts]]) << " to " 
	       << wl[n_elts] << "th list" << endl;
	  for(int l=0; l<s.level; ++l) cout << " ";
	  cout << alist << endl;
	}

	isIn.add(n_elts);

	++n_elts;
      }

      checkEquality(alist, blist);
    }
    
    isIn.fill();
    isIn.setMax(i_elts-1);

    for(int i=0; i<branch_length; i++) {
      s.restore();
      //int j1 = (10+(50-i-1))%n_elts;
      //int j2 = (10+(50-i-2))%n_elts;

      //cout << "add " << (elts[idx2[j1]]) << " and " << (elts[idx2[j2]]) << endl;
      if(Verbosity>MEDIUM) {
	for(int l=0; l<s.level; ++l) cout << " ";
	cout << s.level << " restore" << endl;
	for(int l=0; l<s.level; ++l) cout << " ";
	cout << alist << endl;      
      }
      
      checkEquality(alist, blist);
    }
    }
  }





template< class T >
class RandomRevNumAffectations : public UnitTest {

public:
  
  RandomRevNumAffectations(const int ql=HIGH, const int qt=MEDIUM) 
    : UnitTest(LOW, ql, qt) { }
  ~RandomRevNumAffectations() {}

  virtual void run() {
    if(Verbosity) cout << "Run " << 10 << " random Reversible<int> checks " ;

    Solver s;
    T x = 1000;
    int val = 0;

    ReversibleNum<T> rx(x, &s);

    int x1 = x;
    int x2;

    for(int k=0; k<10; ++k) {
      s.save();
	
      for(int i=0; i<10; ++i) {
	rx += ((i+k)%20)-10;
	for(int j=0; j<5; ++j)
	  val += rx;
      }

      x2 = rx;
      
      for(int i=0; i<9999990; ++i) {
	if(!(i%10))s.save();
	rx += ((i+k)%20)-10;
	for(int j=0; j<5; ++j)
	  val += rx;
      }
      
      for(int i=0; i<9999990; i+=10) {
	s.restore();
      }
      
      if(x2 != rx) {
	cout << "Error: discrepancy between the values! " 
	     << x << " " << rx << endl;
	exit(1);
      }

      s.restore();

      if(x1 != rx) {
	cout << "Error: discrepancy between the values! " 
	     << x << " " << rx << endl;
	exit(1);
      }
    }
  }
  
};


class CostasAllDiffSolutions : public UnitTest {

public:
  
  int size;
  
  CostasAllDiffSolutions(const int sz) 
    : UnitTest(LOW, LOW, LOW) { size=sz; }
  ~CostasAllDiffSolutions() {}

  virtual void run() {
    Vector< Variable > X;
    int i, j;

    for(i=0; i<size; ++i) X.add( Variable(1,size) );
    
    cout << X << endl;

    Solver s;

    s.add( AllDiff(X) );

    cout << s << endl;


    Vector< Variable > difference[size-2];
    Vector< Variable > scope;
    for(i=1; i<size-1; ++i) {
      for(j=i; j<size; ++j) {
	Variable Y(1-size, size-1);

	s.add(Y);
	difference[i-1].add(Y);
	scope.add(Y);
	scope.add(X[j]);
	scope.add(X[j-i]);

	Constraint *con = new PredicateAdd(scope);
	con->initialise();

	s.add(con);
	scope.clear();
      }

      
      Constraint *con = new ConstraintAllDiff(difference[i-1]);
      con->initialise();

      s.add(con);
    }

//     Vector< Variable > scope;
//     for(int i=1; i<size-1; ++i) {
//       //differences at distance i
//       scope.clear();
//       for(int j=0; j<size-i; ++j) {
// 	scope.add(X[j]-X[j+i]);
//       }
//       s.add( AllDiff(scope) );
//     }
  
    cout << s << endl;

    s.initialise();

    s.depth_first_search(X);

    for(int i=0; i<size; ++i)
      cout << setw(3) << X[i].get_solution_value() << " " ;
    cout << endl;
   for(int i=1; i<size; ++i) {
     for(int j=0; j<size-i; ++j) {
       cout << setw(3) << (X[j].get_solution_value()-
			   X[j+i].get_solution_value())
	    << " " ;
      }
     cout << endl;
    }
  }

};


class CostasNotEqualSolutions : public UnitTest {

public:
  
  int size;
  
  CostasNotEqualSolutions(const int sz) 
    : UnitTest(LOW, LOW, LOW) { size=sz; }
  ~CostasNotEqualSolutions() {}

  virtual void run() {
    Vector< Variable > X;
    for(int i=0; i<size; ++i) X.add( Variable(1,size) );
    
    cout << X << endl;

    Solver s;
    int i, j, k;
    for(i=0; i<size; ++i)
      for(j=i+1; j<size; ++j)
	s.add( X[i] != X[j] );

    cout << s << endl;

    Vector< Variable > difference[size-2];
    Vector< Variable > scope;
    for(i=1; i<size-1; ++i) {
      for(j=i; j<size; ++j) {
	Variable Y(1-size, size-1);

	s.add(Y);
	difference[i-1].add(Y);
	scope.add(Y);
	scope.add(X[j]);
	scope.add(X[j-i]);

	Constraint *con = new PredicateAdd(scope);
	con->initialise();

	s.add(con);
	scope.clear();
      }
      for(j=1; j<size-i; ++j)
	for(k=0; k<j; ++k) {
	  scope.add(difference[i-1][j]);
	  scope.add(difference[i-1][k]);
	  
	  Constraint *con = new ConstraintNotEqual(scope);
	  con->initialise();
	  
	  s.add(con);
	  scope.clear();
	}
    }

    cout << s << endl;

    s.initialise();

    s.depth_first_search(X, 
			 new GenericHeuristic< GenericDVO< MinDomain >, MinValue >(&s), 
			 new Geometric(256, 1.1));

    for(int i=0; i<size; ++i)
      cout << setw(3) << X[i].get_solution_value() << " " ;
    cout << endl;
   for(int i=1; i<size; ++i) {
     for(int j=0; j<size-i; ++j) {
       cout << setw(3) << (X[j].get_solution_value()-
			   X[j+i].get_solution_value())
	    << " " ;
      }
     cout << endl;
    }
  }

};


class PigeonsSolutions : public UnitTest {

public:
  
  int size;
  
  PigeonsSolutions(const int sz) 
    : UnitTest(LOW, LOW, LOW) { size=sz; }
  ~PigeonsSolutions() {}

  virtual void run() {
    Vector< Variable > X;
    for(int i=0; i<size; ++i) X.add( Variable(1,size-1) );
    
    cout << X << endl;

    Solver s;
    int i, j;
    for(i=0; i<size; ++i)
      for(j=i+1; j<size; ++j)
	s.add( X[i] != X[j] );

    cout << s << endl;

    s.initialise();

    s.depth_first_search(X, 
			 new GenericHeuristic< NoOrder, MinValue >(&s), 
			 new NoRestart);

    cout << s.statistics << endl;
  }

};


int main(int argc, char *argv[])
{  

  usrand(12345);

  std::vector<UnitTest*> tests;
 
  tests.push_back(new PigeonsSolutions(atoi(argv[1])));
  /*
  tests.push_back(new CostasNotEqualSolutions(atoi(argv[1]))); 
  tests.push_back(new PigeonsSolutions(atoi(argv[1])));
 tests.push_back(new CostasAllDiffSolutions(atoi(argv[1])));
  tests.push_back(new RandomCListRandomEraseAndRestore<4>());
  tests.push_back(new RandomDomainRandomRemoveRangeAndRestore());
  tests.push_back(new RandomDomainRandomSetDomainBitsetAndRestore());
  tests.push_back(new RandomDomainRandomSetDomainAndRestore());
  tests.push_back(new RandomDomainRandomSetMaxAndRestore());
  tests.push_back(new RandomDomainRandomSetMinAndRestore());
  tests.push_back(new RandomDomainRandomRemove());
  tests.push_back(new RandomRevNumAffectations<int>());
  */
  //tests[0]->Verbosity = HIGH;
  //tests[0]->Quality = HIGH;
  //tests[0]->Quantity = EXTREME;
  //tests[1]->Verbosity = EXTREME;

  while(tests.size() > 0) {
    UnitTest *t = tests.back();
    
    double TIME = getRunTime();
    t->run();
    cout << (getRunTime() - TIME) << endl ;
    
    delete t;
    tests.pop_back();
  }

}


