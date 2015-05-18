
#include <assert.h> 
#include <stdlib.h>
#include <fstream>


#include <mistral_solver.hpp>
#include <mistral_variable.hpp>
#include <mistral_search.hpp>



using namespace std;
using namespace Mistral;


class Puzzle {

public:
	int height; // j, y
	int width;  // i, x
	int surface;
	
	BitSet* grid;
	
	Vector<int> squares;
	Vector<int> trail;
	
	int **limit; // used to compute only non-dominated squares
	int **corner_of; // for each cell, the size of its max square
	IntStack **covered_by; // for each cell, the id of the max squares that cover it
	Vector<int> max_y; // for each maximum square, the y coordinate of its top-left corner 
	Vector<int> max_x; // for each maximum square, the x coordinate of its top-left corner 
	
	BitSet set_cover; // set of maximum squares that might be enough to cover the whole shape
	IntStack unique_squares;
	
	IntStack *choices;
	
	
	
	Puzzle() {
		width = 0;
		height = 0;
		grid = NULL;
		limit = NULL;
		corner_of = NULL;
		covered_by = NULL;
		choices = NULL;
	}
	
	virtual ~Puzzle() {
		if(grid) {
			delete [] grid;
			delete [] choices;
			for(int j=0; j<height; ++j) {
				delete [] covered_by[j];
				delete [] corner_of[j];
				delete [] limit[j];
			}
			delete [] covered_by;
			delete [] corner_of;
			delete [] limit;
		}
	}
	
	
	inline int num_squares() { return squares.size/3; }
	
	inline void save() { trail.add(squares.size); }
	
	inline void add_square(const int y, const int x, const int c) {
		
		// // cout << "add square " << y << "." << x << ": " << c << endl;
		//
		// write(cout);
		//
		//
		// int rs = 0;
		// for(int j=0; j<height; ++j) {
		// 	for(int i=0; i<width; ++i) {
		// 		rs += grid[j].contain(i);
		// 	}
		// }
		// assert(rs == surface);
		
		if(c>1) {
			for(int j=y; j<y+c; ++j) {
				
				// cout << (x >> BitSet::EXP) << " " << ((x+c-1) >> BitSet::EXP) << endl;
				//
				// grid[j].print_bits(cout);
				//
				// cout << " xor [" << x << ", " << (x+c-1) << "] = ";
				// cout.flush();
				
				// cout << grid[j] << " xor [" << x << "," << (x+c-1) << "] = ";
				grid[j].xor_with(x,x+c-1);
				// cout << grid[j] << endl;
				
				// grid[j].print_bits(cout);
				//
				// cout << endl;
			}
		} else {
			grid[y].remove(x);
		}
		
		surface -= (c*c);
	
		squares.add(y);
		squares.add(x);
		squares.add(c);
		
		//write(cout);
		
		//  rs = 0;
		// for(int j=0; j<height; ++j) {
		// 	for(int i=0; i<width; ++i) {
		// 		rs += grid[j].contain(i);
		// 	}
		// }
		// assert(rs == surface);
		
	}
	
	inline void restore() {
		//
		// cout << trail << endl;
		// cout << squares << endl;
		//
		int l = trail.pop();
		
		do {
			int c = squares.pop();
			int x = squares.pop();
			int y = squares.pop();
			
				//
			// cout << "undo " << y << "." << x << ": " << c << endl;
			//
			surface += (c*c);
		
			if(c>1) {
				for(int j=y; j<y+c; ++j) {
					grid[j].union_with(x,x+c-1);
				}
			} else {
				grid[y].add(x);
			}		
		} while(squares.size>l);
		
		
		// int rs = 0;
		// for(int j=0; j<height; ++j) {
		// 	for(int i=0; i<width; ++i) {
		// 		rs += grid[j].contain(i);
		// 	}
		// }
		// assert(rs == surface);
	}
	
	
	void init_maximum_squares() {
		corner_of = new int*[height];
		limit = new int*[height];
		choices = new IntStack[height*width+1];
		covered_by = new IntStack*[height];
		for(int j=0; j<height; ++j) {
			corner_of[j] = new int[width];
			limit[j] = new int[width];
			fill(corner_of[j], corner_of[j]+width, -1);
			fill(limit[j], limit[j]+width, 0);
			
			covered_by[j] = new IntStack[width];
			for(int i=0; i<width; ++i) {
					choices[j*width+i].initialise(0,height*width-1,height*width,true);
					covered_by[j][i].initialise(0,height*width-1,height*width,false);
			}
		}
		choices[height*width].initialise(0,height*width-1,height*width,true);
		set_cover.initialise(0,height*width-1,BitSet::empt);
		unique_squares.initialise(0,height*width-1,height*width,false);
	}
	
	
	void clear_maximum_squares() {
		for(int j=0; j<height; ++j) {
			for(int i=0; i<width; ++i) {
				covered_by[j][i].clear();
			}
			fill(limit[j], limit[j]+width, 0);
		}
		while(!max_x.empty()) {
			int y = max_y.pop();
			int x = max_x.pop();
			int c = corner_of[y][x];
			corner_of[y][x] = -1;
		}
		unique_squares.clear();
		set_cover.clear();
	}
	
	int maximise_square(const int j, const int i) {
		int c = 1;
		
		// cout << "  " << c << endl;
		while(((j+c) < height) && ((i+c) < width)) {
			
			
			bool ok = true;
			for(int x=i; x<=i+c; ++x) {
				if(!grid[j+c].contain(x)) {
					ok =false;
					break;
				}
			}
			if(ok) {
				for(int y=j; y<=j+c; ++y) {
					if(!grid[y].contain(i+c)){
						ok = false;
						break;
					}
				}
			}
			if(ok) {
				++c;
			} else {
				break;
			}
			
			// cout << "  " << c << endl;
		}
		return c;
	}
	
	void compute_maximum_squares() {
		
		// cout << "compute max squares" << endl;
		
		int id_sq = 0;
		for(int j=0; j<height; ++j) {
			for(int i=0; i<width; ++i) {
				
				// cout << j << "," << i << endl;
				
				if(grid[j].contain(i)) {
					if(((width-i) > limit[j][i]) && ((height-j) > limit[j][i])) {
						int c = maximise_square(j, i);
						
						// cout << c << " > " << limit[j][i] << endl << endl;
                
						if(c>limit[j][i]) {
							
							max_x.add(i);
							max_y.add(j);
							corner_of[j][i] = c; 
							bool unique = true;
							
							for(int y=j; y<j+c; ++y) {
								for(int x=i; x<i+c; ++x) {
									if(covered_by[y][x].size>0) {
										unique = false;
										for(int k=0; k<covered_by[y][x].size; ++k)
											if(unique_squares.contain(covered_by[y][x][k]))
												unique_squares.remove(covered_by[y][x][k]);
									}
									covered_by[y][x].add(id_sq);
									limit[y][x]=max(limit[y][x], min(j+c-y, i+c-x));
								}
							}
							
							if(unique)
								unique_squares.add(id_sq);
							
							++id_sq;
						}
					}
				}
			}
		}
	}
	

	
	int reduce() {
		int k, y, x;
		for(int u=0; u<unique_squares.size; ++u) {
			k = unique_squares[u];
			//cout << " (" << corner_of[max_y[k]][max_x[k]] << ", " << max_y[k] << ", " << max_x[k] << ")";
			y = max_y[k];
			x = max_x[k];
			add_square(y, x, corner_of[y][x]);
			covered_by[y][x].remove(k);
		}
		// for(int k=max_x.size-1; k>=0; --k) {
		//
		// }

		return unique_squares.size;
	}
	
	int compute_cover() {
		int cover_size = 0;
		for(int j=0; j<height; ++j) {
			for(int i=0; i<width; ++i) {
				if(covered_by[j][i].size==1) {
					int sq = covered_by[j][i][0];
					
					if(!set_cover.contain(sq)) {
						set_cover.add(sq);
						++cover_size;
					}
				}
			}
		}
		return cover_size;
	}
	
	
	void read(const char* filename) {
		
		cout << "read from " << filename << endl;
		
		ifstream infile(filename, ios_base::in);
		char c;
		
		surface = 0;
		
		infile >> height;
		infile >> width;
		infile.ignore(1000, '\n');
		
		grid = new BitSet[height];
		
		for(int j=0; j<height; ++j) {
			grid[j].initialise(0,width-1,BitSet::empt);
			for(int i=0; i<width; ++i) {
				c = infile.get();
				cout << c; 
				if(c != ' ') {
					grid[j].add(i);
					++surface;
				}
			}
			cout << endl;
			infile.ignore(1000, '\n');
		}
		
		infile.close();
		
		init_maximum_squares();
	}
	
	
	ostream& write( ostream& os ) {
		
		// cout << "print puzzle" << endl;
		
		int **buffer = new int*[height];
		for(int j=0; j<height; ++j) {
			buffer[j] = new int[width];
			fill(buffer[j], buffer[j]+width, -2);
		}
		
		int y,x,c;
		for(int k=0; k<squares.size;) {
			y = squares[k++];
			x = squares[k++];
			c = squares[k++];
			
			if(c==1) buffer[y][x] = -1;
			else {
				for(int j=y; j<y+c; ++j) {
					for(int i=x; i<x+c; ++i) {
						buffer[j][i]=((k/3)%10);
					}
				}
			}
		}
		
		for(int j=0; j<height; ++j) {
			for(int i=0; i<width; ++i) {
				if(buffer[j][i]==-1) {
					assert(!grid[j][i]);
					os << "X";
				} else if(buffer[j][i]>-1) {
					os << buffer[j][i];
					assert(!grid[j][i]);
				} else if(grid[j][i]) {
					os << ".";
				} else {
					os << " ";
				}
			}
			os << endl;
		}
		
		// cout << "unique:" ;
		// for(int u=0; u<unique_squares.size; ++u) {
		// 	int k = unique_squares[u];
		// 	cout << " (" << corner_of[max_y[k]][max_x[k]] << ", " << max_y[k] << ", " << max_x[k] << ")";
		// }
		// cout << endl << "multi:";
		// for(int k=0; k<max_x.size; ++k) {
		// 	if(!unique_squares.contain(k))
		// 		cout << " (" << corner_of[max_y[k]][max_x[k]] << ", " << max_y[k] << ", " << max_x[k] << ")";
		// }
		// cout << endl;
		// for(int k=0; k<max_x.size; ++k) {
		// 	cout << "(" << corner_of[max_y[k]][max_x[k]] << ", " << max_y[k] << ", " << max_x[k] << ") ";
		// }
		// cout << endl;
		// for(int j=0; j<height; ++j) {
		// 	for(int i=0; i<width; ++i) {
		// 		if(grid[j].contain(i))
		// 			cout << j << "," << i << ": " << covered_by[j][i] << endl;
		// 	}
		// }
		// cout << endl;
		
		return os;
		
	}
	
	
	void print_choices() {
		for(int k=0; k<=trail.size; ++k) {
			cout << choices[k] << endl;
		}
		
	}

};


// void perm(const int n) {
//
// 	// int limit = 10;
//
// 	Vector<int> trail;
// 	//int *num_branches = new int[n+1];
// 	//std::fill(num_branches, num_branches+n+1, 0);
// 	IntStack sequence;
// 	sequence.initialise(0, n-1, n, false);
//
// 	bool finished=false;
//
// 	while(!finished) {
//
// 		// //if(--limit<0) break;
// 		//
// 		cout << endl;
// 		for(int i=0; i<trail.size; ++i) {
// 			cout << " " << setw(2) << trail[i] ;
// 		}
// 		// cout << endl;
// 		// for(int i=0; i<sequence.size; ++i) {
// 		// 	cout << " " << setw(2) << num_branches[i] ;
// 		// }
// 		// cout << " |";
// 		// for(int i=sequence.size; i<n; ++i) {
// 		// 	cout << " " << setw(2) << num_branches[i] ;
// 		// }
// 		cout << endl;
// 		for(int i=0; i<sequence.size; ++i) {
// 			cout << " " << setw(2) << sequence[i] ;
// 		}
// 		cout << " |";
// 		for(int i=sequence.size; i<n; ++i) {
// 			cout << " " << setw(2) << sequence[i] ;
// 		}
// 		cout << endl;
//
//
// 		if(sequence.size+num_branches[trail.size]>=n) {
// 			if(sequence.size==n) {
// 				cout << "sol: ";
// 				for(int i=0; i<n; ++i) {
// 					cout << " " << setw(2) << sequence[i] ;
// 				}
// 				cout << endl;
// 			}
//
// 			if(trail.empty()) {
// 				finished = true;
// 			} else {
//
// 				cout << endl << "backtrack to level " << trail.size-1 ;
//
// 				num_branches[trail.size] = 0;
// 				sequence.size = trail.pop();
// 				++num_branches[trail.size];
// 			}
// 		} else {
// 			cout << endl << "decision at level " << trail.size << " (" << (num_branches[trail.size]+1) << "th)" ;
//
// 			trail.add(sequence.size);
// 			sequence.add(sequence.back(num_branches[trail.size-1]));
// 		}
// 	}
//
// 	delete [] num_branches;
// }



void dfs(Puzzle& p, const int verbosity) {
	
	//int limit = 2000000;
	
	//usrand(12345);
	
	// int shuffle = new int[p.height*p.width];
	// for(int k=0; k<p.height*p.width; ++k) {
	// 	shuffle[k] = k;
	// }
	// for(int k=0; k<p.height*p.width; ++k) {
	// 	int r = randint(p.height*p.width-k);
	// 	int aux = shuffle[k+r];
	// 	shuffle[k+r] = shuffle[k];
	// 	shuffle[]
	// }
	
	
	int ub = -1, lb;
	int y, x, k, l, r;
	bool solved = false;
	
	//p.init_maximum_squares();
	
	while(!solved) {
		
		//if(--limit<0) break;
		
		p.clear_maximum_squares();
		p.compute_maximum_squares();
		p.reduce();
		
		if(verbosity>2)
			cout << p.trail << endl;
		if(verbosity>1)
			p.write(cout);
		
		lb = p.compute_cover();
		
		if(verbosity>1)
			cout << "#sq=" << p.num_squares() << " + lb=" << lb << " <> " << ub << endl << endl;
		
		if(ub>=0 && (p.num_squares() + lb) >= ub) {
			
			if(p.trail.size) {
				
				assert(p.trail.size < p.height*p.width);
				
				p.choices[p.trail.size].fill();
				p.restore();
				y = p.squares.back(0);
				x = p.squares.back(-1);
				
				if(verbosity>2)
					cout << "backtrack " << y << " " << x << endl;
				
				
				p.choices[p.trail.size].remove(y*p.width+x);
			} else {
				solved = true;
				
				if(verbosity) {
					cout << "solved!" << endl;
				}
			}
		} else if(p.surface) {
			
			r = randint(p.max_x.size);
			//
			
			//cout << endl;
			for(l=0; l<p.max_x.size; ++l) {
			 	k = ((l+r)%(p.max_x.size));
				
				cout << k 
					//<< " " << p.max_x.size << " in " << p.choices[p.trail.size] << " \\ " << p.unique_squares 
					<< endl;
				
				if(!p.unique_squares.contain(k)) {
					
					cout << " not unique " << p.max_y << endl;
					
					assert(k<p.max_x.size);
					assert(k<p.max_y.size);
					
					y = p.max_y[k];
					x = p.max_x[k];
					
					
					cout << y << " " << x << endl;
					
					assert(p.trail.size < p.height*p.width);
					
					
					cout  << (y*p.width+x) << " in "<< p.choices[p.trail.size] << " ?" << endl;
					
					if(p.choices[p.trail.size].contain(y*p.width+x)) {
						
						if(verbosity>2)
							cout << "decision " << y << " " << x << " in " << p.choices[p.trail.size] << " \\ " << p.unique_squares << endl;

						p.save();
						
						p.add_square(y,x,p.corner_of[y][x]);
											
						break;
					}
					
				}
			}
			
			if(l>=p.max_x.size) {
				
				if(p.trail.size) {
					
					assert(p.trail.size < p.height*p.width);
					
					p.choices[p.trail.size].fill();
					p.restore();
					y = p.squares.back(0);
					x = p.squares.back(-1);
			
					if(verbosity>2)
						cout << "backtrack " << y << " " << x << endl;
			
					
			
					p.choices[p.trail.size].remove(y*p.width+x);
				} else {
					solved = true;
				
					if(verbosity) {
						cout << "solved!" << endl;
					}
				}
			
			}
			
		} else {
			
			if(verbosity) {
				cout << "solution #sq=" << p.num_squares() << endl;
				p.write(cout);
			}
			
			ub = p.num_squares();
			
			assert(p.trail.size < p.height*p.width);
			
			p.choices[p.trail.size].fill();
			p.restore();
			y = p.squares.back(0);
			x = p.squares.back(-1);
			
			if(verbosity>2)
				cout << "backtrack " << y << " " << x << endl;
			
			p.choices[p.trail.size].remove(y*p.width+x);
						
		}
	}
}



int main(int argc, char **argv)
{

	Puzzle p;
	int cs;
	
	p.read("ex1.txt");
	//p.init_maximum_squares();
	
	//perm(4);
	
	
	dfs(p, 3);
	
	
	// p.init_maximum_squares();
	// p.clear_maximum_squares();
	// p.compute_maximum_squares();
	// p.reduce();
	// cs = p.compute_cover();
	//
	// p.write(cout);
	// cout << cs << endl;
	//
	// p.save();
	// p.add_square(0,2,2);
	// p.clear_maximum_squares();
	// p.compute_maximum_squares();
	// p.reduce();
	// cs = p.compute_cover();
	//
	// p.write(cout);
	// cout << cs << endl;
	//
	// p.save();
	// p.add_square(2,3,3);
	// p.clear_maximum_squares();
	// p.compute_maximum_squares();
	// p.reduce();
	// cs = p.compute_cover();
	//
	// p.write(cout);
	// cout << cs << endl;
	//
	// p.restore();
	// p.clear_maximum_squares();
	// p.compute_maximum_squares();
	// cs = p.compute_cover();
	//
	// p.write(cout);
	// cout << cs << endl;
	//
	// p.restore();
	// p.clear_maximum_squares();
	// p.compute_maximum_squares();
	// cs = p.compute_cover();
	//
	// p.write(cout);
	// cout << cs << endl;


}



