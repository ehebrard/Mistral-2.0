
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
	
	BitSet grid;
	
	Vector<int> squares;
	Vector<int> rectangles;
	Vector<int> trail;
	
	int **limit; // used to compute only non-dominated squares
	int **corner_of; // for each cell, the size of its max square
	IntStack **covered_by; // for each cell, the id of the max squares that cover it
	Vector<int> max_y; // for each maximum square, the y coordinate of its top-left corner 
	Vector<int> max_x; // for each maximum square, the x coordinate of its top-left corner 
	
	BitSet set_cover; // set of maximum squares that might be enough to cover the whole shape
	IntStack unique_squares;
	
	IntStack *choices;
	
	Vector<int> best_solution;
	
	
	
	Puzzle() {
		width = 0;
		height = 0;
		limit = NULL;
		corner_of = NULL;
		covered_by = NULL;
		choices = NULL;
	}
	
	virtual ~Puzzle() {
		if(limit) {
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
	
	inline void add_max_square(const int y, const int x) {
		
		//cout << "add "<< y << " " << x << endl;
		
		int i, j, c, extendable;
		if(grid.fast_contain(y*width+x)) {
			c = 1;
			
			while(x+c<width && y+c<height) {
				
				// cout << grid[y+c] << " >= [" << x << ", " << (x+c) << "]?" << endl;
				
				if(!grid.includes((y+c)*width+x,(y+c)*width+x+c)) break;
				
				// cout << "yes!" << endl;
					
				for(j=y; j<y+c && grid.fast_contain(j*width+x+c); ++j);// {};
				if(j<y+c) break;
				
				// cout << "ok for the bottom row" << endl;
				
				++c;
			}
		
			// cout << "  --> " << c << endl;
				
			if(c>1) {
				for(j=y; j<y+c; ++j) {
					grid.xor_with(j*width+x,j*width+x+c-1);
				}
			} else {
				grid.fast_remove(y*width+x);
			}
		
			surface -= (c*c);
	
			squares.add(y);
			squares.add(x);
			squares.add(c);
		}
	}
	
	
	inline void add_square(const int y, const int x, const int c) {
		
		 // cout << "add square " << y << "." << x << ": " << c << endl;
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
		//  cout << surface << endl;
		
		if(c>1) {
			for(int j=y; j<y+c; ++j) {
				

				 // grid[j].print_bits(cout);
				 //
				 //  cout << " xor [" << x << ", " << (x+c-1) << "] = ";
				 // 					 cout.flush();
			
				 
				grid.xor_with(j*width+x,j*width+x+c-1);
	
				 // grid[j].print_bits(cout);
				 //  cout << endl;
			}
		} else {
			
			// cout << y << " " << grid[y] << " - " << x << " = ";
			// cout.flush();
			
			// assert(grid[y].contain(x));
			
			grid.fast_remove(y*width+x);
			
			// cout << grid[y] << endl;
		}
		
		surface -= (c*c);
	
		squares.add(y);
		squares.add(x);
		squares.add(c);
		
		//write(cout);
		
		// cout << surface << endl;
		//  rs = 0;
		// for(int j=0; j<height; ++j) {
		// 	for(int i=0; i<width; ++i) {
		// 		rs += grid[j].contain(i);
		// 	}
		// }
		// assert(rs == surface);
		// cout << surface << endl;
		
	}
	
	inline void restore() {
		
		// int rs = 0;
		// for(int j=0; j<height; ++j) {
		// 	for(int i=0; i<width; ++i) {
		// 		rs += grid[j].contain(i);
		// 	}
		// }
		// assert(rs == surface);
		
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
					grid.union_with(j*width+x,j*width+x+c-1);
				}
			} else {
				grid.fast_add(y*width+x);
			}		
		} while(squares.size>l);
		
		
		//  rs = 0;
		// for(int j=0; j<height; ++j) {
		// 	for(int i=0; i<width; ++i) {
		// 		rs += grid[j].contain(i);
		// 	}
		// }
		// assert(rs == surface);
	}
	
	void apply(Vector<int>& somesquares, const bool propagate=true) {
		int y, x;
		if(propagate)
			propagate_and_bound(-1);
		int nsq = squares.size;
		for(int k=0; k<somesquares.size; k+=3) {
			y = somesquares[k];
			x = somesquares[k+1];
			add_max_square(y,x);
			
			if(nsq<squares.size) {
				if(propagate)
					propagate_and_bound(-1);
				nsq = squares.size;
				//write(cout);
			}
		}
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
		rectangles.clear();
	}
	
	int maximise_square(const int j, const int i) {
		int c = 1;
		
		// cout << "  " << c << endl;
		while(((j+c) < height) && ((i+c) < width)) {
			
			
			bool ok = true;
			for(int x=i; x<=i+c; ++x) {
				if(!grid.contain((j+c)*width+x)) {
					ok =false;
					break;
				}
			}
			if(ok) {
				for(int y=j; y<=j+c; ++y) {
					if(!grid.contain(y*width+i+c)){
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
				
				if(grid.contain(j*width+i)) {
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
	
	bool is_enclosed(const int y1, const int x1, const int y2, const int x2) {
		
		bool enclosed = true;
		
		// cout << "check " << y1 << ", " << x1 << ", " << y2 << ", " << x2 << endl;
		
		// check top
		if(y1) {
			for(int i=x1; enclosed && i<x2; ++i) {
				enclosed = (!grid.contain((y1-1)*width+i) || !grid.contain((y1-1)*width+i+1));
			}
		}
		
		if(enclosed) {
			
			// cout << "ok top" << endl;
			
			// check bottom
			if(y2<height-1) {
				for(int i=x1; enclosed && i<x2; ++i) {
					enclosed = (!grid.contain((y2+1)*width+i) || !grid.contain((y2+1)*width+i+1));
				}
			}
		
			if(enclosed) {
				
				// cout << "ok bottom" << endl;
				
				// check left
				if(x1) {
					for(int j=y1; enclosed && j<y2; ++j) {
						enclosed = (!grid.contain(j*width+x1-1) || !grid.contain(j*width+x1-1));
					}
				}
		
				if(enclosed) {
					
					// cout << "ok left" << endl;
					
					// check right
					if(x2<width-1) {
						for(int j=y1; enclosed && j<y2; ++j) {
							//cout << "  (" << j << ", " << (j+1) <<")" <<endl;
							enclosed = (!grid.contain(j*width+x2+1) || !grid.contain(j*width+x2+1));
						}
					}
					
					if(enclosed) {
					
					// if(enclosed) cout << "ok right" << endl;
					
						// check top left corner
						if(y1 && x1) {
							enclosed = (!grid.contain((y1-1)*width+x1-1) || !grid.contain(y1*width+x1-1) || !grid.contain((y1-1)*width+x1));
						}
						
						if(enclosed) {
					
						// if(enclosed) cout << "ok right" << endl;
					
							// check bottom left corner
							if(y2<height-1 && x1) {
								enclosed = (!grid.contain((y2+1)*width+x1-1) || !grid.contain(y2*width+x1-1) || !grid.contain((y2+1)*width+x1));
							}
						
							if(enclosed) {
					
							// if(enclosed) cout << "ok right" << endl;
					
								// check top right corner
								if(y1 && x2<width-1) {
									enclosed = (!grid.contain((y1-1)*width+x2+1) || !grid.contain(y1*width+x2+1) || !grid.contain((y1-1)*width+x2));
								}
						
								if(enclosed) {
					
								// if(enclosed) cout << "ok right" << endl;
					
									// check bottom right corner
									if(y2<height-1 && x2<width-1) {
										enclosed = (!grid.contain((y2+1)*width+x2+1) || !grid.contain(y2*width+x2+1) || !grid.contain((y2+1)*width+x2));
									}
						
						
								}	
							}	
						}	
					}	
				}
			}
		}
		
		return enclosed;
	}
	
	
	void compute_rectangles() {
		int y, x, c, a, b, extended;
		for(int k=0; k<max_x.size; ++k) {
			y = max_y[k];
			x = max_x[k];
			c = corner_of[y][x];
			
			if(c>1) {
			
			//cout << "extend " << y << " " << x << endl;
			
			// try to extend it by adding rows
			extended = true;
			for(a=y+c-1; extended && a<height-1;) {
				for(int i=x; extended && i<x+c; ++i) {
					extended = grid.contain((a+1)*width+i);
				}
				if(extended) ++a;
			}
			extended = true;
			for(b=y; extended && b>0;) {
				for(int i=x; extended && i<x+c; ++i) {
					extended = grid.contain((b-1)*width+i);
				}
				if(extended) --b;
			}
			
			if(
				 b<y || 
			a>=y+c) {
				int l = rectangles.size;
				if(l==0 || rectangles[l-1]<(x+c-1) || rectangles[l-2]<a || rectangles[l-3]>x || rectangles[l-4]>b) {
				
					//cout << "extend by adding rows!" << b << " " << x << " " << a << " " << (x+c-1) << endl;
				
					if(is_enclosed(b,x,a,x+c-1)) {
					rectangles.add(b);
				//rectangles.add(y);
				rectangles.add(x);
				rectangles.add(a);
				rectangles.add(x+c-1);
			}
			}
			}
			
			// try to extend it by adding columns
			extended = true;
			for(a=x+c-1; extended && a<width-1;) {
				for(int j=y; extended && j<y+c; ++j) {
					extended = grid.contain(j*width+a+1);
				}
				if(extended) ++a;
			}
			extended = true;
			for(b=x; extended && b>0;) {
				for(int j=y; extended && j<y+c; ++j) {
					extended = grid.contain(j*width+b-1);
				}
				if(extended) --b;
			}

			
			if(
				 b<x ||
			a>=x+c) {
				
				int l = rectangles.size;
				if(l==0 || rectangles[l-1]<a || rectangles[l-2]<(y+c-1) || rectangles[l-3]>b || rectangles[l-4]>y) {
					
					//cout << "extend by adding columns! " << y << " " << b << " " << (y+c-1) << " " << a << endl;
					
					if(is_enclosed(y,b,y+c-1,a)) {
					rectangles.add(y);
					rectangles.add(b);
				//rectangles.add(x);
				rectangles.add(y+c-1);
				rectangles.add(a);
			}
			}
			}
			
		}
		
	}
	}
	
	
	int propagate_and_bound(const int ub) {
		
		bool ufilter = true;
		bool rfilter = true;
		
		while(ufilter || rfilter) {
			clear_maximum_squares();

			compute_maximum_squares();
			
			if(rfilter) {
				reduce_unique();
			}
			
			if(ufilter) {
				compute_rectangles();
				reduce_rectangle();
			}
			
			ufilter = !unique_squares.empty();
			rfilter = !rectangles.empty();
		}
		
		return (ub>=0 ? compute_cover() : 0);
	}

	
	int reduce_unique() {
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
	
	int reduce_rectangle() {
		int y1, x1, y2, x2;
		for(int u=0; u<rectangles.size; u+=4) {
			y1 = rectangles[u];
			x1 = rectangles[u+1];
			y2 = rectangles[u+2];
			x2 = rectangles[u+3];

			
			while(y1<=y2 && x1<=x2) {
				if(y2-y1 > x2-x1) {
					add_square(y1, x1, (x2-x1+1));
					y1 = y1+x2-x1+1;
				} else {
					add_square(y1, x1, (y2-y1+1));
					x1 = x1+y2-y1+1;
				}
			}

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
		
		grid.initialise(0,height*width-1,BitSet::empt);
		
		for(int j=0; j<height; ++j) {
			for(int i=0; i<width; ++i) {
				c = infile.get();
				//cout << c; 
				if(c != ' ') {
					grid.add(j*width+i);
					++surface;
				}
			}
			//cout << endl;
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
		
		int y,x,c,r=1;
		for(int k=0; k<squares.size;) {
			y = squares[k++];
			x = squares[k++];
			c = squares[k++];
			
			if(c==1) buffer[y][x] = -1;
			else {
				for(int j=y; j<y+c; ++j) {
					for(int i=x; i<x+c; ++i) {
						buffer[j][i]=(r%10);
					}
				}
				++r;
			}
		}
		
		for(int j=0; j<height; ++j) {
			for(int i=0; i<width; ++i) {
				if(buffer[j][i]==-1) {
					assert(!grid[j*width+i]);
					os << "X";
				} else if(buffer[j][i]>-1) {
					os << buffer[j][i] ;
					assert(!grid[j*width+i]);
				} else if(grid[j*width+i]) {
					os << ".";
				} else {
					os << " ";
				}
			}
			os << endl;
		}
		
		// cout << max_x.size << endl;
		// cout << rectangles.size << endl;
		// for(int k=0; k<rectangles.size; k+=4) {
		//
		//
		//
		// 	int y1 = rectangles[k];
		// 	int x1 = rectangles[k+1];
		// 	int y2 = rectangles[k+2];
		// 	int x2 = rectangles[k+3];
		//
		// 	cout << "(" << y1 << "," << x1 << "," << y2 << "," << x2 << ")\n";
		//
		// 	for(int j=y1; j<=y2; ++j) {
		// 		for(int i=x1; i<=x2; ++i) {
		// 			if(buffer[j][i] != -2) {
		// 				cout << j << " " << i << ": " << buffer[j][i] << endl;
		// 				exit(1);
		// 			}
		// 			assert(buffer[j][i] == -2);
		// 			buffer[j][i] = 0;
		// 		}
		// 	}
		//
		// 	for(int j=0; j<height; ++j) {
		// 		for(int i=0; i<width; ++i) {
		// 			if(buffer[j][i]==-1) {
		// 				//assert(!grid[j][i]);
		// 				os << "X";
		// 			} else if(buffer[j][i]>-1) {
		// 				os << buffer[j][i];
		// 				//assert(!grid[j][i]);
		// 			} else if(grid[j][i]) {
		// 				os << ".";
		// 			} else {
		// 				os << " ";
		// 			}
		// 		}
		// 		os << endl;
		// 	}
		//
		//
		// 	for(int j=y1; j<=y2; ++j) {
		// 		for(int i=x1; i<=x2; ++i) {
		// 			buffer[j][i] = -2;
		// 		}
		// 	}
		//
		// }
		// cout << endl;
		// //exit(1);
		
		
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
	
	ostream& print_squares(ostream& os) {
		int y,x,c;
		os << "squares";
		for(int k=0; k<best_solution.size;) {
			y = best_solution[k++];
			x = best_solution[k++];
			c = best_solution[k++];
			
			os << " " << x << " " << y << " " << c;
		}
		os << endl;
		return os;
	}

};


class Nogood {

public:
	// int height;
	// int width;
	int value;

	BitSet grid;

	// Nogood(Puzzle& p, const int ub) {
	// 	// height = p.height;
	// 	// width = p.width;
	//
	// 	value = ub-p.num_squares();
	//
	// 	grid.initialise(p.grid);
	// }
	
	Nogood() {};
	
	void initialise(const int size) {
		value = 0;
		grid.initialise(0, size-1, BitSet::empt);	
	}
	
	virtual ~Nogood() {}
	
	void store(Puzzle& p, const int ub) {
		value = ub-p.num_squares();
		grid.copy(p.grid);
	}
	
	Nogood& operator=(const Nogood& n) {
     grid.copy(n.grid);
		 value = n.value;
		 return *this;
  }


	bool operator==(const Nogood& n) {
      return grid.equal(n.grid);
  }
	
	bool operator==(const Puzzle& p) {
      return grid.equal(p.grid);
  }

	
	// inline bool subsume(Puzzle& p, const int ub) {
	//
	// 	if(p.num_squares() + value >= ub) {
	// 		if(grid.includes(p.grid)) {
	// 		 return true;
	// 		}
	// 	}
	//
	// 	return false;
	// }
	//
	//
	// inline bool subsume(Nogood& n) {
	//
	// 	if(grid.includes(n.grid)) {
	// 	 return true;
	// 	}
	//
	// 	return false;
	// }


	ostream& write( ostream& os, const int height, const int width ) {
		os << value << endl;
		for(int j=0; j<height; ++j) {
			for(int i=0; i<width; ++i) {
				if(grid[j*width+i]) {
					os << ".";
				} else {
					os << " ";
				}
			}
			os << endl;
		}
		
		return os;
	}


};



class nPuzzle {

public:
	int height; // j, y
	int width;  // i, x
	int surface;
	
	BitSet* grid;
	
	Vector<int> squares;
	Vector<int> rectangles;
	Vector<int> trail;
	
	int **limit; // used to compute only non-dominated squares
	int **corner_of; // for each cell, the size of its max square
	IntStack **covered_by; // for each cell, the id of the max squares that cover it
	Vector<int> max_y; // for each maximum square, the y coordinate of its top-left corner 
	Vector<int> max_x; // for each maximum square, the x coordinate of its top-left corner 
	
	BitSet set_cover; // set of maximum squares that might be enough to cover the whole shape
	IntStack unique_squares;
	
	IntStack *choices;
	
	Vector<int> best_solution;
	
	
	
	nPuzzle() {
		width = 0;
		height = 0;
		grid = NULL;
		limit = NULL;
		corner_of = NULL;
		covered_by = NULL;
		choices = NULL;
	}
	
	virtual ~nPuzzle() {
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
	
	inline void add_max_square(const int y, const int x) {
		
		//cout << "add "<< y << " " << x << endl;
		
		int i, j, c, extendable;
		if(grid[y].fast_contain(x)) {
			c = 1;
			
			while(x+c<width && y+c<height) {
				
				// cout << grid[y+c] << " >= [" << x << ", " << (x+c) << "]?" << endl;
				
				if(!grid[y+c].includes(x,x+c)) break;
				
				// cout << "yes!" << endl;
					
				for(j=y; j<y+c && grid[j].fast_contain(x+c); ++j);// {};
				if(j<y+c) break;
				
				// cout << "ok for the bottom row" << endl;
				
				++c;
			}
		
			// cout << "  --> " << c << endl;
				
			if(c>1) {
				for(j=y; j<y+c; ++j) {
					grid[j].xor_with(x,x+c-1);
				}
			} else {
				grid[y].fast_remove(x);
			}
		
			surface -= (c*c);
	
			squares.add(y);
			squares.add(x);
			squares.add(c);
		}
	}
	
	
	inline void add_square(const int y, const int x, const int c) {
		
		 // cout << "add square " << y << "." << x << ": " << c << endl;
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
		//  cout << surface << endl;
		
		if(c>1) {
			for(int j=y; j<y+c; ++j) {
				

				 // grid[j].print_bits(cout);
				 //
				 //  cout << " xor [" << x << ", " << (x+c-1) << "] = ";
				 // 					 cout.flush();
			
				 
				grid[j].xor_with(x,x+c-1);
	
				 // grid[j].print_bits(cout);
				 //  cout << endl;
			}
		} else {
			
			// cout << y << " " << grid[y] << " - " << x << " = ";
			// cout.flush();
			
			// assert(grid[y].contain(x));
			
			grid[y].fast_remove(x);
			
			// cout << grid[y] << endl;
		}
		
		surface -= (c*c);
	
		squares.add(y);
		squares.add(x);
		squares.add(c);
		
		//write(cout);
		
		// cout << surface << endl;
		//  rs = 0;
		// for(int j=0; j<height; ++j) {
		// 	for(int i=0; i<width; ++i) {
		// 		rs += grid[j].contain(i);
		// 	}
		// }
		// assert(rs == surface);
		// cout << surface << endl;
		
	}
	
	inline void restore() {
		
		// int rs = 0;
		// for(int j=0; j<height; ++j) {
		// 	for(int i=0; i<width; ++i) {
		// 		rs += grid[j].contain(i);
		// 	}
		// }
		// assert(rs == surface);
		
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
				grid[y].fast_add(x);
			}		
		} while(squares.size>l);
		
		
		//  rs = 0;
		// for(int j=0; j<height; ++j) {
		// 	for(int i=0; i<width; ++i) {
		// 		rs += grid[j].contain(i);
		// 	}
		// }
		// assert(rs == surface);
	}
	
	void apply(Vector<int>& somesquares, const bool propagate=true) {
		int y, x;
		if(propagate)
			propagate_and_bound(-1);
		int nsq = squares.size;
		for(int k=0; k<somesquares.size; k+=3) {
			y = somesquares[k];
			x = somesquares[k+1];
			add_max_square(y,x);
			
			if(nsq<squares.size) {
				if(propagate)
					propagate_and_bound(-1);
				nsq = squares.size;
				write(cout);
			}
		}
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
		rectangles.clear();
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
	
	bool is_enclosed(const int y1, const int x1, const int y2, const int x2) {
		
		bool enclosed = true;
		
		// cout << "check " << y1 << ", " << x1 << ", " << y2 << ", " << x2 << endl;
		
		// check top
		if(y1) {
			for(int i=x1; enclosed && i<x2; ++i) {
				enclosed = (!grid[y1-1].contain(i) || !grid[y1-1].contain(i+1));
			}
		}
		
		if(enclosed) {
			
			// cout << "ok top" << endl;
			
			// check bottom
			if(y2<height-1) {
				for(int i=x1; enclosed && i<x2; ++i) {
					enclosed = (!grid[y2+1].contain(i) || !grid[y2+1].contain(i+1));
				}
			}
		
			if(enclosed) {
				
				// cout << "ok bottom" << endl;
				
				// check left
				if(x1) {
					for(int j=y1; enclosed && j<y2; ++j) {
						enclosed = (!grid[j].contain(x1-1) || !grid[j+1].contain(x1-1));
					}
				}
		
				if(enclosed) {
					
					// cout << "ok left" << endl;
					
					// check right
					if(x2<width-1) {
						for(int j=y1; enclosed && j<y2; ++j) {
							//cout << "  (" << j << ", " << (j+1) <<")" <<endl;
							enclosed = (!grid[j].contain(x2+1) || !grid[j+1].contain(x2+1));
						}
					}
					
					if(enclosed) {
					
					// if(enclosed) cout << "ok right" << endl;
					
						// check top left corner
						if(y1 && x1) {
							enclosed = (!grid[y1-1].contain(x1-1) || !grid[y1].contain(x1-1) || !grid[y1-1].contain(x1));
						}
						
						if(enclosed) {
					
						// if(enclosed) cout << "ok right" << endl;
					
							// check bottom left corner
							if(y2<height-1 && x1) {
								enclosed = (!grid[y2+1].contain(x1-1) || !grid[y2].contain(x1-1) || !grid[y2+1].contain(x1));
							}
						
							if(enclosed) {
					
							// if(enclosed) cout << "ok right" << endl;
					
								// check top right corner
								if(y1 && x2<width-1) {
									enclosed = (!grid[y1-1].contain(x2+1) || !grid[y1].contain(x2+1) || !grid[y1-1].contain(x2));
								}
						
								if(enclosed) {
					
								// if(enclosed) cout << "ok right" << endl;
					
									// check bottom right corner
									if(y2<height-1 && x2<width-1) {
										enclosed = (!grid[y2+1].contain(x2+1) || !grid[y2].contain(x2+1) || !grid[y2+1].contain(x2));
									}
						
						
								}	
							}	
						}	
					}	
				}
			}
		}
		
		return enclosed;
	}
	
	
	void compute_rectangles() {
		int y, x, c, a, b, extended;
		for(int k=0; k<max_x.size; ++k) {
			y = max_y[k];
			x = max_x[k];
			c = corner_of[y][x];
			
			if(c>1) {
			
			//cout << "extend " << y << " " << x << endl;
			
			// try to extend it by adding rows
			extended = true;
			for(a=y+c-1; extended && a<height-1;) {
				for(int i=x; extended && i<x+c; ++i) {
					extended = grid[a+1].contain(i);
				}
				if(extended) ++a;
			}
			extended = true;
			for(b=y; extended && b>0;) {
				for(int i=x; extended && i<x+c; ++i) {
					extended = grid[b-1].contain(i);
				}
				if(extended) --b;
			}
			
			if(
				 b<y || 
			a>=y+c) {
				int l = rectangles.size;
				if(l==0 || rectangles[l-1]<(x+c-1) || rectangles[l-2]<a || rectangles[l-3]>x || rectangles[l-4]>b) {
				
					//cout << "extend by adding rows!" << b << " " << x << " " << a << " " << (x+c-1) << endl;
				
					if(is_enclosed(b,x,a,x+c-1)) {
					rectangles.add(b);
				//rectangles.add(y);
				rectangles.add(x);
				rectangles.add(a);
				rectangles.add(x+c-1);
			}
			}
			}
			
			// try to extend it by adding columns
			extended = true;
			for(a=x+c-1; extended && a<width-1;) {
				for(int j=y; extended && j<y+c; ++j) {
					extended = grid[j].contain(a+1);
				}
				if(extended) ++a;
			}
			extended = true;
			for(b=x; extended && b>0;) {
				for(int j=y; extended && j<y+c; ++j) {
					extended = grid[j].contain(b-1);
				}
				if(extended) --b;
			}

			
			if(
				 b<x ||
			a>=x+c) {
				
				int l = rectangles.size;
				if(l==0 || rectangles[l-1]<a || rectangles[l-2]<(y+c-1) || rectangles[l-3]>b || rectangles[l-4]>y) {
					
					//cout << "extend by adding columns! " << y << " " << b << " " << (y+c-1) << " " << a << endl;
					
					if(is_enclosed(y,b,y+c-1,a)) {
					rectangles.add(y);
					rectangles.add(b);
				//rectangles.add(x);
				rectangles.add(y+c-1);
				rectangles.add(a);
			}
			}
			}
			
		}
		
	}
	}
	
	
	int propagate_and_bound(const int ub) {
		
		bool ufilter = true;
		bool rfilter = true;
		
		while(ufilter || rfilter) {
			clear_maximum_squares();

			compute_maximum_squares();
			
			if(rfilter) {
				reduce_unique();
			}
			
			if(ufilter) {
				compute_rectangles();
				reduce_rectangle();
			}
			
			ufilter = !unique_squares.empty();
			rfilter = !rectangles.empty();
		}
		
		return (ub>=0 ? compute_cover() : 0);
	}

	
	int reduce_unique() {
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
	
	int reduce_rectangle() {
		int y1, x1, y2, x2;
		for(int u=0; u<rectangles.size; u+=4) {
			y1 = rectangles[u];
			x1 = rectangles[u+1];
			y2 = rectangles[u+2];
			x2 = rectangles[u+3];

			
			while(y1<=y2 && x1<=x2) {
				if(y2-y1 > x2-x1) {
					add_square(y1, x1, (x2-x1+1));
					y1 = y1+x2-x1+1;
				} else {
					add_square(y1, x1, (y2-y1+1));
					x1 = x1+y2-y1+1;
				}
			}

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
		
		int y,x,c,r=1;
		for(int k=0; k<squares.size;) {
			y = squares[k++];
			x = squares[k++];
			c = squares[k++];
			
			if(c==1) buffer[y][x] = -1;
			else {
				for(int j=y; j<y+c; ++j) {
					for(int i=x; i<x+c; ++i) {
						buffer[j][i]=(r%10);
					}
				}
				++r;
			}
		}
		
		for(int j=0; j<height; ++j) {
			for(int i=0; i<width; ++i) {
				if(buffer[j][i]==-1) {
					assert(!grid[j][i]);
					os << "X";
				} else if(buffer[j][i]>-1) {
					os << buffer[j][i] ;
					assert(!grid[j][i]);
				} else if(grid[j][i]) {
					os << ".";
				} else {
					os << " ";
				}
			}
			os << endl;
		}
		
		// cout << max_x.size << endl;
		// cout << rectangles.size << endl;
		// for(int k=0; k<rectangles.size; k+=4) {
		//
		//
		//
		// 	int y1 = rectangles[k];
		// 	int x1 = rectangles[k+1];
		// 	int y2 = rectangles[k+2];
		// 	int x2 = rectangles[k+3];
		//
		// 	cout << "(" << y1 << "," << x1 << "," << y2 << "," << x2 << ")\n";
		//
		// 	for(int j=y1; j<=y2; ++j) {
		// 		for(int i=x1; i<=x2; ++i) {
		// 			if(buffer[j][i] != -2) {
		// 				cout << j << " " << i << ": " << buffer[j][i] << endl;
		// 				exit(1);
		// 			}
		// 			assert(buffer[j][i] == -2);
		// 			buffer[j][i] = 0;
		// 		}
		// 	}
		//
		// 	for(int j=0; j<height; ++j) {
		// 		for(int i=0; i<width; ++i) {
		// 			if(buffer[j][i]==-1) {
		// 				//assert(!grid[j][i]);
		// 				os << "X";
		// 			} else if(buffer[j][i]>-1) {
		// 				os << buffer[j][i];
		// 				//assert(!grid[j][i]);
		// 			} else if(grid[j][i]) {
		// 				os << ".";
		// 			} else {
		// 				os << " ";
		// 			}
		// 		}
		// 		os << endl;
		// 	}
		//
		//
		// 	for(int j=y1; j<=y2; ++j) {
		// 		for(int i=x1; i<=x2; ++i) {
		// 			buffer[j][i] = -2;
		// 		}
		// 	}
		//
		// }
		// cout << endl;
		// //exit(1);
		
		
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
	
	ostream& print_squares(ostream& os) {
		int y,x,c;
		os << "squares";
		for(int k=0; k<best_solution.size;) {
			y = best_solution[k++];
			x = best_solution[k++];
			c = best_solution[k++];
			
			os << " " << x << " " << y << " " << c;
		}
		os << endl;
		return os;
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


class Algo {

public:
	
	Nogood tmp;
	Nogood *base;
	
	int ng_limit;
	int nb_nogoods;
	
	Algo(Puzzle& p, const int nl) : ng_limit(nl) {
		nb_nogoods = 0;
		
		tmp.initialise(p.width*p.height);
		
		base = new Nogood[ng_limit];
		for(int k=0; k<ng_limit; ++k) {
			base[k].initialise(p.width*p.height);
		}
		
	}
	
	virtual ~Algo() { delete [] base; }
	
	
	

bool dfs(Puzzle& p, int& ub, const int limit, const bool randomized, const int verbosity) {
	
	//int limit = 2000000;

	
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
	
	
	int lb;
	int y, x, k, l, r=0;
	bool solved = false;
	
	//p.init_maximum_squares();
	
	int iteration = 0;
	while(!solved) {
		
		if(++iteration >= limit) break;
		
		// p.clear_maximum_squares();
		// p.compute_maximum_squares();
		// p.reduce_unique();
		// p.compute_rectangles();
		
		lb = p.propagate_and_bound(ub);
		
		
		// cout << "check if a nogood applies" << endl;
		for(int k=min(ng_limit,nb_nogoods)-1; k>=0; --k) {
			if( base[k] == p ) {
				
					// cout << "yes: " << endl;
				
				// base[k].write(cout, p.height, p.width);
				
				// cout << endl;
				// p.write(cout);
					
				// cout << endl << p.num_squares() << " + " << base[k].value << " >= " << ub << "?\n" ;
				//
				// exit(1);
				
				if(base[k].value > lb)
					lb = base[k].value;
			}
		}

		
		if(verbosity>3)
			cout << p.trail << endl;
		if(verbosity>2)
			p.write(cout);
		
		// p.reduce_rectangle();
		// p.rectangles.clear();
		
		if(verbosity>2)
			p.write(cout);
		
		
		// lb = p.compute_cover();
		
		if(verbosity>1) {
			for(int l=0; l<p.trail.size; ++l) cout << " ";
			cout << "#sq=" << p.num_squares() << " + lb=" << lb << " <> " << ub << endl;
		}
		
		if(ub>=0 && (p.num_squares() + lb) >= ub) {
			
			if(p.trail.size) {
				
				p.choices[p.trail.size].fill();
				p.restore();
				y = p.squares.back(0);
				x = p.squares.back(-1);
				
				if(verbosity>1) {
			for(int l=0; l<p.trail.size; ++l) cout << " ";
					cout << "backtrack " << y << " " << x << endl;
				}
				
				p.choices[p.trail.size].remove(y*p.width+x);
			} else {
				solved = true;
				
				if(verbosity) {
					cout << "solved!" << endl;
				}
			}
		} else if(p.surface) {
			
			if(randomized && p.max_x.size>1 && p.trail.size<=10)
				//r = randint(p.max_x.size);
				r = randint(2);
			else 
				r = 0;
			//
			
			//cout << endl;
			for(l=0; l<p.max_x.size; ++l) {
			//for(l=p.max_x.size-1; l>=0; --l) {
				if(r>0) {
					k = ((l+r)%(p.max_x.size));
				} else {
					k = l;
				}
	
				if(!p.unique_squares.contain(k)) {
						
					y = p.max_y[k];
					x = p.max_x[k];
					
					if(p.choices[p.trail.size].contain(y*p.width+x)) {
						
						if(verbosity>1) {
			for(int l=0; l<p.trail.size; ++l) cout << " ";
							cout << "decision " << y << " " << x << " in " << p.choices[p.trail.size] << " \\ " << p.unique_squares << endl;
						}
						
						p.save();
						
						p.add_square(y,x,p.corner_of[y][x]);
											
						break;
					}
					
				}
			}
			
			if(l>=p.max_x.size) {
				
				if(p.trail.size) {
					
					//assert(p.trail.size < p.height*p.width);
					
					
					// cout << "backtrack after having exhausted all possible choices for level " << p.trail.size << " creating a new nogood" << endl;
					
					tmp.store(p, ub);
					
					// tmp.write(cout, p.height, p.width);
					//
					// cout << "check if it already exists" << endl;
					//
					bool exists = false;
					for(int k=min(nb_nogoods,ng_limit)-1; !exists && k>=0; --k) {
						exists = (base[k] == tmp);
					}
					
					if(!exists) {
						base[(nb_nogoods++)%ng_limit] = tmp;
					}
					
					
					// //cout << "check if it is subsumed "
					//
					// if(nb_nogood < ng_limit) {
					// 	base[nb_nogood++].store()
					// }
					//
					// cout << endl;
					
					
					p.choices[p.trail.size].fill();
					p.restore();
					y = p.squares.back(0);
					x = p.squares.back(-1);
			
					if(verbosity>1) {
			for(int l=0; l<p.trail.size; ++l) cout << " ";
						cout << "backtrack " << y << " " << x << endl;
					}
					
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
			 			for(int l=0; l<p.trail.size; ++l) cout << " ";
				cout << "solution #sq=" << p.num_squares() << " " << iteration << " iterations " << endl;
				
				if(verbosity>1)
					p.write(cout);
				
				p.best_solution.clear();
				for(int k=0; k<p.squares.size; ++k) {
					p.best_solution.add(p.squares[k]);
				}
				//p.print_squares(cout);
				//cout << endl;
			}
			
			ub = p.num_squares();
			
			p.choices[p.trail.size].fill();
			p.restore();
			y = p.squares.back(0);
			x = p.squares.back(-1);
			
			if(verbosity>1) {
			for(int l=0; l<p.trail.size; ++l) cout << " ";
				cout << "backtrack " << y << " " << x << endl;
			}
			
			p.choices[p.trail.size].remove(y*p.width+x);
						
		}
	}
	
	while(!p.trail.empty()) {
		p.restore();
	}
	
	return solved;
}


void restart(Puzzle &p, const int olimit, const int verbosity) {
	
	int total = 0;
	
	int limit = olimit/4;
	int ub = -1;
	bool randomized = false;
	
	bool solved = false;
	
	while(!solved && total<olimit) {
		if(verbosity)
			cout << "ub = " << ub << ", restart for " << limit << " iterations" << endl; 
		solved = dfs(p, ub, limit, randomized, verbosity);
		
		total += limit;
		
		if(!randomized) {
			randomized = true;
			limit = 100;
		} else {
			limit = (int)((double)limit * 1.2);
			if(limit > (olimit-total)) {
				limit = olimit-total;
			}
		}
	}
	
}


// void ls(Puzzle &p) {
//
// 	Vector<int> solution;
// 	for(int j=0; j<p.height; ++j) {
// 		for(int i=0; i<p.width; ++i) {
// 			if(p.grid[j].contain(i)) {
// 				solution.add(j);
// 				solution.add(i);
// 				solution.add(-1);
// 			}
// 		}
// 	}
//
// 	p.apply(solution);
//
//
//
// }


};

int main(int argc, char **argv)
{
	usrand(12345);


	Puzzle p;
	int cs;
	
	const char* filename = "ex1.txt";
	if(argc>1) {
		filename = argv[1];
	}
	int limit = 100000;
	if(argc>2) {
		limit = atoi(argv[2]);
	}
	
	p.read(filename);
	//p.init_maximum_squares();
	
	//perm(4);
	
	int ub = -1;
	
	Algo a(p, 500);
	
	//a.dfs(p, ub, limit, false, 1);
	
	a.restart(p, limit, 1);
	//

	//cout << p.best_solution.size/3 << endl;
	p.print_squares(cout);
	cout << endl;
	
	// ls(p);
	//
	// cout << p.num_squares() << endl;
	
	//p.write(cout);
	//cout << endl;
	
	
	// p.write(cout);
	// p.save();
	
	
	p.apply(p.best_solution);
	//
	// cout << p.squares.size << endl;
	//
	
	
	p.write(cout);
	cout << endl;
	// p.restore();
	//
	//
	// for(int k=0; k<p.best_solution.size; k+=3) {
	// 	cout << k << " " << p.best_solution[k] << " " << p.best_solution[k+1] << " " << p.best_solution[k+2] << endl;
	// }
	//
	//
	// p.best_solution[262] = 15;
	// p.best_solution[268] = 13;
	//
	// p.best_solution.add(19);
	// p.best_solution.add(13);
	// p.best_solution.add(-1);
	//
	// p.best_solution.add(19);
	// p.best_solution.add(14);
	// p.best_solution.add(-1);
	//
	//
	// p.save();
	// p.apply(p.best_solution);
	//
	// cout << p.squares.size << endl;
	//
	// p.write(cout);
	// p.restore();
	
	
	
	// for(int k=0; k<p.best_solution.size; ++k) {
	//
	// }
	
	
	// p.init_maximum_squares();
	// p.clear_maximum_squares();
	// p.compute_maximum_squares();
	// p.reduce_unique();
	// cs = p.compute_cover();
	//
	// p.write(cout);
	// cout << cs << endl;
	//
	// p.save();
	// p.add_square(0,2,2);
	// p.clear_maximum_squares();
	// p.compute_maximum_squares();
	// p.reduce_unique();
	// cs = p.compute_cover();
	//
	// p.write(cout);
	// cout << cs << endl;
	//
	// p.save();
	// p.add_square(2,3,3);
	// p.clear_maximum_squares();
	// p.compute_maximum_squares();
	// p.reduce_unique();
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



