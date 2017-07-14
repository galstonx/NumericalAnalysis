namespace NumericalAnalysis {




  template<typename T>
  Matrix<T>::Matrix(unsigned long long r,unsigned long long c) : n_rows(r), n_cols(c), data(r,std::vector<T>(c)) {

  }


  template<typename T>
  unsigned long long Matrix<T>::rows() const {
    return n_rows;
  }

  template<typename T>
  unsigned long long Matrix<T>::cols() const {
    return n_cols;
  }


  template<typename T>
  T Matrix<T>::get(unsigned long long i,unsigned long long j) const  {
    return data[i][j];
  }

  template<typename T>
  int Matrix<T>::set(unsigned long long i,unsigned long long j,T val) {
    data[i][j]=val;
    return 0;
  }

  template<typename T>
  int Matrix<T>::swapRows(unsigned long long i,unsigned long long j) {
    std::vector<T> tmp_row(data[i]);
    data[i]=data[j];
    data[j]=tmp_row;

    return 0;
  }


  // replace row 2 with c*row1+row2
  template<typename T>
  int Matrix<T>::rowOperation(T c,unsigned long long row1,unsigned long long row2) {
    for(unsigned long long i=0;i<n_cols;++i) {
      data[row2][i]=data[row2][i]+c*data[row1][i];
    }
    return 0;
  }

  template<typename T>
  int Matrix<T>::rowReduce() {
    unsigned long long p_col=0;
    unsigned long long p_row=0;
    while(p_row<n_rows) {
      T p_val=data[p_row][p_col];
      // pivot is not 0, use it to 0 out column
      if(p_val!=0) {
	for(unsigned long long j=p_row+1;j<n_rows;++j) {
	  // only need to 0 out if nonzero
	  T x=data[j][p_col];
	  if(x!=0) {
	    //rowOperation(-data[j][p_col]/p_val,p_row,j);
	    data[j][p_col]=0;
	    for(unsigned long long k=p_col+1;k<n_cols;++k) {
	      data[j][k]=data[j][k]-x/p_val*data[p_row][k];
	    }
	  }
	}
	// now increment p_row and p_col
	p_row++;
	p_col++;
      }
      // pivot is 0, need to try to swap rows
      else {
	bool swapped=false;
	// check for nonzero entry to swap with
	for(unsigned long long j=p_row+1;j<n_rows && !swapped;++j) {
	  if(data[j][p_col]!=0) {
	    swapRows(p_row,j);
	    swapped=true;
	  }
	}
	// if couldn't swap, increment pivot column
	// either way don't increment p_row
	if(!swapped) p_col++;
      }
    }
    return 0;
  }

  template<typename T>
  int Matrix<T>::GaussBackSolve(Point<T>& ans) {
    int rv=0;
    unsigned long long p_row=n_rows-1;
    unsigned long long p_col=0;
    unsigned long long last_p_col=n_cols-2;
    bool has_solution=false;
    bool done=false;
    while(! done) {
      // first find the pivot
      T p_val=get(p_row,p_col);
      while(p_val==0 && p_col<last_p_col) {
	p_col++;
	p_val=get(p_row,p_col);
      }
      // if no pivot, check if inconsistent
      if( p_col==n_cols-1 ) {
	if(get(p_row,n_cols-1)!=0) { // inconsistent
	  done=true;
	  rv=1;
	}
      }
      else { // found a pivot, update answer
	T x=get(p_row,n_cols-1);
	for(unsigned long long i=p_col+1;i<n_cols-1;++i) {
	  x=x-get(p_row,i)*ans[i];
	}
	x=x/p_val;
	ans[p_col]=x;
      }
      // move on to next row
      if(p_row==0) {
	done=true;
      }
      else {
	p_row--;
	p_col=0;
      }
    }
    return rv;
  }


  template <typename T>
  std::ostream& operator<<(std::ostream& str,const Matrix<T>& matrix) {
    for(unsigned long long i=0;i<matrix.n_rows;++i) {
      for(unsigned long long j=0;j<matrix.n_cols;++j) {
	str << matrix.get(i,j) << " ";
      }
      str << std::endl;
    }
    return str;
  }


}
