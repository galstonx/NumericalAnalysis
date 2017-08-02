
namespace NumericalAnalysis {

  template<typename T>
  int PoissonBasicMethod2d<T>::init_matrix() {
    unsigned long long counter=0;
    for(unsigned long long i=0;i<n_y_pts;++i) {
      for(unsigned long long j=0;j<n_x_pts;++j) {
	p_matrix.set(counter,counter,4);
	if(j!=n_x_pts-1) p_matrix.set(counter,counter+1,-1);
	if(j!=0) p_matrix.set(counter,counter-1,-1);
	if(i!=n_y_pts-1) {
	  p_matrix.set(counter+n_x_pts,counter,-1);
	}
	if(i!=0) {
	  p_matrix.set(counter-n_x_pts,counter,-1);
	}
	++counter;
      }
    }
    return 0;
  }

  template<typename T>
  int PoissonBasicMethod2d<T>::init_nonhomog_vector() {
    unsigned long long counter=0;
    T rv;
    Point<T> pt(2);
    for(unsigned long long i=0;i<n_y_pts;++i) {
      pt[1]=ymin+(i+1)*y_stepsize;
      for(unsigned long long j=0;j<n_x_pts;++j) {
	pt[0]=xmin+(j+1)*x_stepsize;
	(*rf_nonhomog)(pt,rv);
	p_matrix.set(counter,n_x_pts*n_y_pts,-x_stepsize*x_stepsize*rv);
	++counter;
      }
    }
    //bottom bc
    T tmp;
    pt[1]=ymin;
    for(unsigned long long i=0;i<n_x_pts;++i) {
      pt[0]=xmin+i*x_stepsize;
      (*rf_bc)(pt,tmp);
      tmp=tmp+p_matrix.get(i,n_x_pts*n_y_pts);
      p_matrix.set(i,n_x_pts*n_y_pts,tmp);
    }
    //top bc
    pt[1]=ymax;
    for(unsigned long long i=0;i<n_x_pts;++i) {
      pt[0]=xmin+i*x_stepsize;
      (*rf_bc)(pt,tmp);
      tmp=tmp+p_matrix.get((n_y_pts-1)*n_x_pts+i,n_x_pts*n_y_pts);
      p_matrix.set((n_y_pts-1)*n_x_pts+i,n_x_pts*n_y_pts,tmp);
    }
    //left bc
    pt[0]=xmin;
    for(unsigned long long j=0;j<n_y_pts;++j) {
      pt[1]=ymin+j*y_stepsize;
      (*rf_bc)(pt,tmp);
      tmp=tmp+p_matrix.get(j*n_x_pts,n_x_pts*n_y_pts);
      p_matrix.set(j*n_x_pts,n_x_pts*n_y_pts,tmp);
    }
    //right bc
    pt[0]=xmax;
    for(unsigned long long j=0;j<n_y_pts;++j) {
      pt[1]=ymin+j*y_stepsize;
      (*rf_bc)(pt,tmp);
      tmp=tmp+p_matrix.get(j*n_x_pts+n_x_pts-1,n_x_pts*n_y_pts);
      p_matrix.set(j*n_x_pts+n_x_pts-1,n_x_pts*n_y_pts,tmp);
    }
    return 0;
  }

  template<typename T>
  void PoissonBasicMethod2d<T>::printMatrix() {
    std::cout << p_matrix << std::endl;
  }

  template<typename T>
  void PoissonBasicMethod2d<T>::solve(Point<T>& soltn) {
    p_matrix.rowReduce();
    p_matrix.GaussBackSolve(soltn_without_bc);
    
    T tmp;
    Point<T> pt(2);
    // bottom
    pt[1]=ymin;
    for(unsigned long long i=0;i<=n_x_pts+1;++i) {
      pt[0]=xmin+i*x_stepsize;
      (*rf_bc)(pt,tmp);
      soltn[i]=tmp;
    }
    // top
    pt[1]=ymax;
    for(unsigned long long i=0;i<=n_x_pts+1;++i) {
      pt[0]=xmin+i*x_stepsize;
      (*rf_bc)(pt,tmp);
      soltn[(n_x_pts+2)*(n_y_pts+1)+i]=tmp;
    }
    // left
    pt[0]=xmin;
    for(unsigned long long i=1;i<=n_y_pts;++i) {
      pt[1]=ymin+i*y_stepsize;
      (*rf_bc)(pt,tmp);
      soltn[i*(n_x_pts+2)]=tmp;
    }
    // left
    pt[0]=xmax;
    for(unsigned long long i=1;i<=n_y_pts;++i) {
      pt[1]=ymin+i*y_stepsize;
      (*rf_bc)(pt,tmp);
      soltn[i*(n_x_pts+2)+n_x_pts+1]=tmp;
    }
    
    // middle
    for(unsigned long long i=1;i<=n_x_pts;++i) {
      for(unsigned long long j=1;j<=n_y_pts;++j) {
	soltn[j*(n_x_pts+2)+i]=soltn_without_bc[(j-1)*(n_x_pts)+i-1];
      }
    }
  }






















}
