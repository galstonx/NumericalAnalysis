
namespace NumericalAnalysis {

  template<typename T>
  PoissonBasicMethod2d<T>::PoissonBasicMethod2d(const RealFunction<T>* f,const RealFunction<T>* g,T x0,T x1,T y0,T y1,unsigned long long x_steps,unsigned long long y_steps) : p_matrix((x_steps-1)*(y_steps-1),(x_steps-1)*(y_steps-1)+1),soltn_without_bc((x_steps-1)*(y_steps-1)) {
    rf_nonhomog=f;
    rf_bc=g;
    n_x_pts=x_steps-1;
    n_y_pts=y_steps-1;
    xmin=x0;
    xmax=x1;
    ymin=y0;
    ymax=y1;
    x_stepsize=(xmax-xmin)/(n_x_pts+1);
    y_stepsize=(ymax-ymin)/(n_y_pts+1);
    init_matrix();
    init_nonhomog_vector();
  }

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
	rf_nonhomog->eval(pt,rv);
	p_matrix.set(counter,n_x_pts*n_y_pts,-x_stepsize*x_stepsize*rv);
	++counter;
      }
    }
    //bottom bc
    T tmp;
    pt[1]=ymin;
    for(unsigned long long i=0;i<n_x_pts;++i) {
      pt[0]=xmin+i*x_stepsize;
      rf_bc->eval(pt,tmp);
      tmp=tmp+p_matrix.get(i,n_x_pts*n_y_pts);
      p_matrix.set(i,n_x_pts*n_y_pts,tmp);
    }
    //top bc
    pt[1]=ymax;
    for(unsigned long long i=0;i<n_x_pts;++i) {
      pt[0]=xmin+i*x_stepsize;
      rf_bc->eval(pt,tmp);
      tmp=tmp+p_matrix.get((n_y_pts-1)*n_x_pts+i,n_x_pts*n_y_pts);
      p_matrix.set((n_y_pts-1)*n_x_pts+i,n_x_pts*n_y_pts,tmp);
    }
    //left bc
    pt[0]=xmin;
    for(unsigned long long j=0;j<n_y_pts;++j) {
      pt[1]=ymin+j*y_stepsize;
      rf_bc->eval(pt,tmp);
      tmp=tmp+p_matrix.get(j*n_x_pts,n_x_pts*n_y_pts);
      p_matrix.set(j*n_x_pts,n_x_pts*n_y_pts,tmp);
    }
    //right bc
    pt[0]=xmax;
    for(unsigned long long j=0;j<n_y_pts;++j) {
      pt[1]=ymin+j*y_stepsize;
      rf_bc->eval(pt,tmp);
      tmp=tmp+p_matrix.get(j*n_x_pts+n_x_pts-1,n_x_pts*n_y_pts);
      p_matrix.set(j*n_x_pts+n_x_pts-1,n_x_pts*n_y_pts,tmp);
    }
    return 0;
  }

  template<typename T>
  int PoissonBasicMethod2d<T>::printMatrix() {
    std::cout << p_matrix << std::endl;
  }

  template<typename T>
  int PoissonBasicMethod2d<T>::solve(Point<T>& soltn) {
    p_matrix.rowReduce();
    p_matrix.GaussBackSolve(soltn_without_bc);
    
    T tmp;
    Point<T> pt(2);
    // bottom
    pt[1]=ymin;
    for(unsigned long long i=0;i<=n_x_pts+1;++i) {
      pt[0]=xmin+i*x_stepsize;
      rf_bc->eval(pt,tmp);
      soltn[i]=tmp;
    }
    // top
    pt[1]=ymax;
    for(unsigned long long i=0;i<=n_x_pts+1;++i) {
      pt[0]=xmin+i*x_stepsize;
      rf_bc->eval(pt,tmp);
      soltn[(n_x_pts+2)*(n_y_pts+1)+i]=tmp;
    }
    // left
    pt[0]=xmin;
    for(unsigned long long i=1;i<=n_y_pts;++i) {
      pt[1]=ymin+i*y_stepsize;
      rf_bc->eval(pt,tmp);
      soltn[i*(n_x_pts+2)]=tmp;
    }
    // left
    pt[0]=xmax;
    for(unsigned long long i=1;i<=n_y_pts;++i) {
      pt[1]=ymin+i*y_stepsize;
      rf_bc->eval(pt,tmp);
      soltn[i*(n_x_pts+2)+n_x_pts+1]=tmp;
    }
    
    // middle
    for(unsigned long long i=1;i<=n_x_pts;++i) {
      for(unsigned long long j=1;j<=n_y_pts;++j) {
	soltn[j*(n_x_pts+2)+i]=soltn_without_bc[(j-1)*(n_x_pts)+i-1];
      }
    }
    
    return 0;
  }






















}
