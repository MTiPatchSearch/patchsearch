// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <ctime>
using namespace Rcpp;
const double EPS =  1.e-8;
const double MAXITER =  200;

double *random_vect(int n) {
  int i;
  double *v=(double *) calloc(n, sizeof(double));
  for (i=0; i<n; i++)
    v[i]= (double)std::rand()/(RAND_MAX+1.0);
  return(v);
}

double inner(double *x, double *y, int n) {
  int i;
  double sum;

  for (sum=0, i=0; i<n; sum+=x[i]*y[i],i++);
  return sum;

}

void product(double *result, double *A, double *x, int n) {
  int i, j;
  double sum;

  for (i=0; i<n; i++) {
    sum=0;
    for (j=0; j<n; j++)
      sum+=A[i+n*j]*x[j];
    result[i]=sum;
  }
}


int power(double *a, int n, int maxiter, double eps, double *v, double *w0, int start) {
  int niter,i;
  double *y;
  double sum, l, normy, d, w[4];

  y=w0;
  if (start) y=random_vect(n);
  niter=0;
  do {
    normy=sqrt(inner(y,y,n));
    for (i=0; i<n; i++) w[i]=y[i]/normy;
    product(y, a, w, n);
    l=inner(w,y,n);
    niter++;
    for (sum=0,i=0; i<n; i++) {
      d=y[i]-l*w[i];
      sum+=d*d;
    }
    d=sqrt(sum);
  }
  while (d>eps*fabs(l) && niter<maxiter);
  if (start) free(y);
  *v=l;
  for (i=0; i<4; i++) w0[i]=w[i];
 return niter;
}


double best_shift(double *a, int n) {
  double m, M, s;
  double t, sum;
  int i, j;
  t=a[0];
  for (i=1; i<n; i++) t=std::max(t, a[i+n*i]);
  M=t;
  t=a[0];
  for (i=0; i<n; i++) {
    for (sum=0,j=0; j<n; j++)
      if (j!=i) sum+=fabs(a[i+n*j]);
    t=std::min(t, a[i+n*i]-sum);
  }
  m=t;
  s=-0.5*(M+m);
  for (i=0; i<n; i++)
    a[i+n*i]=a[i+n*i]+s;
  return s;
}


int shift_power(double *a, int n, int maxiter, double eps, double *v, double *w, int start) {
  double sh;
  int niter;
  sh=best_shift(a, n);

  niter=power(a, n, maxiter, eps, v, w, start);
  *v=*v-sh;
  return niter;
}

/*
 * calculs des coefficients d'Euler
 *
 */
void euler(double *X,  double *Y, int n, double *v, double *w, int start) {
  double  K[3][3], A[16];
  double  sum;
  int i,j,k;


  for (i=0; i<3;i++)
    for (j=0; j<3; j++) {
      sum=0;
      for (k=0; k<n;k++)
	sum+=Y[i+3*k]*X[j+3*k];
      K[i][j]=sum;
    }

  A[0+4*0]=K[0][0]+K[1][1]+K[2][2];
  A[1+4*1]=K[0][0]-K[1][1]-K[2][2];
  A[2+4*2]=-K[0][0]+K[1][1]-K[2][2];
  A[3+4*3]=-K[0][0]-K[1][1]+K[2][2];
  A[0+4*1]=A[1+4*0]=K[1][2]-K[2][1];
  A[0+4*2]=A[2+4*0]=K[2][0]-K[0][2];
  A[0+4*3]=A[3+4*0]=K[0][1]-K[1][0];
  A[1+4*2]=A[2+4*1]=K[0][1]+K[1][0];
  A[1+4*3]=A[3+4*1]=K[0][2]+K[2][0];
  A[2+4*3]=A[3+4*2]=K[1][2]+K[2][1];
  shift_power(A, 4, MAXITER , EPS, v, w, start);
}

/*
 * rmsd  entre X:3xN et Y:3xN stockés dans des vecteurs
 * par colonnes mises bout à bout.
 * effet de bord: centrage de X et Y
 */

//' Compute RMSD between matrices of coordinates X and Y
//'
//' @param X matrix Nx3
//' @param Y matrix Nx3
//' @export
// [[Rcpp::export]]
double rmsd(NumericMatrix X,  NumericMatrix Y) {
  double Xm[3], Ym[3];
  double x, y, sum, v, r;
  int i,j,n;

  double w[4];
  for (i=0; i<4; i++) w[i]=0;

  n=X.nrow();
  for (i=0; i<3;i++) {
    sum=0;
    for (j=0; j<n;j++)
      sum+=X(j,i);
    Xm[i]=sum/n;
  }

  for (i=0; i<3;i++) {
    sum=0;
    for (j=0; j<n;j++)
      sum+=Y(j,i);
    Ym[i]=sum/n;
  }
  double XX[3*n],YY[3*n];
  sum=0;
  for (i=0; i<3; i++)
    for (j=0; j<n; j++) {
      x=X(j,i)-Xm[i];
      y=Y(j,i)-Ym[i];
      XX[i+3*j]=x;
      YY[i+3*j]=y;
      sum+=x*x+y*y;
    }

  euler(XX, YY, n, &v, w, 1);

  r=sqrt(fabs(sum-2*v)/n);
  return r;
}

//' Compute deviation between matrices of coordinates X index I and Y index J
//'
//' @param X matrix Nx3
//' @param I vector M
//' @param Y matrix Nx3
//' @param J vector M
//' @export
// [[Rcpp::export]]
NumericVector devsub(NumericMatrix X,  IntegerVector I, NumericMatrix Y, IntegerVector J) {
  double Xm[3], Ym[3];
  double x, y, sum, sum2, v;
  int i,j,k,n;

  n=I.size();
  if (n != J.size()) {
    std::cout << "structures X and Y not of same size!" << I.size() <<  J.size() << std::endl;
    //printf("structures X and Y not of same size! %d %d\n", I.size(), J.size());
    return R_NilValue;
  }

  NumericVector result(n);

  double w[4];
  for (i=0; i<4; i++) w[i]=1.0;

  n=I.size();
  for (i=0; i<3;i++) {
    sum=0;
    for (j=0; j<n;j++)
      sum+=X(I(j),i);
    Xm[i]=sum/n;
  }

  for (i=0; i<3;i++) {
    sum=0;
    for (j=0; j<n;j++)
      sum+=Y(J(j),i);
    Ym[i]=sum/n;
  }

  double XX[3*n],YY[3*n];
  sum=0;
  for (i=0; i<3; i++)
    for (j=0; j<n; j++) {
      x=X(I(j),i)-Xm[i];
      y=Y(J(j),i)-Ym[i];
      XX[i+3*j]=x;
      YY[i+3*j]=y;
    }

  euler(XX, YY, n, &v, w, 0);
  double R[9];

    R[0+3*0] = w[0]*w[0]+w[1]*w[1]-w[2]*w[2]-w[3]*w[3];
    R[1+3*1] = w[0]*w[0]-w[1]*w[1]+w[2]*w[2]-w[3]*w[3];
    R[2+3*2] = w[0]*w[0]-w[1]*w[1]-w[2]*w[2]+w[3]*w[3];
    R[0+3*1] = 2*(w[1]*w[2]+w[0]*w[3]);
    R[0+3*2] = 2*(w[1]*w[3]-w[0]*w[2]);
    R[1+3*0] = 2*(w[1]*w[2]-w[0]*w[3]);
    R[1+3*2] = 2*(w[2]*w[3]+w[0]*w[1]);
    R[2+3*0] = 2*(w[1]*w[3]+w[0]*w[2]);
    R[2+3*1] = 2*(w[2]*w[3]-w[0]*w[1]);

 for (k=0; k<I.size(); k++) {
    sum2=0;
    for (i=0; i<3; i++) {
      sum=0;
      for (j=0; j<3; j++) {
	sum+=R[i+3*j]*XX[j+3*k];
      }
      sum2+=(sum-YY[i+3*k])*(sum-YY[i+3*k]);
    }
    result(k)=sqrt(sum2);
  }

  return result;
}

double distloc(NumericMatrix X, int ind1, int ind2) {
  int i;
  double x,y,sum=0;
  for (i=0;i<3;i++) {
    x=X(ind1-1,i);
    y=X(ind2-1,i);
    sum+=(x-y)*(x-y);
  }
  return sqrt(sum);
}

//' Get set of neighbors that could be added to the clique C
//'
//' @param C
//' @param V
//' @param X
//' @param Y
//' @param thresh
//' @param nbnei
//' @export
// [[Rcpp::export]]
IntegerVector getNeighborsc(IntegerVector C, IntegerMatrix V, NumericMatrix X, NumericMatrix Y, double thresh, int nbnei) {
  int i,j, inC,count=0;
  int nC=C.size();
  IntegerVector K;
  IntegerVector I(nC), J(nC);
  int N=X.nrow();

  for (int i=0; i < nC; i++) {
    I(i)=(C(i)-1)%N+1;
    J(i)=(int)(C(i)-1)/N+1;
  }

  double deltadist, d1, d2;
  for (j=0; j<V.nrow(); j++) {
    count=0;
    inC=0;
    for (i=0; i<nC; i++) {
      if (J(i)!=V(j,0) || I(i)!=V(j,1)) {
	d1=distloc(Y,J(i),V(j,0));
	d2=distloc(X,I(i),V(j,1));
	deltadist=fabs(d1-d2);
	//deltadist=fabs(d1-d2)/std::min(d1,d2);
	if (deltadist<=thresh) count++;
      }
      else inC=1;
    }
    if (count>=nbnei && !inC) K.push_back((V(j,0)-1)*N+V(j,1));
  }
  return K;
}
//' enrich clique with neighbors
//'
//' @param C
//' @param K
//' @param X
//' @param Y
//' @param thresh
//' @export
// [[Rcpp::export]]
IntegerVector enrichloopc(IntegerVector C, IntegerVector K, NumericMatrix X, NumericMatrix Y, double thresh) {

  int i1,j1,ik,jk;
  IntegerVector Ind, Jnd;
  double d;
  int N=X.nrow(), last;
  int nK=K.size();
  int nC=C.size();
  int improved=1, moved=0;;
  NumericVector dev(nC+nK);

  IntegerVector I(nC), J(nC);
  for (int i=0; i < nC; i++) {
    I(i)=(C(i)-1)%N;
    J(i)=(int)(C(i)-1)/N;
  }
  last=nC;
  while (improved) {
    //printf("improvement loop...\n");
    improved=0;
      for(int k=0; k < nK; k++) {
	moved=0;
	if (K(k)<0) continue;
	ik=(K(k)-1)%N;
	jk=(int)(K(k)-1)/N;
	I.push_back(ik);
	J.push_back(jk);
	dev=devsub(X,I,Y,J);
	d=dev(last);
	I.erase(last);
	J.erase(last);
	Ind = match(IntegerVector::create(ik),I);
	Jnd = match(IntegerVector::create(jk),J);
	if (Ind(0)==NA_INTEGER && Jnd(0)==NA_INTEGER) {
	  if (d<=thresh) {
	    I.push_back(ik);
	    J.push_back(jk);
	    moved=1;
	  }
	}
	else if (Ind(0)!=NA_INTEGER && Jnd(0)==NA_INTEGER) {
	    i1=Ind(0)-1;
	    if (d<dev(i1)) {
	      J(i1)=jk;
	      moved=1;
	    }
	}
	else if (Ind(0)==NA_INTEGER && Jnd(0)!=NA_INTEGER) {
	    j1=Jnd(0)-1;
	    if (d<dev(j1)) {
	      I(j1)=ik;
	      moved=1;
	    }
	}
	else if (Ind(0)!=NA_INTEGER && Jnd(0)!=NA_INTEGER){
	    i1=Ind(0)-1;
	    j1=Jnd(0)-1;
	    if (d<dev(i1) && d<dev(j1)) {
	      I(i1)=ik; J(i1)=jk;
	      I.erase(j1);J.erase(j1);
	      //I.push_back(ik);J.push_back(jk);
	      moved=1;
	    }
	}
	last=I.size();
	if (moved) {K(k)=-1;improved=1;} // sort of remove k from K
      }
  }
  if (I.size()==0) return R_NilValue;
  IntegerVector Cplus(I.size());
  for (int i=0; i<I.size(); i++) {
    Cplus(i)=J(i)*N+I(i)+1;
  }
  return Cplus;
}

IntegerVector enrichc(IntegerVector C, IntegerVector K, NumericMatrix X, NumericMatrix Y, double thresh) {

  int i1,j1,ik,jk;
  IntegerVector Ind, Jnd;
  double d;
  int N=X.nrow(), last;
  int nK=K.size();
  int nC=C.size();
  NumericVector dev(nC+nK);

  IntegerVector I(nC), J(nC);
  for (int i=0; i < nC; i++) {
    I(i)=(C(i)-1)%N;
    J(i)=(int)(C(i)-1)/N;
  }
  last=nC;
  for(int k=0; k < nK; k++) {
    ik=(K(k)-1)%N;
    jk=(int)(K(k)-1)/N;
    I.push_back(ik);
    J.push_back(jk);
    dev=devsub(X,I,Y,J);
    d=dev(last);
    I.erase(last);
    J.erase(last);
    Ind = match(IntegerVector::create(ik),I);
    Jnd = match(IntegerVector::create(jk),J);
    if (Ind(0)==NA_INTEGER && Jnd(0)==NA_INTEGER) {
      if (sqrt(d)<=thresh) {
	I.push_back(ik);
	J.push_back(jk);
      }
    }
    else if (Ind(0)!=NA_INTEGER && Jnd(0)==NA_INTEGER) {
	i1=Ind(0)-1;
	if (d<dev(i1)) {
	  J(i1)=jk;
	}
    }
    else if (Ind(0)==NA_INTEGER && Jnd(0)!=NA_INTEGER) {
	j1=Jnd(0)-1;
	if (d<dev(j1)) {
	  I(j1)=ik;
	}
    }
    else if (Ind(0)!=NA_INTEGER && Jnd(0)!=NA_INTEGER){
	i1=Ind(0)-1;
	j1=Jnd(0)-1;
	if (d<dev(i1) && d<dev(j1)) {
	  I(i1)=ik; J(i1)=jk;
	  I.erase(j1);J.erase(j1);
	  //I.push_back(ik);J.push_back(jk);
	}
    }
    last=I.size();
  }

  if (I.size()==0) return R_NilValue;
  IntegerVector Cplus(I.size());
  for (int i=0; i<I.size(); i++) {
    Cplus(i)=J(i)*N+I(i)+1;
  }
  return Cplus;
}

// (V(v,0),V(v,1)): link number v (query, target)
//' Constructs vertex table V
//'
//' @param XProp
//' @param Xseq
//' @param YProp
//' @param Yseq
//' @param S
//' @param sthresh
//' @export
// [[Rcpp::export]]
IntegerMatrix vertexc(CharacterVector XProp, IntegerVector Xseq, CharacterVector YProp, IntegerVector Yseq, NumericMatrix S, double sthresh) {
  int i,ip,v=0;
  int N=XProp.size();
  int M=YProp.size();
  IntegerMatrix V(N*M,2);
  for(i=0; i<M; i++)
    for (ip=0; ip<N; ip++)
      if (strcmp(XProp(ip),YProp(i))==0 && S(Xseq(ip),Yseq(i))>=sthresh) { /* S 23x23 */
	V(v,0)=i+1;
	V(v,1)=ip+1;
	v++;
      }
  if (v==0) {
    fprintf(stderr,"no vertices\n");
    return V;
  }
  return V(Range(0,v-1),_);
}


//' Constructs product graph
//'
//' @param X matrix Nx3
//' @param Y matrix Nx3
//' @param V matrix Nvx2
//' @export
// [[Rcpp::export]]
IntegerMatrix buildGraphc(NumericMatrix X, NumericMatrix Y, IntegerMatrix V, IntegerVector localdegree, double thresh, double localthresh) {
  int i,j, nv, e=0;
  double d,d1,d2;
  std::vector<int> Etmp;

  nv=V.nrow();
  printf("buildGraphc::thresh=%lf %lf %d\n", thresh, localthresh, nv);
  for(i=0; i<nv-1; i++) {
    for (j=i+1; j<nv; j++) {
      //printf("%d %d %d %d\n",V(i,0),V(j,0),V(i,1),V(j,1));
      d1=distloc(Y,V(i,0),V(j,0));
      d2=distloc(X,V(i,1),V(j,1));
      d=fabs(d1-d2)/std::min(d1,d2);
      //printf("d=%lf\n",d);
      if (d<=thresh && d1>EPS && d2>EPS) {
      //if (d<=thresh*std::min(d1,d2) && d1>EPS && d2>EPS) {
	      Etmp.push_back(V(i,1));
	      Etmp.push_back(V(j,1));
	      Etmp.push_back(V(i,0));
	      Etmp.push_back(V(j,0));
	      e+=1;
	      if (d1<=localthresh && d2<=localthresh) {
		localdegree(i)+=1;
		localdegree(j)+=1;
	      }
	}
    }
  }
  int n=Etmp.size();
  printf("n=%d,e=%d\n",n,e);
  IntegerMatrix E(e,4);
  for (i=0,j=0; i<n; i+=4,j++) {
    E(j,0)=Etmp[i];
    E(j,1)=Etmp[i+1];
    E(j,2)=Etmp[i+2];
    E(j,3)=Etmp[i+3];
  }
  Etmp.clear();
return E;
}
