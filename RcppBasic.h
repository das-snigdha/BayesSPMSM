using namespace Rcpp;

NumericVector c(double a, double b){
  NumericVector x = NumericVector::create(a, b) ;
  return x ;
}

NumericVector c(double b, NumericVector a){
  NumericVector x(a.size() + 1) ;
  x[0] = b ;
  for(int i=1; i<=a.size(); i++)
    x[i] = a[i] ;
  return x ;
}

NumericVector c(NumericVector a, double b){
  NumericVector x(a.size() + 1) ;
  for(int i=0; i<a.size(); i++)
    x[i] = a[i] ;
  x[a.size()] = b ;
  return x ;
}

NumericVector c(NumericVector a, NumericVector b){
  NumericVector x(a.size() + b.size()) ;
  for(int i=0; i<a.size(); i++)
    x[i] = a[i] ;
  for(int i=0; i<b.size(); i++)
    x[i + a.size()] = b[i] ;
  return x ;
}

IntegerVector c(int a, int b){
  IntegerVector x = IntegerVector::create(a, b) ;
  return x ;
}

IntegerVector c(IntegerVector a, int b){
  IntegerVector x(a.size() + 1) ;
  for(int i=0; i<a.size(); i++)
    x[i] = a[i] ;
  x[a.size()] = b ;
  return x ;
}

IntegerVector c(IntegerVector a, IntegerVector b){
  IntegerVector x(a.size() + b.size()) ;
  for(int i=0; i<a.size(); i++)
    x[i] = a[i] ;
  for(int i=0; i<b.size(); i++)
    x[i + a.size()] = b[i] ;
  return x ;
}

void asMatrix(NumericVector x, NumericMatrix res, bool byrow=false){
  if(x.size() != res.nrow()*res.ncol())
    Rprintf("Error: Dimensions do not match in asMatrix\n") ;
  int ind=0 ;
  if(byrow==false){
    for(int j=0; j<res.ncol(); j++){
      for(int i=0; i<res.nrow(); i++){
        res(i,j) = x[ind] ;
        ind++ ;
      }
    }
  }
  else{
    for(int i=0; i<res.nrow(); i++){
      for(int j=0; j<res.ncol(); j++){
        res(i,j) = x[ind] ;
        ind++ ;
      }
    }
  }
}

void asVector(NumericMatrix A, NumericVector res, bool byrow=false){
  if(res.size() != A.nrow()*A.ncol())
    Rprintf("Error: Dimensions do not match in asVector\n") ;
  int ind=0 ;
  if(byrow==false){
    for(int j=0; j<A.ncol(); j++){
      for(int i=0; i<A.nrow(); i++){
        res[ind] = A(i,j) ;
        ind++ ;
      }
    }
  }
  else{
    for(int i=0; i<A.nrow(); i++){
      for(int j=0; j<A.ncol(); j++){
        res[ind] = A(i,j) ;
        ind++ ;
      }
    }
  }
}

NumericVector seq(double from, double to, int length){
  NumericVector res(length) ;
  double tmp = (to-from) / (double)(length-1) ;
  res[0] = from ;
  res[length-1] = to ;
  for(int i=1; i<length-1; i++)
    res[i] = res[i-1] + tmp ;
  return(res) ;
}

int sum(IntegerVector x){
  int res = 0 ;
  for(int i=0; i<x.size(); i++)
    res += x[i] ;
  return(res) ;
}

double sum(NumericVector x){
  double res = 0.0 ;
  for(int i=0; i<x.size(); i++)
    res += x[i] ;
  return(res) ;
}

void sum(NumericMatrix A, double b, NumericMatrix res){
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      res(i,j) = A(i,j) + b ;
}

void sum(NumericMatrix A, NumericMatrix B, NumericMatrix res){
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      res(i,j) = A(i,j) + B(i,j) ;
}

NumericMatrix sum(NumericMatrix A, double b){
  
  NumericMatrix res(A.nrow(), A.ncol()) ;
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      res(i,j) = A(i,j) + b ;
  return(res) ;
}

NumericMatrix sum(NumericMatrix A, NumericMatrix B){
  
  NumericMatrix res(A.nrow(), A.ncol()) ;
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      res(i,j) = A(i,j) + B(i,j) ;
  return(res) ;
}

NumericVector colSums(NumericMatrix X){
  int n=X.nrow(), p=X.ncol() ;
  NumericVector res(p) ;
  for(int j=0; j<p; j++)
    for(int i=0; i<n; i++)
      res[j] += X(i,j) ;
  return(res) ;
}

NumericVector rowSums(NumericMatrix X){
  int n=X.nrow(), p=X.ncol() ;
  NumericVector res(n) ;
  for(int i=0; i<n; i++)
    for(int j=0; j<p; j++)
      res[i] += X(i,j) ;
  return(res) ;
}

double innerProduct(NumericVector x, NumericVector y){
  double res=0.0 ;
  for(int i=0; i<x.size(); i++)
    res += x[i]*y[i] ;
  return(res) ;
}

double norm(NumericVector x, int p=2, bool power=false){
  double res=0.0 ;
  for(int i=0; i<x.size(); i++)
    res += pow(x[i], p) ;
  if(power==true)
    return(res) ;
  else
    return(pow(res, 1.0/(double)p)) ;
}

void product(NumericMatrix A, NumericMatrix B, NumericMatrix res){
  for(int i=0; i<A.nrow(); i++)	{
    for(int j=0; j<B.ncol(); j++)	{
      res(i,j) = 0.0 ;
      for(int k=0; k<A.ncol(); k++)
        res(i,j) += A(i,k)*B(k,j) ;
    }
  }
}

void product(NumericMatrix A, double x, NumericMatrix res){
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      res(i,j) = x * A(i,j) ;
}

void product(NumericMatrix A, NumericVector x, NumericVector res){
  for(int i=0; i<A.nrow(); i++)	{
    res[i] = 0.0 ;
    for(int j=0; j<A.ncol(); j++)
      res[i] += A(i,j)*x[j] ;
  }
}

NumericVector product(NumericMatrix A, NumericVector x){
  NumericVector res(x.size()) ;
  for(int i=0; i<A.nrow(); i++)	{
    res[i] = 0.0 ;
    for(int j=0; j<A.ncol(); j++)
      res[i] += A(i,j)*x[j] ;
  }
  return(res) ;
}

NumericMatrix product(NumericMatrix A, NumericMatrix B){
  NumericMatrix res(A.nrow(), B.ncol()) ;
  for(int i=0; i<A.nrow(); i++)	{
    for(int j=0; j<B.ncol(); j++)	{
      res(i,j) = 0.0 ;
      for(int k=0; k<A.ncol(); k++)
        res(i,j) += A(i,k)*B(k,j) ;
    }
  }
  return(res) ;
}

void transpose(NumericMatrix A, NumericMatrix t_A){
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      t_A(j,i) = A(i,j) ;
}

NumericMatrix transpose(NumericMatrix A){
  NumericMatrix t_A(A.nrow(), A.ncol()) ;
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      t_A(j,i) = A(i,j) ;
  return(t_A) ;
}

double quadraticForm(NumericVector x, NumericMatrix A) {
  int n=A.nrow() ;
  double res=0.0 ;
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      res += x[i]*x[j]*A(i,j) ;
  return(res) ;
}

void Atx(NumericMatrix A, NumericVector x, NumericVector res){
  for(int i=0; i<A.ncol(); i++)	{
    res[i] = 0.0 ;
    for(int j=0; j<A.nrow(); j++)
      res[i] += A(j,i)*x[j] ;
  }
}

void AtA(NumericMatrix A, NumericMatrix res){
  // res = t(A) %*% A
  for(int i=0; i<A.ncol(); i++){
    for(int j=0; j<A.ncol(); j++){
      res(i,j) = 0.0 ;
      for(int k=0; k<A.nrow(); k++)
        res(i,j) += A(k,i)*A(k,j) ;
    }
  }
}

void AAt(NumericVector A, NumericMatrix res){
  // res = A %*% t(A) 
  for(int i=0; i<A.size(); i++){
    for(int j=0; j<A.size(); j++){
      res(i,j) = A[i]*A[j] ;
    }
  }
}

double trace(NumericMatrix A){
  double res=0.0 ;
  for(int i=0; i<A.nrow(); i++)
    res += A(i,i) ;
  return(res) ;
}

void vec_mem_cpy(NumericVector from, NumericVector to){
  for(int i=0; i<to.size(); i++)
    to[i] = from[i] ;
}

void vec_mem_cpy(NumericVector from, NumericVector to, int n){
  for(int i=0; i<n; i++)
    to[i] = from[i] ;
}

void vec_mem_cpy(IntegerVector from, IntegerVector to){
  for(int i=0; i<to.size(); i++)
    to[i] = from[i] ;
}

void mtx_mem_cpy(NumericMatrix from, NumericMatrix to){
  for(int i=0; i<to.nrow(); i++)
    for(int j=0; j<to.ncol(); j++)
      to(i,j) = from(i,j) ;
}

int sample(IntegerVector x, NumericVector prob, int length, bool normalized=false){
  
  double p=runif(1)[0], cum=0.0, p_sum=0.0 ;
  NumericVector n_prob(length) ;
  
  if(normalized==false)	{
    for(int k=0; k<length; k++)
      p_sum += prob[k] ;
  }else	{
    p_sum = 1.0 ;
  }
  for(int k=0; k<length-1; k++)	{
    n_prob[k] = prob[k]/p_sum ;
    cum += n_prob[k] ;
    if(p < cum)
      return(x[k]) ;
  }
  return(x[length-1]) ;
}

void chol(NumericMatrix A, NumericMatrix res){
  int n = A.nrow() ;
  arma::mat tmp = arma::chol(as<arma::mat>(A)) ;
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      res(i,j) = tmp(i,j) ;
}

NumericMatrix chol(NumericMatrix A){
  int n = A.nrow() ;
  NumericMatrix res(n,n) ;
  arma::mat tmp = arma::chol(as<arma::mat>(A)) ;
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      res(i,j) = tmp(i,j) ;
  return(res) ;
}

void solve(NumericMatrix A, NumericMatrix A_inv){
  arma::mat tmp = as<arma::mat>(A).i() ;
  
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      A_inv(i,j) = tmp(i,j) ;
}

NumericMatrix solve(NumericMatrix A){
  arma::mat res = as<arma::mat>(A).i() ;
  return(wrap(res)) ;
}

NumericMatrix rmvnorm(int n, NumericVector mu, NumericMatrix Sigma){
  arma::mat Y = arma::randn(n,mu.size());
  arma::mat Z = arma::repmat(as<arma::vec>(mu),1,n).t() + Y * arma::chol(as<arma::mat>(Sigma)) ;
  return(wrap(Z)) ;
}

NumericVector rmvnorm(NumericVector mu, NumericMatrix Sigma){
  arma::mat Y = arma::randn(1,mu.size());
  arma::mat Z = arma::repmat(as<arma::vec>(mu),1,1).t() + Y * arma::chol(as<arma::mat>(Sigma)) ;
  return(wrap(Z)) ;
}

void rmvnorm(NumericVector mu, NumericMatrix Sigma, NumericVector res){
  arma::mat Y = arma::randn(1,mu.size());
  arma::mat Z = arma::repmat(as<arma::vec>(mu),1,1).t() + Y * arma::chol(as<arma::mat>(Sigma)) ;
  for(int i=0; i<res.size(); i++)
    res[i] = Z(0,i) ;
}

double dmvnorm(NumericVector x, NumericVector mean, NumericMatrix Sigma, bool logd = false){
  int d = x.size();
  double res ;
  NumericVector z(d) ;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(as<arma::mat>(Sigma))))) ;
  
  res = -0.5*d*log(2.0*M_PI) + arma::sum(log(rooti.diag())) ;
  
  product(wrap(rooti), x-mean, z) ;
  res -= 0.5 * innerProduct(z,z) ;
  
  
  if (logd==false) {
    res=exp(res);
  }
  return(res) ;
}

NumericVector rigamma(int n, double shape, double rate){
  return(1.0/rgamma(n, shape, 1.0/rate)) ;
}

NumericMatrix rDirichlet(int n, NumericVector alpha){
  int K=alpha.size() ;
  double sum_gam ;
  NumericVector gam(K) ;
  NumericMatrix res(n,K) ;
  
  for(int i=0; i<n; i++){
    for(int k=0; k<K; k++)
      gam[k] = rgamma(1, alpha[k], 1.0)[0] ;
    sum_gam = sum(gam) ;
    for(int k=0; k<K; k++)
      res(i,k) = gam[k]/sum_gam ;
  }
  return(res) ;
}

void rWishart(int n, NumericMatrix V, NumericMatrix res){
  int d=V.ncol() ;
  NumericMatrix X=rmvnorm(n, rep(0.0,d), V) ;
  AtA(X, res) ;
}

void rInvWishart(int n, NumericMatrix V, NumericMatrix res){
  NumericMatrix tmp(V.nrow(),V.ncol()) ;
  rWishart(n, V, tmp) ;
  solve(tmp, res) ;
}

NumericMatrix diag_c(int n, double x = 1.0){
  NumericMatrix res(n, n);
  for(int i=0; i<n; i++){
    res(i,i) = x ;
  }
  return(res) ;
}

NumericVector log_vec(NumericVector x){
  int n = x.size() ;
  NumericVector res(n) ;
  
  for(int i = 0; i<n; i++){
    res[i] = log(x[i]) ; 
  }
  
  return(res) ;
}

NumericVector exp_vec(NumericVector x){
  int n = x.size() ;
  NumericVector res(n) ;
  
  for(int i = 0; i<n; i++){
    res[i] = exp(x[i]) ; 
  }
  
  return(res) ;
}