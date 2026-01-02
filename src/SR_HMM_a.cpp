#include <TMB.hpp>
#include <iostream>


// Code taken from Brooke
// Set up Lambert's W function to use to calculate SMSY
// Code taken from https://kaskr.github.io/adcomp/lambert_8cpp_source.html
// Step 1: Code up a plain C version
// Double version of Lambert W function
double LambertW(double x) {
  double logx = log(x);
  double y = (logx > 0 ? logx : 0);
  int niter = 100, i=0;
  for (; i < niter; i++) {
    if ( fabs( logx - log(y) - y) < 1e-9) break;
    y -= (y - exp(logx - y)) / (1 + y);
  }
  if (i == niter) Rf_warning("W: failed convergence");
  return y;
}

TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  LambertW
  ,
  // OUTPUT_DIM
  1,
  // ATOMIC_DOUBLE
  ty[0] = LambertW(tx[0]); // Call the 'double' version
,
// ATOMIC_REVERSE
Type W  = ty[0];                    // Function value from forward pass
Type DW = 1. / (exp(W) * (1. + W)); // Derivative
px[0] = DW * py[0];                 // Reverse mode chain rule
)
  
  // Scalar version
  template<class Type>
  Type LambertW(Type x){
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return LambertW(tx)[0];
  }
  



template <class Type>
Type ddirichlet(vector<Type> x, vector<Type> a, int do_log)
{
  Type phi = a.sum();
  int n = x.size();
  Type ll = lgamma(phi);
  for(int i = 0; i < n; i++) ll +=  -lgamma(a(i)) + (a(i) - 1.0) * log(x(i));
  if(do_log == 1) return ll;
  else return exp(ll);
}

template<class Type>
vector<Type> segment_1(vector<Type> yt, vector<Type> st, matrix<Type> qij,vector<Type> pi1,vector<Type> logalpha,Type beta,Type sigma,int t){
  
  int k_regime = logalpha.size();
  Type small = pow(10,-300);
  vector<Type> sr = log(pi1);
  
  for(int j = 0;j < k_regime;++j){
    Type f_now = logalpha(j) - beta*st(0);
    sr(j) += dnorm(yt(0), f_now, sigma,true);
  }
  
  for(int i = 1;i <= t;++i){
    vector<Type> sr_new = sr;    
    
    for(int j = 0;j < k_regime;++j){
      sr_new(j) = sr(0) +qij(0,j);
    
      for(int jj = 1;jj < k_regime;++jj){
        Type temp = sr(jj) +qij(jj,j);
        sr_new(j) = logspace_add(sr_new(j),temp);
      }
    }
    sr = sr_new;
 
    for(int j = 0;j < k_regime;++j){
      Type f_now = logalpha(j) - beta*st(i);
      sr(j) += dnorm(yt(i), f_now, sigma,true);
    }
  }
  return sr;
}


template<class Type>
Type segment_2(vector<Type> yt, vector<Type> st, matrix<Type> qij,vector<Type>
pi1,vector<Type> logalpha,Type beta,Type sigma,int rt,int t){
  
  int k_regime = logalpha.size();
  int n = yt.size();
  vector<Type> sr = qij.row(rt);
  
  for(int j = 0;j < k_regime;++j){
    Type f_now = logalpha(j) - beta*st(t+1);
    sr(j) += dnorm(yt(t+1), f_now, sigma,true);
  }
 
  for(int i = t+2;i < n;++i){
    vector<Type> sr_new = sr;
    for(int j = 0;j < k_regime;++j){
      sr_new(j) = sr(0) +qij(0,j);
      for(int jj = 1;jj < k_regime;++jj){
        Type temp = sr(jj) +qij(jj,j);
        sr_new(j) = logspace_add(sr_new(j),temp);
      }
    }
    sr = sr_new;
    for(int j = 0;j < k_regime;++j){
      Type f_now = logalpha(j) - beta*st(i);
      sr(j) += dnorm(yt(i), f_now, sigma,true);
    }
  }
  Type seg2 = sr(0);
  for(int j = 1;j < k_regime;++j){
    seg2 = logspace_add(seg2,sr(j));
  }
  return seg2;
}

//main //////////////////////////////////////////////
template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(yt);
  DATA_UPDATE(yt);
  DATA_VECTOR(st);
  DATA_SCALAR(logalpha_u); //upper bound for a
  DATA_SCALAR(logalpha_l); //lower bound for a
  DATA_SCALAR(beta_u);  //upper bound for b
  DATA_SCALAR(beta_l); //lower bound for b
  //DATA_SCALAR(sigma_u); //upper bound for sigma

  DATA_VECTOR(initDist);

  DATA_MATRIX(alpha_dirichlet); //prior inputs for dirichlet 
  
  DATA_INTEGER(priors_flag);
  DATA_INTEGER(stan_flag);
  DATA_SCALAR(logb_p_mean); //mean for logb prior
  DATA_SCALAR(logb_p_sd); //sd for logb prior

  PARAMETER_VECTOR(lalpha);
  PARAMETER(lbeta);
  PARAMETER(logsigma);
  //PARAMETER_VECTOR(pi1_tran);
  PARAMETER_MATRIX(qij_tran);

  int k_regime = lalpha.size();
  

  Type beta = (beta_u-beta_l)/(1+exp(-lbeta))+beta_l;// when lbeta is negative infinity, beta=0; when lbeta is positive infinity, beta=beta_u
  Type Smax = Type(1.0)/beta;
  vector<Type> logalpha(k_regime), Smsy(k_regime), umsy(k_regime);
  
  logalpha(0) = (logalpha_u-logalpha_l)/(1+exp(-lalpha(0)))+logalpha_l;
  
  for(int i = 1;i < k_regime;++i){
    logalpha(i) = logalpha(i-1) + (logalpha_u-logalpha(i-1))/(1+exp(-lalpha(i)));

  } // logalpha(1) from logalpha(0) to logalpha_u

  Type sigma = exp(logsigma);
  vector<Type> pi1(k_regime);
 
  for(int i = 0;i < k_regime;++i){
    pi1(i) = initDist(i);
  }
  
  
  Type small = pow(10,-300);
  matrix<Type> qij(k_regime,k_regime);
  for(int i = 0;i < k_regime;++i){
    for(int j = 0;j < k_regime-1;++j){
      qij(i,j) = exp(qij_tran(i,j));
    }
    qij(i,k_regime-1) = 1;
    vector<Type> qij_row = qij.row(i);
    Type row_sum = qij_row.sum();
    for(int j = 0;j < k_regime;++j){
      qij(i,j) = qij(i,j)/row_sum;
      qij(i,j) = log(qij(i,j)+small);
    }
  } 
  int n = yt.size();
  vector<Type> sr = segment_1(yt, st, qij,pi1,logalpha,beta,sigma,n-1);
  Type nll = sr(0);
  
  for(int j = 1;j < k_regime;++j){
    nll = logspace_add(nll,sr(j));
  }
  nll = -nll;

// predict r_t ///////////////////////////////////
  matrix<Type> r_pred(k_regime,n);
  for(int i = 0;i < n-1;++i){
    sr = segment_1(yt, st, qij,pi1,logalpha,beta,sigma,i);
    for(int j = 0;j < k_regime;++j){
      Type tempt = sr(j) + segment_2(yt, st, qij,pi1,logalpha,beta,sigma,j,i) + nll;
      r_pred(j,i) = exp(tempt);
    }
  }
  sr = segment_1(yt, st, qij,pi1,logalpha,beta,sigma,n-1);
  for(int j = 0;j < k_regime;++j){
    Type tempt = sr(j) + nll;
    r_pred(j,n-1) = exp(tempt);

    Smsy(j) = (1 - LambertW(exp(1-logalpha(j))) ) / beta;
    umsy(j) = (1 - LambertW(exp(1-logalpha(j))) ); 
  }
 
  qij = exp(qij.array());

  //priors
  Type pnll = Type(0.0);
  if(priors_flag == 1){
    
 
    Type logbeta = log(beta);
    //pnll -= dnorm(logbeta,Type(-12.0),Type(3.0),true);
    pnll -= dnorm(logbeta,logb_p_mean,logb_p_sd,true);
    

    //pnll -= dgamma(sigma,Type(2.0),Type(1.0)/Type(3.0),true);
    pnll -= dnorm(sigma,Type(0.0),Type(1.0),true) - log(pnorm(Type(0.0), Type(0.0),Type(1.0)));
    if(stan_flag) pnll -= logsigma;
   
    for(int j = 0;j < k_regime;++j){
      
      pnll -=dnorm(logalpha(j),Type(1.5),Type(2.5),true);
     
      vector<Type> qijtmp = qij.row(j);
      vector<Type> alpha_dirichlettmp = alpha_dirichlet.row(j);
      
      pnll -= ddirichlet(qijtmp,alpha_dirichlettmp,true);   
    }
    
  }

  Type ans= nll + pnll;

REPORT(beta);
REPORT(logalpha);
REPORT(sigma);
REPORT(qij);
REPORT(r_pred); 
REPORT(umsy);
REPORT(Smsy);  
REPORT(Smax);  
REPORT(nll);
REPORT(pnll);   


ADREPORT(logalpha);
ADREPORT(beta);
ADREPORT(sigma);
ADREPORT(qij);
ADREPORT(umsy);
ADREPORT(Smsy);
ADREPORT(Smax); 

return ans;

}
