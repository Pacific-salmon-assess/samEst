#include <TMB.hpp>
#include <iostream>

template<class Type>
vector<Type> segment_1(vector<Type> yt, vector<Type> st, matrix<Type> qij,vector<Type> pi1,vector<Type> alpha,vector<Type> beta,vector<Type> sigma,int t){
  
  int k_regime = beta.size();
  Type small = pow(10,-300);
  vector<Type> sr = log(pi1 + small);
  
  for(int j = 0;j < k_regime;++j){
    Type f_now = alpha(j) - beta(j)*st(0);
    sr(j) += dnorm(yt(0), f_now, sigma(j),true);
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
      Type f_now = alpha(j) - beta(j)*st(i);
      sr(j) += dnorm(yt(i), f_now, sigma(j),true);
    }
  }
  return sr;
}


template<class Type>
Type segment_2(vector<Type> yt, vector<Type> st, matrix<Type> qij,vector<Type>
pi1,vector<Type> alpha,vector<Type> beta,vector<Type> sigma,int rt,int t){
  
  int k_regime = beta.size();
  int n = yt.size();
  vector<Type> sr = qij.row(rt);
  
  for(int j = 0;j < k_regime;++j){
    Type f_now = alpha(j) - beta(j)*st(t+1);
    sr(j) += dnorm(yt(t+1), f_now, sigma(j),true);
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
      Type f_now = alpha(j) - beta(j)*st(i);
      sr(j) += dnorm(yt(i), f_now, sigma(j),true);
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
  DATA_VECTOR(st);
  DATA_SCALAR(alpha_u); //upper bound for b
  DATA_SCALAR(alpha_l); //lower bound for b
  DATA_SCALAR(beta_u);  //upper bound for a
  DATA_SCALAR(sigma_u); //upper bound for sigma

  PARAMETER_VECTOR(lalpha);
  PARAMETER_VECTOR(lbeta);
  PARAMETER_VECTOR(lsigma);
  PARAMETER_VECTOR(pi1_tran);
  PARAMETER_MATRIX(qij_tran);

  int k_regime = lbeta.size();
  // vector<Type> alpha = alpha_tr;
  // for(int i = 1;i < k_regime;++i){
  // alpha(i) = alpha(i-1) + exp(alpha_tr(i));
  // }

  vector<Type> beta = beta_u/(1+exp(-lbeta));// when lbeta is negative infinity, beta=0; when lbeta is positive infinity, beta=beta_u
  vector<Type> alpha(k_regime);
  
  alpha(0) = (alpha_u-alpha_l)/(1+exp(-lalpha(0)))+alpha_l;
  
  for(int i = 1;i < k_regime;++i){
    alpha(i) = alpha(i-1) + (alpha_u-alpha(i-1))/(1+exp(-lalpha(i)));
  } // alpha(1) from alpha(0) to alpha_u

  vector<Type> sigma = sigma_u/(1+exp(-lsigma));
  vector<Type> pi1(k_regime);
 
  for(int i = 0;i < k_regime-1;++i){
    pi1(i) = exp(pi1_tran(i));
  }
  pi1(k_regime-1) = 1;
  pi1 = pi1/(pi1.sum());
  
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
  vector<Type> sr = segment_1(yt, st, qij,pi1,alpha,beta,sigma,n-1);
  Type nll = sr(0);
  
  for(int j = 1;j < k_regime;++j){
    nll = logspace_add(nll,sr(j));
  }
  nll = -nll;

// predict r_t ///////////////////////////////////
  matrix<Type> r_pred(k_regime,n);
  for(int i = 0;i < n-1;++i){
    sr = segment_1(yt, st, qij,pi1,alpha,beta,sigma,i);
    for(int j = 0;j < k_regime;++j){
      Type tempt = sr(j) + segment_2(yt, st, qij,pi1,alpha,beta,sigma,j,i) + nll;
      r_pred(j,i) = exp(tempt);
    }
  }
  sr = segment_1(yt, st, qij,pi1,alpha,beta,sigma,n-1);
  for(int j = 0;j < k_regime;++j){
    Type tempt = sr(j) + nll;
    r_pred(j,n-1) = exp(tempt);
  }
 
 qij = exp(qij.array());

REPORT(beta);
REPORT(alpha);
REPORT(sigma);
REPORT(pi1);
REPORT(qij);
REPORT(r_pred);     


ADREPORT(alpha);
ADREPORT(beta);
ADREPORT(sigma);
ADREPORT(pi1);
ADREPORT(qij);

return nll;

}
