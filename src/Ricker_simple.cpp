#include <TMB.hpp>

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
  

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs_S);    // observed  Spawner
  DATA_VECTOR(obs_logRS);   // observed log recruitment
  DATA_INTEGER(priors);
  
  PARAMETER(alpha);
  PARAMETER(logbeta);
  PARAMETER(logsigobs);
  
  int timeSteps=obs_logRS.size();

  //priors - based on evaluation done with the prior predictive check
  //Type ans= Type(0);
  Type nll= Type(0.0);
  Type pnll = Type(0.0);
  
  //model
  Type beta = exp(logbeta);
  Type sigobs = exp(logsigobs);
  Type Smax  = Type(1.0)/beta;

  if(priors == 1){
    //ans -=dnorm(alpha,Type(0.0),Type(2.5),true);
    pnll -=dgamma(alpha,Type(3.0),Type(1.0),true);
    pnll -=dnorm(logbeta,Type(-12.0),Type(3.0),true);  
    pnll -= dgamma(sigobs,Type(2.0),Type(1.0)/Type(3.0),true);
  }
  //ans -= dnorm(logsigobs,Type(0.0),Type(2.0),true);
  //ans -= dexp(sigobs,Type(2.0),true);
  //ans -= dt(sigobs,Type(3.0),true);
  
  
  vector<Type> pred_logRS(timeSteps), pred_logR(timeSteps), residuals(timeSteps); 
   
  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logRS(i))){
      pred_logRS(i) = alpha - beta * obs_S(i) ; 
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));
      residuals(i) = obs_logRS(i) - pred_logRS(i);
      nll+=-dnorm(obs_logRS(i),pred_logRS(i),sigobs,true);
      
    }
  
  }
  //Type umsy = Type(.5) * alpha - Type(0.07) * (alpha * alpha);
  //Type Smsy = alpha/beta * (Type(0.5) -Type(0.07) * alpha);
  Type umsy = (Type(1) - LambertW(exp(1-alpha)));
  Type Smsy = (Type(1) - LambertW(exp(1-alpha))) / beta;
  Type ans= nll + pnll;


  SIMULATE {
   vector<Type> R_Proj(timeSteps);
   for(int i=0; i<timeSteps; ++i){
     R_Proj(i) = exp(rnorm(pred_logR(i), sigobs));
   }
   REPORT(R_Proj);
  }

  REPORT(pred_logR)
  REPORT(pred_logRS)
  REPORT(residuals)
  REPORT(alpha)  
  REPORT(beta)
  REPORT(sigobs)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  REPORT(nll);
  REPORT(pnll);  

  ADREPORT(alpha);
  ADREPORT(beta);
  ADREPORT(sigobs);
  ADREPORT(umsy);
  ADREPORT(Smsy);
  
  return ans;
}

