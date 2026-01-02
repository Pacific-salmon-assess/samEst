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
  DATA_UPDATE( obs_logRS );
  DATA_INTEGER(priors_flag); //flag indicating wether or not priors should be used
  DATA_INTEGER(stan_flag); //flag indicating wether or not tmbstan is used 
  DATA_SCALAR(sig_p_sd); //sd for sigma prior
  DATA_SCALAR(logb_p_mean); //mean for logb prior
  DATA_SCALAR(logb_p_sd); //sd for logb prior

  //DATA_SCALAR(sig_p_mean);
  
  //lfo quantities
  DATA_SCALAR(y_oos); //log(recruits per spawner) next year
  DATA_SCALAR(x_oos); //spawners in next year

  
  
  PARAMETER(logalpha);
  PARAMETER(logbeta);
  PARAMETER(logsigobs);
  
  int timeSteps=obs_logRS.size();

  //priors - based on evaluation done with the prior predictive check
  
  Type nll= Type(0.0);
  Type pnll = Type(0.0);
  
  //model
  Type beta = exp(logbeta);
  Type sigobs = exp(logsigobs);
  Type Smax  = Type(1.0)/beta;

  if(priors_flag == 1){
    
    pnll -=dnorm(logalpha,Type(1.5),Type(2.5),true);
    
    pnll -= dnorm(logbeta,logb_p_mean,logb_p_sd,true);
    
    
    pnll -= dnorm(sigobs,Type(0.0),sig_p_sd,true) - log(pnorm(Type(0.0), Type(0.0),sig_p_sd));
    if(stan_flag) pnll -= logsigobs; //Jacobian for half normal prior
  }

  
  vector<Type> pred_logRS(timeSteps), pred_logR(timeSteps), residuals(timeSteps), ll(timeSteps) ; 
   
  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logRS(i))){
      pred_logRS(i) = logalpha - beta * obs_S(i) ; 
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));
      residuals(i) = obs_logRS(i) - pred_logRS(i);
      ll(i) = dnorm(obs_logRS(i),pred_logRS(i),sigobs,true);
      nll+=-ll(i);
      
    }
  
  }
  
  Type umsy = (Type(1) - LambertW(exp(1-logalpha)));
  Type Smsy = (Type(1) - LambertW(exp(1-logalpha))) / beta;
  Type ans= nll + pnll;

  Type pred_oos = logalpha - beta * x_oos;
  Type log_lik_oos = dnorm(y_oos,pred_oos,sigobs,true);


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
  REPORT(logalpha)  
  REPORT(beta)
  REPORT(sigobs)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  REPORT(nll);
  REPORT(ll);
  REPORT(pnll); 
  REPORT(log_lik_oos); 

  ADREPORT(logalpha);
  ADREPORT(beta);
  ADREPORT(sigobs);
  ADREPORT(umsy);
  ADREPORT(Smsy);
  
  return ans;
}