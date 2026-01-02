#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


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
Type minus_one_to_one(Type x)
{
  return Type(2) * invlogit(x) - Type(1);
}

 // dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres;
  logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  if(give_log)return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs_S);    // observed  Spawner
  DATA_VECTOR(obs_logRS);   // observed log recruitment
  DATA_UPDATE(obs_logRS);
  DATA_INTEGER(priors_flag); //flag indicating wether or not priors should be used
  DATA_INTEGER(stan_flag); //flag indicating wether or not 

  DATA_SCALAR(sig_p_sd); //sd for sigma prior
  DATA_SCALAR(logb_p_sd); //sd for logb prior
  DATA_SCALAR(logb_p_mean); //mean for logb prior

  //lfo quantities
  DATA_SCALAR(y_oos); //log(recruits per spawner) next year
  DATA_SCALAR(x_oos); //spawners in next year


  PARAMETER(logalpha);
  PARAMETER(logbeta);
  PARAMETER(logsigobs);
  PARAMETER(rho);
  
  
  int timeSteps=obs_logRS.size();

  Type rhoo = minus_one_to_one(rho);

  
  Type beta = exp(logbeta);
  Type sigobs = exp(logsigobs);
  Type Smax  = Type(1.0)/beta;
  

  Type sigAR  = sigobs*sqrt(1-pow(rhoo,2));

  
  //priors - based on evaluation done with the prior predictive check
  //Type ans = Type(0);
  Type nll= Type(0);
  Type pnll = Type(0.0);

  if(priors_flag == 1){
    
    pnll -=dnorm(logalpha,Type(1.5),Type(2.5),true);
    
    pnll -= dnorm(logbeta,logb_p_mean,logb_p_sd,true);
    
    pnll -= dnorm(sigobs,Type(0.0),sig_p_sd,true) - log(pnorm(Type(0.0), Type(0.0),sig_p_sd));
    if(stan_flag) pnll -= logsigobs;
  }
  
  vector<Type> pred_logRS(timeSteps), pred_logR(timeSteps), residuals(timeSteps) ;
 
  pred_logRS(0) = logalpha - beta * obs_S(0) ;
  pred_logR(0) = pred_logRS(0) + log(obs_S(0));
  residuals(0) = obs_logRS(0) - pred_logRS(0);
    
  nll+= -dnorm(obs_logRS(0),pred_logRS(0),sigobs,true);

  for(int i=1;i<timeSteps;i++){
    if(!isNA(obs_logRS(i))){
      
      pred_logRS(i) = logalpha - beta * obs_S(i) ;
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));
      residuals(i) = obs_logRS(i) - pred_logRS(i);     
      nll+=-dnorm(obs_logRS(i),pred_logRS(i) + residuals(i-1) * rhoo ,sigAR,true);      
    } 
  }
  
  

  Type umsy = (Type(1) - LambertW(exp(1-logalpha)));
  Type Smsy = (Type(1) - LambertW(exp(1-logalpha))) / beta;
  
  Type ans = nll + pnll;

  Type pred_oos = logalpha - beta * x_oos+ residuals(timeSteps-1) * rhoo;
  Type log_lik_oos = dnorm(y_oos,pred_oos,sigAR,true);

  REPORT(logalpha)
  REPORT(beta)
  REPORT(rhoo)
  REPORT(pred_logRS)
  REPORT(residuals)
  REPORT(sigobs)
  REPORT(sigAR)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  REPORT(nll);
  REPORT(pnll);  
  REPORT(log_lik_oos);
 
  ADREPORT(logalpha);
  ADREPORT(beta);
  REPORT(rhoo);
  ADREPORT(sigobs);
  ADREPORT(sigAR);
  ADREPORT(umsy);
  ADREPORT(Smsy);
  

  return ans;
}

