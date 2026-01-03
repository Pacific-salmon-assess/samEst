#include <TMB.hpp>

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
  


template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

 // dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres;
  logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  if(give_log)return logres; else return exp(logres);
}


template <class Type>
Type dstudent(Type x, Type mean, Type sigma, Type df, int give_log = 0) {
  // from metRology::dt.scaled()
  // dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - log(sd)
  Type logres = dt((x - mean) / sigma, df, true) - log(sigma);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  DATA_IVECTOR(options_z);
  DATA_VECTOR(obs_S);    // observed  Spawner
  DATA_VECTOR(obs_logRS);   // observed recruitment
  DATA_UPDATE(obs_logRS);   // line to calculate EDF
  
  DATA_INTEGER(priors_flag); //flag indicating wether or not priors should be used
  DATA_INTEGER(stan_flag); //flag indicating wether or not tmbstan is used 
  DATA_SCALAR(sig_p_sd); //sd for sigma prior
  DATA_SCALAR(siga_p_sd); //sd for siga prior
  DATA_SCALAR(sigb_p_sd); //sd for sigb prior


  //DATA_SCALAR(logb_p_mean); //mean for logb prior
  //DATA_SCALAR(logb_p_sd); //sd for logb prior
  DATA_SCALAR(logsmax_p_sd); //sd for logsmax prior
  DATA_SCALAR(logsmax_p_mean); //mean for logsmax prior

  //fixed parameters
  PARAMETER(logalpha); //initial alpha value
  PARAMETER(logSmax); //log Smax from ricker curve
  PARAMETER(logsigobs);
  PARAMETER(logsiga);
  PARAMETER(logsigb);

  //random effects
 
  PARAMETER_VECTOR(epslogalpha_t);
  PARAMETER_VECTOR(epslogsmax_t);
    
  int timeSteps=obs_logRS.size();

 
  Type sigobs = exp(logsigobs);
  Type siga = exp(logsiga);
  Type sigb = exp(logsigb);


  vector<Type> pred_logR(timeSteps), pred_logRS(timeSteps),umsy(timeSteps), Smsy(timeSteps) ;
  vector<Type> residuals(timeSteps),Srep(timeSteps), ll(timeSteps), logalpha_t(timeSteps); 
  vector<Type>  beta_t(timeSteps), Smax_t(timeSteps); 
 
  //priors on parameters
 
  //Type ans= Type(0);
  Type nll = Type(0.0);
  Type renll = Type(0.0);
  Type pnll = Type(0.0);
  
  if(priors_flag == 1){
    
    pnll -=dnorm(logalpha,Type(1.5),Type(2.5),true);
    pnll -= dnorm(logSmax,logsmax_p_mean,logsmax_p_sd,true);
       
    pnll  -= dnorm(sigobs,Type(0.0),sig_p_sd,true) - log(pnorm(Type(0.0), Type(0.0),sig_p_sd));
    if(stan_flag){
      pnll -= logsigobs;
    } 
    if( options_z(0)==1 ){
      pnll  -= dnorm(siga,Type(0.0),siga_p_sd,true) - log(pnorm(Type(0.0), Type(0.0),siga_p_sd));
      if(stan_flag){
        pnll -= logsiga;
      } 
    }
    if( options_z(1)==1 ){
      pnll  -= dnorm(sigb,Type(0.0),sigb_p_sd,true) - log(pnorm(Type(0.0), Type(0.0),sigb_p_sd));    
      if(stan_flag){
        pnll -= logsigb;
      } 
    }
  }

  for( int t=0; t<obs_logRS.size(); t++ ){
    if( options_z(0)==1 ){
      if(t==0) renll -= dnorm( epslogalpha_t(t), Type(0.0), siga, true );
      if(t>=1) renll -= dnorm( epslogalpha_t(t), epslogalpha_t(t-1), siga, true );
    }
    if( options_z(1)==1 ){
      if(t==0) renll -= dnorm( epslogsmax_t(t), Type(0.0), sigb, true );
      if(t>=1) renll -= dnorm( epslogsmax_t(t), epslogsmax_t(t-1), sigb, true );
    }
  }



  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logRS(i))){
      logalpha_t(i) = logalpha + epslogalpha_t(i);
      Smax_t(i) = exp(logSmax+epslogsmax_t(i));
      beta_t(i) = Type(1.0)/Smax_t(i);
      pred_logRS(i) = logalpha_t(i) - beta_t(i) * obs_S(i) ;
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));
      
      Srep(i) = logalpha_t(i)/beta_t(i);

      Smsy(i) = (Type(1) - LambertW(exp(1-logalpha_t(i))) ) / beta_t(i);
      umsy(i) = (Type(1) - LambertW(exp(1-logalpha_t(i))) ); 

      residuals(i) = obs_logRS(i) - pred_logRS(i);
      ll(i) = dnorm(obs_logRS(i),pred_logRS(i),sigobs,true);
      nll+=-ll(i);
    }
  }
  
  Type Smax = exp(logSmax);
  Type beta = Type(1.0)/Smax;

  Type ans= nll + renll + pnll;

  REPORT(logalpha)
  REPORT(logSmax)
  REPORT(beta)
  REPORT(Smax)
  REPORT(logalpha_t)
  REPORT(beta_t)
  REPORT(sigobs)
  REPORT(siga)
  REPORT(sigb)
  REPORT(pred_logRS)
  REPORT(pred_logR)
  REPORT(residuals)
  REPORT(Smax_t)
  REPORT(umsy)
  REPORT(Smsy)
  REPORT(Srep)
  REPORT(nll);
  REPORT(ll);
  REPORT(renll);
  REPORT(pnll); 
  REPORT(epslogsmax_t); 
  REPORT(epslogalpha_t);
   
  return ans;
}

