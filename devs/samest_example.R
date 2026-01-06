#install.packages("remotes") 
remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')

library(samEst);library(ggplot2)

data(harck) #load example dataset - harrison lake chinook

#note - data fed to next functions must have spawners (S), recruits (R), broodyear (by), and annual productivity (logRS) with these dataframe names
df=data.frame(S=harck$S,R=harck$R,by=harck$by,logRS=harck$logRS)

#There are both TMB and Bayesian (stan) options for each function - Bayesian requires rstan

#To estimate a traditional (static) Ricker model, with 1-year residual autocorrelation

static_mod=samEst::ricker_TMB(data=df,AC=T,priors_flag=0)
samEst::sr_plot(df=df,mod=static_mod,type='static',form='tmb',title='Harrison chinook')
#plot function to visualize - using the different model/estimation forms
#list of parameters by $
static_mod$logalpha

#priors_flag = no priors on any parameters, can remove this (defaults to 1) to have reasonable priors applied for alpha,smax,sigma, priors are required for bayesian version
static_mod_st=samEst::ricker_stan(data=df,ac=T)
samEst::sr_plot(df=df,mod=static_mod_st,type='static',form='stan',title='Harrison chinook')
#stan form of plot function

#Time-varying productivity model variant
tv_mod_tmb=samEst::ricker_rw_TMB(data=df,tv.par='a')
tv_mod_tmb$logalpha #maximum likelihood estimates of log(alpha) by year
samEst::sr_plot(df=df,mod=tv_mod_tmb,type='rw',par='a',form='tmb',title='Harrison chinook')

tv_mod_st=samEst::ricker_rw_stan(data=df,par='a') 
samEst::sr_plot(df=df,mod=tv_mod_st,type='rw',par='a',form='stan',title='Harrison chinook')

