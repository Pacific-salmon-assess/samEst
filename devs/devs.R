#============================================
#commands used for packages development and testing
#Catarina wor
#August 2022
#============================================


devtools::document()


devtools::load_all()




p_tmbstan<-ricker_TMBstan(data=harck)
pac_tmbstan<-ricker_TMBstan(data=harck, AC=TRUE)
ptva_tmbstan<- ricker_rw_TMBstan(data=harck,tv.par="a")
ptvb_tmbstan<- ricker_rw_TMBstan(data=harck,tv.par="b")
ptvab_tmbstan<- ricker_rw_TMBstan(data=harck,tv.par="both")

phmma_tmbstan<- ricker_rw_TMBstan(data=harck,tv.par="a")
phmmb_tmbstan<- ricker_rw_TMBstan(data=harck,tv.par="b")
phmmab_tmbstan<- ricker_rw_TMBstan(data=harck,tv.par="both")

#render the vignette
devtools::build_rmd("vignettes/samEst.Rmd")
usethis::use_logo("man/figures/samEst_large.png")


#lfo tmb test
#combine simple and AR

p <- ricker_TMB_deprecated(data=harck, priors_flag=1)
p2 <- ricker_TMB(data=harck,Smax_mean=220000,Smax_sd=230000, priors_flag=1)




p$Smax
p2$Smax


p$logalpha   
p2$logalpha

p$beta
p2$beta


    
p$sig      
p2$sigma     
p2$sigma_noar

p$Smsy
p2$Smsy

p$umsy
p2$umsy
           
    
p$rho   
p2$rho

pac<-ricker_TMB(data=harck, AC=TRUE)
pac2<-ricker_TMB2(data=harck, AC=TRUE)



pac$logalpha   
pac2$logalpha

pac$beta
pac2$beta

pac$Smax
pac2$Smax
    
pac$sigar      
pac2$sigma     
pac$sig      
pac2$sigma_noar

pac$Smsy
pac2$Smsy

pac$umsy
pac2$umsy
           
    
pac$rho   
pac2$rho




ptva <- ricker_rw_TMB(data=harck,tv.par="a", Smax_mean=250000,Smax_sd=200000)

ptva2 <- ricker_rw_TMB_logb(data=harck,tv.par="a")

ptva$logalpha
ptva2$logalpha

ptva$Smax
ptva2$Smax



ptvb <- ricker_rw_TMB(data=harck,tv.par="b",sig_p_sd=1)
ptvb2 <- ricker_rw_TMB_logb(data=harck,tv.par="b",sig_p_sd=1)

ptvb$logalpha
ptvb2$logalpha

ptvb$Smax
ptvb2$Smax


ptvab <- ricker_rw_TMB(data=harck,tv.par="both")
ptvab2 <- ricker_rw_TMB_logb(data=harck,tv.par="both")

ptvab$logalpha
ptvab2$logalpha

ptvab$Smax
ptvab2$Smax

ptvab$tmb_obj$report()


################################################
#hmmm


phmm <- ricker_hmm_TMB(data=harck, tv.par='both')
phmm_logb <- ricker_hmm_TMB2_logb(data=harck, tv.par='both')
phmm2 <- ricker_hmm_TMB2(data=harck, tv.par='both')


phmm$logalpha
phmm2$logalpha
phmm_logb$logalpha

phmm$beta
phmm2$beta
phmm_logb$beta  

phmm$Smax
phmm2$Smax 
phmm_logb$Smax   

phmm$sigma    
phmm2$sigma     

phmm$qij  
phmm2$qij   

phmm$Smsy  
phmm2$Smsy    

phmm$Smax
phmm2$Smax      

phmm$umsy     
phmm2$umsy

phmm$probregime-phmm2$probregime 

phmm$regime 
phmm2$regime


phmma <- ricker_hmm_TMB(data=harck, tv.par='a')
phmma2 <- ricker_hmm_TMB2(data=harck, tv.par='a')
phmma_logb <- ricker_hmm_TMB2_logb(data=harck, tv.par='a')

phmma$logalpha
phmma2$logalpha
$logalpha



phmma$beta
phmma2$beta  
phmma_logb$beta
 
phmma$Smax
phmma2$Smax      
phmma_logb$Smax    

phmma$sigma    
phmma2$sigma     

phmma$qij  -phmma2$qij   

phmma$Smsy  
phmma2$Smsy    

phmma$Smax
phmma2$Smax      

phmma$umsy     
phmma2$umsy

phmma$probregime-phmma2$probregime 

phmma$regime 
phmma2$regime


phmmb <- ricker_hmm_TMB(data=harck, tv.par='b',priors_flag =1)
phmmb2 <- ricker_hmm_TMB2(data=harck, tv.par='b',priors_flag = 1,Smax_limits=c(100,10000000))
phmmb_logb <- ricker_hmm_TMB2_logb(data=harck, tv.par='b',priors_flag =1)

phmmb$conv_problem
phmmb2$conv_problem


phmmb$logalpha
phmmb2$logalpha
phmmb_logb$logalpha

phmmb$beta
phmmb2$beta
phmmb_logb$beta


phmmb$Smax
phmmb2$Smax      
phmmb_logb$Smax 

phmmb$sigma    
phmmb2$sigma     

phmmb$qij  -phmmb2$qij   

phmmb$Smsy  
phmmb2$Smsy    

phmmb$umsy     
phmmb2$umsy

phmmb$probregime
phmmb2$probregime 

phmmb$regime -phmmb2$regime

phmmb$tmb_obj$report()$pnll
phmmb2$tmb_obj$report()$pnll

phmmb$tmb_obj$report()$nll
phmmb2$tmb_obj$report()$nll

#===============================
#old tests
p <- ricker_TMB(data=harck)
pac<-ricker_TMB(data=harck, AC=TRUE)
ptva<- ricker_rw_TMB(data=harck,tv.par="a")
ptvb <- ricker_rw_TMB(data=harck,tv.par="b",sig_p_sd=1)
ptvab <- ricker_rw_TMB(data=harck,tv.par="both",sig_p_sd=.5)
phmma <- ricker_hmm_TMB2(data=harck, tv.par='a')
phmmb <- ricker_hmm_TMB2(data=harck, tv.par='b')
phmm <- ricker_hmm_TMB2(data=harck, tv.par='both')


names(p)
names(pac)
names(ptva)
names(ptvb)
names(ptvab)
names(phmma)
names(phmmb)
names(phmm)




#reference points
smsyCalc(p$loga,p$beta)
umsyCalc(p$loga)
sgenCalcDirect(p$loga,p$beta)


#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)
#library(samEst)

#read in data 
#use Harrison as an example


#sr <- read.csv("C:/Users/worc/Documents/timevarproject/simeval/data/samsimHarCk/HARSR.csv")

#head(sr)


#har<-data.frame(by=sr$Brood.Year,
#	S=sr$Sum.Total.Spawners,
#	R=sr$AEQ_Recruitment..age.2.5.,
#	logRS=log(sr$AEQ_Recruitment..age.2.5./sr$Sum.Total.Spawners))


#harck<-har[!is.na(har$S),]

#setwd("C:\\Users\\worc\\Documents\\timevarproject\\samEst\\data")
#usethis::use_data(harck)
#library(samEst)
data(harck)
plot(harck$S,harck$R)


##testing functions

library(samEst)
sgen<-samEst::sGenCalc(p$alpha,p$beta,p$Smsy)



p$Smax
p$alpha

fullLL<-ricker_kf_TMB(data=harck,fullLL=TRUE)
reLL<-ricker_kf_TMB(data=harck)

fullLL$tmb_obj$fn()
reLL$tmb_obj$fn()

2*fullLL$tmb_obj$fn()+2*3

2*reLL$tmb_obj$fn()+2*3

pkf<-ricker_TMB(data=harck, AC=TRUE)


ip_logb_mean<-log(1/(max(harck$S)*.5))
ip_logb_sd<-sqrt(log(1+1))

p_ip <- ricker_TMB(data=harck,logb_p_mean=ip_logb_mean,logb_p_sd=ip_logb_sd)
p_ip$Smax
p_ip$alpha

simple_mod <- samEst::compile_code(type='static', ac=FALSE, par='n',lambertW = FALSE)
b <- ricker_stan(data=harck,iter = 800, mod=simple_mod)

pnp <- ricker_TMB(data=harck,prior=0)

pb <- ricker_stan(data=harck,iter = 2000)

mymod <- compile_code(type='static',ac=TRUE,par='n',caphigh=FALSE)
pb2 <- ricker_stan(data=harck,iter = 2000,AC=TRUE, mod = mymod)

names(pb)
p$alpha
pb$alpha

p$beta
pb$beta


pac<-ricker_TMB(data=harck, AC=TRUE)
pac[1:10]

lfostatic<-tmb_mod_lfo_cv(data=harck,model='static')
lfoac <- tmb_mod_lfo_cv(data=harck,model='staticAC')


sum(lfostatic)
sum(lfoac)
names(p)

ptva<- ricker_rw_TMB(data=harck,tv.par="a",sig_p_sd=1,AICc_type="conditional")
ptva$alpha
ptva$grad_i
ptvamaic<- ricker_rw_TMB(data=harck,tv.par="a",sig_p_sd=1, AICc_type="marginal")

ptvath_centered<- ricker_rw_TMB_centered(data=harck,tv.par="a",sig_p_sd=1,AICc_type="conditional")
ptvath_centered$logalpha
ptvath_centered$beta
ptvath_centered$grad_i

ptvbth_centered<- ricker_rw_TMB_centered(data=harck,tv.par="b",sig_p_sd=1,AICc_type="conditional")
ptvbth_centered$logalpha
ptvbth_centered$beta
ptvbth_centered$grad_i
names(ptvath_centered)




data.frame(
  value=c(p$AICc
ptvamaic$AICc
ptva$AICc), 

ptva_ip<- ricker_rw_TMB(data=harck,tv.par="a",sig_p_sd=1,logb_p_mean=ip_logb_mean,logb_p_sd=ip_logb_sd)

ptva[1:6]
ptva_ip[1:6]
pac$tmb_obj$report()$nll

pkfa<- ricker_kf_TMB(data=harck)
pkfa[1:6]

pkfa2<- ricker_kf_TMB(data=harck,fullLL=T)
pkfa2[1:6]

lfoalpha <- tmb_mod_lfo_cv(data=harck,tv.par='a', siglfo="obs")
sum(lfoalpha$lastparam)
sum(lfoalpha$last3paramavg)
sum(lfoalpha$last5paramavg)


phmm <- ricker_hmm_TMB(data=harck, tv.par='both')
phmm[1:10]

dirpr<-matrix(c(4,1,1,4),2,2)
phmm_dirpr <- ricker_hmm_TMB(data=harck, tv.par='both',dirichlet_prior=dirpr)
phmm_dirpr[1:10]


pbhmm <- ricker_hmm_stan(data, par='b')



phmmb <- ricker_hmm_TMB(data=harck, tv.par='b')
phmmb[1:8]


ptvb <- ricker_rw_TMB(data=harck,tv.par="b")

ptvab <- ricker_rw_TMB(data=harck,tv.par="both")



lfohmm <- tmb_mod_lfo_cv(data=harck,tv.par='HMM')

sum(lfohmm$regime_pick,na.rm=T)
sum(lfohmm$regime_average)


phmma <- ricker_hmm_TMB(data=harck, tv.par='a')

phmma <- ricker_hmm_TMB(data=harck, tv.par='a')


phmma[1:5]

lfohmma <- tmb_mod_lfo_cv(data=harck,tv.par='HMM_a')

sum(lfohmma$regime_pick)
sum(lfohmma$regime_average)

lfohmmb <- tmb_mod_lfo_cv(data=harck,tv.par='HMM_b')

sum(lfohmmb$regime_pick,na.rm=T)
sum(lfohmmb$regime_average)


phmmb <- ricker_hmm_TMB(data=harck, tv.par='b')
phmmb[1:5]

ptva$EDF
ptvb$EDF
ptvab$EDF


#stan functions

#lfo tmb testing
lfostatic<-tmb_mod_lfo_cv(data=harck,model='static', L=round((2/3)*nrow(harck)))
lfoac <- tmb_mod_lfo_cv(data=harck,model='staticAC', L=round((2/3)*nrow(harck)))
lfoalpha <- tmb_mod_lfo_cv(data=harck,model='rw_a', siglfo="obs", L=round((2/3)*nrow(harck)))
lfobeta <- tmb_mod_lfo_cv(data=harck,model='rw_b', siglfo="obs", L=round((2/3)*nrow(harck)))
lfoalphabeta <- tmb_mod_lfo_cv(data=harck,model='rw_both', siglfo="obs", L=round((2/3)*nrow(harck)))
lfohmma <- tmb_mod_lfo_cv(data=harck,model='HMM_a', L=round((2/3)*nrow(harck)))
lfohmmb <- tmb_mod_lfo_cv(data=harck,model='HMM_b', L=round((2/3)*nrow(harck)))



p <- ricker_TMB(data=harck)
pac<-ricker_TMB(data=harck, AC=TRUE)
ptva<- ricker_rw_TMB(data=harck,tv.par="a")
ptvb <- ricker_rw_TMB(data=harck,tv.par="b",sig_p_sd=1)
ptvab <- ricker_rw_TMB(data=harck,tv.par="both",sig_p_sd=.5)
phmma <- ricker_hmm_TMB(data=harck, tv.par='a')
phmmb <- ricker_hmm_TMB(data=harck, tv.par='b')
phmm <- ricker_hmm_TMB(data=harck, tv.par='both')


phmma$tmb_obj$report()$nll+
phmma$tmb_obj$report()$pnll

phmma$tmb_obj$fn()

c(
p$AICc,
pac$AICc,
ptva$AICc,
ptvb$AICc, 
ptvab$AICc,
phmma$AICc,
phmmb$AICc,
phmm$AICc)

c(
p$BIC,
pac$BIC,
ptva$BIC,
ptvb$BIC, 
ptvab$BIC,
phmma$BIC,
phmmb$BIC,
phmm$BIC)




#=====================================================
#RTMB


devtools::document()
devtools::load_all()
#library(RTMB)


p <- ricker_TMB(data=harck)
p2 <- ricker_RTMB(data=harck)

p$alpha
p2$logalpha

pac <- ricker_TMB(data=harck,AC=TRUE,priors_flag=0)
pac$rho
pac$alpha
pac$beta
pac$Smax
pac$tmb_obj$fn()
names(pac)

devtools::document()
devtools::load_all()
pac2 <- ricker_RTMB(data=harck,AC=TRUE,priors_flag=0)
pac2$rho
pac2$logalpha
pac2$beta
pac2$Smax
pac2$sigAR

pac2$obj$fn()


data=harck
dat<- list(
    obs_S = as.numeric(data$S),
    obs_logRS = data$logRS,
    by = data$by,
    priors_flag=1,
    stan_flag=0,
    sig_p_sd=1,
    logb_p_mean=-12,
    logb_p_sd=3
  )
  

  magS <- log10_ceiling(max(data$S))
  initlm <- lm(obs_logRS~obs_S, data=dat)

 
param <- list(
    logalpha   = initlm$coefficients[[1]],
    logbeta = ifelse(initlm$coefficients[[2]]>0,log(magS),log(-initlm$coefficients[[2]])),
    logsigobs = log(1),
    ar1_phi=3.0)


library(RTMB)
ricker_ac_RTMB_fn(param,dat)


tmbfn<-function(param){ricker_ac_RTMB_fn(param,dat=dat)}
obj <- MakeADFun(tmbfn, param, silent=TRUE)


opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)


p$alpha
obj$report()$logalpha

p$Smax
obj$report()$Smax



 final_grads <- sdr$gradient.fixed
  bad_eig <- FALSE
  conv_problem<-FALSE
  if (!is.null(sdr$pdHess)) {
    if (!sdr$pdHess) {
      warning("The model may not have converged: ",
        "non-positive-definite Hessian matrix.", call. = FALSE)
      conv_problem<-TRUE
    } else {
      eigval <- try(1 / eigen(sdr$cov.fixed)$values, silent = TRUE)
      if (is(eigval, "try-error") || (min(eigval) < .Machine$double.eps * 10)) {
        warning("The model may not have converged: ",
          "extreme or very small eigen values detected.", call. = FALSE)
        bad_eig <- TRUE
        conv_problem<-TRUE
      }
      if (any(final_grads > 0.01)){
        warning("The model may not have converged. ",
          "Maximum final gradient: ", max(final_grads), ".", call. = FALSE)
        conv_problem<-TRUE



(Smsy <- (1 - gsl::lambert_W0(exp(1 - p$alpha))) /p$beta)
(Smsy2 <- (1 - samest_lambertW(exp(1 - p$alpha))) /p$beta)




pl <- as.list(sdr, "Est")
plsd <- as.list(sdr, "Std")



p<-ricker_TMB(data=df,logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)

ptva <- ricker_rw_TMB(data=df,tv.par='a',logb_p_mean=logbeta_pr,
    logb_p_sd=logbeta_pr_sig, deltaEDF=0.01)

ptvb <- ricker_rw_TMB(data=df, tv.par='b',sigb_p_sd=1,
  logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, deltaEDF=0.01)



p$AICc
ptva$AICc
ptvb$AICc

ptva$EDF
ptvb$EDF





