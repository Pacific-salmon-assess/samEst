#============================================
#commands used for packages development and testing
#Catarina wor
#August 2022
#============================================



devtools::document()
devtools::load_all()


#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)


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

p <- ricker_TMB(data=harck)
p$Smax
p$alpha


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


ptva<- ricker_rw_TMB(data=harck,tv.par="a",sig_p_sd=1)
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




#stan functions

#lfo tmb testing
lfostatic<-tmb_mod_lfo_cv(data=harck,model='static', L=round((2/3)*nrow(harck)))
lfoac <- tmb_mod_lfo_cv(data=harck,model='staticAC', L=round((2/3)*nrow(harck)))
lfoalpha <- tmb_mod_lfo_cv(data=harck,model='rw_a', siglfo="obs", L=round((2/3)*nrow(harck)))
lfobeta <- tmb_mod_lfo_cv(data=harck,model='rw_b', siglfo="obs", L=round((2/3)*nrow(harck)))
lfoalphabeta <- tmb_mod_lfo_cv(data=harck,model='rw_both', siglfo="obs", L=round((2/3)*nrow(harck)))
lfohmma <- tmb_mod_lfo_cv(data=harck,model='HMM_a', L=round((2/3)*nrow(harck)))
lfohmmb <- tmb_mod_lfo_cv(data=harck,model='HMM_b', L=round((2/3)*nrow(harck)))
lfohmm <- tmb_mod_lfo_cv(data=harck,model='HMM', L=round((2/3)*nrow(harck)))


p <- ricker_TMB(data=harck)
pac<-ricker_TMB(data=harck, AC=TRUE)
ptva<- ricker_rw_TMB(data=harck,tv.par="a")
ptvb <- ricker_rw_TMB(data=harck,tv.par="b",sig_p_sd=1)
ptvab <- ricker_rw_TMB(data=harck,tv.par="both",sig_p_sd=.5)
phmma <- ricker_hmm_TMB(data=harck, tv.par='a')
phmmb <- ricker_hmm_TMB(data=harck, tv.par='b')
phmm <- ricker_hmm_TMB(data=harck, tv.par='both')




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